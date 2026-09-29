//! Multi-scalar multiplication over G1 — Pippenger's bucket method with
//! the modern low-level treatment:
//!
//! * **Jacobian bucket accumulation** — `madd-2007-bl` (7M+4S) mixed
//!   additions; fresh buckets are affine copies with `Z = 1`.
//! * **Precomputed window digits** — every scalar is decomposed once into
//!   `c`-bit windows; the per-point `to_repr` cost is paid a single time.
//! * **Classic all-levels bucket reduction** — `Σ_j j·B_j` is computed by
//!   the running-prefix trick *at every bucket level*; the implicit
//!   multiplication by the bucket index comes precisely from those
//!   repeated additions, so empty levels must still advance the sum
//!   whenever the running prefix is non-identity.
//! * **Adaptive window schedule** — `c ≈ log₂n − log₂log₂n + 1` keeps the
//!   bucket space `2^c` on the order of the batch size, the classic
//!   Pippenger optimum (clamped to [4, 13]).
//! * **Persistent bucket arena** — the bucket array is thread-local and
//!   reused across calls; only touched entries are reset.
//! * **Thread-parallel large batches** — static chunking over the points,
//!   each chunk running the full window schedule.

use zoda_bls::curve::Jacobian;
use zoda_bls::g1::{G1Affine, G1Config, G1Projective};
use zoda_math::Fr;

pub type G1J = Jacobian<G1Config>;

/// Window size for a batch of `n` points: the classic Pippenger optimum.
/// `2^c` stays on the order of `n`, balancing accumulation (n adds per
/// window) against reduction (≈2·2^c adds per window).
fn window(n: usize) -> usize {
    if n < 16 {
        4
    } else {
        let l = (n as f64).log2();
        let c = (l - l.log2() + 1.0) as usize;
        c.clamp(4, 13)
    }
}

thread_local! {
    /// Persistent per-thread bucket arena. Invariant: every entry is the
    /// identity (Z = 0) whenever it is not being used by an active window —
    /// touched entries are reset at the end of each window.
    static BUCKETS: std::cell::RefCell<Vec<G1J>> =
        std::cell::RefCell::new(Vec::new());
}

/// Single-threaded Pippenger MSM: Σ sᵢ·Pᵢ.
pub fn msm(points: &[G1Affine], scalars: &[Fr]) -> G1Projective {
    msm_jacobian(points, scalars).to_homogeneous()
}

/// The MSM core, staying in Jacobian coordinates (for callers that keep
/// accumulating, e.g. the FK20 G1-FFT pipeline).
pub fn msm_jacobian(points: &[G1Affine], scalars: &[Fr]) -> G1J {
    assert_eq!(points.len(), scalars.len());
    let n = points.len();
    if n == 0 {
        return G1J::identity();
    }
    let c = window(n);
    let nwin = (256 + c - 1) / c;

    // ---- digit decomposition (window-major layout for locality) ----
    let mut digits = vec![0u16; nwin * n];
    for (i, sc) in scalars.iter().enumerate() {
        let r = sc.to_repr(); // canonical little-endian limbs
        for w in 0..nwin {
            let lo = w * c;
            let mut d: u32 = 0;
            for b in 0..c {
                let bit = lo + b;
                if bit >= 256 {
                    break;
                }
                if (r[bit / 64] >> (bit % 64)) & 1 == 1 {
                    d |= 1 << b;
                }
            }
            digits[w * n + i] = d as u16;
        }
    }

    BUCKETS.with(|arena| {
        let mut arena = arena.borrow_mut();
        let nbuckets = 1usize << c;
        if arena.len() < nbuckets {
            arena.resize(nbuckets, G1J::identity());
        }
        let buckets = &mut arena[..nbuckets];

        let mut result = G1J::identity();
        let mut touched: Vec<usize> = Vec::new();
        // process windows from the most significant down
        for w in (0..nwin).rev() {
            if !result.is_identity() {
                for _ in 0..c {
                    result = result.double();
                }
            }
            let dm = &digits[w * n..(w + 1) * n];
            touched.clear();
            for (i, &d) in dm.iter().enumerate() {
                if d == 0 {
                    continue;
                }
                if points[i].infinity {
                    continue; // the identity contributes nothing
                }
                let d = d as usize;
                let b = &mut buckets[d];
                if b.z.is_zero() {
                    // fresh bucket: copy the affine point with Z = 1 (free)
                    *b = G1J {
                        x: points[i].x,
                        y: points[i].y,
                        z: zoda_bls::Fp::one(),
                    };
                    touched.push(d);
                } else {
                    *b = b.add_mixed(&points[i]);
                }
            }
            if touched.is_empty() {
                continue;
            }
            // Classic bucket reduction: for j = 2^c-1 .. 1,
            //   running += B_j ;  sum += running
            // so Σ_j j·B_j emerges from the repeated prefix additions.
            // Levels where both the bucket and the running prefix are the
            // identity contribute nothing and are skipped.
            let mut running = G1J::identity();
            let mut sum = G1J::identity();
            for j in (1..nbuckets).rev() {
                let b = buckets[j];
                let bucket_empty = b.z.is_zero();
                if bucket_empty && running.is_identity() {
                    continue;
                }
                if !bucket_empty {
                    if running.is_identity() {
                        running = b;
                    } else if b.z == zoda_bls::Fp::one() {
                        // pristine bucket: (x, y) is directly affine
                        running = running.add_mixed(&G1Affine {
                            x: b.x,
                            y: b.y,
                            infinity: false,
                        });
                    } else {
                        running = running.add(&b);
                    }
                    buckets[j] = G1J::identity(); // reset for the next window
                }
                if !running.is_identity() {
                    if sum.is_identity() {
                        sum = running;
                    } else {
                        sum = sum.add(&running);
                    }
                }
            }
            result = if result.is_identity() {
                sum
            } else {
                result.add(&sum)
            };
        }
        result
    })
}

/// Multi-threaded MSM for large batches (falls back to single-threaded for
/// small ones). Uses std::thread::scope with static point chunking; each
/// chunk runs the full window schedule at full height.
pub fn msm_parallel(points: &[G1Affine], scalars: &[Fr]) -> G1Projective {
    msm_parallel_jacobian(points, scalars).to_homogeneous()
}

/// Parallel MSM core returning Jacobian.
pub fn msm_parallel_jacobian(points: &[G1Affine], scalars: &[Fr]) -> G1J {
    let n = points.len();
    if n < 512 {
        return msm_jacobian(points, scalars);
    }
    let cpus = std::thread::available_parallelism()
        .map(|c| c.get())
        .unwrap_or(1);
    let chunks = cpus.min(n / 256).max(1);
    if chunks <= 1 {
        return msm_jacobian(points, scalars);
    }
    let per = (n + chunks - 1) / chunks;
    let mut results: Vec<G1J> = Vec::with_capacity(chunks);
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(chunks);
        for ci in 0..chunks {
            let start = ci * per;
            if start >= n {
                break;
            }
            let end = (start + per).min(n);
            let pts = &points[start..end];
            let scs = &scalars[start..end];
            handles.push(scope.spawn(move || msm_jacobian(pts, scs)));
        }
        for h in handles {
            results.push(h.join().expect("msm thread panicked"));
        }
    });
    let mut acc = G1J::identity();
    for r in results {
        if acc.is_identity() {
            acc = r;
        } else if !r.is_identity() {
            acc = acc.add(&r);
        }
    }
    acc
}

/// Weighted G1 sum Σ sᵢ·Pᵢ given points already in affine form (thin
/// wrapper used by the KZG batch verifier).
pub fn g1_lincomb(points: &[G1Affine], scalars: &[Fr]) -> G1Affine {
    msm_parallel_jacobian(points, scalars).to_affine()
}

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_math::PrimeField;

    #[test]
    fn msm_matches_naive() {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"msm-test-seed-000000000000000000");
        let g = G1Projective::generator();
        for n in [1usize, 2, 3, 7, 33, 64, 128, 300, 1000, 4096] {
            let pts: Vec<G1Affine> = (0..n)
                .map(|i| {
                    g.mul_limbs(&[(i * 7919 + 3) as u64, i as u64 * 31, 7, 0])
                        .to_affine()
                })
                .collect();
            let scalars: Vec<zoda_math::Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
            let got = msm(&pts, &scalars);
            // naive
            let mut want = G1Projective::identity();
            for (p, s) in pts.iter().zip(scalars.iter()) {
                want = want.add(&p.to_projective().mul_fr(s));
            }
            assert_eq!(
                got.to_affine(),
                want.to_affine(),
                "msm mismatch for n={}",
                n
            );
        }
    }

    #[test]
    fn msm_edge_cases() {
        let g = G1Projective::generator();
        let p1 = g.mul_limbs(&[7, 0, 0, 0]).to_affine();
        let p2 = g.mul_limbs(&[11, 0, 0, 0]).to_affine();
        // zero-length
        assert!(msm(&[], &[]).is_identity());
        // zero scalar
        let zero = zoda_math::Fr::zero();
        assert!(msm(&[p1], &[zero]).is_identity());
        // all-ones small batch
        let one = zoda_math::Fr::ONE;
        let two = zoda_math::Fr::from_u64(2);
        let got = msm(&[p1, p2], &[one, one]).to_affine();
        let want = p1.to_projective().add(&p2.to_projective()).to_affine();
        assert_eq!(got, want);
        // duplicated points (bucket collision path): 1·P + 2·P + 1·Q = 3P + Q
        let got = msm(&[p1, p1, p2], &[one, two, one]).to_affine();
        let want = p1
            .to_projective()
            .add(&p1.to_projective())
            .add(&p1.to_projective())
            .add(&p2.to_projective())
            .to_affine();
        assert_eq!(got, want);
        // a point and its negation cancel
        let got = msm(&[p1, p1.neg()], &[one, one]).to_affine();
        assert!(got.is_identity());
        // mostly-zero scalars (the single-nonzero path)
        let mut rng = zoda_math::ZodaRng::from_seed(*b"msm-zero-scalar-0000000000000000");
        let n = 256;
        let pts: Vec<G1Affine> = (0..n)
            .map(|i| g.mul_limbs(&[i as u64, (i * 13) as u64, 3, 0]).to_affine())
            .collect();
        for nz in [1usize, 2, 5] {
            let mut scalars = vec![zoda_math::Fr::zero(); n];
            for s in scalars.iter_mut().take(nz) {
                *s = rng.next_fr(false);
            }
            let got = msm(&pts, &scalars).to_affine();
            let mut want = G1Projective::identity();
            for (p, s) in pts.iter().zip(scalars.iter()) {
                want = want.add(&p.to_projective().mul_fr(s));
            }
            assert_eq!(got, want.to_affine(), "nz={}", nz);
        }
        // parallel path agrees with the serial path
        let mut rng = zoda_math::ZodaRng::from_seed(*b"msm-parl-seed-000000000000000000");
        let n = 2048;
        let pts: Vec<G1Affine> = (0..n)
            .map(|i| g.mul_limbs(&[i as u64, i as u64 * 31, 0, 0]).to_affine())
            .collect();
        let scalars: Vec<zoda_math::Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
        assert_eq!(
            msm(&pts, &scalars).to_affine(),
            msm_parallel(&pts, &scalars).to_affine()
        );
        // repeated calls reuse the arena correctly (stale-bucket check)
        let again = msm(&pts, &scalars).to_affine();
        assert_eq!(again, msm(&pts, &scalars).to_affine());
    }
}
