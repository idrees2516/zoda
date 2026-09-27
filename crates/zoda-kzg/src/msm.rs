//! Pippenger multi-scalar multiplication over G1 with a simple
//! thread-parallel path for large batches.

use zoda_bls::g1::{G1Affine, G1Projective};
use zoda_bls::Fp;

/// Window size lookup for a given batch size (following the standard
/// Pippenger heuristic).
fn window(n: usize) -> usize {
    if n < 32 {
        8
    } else {
        (f64::ln(n as f64) / 0.693 + 2.0) as usize
    }
}

/// Single-threaded Pippenger MSM: Σ sᵢ·Pᵢ.
pub fn msm(points: &[G1Affine], scalars: &[zoda_math::Fr]) -> G1Projective {
    assert_eq!(points.len(), scalars.len());
    let n = points.len();
    if n == 0 {
        return G1Projective::identity();
    }
    let w = window(n).min(16);
    let bits = 255; // scalar field width (top bit always 0)
    let windows = (bits + w - 1) / w;

    // bucket accumulation per window
    let mut result = G1Projective::identity();
    for wi in (0..windows).rev() {
        // result = result * 2^w
        for _ in 0..w {
            result = result.double();
        }
        let buckets = 1usize << w;
        let mut acc = vec![G1Projective::identity(); buckets];
        for (pt, sc) in points.iter().zip(scalars.iter()) {
            let limbs = sc.to_repr();
            // extract the wi-th window
            let bit_lo = wi * w;
            let mut idx = 0usize;
            for b in 0..w {
                let bit = bit_lo + b;
                if bit >= 256 {
                    break;
                }
                let limb = bit / 64;
                let off = bit % 64;
                if (limbs[limb] >> off) & 1 == 1 {
                    idx |= 1 << b;
                }
            }
            if idx != 0 {
                acc[idx] = acc[idx].add_mixed(pt);
            }
        }
        // sum buckets: running = Σ_{j} j·bucket[j] via the standard trick
        let mut running = G1Projective::identity();
        let mut sum = G1Projective::identity();
        for j in (1..buckets).rev() {
            running = running.add(&acc[j]);
            sum = sum.add(&running);
        }
        result = result.add(&sum);
    }
    result
}

/// Multi-threaded MSM for large batches (falls back to single-threaded for
/// small ones). Uses std::thread::scope with static chunking.
pub fn msm_parallel(points: &[G1Affine], scalars: &[zoda_math::Fr]) -> G1Projective {
    let n = points.len();
    if n < 512 {
        return msm(points, scalars);
    }
    let cpus = std::thread::available_parallelism()
        .map(|c| c.get())
        .unwrap_or(1);
    let chunks = cpus.min(n / 256).max(1);
    let per = (n + chunks - 1) / chunks;
    let mut results = Vec::with_capacity(chunks);
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
            handles.push(scope.spawn(move || msm(pts, scs)));
        }
        for h in handles {
            results.push(h.join().expect("msm thread panicked"));
        }
    });
    let mut acc = G1Projective::identity();
    for r in results {
        acc = acc.add(&r);
    }
    acc
}

/// Weighted G1 sum Σ sᵢ·Pᵢ given points already in affine form (thin
/// wrapper used by the KZG batch verifier).
pub fn g1_lincomb(points: &[G1Affine], scalars: &[zoda_math::Fr]) -> G1Affine {
    msm_parallel(points, scalars).to_affine()
}

/// Convert an Fp point field into an affine point (validation helper).
pub fn fp_to_affine(_: &Fp) -> G1Affine {
    unreachable!("unused")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn msm_matches_naive() {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"msm-test-seed-000000000000000000");
        let g = G1Projective::generator();
        for n in [1usize, 2, 7, 33, 128, 300] {
            let pts: Vec<G1Affine> = (0..n)
                .map(|i| g.mul_limbs(&[(i * 7919 + 3) as u64, 0, 0, 0]).to_affine())
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
}
