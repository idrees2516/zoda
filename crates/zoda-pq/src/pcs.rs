//! The lattice polynomial commitment scheme: BDLOP Module-LWE commitments
//! + Lyubashevsky-style linear openings with **narrow ring challenges**
//! (coefficients in {−1, 0, 1}, bounded weight), made non-interactive via
//! Fiat–Shamir.
//!
//! * Commit: `C = A·r + B·m ∈ R_q^{k̄}` (r short CBD, m = coefficient
//!   chunks). Hiding: Module-LWE. Binding: Module-SIS.
//! * Open at a ring point `ζ` (a scalar `z` embeds as the constant ring
//!   element, recovering the classical evaluation): the prover answers
//!   with `(C₁, v₁, r₂, m₂)` where `r₂ = r₁ + x·r`, `m₂ = m₁ + x·m` for
//!   the Fiat–Shamir challenge `x`, rejection-sampled so `r₂` stays short.
//! * Verify: `A·r₂ + B·m₂ = C₁ + x·C` and `eval(m₂, ζ) = v₁ + x·v`.
//!
//! The evaluation functional `m ↦ Σⱼ mⱼ(ζ)·ζ^{jn}` is R-linear, which is
//! exactly what makes the second verification equation valid.

use crate::params::PqParams;
use crate::rq::{fq_to_poly, Rq, RqVector};
use zoda_math::sha256::sha256;
use zoda_math::{Fq, PrimeField, ZodaRng};

/// Challenge weight (number of ±1 coefficients).
pub const CHALLENGE_WEIGHT: usize = 32;
/// Rejection-sampling infinity-norm bound for `r₂` (η·ω + margin).
const BETA_MARGIN: u32 = 64;

/// A public commitment.
#[derive(Clone, PartialEq)]
pub struct LatticeCommitment {
    /// `C = A·r + B·m`
    pub c: Vec<Rq>,
    /// Masking commitment `C₁ = A·r₁ + B·m₁` (published up-front).
    pub c1: Vec<Rq>,
}

/// The opening proof.
#[derive(Clone, PartialEq)]
pub struct LatticeProof {
    pub c1: Vec<Rq>,
    pub v1: Rq,
    pub r2: RqVector,
    pub m2: RqVector,
}

/// Everything the prover needs to open (carried opaquely).
#[derive(Clone)]
pub struct OpeningKey {
    r: RqVector,
    m: Vec<Rq>,
    /// the committed polynomial (for re-derivation)
    pub poly: Vec<Fq>,
}

/// Deterministic expansion of a public matrix from the seed.
fn expand_matrix(seed: &[u8; 32], domain: u8, rows: usize, cols: usize, n: usize) -> Vec<Vec<Rq>> {
    let mut out = Vec::with_capacity(rows);
    for i in 0..rows {
        let mut row = Vec::with_capacity(cols);
        for j in 0..cols {
            let mut coeffs = Vec::with_capacity(n);
            let mut counter = 0u32;
            while coeffs.len() < n {
                let mut input = Vec::with_capacity(32 + 8);
                input.extend_from_slice(seed);
                input.push(domain);
                input.extend_from_slice(&(i as u16).to_le_bytes());
                input.extend_from_slice(&(j as u16).to_le_bytes());
                input.extend_from_slice(&counter.to_le_bytes());
                for chunk in sha256(&input).chunks_exact(3) {
                    let v = chunk[0] as u32 | (chunk[1] as u32) << 8 | (chunk[2] as u32) << 16;
                    if v < crate::rq::Q {
                        coeffs.push(Fq(v));
                        if coeffs.len() == n {
                            break;
                        }
                    }
                }
                counter += 1;
            }
            row.push(Rq(coeffs));
        }
        out.push(row);
    }
    out
}

fn cbd(n: usize, eta: u32, rng: &mut ZodaRng) -> Rq {
    let mut coeffs = Vec::with_capacity(n);
    for _ in 0..n {
        let mut v = 0i64;
        for _ in 0..eta {
            v += rng.next_below(2) as i64 - rng.next_below(2) as i64;
        }
        let v = if v < 0 {
            Fq((crate::rq::Q as i64 + v) as u32)
        } else {
            Fq(v as u32)
        };
        coeffs.push(v);
    }
    Rq(coeffs)
}

fn sample_short(k: usize, n: usize, eta: u32, rng: &mut ZodaRng) -> RqVector {
    RqVector((0..k).map(|_| cbd(n, eta, rng)).collect())
}

fn inf_norm(v: &RqVector) -> u32 {
    let mut m = 0u32;
    for r in &v.0 {
        for c in &r.0 {
            let d = c.0.min(crate::rq::Q - c.0);
            m = m.max(d);
        }
    }
    m
}

fn mat_vec_mul(a: &[Vec<Rq>], v: &[Rq]) -> Vec<Rq> {
    a.iter()
        .map(|row| {
            let mut acc = Rq(vec![Fq(0); v[0].len()]);
            for (ai, vi) in row.iter().zip(v.iter()) {
                acc = acc.add(&ai.mul(vi));
            }
            acc
        })
        .collect()
}

fn public_seed() -> [u8; 32] {
    sha256(b"zoda-pq-lattice-pcs-public-matrix-seed-v1")
}

/// Sample a narrow challenge: coefficients in {−1,0,1}, weight ≤ ω.
fn sample_challenge(n: usize, transcript: &[u8]) -> Rq {
    let mut coeffs = vec![Fq(0); n];
    let mut placed = 0usize;
    let mut counter = 0u32;
    while placed < CHALLENGE_WEIGHT {
        let mut input = Vec::with_capacity(transcript.len() + 4);
        input.extend_from_slice(transcript);
        input.extend_from_slice(&counter.to_le_bytes());
        let h = sha256(&input);
        // 2 bytes: position, 1 byte: sign
        for chunk in h.chunks_exact(3) {
            if placed >= CHALLENGE_WEIGHT {
                break;
            }
            let pos = (chunk[0] as u32 | (chunk[1] as u32) << 8) as usize % n;
            if coeffs[pos].is_zero() {
                coeffs[pos] = if chunk[2] & 1 == 0 {
                    Fq(1)
                } else {
                    Fq(crate::rq::Q - 1)
                };
                placed += 1;
            }
        }
        counter += 1;
        if counter > 128 {
            break; // safety valve (extremely unlikely)
        }
    }
    Rq(coeffs)
}

/// Evaluate the chunked polynomial at a ring point ζ:
/// `P(ζ) = Σⱼ mⱼ(ζ)·ζ^{jn}` (Horner over chunks with the running power
/// ζ^n).
pub fn eval_chunks(m: &[Rq], zeta: &Rq) -> Rq {
    let n = zeta.len();
    // zeta^n by repeated squaring (n is a power of two)
    let mut zpow = zeta.clone();
    let mut p = n;
    while p > 1 {
        zpow = zpow.mul(&zpow);
        p /= 2;
    }
    let eval_ring = |p: &Rq| -> Rq {
        let mut acc = Rq(vec![Fq(0); n]);
        for c in p.0.iter().rev() {
            let mut term = zeta.mul(&acc);
            term.0[0] = term.0[0] + *c;
            acc = term;
        }
        acc
    };
    let mut it = m.iter().rev();
    let mut acc = eval_ring(it.next().expect("nonempty"));
    for chunk in it {
        acc = acc.mul(&zpow);
        acc = acc.add(&eval_ring(chunk));
    }
    acc
}

/// An admissible evaluation point: ζ = t·X with t an n-th root of unity
/// mod q, so that ζ^n = −1 and the evaluation functional is R-linear
/// (a ring homomorphism on R_q — required by the σ-protocol's second
/// verification equation).
pub fn ring_point(t: u64, n: usize) -> Rq {
    let mut c = vec![Fq(0); n];
    c[1] = Fq((t % crate::rq::Q as u64) as u32);
    Rq(c)
}

/// A random admissible ring point from a seed: t = g^((q-1)/n · k).
pub fn ring_point_from_seed(seed: &[u8], n: usize) -> Rq {
    let h = sha256(seed);
    let k = u64::from_le_bytes(h[..8].try_into().unwrap()) % n as u64;
    let exp = ((crate::rq::Q as u64 - 1) / n as u64) * k;
    let mut t: u64 = 1;
    let mut base: u64 = 10;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            t = (t * base) % crate::rq::Q as u64;
        }
        base = (base * base) % crate::rq::Q as u64;
        e >>= 1;
    }
    ring_point(t, n)
}

/// Verify ζ^n == −1 (admissibility).
pub fn is_admissible_ring_point(zeta: &Rq) -> bool {
    let n = zeta.len();
    let mut zpow = zeta.clone();
    let mut p = n;
    while p > 1 {
        zpow = zpow.mul(&zpow);
        p /= 2;
    }
    let mut expect = vec![Fq(0); n];
    expect[0] = -Fq(1);
    zpow == Rq(expect)
}

pub struct LatticePcs;

impl LatticePcs {
    /// Commit to a coefficient-form polynomial over Z_q.
    pub fn commit(params: &PqParams, poly: &[Fq], rng: &mut ZodaRng) -> (LatticeCommitment, OpeningKey) {
        let n = params.n;
        let seed = public_seed();
        let m = fq_to_poly(poly, n);
        let k = m.len();
        let a = expand_matrix(&seed, 0xA1, params.k_bar, params.m_bar, n);
        let b = expand_matrix(&seed, 0xB2, params.k_bar, k, n);

        let r = sample_short(params.m_bar, n, params.eta, rng);
        let ar = mat_vec_mul(&a, &r.0);
        let bm = mat_vec_mul(&b, &m);
        let c: Vec<Rq> = ar.iter().zip(bm.iter()).map(|(x, y)| x.add(y)).collect();

        let m1 = RqVector::uniform(k, n, rng);
        let r1 = sample_short(params.m_bar, n, params.eta, rng);
        let ar1 = mat_vec_mul(&a, &r1.0);
        let bm1 = mat_vec_mul(&b, &m1.0);
        let c1: Vec<Rq> = ar1.iter().zip(bm1.iter()).map(|(x, y)| x.add(y)).collect();

        (
            LatticeCommitment { c, c1 },
            OpeningKey { r, m, poly: poly.to_vec() },
        )
    }

    /// Plain scalar evaluation of the polynomial (helper, mod q).
    pub fn evaluate(poly: &[Fq], z: u64) -> u32 {
        let mut acc: u64 = 0;
        for c in poly.iter().rev() {
            acc = (acc * (z % crate::rq::Q as u64) + c.0 as u64) % crate::rq::Q as u64;
        }
        acc as u32
    }

    /// The ring-element value of the committed polynomial at the ring
    /// point ζ (helper used by tests and the docs).
    pub fn eval_at_ring(poly: &[Fq], zeta: &Rq) -> Rq {
        let n = zeta.len();
        eval_chunks(&fq_to_poly(poly, n), zeta)
    }

    /// Fiat–Shamir transcript for a statement (c, c1, ζ, v, and — for the
    /// prover — nothing secret).
    fn transcript(c: &[Rq], c1: &[Rq], zeta: &Rq, v: &Rq) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"zoda-pq-fs-v1");
        for ri in c {
            buf.extend_from_slice(&ri.to_bytes());
        }
        for ri in c1 {
            buf.extend_from_slice(&ri.to_bytes());
        }
        buf.extend_from_slice(&zeta.to_bytes());
        buf.extend_from_slice(&v.to_bytes());
        buf
    }

    /// Produce a non-interactive opening proof at the ring point ζ.
    pub fn open(
        params: &PqParams,
        key: &OpeningKey,
        zeta: &Rq,
        v: &Rq,
        rng: &mut ZodaRng,
    ) -> LatticeProof {
        let n = params.n;
        let k = key.m.len();
        let seed = public_seed();
        let a = expand_matrix(&seed, 0xA1, params.k_bar, params.m_bar, n);
        let b = expand_matrix(&seed, 0xB2, params.k_bar, k, n);

        // Recompute C from the key so the transcript binds the real
        // commitment value.
        let ar = mat_vec_mul(&a, &key.r.0);
        let bm = mat_vec_mul(&b, &key.m);
        let c: Vec<Rq> = ar.iter().zip(bm.iter()).map(|(x, y)| x.add(y)).collect();

        // Masking material, resampled by rejection sampling until the
        // response r2 = r1 + x·r is short enough.
        let beta = params.eta * CHALLENGE_WEIGHT as u32 + BETA_MARGIN;
        let (r2, m1, c1, v1, x) = loop {
            let m1 = RqVector::uniform(k, n, rng);
            let r1 = sample_short(params.m_bar, n, params.eta, rng);
            let ar1 = mat_vec_mul(&a, &r1.0);
            let bm1 = mat_vec_mul(&b, &m1.0);
            let c1: Vec<Rq> = ar1
                .iter()
                .zip(bm1.iter())
                .map(|(p, q)| p.add(q))
                .collect();
            let v1 = eval_chunks(&m1.0, zeta);
            let tr = Self::transcript(&c, &c1, zeta, v);
            let x = sample_challenge(n, &tr);
            let r2 = RqVector(
                r1.0
                    .iter()
                    .zip(key.r.0.iter())
                    .map(|(r1i, ri)| r1i.add(&x.mul(ri)))
                    .collect(),
            );
            if inf_norm(&r2) <= beta {
                break (r2, m1, c1, v1, x);
            }
        };

        // m2 = m1 + x·m
        let m2 = RqVector(
            m1.0
                .iter()
                .zip(key.m.iter())
                .map(|(m1i, mi)| m1i.add(&x.mul(mi)))
                .collect(),
        );

        LatticeProof { c1, v1, r2, m2 }
    }

    /// Verify an opening at the ring point ζ.
    pub fn verify(
        params: &PqParams,
        com: &LatticeCommitment,
        zeta: &Rq,
        v: &Rq,
        proof: &LatticeProof,
    ) -> bool {
        let n = params.n;
        let m_bar = params.m_bar;
        let k = proof.m2.len();
        if proof.r2.len() != m_bar || com.c.len() != params.k_bar {
            return false;
        }
        let seed = public_seed();
        let a = expand_matrix(&seed, 0xA1, params.k_bar, m_bar, n);
        let b = expand_matrix(&seed, 0xB2, params.k_bar, k, n);

        // strong Fiat–Shamir: bind to (C, C1, ζ, v)
        let tr = Self::transcript(&com.c, &proof.c1, zeta, v);
        let x = sample_challenge(n, &tr);

        // Equation 1: A·r2 + B·m2 == C1 + x·C
        let ar2 = mat_vec_mul(&a, &proof.r2.0);
        let bm2 = mat_vec_mul(&b, &proof.m2.0);
        let lhs: Vec<Rq> = ar2.iter().zip(bm2.iter()).map(|(p, q)| p.add(q)).collect();
        let rhs: Vec<Rq> = proof
            .c1
            .iter()
            .zip(com.c.iter())
            .map(|(c1i, ci)| c1i.add(&x.mul(ci)))
            .collect();
        if lhs != rhs {
            return false;
        }

        // shortness of r2 (implicit MSIS binding)
        let beta = params.eta * CHALLENGE_WEIGHT as u32 + BETA_MARGIN;
        if inf_norm(&proof.r2) > beta {
            return false;
        }

        // Equation 2: eval(m2, ζ) == v1 + x(ζ)·v
        // (the evaluation functional is an R-module map only up to the
        // character x ↦ x(ζ): L(x·m) = x(ζ)·L(m) for admissible ζ)
        let lhs2 = eval_chunks(&proof.m2.0, zeta);
        let x_at_zeta = {
            // evaluate the challenge ring element at ζ
            let mut acc = Rq(vec![Fq(0); n]);
            for c in x.0.iter().rev() {
                let mut term = zeta.mul(&acc);
                term.0[0] = term.0[0] + *c;
                acc = term;
            }
            acc
        };
        let rhs2 = proof.v1.add(&x_at_zeta.mul(v));
        lhs2 == rhs2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expand_matrix_deterministic() {
        let a = expand_matrix(&public_seed(), 0xA1, 2, 2, 8);
        let b = expand_matrix(&public_seed(), 0xA1, 2, 2, 8);
        assert_eq!(a, b);
        let c = expand_matrix(&public_seed(), 0xB2, 2, 2, 8);
        assert_ne!(a, c);
    }

    #[test]
    fn challenge_is_narrow() {
        let x = sample_challenge(64, b"transcript-test");
        let w: usize = x
            .0
            .iter()
            .filter(|c| c.0 == 1 || c.0 == crate::rq::Q - 1)
            .count();
        assert!(w <= CHALLENGE_WEIGHT && w > 0, "weight {}", w);
        // deterministic
        let y = sample_challenge(64, b"transcript-test");
        assert_eq!(x, y);
    }

    #[test]
    fn ring_points_are_admissible() {
        let zeta = ring_point_from_seed(b"test-ring-point", 64);
        assert!(is_admissible_ring_point(&zeta));
        assert!(is_admissible_ring_point(&ring_point(1, 64)));
        let mut c = vec![Fq(0); 64];
        c[0] = Fq(5);
        assert!(!is_admissible_ring_point(&Rq(c)));
    }

    #[test]
    fn evaluation_is_ring_hom_at_admissible_points() {
        let mut rng = ZodaRng::from_seed(*b"pq-eval-test-seed-00000000000000");
        let n = 32;
        let zeta = ring_point_from_seed(b"hom-test", n);
        let a = Rq((0..n).map(|_| Fq(rng.next_below(crate::rq::Q as u64) as u32)).collect());
        let b = Rq((0..n).map(|_| Fq(rng.next_below(crate::rq::Q as u64) as u32)).collect());
        let ab = a.mul(&b);
        let lhs = eval_chunks(&[ab], &zeta);
        let rhs = eval_chunks(&[a], &zeta).mul(&eval_chunks(&[b], &zeta));
        assert_eq!(lhs, rhs);
    }
}
