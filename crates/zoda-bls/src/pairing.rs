//! The optimal ate pairing on BLS12-381.
//!
//! * Miller loop: Beuchat et al. algorithms 26/27 (eprint 2010/354),
//!   64 iterations over the BLS parameter x = -0xd201000000010000.
//! * Final exponentiation: the Fuentes–Castañeda–Rodríguez-Henríquez
//!   "fountain" chain (cyclotomic exponentiation + Frobenius), matching the
//!   structure used by zkcrypto/bls12_381 and blst.
//!
//! Frobenius coefficients are **derived and verified at initialisation**
//! from the tower definition (ξ = 1 + u), not transcribed from tables:
//! Γ₁[k] = ξ^((pᵏ−1)/3) and γ₁ = w^(p−1) composed by Fp6-Frobenius —
//! eliminating an entire class of transcription bugs.

use crate::fp::Fp;
use crate::fp2::Fp2;
use crate::fp6::Fp6;
use crate::fp12::Fp12;
use crate::g1::G1Affine;
use crate::g2::{G2Affine, G2Projective};
use std::sync::OnceLock;

/// The BLS12-381 parameter x (positive magnitude; the true x is negative).
pub const BLS_X: u64 = 0xd201_0000_0001_0000;
pub const BLS_X_IS_NEGATIVE: bool = true;

// ---------------------------------------------------------------------------
// Little-endian big-int helpers (used only during one-time initialisation)
// ---------------------------------------------------------------------------

fn mul_limbs_big(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64; a.len() + b.len()];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        let mut carry = 0u128;
        for (j, &y) in b.iter().enumerate() {
            let t = out[i + j] as u128 + (x as u128) * (y as u128) + carry;
            out[i + j] = t as u64;
            carry = t >> 64;
        }
        let mut idx = i + b.len();
        while carry > 0 && idx < out.len() {
            let t = out[idx] as u128 + carry;
            out[idx] = t as u64;
            carry = t >> 64;
            idx += 1;
        }
    }
    while out.len() > 0 && *out.last().unwrap() == 0 {
        out.pop();
    }
    out
}

fn add_limbs_big(a: &[u64], b: &[u64]) -> Vec<u64> {
    let n = a.len().max(b.len());
    let mut out = vec![0u64; n];
    let mut carry = 0u128;
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(0) as u128;
        let y = b.get(i).copied().unwrap_or(0) as u128;
        let t = x + y + carry;
        out[i] = t as u64;
        carry = t >> 64;
    }
    if carry > 0 {
        out.push(carry as u64);
    }
    out
}

/// Divide a little-endian limb integer by 3 (exact).
fn div_by_3(v: &[u64]) -> Vec<u64> {
    let mut out = vec![0u64; v.len()];
    let mut rem: u128 = 0;
    for i in (0..v.len()).rev() {
        let cur = (rem << 64) | v[i] as u128;
        out[i] = (cur / 3) as u64;
        rem = cur % 3;
    }
    while out.len() > 0 && *out.last().unwrap() == 0 {
        out.pop();
    }
    out
}

// ---------------------------------------------------------------------------
// Runtime-derived Frobenius coefficients
// ---------------------------------------------------------------------------

struct FrobeniusCoeffs {
    /// Γ₁[k] = ξ^((p^k − 1)/3) ∈ Fp2 (Fp6-level coefficients).
    gamma: [Fp2; 7],
    /// γ_k = w^(p^k − 1) ∈ Fp6 (Fp12-level coefficients).
    wgamma: [Fp6; 7],
}

fn frob() -> &'static FrobeniusCoeffs {
    static FROB: OnceLock<FrobeniusCoeffs> = OnceLock::new();
    FROB.get_or_init(|| {
        let p: Vec<u64> = Fp::MODULUS.to_vec();
        let p_minus_1: Vec<u64> = {
            // p − 1: subtract 1 once with borrow propagation
            let mut e = p.clone();
            let mut borrow: u128 = 1;
            for limb in e.iter_mut() {
                let t = (*limb as u128).wrapping_sub(borrow);
                *limb = t as u64;
                borrow = (t >> 64) & 1;
            }
            e
        };

        let mut out = FrobeniusCoeffs {
            gamma: [Fp2::zero(); 7],
            wgamma: [Fp6::zero(); 7],
        };

        // Γ₁[k] = ξ^((p^k − 1)/3); p^k − 1 = (p^(k−1) − 1)·p + (p − 1)
        let mut pk_minus1: Vec<u64> = vec![]; // p^0 − 1 = 0
        for k in 1..=6usize {
            pk_minus1 = add_limbs_big(&mul_limbs_big(&pk_minus1, &p), &p_minus_1);
            let e = div_by_3(&pk_minus1);
            out.gamma[k] = Fp2::xi().pow_limbs(&e);
        }

        // γ₁ = w^(p − 1) via one direct exponentiation in Fp12.
        // w is the generator of Fp12 = Fp6[w]/(w² − v) over Fp6:
        // w = (0, 1) with the Fp6 identity as its w-coefficient.
        let w = Fp12 {
            c0: Fp6::zero(),
            c1: Fp6::one(),
        };
        let gamma1 = w.pow_limbs_fp12(&p_minus_1);
        debug_assert!(
            gamma1.c1.is_zero(),
            "w^(p-1) must lie in the Fp6 (c0) part"
        );
        out.wgamma[1] = gamma1.c0;

        // γ_{k+1} = (γ_k)^p · γ₁: from w^(p^(k+1)−1) = (w^(p^k−1))^p · w^(p−1).
        // The first factor is an Fp6 Frobenius of power ONE (Γ₁[1]).
        // (constants are passed explicitly: reading the global here would
        // deadlock the OnceLock we are currently initialising)
        for k in 1..6usize {
            let gamma1 = out.gamma[1];
            let frobbed = frob6_with_gamma(out.wgamma[k], 1, &gamma1);
            out.wgamma[k + 1] = frobbed * out.wgamma[1];
        }

        // Self-verification: γ_k² == Γ₁[k] (as an element of Fp6 ⊂ Fp12).
        for k in 1..=6usize {
            let check = out.wgamma[k] * out.wgamma[k];
            assert_eq!(
                check,
                Fp6 {
                    c0: out.gamma[k],
                    c1: Fp2::zero(),
                    c2: Fp2::zero()
                },
                "frobenius constant self-check failed for k={}",
                k
            );
        }
        out
    })
}

/// Fp6 Frobenius^k with an explicit Γ₁[k] (no global read — safe during
/// initialisation).
fn frob6_with_gamma(x: Fp6, k: usize, gamma: &Fp2) -> Fp6 {
    let conj = |c: Fp2| -> Fp2 {
        if k % 2 == 1 {
            c.conjugate()
        } else {
            c
        }
    };
    Fp6 {
        c0: conj(x.c0),
        c1: conj(x.c1) * *gamma,
        c2: conj(x.c2) * *gamma * *gamma,
    }
}

impl Fp12 {
    /// Frobenius^k: x → x^(p^k).
    pub fn frobenius_k(&self, k: usize) -> Fp12 {
        let f = frob();
        let gamma = f.gamma[k];
        Fp12 {
            c0: frob6_with_gamma(self.c0, k, &gamma),
            c1: frob6_with_gamma(self.c1, k, &gamma) * f.wgamma[k],
        }
    }

    /// Exponentiation by arbitrary little-endian limbs.
    pub fn pow_limbs_fp12(&self, exp: &[u64]) -> Fp12 {
        let mut acc = Fp12::one();
        let mut started = false;
        for limb in exp.iter().rev() {
            for b in (0..64).rev() {
                if started {
                    acc = acc.square();
                }
                if (limb >> b) & 1 == 1 {
                    if started {
                        acc = acc * *self;
                    } else {
                        acc = *self;
                        started = true;
                    }
                }
            }
        }
        acc
    }
}

// ---------------------------------------------------------------------------
// Miller loop (Beuchat et al., eprint 2010/354)
// ---------------------------------------------------------------------------

/// Precomputed G2 data for the Miller loop.
pub struct G2Prepared {
    pub infinity: bool,
    /// 68 line-coefficient triples (a, b, c) ∈ Fp2³
    pub coeffs: Vec<(Fp2, Fp2, Fp2)>,
}

impl From<G2Affine> for G2Prepared {
    fn from(q: G2Affine) -> G2Prepared {
        let mut r = G2Projective {
            x: q.x,
            y: q.y,
            z: Fp2::one(),
        };
        let mut coeffs = Vec::with_capacity(68);
        let mut found_one = false;
        for i in (0..64).rev() {
            let bit = ((BLS_X >> 1) >> i) & 1 == 1;
            if !found_one {
                found_one = bit;
                continue;
            }
            coeffs.push(doubling_step(&mut r));
            if bit {
                coeffs.push(addition_step(&mut r, &q));
            }
        }
        coeffs.push(doubling_step(&mut r));
        assert_eq!(coeffs.len(), 68);
        G2Prepared {
            infinity: q.infinity,
            coeffs,
        }
    }
}

/// Algorithm 26 (point doubling step) — line coefficients.
fn doubling_step(r: &mut G2Projective) -> (Fp2, Fp2, Fp2) {
    let tmp0 = r.x.square();
    let tmp1 = r.y.square();
    let tmp2 = tmp1.square();
    let tmp3 = (tmp1 + r.x).square() - tmp0 - tmp2;
    let tmp3 = tmp3 + tmp3;
    let tmp4 = tmp0 + tmp0 + tmp0;
    let tmp6 = r.x + tmp4;
    let tmp5 = tmp4.square();
    let zsquared = r.z.square();
    r.x = tmp5 - tmp3 - tmp3;
    r.z = (r.z + r.y).square() - tmp1 - zsquared;
    r.y = (tmp3 - r.x) * tmp4;
    let tmp2 = tmp2 + tmp2;
    let tmp2 = tmp2 + tmp2;
    let tmp2 = tmp2 + tmp2;
    r.y = r.y - tmp2;
    let tmp3 = tmp4 * zsquared;
    let tmp3 = tmp3 + tmp3;
    let tmp3 = -tmp3;
    let tmp6 = tmp6.square() - tmp0 - tmp5;
    let tmp1 = tmp1 + tmp1;
    let tmp1 = tmp1 + tmp1;
    let tmp6 = tmp6 - tmp1;
    let tmp0 = r.z * zsquared;
    let tmp0 = tmp0 + tmp0;
    (tmp0, tmp3, tmp6)
}

/// Algorithm 27 (mixed addition step) — line coefficients.
fn addition_step(r: &mut G2Projective, q: &G2Affine) -> (Fp2, Fp2, Fp2) {
    let zsquared = r.z.square();
    let ysquared = q.y.square();
    let t0 = zsquared * q.x;
    let t1 = ((q.y + r.z).square() - ysquared - zsquared) * zsquared;
    let t2 = t0 - r.x;
    let t3 = t2.square();
    let t4 = t3 + t3;
    let t4 = t4 + t4;
    let t5 = t4 * t2;
    let t6 = t1 - r.y - r.y;
    let t9 = t6 * q.x;
    let t7 = t4 * r.x;
    r.x = t6.square() - t5 - t7 - t7;
    r.z = (r.z + t2).square() - zsquared - t3;
    let t10 = q.y + r.z;
    let t8 = (t7 - r.x) * t6;
    let t0 = r.y * t5;
    let t0 = t0 + t0;
    r.y = t8 - t0;
    let t10 = t10.square() - ysquared;
    let ztsquared = r.z.square();
    let t10 = t10 - ztsquared;
    let t9 = t9 + t9 - t10;
    let t10 = r.z + r.z;
    let t6 = -t6;
    let t1 = t6 + t6;
    (t10, t1, t9)
}

/// Apply a line evaluation: sparse Fp12 multiply.
fn ell(f: Fp12, coeffs: &(Fp2, Fp2, Fp2), p: &G1Affine) -> Fp12 {
    let mut c0 = coeffs.0;
    let mut c1 = coeffs.1;
    let py = p.y;
    let px = p.x;
    c0 = Fp2 {
        c0: c0.c0 * py,
        c1: c0.c1 * py,
    };
    c1 = Fp2 {
        c0: c1.c0 * px,
        c1: c1.c1 * px,
    };
    f.mul_by_014(&coeffs.2, &c1, &c0)
}

/// Σᵢ e(pᵢ, qᵢ) — the Miller-loop aggregation primitive (one final
/// exponentiation amortises over any number of terms).
pub fn multi_miller_loop(terms: &[(&G1Affine, &G2Prepared)]) -> Fp12 {
    let mut f = Fp12::one();
    let mut index = 0usize;
    let mut found_one = false;
    for i in (0..64).rev() {
        let bit = ((BLS_X >> 1) >> i) & 1 == 1;
        if !found_one {
            found_one = bit;
            continue;
        }
        for (p, q) in terms {
            if !(p.infinity || q.infinity) {
                f = ell(f, &q.coeffs[index], p);
            }
        }
        index += 1;
        if bit {
            for (p, q) in terms {
                if !(p.infinity || q.infinity) {
                    f = ell(f, &q.coeffs[index], p);
                }
            }
            index += 1;
        }
        f = f.square();
    }
    // final doubling step
    for (p, q) in terms {
        if !(p.infinity || q.infinity) {
            f = ell(f, &q.coeffs[index], p);
        }
    }
    if BLS_X_IS_NEGATIVE {
        f = f.conjugate();
    }
    f
}

// ---------------------------------------------------------------------------
// Final exponentiation (Fuentes–Castañeda–Rodríguez-Henríquez "fountain")
// ---------------------------------------------------------------------------

/// Cyclotomic exponentiation by |x|, conjugated for the negative parameter.
fn cyclotomic_exp_x(f: Fp12) -> Fp12 {
    let mut tmp = Fp12::one();
    let mut found_one = false;
    for i in (0..64).rev() {
        let bit = (BLS_X >> i) & 1 == 1;
        if found_one {
            tmp = tmp.cyclotomic_square();
        } else {
            found_one = bit;
        }
        if bit {
            tmp = tmp * f;
        }
    }
    tmp.conjugate()
}

/// The final exponentiation: f → f^((p¹² − 1)/r).
pub fn final_exponentiation(f: &Fp12) -> Option<Fp12> {
    if f.is_zero() {
        return None;
    }
    // Easy part: f^((p^6 − 1)·(p^2 + 1))
    let t0 = f.frobenius_k(6);
    let t1 = f.invert()?;
    let mut t2 = t0 * t1;
    let t1 = t2;
    t2 = t2.frobenius_k(2);
    t2 = t2 * t1;
    let f = t2;

    // Hard part — fountain chain.
    let t1 = f.cyclotomic_square().conjugate();
    let mut t3 = cyclotomic_exp_x(f);
    let mut t4 = t1 * t3;
    let t1 = cyclotomic_exp_x(t4);
    t4 = t4.conjugate();
    let mut f = f * t4;
    t4 = t3.cyclotomic_square();
    let t0 = cyclotomic_exp_x(t1);
    t3 = t3 * t0;
    t3 = t3.frobenius_k(2);
    f = f * t3;
    t4 = t4 * cyclotomic_exp_x(t0);
    f = f * cyclotomic_exp_x(t4);
    t4 = t4 * t2.conjugate();
    t2 = t2 * t1;
    t2 = t2.frobenius_k(3);
    f = f * t2;
    t4 = t4.frobenius_k(1);
    f = f * t4;

    Some(f)
}

/// The full pairing e(p, q).
pub fn pairing(p: &G1Affine, q: &G2Affine) -> Option<Fp12> {
    if p.infinity || q.infinity {
        return Some(Fp12::one());
    }
    let qp = G2Prepared::from(*q);
    let ml = multi_miller_loop(&[(p, &qp)]);
    final_exponentiation(&ml)
}

/// Check Π e(pᵢ, qᵢ) == 1 with a single final exponentiation.
pub fn pairing_check(terms: &[(&G1Affine, &G2Prepared)]) -> bool {
    let ml = multi_miller_loop(terms);
    match final_exponentiation(&ml) {
        Some(f) => f.is_one(),
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::g1::G1Projective;
    use crate::g2::G2Projective;

    #[test]
    fn frobenius_constants_self_check() {
        let f = frob();
        for k in 1..=6 {
            assert!(!f.gamma[k].is_zero());
            assert!(!f.wgamma[k].is_zero());
        }
    }

    #[test]
    fn frobenius_matches_exponentiation() {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"frob-test-seed-00000000000000000");
        let mut b = [0u8; 32];
        rng.next_bytes(&mut b);
        let x = Fp12 {
            c0: Fp6 {
                c0: Fp2::new(Fp::from_le_bytes32(&b), Fp::from_le_bytes32(&b)),
                c1: Fp2::one(),
                c2: Fp2::zero(),
            },
            c1: Fp6 {
                c0: Fp2::one(),
                c1: Fp2::zero(),
                c2: Fp2::one(),
            },
        };
        let p: Vec<u64> = Fp::MODULUS.to_vec();
        let mut e: Vec<u64> = vec![1];
        for k in 1..=6usize {
            e = mul_limbs_big(&e, &p);
            let by_exp = x.pow_limbs_fp12(&e);
            let by_frob = x.frobenius_k(k);
            assert_eq!(by_exp, by_frob, "frobenius mismatch k={}", k);
        }
    }

    #[test]
    fn pairing_bilinear() {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"pairing-test-seed-00000000000000");
        for _ in 0..3 {
            let a = rng.next_fr(true);
            let b = rng.next_fr(true);
            // e(ab·P, Q) == e(P, ab·Q) == e(aP, bQ)
            let ab = a * b;
            let g1_ab = G1Projective::generator().mul_fr(&ab).to_affine();
            let g2_ab = G2Projective::generator().mul_fr(&ab).to_affine();
            let g1_a = G1Projective::generator().mul_fr(&a).to_affine();
            let g2_b = G2Projective::generator().mul_fr(&b).to_affine();
            let e1 = pairing(&g1_ab, &G2Affine::generator()).unwrap();
            let e2 = pairing(&G1Affine::generator(), &g2_ab).unwrap();
            let e3 = pairing(&g1_a, &g2_b).unwrap();
            assert_eq!(e1, e2, "e(abP,Q) != e(P,abQ)");
            assert_eq!(e1, e3, "e(abP,Q) != e(aP,bQ)");
            // e((a+c)P, Q) = e(aP,Q)·e(cP,Q)
            let c = rng.next_fr(true);
            let g1_apc = G1Projective::generator().mul_fr(&(a + c)).to_affine();
            let g1_c = G1Projective::generator().mul_fr(&c).to_affine();
            let e_sum = pairing(&g1_apc, &G2Affine::generator()).unwrap();
            let e_a = pairing(&g1_a, &G2Affine::generator()).unwrap();
            let e_c = pairing(&g1_c, &G2Affine::generator()).unwrap();
            assert_eq!(e_sum, e_a * e_c, "additivity failed");
        }
    }

    #[test]
    fn pairing_nondegenerate() {
        let e = pairing(&G1Affine::generator(), &G2Affine::generator()).unwrap();
        assert!(!e.is_one());
        // e^r == 1
        let e_r = e.pow_limbs_fp12(&zoda_math::Fr::MODULUS);
        assert!(e_r.is_one());
    }

    #[test]
    fn pairing_check_works() {
        let g1 = G1Affine::generator();
        let g2 = G2Prepared::from(G2Affine::generator());
        // e(g1, g2) ≠ 1
        assert!(!pairing_check(&[(&g1, &g2)]));
        // e(g1, g2) · e(-g1, g2) == 1
        let neg_g1 = g1.neg();
        assert!(pairing_check(&[(&g1, &g2), (&neg_g1, &g2)]));
        // identity terms contribute nothing
        assert!(pairing_check(&[
            (&g1, &g2),
            (&G1Affine::identity(), &g2),
            (&neg_g1, &g2)
        ]));
    }
}
