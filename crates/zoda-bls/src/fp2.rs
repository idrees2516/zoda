//! `Fp2 = Fp[u] / (u² + 1)` — the quadratic extension of BLS12-381's base
//! field (Karatsuba arithmetic, norm-based inversion and square roots).

use crate::fp::Fp;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct Fp2 {
    pub c0: Fp,
    pub c1: Fp,
}

impl Fp2 {
    pub const fn zero() -> Fp2 {
        Fp2 {
            c0: Fp::zero(),
            c1: Fp::zero(),
        }
    }
    pub const fn one() -> Fp2 {
        Fp2 {
            c0: Fp::one(),
            c1: Fp::zero(),
        }
    }
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// The cubic non-residue ξ = 1 + u used to build Fp6 (standard
    /// BLS12-381 tower: Fp6 = Fp2[v]/(v³ − ξ), Fp12 = Fp6[w]/(w² − v)).
    pub fn xi() -> Fp2 {
        Fp2 {
            c0: Fp::one(),
            c1: Fp::one(),
        }
    }

    pub fn new(c0: Fp, c1: Fp) -> Fp2 {
        Fp2 { c0, c1 }
    }

    /// Conjugate (Frobenius on Fp2).
    #[inline]
    pub fn conjugate(&self) -> Fp2 {
        Fp2 {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    /// Norm a·ā = c0² + c1² ∈ Fp.
    #[inline]
    pub fn norm(&self) -> Fp {
        (self.c0 * self.c0 + self.c1 * self.c1)
    }

    #[inline]
    pub fn square(&self) -> Fp2 {
        // (a + bu)² = (a² − b²) + 2ab·u   (u² = −1)
        let a = self.c0;
        let b = self.c1;
        let ab = a * b;
        Fp2 {
            c0: (a + b) * (a - b),
            c1: ab + ab,
        }
    }

    pub fn invert(&self) -> Option<Fp2> {
        let n = self.norm().invert()?;
        Some(Fp2 {
            c0: self.c0 * n,
            c1: -(self.c1 * n),
        })
    }

    /// Square root in Fp2 via the norm method (p ≡ 3 mod 4).
    pub fn sqrt(&self) -> Option<Fp2> {
        if self.is_zero() {
            return Some(Fp2::zero());
        }
        if self.c1.is_zero() {
            // pure real: sqrt in Fp or i·sqrt(−x)
            return match self.c0.sqrt() {
                Some(r) => Some(Fp2 { c0: r, c1: Fp::zero() }),
                None => {
                    // -c0 must be a square (−1 is a non-residue)
                    let r = (-self.c0).sqrt()?;
                    Some(Fp2 {
                        c0: Fp::zero(),
                        c1: r,
                    })
                }
            };
        }
        // x = α + βu. t = sqrt(norm(x)) ∈ Fp; then u = sqrt((α ± t)/2), v = β/(2u)
        let alpha = self.c0;
        let beta = self.c1;
        let norm = alpha * alpha + beta * beta;
        let t = norm.sqrt()?;
        // exactly one of (α+t)/2, (α−t)/2 is a square in Fp
        let half = Fp::from_u64(2).invert().unwrap();
        let mut u = None;
        for sign in [1i8, -1i8] {
            let cand = (alpha + if sign > 0 { t } else { -t }) * half;
            if let Some(r) = cand.sqrt() {
                u = Some(r);
                break;
            }
        }
        let u = u?;
        if u.is_zero() {
            return None;
        }
        let v = beta * (Fp::from_u64(2) * u).invert().unwrap();
        let candidate = Fp2 { c0: u, c1: v };
        if candidate.square() == *self {
            Some(candidate)
        } else {
            None
        }
    }

    /// Exponentiation by an arbitrary-length LE limb exponent.
    pub fn pow_limbs(&self, exp: &[u64]) -> Fp2 {
        let mut acc = Fp2::one();
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

    /// sgn0 per RFC 9380 (sign of c0, then of c1 if c0 is zero).
    pub fn sgn0(&self) -> bool {
        // sign_0(x) = x mod 2 (parity of the standard representation)
        let c0_parity = {
            let r = self.c0.to_repr();
            r[0] & 1 == 1
        };
        let c1_parity = {
            let r = self.c1.to_repr();
            r[0] & 1 == 1
        };
        let c0_zero = self.c0.is_zero();
        c0_parity || (c0_zero && c1_parity)
    }

    /// Canonical big-endian byte encoding (c1 || c0), 96 bytes.
    pub fn to_be_bytes(&self) -> [u8; 96] {
        let mut out = [0u8; 96];
        let c1 = self.c1.to_repr();
        let c0 = self.c0.to_repr();
        // repr limbs are little-endian; emit most significant limb first
        for i in 0..6 {
            out[i * 8..i * 8 + 8].copy_from_slice(&c1[5 - i].to_be_bytes());
        }
        for i in 0..6 {
            out[48 + i * 8..48 + i * 8 + 8].copy_from_slice(&c0[5 - i].to_be_bytes());
        }
        out
    }

    /// Multiply by (u + 1) (used for Fp6 non-residue power computations).
    #[inline]
    pub fn mul_by_1_plus_u(&self) -> Fp2 {
        // (a + bu)(1 + u) = (a − b) + (a + b)u
        Fp2 {
            c0: self.c0 - self.c1,
            c1: self.c0 + self.c1,
        }
    }
}

impl Add for Fp2 {
    type Output = Fp2;
    #[inline(always)]
    fn add(self, rhs: Fp2) -> Fp2 {
        Fp2 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}
impl Sub for Fp2 {
    type Output = Fp2;
    #[inline(always)]
    fn sub(self, rhs: Fp2) -> Fp2 {
        Fp2 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}
impl Neg for Fp2 {
    type Output = Fp2;
    #[inline(always)]
    fn neg(self) -> Fp2 {
        Fp2 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}
impl Mul for Fp2 {
    type Output = Fp2;
    #[inline(always)]
    fn mul(self, rhs: Fp2) -> Fp2 {
        // Karatsuba: c0 = a0b0 − a1b1, c1 = (a0+a1)(b0+b1) − a0b0 − a1b1
        let v0 = self.c0 * rhs.c0;
        let v1 = self.c1 * rhs.c1;
        let t = (self.c0 + self.c1) * (rhs.c0 + rhs.c1);
        Fp2 {
            c0: v0 - v1,
            c1: t - v0 - v1,
        }
    }
}
impl AddAssign for Fp2 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Fp2) {
        *self = *self + rhs;
    }
}
impl SubAssign for Fp2 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Fp2) {
        *self = *self - rhs;
    }
}
impl MulAssign for Fp2 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Fp2) {
        *self = *self * rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand_fp2(i: u64) -> Fp2 {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"fp2-test-seed-000000000000000000");
        let _ = i;
        let mut b0 = [0u8; 32];
        let mut b1 = [0u8; 32];
        rng.next_bytes(&mut b0);
        rng.next_bytes(&mut b1);
        Fp2::new(Fp::from_le_bytes32(&b0), Fp::from_le_bytes32(&b1))
    }

    #[test]
    fn fp2_axioms() {
        for _ in 0..100 {
            let a = rand_fp2(0);
            let b = rand_fp2(1);
            let c = rand_fp2(2);
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a.square(), a * a);
            if !a.is_zero() {
                assert_eq!(a * a.invert().unwrap(), Fp2::one());
            }
        }
    }

    #[test]
    fn fp2_sqrt() {
        for i in 0..25 {
            let a = rand_fp2(i);
            let sq = a.square();
            let r = sq.sqrt().expect("squares are QR");
            assert!(r == a || r == -a, "sqrt failed for i={}", i);
        }
        // i (u) has sqrt since −1 is a square in Fp2
        let u = Fp2::new(Fp::zero(), Fp::one());
        let r = u.sqrt().unwrap();
        assert_eq!(r.square(), u);
    }

    #[test]
    fn fp2_conj_norm() {
        let a = rand_fp2(0);
        assert_eq!(a * a.conjugate(), Fp2::new(a.norm(), Fp::zero()));
    }
}
