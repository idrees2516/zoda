//! `Fp12 = Fp6[w] / (w² − v)` — the top of the BLS12-381 tower, target group
//! of the pairing. Includes cyclotomic squaring (Granger–Scott), the sparse
//! `mul_by_014` used by the Miller loop and Frobenius maps with
//! runtime-derived constants.

use crate::fp2::Fp2;
use crate::fp6::Fp6;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
}

impl Fp12 {
    pub const fn one() -> Fp12 {
        Fp12 {
            c0: Fp6::one(),
            c1: Fp6::zero(),
        }
    }
    pub fn is_one(&self) -> bool {
        self.c0 == Fp6::one() && self.c1.is_zero()
    }
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    #[inline]
    pub fn square(&self) -> Fp12 {
        // (a + b w)² = (a² + b² v) + 2ab w
        let a = self.c0;
        let b = self.c1;
        let ab = a * b;
        Fp12 {
            c0: (a + b) * (a + b.mul_by_nonresidue()) - ab - ab.mul_by_nonresidue(),
            c1: ab + ab,
        }
    }

    /// Cyclotomic squaring — valid only on elements of the cyclotomic
    /// subgroup `GΦ₆(Fp6)` (pairing outputs); roughly 2× faster than the
    /// generic square. Algorithm 5.5.4 of the Guide to Pairing-Based
    /// Cryptography, as described by Granger–Scott (eprint 2009/565).
    pub fn cyclotomic_square(&self) -> Fp12 {
        fn fp4_square(a: Fp2, b: Fp2) -> (Fp2, Fp2) {
            let t0 = a * a;
            let t1 = b * b;
            let mut t2 = t1.mul_by_1_plus_u(); // t1 * (1+u): ξ-multiplication
            let c0 = t2 + t0;
            t2 = a + b;
            t2 = t2 * t2;
            t2 = t2 - t0;
            let c1 = t2 - t1;
            (c0, c1)
        }
        // Unpack into 6 Fp2 elements.
        let mut z0 = self.c0.c0;
        let mut z4 = self.c0.c1;
        let mut z3 = self.c0.c2;
        let mut z2 = self.c1.c0;
        let mut z1 = self.c1.c1;
        let mut z5 = self.c1.c2;

        let (t0, t1) = fp4_square(z0, z1);
        z0 = t0 - self.c0.c0;
        z0 = z0 + z0 + t0;
        z1 = t1 + self.c1.c1;
        z1 = z1 + z1 + t1;

        let (t0, t1) = fp4_square(z2, z3);
        let (t2, t3) = fp4_square(z4, z5);

        z4 = t0 - self.c0.c1;
        z4 = z4 + z4 + t0;
        z5 = t1 + self.c1.c2;
        z5 = z5 + z5 + t1;

        let t0 = t3.mul_by_1_plus_u();
        z2 = t0 + self.c1.c0;
        z2 = z2 + z2 + t0;
        z3 = t2 - self.c0.c2;
        z3 = z3 + z3 + t2;

        Fp12 {
            c0: Fp6 {
                c0: z0,
                c1: z4,
                c2: z3,
            },
            c1: Fp6 {
                c0: z2,
                c1: z1,
                c2: z5,
            },
        }
    }

    pub fn invert(&self) -> Option<Fp12> {
        // (a + bw)⁻¹ = (a − bw)/(a² − b²v)
        let a = self.c0;
        let b = self.c1;
        let den = a * a - b.mul_by_nonresidue() * b;
        let inv = den.invert()?;
        Some(Fp12 {
            c0: a * inv,
            c1: -(b * inv),
        })
    }

    /// Conjugation = Frobenius⁶ — the cheap inverse on the cyclotomic
    /// subgroup.
    #[inline]
    pub fn conjugate(&self) -> Fp12 {
        Fp12 {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    /// Sparse multiplication by an element with Fp2 coefficients at
    /// positions 0, 1 and 4 of the six-coefficient vector
    /// (c0.c0, c0.c1, c0.c2, c1.c0, c1.c1, c1.c2).
    #[inline]
    pub fn mul_by_014(&self, c0: &Fp2, c1: &Fp2, c4: &Fp2) -> Fp12 {
        let aa = self.c0.mul_by_01(c0, c1);
        let bb = self.c1.mul_by_1(c4);
        let o = *c1 + *c4;
        let mut c1_new = self.c1 + self.c0;
        c1_new = c1_new.mul_by_01(c0, &o);
        c1_new = c1_new - aa - bb;
        let mut c0_new = bb.mul_by_nonresidue();
        c0_new = c0_new + aa;
        Fp12 {
            c0: c0_new,
            c1: c1_new,
        }
    }
}

impl Add for Fp12 {
    type Output = Fp12;
    #[inline(always)]
    fn add(self, rhs: Fp12) -> Fp12 {
        Fp12 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}
impl Sub for Fp12 {
    type Output = Fp12;
    #[inline(always)]
    fn sub(self, rhs: Fp12) -> Fp12 {
        Fp12 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}
impl Neg for Fp12 {
    type Output = Fp12;
    #[inline(always)]
    fn neg(self) -> Fp12 {
        Fp12 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}
impl Mul for Fp12 {
    type Output = Fp12;
    #[inline]
    fn mul(self, rhs: Fp12) -> Fp12 {
        // Karatsuba over the quadratic extension with w² = v:
        let a = self.c0;
        let b = self.c1;
        let c = rhs.c0;
        let d = rhs.c1;
        let ac = a * c;
        let bd = b * d;
        let mid = (a + b) * (c + d);
        Fp12 {
            c0: ac + bd.mul_by_nonresidue(),
            c1: mid - ac - bd,
        }
    }
}
impl AddAssign for Fp12 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Fp12) {
        *self = *self + rhs;
    }
}
impl SubAssign for Fp12 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Fp12) {
        *self = *self - rhs;
    }
}
impl MulAssign for Fp12 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Fp12) {
        *self = *self * rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp::Fp;

    fn rand_fp2(rng: &mut zoda_math::ZodaRng) -> Fp2 {
        let mut b0 = [0u8; 32];
        let mut b1 = [0u8; 32];
        rng.next_bytes(&mut b0);
        rng.next_bytes(&mut b1);
        Fp2::new(Fp::from_le_bytes32(&b0), Fp::from_le_bytes32(&b1))
    }

    fn rand_fp12() -> Fp12 {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"fp12-test-seed-00000000000000000");
        Fp12 {
            c0: Fp6 {
                c0: rand_fp2(&mut rng),
                c1: rand_fp2(&mut rng),
                c2: rand_fp2(&mut rng),
            },
            c1: Fp6 {
                c0: rand_fp2(&mut rng),
                c1: rand_fp2(&mut rng),
                c2: rand_fp2(&mut rng),
            },
        }
    }

    #[test]
    fn fp12_axioms() {
        for _ in 0..25 {
            let a = rand_fp12();
            let b = rand_fp12();
            let c = rand_fp12();
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a.square(), a * a);
            if !a.is_zero() {
                assert_eq!(a * a.invert().unwrap(), Fp12::one());
            }
            assert_eq!(a * a.conjugate(), {
                let norm6 = a.c0 * a.c0 - a.c1.mul_by_nonresidue() * a.c1;
                Fp12 { c0: norm6, c1: Fp6::zero() }
            });
        }
    }

    #[test]
    fn fp12_sparse_014() {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"fp12-sparse-seed-000000000000000");
        let a = rand_fp12();
        let b0 = rand_fp2(&mut rng);
        let b1 = rand_fp2(&mut rng);
        let b4 = rand_fp2(&mut rng);
        let full = a
            * Fp12 {
                c0: Fp6 {
                    c0: b0,
                    c1: b1,
                    c2: Fp2::zero(),
                },
                c1: Fp6 {
                    c0: Fp2::zero(),
                    c1: b4,
                    c2: Fp2::zero(),
                },
            };
        assert_eq!(a.mul_by_014(&b0, &b1, &b4), full);
    }
}
