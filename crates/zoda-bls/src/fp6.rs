//! `Fp6 = Fp2[v] / (v³ − ξ)` with ξ = 9 + u — the middle layer of the
//! BLS12-381 tower (D-type/psuedo-M-type multiplication, norm inversion).

use crate::fp2::Fp2;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct Fp6 {
    pub c0: Fp2,
    pub c1: Fp2,
    pub c2: Fp2,
}

impl Fp6 {
    pub const fn zero() -> Fp6 {
        Fp6 {
            c0: Fp2::zero(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }
    pub const fn one() -> Fp6 {
        Fp6 {
            c0: Fp2::one(),
            c1: Fp2::zero(),
            c2: Fp2::zero(),
        }
    }
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    /// Multiply by the non-residue v (v³ = ξ).
    #[inline]
    pub fn mul_by_nonresidue(&self) -> Fp6 {
        // (c0 + c1 v + c2 v²)·v = c0 v + c1 v² + c2 ξ
        Fp6 {
            c0: self.c2 * Fp2::xi(),
            c1: self.c0,
            c2: self.c1,
        }
    }

    #[inline]
    pub fn square(&self) -> Fp6 {
        // Schoolbook over the cubic basis with v³ = ξ:
        // (a + b v + c v²)² = (a² + 2bc ξ) + (2ab + c²ξ) v + (b² + 2ac) v²
        let a = self.c0;
        let b = self.c1;
        let c = self.c2;
        let a2 = a * a;
        let b2 = b * b;
        let c2 = c * c;
        let bc = b * c;
        let ab = a * b;
        let ac = a * c;
        Fp6 {
            c0: a2 + bc * Fp2::xi() + bc * Fp2::xi(),
            c1: ab + ab + c2 * Fp2::xi(),
            c2: b2 + ac + ac,
        }
    }

    pub fn invert(&self) -> Option<Fp6> {
        // Norm-based inverse for the cubic extension with v³ = ξ.
        let a = self.c0;
        let b = self.c1;
        let c = self.c2;
        let xi = Fp2::xi();
        // t0 = a² − bcξ ; t1 = c²ξ − ab ; t2 = b² − ac
        let t0 = a * a - (b * c) * xi;
        let t1 = (c * c) * xi - a * b;
        let t2 = b * b - a * c;
        // denominator = a·t0 + ξ(b·t2 + c·t1)
        let den = a * t0 + xi * (b * t2 + c * t1);
        let inv = den.invert()?;
        Some(Fp6 {
            c0: t0 * inv,
            c1: t1 * inv,
            c2: t2 * inv,
        })
    }

    /// Sparse multiplication by (c0, 0, 0).
    #[inline]
    pub fn mul_by_0(&self, c0: &Fp2) -> Fp6 {
        Fp6 {
            c0: self.c0 * *c0,
            c1: self.c1 * *c0,
            c2: self.c2 * *c0,
        }
    }

    /// Sparse multiplication by (0, c1, 0).
    #[inline]
    pub fn mul_by_1(&self, c1: &Fp2) -> Fp6 {
        // (b v)·(a + b' v + c' v²) = b a v + b b' v² + b c' ξ
        Fp6 {
            c0: self.c2 * *c1 * Fp2::xi(),
            c1: self.c0 * *c1,
            c2: self.c1 * *c1,
        }
    }

    /// Sparse multiplication by (c0, c1, 0).
    #[inline]
    pub fn mul_by_01(&self, c0: &Fp2, c1: &Fp2) -> Fp6 {
        // (a + b v + c v²)(x + y v) = (ax + cyξ) + (ay + bx) v + (by + cx) v²
        let a = self.c0;
        let b = self.c1;
        let c = self.c2;
        let x = *c0;
        let y = *c1;
        let ax = a * x;
        let by = b * y;
        Fp6 {
            c0: ax + (c * y) * Fp2::xi(),
            c1: (a + b) * (x + y) - ax - by,
            c2: by + c * x,
        }
    }
}

impl Add for Fp6 {
    type Output = Fp6;
    #[inline(always)]
    fn add(self, rhs: Fp6) -> Fp6 {
        Fp6 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
            c2: self.c2 + rhs.c2,
        }
    }
}
impl Sub for Fp6 {
    type Output = Fp6;
    #[inline(always)]
    fn sub(self, rhs: Fp6) -> Fp6 {
        Fp6 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
            c2: self.c2 - rhs.c2,
        }
    }
}
impl Neg for Fp6 {
    type Output = Fp6;
    #[inline(always)]
    fn neg(self) -> Fp6 {
        Fp6 {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }
}
impl Mul for Fp6 {
    type Output = Fp6;
    #[inline]
    fn mul(self, rhs: Fp6) -> Fp6 {
        // Schoolbook with Karatsuba over the three coefficients:
        // a = (a0,a1,a2), b = (b0,b1,b2)
        // c0 = a0b0 + ξ(a1b2 + a2b1)
        // c1 = a0b1 + a1b0 + ξ a2b2
        // c2 = a0b2 + a1b1 + a2b0
        let (a0, a1, a2) = (self.c0, self.c1, self.c2);
        let (b0, b1, b2) = (rhs.c0, rhs.c1, rhs.c2);
        let a0b0 = a0 * b0;
        let a1b1 = a1 * b1;
        let a2b2 = a2 * b2;
        // Karatsuka-style symmetric combinations
        let s01 = (a0 + a1) * (b0 + b1) - a0b0 - a1b1;
        let s02 = (a0 + a2) * (b0 + b2) - a0b0 - a2b2;
        let s12 = (a1 + a2) * (b1 + b2) - a1b1 - a2b2;
        let xi = Fp2::xi();
        Fp6 {
            c0: a0b0 + s12 * xi,
            c1: s01 + a2b2 * xi,
            c2: s02 + a1b1,
        }
    }
}
impl AddAssign for Fp6 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Fp6) {
        *self = *self + rhs;
    }
}
impl SubAssign for Fp6 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Fp6) {
        *self = *self - rhs;
    }
}
impl MulAssign for Fp6 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Fp6) {
        *self = *self * rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fp::Fp;

    fn rand_fp6() -> Fp6 {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"fp6-test-seed-000000000000000000");
        let mut r = |i: usize| {
            let mut b = [0u8; 32];
            rng.next_bytes(&mut b);
            let _ = i;
            Fp::from_le_bytes32(&b)
        };
        Fp6 {
            c0: Fp2::new(r(0), r(1)),
            c1: Fp2::new(r(2), r(3)),
            c2: Fp2::new(r(4), r(5)),
        }
    }

    #[test]
    fn fp6_axioms() {
        for _ in 0..50 {
            let a = rand_fp6();
            let b = rand_fp6();
            let c = rand_fp6();
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a.square(), a * a);
            if !a.is_zero() {
                assert_eq!(a * a.invert().unwrap(), Fp6::one());
            }
        }
    }

    #[test]
    fn fp6_nonresidue() {
        // v³ = ξ
        let v = Fp6 {
            c0: Fp2::zero(),
            c1: Fp2::one(),
            c2: Fp2::zero(),
        };
        // v·nr = v², v²·nr = v³ = ξ
        let v3 = v.mul_by_nonresidue().mul_by_nonresidue();
        assert_eq!(v3.c0, Fp2::xi());
        assert!(v3.c1.is_zero() && v3.c2.is_zero());
    }

    #[test]
    fn fp6_sparse_muls() {
        let a = rand_fp6();
        let c0 = a.c1;
        let c1 = a.c2;
        let full = a
            * Fp6 {
                c0,
                c1,
                c2: Fp2::zero(),
            };
        assert_eq!(a.mul_by_01(&c0, &c1), full);
        let only0 = a * Fp6 { c0, c1: Fp2::zero(), c2: Fp2::zero() };
        assert_eq!(a.mul_by_0(&c0), only0);
        let only1 = a * Fp6 { c0: Fp2::zero(), c1, c2: Fp2::zero() };
        assert_eq!(a.mul_by_1(&c1), only1);
    }
}
