//! Generic short-Weierstrass curve machinery (projective + affine, complete
//! addition formulas, ZCash-format compression) shared by G1 and G2.
//!
//! Addition/doubling use the complete Renes–Costello–Batina formulas for
//! curves with a = 0 (eprint 2015/1060, Algorithms 7–9) — no edge cases.

use crate::fp::Fp;
use crate::fp2::Fp2;

/// The minimal field interface needed by the generic curve code.
pub trait ExtField: Copy + PartialEq + 'static {
    fn f_zero() -> Self;
    fn f_one() -> Self;
    fn f_add(self, o: Self) -> Self;
    fn f_sub(self, o: Self) -> Self;
    fn f_mul(self, o: Self) -> Self;
    fn f_square(self) -> Self;
    fn f_neg(self) -> Self;
    fn f_is_zero(&self) -> bool;
    fn f_lex_largest(&self) -> bool;
    fn f_sqrt(x: &Self) -> Option<Self>;
    fn f_invert(&self) -> Option<Self>;
}

impl ExtField for Fp {
    fn f_zero() -> Self { Fp::zero() }
    fn f_one() -> Self { Fp::one() }
    fn f_add(self, o: Self) -> Self { self + o }
    fn f_sub(self, o: Self) -> Self { self - o }
    fn f_mul(self, o: Self) -> Self { self * o }
    fn f_square(self) -> Self { self.square() }
    fn f_neg(self) -> Self { -self }
    fn f_is_zero(&self) -> bool { self.is_zero() }
    fn f_lex_largest(&self) -> bool { self.lexicographically_largest() }
    fn f_sqrt(x: &Self) -> Option<Self> { x.sqrt() }
    fn f_invert(&self) -> Option<Self> { self.invert() }
}

impl ExtField for Fp2 {
    fn f_zero() -> Self { Fp2::zero() }
    fn f_one() -> Self { Fp2::one() }
    fn f_add(self, o: Self) -> Self { self + o }
    fn f_sub(self, o: Self) -> Self { self - o }
    fn f_mul(self, o: Self) -> Self { self * o }
    fn f_square(self) -> Self { self.square() }
    fn f_neg(self) -> Self { -self }
    fn f_is_zero(&self) -> bool { self.is_zero() }
    fn f_lex_largest(&self) -> bool {
        self.c1.lexicographically_largest() || (self.c1.is_zero() && self.c0.lexicographically_largest())
    }
    fn f_sqrt(x: &Self) -> Option<Self> { x.sqrt() }
    fn f_invert(&self) -> Option<Self> { self.invert() }
}

/// Parameters of a short-Weierstrass curve y² = x³ + b over an extension
/// field, plus its canonical byte encoding rules.
pub trait CurveConfig: 'static + Copy + PartialEq {
    type Base: ExtField;
    /// The curve coefficient b.
    fn b() -> Self::Base;
    /// 3·b (used by the complete formulas).
    fn b3() -> Self::Base;
    /// Fixed affine generator (x, y).
    fn generator_xy() -> (Self::Base, Self::Base);
    /// Number of bytes per base-field element in compressed form.
    fn elem_bytes() -> usize;

    /// Canonical big-endian bytes of a base element.
    fn base_to_be_bytes(x: &Self::Base, out: &mut [u8]);
    /// Parse a canonical big-endian element; `None` if ≥ p or malformed.
    fn base_from_be_bytes(bytes: &[u8]) -> Option<Self::Base>;
}

/// Affine point (or the identity, flagged).
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct Affine<C: CurveConfig> {
    pub x: C::Base,
    pub y: C::Base,
    pub infinity: bool,
}

/// Projective point (x : y : z); identity is z = 0.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct Projective<C: CurveConfig> {
    pub x: C::Base,
    pub y: C::Base,
    pub z: C::Base,
}

impl<C: CurveConfig> Projective<C> {
    pub fn identity() -> Self {
        Projective {
            x: C::Base::f_zero(),
            y: C::Base::f_one(),
            z: C::Base::f_zero(),
        }
    }

    pub fn is_identity(&self) -> bool {
        self.z.f_is_zero()
    }

    pub fn generator() -> Self {
        let (x, y) = C::generator_xy();
        Projective { x, y, z: C::Base::f_one() }
    }

    /// RCB Algorithm 9 (complete doubling, a = 0).
    pub fn double(&self) -> Self {
        let t0 = self.y.f_square();
        let mut z3 = t0.f_add(t0);
        z3 = z3.f_add(z3);
        z3 = z3.f_add(z3);
        let t1 = self.y.f_mul(self.z);
        let t2 = self.z.f_square();
        let t2 = C::b3().f_mul(t2);
        let x3 = t2.f_mul(z3);
        let y3 = t0.f_add(t2);
        z3 = t1.f_mul(z3);
        let t1 = t2.f_add(t2);
        let t2 = t1.f_add(t2);
        let t0 = t0.f_sub(t2);
        let y3 = t0.f_mul(y3);
        let y3 = x3.f_add(y3);
        let t1 = self.x.f_mul(self.y);
        let x3 = t0.f_mul(t1);
        let x3 = x3.f_add(x3);
        let out = Projective { x: x3, y: y3, z: z3 };
        if out.is_identity() && !self.is_identity() {
            // complete formulas handle identity; select for safety
            return Self::identity();
        }
        out
    }

    /// RCB Algorithm 7 (complete addition, a = 0).
    pub fn add(&self, rhs: &Self) -> Self {
        let t0 = self.x.f_mul(rhs.x);
        let t1 = self.y.f_mul(rhs.y);
        let t2 = self.z.f_mul(rhs.z);
        let t3 = self.x.f_add(self.y);
        let t4 = rhs.x.f_add(rhs.y);
        let t3 = t3.f_mul(t4);
        let t4 = t0.f_add(t1);
        let t3 = t3.f_sub(t4);
        let t4 = self.y.f_add(self.z);
        let x3 = rhs.y.f_add(rhs.z);
        let t4 = t4.f_mul(x3);
        let x3 = t1.f_add(t2);
        let t4 = t4.f_sub(x3);
        let x3 = self.x.f_add(self.z);
        let y3 = rhs.x.f_add(rhs.z);
        let x3 = x3.f_mul(y3);
        let y3 = t0.f_add(t2);
        let y3 = x3.f_sub(y3);
        let x3 = t0.f_add(t0);
        let t0 = x3.f_add(t0);
        let t2 = C::b3().f_mul(t2);
        let z3 = t1.f_add(t2);
        let t1 = t1.f_sub(t2);
        let y3 = C::b3().f_mul(y3);
        let x3 = t4.f_mul(y3);
        let t2 = t3.f_mul(t1);
        let x3 = t2.f_sub(x3);
        let y3 = y3.f_mul(t0);
        let t1 = t1.f_mul(z3);
        let y3 = t1.f_add(y3);
        let t0 = t0.f_mul(t3);
        let z3 = z3.f_mul(t4);
        let z3 = z3.f_add(t0);
        Projective { x: x3, y: y3, z: z3 }
    }

    /// RCB Algorithm 8 (complete mixed addition, a = 0) — affine rhs.
    pub fn add_mixed(&self, rhs: &Affine<C>) -> Self {
        if rhs.infinity {
            return *self;
        }
        let t0 = self.x.f_mul(rhs.x);
        let t1 = self.y.f_mul(rhs.y);
        let t3 = rhs.x.f_add(rhs.y);
        let t4 = self.x.f_add(self.y);
        let t3 = t3.f_mul(t4);
        let t4 = t0.f_add(t1);
        let t3 = t3.f_sub(t4);
        let t4 = rhs.y.f_mul(self.z);
        let t4 = t4.f_add(self.y);
        let y3 = rhs.x.f_mul(self.z);
        let y3 = y3.f_add(self.x);
        let x3 = t0.f_add(t0);
        let t0 = x3.f_add(t0);
        let t2 = C::b3().f_mul(self.z);
        let z3 = t1.f_add(t2);
        let t1 = t1.f_sub(t2);
        let y3 = C::b3().f_mul(y3);
        let x3 = t4.f_mul(y3);
        let t2 = t3.f_mul(t1);
        let x3 = t2.f_sub(x3);
        let y3 = y3.f_mul(t0);
        let t1 = t1.f_mul(z3);
        let y3 = t1.f_add(y3);
        let t0 = t0.f_mul(t3);
        let z3 = z3.f_mul(t4);
        let z3 = z3.f_add(t0);
        Projective { x: x3, y: y3, z: z3 }
    }

    pub fn neg(&self) -> Self {
        Projective {
            x: self.x,
            y: self.y.f_neg(),
            z: self.z,
        }
    }

    /// Convert to affine coordinates.
    ///
    /// The complete formulas use **homogeneous** projective coordinates:
    /// the affine point is (X/Z, Y/Z).
    pub fn to_affine(&self) -> Affine<C> {
        if self.is_identity() {
            return Affine {
                x: C::Base::f_zero(),
                y: C::Base::f_one(),
                infinity: true,
            };
        }
        let zinv = self.z.f_invert().expect("nonzero z has an inverse");
        Affine {
            x: self.x.f_mul(zinv),
            y: self.y.f_mul(zinv),
            infinity: false,
        }
    }

    /// Windowed (w = 4) scalar multiplication by a little-endian limb
    /// scalar.
    pub fn mul_limbs(&self, scalar: &[u64]) -> Self {
        // Build 4-bit window table of [0P, 1P, ..., 15P]
        let mut table = Vec::with_capacity(16);
        table.push(Self::identity());
        table.push(*self);
        for i in 2..16 {
            table.push(table[i - 1].add(self));
        }
        let mut acc = Self::identity();
        // Iterate over 4-bit nibbles, most significant first.
        let total_bits = scalar.len() * 64;
        let mut bit = total_bits;
        while bit > 0 {
            let start = bit.saturating_sub(4);
            // extract nibble covering bits [start, bit)
            let mut nib = 0u32;
            for b in start..bit {
                let limb = b / 64;
                let off = b % 64;
                if limb < scalar.len() && (scalar[limb] >> off) & 1 == 1 {
                    nib |= 1 << (b - start);
                }
            }
            // four doublings (skip leading zeros implicitly by add of zero)
            if !acc.is_identity() {
                acc = acc.double().double().double().double();
            }
            if nib != 0 {
                acc = acc.add(&table[nib as usize]);
            }
            bit = start;
        }
        acc
    }
}

impl<C: CurveConfig> core::ops::Add for Projective<C> {
    type Output = Projective<C>;
    #[inline]
    fn add(self, rhs: Projective<C>) -> Projective<C> {
        Projective::add(&self, &rhs)
    }
}
impl<C: CurveConfig> core::ops::Add<&Projective<C>> for &Projective<C> {
    type Output = Projective<C>;
    #[inline]
    fn add(self, rhs: &Projective<C>) -> Projective<C> {
        Projective::add(self, rhs)
    }
}
impl<C: CurveConfig> core::ops::Sub for Projective<C> {
    type Output = Projective<C>;
    #[inline]
    fn sub(self, rhs: Projective<C>) -> Projective<C> {
        self.add(&rhs.neg())
    }
}

impl<C: CurveConfig> core::ops::Neg for Affine<C> {
    type Output = Affine<C>;
    #[inline]
    fn neg(self) -> Affine<C> {
        Affine::neg(&self)
    }
}

impl<C: CurveConfig> Affine<C> {
    pub fn identity() -> Self {
        Affine {
            x: C::Base::f_zero(),
            y: C::Base::f_one(),
            infinity: true,
        }
    }
    pub fn generator() -> Self {
        let (x, y) = C::generator_xy();
        Affine { x, y, infinity: false }
    }
    pub fn is_identity(&self) -> bool {
        self.infinity
    }
    pub fn neg(&self) -> Self {
        if self.infinity {
            *self
        } else {
            Affine {
                x: self.x,
                y: self.y.f_neg(),
                infinity: false,
            }
        }
    }
    /// Check y² = x³ + b.
    pub fn is_on_curve(&self) -> bool {
        if self.infinity {
            return true;
        }
        let y2 = self.y.f_square();
        let x3b = self.x.f_square().f_mul(self.x).f_add(C::b());
        y2 == x3b
    }
    pub fn to_projective(&self) -> Projective<C> {
        if self.infinity {
            return Projective::identity();
        }
        Projective {
            x: self.x,
            y: self.y,
            z: C::Base::f_one(),
        }
    }

    /// Compressed encoding (ZCash format):
    /// bit 7 = compressed, bit 6 = infinity, bit 5 = y is lex-largest.
    pub fn to_compressed(&self) -> Vec<u8> {
        let n = C::elem_bytes();
        let mut out = vec![0u8; n];
        let x = if self.infinity {
            C::Base::f_zero()
        } else {
            self.x
        };
        C::base_to_be_bytes(&x, &mut out[..n]);
        out[0] |= 0x80;
        if self.infinity {
            out[0] |= 0x40;
        } else if self.y.f_lex_largest() {
            out[0] |= 0x20;
        }
        out
    }

    /// Decompress a compressed point encoding.
    pub fn from_compressed(bytes: &[u8]) -> Option<Self> {
        let n = C::elem_bytes();
        if bytes.len() != n {
            return None;
        }
        if bytes[0] & 0x80 == 0 {
            return None; // not compressed
        }
        let infinity = bytes[0] & 0x40 != 0;
        let largest = bytes[0] & 0x20 != 0;
        let mut xbytes = vec![0u8; n];
        xbytes.copy_from_slice(&bytes[..n]);
        xbytes[0] &= 0x1f;
        if infinity {
            // all remaining bits must be zero
            if xbytes.iter().any(|&b| b != 0) {
                return None;
            }
            return Some(Self::identity());
        }
        let x = C::base_from_be_bytes(&xbytes)?;
        // y² = x³ + b
        let y2 = x.f_square().f_mul(x).f_add(C::b());
        let y = <C::Base as ExtField>::f_sqrt(&y2)?;
        let y = if y.f_lex_largest() == largest { y } else { y.f_neg() };
        Some(Affine { x, y, infinity: false })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fp_extfield_impls() {
        let a = Fp::from_u64(5);
        let b = Fp::from_u64(7);
        assert_eq!(a.f_mul(b), Fp::from_u64(35));
        let c = Fp2::new(Fp::from_u64(1), Fp::from_u64(2));
        assert_eq!(c.f_square(), c.f_mul(c));
    }
}
