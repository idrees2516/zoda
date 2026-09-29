//! The 23-bit NTT-friendly prime field used by the lattice constructions:
//! `q = 8380417 = 2^23 - 2^13 + 1` (the Dilithium modulus).
//!
//! Plain `u32` values reduced with u64 intermediates; the small modulus
//! keeps every operation in fast native arithmetic and admits negacyclic
//! NTTs up to length 8192 (`2^13 | q - 1`).

use crate::PrimeField;
use core::hash::{Hash, Hasher};

pub const Q: u32 = 8_380_417;
/// Two-adicity of `q - 1`.
pub const Q_TWO_ADICITY: u32 = 13;
/// Multiplicative generator of Z_q^*.
pub const Q_GENERATOR: u32 = 10;

#[derive(Copy, Clone, PartialEq, Eq, Default)]
#[repr(transparent)]
pub struct Fq(pub u32);

impl core::fmt::Debug for Fq {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Fq({})", self.0)
    }
}

impl Fq {
    #[inline(always)]
    pub fn add_impl(a: u32, b: u32) -> u32 {
        let s = a + b;
        if s >= Q {
            s - Q
        } else {
            s
        }
    }
    #[inline(always)]
    pub fn sub_impl(a: u32, b: u32) -> u32 {
        a.wrapping_sub(b).wrapping_add(Q) % Q
    }
    #[inline(always)]
    pub fn mul_impl(a: u32, b: u32) -> u32 {
        ((a as u64 * b as u64) % Q as u64) as u32
    }

    /// Primitive `2^k`-th root of unity (k ≤ 13).
    pub fn root_of_unity(k: u32) -> Fq {
        debug_assert!(k > 0 && k <= Q_TWO_ADICITY);
        let e = ((Q - 1) >> k) as u64;
        let mut res: u64 = 1;
        let mut base: u64 = Q_GENERATOR as u64;
        let mut ee = e;
        while ee != 0 {
            if ee & 1 == 1 {
                res = (res * base) % Q as u64;
            }
            base = (base * base) % Q as u64;
            ee >>= 1;
        }
        Fq(res as u32)
    }
}

impl core::ops::Add for Fq {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        Fq(Self::add_impl(self.0, rhs.0))
    }
}
impl core::ops::Sub for Fq {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        Fq(Self::sub_impl(self.0, rhs.0))
    }
}
impl core::ops::Mul for Fq {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        Fq(Self::mul_impl(self.0, rhs.0))
    }
}
impl core::ops::AddAssign for Fq {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl core::ops::SubAssign for Fq {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}
impl core::ops::MulAssign for Fq {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl core::ops::Neg for Fq {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        if self.0 == 0 {
            self
        } else {
            Fq(Q - self.0)
        }
    }
}

impl Hash for Fq {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl PrimeField for Fq {
    const LIMBS: usize = 1;
    const ONE: Self = Fq(1);
    const GENERATOR: Self = Fq(Q_GENERATOR);
    const TWO_ADICITY: u32 = Q_TWO_ADICITY;

    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        core::ops::Add::add(self, rhs)
    }
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        core::ops::Sub::sub(self, rhs)
    }
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        core::ops::Mul::mul(self, rhs)
    }
    #[inline(always)]
    fn neg(self) -> Self {
        core::ops::Neg::neg(self)
    }
    fn invert(self) -> Option<Self> {
        if self.0 == 0 {
            return None;
        }
        let mut res: u64 = 1;
        let mut base = self.0 as u64;
        let mut e = (Q - 2) as u64;
        while e != 0 {
            if e & 1 == 1 {
                res = (res * base) % Q as u64;
            }
            base = (base * base) % Q as u64;
            e >>= 1;
        }
        Some(Fq(res as u32))
    }
    fn from_u64(v: u64) -> Self {
        Fq((v % Q as u64) as u32)
    }
    fn to_le_bytes(self) -> [u8; 32] {
        let mut out = [0u8; 32];
        out[..4].copy_from_slice(&self.0.to_le_bytes());
        out
    }
    fn from_le_bytes_mod_order(bytes: &[u8; 32]) -> Self {
        let mut w = [0u8; 4];
        w.copy_from_slice(&bytes[..4]);
        Fq(u32::from_le_bytes(w) % Q)
    }
    fn from_be_bytes_mod_order(bytes: &[u8]) -> Self {
        let mut v: u64 = 0;
        for b in bytes.iter().take(8) {
            v = (v << 8) | *b as u64;
        }
        Fq((v % Q as u64) as u32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fq_axioms() {
        let mut rng = crate::ZodaRng::from_seed(*b"fq-test-seed-bb00000000000000000");
        for _ in 0..1000 {
            let a = Fq(rng.next_u32() % Q);
            let b = Fq(rng.next_u32() % Q);
            let c = Fq(rng.next_u32() % Q);
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a * a.invert().unwrap(), Fq(1));
        }
    }

    #[test]
    fn roots() {
        for k in 1..=13u32 {
            let w = Fq::root_of_unity(k);
            // order exactly 2^k
            let mut x = w;
            let mut e = 1u64;
            while e < (1u64 << k) {
                x = x * x;
                e *= 2;
            }
            assert_eq!(x, Fq(1), "w^(2^{}) != 1", k);
            if k > 1 {
                assert_ne!(w.pow_u64(1 << (k - 1)), Fq(1), "not primitive for k={}", k);
            }
        }
    }
}
