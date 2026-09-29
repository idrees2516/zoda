//! The Goldilocks field `p = 2^64 - 2^32 + 1`.
//!
//! Native `u64` arithmetic with the special-form reduction makes this the
//! fastest field in the stack — used for standalone ZODA tensor codes and
//! benchmarks where BLS12-381 compatibility is not required.

use crate::PrimeField;
use core::hash::{Hash, Hasher};

#[derive(Copy, Clone, PartialEq, Eq, Default)]
#[repr(transparent)]
pub struct Goldilocks(pub u64);

pub const GOLDILOCKS_MODULUS: u64 = 0xffff_ffff_0000_0001;

impl Goldilocks {
    /// Reduce the 128-bit value `t` modulo p using `2^64 ≡ 2^32 - 1`.
    ///
    /// Iteratively folds the high half into the low half; preserves the
    /// residue class exactly and converges in ≤ 4 steps, then applies a
    /// single conditional subtraction. This is the correctness-first
    /// formulation of the standard Goldilocks reduction.
    #[inline(always)]
    pub const fn reduce128(t: u128) -> u64 {
        let mut v = t;
        while (v >> 64) != 0 {
            let lo = (v as u64) as u128;
            let hi = v >> 64;
            v = lo + hi * 0xffff_ffff;
        }
        let mut r = v as u64;
        if r >= GOLDILOCKS_MODULUS {
            r -= GOLDILOCKS_MODULUS;
        }
        r
    }

    #[inline(always)]
    pub const fn add_impl(a: u64, b: u64) -> u64 {
        let (s, over) = a.overflowing_add(b);
        let mut r = s;
        if over {
            r = r.wrapping_add(0xffff_ffff); // 2^64 ≡ 2^32-1
        }
        if r >= GOLDILOCKS_MODULUS {
            r -= GOLDILOCKS_MODULUS;
        }
        r
    }

    #[inline(always)]
    pub const fn sub_impl(a: u64, b: u64) -> u64 {
        let (d, under) = a.overflowing_sub(b);
        if under {
            d.wrapping_sub(0xffff_ffff) // -2^64 ≡ -(2^32-1)
        } else {
            d
        }
    }
}

impl core::fmt::Debug for Goldilocks {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Goldilocks({:#x})", self.0)
    }
}

impl core::ops::Add for Goldilocks {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        Goldilocks(Self::add_impl(self.0, rhs.0))
    }
}
impl core::ops::Sub for Goldilocks {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        Goldilocks(Self::sub_impl(self.0, rhs.0))
    }
}
impl core::ops::Mul for Goldilocks {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        Goldilocks(Self::reduce128((self.0 as u128) * (rhs.0 as u128)))
    }
}
impl core::ops::AddAssign for Goldilocks {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl core::ops::SubAssign for Goldilocks {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}
impl core::ops::MulAssign for Goldilocks {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl core::ops::Neg for Goldilocks {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        if self.0 == 0 {
            self
        } else {
            Goldilocks(GOLDILOCKS_MODULUS - self.0)
        }
    }
}

impl Hash for Goldilocks {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}

impl PrimeField for Goldilocks {
    const LIMBS: usize = 1;
    const ONE: Self = Goldilocks(1);
    const GENERATOR: Self = Goldilocks(7); // 7 is a Goldilocks generator
    const TWO_ADICITY: u32 = 32;

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
        Some(Goldilocks(pow_mod(self.0, GOLDILOCKS_MODULUS - 2)))
    }
    fn from_u64(v: u64) -> Self {
        Goldilocks(v % GOLDILOCKS_MODULUS)
    }
    fn to_le_bytes(self) -> [u8; 32] {
        let mut out = [0u8; 32];
        out[..8].copy_from_slice(&self.0.to_le_bytes());
        out
    }
    fn from_le_bytes_mod_order(bytes: &[u8; 32]) -> Self {
        let mut w = [0u8; 8];
        w.copy_from_slice(&bytes[..8]);
        Goldilocks(u64::from_le_bytes(w) % GOLDILOCKS_MODULUS)
    }
    fn from_be_bytes_mod_order(bytes: &[u8]) -> Self {
        let mut le = [0u8; 32];
        for (i, b) in bytes.iter().rev().take(32).enumerate() {
            le[i] = *b;
        }
        Self::from_le_bytes_mod_order(&le)
    }
}

/// Modular exponentiation in the Goldilocks field.
pub fn pow_mod(mut base: u64, mut e: u64) -> u64 {
    let mut res: u64 = 1;
    while e != 0 {
        if e & 1 == 1 {
            res = Goldilocks::reduce128((res as u128) * (base as u128));
        }
        base = Goldilocks::reduce128((base as u128) * (base as u128));
        e >>= 1;
    }
    res
}

impl Goldilocks {
    /// Primitive root of unity of order `2^k` (k ≤ 32).
    pub fn root_of_unity(k: u32) -> Goldilocks {
        debug_assert!(k > 0 && k <= 32);
        // (p-1) = 2^32 * odd with odd = 0xffffffff
        let e = 0xffff_ffffu64 << (32 - k);
        Goldilocks(pow_mod(7, e))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn goldilocks_axioms() {
        let mut rng = crate::ZodaRng::from_seed(*b"goldilocks-test-seed-aaaaaaaaaaa");
        for _ in 0..500 {
            let a = Goldilocks(rng.next_u64());
            let b = Goldilocks(rng.next_u64());
            let c = Goldilocks(rng.next_u64());
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert!((a * b).0 < GOLDILOCKS_MODULUS);
            // ground-truth check against u128 modular arithmetic
            let t = (a.0 as u128 * b.0 as u128) % GOLDILOCKS_MODULUS as u128;
            assert_eq!((a * b).0 as u128, t);
            if !a.is_zero() {
                assert_eq!(a * a.invert().unwrap(), Goldilocks(1));
            }
        }
    }

    #[test]
    fn roots() {
        for k in 1..=16u32 {
            let w = Goldilocks::root_of_unity(k);
            assert_eq!(w.pow_u64(1 << k), Goldilocks(1));
            assert_ne!(w.pow_u64(1 << (k - 1)), Goldilocks(1));
        }
    }
}
