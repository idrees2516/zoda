//! Macro-generated Montgomery-form prime fields.
//!
//! The `mont_field!` macro emits a fully reduced, branch-free prime-field
//! type with CIOS (coarsely integrated operand scanning) multiplication and
//! Fermat inversion computed directly in the Montgomery domain.
//!
//! Both the 4-limb BLS12-381 scalar field `Fr` and the 6-limb base field
//! `Fp` (in `zoda-bls`) are generated from this macro, guaranteeing a single
//! audited arithmetic core.

/// Emit the core of a Montgomery-form prime field type (no trait impls).
///
/// * `$name`  – type name
/// * `$limbs` – limb count (4 or 6)
/// * `$mod`   – modulus limbs, least significant first
/// * `$inv`   – `-m^-1 mod 2^64`
/// * `$r1`    – `R mod m` (the Montgomery representation of 1)
/// * `$r2`    – `R^2 mod m` (used to enter Montgomery form)
/// * `$gen`   – generator limbs (standard integer form)
/// * `$adic`  – two-adicity of `m - 1`
#[macro_export]
macro_rules! mont_field {
    ($name:ident, $limbs:expr, $mod:expr, $inv:expr, $r1:expr, $r2:expr, $gen:expr, $adic:expr) => {
        /// A prime-field element in Montgomery form.
        #[derive(Copy, Clone, PartialEq, Eq, Hash, Default)]
        #[repr(transparent)]
        pub struct $name(pub [u64; $limbs]);

        impl $name {
            pub const MODULUS: [u64; $limbs] = $mod;
            pub const INV: u64 = $inv;
            pub const R1: [u64; $limbs] = $r1;
            pub const R2: [u64; $limbs] = $r2;
            pub const TWO_ADICITY: u32 = $adic;
            pub const GENERATOR_RAW: [u64; $limbs] = $gen;

            #[inline(always)]
            const fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
                let t = a as u128 + b as u128 + carry as u128;
                (t as u64, (t >> 64) as u64)
            }
            #[inline(always)]
            const fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
                let t = (a as u128).wrapping_sub(b as u128 + borrow as u128);
                (t as u64, ((t >> 64) as u64) & 1)
            }
            const fn geq(a: &[u64; $limbs], b: &[u64; $limbs]) -> bool {
                let mut i = $limbs;
                while i > 0 {
                    i -= 1;
                    if a[i] > b[i] {
                        return true;
                    }
                    if a[i] < b[i] {
                        return false;
                    }
                }
                true
            }

            /// Montgomery multiplication (CIOS). Result fully reduced.
            #[inline(always)]
            pub const fn mont_mul(a: [u64; $limbs], b: [u64; $limbs]) -> [u64; $limbs] {
                let mut t = [0u64; $limbs + 2];
                let mut i = 0;
                while i < $limbs {
                    let mut carry = 0u128;
                    let mut j = 0;
                    while j < $limbs {
                        let s = t[j] as u128 + a[j] as u128 * b[i] as u128 + carry;
                        t[j] = s as u64;
                        carry = s >> 64;
                        j += 1;
                    }
                    let s = t[$limbs] as u128 + carry;
                    t[$limbs] = s as u64;
                    t[$limbs + 1] = (s >> 64) as u64;

                    let m = t[0].wrapping_mul(Self::INV);
                    let s = t[0] as u128 + m as u128 * Self::MODULUS[0] as u128;
                    let mut carry = s >> 64;
                    let mut j = 1;
                    while j < $limbs {
                        let s = t[j] as u128 + m as u128 * Self::MODULUS[j] as u128 + carry;
                        t[j - 1] = s as u64;
                        carry = s >> 64;
                        j += 1;
                    }
                    let s = t[$limbs] as u128 + carry;
                    t[$limbs - 1] = s as u64;
                    carry = s >> 64;
                    t[$limbs] = t[$limbs + 1] + carry as u64;
                    t[$limbs + 1] = 0;
                    i += 1;
                }

                let mut res = [0u64; $limbs];
                let mut i = 0;
                while i < $limbs {
                    res[i] = t[i];
                    i += 1;
                }
                if t[$limbs] != 0 || Self::geq(&res, &Self::MODULUS) {
                    let mut borrow = 0u64;
                    let mut i = 0;
                    while i < $limbs {
                        let (d, b) = Self::sbb(res[i], Self::MODULUS[i], borrow);
                        res[i] = d;
                        borrow = b;
                        i += 1;
                    }
                }
                res
            }

            const fn one_raw() -> [u64; $limbs] {
                let mut a = [0u64; $limbs];
                a[0] = 1;
                a
            }

            /// Standard integer limbs → Montgomery form.
            #[inline(always)]
            pub const fn to_mont(a: [u64; $limbs]) -> [u64; $limbs] {
                Self::mont_mul(a, Self::R2)
            }
            /// Montgomery form → standard integer limbs.
            #[inline(always)]
            pub const fn from_mont(a: [u64; $limbs]) -> [u64; $limbs] {
                Self::mont_mul(a, Self::one_raw())
            }
            /// Standard integer limbs of this value.
            #[inline]
            pub const fn to_repr(self) -> [u64; $limbs] {
                Self::from_mont(self.0)
            }
            /// Build from standard integer limbs.
            #[inline]
            pub const fn from_repr_limbs(a: [u64; $limbs]) -> Self {
                Self(Self::to_mont(a))
            }
            /// Modular addition on Montgomery forms.
            #[inline(always)]
            pub const fn add_limbs(a: [u64; $limbs], b: [u64; $limbs]) -> [u64; $limbs] {
                let mut out = [0u64; $limbs];
                let mut carry = 0u64;
                let mut i = 0;
                while i < $limbs {
                    let (s, c) = Self::adc(a[i], b[i], carry);
                    out[i] = s;
                    carry = c;
                    i += 1;
                }
                if carry != 0 || Self::geq(&out, &Self::MODULUS) {
                    let mut borrow = 0u64;
                    let mut i = 0;
                    while i < $limbs {
                        let (d, b) = Self::sbb(out[i], Self::MODULUS[i], borrow);
                        out[i] = d;
                        borrow = b;
                        i += 1;
                    }
                }
                out
            }
            /// Modular subtraction on Montgomery forms.
            #[inline(always)]
            pub const fn sub_limbs(a: [u64; $limbs], b: [u64; $limbs]) -> [u64; $limbs] {
                let mut out = [0u64; $limbs];
                let mut borrow = 0u64;
                let mut i = 0;
                while i < $limbs {
                    let (d, b) = Self::sbb(a[i], b[i], borrow);
                    out[i] = d;
                    borrow = b;
                    i += 1;
                }
                if borrow != 0 {
                    let mut carry = 0u64;
                    let mut i = 0;
                    while i < $limbs {
                        let (s, c) = Self::adc(out[i], Self::MODULUS[i], carry);
                        out[i] = s;
                        carry = c;
                        i += 1;
                    }
                }
                out
            }
            /// Negation on Montgomery forms.
            #[inline(always)]
            pub const fn neg_limbs(a: [u64; $limbs]) -> [u64; $limbs] {
                if Self::is_zero_limbs(a) {
                    return a;
                }
                Self::sub_limbs(Self::MODULUS, a)
            }
            const fn is_zero_limbs(a: [u64; $limbs]) -> bool {
                let mut i = 0;
                while i < $limbs {
                    if a[i] != 0 {
                        return false;
                    }
                    i += 1;
                }
                true
            }

            /// Fermat inversion (`a^(p-2)`) computed in the Montgomery
            /// domain. Returns `None` for zero.
            pub fn invert_limbs(a: [u64; $limbs]) -> Option<[u64; $limbs]> {
                if Self::is_zero_limbs(a) {
                    return None;
                }
                // exponent = MODULUS - 2
                let mut e = Self::MODULUS;
                let (v, _b) = Self::sbb(e[0], 2, 0);
                e[0] = v;
                let mut acc = Self::R1;
                let mut base = a;
                // find top set bit
                let mut top = $limbs - 1;
                while top > 0 && e[top] == 0 {
                    top -= 1;
                }
                let bits = 64 - e[top].leading_zeros();
                let mut idx = top;
                let mut started = false;
                while idx < $limbs {
                    let nbits = if idx == top { bits } else { 64 };
                    let mut b = nbits;
                    while b > 0 {
                        b -= 1;
                        if started {
                            acc = Self::mont_mul(acc, acc);
                        }
                        if (e[idx] >> b) & 1 == 1 {
                            if started {
                                acc = Self::mont_mul(acc, base);
                            } else {
                                acc = base;
                                started = true;
                            }
                        }
                    }
                    if idx == 0 {
                        break;
                    }
                    idx -= 1;
                }
                Some(acc)
            }
        }

        impl core::fmt::Debug for $name {
            fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                let v = self.to_repr();
                write!(f, "{}(0x{:016x}", stringify!($name), v[$limbs - 1])?;
                let mut i = $limbs - 1;
                while i > 0 {
                    i -= 1;
                    write!(f, "_{:016x}", v[i])?;
                }
                write!(f, ")")
            }
        }

        impl core::ops::Add for $name {
            type Output = Self;
            #[inline(always)]
            fn add(self, rhs: Self) -> Self {
                Self(Self::add_limbs(self.0, rhs.0))
            }
        }
        impl core::ops::AddAssign for $name {
            #[inline(always)]
            fn add_assign(&mut self, rhs: Self) {
                *self = *self + rhs;
            }
        }
        impl core::ops::Sub for $name {
            type Output = Self;
            #[inline(always)]
            fn sub(self, rhs: Self) -> Self {
                Self(Self::sub_limbs(self.0, rhs.0))
            }
        }
        impl core::ops::SubAssign for $name {
            #[inline(always)]
            fn sub_assign(&mut self, rhs: Self) {
                *self = *self - rhs;
            }
        }
        impl core::ops::Mul for $name {
            type Output = Self;
            #[inline(always)]
            fn mul(self, rhs: Self) -> Self {
                Self(Self::mont_mul(self.0, rhs.0))
            }
        }
        impl core::ops::MulAssign for $name {
            #[inline(always)]
            fn mul_assign(&mut self, rhs: Self) {
                *self = *self * rhs;
            }
        }
        impl core::ops::Neg for $name {
            type Output = Self;
            #[inline(always)]
            fn neg(self) -> Self {
                Self(Self::neg_limbs(self.0))
            }
        }

        impl $name {
            /// Reduce ≤ 256 bits given as 4 little-endian u64 words,
            /// returning the Montgomery form.
            pub const fn from_words_mod_order(words: [u64; 4]) -> Self {
                // two64_m = mont form of 2^64
                let mut two64 = [0u64; $limbs];
                two64[1] = 1;
                let two64_m = Self::to_mont(two64);
                let mut acc = [0u64; $limbs];
                let mut i = 4;
                while i > 0 {
                    i -= 1;
                    acc = Self::mont_mul(acc, two64_m);
                    let mut limb = [0u64; $limbs];
                    limb[0] = words[i];
                    let limb_m = Self::to_mont(limb);
                    acc = Self::add_limbs(acc, limb_m);
                }
                Self(acc)
            }
        }
    };
}

/// Emit the `PrimeField` trait implementation for a `mont_field!`-generated
/// type whose modulus fits in 32 bytes (4 limbs).
#[macro_export]
macro_rules! mont_prime_field {
    ($name:ident, $limbs:expr) => {
        impl $crate::PrimeField for $name {
            const LIMBS: usize = $limbs;
            const ONE: Self = Self(Self::R1);
            const GENERATOR: Self = Self(Self::to_mont(Self::GENERATOR_RAW));
            const TWO_ADICITY: u32 = Self::TWO_ADICITY;

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
            #[inline]
            fn invert(self) -> Option<Self> {
                Self::invert_limbs(self.0).map(Self)
            }
            #[inline]
            fn from_u64(v: u64) -> Self {
                let mut l = [0u64; $limbs];
                l[0] = v;
                Self::from_repr_limbs(l)
            }
            #[inline]
            fn to_le_bytes(self) -> [u8; 32] {
                let mut out = [0u8; 32];
                let r = self.to_repr();
                let mut i = 0;
                while i < $limbs {
                    out[i * 8..i * 8 + 8].copy_from_slice(&r[i].to_le_bytes());
                    i += 1;
                }
                out
            }
            #[inline]
            fn from_le_bytes_mod_order(bytes: &[u8; 32]) -> Self {
                let mut words = [0u64; 4];
                let mut i = 0;
                while i < 4 {
                    let mut w = [0u8; 8];
                    w.copy_from_slice(&bytes[i * 8..i * 8 + 8]);
                    words[i] = u64::from_le_bytes(w);
                    i += 1;
                }
                Self::from_words_mod_order(words)
            }
            #[inline]
            fn from_be_bytes_mod_order(bytes: &[u8]) -> Self {
                let mut le = [0u8; 32];
                let n = bytes.len().min(32);
                let mut i = 0;
                while i < n {
                    le[i] = bytes[n - 1 - i];
                    i += 1;
                }
                Self::from_le_bytes_mod_order(&le)
            }
        }
    };
}
