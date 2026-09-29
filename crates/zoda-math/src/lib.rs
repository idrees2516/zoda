//! # zoda-math
//!
//! Foundation library for the zoda data-availability stack: prime fields
//! (BLS12-381 `Fr`, Goldilocks, and the lattice modulus `Fq`), fast
//! number-theoretic transforms, subproduct-tree polynomial arithmetic,
//! SHA-256 / Keccak-256, Merkle trees, a deterministic CSPRNG and the
//! probability calculators used by sampling.

pub mod fr;
pub mod fq;
pub mod goldilocks;
pub mod keccak;
pub mod merkle;
pub mod mont;
pub mod ntt;
pub mod poly;
pub mod rng;
pub mod sha256;
pub mod stats;
pub mod u256;

pub use fq::Fq;
pub use fr::Fr;
pub use goldilocks::Goldilocks;
pub use ntt::FftDomain;
pub use rng::ZodaRng;

use core::fmt::Debug;
use core::hash::Hash;

/// The common interface implemented by every prime field in zoda.
///
/// Values are always kept fully reduced; every operation is branch-free with
/// respect to the data being processed.
pub trait PrimeField:
    Copy
    + Clone
    + PartialEq
    + Eq
    + Debug
    + Hash
    + Default
    + Send
    + Sync
    + 'static
    + core::ops::Add<Output = Self>
    + core::ops::Sub<Output = Self>
    + core::ops::Mul<Output = Self>
    + core::ops::Neg<Output = Self>
    + core::ops::AddAssign
    + core::ops::SubAssign
    + core::ops::MulAssign
{
    /// Canonical little-endian limb representation length.
    const LIMBS: usize;
    /// The multiplicative identity.
    const ONE: Self;
    /// A multiplicative generator of a large-order subgroup from which roots
    /// of unity can be derived.
    const GENERATOR: Self;
    /// The exponent of the largest power of two dividing `p - 1`.
    const TWO_ADICITY: u32;

    /// Field addition.
    fn add(self, rhs: Self) -> Self;
    /// Field subtraction.
    fn sub(self, rhs: Self) -> Self;
    /// Field multiplication.
    fn mul(self, rhs: Self) -> Self;
    /// Additive inverse.
    fn neg(self) -> Self;
    /// Multiplicative inverse; `None` for zero.
    fn invert(self) -> Option<Self>;
    /// Square.
    #[inline]
    fn square(self) -> Self {
        self * self
    }
    /// Exponentiation by a 64-bit exponent.
    fn pow_u64(self, exp: u64) -> Self {
        let mut res = Self::ONE;
        let mut e = exp;
        let mut base = self;
        while e != 0 {
            if e & 1 == 1 {
                res = res * base;
            }
            base = base.square();
            e >>= 1;
        }
        res
    }
    /// Exponentiation by a little-endian 4×64-bit exponent.
    fn pow_word(&self, exp: &[u64; 4]) -> Self {
        let mut res = Self::ONE;
        for i in (0..4).rev() {
            if exp[i] == 0 && i == 0 {
                continue;
            }
            for b in (0..64).rev() {
                res = res.square();
                if (exp[i] >> b) & 1 == 1 {
                    res = res * *self;
                }
            }
        }
        res
    }
    /// True if this is the additive identity.
    #[inline]
    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
    /// The additive identity.
    #[inline]
    fn zero() -> Self {
        Self::default()
    }
    /// Construct from a small unsigned integer.
    fn from_u64(v: u64) -> Self;
    /// Canonical 32-byte little-endian encoding of the reduced value.
    fn to_le_bytes(self) -> [u8; 32];
    /// Reduce 32 little-endian bytes modulo the field order.
    fn from_le_bytes_mod_order(bytes: &[u8; 32]) -> Self;
    /// Interpret a big-endian integer given as bytes, reduced mod p.
    fn from_be_bytes_mod_order(bytes: &[u8]) -> Self;
    /// Serialise to hex (canonical LE bytes) — for logs and tests.
    fn to_hex(self) -> String {
        let bytes = self.to_le_bytes();
        let mut s = String::with_capacity(66);
        s.push_str("0x");
        for b in bytes.iter().rev() {
            s.push_str(&format!("{:02x}", b));
        }
        s
    }
}

/// Montgomery batch inversion (single inversion + 3(n-1) multiplications).
///
/// Zero entries are mapped to zero.
pub fn batch_invert<F: PrimeField>(v: &mut [F]) {
    if v.is_empty() {
        return;
    }
    let mut acc = F::ONE;
    let mut prefix = vec![F::ONE; v.len()];
    for i in 0..v.len() {
        prefix[i] = acc;
        if !v[i].is_zero() {
            acc = acc * v[i];
        }
    }
    let mut inv = match acc.invert() {
        Some(x) => x,
        None => return, // everything was zero
    };
    for i in (0..v.len()).rev() {
        if v[i].is_zero() {
            continue;
        }
        let vi = v[i];
        v[i] = inv * prefix[i];
        inv = inv * vi;
    }
}
