//! The BLS12-381 scalar field `Fr` — the field over which all KZG polynomial
//! arithmetic and Ethereum blob data live.

use crate::{mont_field, mont_prime_field, PrimeField};

mont_field!(
    Fr,
    4,
    [
        0xffffffff00000001,
        0x53bda402fffe5bfe,
        0x3339d80809a1d805,
        0x73eda753299d7d48
    ],
    0xfffffffeffffffff,
    [
        0x00000001fffffffe,
        0x5884b7fa00034802,
        0x998c4fefecbc4ff5,
        0x1824b159acc5056f
    ],
    [
        0xc999e990f3f29c6d,
        0x2b6cedcb87925c23,
        0x05d314967254398f,
        0x0748d9d99f59ff11
    ],
    [7, 0, 0, 0],
    32
);

mont_prime_field!(Fr, 4);

/// The scalar-field modulus limbs (little-endian), for canonicality
/// checks on raw encodings.
pub const FR_MODULUS: [u64; 4] = [
    0xffffffff00000001,
    0x53bda402fffe5bfe,
    0x3339d80809a1d805,
    0x73eda753299d7d48
];

impl Fr {
    /// The primitive root of unity `7^((r-1)/2^k)` of order `2^k`.
    ///
    /// `k` must be within the two-adicity (32 for BLS12-381's scalar field).
    pub fn root_of_unity(k: u32) -> Fr {
        debug_assert!(k > 0 && k <= 32);
        // r - 1, least-significant limbs first
        let mut e = [
            0xffffffff00000000u64,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ];
        // e >>= k
        let mut carry = 0u64;
        for limb in e.iter_mut().rev() {
            let nc = *limb << (64 - k);
            *limb = (*limb >> k) | carry;
            carry = nc;
        }
        Fr::GENERATOR.pow_word(&e)
    }

    /// Hash a byte string into `Fr` the way Ethereum's KZG spec does:
    /// `int(SHA256(data)) mod r`.
    pub fn hash_to_fr(data: &[u8]) -> Fr {
        let h = crate::sha256::sha256(data);
        Fr::from_be_bytes_mod_order(&h)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn from_int(v: u64) -> Fr {
        Fr::from_u64(v)
    }

    #[test]
    fn mont_rt_matches_const_path() {
        // The ADX/BMI2 clone must be semantically identical to the
        // portable const path (same algorithm, different codegen).
        let mut rng = crate::rng::ZodaRng::from_seed(*b"mont-rt-check-000000000000000000");
        for _ in 0..500 {
            let mut a = [0u64; 4];
            let mut b = [0u64; 4];
            for i in 0..4 {
                a[i] = rng.next_u64();
                b[i] = rng.next_u64();
            }
            // keep the top limb in range so values stay < 2^256 (they may
            // exceed r — the field functions handle unreduced operands)
            a[3] >>= 2;
            b[3] >>= 2;
            assert_eq!(Fr::mont_mul_rt(a, b), Fr::mont_mul(a, b));
        }
    }

    #[test]
    fn wide_mul_and_redc_match_mont_mul() {
        // mul_wide followed by mont_reduce_wide is the lazy-reduction
        // pipeline: it must agree exactly with the fused CIOS product.
        let mut rng = crate::rng::ZodaRng::from_seed(*b"mont-wide-chk-000000000000000000");
        for _ in 0..500 {
            let mut a = [0u64; 4];
            let mut b = [0u64; 4];
            for i in 0..4 {
                a[i] = rng.next_u64();
                b[i] = rng.next_u64();
            }
            // REDC validity needs T = a·b < r·R; limbs < r guarantee that,
            // and arbitrary < 2^256 products also satisfy it (r ~ 2^255,
            // R = 2^256 -> r·R ~ 2^511 > any 2^512 product? no: bound is
            // a·b < r·R; a,b < 2^256 gives a·b < 2^512 while r·R ~ 2^511,
            // so clamp the top limbs to stay in the valid range).
            a[3] &= 0x3fffffffffffffff;
            b[3] &= 0x3fffffffffffffff;
            let wide = Fr::mul_wide_rt(&a, &b);
            assert_eq!(Fr::mont_reduce_wide_rt(&wide), Fr::mont_mul(a, b));
            assert_eq!(Fr::mont_reduce_wide(&wide), Fr::mont_mul(a, b));
        }
    }

    #[test]
    fn basic_arithmetic() {
        let a = from_int(7);
        let b = from_int(5);
        assert_eq!(a + b, from_int(12));
        assert_eq!(a - b, from_int(2));
        assert_eq!(a * b, from_int(35));
        assert_eq!(-a + a, Fr::default());
        assert_eq!(a * Fr::ONE, a);
    }

    #[test]
    fn wraparound() {
        // r - 1 + 2 = 1
        let rm1 = -Fr::ONE;
        assert_eq!(rm1 + from_int(2), from_int(1));
    }

    #[test]
    fn inversion() {
        let a = from_int(123456789);
        let inv = a.invert().unwrap();
        assert_eq!(inv * a, Fr::ONE);
        assert!(Fr::default().invert().is_none());
    }

    #[test]
    fn roots_of_unity() {
        for k in 1..=13u32 {
            let w = Fr::root_of_unity(k);
            assert_eq!(w.pow_u64(1 << k), Fr::ONE, "order wrong for k={}", k);
            if k > 1 {
                assert_ne!(w.pow_u64(1 << (k - 1)), Fr::ONE, "not primitive k={}", k);
            }
        }
    }

    #[test]
    fn pow_and_bytes() {
        let a = from_int(3);
        assert_eq!(a.pow_u64(5), from_int(243));
        let bytes = from_int(0x0123).to_le_bytes();
        let b = Fr::from_le_bytes_mod_order(&bytes);
        assert_eq!(a + b, from_int(3) + from_int(0x123));
        // big-endian conversion round trip
        let mut be = [0u8; 32];
        be[31] = 0x42;
        assert_eq!(Fr::from_be_bytes_mod_order(&be), from_int(0x42));
    }

    #[test]
    fn mont_mul_matches_field_axioms() {
        let mut rng = crate::ZodaRng::from_seed(*b"zoda-math test seed 000000001000");
        for _ in 0..200 {
            let a = rng.next_fr(false);
            let b = rng.next_fr(false);
            let c = rng.next_fr(false);
            // distributivity / associativity — a full field-axiom check of
            // the Montgomery pipeline (any bug in reduction breaks these).
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!((a - b) + b, a);
            if !a.is_zero() {
                assert_eq!(a * a.invert().unwrap(), Fr::ONE);
            }
            // repr round trip
            let back = Fr::from_repr_limbs(a.to_repr());
            assert_eq!(back, a);
        }
    }

    #[test]
    fn root_of_unity_is_primitive() {
        // Cross-check the 8192-th root against the spec derivation:
        // PRIMITIVE_ROOT_OF_UNITY = 7, root = 7^((r-1)/8192)
        let w8192 = Fr::root_of_unity(13);
        assert_eq!(w8192.pow_u64(8192), Fr::ONE);
        assert_ne!(w8192.pow_u64(4096), Fr::ONE);
        // squares give the 4096-th root
        assert_eq!(Fr::root_of_unity(12), w8192 * w8192);
    }
}
