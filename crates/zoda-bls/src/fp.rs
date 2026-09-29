//! BLS12-381 base field `Fp` (381 bits, 6×u64 limbs, Montgomery form).

use zoda_math::mont_field;

mont_field!(
    Fp,
    6,
    [
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a
    ],
    0x89f3fffcfffcfffd,
    [
        0x760900000002fffd,
        0xebf4000bc40c0002,
        0x5f48985753c758ba,
        0x77ce585370525745,
        0x5c071a97a256ec6d,
        0x15f65ec3fa80e493
    ],
    [
        0xf4df1f341c341746,
        0x0a76e6a609d104f1,
        0x8de5476c4c95b6d5,
        0x67eb88a9939d83c0,
        0x9a793e85b519952d,
        0x11988fe592cae3aa
    ],
    [
        0x07d0_0000_0000_000f, // placeholder generator (unused: no FFTs over Fp)
        0,
        0,
        0,
        0,
        0
    ],
    1
);

impl Fp {
    /// Additive identity.
    #[inline]
    pub const fn zero() -> Fp {
        Fp([0; 6])
    }
    /// Multiplicative identity (Montgomery form of 1).
    #[inline]
    pub const fn one() -> Fp {
        Fp(Self::R1)
    }
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0 == [0u64; 6]
    }
    #[inline]
    pub fn invert(self) -> Option<Fp> {
        Self::invert_limbs(self.0).map(Fp)
    }
    #[inline]
    pub fn square(&self) -> Fp {
        *self * *self
    }
    /// Construct from a u64 (standard integer).
    #[inline]
    pub fn from_u64(v: u64) -> Fp {
        let mut l = [0u64; 6];
        l[0] = v;
        Fp::from_repr_limbs(l)
    }
    /// Construct from canonical limbs (standard integer form).
    #[inline]
    pub const fn from_limbs(l: [u64; 6]) -> Fp {
        Fp::from_repr_limbs(l)
    }
    /// Exponentiation by an arbitrary-length little-endian limb exponent.
    pub fn pow_limbs(&self, exp: &[u64]) -> Fp {
        let mut acc = Fp::one();
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

    /// Exponentiation by a 32-bit exponent.
    pub fn pow_u32(&self, mut e: u32) -> Fp {
        let mut acc = Fp::one();
        let mut base = *self;
        while e > 0 {
            if e & 1 == 1 {
                acc = acc * base;
            }
            base = base.square();
            e >>= 1;
        }
        acc
    }

    /// Square root in `Fp` (`p ≡ 3 (mod 4)`): returns `a^((p+1)/4)` when it
    /// is a genuine square, else `None`.
    pub fn sqrt(&self) -> Option<Fp> {
        // exponent = (p + 1) / 4
        let mut e = Self::MODULUS;
        // p + 1
        let mut carry = 1u64;
        for limb in e.iter_mut() {
            let t = (*limb as u128) + (carry as u128);
            *limb = t as u64;
            carry = (t >> 64) as u64;
        }
        // >> 2
        let mut c = 0u64;
        for limb in e.iter_mut().rev() {
            let nc = *limb << 62;
            *limb = (*limb >> 2) | c;
            c = nc;
        }
        let r = self.pow_limbs(&e);
        if r.square() == *self {
            Some(r)
        } else {
            None
        }
    }

    /// True iff `x > (p-1)/2` — the canonical "largest" sign used by point
    /// compression. Operates on the standard (non-Montgomery) form.
    pub fn lexicographically_largest(&self) -> bool {
        let r = self.to_repr();
        // x largest ⇔ 2x >= p  (as 384-bit integers)
        let mut carry = 0u128;
        let mut out = [0u64; 7];
        for i in 0..6 {
            let t = (r[i] as u128) * 2 + carry;
            out[i] = t as u64;
            carry = t >> 64;
        }
        out[6] = carry as u64;
        if out[6] != 0 {
            return true;
        }
        for i in (0..6).rev() {
            if out[i] > Self::MODULUS[i] {
                return true;
            }
            if out[i] < Self::MODULUS[i] {
                return false;
            }
        }
        true // 2x == p — impossible for reduced x (p odd), but be safe
    }

    /// Random-ish Fp from 32 bytes (reduced; test helper).
    pub fn from_le_bytes32(b: &[u8; 32]) -> Fp {
        let mut words = [0u64; 4];
        for i in 0..4 {
            let mut w = [0u8; 8];
            w.copy_from_slice(&b[i * 8..i * 8 + 8]);
            words[i] = u64::from_le_bytes(w);
        }
        Fp::from_words_mod_order(words)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand_fp(seed: u64) -> Fp {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"fp-test-seed-0000000000000000000");
        let _ = seed;
        let mut b = [0u8; 32];
        rng.next_bytes(&mut b);
        Fp::from_le_bytes32(&b)
    }

    #[test]
    fn fp_arithmetic() {
        for _ in 0..100 {
            let a = rand_fp(0);
            let b = rand_fp(1);
            let c = rand_fp(2);
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            if !a.is_zero() {
                assert_eq!(a * a.invert().unwrap(), Fp::one());
            }
        }
        assert_eq!(Fp::zero().is_zero(), true);
        assert_eq!(Fp::one() * Fp::one(), Fp::one());
    }

    #[test]
    fn fp_sqrt() {
        // 4 is a square (2^2)
        let four = Fp::from_u64(4);
        let r = four.sqrt().unwrap();
        assert!(r == Fp::from_u64(2) || r == -Fp::from_u64(2));
        // random squares are squares
        for i in 0..25 {
            let a = rand_fp(i);
            let sq = a.square();
            let r = sq.sqrt().expect("squares are QRs");
            assert!(r == a || r == -a);
        }
    }

    #[test]
    fn fp_lex_order() {
        // 1 < p/2 → not largest; p-1 > p/2 → largest
        assert!(!Fp::from_u64(1).lexicographically_largest());
        assert!((-Fp::one()).lexicographically_largest());
    }
}
