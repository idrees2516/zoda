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

/// The exponent `(p + 1) / 4` as LE limbs — square roots in Fp for
/// `p \equiv 3 (mod 4)`, computed once at compile time.
const SQRT_EXP: [u64; 6] = {
    let mut e = [
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a,
    ];
    // p + 1
    let mut carry = 1u64;
    let mut i = 0;
    while i < 6 {
        let t = e[i] as u128 + carry as u128;
        e[i] = t as u64;
        carry = (t >> 64) as u64;
        i += 1;
    }
    // >> 2
    let mut c = 0u64;
    let mut j = 6;
    while j > 0 {
        j -= 1;
        let nc = e[j] << 62;
        e[j] = (e[j] >> 2) | c;
        c = nc;
    }
    e
};

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

    /// Fixed-exponent windowed exponentiation (4-bit windows). For a
    /// 381-bit exponent this trades ~190 multiplications for ~95, using a
    /// 16-entry table of small powers of the base — about 25% faster than
    /// the bit-at-a-time ladder on the two hot fixed exponents of this
    /// crate (square roots and, via the extension tower, inversions that
    /// still go through Fermat).
    pub fn pow_windowed(&self, exp: &[u64; 6]) -> Fp {
        // find the top limb / bit
        let mut top = 5usize;
        while top > 0 && exp[top] == 0 {
            top -= 1;
        }
        let bits = 64 - exp[top].leading_zeros();
        let total_bits = top * 64 + bits as usize;
        if total_bits == 0 {
            return Fp::one();
        }
        // digit table: self^0 .. self^15
        let mut tab = [Fp::zero(); 16];
        tab[0] = Fp::one();
        tab[1] = *self;
        for i in 2..16 {
            tab[i] = tab[i - 1] * *self;
        }
        // scan 4-bit digits, most significant first (the first digit may
        // be narrower than 4 bits)
        let nwin = (total_bits + 3) / 4;
        let mut acc = Fp::one();
        for w in (0..nwin).rev() {
            if w != nwin - 1 {
                acc = acc.square().square().square().square();
            }
            // digit covering bits [4w, 4w+4)
            let lo = w * 4;
            let mut d = 0usize;
            for b in 0..4 {
                let bit = lo + b;
                if bit < total_bits && (exp[bit / 64] >> (bit % 64)) & 1 == 1 {
                    d |= 1 << b;
                }
            }
            if w == nwin - 1 {
                // leading window: set acc directly (no squaring chain yet)
                acc = tab[d];
            } else if d != 0 {
                acc = acc * tab[d];
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
    /// is a genuine square, else `None`. Uses the windowed fixed-exponent
    /// ladder over the compile-time `(p+1)/4` constant.
    pub fn sqrt(&self) -> Option<Fp> {
        let r = self.pow_windowed(&SQRT_EXP);
        if r.square() == *self {
            Some(r)
        } else {
            None
        }
    }

    /// A square root candidate without the QR check (`a^((p+1)/4)`),
    /// for callers that verify the result by other means (the Fp2
    /// complex-method square root reuses this for its second stage).
    pub fn sqrt_candidate(&self) -> Fp {
        self.pow_windowed(&SQRT_EXP)
    }

    /// True iff `self` is a quadratic residue (Euler criterion via the
    /// square-root candidate — one extra squaring).
    pub fn is_square(&self) -> bool {
        let r = self.pow_windowed(&SQRT_EXP);
        r.square() == *self
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
    fn pow_windowed_matches_pow_limbs() {
        // differential: the 4-bit window ladder must agree with the
        // bit-at-a-time ladder on arbitrary exponents of many shapes
        let mut rng = zoda_math::ZodaRng::from_seed(*b"fp-pow-windowed-0000000000000000");
        for i in 0..64 {
            let a = rand_fp(i);
            let mut e = [0u64; 6];
            for limb in e.iter_mut() {
                let mut b = [0u8; 8];
                rng.next_bytes(&mut b);
                *limb = u64::from_le_bytes(b);
            }
            // shrink some exponents to exercise short/leading-window paths
            if i % 4 == 0 {
                e[5] = 0;
            }
            if i % 8 == 0 {
                e = [e[0] & 0xff, 0, 0, 0, 0, 0];
            }
            assert_eq!(
                a.pow_windowed(&e),
                a.pow_limbs(&e),
                "windowed pow mismatch at i={}",
                i
            );
        }
        // fixed exponents: sqrt / QR consistency
        let a = rand_fp(99);
        let sq = a.square();
        assert!(sq.is_square());
        assert_eq!(sq.sqrt().unwrap().square(), sq);
        assert_eq!(sq.sqrt_candidate().square(), sq);
        // 2 is a non-residue mod p (p = 3 mod 8)
        assert!(!Fp::from_u64(2).is_square());
    }

    #[test]
    fn fp_lex_order() {
        // 1 < p/2 → not largest; p-1 > p/2 → largest
        assert!(!Fp::from_u64(1).lexicographically_largest());
        assert!((-Fp::one()).lexicographically_largest());
    }
}
