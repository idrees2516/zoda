//! `Fp2 = Fp[u] / (u² + 1)` — the quadratic extension of BLS12-381's base
//! field (Karatsuba arithmetic with **lazy reduction** — three wide
//! products and two Montgomery reductions per multiplication — norm-based
//! inversion and square roots).

use crate::fp::Fp;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// `4·p²` as a 12-limb integer: the lazy-reduction positivity offset.
/// `p² ≡ 0 (mod p)`, so adding this constant keeps every wide intermediate
/// strictly positive (no borrow underflow) while preserving the residue
/// class. Bounds: every Karatsuba intermediate stays below `8·p² < p·R`, the
/// valid range for a single wide Montgomery reduction to a fully reduced
/// result.
const FOUR_P2: [u64; 12] = {
    let p2 = Fp::mul_wide(&Fp::MODULUS, &Fp::MODULUS);
    let mut out = [0u64; 12];
    let mut carry = 0u128;
    let mut i = 0;
    while i < 12 {
        let v = ((p2[i] as u128) << 2) | carry;
        out[i] = v as u64;
        carry = v >> 64;
        i += 1;
    }
    out
};

/// 6-limb addition (inputs are reduced Montgomery limbs; the sum of two
/// values below `p` stays below `2p < 2^384`, so no carry can escape).
#[inline(always)]
fn add6(a: [u64; 6], b: [u64; 6]) -> [u64; 6] {
    let mut out = [0u64; 6];
    let mut carry = 0u128;
    for i in 0..6 {
        let t = a[i] as u128 + b[i] as u128 + carry;
        out[i] = t as u64;
        carry = t >> 64;
    }
    debug_assert_eq!(carry, 0);
    out
}

/// 6-limb subtraction on Montgomery limbs, mod-`p` transparent: if
/// `a < b` the result is `a − b + p` (a multiple-of-`p` adjustment that
/// keeps the value non-negative and the residue class unchanged).
#[inline(always)]
fn sub6_mod_p(a: [u64; 6], b: [u64; 6]) -> [u64; 6] {
    let mut out = [0u64; 6];
    let mut borrow = 0u128;
    for i in 0..6 {
        let t = (a[i] as u128).wrapping_sub(b[i] as u128 + borrow);
        out[i] = t as u64;
        borrow = (t >> 64) & 1;
    }
    if borrow != 0 {
        let mut carry = 0u128;
        for i in 0..6 {
            let t = out[i] as u128 + Fp::MODULUS[i] as u128 + carry;
            out[i] = t as u64;
            carry = t >> 64;
        }
    }
    out
}

/// 12-limb addition (lazy intermediates stay below `8p² < 2^765` — no carry).
#[inline(always)]
fn add12(a: &[u64; 12], b: &[u64; 12]) -> [u64; 12] {
    let mut out = [0u64; 12];
    let mut carry = 0u128;
    for i in 0..12 {
        let t = a[i] as u128 + b[i] as u128 + carry;
        out[i] = t as u64;
        carry = t >> 64;
    }
    debug_assert_eq!(carry, 0);
    out
}

/// 12-limb subtraction; `a ≥ b` must hold (guaranteed by the `4p²` offset).
#[inline(always)]
fn sub12(a: &[u64; 12], b: &[u64; 12]) -> [u64; 12] {
    let mut out = [0u64; 12];
    let mut borrow = 0u128;
    for i in 0..12 {
        let t = (a[i] as u128).wrapping_sub(b[i] as u128 + borrow);
        out[i] = t as u64;
        borrow = (t >> 64) & 1;
    }
    debug_assert_eq!(borrow, 0);
    out
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct Fp2 {
    pub c0: Fp,
    pub c1: Fp,
}

/// `1/2` in Fp (Montgomery form), computed at compile time as
/// `(p+1)/2` — valid because `2·(p+1)/2 = p+1 ≡ 1 (mod p)`.
const INV_TWO: Fp = Fp(Fp::to_mont({
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
    // >> 1
    let mut c = 0u64;
    let mut j = 6;
    while j > 0 {
        j -= 1;
        let nc = e[j] << 63;
        e[j] = (e[j] >> 1) | c;
        c = nc;
    }
    e
}));

/// Batch inversion over Fp (Montgomery's trick — one Fermat inversion
/// plus 3(n−1) multiplications), preserving zeros as zeros.
fn batch_invert_fp(vals: &[Fp]) -> Vec<Fp> {
    // a local copy of zoda_math::batch_invert specialised to Fp's Copy type;
    // the generic zoda_math helper is written against Fr-style fields
    let mut prods = Vec::with_capacity(vals.len());
    let mut acc = Fp::one();
    for v in vals {
        prods.push(acc);
        if !v.is_zero() {
            acc = acc * *v;
        }
    }
    // acc holds the running product of non-zero entries; invert once
    let mut inv = match acc.invert() {
        Some(i) => i,
        None => Fp::one(), // all-zero slice
    };
    let mut out = vec![Fp::zero(); vals.len()];
    for i in (0..vals.len()).rev() {
        if vals[i].is_zero() {
            continue; // inv unchanged: zero never entered the product
        }
        out[i] = inv * prods[i];
        inv = inv * vals[i];
    }
    out
}

/// One element's progress through the batch complex-method square root.
enum SqrtStage {
    /// The root is complete (zero or pure-real input).
    Done(Fp2),
    /// Real part fixed; imaginary part is `β/(2u)` (inversion deferred).
    Pair(Fp, Fp),
    /// No square root exists.
    NonSquare,
}

impl Fp2 {
    pub const fn zero() -> Fp2 {
        Fp2 {
            c0: Fp::zero(),
            c1: Fp::zero(),
        }
    }
    pub const fn one() -> Fp2 {
        Fp2 {
            c0: Fp::one(),
            c1: Fp::zero(),
        }
    }
    pub fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// The cubic non-residue ξ = 1 + u used to build Fp6 (standard
    /// BLS12-381 tower: Fp6 = Fp2[v]/(v³ − ξ), Fp12 = Fp6[w]/(w² − v)).
    pub fn xi() -> Fp2 {
        Fp2 {
            c0: Fp::one(),
            c1: Fp::one(),
        }
    }

    pub fn new(c0: Fp, c1: Fp) -> Fp2 {
        Fp2 { c0, c1 }
    }

    /// Conjugate (Frobenius on Fp2).
    #[inline]
    pub fn conjugate(&self) -> Fp2 {
        Fp2 {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    /// Norm a·ā = c0² + c1² ∈ Fp.
    #[inline]
    pub fn norm(&self) -> Fp {
        self.c0 * self.c0 + self.c1 * self.c1
    }

    #[inline]
    pub fn square(&self) -> Fp2 {
        // Complex squaring with lazy reduction:
        //   c0 = (a+b)(a−b) = a² − b²,   c1 = 2ab
        // two wide products + two Montgomery reductions (down from two
        // fully reduced products plus additions).
        let a = self.c0.0;
        let b = self.c1.0;
        let s = add6(a, b);
        let d = sub6_mod_p(a, b);
        let t0 = Fp::mul_wide_rt(&s, &d); // (a+b)(a−b), < 2p²
        let t1 = Fp::mul_wide_rt(&a, &b); // ab, < p²
        // 2ab: shift the wide product left by one bit (still < 2p²)
        let mut t1s = [0u64; 12];
        let mut carry = 0u128;
        for i in 0..12 {
            let v = ((t1[i] as u128) << 1) | carry;
            t1s[i] = v as u64;
            carry = v >> 64;
        }
        debug_assert_eq!(carry, 0);
        Fp2 {
            c0: Fp(Fp::mont_reduce_wide_rt(&t0)),
            c1: Fp(Fp::mont_reduce_wide_rt(&t1s)),
        }
    }

    pub fn invert(&self) -> Option<Fp2> {
        let n = self.norm().invert()?;
        Some(Fp2 {
            c0: self.c0 * n,
            c1: -(self.c1 * n),
        })
    }

    /// Batch inversion over Fp2 (Montgomery's trick): one norm-based
    /// inversion — a single windowed Fp exponentiation — plus 3(n−1) Fp2
    /// multiplications, amortising the dominant cost across the whole
    /// slice. Zeros pass through as zeros (they never enter the running
    /// product, so the shared inverse stays valid for the other entries).
    pub fn invert_batch(vals: &[Fp2]) -> Vec<Fp2> {
        let mut prods = Vec::with_capacity(vals.len());
        let mut acc = Fp2::one();
        for v in vals {
            prods.push(acc);
            if !v.is_zero() {
                acc = acc * *v;
            }
        }
        let mut inv = match acc.invert() {
            Some(i) => i,
            None => Fp2::one(), // all-zero slice
        };
        let mut out = vec![Fp2::zero(); vals.len()];
        for i in (0..vals.len()).rev() {
            if vals[i].is_zero() {
                continue;
            }
            out[i] = inv * prods[i];
            inv = inv * vals[i];
        }
        out
    }

    /// Square root in Fp2 via the complex method for `p \equiv 3 (mod 4)`
    /// — two Fp exponentiations by the compile-time constant `(p+1)/4`
    /// (norm, then real part) instead of one 762-bit Fp2 exponentiation:
    /// roughly 4× faster, and the first stage doubles as the quadratic-
    /// residue test (`a` is a square in Fp2 ⇔ its norm is a square in Fp,
    /// because the norm induces an isomorphism Fp2*/(Fp2*)² → Fp*/(Fp*)²).
    ///
    /// Returns `None` exactly when `self` is a non-square (or zero with a
    /// vanishing component arrangement that has no root).
    pub fn sqrt(&self) -> Option<Fp2> {
        if self.is_zero() {
            return Some(Fp2::zero());
        }
        if self.c1.is_zero() {
            // pure real: sqrt in Fp or i·sqrt(−x)
            return match self.c0.sqrt() {
                Some(r) => Some(Fp2 { c0: r, c1: Fp::zero() }),
                None => {
                    // -c0 must be a square (−1 is a non-residue in Fp)
                    let r = (-self.c0).sqrt()?;
                    Some(Fp2 {
                        c0: Fp::zero(),
                        c1: r,
                    })
                }
            };
        }
        self.sqrt_complex()
    }

    /// The complex-method core for `c1 ≠ 0`. Stages:
    ///   1. `t = sqrt(c0² + c1²)` — `None` ⇒ the norm is a non-residue ⇒
    ///      `self` is a non-square in Fp2 (isomorphism on square classes).
    ///   2. `u = sqrt((c0 ± t)/2)` — exactly one sign is a residue when
    ///      `self` is a square (the two candidate real parts multiply to
    ///      `−(c1/2)²`, and −1 is a non-residue in Fp).
    ///   3. `v = c1 / (2u)`.
    fn sqrt_complex(&self) -> Option<Fp2> {
        let alpha = self.c0;
        let beta = self.c1;
        let norm = alpha * alpha + beta * beta;
        let t = norm.sqrt()?;
        let half = INV_TWO;
        let mut u = None;
        for sign in [1i8, -1i8] {
            let cand = (alpha + if sign > 0 { t } else { -t }) * half;
            if let Some(r) = cand.sqrt() {
                if !r.is_zero() {
                    u = Some(r);
                    break;
                }
            }
        }
        let u = u?;
        let v = beta * (u + u).invert()?;
        let candidate = Fp2 { c0: u, c1: v };
        if candidate.square() == *self {
            Some(candidate)
        } else {
            None
        }
    }

    /// Batch square root over Fp2 (complex method): one batch inversion
    /// of the `2u` denominators across the whole slice replaces a Fermat
    /// inversion per element. `None` entries mark non-squares. This is the
    /// batch hash-to-curve workhorse.
    pub fn sqrt_batch(vals: &[Fp2]) -> Vec<Option<Fp2>> {
        // stage 1+2 for every element; complete results (zero and
        // pure-real inputs) bypass the shared inversion entirely
        let mut stages: Vec<SqrtStage> = Vec::with_capacity(vals.len());
        let mut denoms: Vec<Fp> = Vec::with_capacity(vals.len());
        for v in vals {
            match v.sqrt_stage() {
                SqrtStage::Done(z) => stages.push(SqrtStage::Done(z)),
                SqrtStage::Pair(beta, u) => {
                    denoms.push(u + u);
                    stages.push(SqrtStage::Pair(beta, u));
                }
                SqrtStage::NonSquare => stages.push(SqrtStage::NonSquare),
            }
        }
        // one batch inversion of all denominators
        let invs = batch_invert_fp(&denoms);
        let mut out = Vec::with_capacity(vals.len());
        let mut k = 0usize;
        for (v, st) in vals.iter().zip(stages.iter()) {
            match st {
                SqrtStage::Done(z) => out.push(Some(*z)),
                SqrtStage::NonSquare => out.push(None),
                SqrtStage::Pair(beta, u) => {
                    let c = Fp2 {
                        c0: *u,
                        c1: *beta * invs[k],
                    };
                    k += 1;
                    debug_assert_eq!(c.square(), *v);
                    out.push(Some(c));
                }
            }
        }
        out
    }

    /// Stage 1+2 of the complex method:
    /// * `Done(z)` — the square root is already complete (zero input, or a
    ///   pure-real input whose root is real or pure-imaginary);
    /// * `Pair(β, u)` — the real part `u` is fixed, the imaginary part is
    ///   `β/(2u)` (the caller batches the inversion);
    /// * `NonSquare` — no root exists.
    fn sqrt_stage(&self) -> SqrtStage {
        if self.is_zero() {
            return SqrtStage::Done(Fp2::zero());
        }
        if self.c1.is_zero() {
            if let Some(r) = self.c0.sqrt() {
                return SqrtStage::Done(Fp2 { c0: r, c1: Fp::zero() });
            }
            if let Some(r) = (-self.c0).sqrt() {
                return SqrtStage::Done(Fp2 { c0: Fp::zero(), c1: r });
            }
            return SqrtStage::NonSquare;
        }
        let alpha = self.c0;
        let beta = self.c1;
        let norm = alpha * alpha + beta * beta;
        let t = match norm.sqrt() {
            Some(t) => t,
            None => return SqrtStage::NonSquare,
        };
        for sign in [1i8, -1i8] {
            let cand = (alpha + if sign > 0 { t } else { -t }) * INV_TWO;
            if let Some(r) = cand.sqrt() {
                if !r.is_zero() {
                    return SqrtStage::Pair(beta, r);
                }
            }
        }
        SqrtStage::NonSquare
    }

    /// Exponentiation by an arbitrary-length LE limb exponent.
    pub fn pow_limbs(&self, exp: &[u64]) -> Fp2 {
        let mut acc = Fp2::one();
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

    /// sgn0 per RFC 9380 (sign of c0, then of c1 if c0 is zero).
    pub fn sgn0(&self) -> bool {
        // sign_0(x) = x mod 2 (parity of the standard representation)
        let c0_parity = {
            let r = self.c0.to_repr();
            r[0] & 1 == 1
        };
        let c1_parity = {
            let r = self.c1.to_repr();
            r[0] & 1 == 1
        };
        let c0_zero = self.c0.is_zero();
        c0_parity || (c0_zero && c1_parity)
    }

    /// Canonical big-endian byte encoding (c1 || c0), 96 bytes.
    pub fn to_be_bytes(&self) -> [u8; 96] {
        let mut out = [0u8; 96];
        let c1 = self.c1.to_repr();
        let c0 = self.c0.to_repr();
        // repr limbs are little-endian; emit most significant limb first
        for i in 0..6 {
            out[i * 8..i * 8 + 8].copy_from_slice(&c1[5 - i].to_be_bytes());
        }
        for i in 0..6 {
            out[48 + i * 8..48 + i * 8 + 8].copy_from_slice(&c0[5 - i].to_be_bytes());
        }
        out
    }

    /// Multiply by (u + 1) (used for Fp6 non-residue power computations).
    #[inline]
    pub fn mul_by_1_plus_u(&self) -> Fp2 {
        // (a + bu)(1 + u) = (a − b) + (a + b)u
        Fp2 {
            c0: self.c0 - self.c1,
            c1: self.c0 + self.c1,
        }
    }
}

impl Add for Fp2 {
    type Output = Fp2;
    #[inline(always)]
    fn add(self, rhs: Fp2) -> Fp2 {
        Fp2 {
            c0: self.c0 + rhs.c0,
            c1: self.c1 + rhs.c1,
        }
    }
}
impl Sub for Fp2 {
    type Output = Fp2;
    #[inline(always)]
    fn sub(self, rhs: Fp2) -> Fp2 {
        Fp2 {
            c0: self.c0 - rhs.c0,
            c1: self.c1 - rhs.c1,
        }
    }
}
impl Neg for Fp2 {
    type Output = Fp2;
    #[inline(always)]
    fn neg(self) -> Fp2 {
        Fp2 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}
impl Mul for Fp2 {
    type Output = Fp2;
    #[inline(always)]
    fn mul(self, rhs: Fp2) -> Fp2 {
        // Lazy-reduction Karatsuba (u² = −1):
        //   c0 = a0b0 − a1b1,  c1 = (a0+a1)(b0+b1) − a0b0 − a1b1
        // Three wide (unreduced) products and two Montgomery reductions
        // instead of three fully reduced multiplications — one third of the
        // reduction work removed from the hottest extension-field op.
        //
        // Positivity: T − U ≥ −p², V − T − U ≥ −2p², and the `4p²` offset
        // lifts every wide value into (0, 8p²) ⊂ (0, p·R) — the REDC-valid
        // range — without touching the residue class (p² ≡ 0 mod p).
        let a0 = self.c0.0;
        let a1 = self.c1.0;
        let b0 = rhs.c0.0;
        let b1 = rhs.c1.0;
        let sa = add6(a0, a1);
        let sb = add6(b0, b1);
        let t = Fp::mul_wide_rt(&a0, &b0);
        let u = Fp::mul_wide_rt(&a1, &b1);
        let v = Fp::mul_wide_rt(&sa, &sb);
        // c0 = T + 4p² − U
        let c0w = sub12(&add12(&t, &FOUR_P2), &u);
        // c1 = V + 4p² − T − U
        let c1w = sub12(&sub12(&add12(&v, &FOUR_P2), &t), &u);
        Fp2 {
            c0: Fp(Fp::mont_reduce_wide_rt(&c0w)),
            c1: Fp(Fp::mont_reduce_wide_rt(&c1w)),
        }
    }
}
impl AddAssign for Fp2 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Fp2) {
        *self = *self + rhs;
    }
}
impl SubAssign for Fp2 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: Fp2) {
        *self = *self - rhs;
    }
}
impl MulAssign for Fp2 {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: Fp2) {
        *self = *self * rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The pre-lazy reference: fully reduced Karatsuba multiplication.
    fn fp2_mul_reference(a: Fp2, b: Fp2) -> Fp2 {
        let v0 = a.c0 * b.c0;
        let v1 = a.c1 * b.c1;
        let t = (a.c0 + a.c1) * (b.c0 + b.c1);
        Fp2 { c0: v0 - v1, c1: t - v0 - v1 }
    }

    /// The pre-lazy reference squaring.
    fn fp2_square_reference(a: Fp2) -> Fp2 {
        let ab = a.c0 * a.c1;
        Fp2 { c0: (a.c0 + a.c1) * (a.c0 - a.c1), c1: ab + ab }
    }

    #[test]
    fn lazy_mul_matches_reference() {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"fp2-lazy-mul-0000000000000000000");
        let next = |rng: &mut zoda_math::ZodaRng| -> Fp2 {
            let mut b0 = [0u8; 32];
            let mut b1 = [0u8; 32];
            rng.next_bytes(&mut b0);
            rng.next_bytes(&mut b1);
            Fp2::new(Fp::from_le_bytes32(&b0), Fp::from_le_bytes32(&b1))
        };
        for i in 0..2000 {
            let a = next(&mut rng);
            let b = next(&mut rng);
            assert_eq!(a * b, fp2_mul_reference(a, b), "lazy mul mismatch at {}", i);
            assert_eq!(a.square(), fp2_square_reference(a), "lazy sqr mismatch at {}", i);
            assert_eq!(a * a, fp2_mul_reference(a, a), "self-mul mismatch at {}", i);
        }
        // edge cases: zero, one, u, negatives, one operand zero
        let z = Fp2::zero();
        let o = Fp2::one();
        let u = Fp2::new(Fp::zero(), Fp::one());
        assert_eq!(z * z, z);
        assert_eq!(z * u, z);
        assert_eq!(o * u, u);
        assert_eq!(u * u, -o);
        assert_eq!(u.square(), -o);
        let a = rand_fp2(0);
        assert_eq!(a * o, a);
        assert_eq!((-a) * o, -a);
        assert_eq!(a * z, z);
    }

    #[test]
    fn lazy_offsets_are_p_transparent() {
        // FOUR_P2 must be exactly 4·p² (checked by construction, but pin
        // the top limbs against p² << 2 computed independently)
        let p2 = Fp::mul_wide(&Fp::MODULUS, &Fp::MODULUS);
        let mut expect = [0u64; 12];
        let mut carry = 0u128;
        for i in 0..12 {
            let v = ((p2[i] as u128) << 2) | carry;
            expect[i] = v as u64;
            carry = v >> 64;
        }
        assert_eq!(carry, 0);
        assert_eq!(FOUR_P2, expect);
    }

    fn rand_fp2(i: u64) -> Fp2 {
        let mut rng = zoda_math::ZodaRng::from_seed(*b"fp2-test-seed-000000000000000000");
        let _ = i;
        let mut b0 = [0u8; 32];
        let mut b1 = [0u8; 32];
        rng.next_bytes(&mut b0);
        rng.next_bytes(&mut b1);
        Fp2::new(Fp::from_le_bytes32(&b0), Fp::from_le_bytes32(&b1))
    }

    #[test]
    fn fp2_axioms() {
        for _ in 0..100 {
            let a = rand_fp2(0);
            let b = rand_fp2(1);
            let c = rand_fp2(2);
            assert_eq!((a * b) * c, a * (b * c));
            assert_eq!(a * (b + c), a * b + a * c);
            assert_eq!(a.square(), a * a);
            if !a.is_zero() {
                assert_eq!(a * a.invert().unwrap(), Fp2::one());
            }
        }
    }

    #[test]
    fn fp2_sqrt() {
        for i in 0..25 {
            let a = rand_fp2(i);
            let sq = a.square();
            let r = sq.sqrt().expect("squares are QR");
            assert!(r == a || r == -a, "sqrt failed for i={}", i);
        }
        // i (u) has sqrt since −1 is a square in Fp2
        let u = Fp2::new(Fp::zero(), Fp::one());
        let r = u.sqrt().unwrap();
        assert_eq!(r.square(), u);
    }

    #[test]
    fn fp2_invert_batch_matches_single() {
        let mut vals: Vec<Fp2> = Vec::new();
        for i in 0..32 {
            vals.push(rand_fp2(i));
        }
        vals.push(Fp2::zero()); // zeros must pass through untouched
        vals.push(Fp2::one());
        let batch = Fp2::invert_batch(&vals);
        for (v, b) in vals.iter().zip(batch.iter()) {
            if v.is_zero() {
                assert!(b.is_zero());
            } else {
                assert_eq!(*b, v.invert().unwrap());
            }
        }
        // all-zero slice must not panic
        assert!(Fp2::invert_batch(&[Fp2::zero(), Fp2::zero()])
            .iter()
            .all(|x| x.is_zero()));
    }

    #[test]
    fn fp2_sqrt_batch_matches_single() {
        // differential: the batch complex-method square root must agree
        // with the single-element path on squares, non-squares and the
        // pure-real / pure-imaginary edge cases
        let mut vals: Vec<Fp2> = Vec::new();
        for i in 0..40 {
            let a = rand_fp2(i);
            vals.push(a.square()); // square
            vals.push(a * Fp2::xi()); // ξ = 1+u is a non-square ⇒ non-square
            vals.push(a);
        }
        // pure-real and pure-imaginary entries
        vals.push(Fp2::new(Fp::from_u64(4), Fp::zero()));
        vals.push(Fp2::new(-Fp::from_u64(4), Fp::zero()));
        vals.push(Fp2::new(Fp::zero(), Fp::from_u64(4)));
        vals.push(Fp2::zero());
        let batch = Fp2::sqrt_batch(&vals);
        for (v, b) in vals.iter().zip(batch.iter()) {
            match (v.sqrt(), b) {
                (Some(s), Some(sb)) => {
                    assert_eq!(s, *sb, "batch/single sqrt disagree");
                    assert_eq!(s.square(), *v);
                }
                (None, None) => {}
                (s, b2) => panic!("QR disagreement: single={:?} batch={:?}", s.is_some(), b2.is_some()),
            }
        }
    }

    #[test]
    fn fp2_conj_norm() {
        let a = rand_fp2(0);
        assert_eq!(a * a.conjugate(), Fp2::new(a.norm(), Fp::zero()));
    }
}
