//! GLV (Galbraith–Lin–Scott) endomorphisms and derived speedups for
//! BLS12-381.
//!
//! * **G1**: the j = 0 automorphism φ(x, y) = (β·x, y) with β a primitive
//!   cube root of unity in Fp acts on the r-torsion as [λ] with
//!   λ² + λ + 1 ≡ 0 (mod r). λ is *derived* at runtime as −p² mod r (p² is
//!   a primitive 6th root of unity mod r because r | Φ₁₂(p)) and matched
//!   against β through the canonical generator — no transcribed tables.
//! * **G2**: the untwist–Frobenius–twist endomorphism ψ(x, y) =
//!   (α·x̄, β₂·ȳ) (p ≡ 3 mod 4, so Frobenius on Fp2 is conjugation) acts as
//!   [λ₂] on G2. (α, β₂, λ₂) are *derived* by scanning the order-6 roots of
//!   unity mod r and verifying that the coordinate ratios transfer across
//!   independent points, plus the curve constraint α³ = β₂² = ξ/ξ̄.
//! * **Lattice decomposition**: each scalar k splits as
//!   k ≡ k₀ + k₁·λ (mod r) with |kᵢ| ≲ 2¹³¹ via a Lagrange–Gauss reduced
//!   basis of {(x, y) : x + λy ≡ 0 (mod r)} and Babai nearest-plane
//!   rounding. `Jacobian::mul_two_scalar` then executes the two half-width
//!   scalars in one interleaved scan — ~1.5–1.6× over a full-width
//!   multiplication.
//! * **Cofactor clearing** on G2 replaces the naive 640-bit [h_eff]·P with
//!   the endomorphism chain [x₀²−x₀+1]ψ(P) + [x₀−1]ψ²(P) (x₀ = −BLS_X),
//!   exactly matching the RFC 9380 h_eff action (verified against the
//!   naive ground truth and the RFC test vectors).
//!
//! ⚠️ **Scope**: `mul_*_public` is *variable-time* and mathematically valid
//! only on the r-torsion (φ = [λ] fails off G1/G2). It is used exclusively
//! for public points already subgroup-verified (SRS elements after
//! `validate_kzg_g1`, hash-to-curve outputs, checked pubkeys/signatures)
//! and never for secret scalars — signing keeps the plain windowed path.

use crate::curve::{Jacobian, SignedMag};
use crate::fp::Fp;
use crate::fp2::Fp2;
use crate::g1::{G1Affine, G1Projective};
use crate::g2::{G2Affine, G2Config, G2Projective};
use crate::pairing::BLS_X;
use std::sync::OnceLock;
use zoda_math::{Fr, PrimeField};

// ===========================================================================
// Small signed bignum (10 limbs = 640 bits, enough for every intermediate:
// lattice dot products ≤ 2^513, cofactor scalars, Babai products)
// ===========================================================================

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
struct Sb {
    neg: bool,
    m: [u64; 10],
}

impl Sb {
    fn zero() -> Self {
        Sb { neg: false, m: [0; 10] }
    }
    fn one() -> Self {
        let mut s = Self::zero();
        s.m[0] = 1;
        s
    }
    fn from_u64(v: u64) -> Self {
        let mut s = Self::zero();
        s.m[0] = v;
        s
    }
    fn from_limbs(src: &[u64]) -> Self {
        let mut s = Self::zero();
        for (i, &l) in src.iter().take(10).enumerate() {
            s.m[i] = l;
        }
        s
    }
    fn is_zero(&self) -> bool {
        self.m == [0u64; 10]
    }
    /// Significant bits of the magnitude.
    fn bits(&self) -> usize {
        for i in (0..10).rev() {
            if self.m[i] != 0 {
                return i * 64 + (64 - self.m[i].leading_zeros() as usize);
            }
        }
        0
    }
    fn bit(&self, i: usize) -> bool {
        i < 640 && (self.m[i / 64] >> (i % 64)) & 1 == 1
    }
    fn set_bit(&mut self, i: usize) {
        self.m[i / 64] |= 1 << (i % 64);
    }
    /// Magnitude comparison.
    fn cmp_mag(&self, o: &Sb) -> core::cmp::Ordering {
        use core::cmp::Ordering::*;
        for i in (0..10).rev() {
            if self.m[i] > o.m[i] {
                return Greater;
            }
            if self.m[i] < o.m[i] {
                return Less;
            }
        }
        Equal
    }
    /// self += 1 on the magnitude (ignores sign; used on non-negative q).
    fn add_mag_one(&mut self) {
        let mut carry = 1u128;
        for l in self.m.iter_mut() {
            let t = *l as u128 + carry;
            *l = t as u64;
            carry = t >> 64;
            if carry == 0 {
                break;
            }
        }
    }
    /// self <<= 1 (magnitude only).
    fn shl1_mag(&mut self) {
        let mut carry = 0u64;
        for l in self.m.iter_mut() {
            let nc = *l >> 63;
            *l = (*l << 1) | carry;
            carry = nc;
        }
    }
    /// Unsigned magnitude addition (signs ignored).
    fn add_mag(a: &Sb, b: &Sb) -> Sb {
        let mut out = Sb::zero();
        let mut carry = 0u128;
        for i in 0..10 {
            let t = a.m[i] as u128 + b.m[i] as u128 + carry;
            out.m[i] = t as u64;
            carry = t >> 64;
        }
        out
    }
    /// Unsigned magnitude subtraction (a ≥ b required).
    fn sub_mag(a: &Sb, b: &Sb) -> Sb {
        let mut out = Sb::zero();
        let mut borrow = 0u128;
        for i in 0..10 {
            let t = (a.m[i] as u128).wrapping_sub(b.m[i] as u128 + borrow);
            out.m[i] = t as u64;
            borrow = (t >> 64) & 1;
        }
        debug_assert_eq!(borrow, 0);
        out
    }
    /// Signed addition.
    fn add(&self, o: &Sb) -> Sb {
        if self.neg == o.neg {
            let mut r = Sb::add_mag(self, o);
            r.neg = self.neg;
            r
        } else {
            match self.cmp_mag(o) {
                core::cmp::Ordering::Equal => Sb::zero(),
                core::cmp::Ordering::Greater => {
                    let mut r = Sb::sub_mag(self, o);
                    r.neg = self.neg;
                    r
                }
                core::cmp::Ordering::Less => {
                    let mut r = Sb::sub_mag(o, self);
                    r.neg = o.neg;
                    r
                }
            }
        }
    }
    fn neg(&self) -> Sb {
        if self.is_zero() {
            *self
        } else {
            let mut r = *self;
            r.neg = !r.neg;
            r
        }
    }
    /// Signed subtraction self − o.
    fn sub(&self, o: &Sb) -> Sb {
        self.add(&o.neg())
    }
    /// Signed schoolbook multiplication.
    fn mul(&self, o: &Sb) -> Sb {
        let mut out = Sb::zero();
        for i in 0..10 {
            if self.m[i] == 0 {
                continue;
            }
            let mut carry = 0u128;
            for j in 0..10 - i {
                let t = out.m[i + j] as u128 + (self.m[i] as u128) * (o.m[j] as u128) + carry;
                out.m[i + j] = t as u64;
                carry = t >> 64;
            }
            // carry beyond the top is impossible: products of our bounded
            // inputs (< 2^520) always fit in 640 bits
            debug_assert_eq!(carry, 0);
        }
        out.neg = self.neg != o.neg && !out.is_zero();
        out
    }
}

/// Floor division of magnitudes: (|num| / |den|, |num| mod |den|).
fn div_rem_mag(num: &Sb, den: &Sb) -> (Sb, Sb) {
    debug_assert!(!den.is_zero());
    let mut q = Sb::zero();
    let mut r = Sb::zero();
    let nbits = num.bits();
    let mut i = nbits;
    while i > 0 {
        i -= 1;
        r.shl1_mag();
        if num.bit(i) {
            r.m[0] |= 1;
        }
        if r.cmp_mag(den) != core::cmp::Ordering::Less {
            r = Sb::sub_mag(&r, den);
            q.set_bit(i);
        }
    }
    (q, r)
}

/// round(num / den) for den > 0 (ties away from zero).
fn round_div(num: &Sb, den: &Sb) -> Sb {
    debug_assert!(!den.neg);
    let (mut q, mut r) = div_rem_mag(num, den);
    r.shl1_mag();
    if r.cmp_mag(den) != core::cmp::Ordering::Less {
        q.add_mag_one();
    }
    q.neg = num.neg && !q.is_zero();
    q
}

// ===========================================================================
// Lattice machinery
// ===========================================================================

type Vec2 = (Sb, Sb);

fn dot(a: &Vec2, b: &Vec2) -> Sb {
    a.0.mul(&b.0).add(&a.1.mul(&b.1))
}

/// Lagrange–Gauss reduction of a 2D lattice basis (v1, v2) to a nearly
/// orthogonal basis with short vectors (norm ≈ √det·O(1)).
fn gauss_reduce(mut v1: Vec2, mut v2: Vec2) -> (Vec2, Vec2) {
    let mut guard = 0;
    loop {
        let n1 = dot(&v1, &v1);
        let n2 = dot(&v2, &v2);
        if n1.cmp_mag(&n2) == core::cmp::Ordering::Greater {
            core::mem::swap(&mut v1, &mut v2);
        }
        let num = dot(&v1, &v2);
        let den = dot(&v1, &v1);
        let mu = round_div(&num, &den);
        if mu.is_zero() {
            break;
        }
        v2.0 = v2.0.sub(&mu.mul(&v1.0));
        v2.1 = v2.1.sub(&mu.mul(&v1.1));
        guard += 1;
        if guard > 200 {
            break; // defensive; converges in ~O(log) steps in practice
        }
    }
    (v1, v2)
}

/// One GLV instance: an eigenvalue λ (acting as [λ] via a cheap coordinate
/// map) plus a reduced decomposition lattice.
struct Glv {
    lambda: Fr,
    /// |v_i|² of the reduced basis (Babai denominators).
    n1: Sb,
    n2: Sb,
    basis: (Vec2, Vec2),
}

impl Glv {
    fn new(lambda: Fr) -> Glv {
        let lam = Sb::from_limbs(&lambda.to_repr());
        let r = Sb::from_limbs(&Fr::MODULUS);
        // L = {(x, y) : x + λy ≡ 0 (mod r)} — start from
        // (−λ mod r, 1) and (r, 0), then reduce.
        let v1 = (r.sub(&lam), Sb::one());
        let v2 = (r, Sb::zero());
        let (b1, b2) = gauss_reduce(v1, v2);
        Glv {
            lambda,
            n1: dot(&b1, &b1),
            n2: dot(&b2, &b2),
            basis: (b1, b2),
        }
    }

    /// Babai nearest-plane decomposition of k ≡ k₀ + k₁·λ (mod r).
    fn decompose(&self, k: &Fr) -> (SignedMag, SignedMag) {
        let mut tx = Sb::from_limbs(&k.to_repr());
        let mut ty = Sb::zero();
        // subtract round(<t, v_i>/<v_i, v_i>)·v_i for i = 1, 2
        for (v, n) in [(&self.basis.0, &self.n1), (&self.basis.1, &self.n2)] {
            let num = tx.mul(&v.0).add(&ty.mul(&v.1));
            let q = round_div(&num, n);
            tx = tx.sub(&q.mul(&v.0));
            ty = ty.sub(&q.mul(&v.1));
        }
        debug_assert!(tx.bits() <= 192 && ty.bits() <= 192);
        (
            SignedMag { neg: tx.neg, m: [tx.m[0], tx.m[1], tx.m[2]] },
            SignedMag { neg: ty.neg, m: [ty.m[0], ty.m[1], ty.m[2]] },
        )
    }
}

// ===========================================================================
// G1: β-endomorphism
// ===========================================================================

struct G1Glv {
    glv: Glv,
    beta: Fp,
}


/// p mod r computed exactly with the bignum (the generic byte-based
/// field entry points truncate to 32 bytes, which p exceeds).
fn p_mod_r() -> Fr {
    let (_, rem) = div_rem_mag(
        &Sb::from_limbs(&Fp::MODULUS),
        &Sb::from_limbs(&Fr::MODULUS),
    );
    let mut l = [0u64; 4];
    l.copy_from_slice(&rem.m[..4]);
    Fr::from_repr_limbs(l)
}

fn g1_glv() -> &'static G1Glv {
    static G: OnceLock<G1Glv> = OnceLock::new();
    G.get_or_init(|| {
        // ---- β: primitive cube root of unity in Fp ----
        // (p − 1)/3 via bignum division
        let p_minus_1 = Sb::from_limbs(&Fp::MODULUS).sub(&Sb::one());
        let (exp, _) = div_rem_mag(&p_minus_1, &Sb::from_u64(3));
        let mut exp_limbs = [0u64; 10];
        exp_limbs.copy_from_slice(&exp.m);
        let mut beta = Fp::one();
        for g in [2u64, 3, 5, 7, 11, 13, 17] {
            let cand = Fp::from_u64(g).pow_limbs(&exp_limbs);
            if cand != Fp::one() && cand * cand * cand == Fp::one() {
                beta = cand;
                break;
            }
        }
        assert!(beta != Fp::one(), "no cube root of unity found");

        // ---- λ: root of λ² + λ + 1 ≡ 0 (mod r), matched to β ----
        // p² is a primitive 6th root of unity mod r (r | Φ₁₂(p)), so
        // −p² and −p⁴ = 1 − p² are the two primitive cube roots.
        let pmr = p_mod_r();
        let u = pmr * pmr; // p² mod r
        let cands = [-u, u - Fr::ONE];
        let gen = G1Affine::generator();
        let mut found: Option<(Fr, Fp)> = None;
        for l in cands {
            assert_eq!(l * l + l + Fr::ONE, Fr::zero(), "λ must satisfy λ²+λ+1=0");
            let r = G1Projective::generator()
                .mul_limbs(&l.to_repr())
                .to_affine();
            for b in [beta, beta * beta] {
                if r.x == b * gen.x && r.y == gen.y {
                    found = Some((l, b));
                }
            }
        }
        let (lambda, beta) = found.expect("no (λ, β) pair matches the G1 generator");

        G1Glv {
            glv: Glv::new(lambda),
            beta,
        }
    })
}

/// φ(x, y) = (β·x, y) — the j = 0 automorphism, one Fp multiplication.
pub fn g1_phi(p: &G1Projective) -> G1Projective {
    let beta = g1_glv().beta;
    G1Projective { x: beta * p.x, y: p.y, z: p.z }
}


/// GLV-accelerated scalar multiplication on a **public, subgroup-verified**
/// G1 point: `[k]P = [k₀]P + [k₁]φ(P)` over half-width signed scalars.
/// The φ-multiples table is *derived* from the P-table by the coordinate
/// map (φ([i]P) = [i]φ(P)), so only one table of curve additions is built.
pub fn mul_g1_public(p: &G1Projective, k: &Fr) -> G1Projective {
    let g = g1_glv();
    let (k0, k1) = g.glv.decompose(k);
    let base = if k0.neg { Jacobian::from_homogeneous(p).neg_j() } else { Jacobian::from_homogeneous(p) };
    let tab_p = base.window_table_32();
    // φ-table: (β·x, y) per entry, sign of k1 folded by negating y
    let beta = g.beta;
    // tab_p carries k0's sign fold; φ(tab_p) inherits it, but the second
    // scalar needs k1's fold — negate exactly when the two differ.
    let flip = k0.neg != k1.neg;
    let mut tab_phi = [crate::g1::G1Affine::identity(); 32];
    for i in 0..32 {
        tab_phi[i] = crate::g1::G1Affine {
            x: beta * tab_p[i].x,
            y: if flip { -tab_p[i].y } else { tab_p[i].y },
            infinity: tab_p[i].infinity,
        };
    }
    Jacobian::mul_two_tables(&tab_p, &tab_phi, &k0, &k1).to_homogeneous()
}

/// λ such that φ acts as [λ] on G1.
pub fn g1_lambda() -> Fr {
    g1_glv().glv.lambda
}

// ===========================================================================
// G2: untwist–Frobenius–twist endomorphism ψ
// ===========================================================================

struct G2Glv {
    glv: Glv,
    psi_alpha: Fp2,
    psi_beta: Fp2,
}

fn g2_glv() -> &'static G2Glv {
    static G: OnceLock<G2Glv> = OnceLock::new();
    G.get_or_init(|| {
        // λ₂ candidates: the odd powers of p mod r (p has order 12;
        // the eigenvalue of the untwist–Frobenius–twist map is p itself).
        let pmr = p_mod_r();
        let pmr3 = pmr * pmr * pmr;
        let cands = [pmr, -pmr, pmr3, -pmr3];

        // two independent G2 points for the ratio-transfer test
        let q0 = G2Projective::generator()
            .mul_limbs(&[0x9e3779b97f4a7c15, 0x1234567890abcdef])
            .to_affine();
        let q1 = G2Projective::generator()
            .mul_limbs(&[0xf00dcafeba5ebeef, 0x0123456789abcdef])
            .to_affine();

        // curve constraint: α³ = β² = δ = ξ/ξ̄
        let xi = Fp2::xi();
        let delta = xi * xi.conjugate().invert().expect("ξ invertible");

        let mut found: Option<(Fr, Fp2, Fp2)> = None;
        for l in cands {
            let r0 = q0.to_projective().mul_limbs(&l.to_repr()).to_affine();
            let a0 = r0.x * q0.x.conjugate().invert().unwrap();
            let b0 = r0.y * q0.y.conjugate().invert().unwrap();
            // the ratios must transfer to an independent point
            let r1 = q1.to_projective().mul_limbs(&l.to_repr()).to_affine();
            if r1.x != a0 * q1.x.conjugate() || r1.y != b0 * q1.y.conjugate() {
                continue;
            }
            // and satisfy the curve constraint
            if !(a0 * a0 * a0 == delta && b0 * b0 == delta) {
                continue;
            }
            found = Some((l, a0, b0));
            break;
        }
        let (lambda2, alpha, beta2) = found.expect("no (λ₂, α, β₂) matches the G2 endomorphism");
        // λ₂ is a 12th root of unity: λ⁴ − λ² + 1 = 0 (the minimal
        // polynomial of p mod r).
        let l2 = lambda2 * lambda2;
        assert!(
            l2 * l2 - l2 + Fr::ONE == Fr::zero(),
            "λ₂ must satisfy λ⁴ − λ² + 1 = 0"
        );

        G2Glv {
            glv: Glv::new(lambda2),
            psi_alpha: alpha,
            psi_beta: beta2,
        }
    })
}

/// ψ(x, y) = (α·x̄, β₂·ȳ) — the untwist–Frobenius–twist endomorphism of the
/// degree-6 twist (Frobenius on Fp2 is conjugation since p ≡ 3 mod 4).
pub fn g2_psi(p: &G2Projective) -> G2Projective {
    let g = g2_glv();
    G2Projective {
        x: g.psi_alpha * p.x.conjugate(),
        y: g.psi_beta * p.y.conjugate(),
        z: p.z.conjugate(),
    }
}

/// GLV-accelerated scalar multiplication on a **public, subgroup-verified**
/// G2 point. The ψ-multiples table is derived from the P-table by the
/// coordinate map (ψ([i]P) = [i]ψ(P)): two Fp2 multiplies per entry
/// instead of a full table of Fp2 curve additions.
pub fn mul_g2_public(p: &G2Projective, k: &Fr) -> G2Projective {
    let g = g2_glv();
    let (k0, k1) = g.glv.decompose(k);
    let base = if k0.neg {
        Jacobian::<G2Config>::from_homogeneous(p).neg_j()
    } else {
        Jacobian::<G2Config>::from_homogeneous(p)
    };
    let tab_p = base.window_table_32();
    let (alpha, beta2) = (g.psi_alpha, g.psi_beta);
    let flip = k0.neg != k1.neg;
    let mut tab_psi = [crate::g2::G2Affine::identity(); 32];
    for i in 0..32 {
        let psi_x = alpha * tab_p[i].x.conjugate();
        let psi_y = beta2 * tab_p[i].y.conjugate();
        tab_psi[i] = crate::g2::G2Affine {
            x: psi_x,
            y: if flip { -psi_y } else { psi_y },
            infinity: tab_p[i].infinity,
        };
    }
    Jacobian::mul_two_tables(&tab_p, &tab_psi, &k0, &k1).to_homogeneous()
}

/// λ₂ such that ψ acts as [λ₂] on G2.
pub fn g2_lambda() -> Fr {
    g2_glv().glv.lambda
}

// ===========================================================================
// Cofactor clearing via the ψ-chain
// ===========================================================================

/// Clear the G2 cofactor exactly as `h_eff · P` (RFC 9380 §8.8.2) using
/// the Budroni–Pintore endomorphism chain (RFC 9380 Appendix G.3):
///
/// ```text
/// c1 = −BLS_X (the signed BLS parameter)
/// t1 = c1·P          t2 = ψ(P)
/// t3 = ψ2(2P) − t2   t2 = c1·(t1 + t2)
/// t3 = t3 + t2 − t1 − P
/// ```
///
/// Two 64-bit scalar multiplications, one doubling, a handful of
/// additions, and the cheap ψ/ψ2 coordinate maps — instead of the naive
/// 640-bit `h_eff` multiplication (~4× faster).
pub fn clear_cofactor_g2(p: &G2Projective) -> G2Projective {
    let c1_mul = |pt: &G2Projective| -> G2Projective { pt.mul_limbs(&[BLS_X]).neg() };
    let t1 = c1_mul(p); // t1 = c1·P
    let t2 = g2_psi(p); // t2 = ψ(P)
    let mut t3 = g2_psi2(&p.double()); // ψ2(2P)
    t3 = t3.add(&t2.neg()); // t3 = ψ2(2P) − t2
    let t2 = c1_mul(&t1.add(&t2)); // t2 = c1·(t1 + t2)
    t3 = t3.add(&t2).add(&t1.neg()).add(&p.neg());
    t3
}

/// ψ2 = ψ∘ψ (RFC 9380 Appendix G.3): `(x, y) ↦ (c₁₂·x, −y)` with
/// `c₁₂ = 2^−((p−1)/3)` — two Fp2 multiplications by a constant... one
/// constant multiply and a negation, exacting the double ψ for free.
pub fn g2_psi2(p: &G2Projective) -> G2Projective {
    let c = psi2_const();
    G2Projective {
        x: c * p.x,
        y: -p.y,
        z: p.z,
    }
}

fn psi2_const() -> Fp2 {
    static C: OnceLock<Fp2> = OnceLock::new();
    *C.get_or_init(|| {
        // exponent (p − 1)/3
        let p_minus_1 = Sb::from_limbs(&Fp::MODULUS).sub(&Sb::one());
        let (exp, _) = div_rem_mag(&p_minus_1, &Sb::from_u64(3));
        let mut e = [0u64; 10];
        e.copy_from_slice(&exp.m);
        let two = Fp2::new(Fp::from_u64(2), Fp::zero());
        let c = two.pow_limbs(&e).invert().expect("2^((p-1)/3) invertible");
        // c³ = 2^−(p−1) = 1 by Fermat
        assert_eq!(c * c * c, Fp2::one(), "psi2 constant self-check");
        c
    })
}

// ===========================================================================
// Subgroup checks
// ===========================================================================

/// r-torsion membership [r]P == O (plain path — valid on any curve point).
pub fn subgroup_check_g1(p: &G1Affine) -> bool {
    if p.infinity {
        return true;
    }
    p.to_projective().mul_limbs(&Fr::MODULUS).is_identity()
}

/// r-torsion membership for G2 (plain path).
pub fn subgroup_check_g2(p: &G2Affine) -> bool {
    if p.infinity {
        return true;
    }
    p.to_projective().mul_limbs(&Fr::MODULUS).is_identity()
}

/// The G1 cofactor `h₁ = #E(Fp)/r = (p + |x₀|)/r` as little-endian limbs
/// (~2^126, three limbs). #E(Fp) = p + 1 − t with t = x₀ + 1 and
/// x₀ = −BLS_X, so #E = p + BLS_X. Verified exactly at first use.
pub fn g1_cofactor_h1() -> [u64; 3] {
    static H1: OnceLock<[u64; 3]> = OnceLock::new();
    *H1.get_or_init(|| {
        // n = p + BLS_X
        let n = Sb::from_limbs(&Fp::MODULUS).add(&Sb::from_u64(BLS_X));
        let (h, rem) = div_rem_mag(&n, &Sb::from_limbs(&Fr::MODULUS));
        assert!(rem.is_zero(), "#E(Fp) must be divisible by r");
        assert!(h.bits() <= 192, "h1 unexpectedly large");
        // self-check: h₁ · r == p + BLS_X
        assert!(h.mul(&Sb::from_limbs(&Fr::MODULUS)) == n, "h1 * r != #E");
        [h.m[0], h.m[1], h.m[2]]
    })
}

/// Batch r-torsion membership for G2. Random-combination batching is
/// **unsound** on E'(Fp2): a small-order point `P` (e.g. of order 3 —
/// small factors divide the cofactor h₂ ≈ 2^507) vanishes from the
/// combination whenever `ord(P) | c`, which an adversary can grind in a
/// handful of transcript attempts. The `c ≡ 1 (mod h)` trick that fixes
/// this for G1 (see `zoda_kzg::eip4844::batch_subgroup_check_g1`) would
/// need ~635-bit scalars here — costlier than the direct check. So this
/// performs the plain per-point `[r]P == O` verification (the amortised
/// pairing cost of a batch still dominates overall).
pub fn batch_subgroup_check_g2(points: &[G2Affine]) -> bool {
    points.iter().all(|p| subgroup_check_g2(p))
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::CurveConfig;
    use zoda_math::ZodaRng;

    fn random_g1(s: u64) -> G1Projective {
        G1Projective::generator().mul_limbs(&[s, s.wrapping_mul(0x9e3779b97f4a7c15), 0, 0])
    }

    fn random_g2(s: u64) -> G2Projective {
        G2Projective::generator().mul_limbs(&[s, s.wrapping_mul(0x51ed270b), 0, 0])
    }

    /// A uniformly random point of E'(Fp2) (almost surely NOT in G2).
    fn random_curve_point_g2(rng: &mut ZodaRng) -> G2Projective {
        loop {
            let mut b0 = [0u8; 32];
            let mut b1 = [0u8; 32];
            rng.next_bytes(&mut b0);
            rng.next_bytes(&mut b1);
            let x = Fp2::new(Fp::from_le_bytes32(&b0), Fp::from_le_bytes32(&b1));
            let y2 = x.square() * x + G2Config::b();
            if let Some(y) = y2.sqrt() {
                return G2Projective { x, y, z: Fp2::one() };
            }
        }
    }

    #[test]
    fn g1_endomorphism_action() {
        // [λ]P == φ(P) for random P ∈ G1
        for i in 1..6u64 {
            let p = random_g1(i * 0x1234 + 7).to_affine();
            let lam = g1_glv().glv.lambda;
            // scalar-mul P by λ directly (P is in G1)
            let scaled = p.to_projective().mul_limbs(&lam.to_repr()).to_affine();
            let phi = g1_phi(&p.to_projective()).to_affine();
            assert_eq!(scaled, phi, "G1 endomorphism mismatch at i={}", i);
        }
    }

    #[test]
    fn g2_endomorphism_action() {
        for i in 1..4u64 {
            let p = random_g2(i * 991 + 3).to_affine();
            let lam = g2_glv().glv.lambda;
            let scaled = p.to_projective().mul_limbs(&lam.to_repr()).to_affine();
            let psi = g2_psi(&p.to_projective()).to_affine();
            assert_eq!(scaled, psi, "G2 endomorphism mismatch at i={}", i);
        }
    }

    #[test]
    fn glv_decompose_bounds() {
        let mut rng = ZodaRng::from_seed(*b"glv-decompose-000000000000000000");
        for _ in 0..50 {
            let k = rng.next_fr(true);
            let (k0, k1) = g1_glv().glv.decompose(&k);
            assert!(k0.bits() <= 134, "k0 too large: {}", k0.bits());
            assert!(k1.bits() <= 134, "k1 too large: {}", k1.bits());
            // congruence: [k]g == [k0]g + [k1]φ(g) (checked via mul_g1_public)
        }
    }

    #[test]
    fn glv_mul_matches_naive_g1() {
        let mut rng = ZodaRng::from_seed(*b"glv-g1-naive-0000000000000000000");
        let p = random_g1(0x5eed).to_affine();
        let mut scalars: Vec<Fr> = vec![Fr::zero(), Fr::ONE, -Fr::ONE, Fr::from_u64(2)];
        for _ in 0..12 {
            scalars.push(rng.next_fr(true));
        }
        for k in &scalars {
            let want = p.to_projective().mul_limbs(&k.to_repr()).to_affine();
            let got = mul_g1_public(&p.to_projective(), k).to_affine();
            assert_eq!(got, want, "G1 GLV mismatch for {:?}", k);
        }
    }

    #[test]
    fn glv_mul_matches_naive_g2() {
        let mut rng = ZodaRng::from_seed(*b"glv-g2-naive-0000000000000000000");
        let p = random_g2(0x600d).to_affine();
        let mut scalars: Vec<Fr> = vec![Fr::zero(), Fr::ONE, -Fr::ONE];
        for _ in 0..6 {
            scalars.push(rng.next_fr(true));
        }
        for k in &scalars {
            let want = p.to_projective().mul_limbs(&k.to_repr()).to_affine();
            let got = mul_g2_public(&p.to_projective(), k).to_affine();
            assert_eq!(got, want, "G2 GLV mismatch for {:?}", k);
        }
    }

    #[test]
    fn psi2_matches_double_psi() {
        let mut rng = ZodaRng::from_seed(*b"psi2-check-000000000000000000000");
        for _ in 0..3 {
            let p = random_curve_point_g2(&mut rng);
            let once = g2_psi(&p);
            let twice = g2_psi(&once);
            assert_eq!(g2_psi2(&p).to_affine(), twice.to_affine());
        }
        // also on a G2 point
        let q = random_g2(0xabc);
        let once = g2_psi(&q);
        assert_eq!(g2_psi2(&q).to_affine(), g2_psi(&once).to_affine());
    }

    #[test]
    fn clear_cofactor_matches_naive_h_eff() {
        // Ground truth: the naive [G2_H_EFF]·P from the hash module.
        let naive = |p: &G2Projective| -> G2Projective {
            p.mul_limbs(&crate::hash::G2_H_EFF)
        };
        let mut rng = ZodaRng::from_seed(*b"cofactor-chain-00000000000000000");
        // random full-curve points (not in G2)
        for _ in 0..3 {
            let p = random_curve_point_g2(&mut rng);
            assert_eq!(clear_cofactor_g2(&p).to_affine(), naive(&p).to_affine());
        }
        // points already in G2: clearing must act as the identity there
        for i in 1..3u64 {
            let p = random_g2(i * 77 + 5);
            assert_eq!(clear_cofactor_g2(&p).to_affine(), naive(&p).to_affine());
        }
    }

    #[test]
    fn batch_subgroup_check_g2_accepts_and_rejects() {
        let good: Vec<G2Affine> = (1..=4u64)
            .map(|i| random_g2(i * 31 + 1).to_affine())
            .collect();
        assert!(batch_subgroup_check_g2(&good));
        // a random curve point is almost surely not in G2
        let mut rng = ZodaRng::from_seed(*b"batch-subgroup-00000000000000000");
        let mut bad = good.clone();
        bad.push(random_curve_point_g2(&mut rng).to_affine());
        assert!(!batch_subgroup_check_g2(&bad));
        // empty is trivially fine
        assert!(batch_subgroup_check_g2(&[]));
    }

    #[test]
    fn signed_bignum_axioms() {
        let mut rng = ZodaRng::from_seed(*b"sbignum-axiom-000000000000000000");
        let rand_sb = |rng: &mut ZodaRng| -> Sb {
            let mut bytes = [0u8; 80];
            rng.next_bytes(&mut bytes);
            let mut s = Sb::zero();
            for i in 0..10 {
                let mut w = [0u8; 8];
                w.copy_from_slice(&bytes[i * 8..i * 8 + 8]);
                s.m[i] = u64::from_le_bytes(w);
            }
            // keep values < 2^320 so products stay in range
            for i in 5..10 {
                s.m[i] = 0;
            }
            s.neg = rng.next_u64() & 1 == 1;
            s
        };
        for _ in 0..50 {
            let a = rand_sb(&mut rng);
            let b = rand_sb(&mut rng);
            let c = rand_sb(&mut rng);
            // add/sub associativity
            assert_eq!(a.add(&b).add(&c), a.add(&b.add(&c)));
            assert_eq!(a.sub(&b).add(&b), a);
            // distributivity of mul over add
            assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
            // round_div on exact multiples
            let d = a.mul(&b);
            if !b.is_zero() && !d.is_zero() {
                let q = round_div(&d, &b);
                assert!(q.is_zero() || q.bits() <= a.bits() + 1);
            }
        }
        // div_rem_mag exactness
        let a = Sb::from_u64(12345678901234567890);
        let b = Sb::from_u64(987654321);
        let (q, r) = div_rem_mag(&a, &b);
        assert_eq!(q.mul(&b).add(&r), a);
        // round_div ties and signs
        let two = Sb::from_u64(2);
        let five = Sb::from_u64(5);
        let three = Sb::from_u64(3);
        assert_eq!(round_div(&five, &two).m[0], 3);
        let neg_five = Sb { neg: true, ..five };
        assert!(round_div(&neg_five, &two).neg);
    }
}
