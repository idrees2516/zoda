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

    /// Windowed (w = 5) scalar multiplication by a little-endian limb
    /// scalar — routed through the Jacobian engine (2M+5S doublings,
    /// 7M+4S mixed additions) and converted back to homogeneous form.
    pub fn mul_limbs(&self, scalar: &[u64]) -> Self {
        Jacobian::from_homogeneous(self)
            .mul_limbs(scalar)
            .to_homogeneous()
    }
}

/// Jacobian point `(X : Y : Z)` representing the affine point
/// `(X/Z², Y/Z³)`. The identity is `Z = 0` (X, Y carry no meaning then).
///
/// Jacobian coordinates carry the fastest known complete-ish formulas for
/// short-Weierstrass curves with `a = 0` (the EFD "2007-bl" family):
/// * doubling — 2M + 5S (`dbl-2007-bl`),
/// * mixed addition — 7M + 4S (`madd-2007-bl`),
/// * general addition — 11M + 5S (`add-2007-bl`),
/// against 8M+3S / 11M / 12M for the complete homogeneous RCB formulas.
/// The rare degenerate cases (`H = 0`) are detected explicitly and routed
/// to doubling / the identity, so the formulas are safe for any inputs.
#[derive(Copy, Clone, Debug)]
pub struct Jacobian<C: CurveConfig> {
    pub x: C::Base,
    pub y: C::Base,
    pub z: C::Base,
}

impl<C: CurveConfig> Jacobian<C> {
    /// The identity (Z = 0).
    #[inline]
    pub fn identity() -> Self {
        Jacobian {
            x: C::Base::f_zero(),
            y: C::Base::f_one(),
            z: C::Base::f_zero(),
        }
    }

    #[inline]
    pub fn is_identity(&self) -> bool {
        self.z.f_is_zero()
    }

    /// Lift an affine point (Z = 1).
    #[inline]
    pub fn from_affine(a: &Affine<C>) -> Self {
        if a.infinity {
            Self::identity()
        } else {
            Jacobian { x: a.x, y: a.y, z: C::Base::f_one() }
        }
    }

    /// Project a homogeneous point `(X : Y : Z)` into Jacobian:
    /// `(X·Z : Y·Z² : Z)` — 3M + 1S, no inversion.
    #[inline]
    pub fn from_homogeneous(p: &Projective<C>) -> Self {
        if p.is_identity() {
            return Self::identity();
        }
        let zz = p.z.f_square();
        Jacobian {
            x: p.x.f_mul(p.z),
            y: p.y.f_mul(zz),
            z: p.z,
        }
    }

    /// Project this Jacobian point into homogeneous coordinates:
    /// `(X·Z : Y : Z³)` — 1M + 1S + 1M, no inversion. (Both triples
    /// represent `(X/Z², Y/Z³)`; the weights `(2,3,1)` map to `(1,1,1)`
    /// by multiplying through by `Z`.)
    #[inline]
    pub fn to_homogeneous(&self) -> Projective<C> {
        if self.is_identity() {
            return Projective::identity();
        }
        let zz = self.z.f_square();
        let zzz = zz.f_mul(self.z);
        Projective {
            x: self.x.f_mul(self.z),
            y: self.y,
            z: zzz,
        }
    }

    /// `dbl-2007-bl` — Jacobian doubling for a = 0 (2M + 5S).
    #[inline]
    pub fn double(&self) -> Self {
        if self.is_identity() {
            return *self;
        }
        let a = self.x.f_square(); // A = X1²
        let b = self.y.f_square(); // B = Y1²
        let c = b.f_square(); // C = B²
        let xb = self.x.f_add(b);
        let d0 = xb.f_square().f_sub(a).f_sub(c); // (X1+B)² − A − C
        let d = d0.f_add(d0); // D = 2·d0
        let e = a.f_add(a).f_add(a); // E = 3A
        let f = e.f_square(); // F = E²
        let c8 = {
            let c2 = c.f_add(c);
            let c4 = c2.f_add(c2);
            c4.f_add(c4) // 8C
        };
        let x3 = {
            let d2 = d.f_add(d);
            f.f_sub(d2) // F − 2D
        };
        let y3 = e.f_mul(d.f_sub(x3)).f_sub(c8);
        let z3 = {
            let yz = self.y.f_mul(self.z);
            yz.f_add(yz) // 2·Y1·Z1
        };
        Jacobian { x: x3, y: y3, z: z3 }
    }

    /// `madd-2007-bl` — Jacobian + affine mixed addition (7M + 4S).
    /// Explicitly handles P ± Q edge cases.
    #[inline]
    pub fn add_mixed(&self, rhs: &Affine<C>) -> Self {
        if self.is_identity() {
            return Self::from_affine(rhs);
        }
        if rhs.infinity {
            return *self;
        }
        let z1z1 = self.z.f_square(); // Z1²
        let u2 = rhs.x.f_mul(z1z1); // U2 = x2·Z1²
        let z1z1z1 = z1z1.f_mul(self.z);
        let s2 = rhs.y.f_mul(z1z1z1); // S2 = y2·Z1³
        let h = u2.f_sub(self.x); // H = U2 − X1
        let s2my = s2.f_sub(self.y); // S2 − Y1
        if h.f_is_zero() {
            if s2my.f_is_zero() {
                return self.double(); // P + P
            }
            return Self::identity(); // P + (−P)
        }
        let hh = h.f_square(); // H²
        let hh2 = hh.f_add(hh); // 2H²
        let i = hh2.f_add(hh2); // I = 4H²
        let j = h.f_mul(i); // J = H·I
        let r = s2my.f_add(s2my); // r = 2(S2 − Y1)
        let v = self.x.f_mul(i); // V = X1·I
        let rr = r.f_square();
        let x3 = rr.f_sub(j).f_sub(v).f_sub(v); // r² − J − 2V
        let y3 = {
            let t = v.f_sub(x3);
            let t2 = self.y.f_mul(j); // Y1·J
            let t4 = t2.f_add(t2); // 2·Y1·J
            r.f_mul(t).f_sub(t4)
        };
        let z3 = {
            let zh = self.z.f_add(h);
            zh.f_square().f_sub(z1z1).f_sub(hh) // (Z1+H)² − Z1Z1 − H²
        };
        Jacobian { x: x3, y: y3, z: z3 }
    }

    /// `add-2007-bl` — general Jacobian addition (11M + 5S), with explicit
    /// P ± Q edge-case handling.
    #[inline]
    pub fn add(&self, rhs: &Self) -> Self {
        if self.is_identity() {
            return *rhs;
        }
        if rhs.is_identity() {
            return *self;
        }
        let z1z1 = self.z.f_square();
        let z2z2 = rhs.z.f_square();
        let u1 = self.x.f_mul(z2z2);
        let u2 = rhs.x.f_mul(z1z1);
        let s1 = self.y.f_mul(rhs.z).f_mul(z2z2);
        let s2 = rhs.y.f_mul(self.z).f_mul(z1z1);
        let h = u2.f_sub(u1);
        let s2ms1 = s2.f_sub(s1);
        if h.f_is_zero() {
            if s2ms1.f_is_zero() {
                return self.double();
            }
            return Self::identity();
        }
        let h2 = h.f_add(h); // 2H
        let i = h2.f_square(); // I = (2H)²
        let j = h.f_mul(i);
        let r = s2ms1.f_add(s2ms1); // r = 2(S2 − S1)
        let v = u1.f_mul(i);
        let rr = r.f_square();
        let x3 = rr.f_sub(j).f_sub(v).f_sub(v);
        let y3 = {
            let t = v.f_sub(x3);
            let t2 = s1.f_mul(j);
            let t4 = t2.f_add(t2);
            r.f_mul(t).f_sub(t4)
        };
        let z3 = {
            let zz = self.z.f_add(rhs.z);
            zz.f_square().f_sub(z1z1).f_sub(z2z2).f_mul(h)
        };
        Jacobian { x: x3, y: y3, z: z3 }
    }

    /// Normalize with a single field inversion.
    #[inline]
    pub fn to_affine(&self) -> Affine<C> {
        if self.is_identity() {
            return Affine::identity();
        }
        match self.z.f_invert() {
            Some(zinv) => {
                let z2 = zinv.f_square();
                let z3 = z2.f_mul(zinv);
                Affine {
                    x: self.x.f_mul(z2),
                    y: self.y.f_mul(z3),
                    infinity: false,
                }
            }
            None => Affine::identity(),
        }
    }

    /// Batch-normalize many Jacobian points with a single field inversion
    /// (Montgomery's trick). Identity inputs (Z = 0) map to the affine
    /// identity and never enter the inversion chain.
    pub fn batch_normalize(points: &[Self]) -> Vec<Affine<C>> {
        let n = points.len();
        let mut out = vec![Affine::identity(); n];
        // prefix products of the nonzero Z's (placeholders for zero Z's)
        let mut prods = Vec::with_capacity(n);
        let mut acc = C::Base::f_one();
        for p in points {
            if p.z.f_is_zero() {
                prods.push(C::Base::f_one());
            } else {
                prods.push(acc);
                acc = acc.f_mul(p.z);
            }
        }
        // acc = product of all nonzero Z's; invert once
        let mut acc = match acc.f_invert() {
            Some(inv) => inv,
            None => return points.iter().map(|p| p.to_affine()).collect(),
        };
        for i in (0..n).rev() {
            let p = points[i];
            if p.z.f_is_zero() {
                continue; // out[i] already the identity
            }
            let zinv = acc.f_mul(prods[i]); // Z_i^{-1}
            acc = acc.f_mul(p.z); // strip Z_i from the running suffix inverse
            let z2 = zinv.f_square();
            let z3 = z2.f_mul(zinv);
            out[i] = Affine {
                x: p.x.f_mul(z2),
                y: p.y.f_mul(z3),
                infinity: false,
            };
        }
        out
    }

    /// Windowed (w = 5) scalar multiplication by a little-endian limb
    /// scalar, in Jacobian coordinates: ~256 doublings + ~51 mixed
    /// additions + a single batch-normalized table — roughly 2x faster
    /// than the homogeneous windowed method.
    pub fn mul_limbs(&self, scalar: &[u64]) -> Self {
        if self.is_identity() {
            return *self;
        }
        const W: usize = 5;
        const T: usize = 1 << W; // 32
        // table[i] = i·P for i in 0..32, built purely projectively
        // (even entries by doubling, odd entries by adding P once), then
        // normalized with a single field inversion.
        let mut tab_proj = [Self::identity(); T];
        tab_proj[1] = *self;
        for i in 2..T {
            if i % 2 == 0 {
                tab_proj[i] = tab_proj[i / 2].double();
            } else {
                tab_proj[i] = tab_proj[i - 1].add(self);
            }
        }
        let tab = Self::batch_normalize(&tab_proj);
        let total_bits = scalar.len() * 64;
        // pad the scan to a whole number of W-bit windows (the extractor
        // treats bits beyond the scalar as zero); W=5 does not divide 256,
        // so starting at `total_bits` would misalign the final window.
        let nwin = (total_bits + W - 1) / W;
        let mut acc = Self::identity();
        let mut bit = nwin * W;
        while bit > 0 {
            let start = bit.saturating_sub(W);
            let mut nib = 0usize;
            for b in start..bit {
                let limb = b / 64;
                let off = b % 64;
                if limb < scalar.len() && (scalar[limb] >> off) & 1 == 1 {
                    nib |= 1 << (b - start);
                }
            }
            if !acc.is_identity() {
                for _ in 0..W {
                    acc = acc.double();
                }
            }
            if nib != 0 {
                acc = acc.add_mixed(&tab[nib]);
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

    type TJ = Jacobian<crate::g1::G1Config>;
    type TP = Projective<crate::g1::G1Config>;

    #[test]
    fn jacobian_matches_homogeneous() {
        let g = TP::generator();
        let gaff = crate::g1::G1Affine::generator();
        let gj = TJ::from_affine(&gaff);
        let mut pj = gj;
        let mut ph = g;
        for i in 0..40 {
            // doubling chain
            pj = pj.double();
            ph = ph.double();
            assert_eq!(pj.to_affine(), ph.to_affine(), "double step {}", i);
            assert!(pj.to_affine().is_on_curve());
        }
        // additions: random-ish combos
        let a = gj.mul_limbs(&[7, 0, 0, 0]);
        let b = gj.mul_limbs(&[11, 0, 0, 0]);
        let c = gj.mul_limbs(&[18, 0, 0, 0]);
        assert_eq!(a.add(&b).to_affine(), c.to_affine());
        // mixed additions
        let a_aff = a.to_affine();
        assert_eq!(gj.add_mixed(&a_aff).to_affine(), gj.add(&a).to_affine());
        // edge cases: P + P, P + (-P), identity
        assert_eq!(gj.add(&gj).to_affine(), gj.double().to_affine());
        assert_eq!(gj.add_mixed(&gj.to_affine()).to_affine(), gj.double().to_affine());
        assert!(gj.add(&gj.neg_j()).is_identity());
        assert!(gj.add_mixed(&gj.to_affine().neg()).is_identity());
        // identity absorbs
        assert_eq!(TJ::identity().add(&gj).to_affine(), gj.to_affine());
        assert_eq!(gj.add(&TJ::identity()).to_affine(), gj.to_affine());
        assert_eq!(TJ::identity().add_mixed(&gj.to_affine()).to_affine(), gj.to_affine());
        assert_eq!(TJ::identity().double().to_affine(), TJ::identity().to_affine());
        // cross-coordinate-system consistency
        assert_eq!(gj.to_homogeneous().to_affine(), gj.to_affine());
        assert_eq!(TJ::from_homogeneous(&g).to_affine(), g.to_affine());
    }

    #[test]
    fn jacobian_scalar_mul_matches() {
        // independent naive reference: MSB double-and-add over affine points
        fn naive_mul(p: &crate::g1::G1Affine, scalar: &[u64]) -> crate::g1::G1Affine {
            let mut acc = crate::g1::G1Projective::identity();
            let total = scalar.len() * 64;
            for b in (0..total).rev() {
                acc = acc.double();
                if (scalar[b / 64] >> (b % 64)) & 1 == 1 {
                    acc = acc.add_mixed(p);
                }
            }
            acc.to_affine()
        }
        let g = crate::g1::G1Affine::generator();
        let gj = TJ::from_affine(&g);
        for scalar in [
            [0u64, 0, 0, 0],
            [1, 0, 0, 0],
            [2, 0, 0, 0],
            [3, 0, 0, 0],
            [0xdeadbeef, 0x12345678, 0x9abcdef0, 0x0fedcba9],
            [0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff],
        ] {
            let want = naive_mul(&g, &scalar);
            let got = gj.mul_limbs(&scalar).to_affine();
            assert_eq!(got, want, "scalar {:x?}", scalar);
            assert!(got.is_on_curve() || got.is_identity());
        }
    }

    #[test]
    fn jacobian_batch_normalize() {
        let g = crate::g1::G1Affine::generator();
        let gj = TJ::from_affine(&g);
        let pts: Vec<TJ> = (1..40u64)
            .map(|i| gj.mul_limbs(&[i, 0, 0, 0]))
            .chain([TJ::identity(), gj.mul_limbs(&[12345, 0, 0, 0])])
            .collect();
        let affine = TJ::batch_normalize(&pts);
        for (p, a) in pts.iter().zip(affine.iter()) {
            assert_eq!(&p.to_affine(), a);
        }
    }
}

impl<C: CurveConfig> Jacobian<C> {
    /// Negate (Y → −Y); the identity maps to itself.
    #[inline]
    pub fn neg_j(&self) -> Self {
        Jacobian {
            x: self.x,
            y: self.y.f_neg(),
            z: self.z,
        }
    }
}

/// A signed magnitude that fits in three 64-bit limbs (≤ 192 bits) — the
/// half-width scalar representation produced by GLV lattice decomposition.
#[derive(Clone, Copy, Debug)]
pub struct SignedMag {
    pub neg: bool,
    pub m: [u64; 3],
}

impl SignedMag {
    pub fn zero() -> Self {
        SignedMag { neg: false, m: [0; 3] }
    }
    pub fn is_zero(&self) -> bool {
        self.m == [0u64; 3]
    }
    /// Number of significant bits in the magnitude.
    pub fn bits(&self) -> usize {
        for i in (0..3).rev() {
            if self.m[i] != 0 {
                return i * 64 + (64 - self.m[i].leading_zeros() as usize);
            }
        }
        0
    }
    /// The W-bit window starting at bit `start` (bits above 192 read as 0).
    fn window_at(&self, start: usize, w: usize) -> usize {
        let mut v = 0usize;
        for b in start..start + w {
            if b < 192 && (self.m[b / 64] >> (b % 64)) & 1 == 1 {
                v |= 1 << (b - start);
            }
        }
        v
    }
}

impl<C: CurveConfig> Jacobian<C> {
    /// Build the w = 5 window table `[0·P, 1·P, …, 31·P]` in affine form
    /// with a single field inversion (batch normalization).
    pub fn window_table_32(&self) -> [Affine<C>; 32] {
        const T: usize = 32;
        let mut tab_proj = [Self::identity(); T];
        tab_proj[1] = *self;
        for i in 2..T {
            tab_proj[i] = if i % 2 == 0 {
                tab_proj[i / 2].double()
            } else {
                tab_proj[i - 1].add(self)
            };
        }
        let mut out = [Affine::identity(); 32];
        let norm = Self::batch_normalize(&tab_proj);
        out[..T].copy_from_slice(&norm);
        out
    }
}

impl<C: CurveConfig> Jacobian<C> {
    /// Interleaved w = 5 scan over two PRECOMPUTED affine tables:
    /// `[k0]P + [k1]Q` where `tab_p`/`tab_q` hold `[0..32)` multiples.
    /// Signs of the half-scalars must already be folded into the tables.
    pub fn mul_two_tables(
        tab_p: &[Affine<C>; 32],
        tab_q: &[Affine<C>; 32],
        k0: &SignedMag,
        k1: &SignedMag,
    ) -> Self {
        const W: usize = 5;
        let top = k0.bits().max(k1.bits());
        if top == 0 {
            return Self::identity();
        }
        let nwin = (top + W - 1) / W;
        let mut acc = Self::identity();
        let mut bit = nwin * W;
        while bit > 0 {
            let start = bit.saturating_sub(W);
            if !acc.is_identity() {
                for _ in 0..W {
                    acc = acc.double();
                }
            }
            let nib0 = k0.window_at(start, W);
            let nib1 = k1.window_at(start, W);
            if nib0 != 0 {
                acc = acc.add_mixed(&tab_p[nib0]);
            }
            if nib1 != 0 {
                acc = acc.add_mixed(&tab_q[nib1]);
            }
            bit = start;
        }
        acc
    }
}

impl<C: CurveConfig> Jacobian<C> {
    /// Interleaved windowed two-scalar multiplication `[k0]P + [k1]Q` over
    /// half-width (≤ 192-bit) signed scalars — the GLV execution engine.
    ///
    /// Signs are folded into the table bases (`[−m]P = [m](−P)`); the two
    /// w = 5 window tables are built together and batch-normalized with a
    /// single field inversion; the shared scan halves the number of
    /// point doublings compared to a full-width scalar multiplication.
    pub fn mul_two_scalar(&self, q: &Self, k0: &SignedMag, k1: &SignedMag) -> Self {
        let tab_p = if k0.neg { self.neg_j() } else { *self }.window_table_32();
        let tab_q = if k1.neg { q.neg_j() } else { *q }.window_table_32();
        Self::mul_two_tables(&tab_p, &tab_q, k0, k1)
    }
}
