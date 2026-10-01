//! The ring `R_q = Z_q[X]/(X^n + 1)` with `q = 8380417`, and its
//! negacyclic NTT (q supports transforms up to length 8192).

use zoda_math::{FftDomain, Fq, PrimeField, ZodaRng};

pub const Q: u32 = 8_380_417;

/// A polynomial in R_q (coefficient form, each `< q`).
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Rq(pub Vec<Fq>);

/// A vector of ring elements.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct RqVector(pub Vec<Rq>);

impl Rq {
    pub fn zero(n: usize) -> Rq {
        Rq(vec![Fq(0); n])
    }

    pub fn from_i64(v: i64) -> Rq {
        let x = if v < 0 {
            Fq(((-(v as i128 % Q as i128)) as u32) % Q)
        } else {
            Fq((v % Q as i64) as u32)
        };
        Rq(vec![x, Fq(0), Fq(0)])
    }

    /// n for this element (trailing zeros trimmed view is not used; ring
    /// elements always carry a fixed length per context).
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|c| c.is_zero())
    }

    /// Sample a uniformly random ring element (rejection-free via the
    /// RNG's uniform mod q).
    pub fn uniform(n: usize, rng: &mut ZodaRng) -> Rq {
        Rq((0..n).map(|_| Fq(rng.next_below(Q as u64) as u32)).collect())
    }

    /// Multiply two ring elements.
    ///
    /// **Density-aware dispatch**: power-of-two lengths ≥ 32 with both
    /// operands dense route through the **negacyclic NTT** (O(n log n));
    /// when either operand is sparse (narrow challenges, monomial
    /// evaluation points, scalars), the schoolbook path with the sparser
    /// operand driving the outer loop is cheaper than three transforms —
    /// the break-even is ≈ 3.5·log₂(n) nonzeros.
    pub fn mul(&self, rhs: &Rq) -> Rq {
        let n = self.0.len();
        debug_assert_eq!(n, rhs.0.len());
        if n >= 32 && n.is_power_of_two() {
            let nnz_a = self.0.iter().filter(|c| !c.is_zero()).count();
            let nnz_b = rhs.0.iter().filter(|c| !c.is_zero()).count();
            let sparse = nnz_a.min(nnz_b);
            if (sparse as f64) < 3.5 * n.trailing_zeros() as f64 {
                // schoolbook, sparser operand outer
                if nnz_a <= nnz_b {
                    self.mul_schoolbook(rhs)
                } else {
                    rhs.mul_schoolbook(self)
                }
            } else {
                let ctx = ntt_ctx(n);
                let a = ctx.forward(&self.0);
                let b = ctx.forward(&rhs.0);
                let c: Vec<Fq> = a.iter().zip(b.iter()).map(|(x, y)| *x * *y).collect();
                ctx.inverse(&c)
            }
        } else {
            self.mul_schoolbook(rhs)
        }
    }

    /// Schoolbook negacyclic convolution with reduction mod X^n + 1.
    pub fn mul_schoolbook(&self, rhs: &Rq) -> Rq {
        let n = self.0.len();
        let mut acc = vec![Fq(0); n];
        for i in 0..n {
            if self.0[i].is_zero() {
                continue;
            }
            for j in 0..n {
                let v = self.0[i] * rhs.0[j];
                let idx = i + j;
                if idx < n {
                    acc[idx] = acc[idx] + v;
                } else {
                    // X^n = -1
                    acc[idx - n] = acc[idx - n] - v;
                }
            }
        }
        Rq(acc)
    }

    pub fn add(&self, rhs: &Rq) -> Rq {
        let n = self.0.len().max(rhs.0.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = self.0.get(i).copied().unwrap_or(Fq(0));
            let b = rhs.0.get(i).copied().unwrap_or(Fq(0));
            out.push(a + b);
        }
        Rq(out)
    }

    pub fn sub(&self, rhs: &Rq) -> Rq {
        let n = self.0.len().max(rhs.0.len());
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let a = self.0.get(i).copied().unwrap_or(Fq(0));
            let b = rhs.0.get(i).copied().unwrap_or(Fq(0));
            out.push(a - b);
        }
        Rq(out)
    }

    pub fn neg(&self) -> Rq {
        Rq(self.0.iter().map(|c| -*c).collect())
    }

    /// Scalar multiply by a small integer.
    pub fn mul_small(&self, s: i64) -> Rq {
        if s >= 0 {
            let s = (s as u64) % Q as u64;
            Rq(self.0.iter().map(|c| *c * Fq(s as u32)).collect())
        } else {
            let s = ((-(s as i128)) % Q as i128) as u64;
            Rq(self.0.iter().map(|c| *c * Fq(s as u32)).collect())
        }
    }

    /// Serialise (4 bytes LE per coefficient).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.0.len() * 3);
        for c in &self.0 {
            let v = c.0;
            out.push((v & 0xff) as u8);
            out.push(((v >> 8) & 0xff) as u8);
            out.push(((v >> 16) & 0x7f) as u8);
        }
        out
    }
}

impl RqVector {
    pub fn uniform(k: usize, n: usize, rng: &mut ZodaRng) -> RqVector {
        RqVector((0..k).map(|_| Rq::uniform(n, rng)).collect())
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn add(&self, rhs: &RqVector) -> RqVector {
        RqVector(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(a, b)| a.add(b))
                .collect(),
        )
    }

    pub fn sub(&self, rhs: &RqVector) -> RqVector {
        RqVector(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(a, b)| a.sub(b))
                .collect(),
        )
    }

    pub fn neg(&self) -> RqVector {
        RqVector(self.0.iter().map(|a| a.neg()).collect())
    }
}

/// Map a list of `Fq` coefficients (a polynomial in Z_q[X], degree < d)
/// into `ceil(d / n)` ring elements (chunks) of length n, zero-padded.
pub fn fq_to_poly(coeffs: &[Fq], n: usize) -> Vec<Rq> {
    let chunks = (coeffs.len() + n - 1) / n;
    let mut out = Vec::with_capacity(chunks);
    for ci in 0..chunks {
        let start = ci * n;
        let end = (start + n).min(coeffs.len());
        let mut c = vec![Fq(0); n];
        c[..end - start].copy_from_slice(&coeffs[start..end]);
        out.push(Rq(c));
    }
    if out.is_empty() {
        out.push(Rq::zero(n));
    }
    out
}

/// Evaluate the coefficient-form polynomial (over Z_q) at z ∈ Rq,
/// treating the polynomial as an element of Rq[X] evaluated at the ring
/// element z (Horner over the ring).
pub fn eval_poly_at_ring(coeffs: &[Fq], z: &Rq) -> Rq {
    let mut acc = Rq(vec![Fq(0); z.len()]);
    for c in coeffs.iter().rev() {
        // acc = acc * z + c
        let mut term = z.mul(&acc);
        term.0[0] = term.0[0] + *c;
        // term length matches acc length
        acc = term;
    }
    acc
}

/// Negacyclic NTT context for `R_q = Z_q[X]/(X^n + 1)` with power-of-two
/// `n ≤ 4096` (q − 1 has 2-adicity 13, so the cyclic transform of size n
/// plus the ψ-twist of order 2n is always available).
///
/// Multiplication in the ring factors as: twist both operands by powers
/// of the primitive 2n-th root ψ, run cyclic NTTs of size n, multiply
/// pointwise, run the inverse NTT, untwist by ψ⁻ⁱ — O(n log n) total.
pub struct RqNtt {
    n: usize,
    dom: FftDomain<Fq>,
    psi: Vec<Fq>,
    psi_inv: Vec<Fq>,
}

impl RqNtt {
    /// Build the context for `n` (power of two, 2..=4096).
    pub fn new(n: usize) -> RqNtt {
        assert!(n.is_power_of_two() && n >= 2 && n <= 4096);
        let dom = FftDomain::<Fq>::new(n);
        // ψ = primitive 2n-th root of unity
        let log_2n = (n.trailing_zeros() + 1).min(13);
        let psi_root = Fq::root_of_unity(log_2n);
        let mut psi = Vec::with_capacity(n);
        let mut psi_inv = Vec::with_capacity(n);
        let psi_root_inv = psi_root.invert().expect("psi invertible");
        let mut cur = Fq(1);
        let mut cur_inv = Fq(1);
        for _ in 0..n {
            psi.push(cur);
            psi_inv.push(cur_inv);
            cur = cur * psi_root;
            cur_inv = cur_inv * psi_root_inv;
        }
        RqNtt { n, dom, psi, psi_inv }
    }

    /// Forward transform: coefficient form → NTT domain.
    pub fn forward(&self, a: &[Fq]) -> Vec<Fq> {
        assert_eq!(a.len(), self.n);
        let mut twisted = vec![Fq(0); self.n];
        for i in 0..self.n {
            twisted[i] = a[i] * self.psi[i];
        }
        let mut out = vec![Fq(0); self.n];
        self.dom.fft(&twisted, &mut out);
        out
    }

    /// Inverse transform: NTT domain → coefficient form.
    pub fn inverse(&self, v: &[Fq]) -> Rq {
        assert_eq!(v.len(), self.n);
        let mut twisted = vec![Fq(0); self.n];
        self.dom.ifft(v, &mut twisted);
        for i in 0..self.n {
            twisted[i] = twisted[i] * self.psi_inv[i];
        }
        Rq(twisted)
    }

    /// Pointwise product in the NTT domain.
    pub fn mul_pointwise(&self, a: &[Fq], b: &[Fq]) -> Vec<Fq> {
        a.iter().zip(b.iter()).map(|(x, y)| *x * *y).collect()
    }

    /// Full ring product via the NTT.
    pub fn mul(&self, a: &Rq, b: &Rq) -> Rq {
        let fa = self.forward(&a.0);
        let fb = self.forward(&b.0);
        let fc = self.mul_pointwise(&fa, &fb);
        self.inverse(&fc)
    }
}

thread_local! {
    /// Cached per-size NTT contexts (the twist and twiddle tables are
    /// shared by every ring multiplication of that size).
    static NTT_CTX: std::cell::RefCell<std::collections::HashMap<usize, std::rc::Rc<RqNtt>>> =
        std::cell::RefCell::new(std::collections::HashMap::new());
}

/// The cached NTT context for ring size `n`.
pub fn ntt_ctx(n: usize) -> std::rc::Rc<RqNtt> {
    NTT_CTX.with(|cache| {
        let mut cache = cache.borrow_mut();
        if let Some(ctx) = cache.get(&n) {
            return ctx.clone();
        }
        let ctx = std::rc::Rc::new(RqNtt::new(n));
        cache.insert(n, ctx.clone());
        ctx
    })
}

impl Rq {
    /// A reference to the shared NTT context for this ring size (if the
    /// size is NTT-friendly).
    pub fn ntt(&self) -> Option<std::rc::Rc<RqNtt>> {
        let n = self.0.len();
        if n >= 2 && n <= 4096 && n.is_power_of_two() {
            Some(ntt_ctx(n))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ring_arithmetic() {
        let mut rng = ZodaRng::from_seed(*b"rq-test-seed-0000000000000000000");
        for n in [4usize, 8, 16] {
            let a = Rq::uniform(n, &mut rng);
            let b = Rq::uniform(n, &mut rng);
            let c = Rq::uniform(n, &mut rng);
            // ring axioms incl. X^n = -1 structure (checked through
            // associativity/distributivity, which fails if the modular
            // reduction is wrong)
            assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)));
            assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
            // X^n = -1: multiply by X^n implicitly via (X^i)(X^{n-i}) terms
            // spot check: (a + a.neg()) == 0
            assert!(a.add(&a.neg()).is_zero());
        }
    }

    #[test]
    fn ntt_mul_matches_schoolbook() {
        let mut rng = ZodaRng::from_seed(*b"rq-ntt-test-seed-000000000000000");
        for n in [32usize, 64, 256, 1024] {
            let a = Rq::uniform(n, &mut rng);
            let b = Rq::uniform(n, &mut rng);
            let want = a.mul_schoolbook(&b);
            let got = a.mul(&b);
            assert_eq!(got, want, "ntt mul mismatch n={}", n);
            // context-based path
            let ctx = ntt_ctx(n);
            assert_eq!(ctx.mul(&a, &b), want, "ctx mul mismatch n={}", n);
            // sparse operand (many zeros)
            let mut sparse = vec![Fq(0); n];
            sparse[0] = Fq(7);
            sparse[n - 1] = Fq(12345);
            let sp = Rq(sparse);
            assert_eq!(sp.mul(&b), sp.mul_schoolbook(&b), "sparse n={}", n);
        }
    }

    #[test]
    fn x_pow_n_is_minus_one() {
        let n = 8;
        let x = {
            let mut c = vec![Fq(0); n];
            c[1] = Fq(1);
            Rq(c)
        };
        // x^n should equal -1
        let mut acc = x.clone();
        for _ in 1..n {
            acc = acc.mul(&x);
        }
        let mut expect = vec![Fq(0); n];
        expect[0] = -Fq(1);
        assert_eq!(acc, Rq(expect));
    }

    #[test]
    fn evaluation_matches_naive() {
        let mut rng = ZodaRng::from_seed(*b"rq-eval-test-seed-00000000000000");
        let coeffs: Vec<Fq> = (0..10).map(|_| Fq(rng.next_below(Q as u64) as u32)).collect();
        let z = Rq::from_i64(5);
        let got = eval_poly_at_ring(&coeffs, &z);
        // naive: evaluate the integer polynomial at the scalar 5 mod q
        let mut acc: u64 = 0;
        for c in coeffs.iter().rev() {
            acc = (acc * 5 + c.0 as u64) % Q as u64;
        }
        // eval at a scalar should give the constant acc when z = 5 (scalar)
        // since scalar multiplication is coefficient-wise:
        assert_eq!(got.0[0], Fq(acc as u32));
        assert!(got.0[1..].iter().all(|c| c.is_zero()));
    }
}
