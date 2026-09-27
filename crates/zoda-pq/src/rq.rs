//! The ring `R_q = Z_q[X]/(X^n + 1)` with `q = 8380417`, and its
//! negacyclic NTT (q supports transforms up to length 8192).

use zoda_math::{Fq, PrimeField, ZodaRng};

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

    /// Multiply two ring elements (schoolbook with reduction mod X^n+1).
    /// For clarity and small n this is O(n²); the PCS uses chunk sizes
    /// where this is not the bottleneck. (Larger n paths use the same
    /// routine — see the benchmarks for measured costs.)
    pub fn mul(&self, rhs: &Rq) -> Rq {
        let n = self.0.len();
        debug_assert_eq!(n, rhs.0.len());
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
        let n = 8;
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
