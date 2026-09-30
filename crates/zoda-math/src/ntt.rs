//! Iterative radix-2 number-theoretic transforms with per-stage twiddle
//! tables and bit-reversal-free natural-order input/output.
//!
//! The transform matches the c-kzg-4844 convention exactly:
//!
//! ```text
//! fft(x, n)[i] = Σ_j x[j] · ω_n^{i·j},  ω_n = ω_{2^13}^{8192/n}
//! ```
//!
//! with `ω_{2^13} = PRIMITIVE_ROOT_OF_UNITY^((r-1)/8192)` for `Fr`, so
//! evaluations interoperate with Ethereum consensus structures bit-for-bit.
//!
//! Optimisations over the textbook transform:
//! * per-stage contiguous twiddle tables (sequential memory access),
//! * one precomputed bit-reversal permutation,
//! * `fft_in_place` variants with zero allocation for steady-state paths,
//! * coset variants used by EIP-7594 recovery.

use crate::PrimeField;

/// A field over which FFTs can be computed (a PrimeField with a known
/// two-adic structure and a generator).
pub trait FftField: PrimeField {
    /// Primitive root of unity of order `2^k`.
    fn fft_root_of_unity(k: u32) -> Self;
}

impl FftField for crate::Fr {
    fn fft_root_of_unity(k: u32) -> Self {
        Self::root_of_unity(k)
    }
}
impl FftField for crate::Goldilocks {
    fn fft_root_of_unity(k: u32) -> Self {
        Self::root_of_unity(k)
    }
}
impl FftField for crate::Fq {
    fn fft_root_of_unity(k: u32) -> Self {
        Self::root_of_unity(k)
    }
}

/// A precomputed FFT domain of size `n = 2^log_n`.
pub struct FftDomain<F: FftField> {
    /// transform size
    pub n: usize,
    #[allow(dead_code)]
    log_n: u32,
    /// per-stage twiddles: `stage_twiddles[s][j] = ω_{2^{s+1}}^j`
    /// laid out contiguously in one vector; stage s occupies
    /// `[2^s - 1 .. 2^{s+1} - 1)`.
    fwd: Vec<F>,
    inv: Vec<F>,
    /// n^{-1}
    inv_n: F,
    /// bit-reversal permutation of length n
    brp: Vec<u32>,
    /// natural-order roots of unity ω_n^j (size n+1, last entry = 1)
    roots: Vec<F>,
}

impl<F: FftField> FftDomain<F> {
    /// Build the domain for `n` (must be a power of two).
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two() && n >= 2, "n must be a power of two ≥ 2");
        let log_n = n.trailing_zeros();
        // total twiddle slots = 1 + 2 + 4 + ... + n/2 = n - 1
        let mut fwd = Vec::with_capacity(n - 1);
        for s in 0..log_n {
            let m = 1usize << (s + 1);
            let w = F::fft_root_of_unity(s + 1);
            for j in 0..m / 2 {
                fwd.push(w.pow_u64(j as u64));
            }
        }
        let mut inv = Vec::with_capacity(n - 1);
        for s in 0..log_n {
            let m = 1usize << (s + 1);
            let w = F::fft_root_of_unity(s + 1).invert().unwrap();
            for j in 0..m / 2 {
                inv.push(w.pow_u64(j as u64));
            }
        }
        let inv_n = F::from_u64(n as u64).invert().unwrap();
        let brp = bit_reversal_permutation(n);
        // natural order roots of unity, n+1 entries (roots[n] = 1)
        let wn = F::fft_root_of_unity(log_n);
        let mut roots = Vec::with_capacity(n + 1);
        let mut cur = F::ONE;
        for _ in 0..n {
            roots.push(cur);
            cur = cur * wn;
        }
        roots.push(F::ONE);
        Self {
            n,
            log_n,
            fwd,
            inv,
            inv_n,
            brp,
            roots,
        }
    }

    /// The transform size.
    #[inline]
    pub fn size(&self) -> usize {
        self.n
    }

    /// Natural-order roots of unity `[1, ω, ω², …, ω^{n-1}, 1]` (n+1 entries).
    #[inline]
    pub fn roots_of_unity(&self) -> &[F] {
        &self.roots
    }

    /// Bit-reversal permutation of `[0, n)`.
    #[inline]
    pub fn brp_indices(&self) -> &[u32] {
        &self.brp
    }

    /// Forward transform: `out[i] = Σ_j in[j]·ω^{ij}`. Both slices must have
    /// length `n`; they may not alias.
    pub fn fft(&self, input: &[F], out: &mut [F]) {
        debug_assert_eq!(input.len(), self.n);
        debug_assert_eq!(out.len(), self.n);
        out[..self.n].copy_from_slice(&input[..self.n]);
        self.fft_in_place(out);
    }

    /// Forward transform in place. `a` must have length exactly `n`.
    pub fn fft_in_place(&self, a: &mut [F]) {
        assert_eq!(a.len(), self.n, "fft_in_place length mismatch");
        // bit-reversal permutation
        for i in 0..self.n {
            let j = self.brp[i] as usize;
            if j > i {
                a.swap(i, j);
            }
        }
        // butterflies
        let mut len = 2usize;
        let mut tw_start = 0usize;
        while len <= self.n {
            let half = len / 2;
            let twiddles = &self.fwd[tw_start..tw_start + half];
            let mut start = 0;
            while start < self.n {
                let (lo, hi) = a.split_at_mut(start + half);
                let (lo_blk, hi_blk) = {
                    let (_, lo_blk) = lo.split_at_mut(start);
                    let (_, hi_blk) = hi.split_at_mut(0);
                    (lo_blk, hi_blk)
                };
                for j in 0..half {
                    let w = twiddles[j];
                    let u = lo_blk[j];
                    let v = hi_blk[j] * w;
                    lo_blk[j] = u + v;
                    hi_blk[j] = u - v;
                }
                start += len;
            }
            tw_start += half;
            len <<= 1;
        }
    }

    /// Forward transform with zero-padding: `input` may be shorter than
    /// the domain size (trailing coefficients are zero).
    pub fn fft_padded(&self, input: &[F], out: &mut [F]) {
        debug_assert!(input.len() <= self.n);
        let mut tmp = vec![F::zero(); self.n];
        tmp[..input.len().min(self.n)].copy_from_slice(&input[..input.len().min(self.n)]);
        self.fft_in_place(&mut tmp);
        out[..self.n].copy_from_slice(&tmp);
    }

    /// Inverse transform (scaled by `1/n`).
    pub fn ifft(&self, input: &[F], out: &mut [F]) {
        out[..self.n].copy_from_slice(&input[..self.n]);
        self.ifft_in_place(out);
    }

    /// Inverse transform in place (scaled by `1/n`).
    pub fn ifft_in_place(&self, a: &mut [F]) {
        assert_eq!(a.len(), self.n, "ifft_in_place length mismatch");
        for i in 0..self.n {
            let j = self.brp[i] as usize;
            if j > i {
                a.swap(i, j);
            }
        }
        let mut len = 2usize;
        let mut tw_start = 0usize;
        while len <= self.n {
            let half = len / 2;
            let twiddles = &self.inv[tw_start..tw_start + half];
            let mut start = 0;
            while start < self.n {
                let (lo, hi) = a.split_at_mut(start + half);
                let (lo_blk, hi_blk) = {
                    let (_, lo_blk) = lo.split_at_mut(start);
                    let (_, hi_blk) = hi.split_at_mut(0);
                    (lo_blk, hi_blk)
                };
                for j in 0..half {
                    let w = twiddles[j];
                    let u = lo_blk[j];
                    let v = hi_blk[j] * w;
                    lo_blk[j] = u + v;
                    hi_blk[j] = u - v;
                }
                start += len;
            }
            tw_start += half;
            len <<= 1;
        }
        for x in a.iter_mut() {
            *x = *x * self.inv_n;
        }
    }

    /// Evaluate the coefficient-form polynomial `coeffs` (length ≤ n) on the
    /// coset `shift · ω^i`: `out[i] = P(shift·ω^i)`.
    pub fn coset_fft(&self, coeffs: &[F], out: &mut [F], shift: F) {
        debug_assert!(coeffs.len() <= self.n);
        let mut tmp = vec![F::zero(); self.n];
        tmp[..coeffs.len()].copy_from_slice(coeffs);
        // prescale: P(shift·ω^i) = Σ c_k shift^k ω^{ik}
        let mut sh = F::ONE;
        for c in tmp.iter_mut().take(coeffs.len()) {
            *c = *c * sh;
            sh = sh * shift;
        }
        self.fft_in_place(&mut tmp);
        out[..self.n].copy_from_slice(&tmp);
    }

    /// Inverse of [`coset_fft`] (scaled by `1/n`).
    pub fn coset_ifft(&self, evals: &[F], out_coeffs: &mut [F], shift: F) {
        let mut tmp = evals[..self.n].to_vec();
        self.ifft_in_place(&mut tmp);
        let inv_shift = shift.invert().unwrap();
        let mut sh = F::ONE;
        for c in tmp.iter_mut() {
            *c = *c * sh;
            sh = sh * inv_shift;
        }
        out_coeffs.copy_from_slice(&tmp);
    }
}

/// Bit-reversal permutation of length `n` (a power of two):
/// `out[i] = reverse_bits(i, n)`.
pub fn bit_reversal_permutation_typed<F: Copy>(v: &mut [F]) {
    let n = v.len();
    debug_assert!(n.is_power_of_two());
    if n <= 1 {
        return;
    }
    for i in 0..n {
        let j = reverse_bits_limited(n, i);
        if j > i {
            v.swap(i, j);
        }
    }
}

/// Reverse the low `log2(order)` bits of `n`.
#[inline]
pub fn reverse_bits_limited(order: usize, n: usize) -> usize {
    debug_assert!(order.is_power_of_two());
    let width = order.trailing_zeros();
    n.reverse_bits() >> (usize::BITS - width)
}

fn bit_reversal_permutation(n: usize) -> Vec<u32> {
    (0..n).map(|i| reverse_bits_limited(n, i) as u32).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Fr;

    #[test]
    fn fft_matches_naive_dft() {
        let mut rng = crate::ZodaRng::from_seed(*b"fft-test-seed-000000000000000000");
        for log in 1..=9u32 {
            let n = 1usize << log;
            let d = FftDomain::<Fr>::new(n);
            let a: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
            let mut out = vec![Fr::zero(); n];
            d.fft(&a, &mut out);
            // naive DFT with the domain's own roots
            let roots = &d.roots_of_unity()[..n];
            for i in 0..n {
                let mut acc = Fr::zero();
                let mut w = Fr::ONE;
                for j in 0..n {
                    // ω^{ij} = (ω^i)^j
                    acc = acc + a[j] * w;
                    w = w * roots[i];
                }
                assert_eq!(acc, out[i], "DFT mismatch n={} i={}", n, i);
            }
            // inverse round trip
            let mut back = vec![Fr::zero(); n];
            d.ifft(&out, &mut back);
            assert_eq!(back, a);
        }
    }

    #[test]
    fn coset_fft_matches_definition() {
        let mut rng = crate::ZodaRng::from_seed(*b"coset-test-seed-0000000000000000");
        let n = 128;
        let d = FftDomain::<Fr>::new(n);
        let shift = rng.next_fr(false);
        let coeffs: Vec<Fr> = (0..n / 2).map(|_| rng.next_fr(false)).collect();
        let mut out = vec![Fr::zero(); n];
        d.coset_fft(&coeffs, &mut out, shift);
        let roots = &d.roots_of_unity()[..n];
        for i in 0..n {
            // P(shift * ω^i) by Horner
            let x = shift * roots[i];
            let mut acc = Fr::zero();
            for c in coeffs.iter().rev() {
                acc = acc * x + *c;
            }
            assert_eq!(acc, out[i], "coset mismatch i={}", i);
        }
        let mut back = vec![Fr::zero(); n];
        d.coset_ifft(&out, &mut back, shift);
        assert_eq!(&back[..coeffs.len()], &coeffs[..]);
    }

    #[test]
    fn goldilocks_and_fq_domains() {
        let d = FftDomain::<crate::Goldilocks>::new(64);
        let a: Vec<crate::Goldilocks> = (0..64)
            .map(|i| crate::Goldilocks((i as u64 * 7919 + 3) % crate::goldilocks::GOLDILOCKS_MODULUS))
            .collect();
        let mut out = vec![crate::Goldilocks(0); 64];
        d.fft(&a, &mut out);
        // spot check: out[0] = sum
        let sum = a.iter().fold(crate::Goldilocks(0), |x, y| x + *y);
        assert_eq!(out[0], sum);
        let dq = FftDomain::<crate::Fq>::new(64);
        let b: Vec<crate::Fq> = (0..64).map(|i| crate::Fq((i * 97 + 5) % crate::fq::Q)).collect();
        let mut oq = vec![crate::Fq(0); 64];
        dq.fft(&b, &mut oq);
        let sumq = b.iter().fold(crate::Fq(0), |x, y| x + *y);
        assert_eq!(oq[0], sumq);
    }
}
