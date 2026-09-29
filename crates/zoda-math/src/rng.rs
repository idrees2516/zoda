//! Deterministic SHA-256-based CSPRNG (counter mode) plus OS seeding.
//!
//! Deterministic seeds make the entire protocol testable and reproducible;
//! `from_os` mixes in `/dev/urandom` for production key material.

use crate::sha256::sha256;
use crate::PrimeField;
use crate::Fr;

#[derive(Clone)]
pub struct ZodaRng {
    state: [u8; 32],
    counter: u64,
    block: [u8; 32],
    used: usize,
}

impl ZodaRng {
    /// Deterministic RNG from a 32-byte seed.
    pub fn from_seed(seed: [u8; 32]) -> Self {
        let state = sha256(&seed);
        ZodaRng {
            state,
            counter: 0,
            block: [0; 32],
            used: 32, // force refill
        }
    }

    /// Non-deterministic RNG seeded from the operating system.
    #[cfg(feature = "std")]
    pub fn from_os() -> Self {
        use std::io::Read;
        let mut seed = [0u8; 32];
        if let Ok(mut f) = std::fs::File::open("/dev/urandom") {
            let _ = f.read_exact(&mut seed);
        } else {
            // fallback: time + address entropy
            let t = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0);
            seed[..16].copy_from_slice(&t.to_le_bytes());
        }
        Self::from_seed(seed)
    }

    fn refill(&mut self) {
        let mut input = Vec::with_capacity(40);
        input.extend_from_slice(&self.state);
        input.extend_from_slice(&self.counter.to_le_bytes());
        self.block = sha256(&input);
        self.counter = self.counter.wrapping_add(1);
        self.used = 0;
    }

    pub fn next_bytes(&mut self, out: &mut [u8]) {
        let mut i = 0;
        while i < out.len() {
            if self.used == 32 {
                self.refill();
            }
            out[i] = self.block[self.used];
            self.used += 1;
            i += 1;
        }
    }

    pub fn next_u32(&mut self) -> u32 {
        let mut b = [0u8; 4];
        self.next_bytes(&mut b);
        u32::from_le_bytes(b)
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut b = [0u8; 8];
        self.next_bytes(&mut b);
        u64::from_le_bytes(b)
    }

    /// Uniform `Fr` element; optionally rejecting zero.
    pub fn next_fr(&mut self, reject_zero: bool) -> Fr {
        loop {
            let mut b = [0u8; 32];
            self.next_bytes(&mut b);
            // interpret LE bytes mod r (slight bias < 2^-101, acceptable
            // for all zoda uses; noted in the security docs)
            let x = Fr::from_le_bytes_mod_order(&b);
            if !x.is_zero() || !reject_zero {
                return x;
            }
        }
    }

    /// Uniform field element of any `PrimeField` (LE bytes mod order).
    pub fn next_field<F: crate::PrimeField>(&mut self) -> F {
        let mut b = [0u8; 32];
        self.next_bytes(&mut b);
        F::from_le_bytes_mod_order(&b)
    }

    /// Uniform integer in `[0, bound)` (rejection sampling, unbiased).
    pub fn next_below(&mut self, bound: u64) -> u64 {
        debug_assert!(bound > 0);
        if bound == 1 {
            return 0;
        }
        // number of valid values per 64-bit draw
        let buckets = u64::MAX / bound;
        let limit = buckets * bound;
        loop {
            let x = self.next_u64();
            if x < limit {
                return x % bound;
            }
        }
    }

    /// Uniform integer in `[lo, hi)`; panics if `hi <= lo`.
    pub fn next_range(&mut self, lo: usize, hi: usize) -> usize {
        assert!(hi > lo);
        lo + self.next_below((hi - lo) as u64) as usize
    }

    /// Fisher-Yates shuffle of `0..n`.
    pub fn shuffle_indices(&mut self, n: usize) -> Vec<usize> {
        let mut v: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = self.next_below((i + 1) as u64) as usize;
            v.swap(i, j);
        }
        v
    }

    /// Uniform `Fq` (lattice field) element.
    pub fn next_fq(&mut self) -> crate::Fq {
        crate::Fq(self.next_below(crate::fq::Q as u64) as u32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic() {
        let mut a = ZodaRng::from_seed(*b"determinism-test-seed-0000000001");
        let mut b = ZodaRng::from_seed(*b"determinism-test-seed-0000000001");
        let mut c = ZodaRng::from_seed(*b"determinism-test-seed-0000000002");
        let va: Vec<u64> = (0..10).map(|_| a.next_u64()).collect();
        let vb: Vec<u64> = (0..10).map(|_| b.next_u64()).collect();
        let vc: Vec<u64> = (0..10).map(|_| c.next_u64()).collect();
        assert_eq!(va, vb);
        assert_ne!(va, vc);
    }

    #[test]
    fn ranges_and_shuffles() {
        let mut r = ZodaRng::from_seed(*b"range-test-seed-0000000000000000");
        for _ in 0..1000 {
            let x = r.next_below(7);
            assert!(x < 7);
        }
        for _ in 0..100 {
            let x = r.next_range(3, 9);
            assert!((3..9).contains(&x));
        }
        let s = r.shuffle_indices(64);
        let mut sorted = s.clone();
        sorted.sort();
        assert_eq!(sorted, (0..64).collect::<Vec<_>>());
    }
}
