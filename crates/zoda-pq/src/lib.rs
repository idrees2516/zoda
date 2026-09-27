//! # zoda-pq — lattice-based post-quantum polynomial commitments
//!
//! A KZG-shaped polynomial commitment built from lattice assumptions:
//! **Module-LWE / Module-SIS** in `R_q = Z_q[X]/(X^n + 1)` with
//! `q = 8380417`, instead of pairings and a trusted setup.
//!
//! ## Construction
//!
//! * **Commit** (BDLOP): `C = A·r + B·m` with `A, B` expanded from a
//!   public seed, `r` short CBD randomness and `m` the coefficient
//!   chunks of the polynomial. Hiding rests on Module-LWE, binding on
//!   Module-SIS; no trusted setup, no pairing — plausibly post-quantum
//!   secure.
//! * **Open** (Lyubashevsky-style σ-protocol): the prover publishes a
//!   masking commitment `C₁ = A·r₁ + B·m₁` and `v₁ = eval(m₁, ζ)`,
//!   derives a narrow ring challenge `x` (coefficients in {−1,0,1},
//!   weight ≤ 32) by Fiat–Shamir, and responds with
//!   `r₂ = r₁ + x·r`, `m₂ = m₁ + x·m`. The verifier checks
//!   `A·r₂ + B·m₂ = C₁ + x·C` and `eval(m₂, ζ) = v₁ + x·v`.
//! * Evaluations are at ring points `ζ ∈ R_q`; a scalar `z` embeds as
//!   the constant ring element, recovering the classical evaluation.
//!
//! ## Parameter sets (NIST levels)
//!
//! See `params.rs` and `docs/POST_QUANTUM.md`. The parameter sizes are
//! first-estimate configurations targeting MSIS/MLWE core-SVP hardness
//! at each level; the implementation has **not** been independently
//! audited and should be treated as a research prototype.

pub mod params;
pub mod pcs;
pub mod rq;

pub use params::PqParams;
pub use pcs::{
    is_admissible_ring_point, ring_point, ring_point_from_seed, LatticeCommitment, LatticePcs,
    LatticeProof, OpeningKey, CHALLENGE_WEIGHT,
};
pub use rq::{fq_to_poly, Rq, RqVector};

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_math::{Fq, PrimeField, ZodaRng};

    #[test]
    fn commitment_roundtrip_fast() {
        let params = PqParams::fast();
        let mut rng = ZodaRng::from_seed(*b"pq-pcs-test-seed-000000000000000");
        let poly: Vec<Fq> = (0..64).map(|_| Fq(rng.next_below(rq::Q as u64) as u32)).collect();
        let (com, key) = LatticePcs::commit(&params, &poly, &mut rng);

        let zeta = ring_point_from_seed(b"test-open-point", params.n);
        let v = LatticePcs::eval_at_ring(&poly, &zeta);
        let proof = LatticePcs::open(&params, &key, &zeta, &v, &mut rng);
        assert!(LatticePcs::verify(&params, &com, &zeta, &v, &proof));

        // wrong point fails
        let zeta2 = ring_point_from_seed(b"test-open-point-2", params.n);
        assert!(!LatticePcs::verify(&params, &com, &zeta2, &v, &proof));

        // wrong value fails
        let mut v_bad = v.clone();
        v_bad.0[0] = v_bad.0[0] + Fq(1);
        assert!(!LatticePcs::verify(&params, &com, &zeta, &v_bad, &proof));

        // tampered proof fails
        let mut bad = proof.clone();
        bad.m2.0[0] = bad.m2.0[0].add(&Rq(vec![Fq(1)]));
        assert!(!LatticePcs::verify(&params, &com, &zeta, &v, &bad));
    }

    #[test]
    fn hiding_and_determinism() {
        let params = PqParams::fast();
        let mut rng = ZodaRng::from_seed(*b"pq-pcs-hide-seed-000000000000000");
        let poly: Vec<Fq> = (0..32).map(|_| Fq(rng.next_below(rq::Q as u64) as u32)).collect();
        let (c1, _) = LatticePcs::commit(&params, &poly, &mut rng);
        let (c2, _) = LatticePcs::commit(&params, &poly, &mut rng);
        assert_ne!(c1.c, c2.c);
        let mut r1 = ZodaRng::from_seed(*b"pq-pcs-det-seed-0000000000000000");
        let mut r2 = ZodaRng::from_seed(*b"pq-pcs-det-seed-0000000000000000");
        let (d1, _) = LatticePcs::commit(&params, &poly, &mut r1);
        let (d2, _) = LatticePcs::commit(&params, &poly, &mut r2);
        assert_eq!(d1.c, d2.c);
    }

    #[test]
    fn high_degree_chunks() {
        let params = PqParams::fast();
        let mut rng = ZodaRng::from_seed(*b"pq-pcs-hd-seed-00000000000000000");
        let poly: Vec<Fq> =
            (0..5 * params.n + 17).map(|_| Fq(rng.next_below(rq::Q as u64) as u32)).collect();
        let (com, key) = LatticePcs::commit(&params, &poly, &mut rng);
        let zeta = ring_point_from_seed(b"test-open-point-hd", params.n);
        let v = LatticePcs::eval_at_ring(&poly, &zeta);
        let proof = LatticePcs::open(&params, &key, &zeta, &v, &mut rng);
        assert!(LatticePcs::verify(&params, &com, &zeta, &v, &proof));
    }
}
