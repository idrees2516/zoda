//! # zoda-bls
//!
//! A dependency-free BLS12-381 library for the zoda data-availability
//! stack: the 381-bit base field, the Fp2/Fp6/Fp12 tower, G1/G2 groups
//! with ZCash-format point compression, the optimal ate pairing
//! (Beuchat–López-Tehraníja... Miller loop + Fuentes et al. final
//! exponentiation), RFC 9380 hash-to-curve and IETF BLS signatures.
//!
//! ## Design notes
//!
//! * All tower and Frobenius constants are either standard small values or
//!   **derived and self-verified at initialisation** (see
//!   [`pairing::frob`] internals) — no hand-copied 381-bit tables beyond
//!   the canonical curve generators.
//! * Point arithmetic uses the complete Renes–Costello–Batina formulas —
//!   no edge cases on identity/double/negation inputs.
//! * Hash-to-curve constants are mechanically extracted from the
//!   zkcrypto/bls12_381 reference tables (see `scripts/extract_h2c_consts.py`)
//!   and validated against the RFC 9380 test vectors.

pub mod curve;
pub mod fp;
pub mod fp2;
pub mod fp6;
pub mod fp12;
pub mod g1;
pub mod g2;
pub mod endomorphism;
pub mod hash;
pub mod h2c_consts;
pub mod pairing;
pub mod sig;

pub use fp::Fp;
pub use fp2::Fp2;
pub use fp6::Fp6;
pub use fp12::Fp12;
pub use g1::{G1Affine, G1Projective};
pub use g2::{G2Affine, G2Projective};
pub use pairing::{pairing, pairing_check, G2Prepared};
pub use sig::{
    aggregate, aggregate_pks, verify, verify_aggregate, verify_batch, verify_batch_strict,
    verify_strict, PublicKey, SecretKey, Signature,
};
