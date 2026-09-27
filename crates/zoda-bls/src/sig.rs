//! BLS signatures on BLS12-381 (IETF draft-irtf-cfrg-bls-signature-04,
//! ciphersuite `BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_` — the
//! Ethereum-flavoured scheme: G1 public keys, G2 signatures).
//!
//! * `SkToPk` — pk = sk · g1
//! * `CoreSign` — σ = sk · H(msg)
//! * `CoreVerify` — e(pk, H(msg)) · e(-g1, σ) == 1
//!
//! `verify` performs on-curve validation of pk and σ; `verify_strict`
//! additionally enforces r-subgroup membership (defends against small
//! subgroup attacks when keys are gathered from untrusted aggregations).

use crate::g1::{G1Affine, G1Projective};
use crate::g2::{G2Affine, G2Projective};
use crate::hash::hash_to_curve_g2;
use crate::pairing::{pairing_check, G2Prepared};
use zoda_math::Fr;

pub const DST: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";

/// A BLS secret key (a scalar in Fr*, stored reduced).
#[derive(Clone, Copy)]
pub struct SecretKey(pub Fr);

/// A BLS public key (G1 affine).
pub type PublicKey = G1Affine;
/// A BLS signature (G2 affine).
pub type Signature = G2Affine;

impl SecretKey {
    /// Deterministic key derivation from 32 bytes of seed material
    /// (`HashToScalar`-style: SHA-256 based key derivation per the draft's
    /// `KeyGen` information model).
    pub fn from_seed(seed: &[u8]) -> SecretKey {
        // I2OSP(0,1) || seed mapped to Fr
        let mut input = Vec::with_capacity(seed.len() + 1);
        input.push(0u8);
        input.extend_from_slice(seed);
        SecretKey(Fr::hash_to_fr(&input))
    }

    /// Build directly from a scalar.
    pub fn from_scalar(fr: Fr) -> SecretKey {
        SecretKey(fr)
    }

    /// Derive the matching public key.
    pub fn public_key(&self) -> PublicKey {
        G1Projective::generator().mul_fr(&self.0).to_affine()
    }

    /// Sign a message.
    pub fn sign(&self, msg: &[u8]) -> Signature {
        let h = hash_to_curve_g2(msg, DST);
        h.mul_fr(&self.0).to_affine()
    }
}

/// Verify a signature against a public key (on-curve checks included).
pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    if pk.infinity || sig.infinity {
        return false;
    }
    if !pk.is_on_curve() || !sig.is_on_curve() {
        return false;
    }
    let h = hash_to_curve_g2(msg, DST).to_affine();
    let neg_g1 = G1Affine::generator().neg();
    // e(pk, H) · e(-g1, σ) == 1
    let h_prepared = G2Prepared::from(h);
    let sig_prepared = G2Prepared::from(*sig);
    pairing_check(&[(pk, &h_prepared), (&neg_g1, &sig_prepared)])
}

/// Verify with full subgroup checks (r-torsion membership).
pub fn verify_strict(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    if pk.infinity || sig.infinity {
        return false;
    }
    if !pk.is_on_curve() || !sig.is_on_curve() {
        return false;
    }
    // pk must be in the r-subgroup of G1: pk·r == O
    let r_limbs = Fr::MODULUS;
    if !pk.to_projective().mul_limbs(&r_limbs).is_identity() {
        return false;
    }
    // σ must be in the r-subgroup of G2: σ·r == O
    if !sig.to_projective().mul_limbs(&r_limbs).is_identity() {
        return false;
    }
    verify(pk, msg, sig)
}

/// Aggregate multiple signatures over the same message.
pub fn aggregate(signatures: &[Signature]) -> Option<Signature> {
    if signatures.is_empty() {
        return None;
    }
    let mut acc = G2Projective::identity();
    for s in signatures {
        if s.infinity {
            return None;
        }
        acc = acc + s.to_projective();
    }
    Some(acc.to_affine())
}

/// Verify an aggregate signature from `pks` all over the same `msg`.
pub fn verify_aggregate(pks: &[PublicKey], msg: &[u8], sig: &Signature) -> bool {
    if pks.is_empty() {
        return false;
    }
    let mut agg = G1Projective::identity();
    for pk in pks {
        if pk.infinity || !pk.is_on_curve() {
            return false;
        }
        agg = agg + pk.to_projective();
    }
    verify(&agg.to_affine(), msg, sig)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sign_verify_roundtrip() {
        let sk = SecretKey::from_seed(b"zoda bls test key material 0001");
        let pk = sk.public_key();
        let msg = b"attack at dawn";
        let sig = sk.sign(msg);
        assert!(verify(&pk, msg, &sig));
        assert!(verify_strict(&pk, msg, &sig));
        // wrong message fails
        assert!(!verify(&pk, b"attack at dusk", &sig));
        // wrong key fails
        let sk2 = SecretKey::from_seed(b"zoda bls test key material 0002");
        assert!(!verify(&sk2.public_key(), msg, &sig));
    }

    #[test]
    fn aggregation() {
        let sks: Vec<_> = (0..3)
            .map(|i| SecretKey::from_seed(format!("zoda agg key {}", i).as_bytes()))
            .collect();
        let msg = b"shared message";
        let sigs: Vec<_> = sks.iter().map(|sk| sk.sign(msg)).collect();
        let agg = aggregate(&sigs).unwrap();
        let pks: Vec<_> = sks.iter().map(|sk| sk.public_key()).collect();
        assert!(verify_aggregate(&pks, msg, &agg));
        // mixing in a signature over a different message breaks it
        let other = sks[0].sign(b"other message");
        let bad_agg = aggregate(&[sigs[0], sigs[1], other]).unwrap();
        assert!(!verify_aggregate(&pks, msg, &bad_agg));
    }

    #[test]
    fn deterministic_keys() {
        let a = SecretKey::from_seed(b"same seed");
        let b = SecretKey::from_seed(b"same seed");
        assert_eq!(a.0, b.0);
        assert_eq!(a.public_key(), b.public_key());
    }
}
