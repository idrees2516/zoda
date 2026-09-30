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

/// Aggregate public keys (multiset sum, same-message aggregation helper).
pub fn aggregate_pks(pks: &[PublicKey]) -> Option<PublicKey> {
    if pks.is_empty() {
        return None;
    }
    let mut acc = G1Projective::identity();
    for pk in pks {
        if pk.infinity {
            return None;
        }
        acc = acc + pk.to_projective();
    }
    Some(acc.to_affine())
}

// ---------------------------------------------------------------------------
// Batch verification
// ---------------------------------------------------------------------------

/// Random per-item coefficients derived from a transcript of every input
/// (unpredictable before the full batch is fixed — grinding a bad batch
/// that passes would take ~2^128 attempts per 128-bit coefficient).
/// Returns two independent challenge sets (pairing, subgroup) so the two
/// checks cannot be correlated.
fn batch_challenges(
    items: &[(&PublicKey, &[u8], &Signature)],
) -> (Vec<zoda_math::Fr>, Vec<zoda_math::Fr>) {
    let mut transcript = Vec::new();
    for (pk, msg, sig) in items {
        transcript.extend_from_slice(&pk.to_compressed());
        transcript.extend_from_slice(&(msg.len() as u64).to_be_bytes());
        transcript.extend_from_slice(msg);
        transcript.extend_from_slice(&sig.to_compressed());
    }
    transcript.extend_from_slice(&(items.len() as u64).to_be_bytes());
    let seed = zoda_math::sha256::sha256(&transcript);
    let mut rng = zoda_math::ZodaRng::from_seed(seed);
    let pairing = (0..items.len()).map(|_| rng.next_fr(true)).collect();
    let subgroup = (0..items.len()).map(|_| rng.next_fr(true)).collect();
    (pairing, subgroup)
}

/// Verify a batch of independent (pk, msg, sig) triples with one final
/// exponentiation amortised over all items:
///
/// `Π e(rᵢ·pkᵢ, H(mᵢ)) · e(−g1, Σ rᵢ·σᵢ) == 1`
///
/// for transcript-derived random `rᵢ`. `n` Miller loops replace `2n`, and
/// the single final exponentiation replaces `n` — roughly 2× per item over
/// independent `verify` calls, with the same acceptance semantics
/// (on-curve checks; h-torsion components pair to 1, exactly as in single
/// verification, because the per-item scalar multiplications use the plain
/// path on untrusted points).
pub fn verify_batch(items: &[(&PublicKey, &[u8], &Signature)]) -> bool {
    if items.is_empty() {
        return true;
    }
    let challenges = batch_challenges(items);
    let mut terms: Vec<(G1Affine, G2Prepared)> = Vec::with_capacity(items.len() + 1);
    let mut sig_sum = G2Projective::identity();
    for ((pk, msg, sig), r) in items.iter().zip(challenges.0.iter()) {
        if pk.infinity || sig.infinity {
            return false;
        }
        if !pk.is_on_curve() || !sig.is_on_curve() {
            return false;
        }
        // plain (non-GLV) multiplication: pk/σ are untrusted, and the
        // h-torsion behaviour must match single verification exactly
        let rp = pk.to_projective().mul_limbs(&r.to_repr()).to_affine();
        terms.push((rp, G2Prepared::from(hash_to_curve_g2(msg, DST).to_affine())));
        let rs = sig.to_projective().mul_limbs(&r.to_repr());
        sig_sum = sig_sum + rs;
    }
    let neg_g1 = G1Affine::generator().neg();
    terms.push((neg_g1, G2Prepared::from(sig_sum.to_affine())));
    let term_refs: Vec<(&G1Affine, &G2Prepared)> = terms.iter().map(|(a, b)| (a, b)).collect();
    pairing_check(&term_refs)
}

/// Batch verification with full subgroup soundness: one Pippenger-based G1
/// batch subgroup check over the public keys, one combined G2 check over
/// the signatures, then the same product-pairing check as [`verify_batch`]
/// — with GLV-accelerated scalar multiplications, since every point has
/// been subgroup-verified first.
pub fn verify_batch_strict(items: &[(&PublicKey, &[u8], &Signature)]) -> bool {
    if items.is_empty() {
        return true;
    }
    let pks: Vec<PublicKey> = items.iter().map(|(pk, _, _)| **pk).collect();
    let sigs: Vec<Signature> = items.iter().map(|(_, _, sig)| **sig).collect();
    if !crate::endomorphism::batch_subgroup_check_g2(&sigs) {
        return false;
    }
    // individual on-curve and per-point G1 subgroup check via the batch
    // combination (points must be on the curve first — decompression
    // guarantees it, but re-check defensively for in-curve points)
    for pk in &pks {
        if pk.infinity || !pk.is_on_curve() {
            return false;
        }
    }
    // batch the G1 subgroup checks: S = Σ cᵢ·pkᵢ with INDEPENDENT
    // challenges and PLAIN multiplications (correct semantics on any
    // point of E(Fp) — GLV would be unsound pre-verification), then a
    // single [r]S == O check
    let challenges = batch_challenges(items);
    {
        let mut s = G1Projective::identity();
        for (pk, c) in pks.iter().zip(challenges.1.iter()) {
            if pk.infinity {
                continue;
            }
            s = s + pk.to_projective().mul_limbs(&c.to_repr());
        }
        if !s.mul_limbs(&zoda_math::Fr::MODULUS).is_identity() {
            return false;
        }
    }
    // every point is now subgroup-verified: the pairing product may use
    // the GLV fast path
    let mut terms: Vec<(G1Affine, G2Prepared)> = Vec::with_capacity(items.len() + 1);
    let mut sig_sum = G2Projective::identity();
    for ((pk, msg, sig), r) in items.iter().zip(challenges.0.iter()) {
        if sig.infinity {
            return false;
        }
        let rp =
            crate::endomorphism::mul_g1_public(&pk.to_projective(), r).to_affine();
        terms.push((rp, G2Prepared::from(hash_to_curve_g2(msg, DST).to_affine())));
        let rs = crate::endomorphism::mul_g2_public(&sig.to_projective(), r);
        sig_sum = sig_sum + rs;
    }
    let neg_g1 = G1Affine::generator().neg();
    terms.push((neg_g1, G2Prepared::from(sig_sum.to_affine())));
    let term_refs: Vec<(&G1Affine, &G2Prepared)> = terms.iter().map(|(a, b)| (a, b)).collect();
    pairing_check(&term_refs)
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
    use crate::fp::Fp;
    use crate::fp2::Fp2;
    use crate::g2::G2Config;
    use crate::curve::CurveConfig;

    fn random_bad_g2_point(seed: u64) -> Signature {
        // a point on E'(Fp2) outside G2 (on-curve, wrong subgroup)
        let mut rng = zoda_math::ZodaRng::from_seed(*b"sig-batch-bad-000000000000000000");
        let _ = seed;
        loop {
            let mut b0 = [0u8; 32];
            let mut b1 = [0u8; 32];
            rng.next_bytes(&mut b0);
            rng.next_bytes(&mut b1);
            let x = Fp2::new(Fp::from_le_bytes32(&b0), Fp::from_le_bytes32(&b1));
            let y2 = x * x * x + G2Config::b();
            if let Some(y) = y2.sqrt() {
                return G2Affine { x, y, infinity: false };
            }
        }
    }

    #[test]
    fn batch_verify_accepts_valid_and_rejects_bad() {
        let sks: Vec<_> = (0..8)
            .map(|i| SecretKey::from_seed(format!("zoda batch key {}", i).as_bytes()))
            .collect();
        let msgs: Vec<Vec<u8>> = (0..8)
            .map(|i| format!("batch message {}", i).into_bytes())
            .collect();
        let items: Vec<(PublicKey, Vec<u8>, Signature)> = sks
            .iter()
            .zip(msgs.iter())
            .map(|(sk, m)| {
                let sig = sk.sign(m);
                (sk.public_key(), m.clone(), sig)
            })
            .collect();
        let refs: Vec<(&PublicKey, &[u8], &Signature)> =
            items.iter().map(|(p, m, s)| (p, m.as_slice(), s)).collect();
        assert!(verify_batch(&refs));
        assert!(verify_batch_strict(&refs));
        // one tampered signature breaks both
        let mut bad = items.clone();
        bad[3].2 = sks[0].sign(&msgs[0]); // wrong message's signature
        let refs: Vec<(&PublicKey, &[u8], &Signature)> =
            bad.iter().map(|(p, m, s)| (p, m.as_slice(), s)).collect();
        assert!(!verify_batch(&refs));
        assert!(!verify_batch_strict(&refs));
        // a wrong public key breaks both
        let mut bad = items.clone();
        bad[5].0 = sks[7].public_key();
        let refs: Vec<(&PublicKey, &[u8], &Signature)> =
            bad.iter().map(|(p, m, s)| (p, m.as_slice(), s)).collect();
        assert!(!verify_batch(&refs));
        assert!(!verify_batch_strict(&refs));
        // a non-G2 signature (on-curve, wrong subgroup) is caught by strict
        let mut bad = items.clone();
        bad[2].2 = random_bad_g2_point(1);
        let refs: Vec<(&PublicKey, &[u8], &Signature)> =
            bad.iter().map(|(p, m, s)| (p, m.as_slice(), s)).collect();
        assert!(!verify_batch_strict(&refs));
        // empty batch is trivially true
        assert!(verify_batch(&[]));
        assert!(verify_batch_strict(&[]));
    }

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
