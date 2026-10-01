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
use crate::hash::{hash_to_curve_g2, hash_to_curve_g2_batch};
use crate::pairing::{final_exponentiation, multi_miller_loop, pairing_check, G2Prepared};
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

/// Worker count for a batch verification: the available parallelism,
/// engaged only when the batch is large enough that per-thread setup is
/// amortised (≥ 8 items), capped by the item count and by 8.
fn batch_threads(n: usize) -> usize {
    if n < 8 {
        return 1;
    }
    let cores = std::thread::available_parallelism()
        .map(|v| v.get())
        .unwrap_or(1);
    cores.min(n).min(8)
}

/// The per-chunk workhorse of [`verify_batch`]: batched hash-to-curve,
/// the plain scalar multiplications (pk/σ are untrusted — the h-torsion
/// semantics of single verification must be preserved exactly), G2
/// preparation of the hashed points, and the chunk's share of the
/// Miller-loop accumulation. Returns the chunk's Miller product and its
/// partial signature sum.
#[allow(clippy::type_complexity)]
fn verify_batch_chunk(
    items: &[(&PublicKey, &[u8], &Signature)],
    challenges: &[Fr],
    strict: bool,
) -> (crate::fp12::Fp12, G2Projective) {
    let msgs: Vec<&[u8]> = items.iter().map(|(_, m, _)| *m).collect();
    let hs = hash_to_curve_g2_batch(&msgs, DST);
    let mut terms: Vec<(G1Affine, G2Prepared)> = Vec::with_capacity(items.len());
    let mut sig_sum = G2Projective::identity();
    for ((pk, _, sig), r) in items.iter().zip(challenges.iter()) {
        let rp = if strict {
            crate::endomorphism::mul_g1_public(&pk.to_projective(), r)
        } else {
            pk.to_projective().mul_limbs(&r.to_repr())
        }
        .to_affine();
        terms.push((rp, G2Prepared::from(hs[terms.len()].to_affine())));
        let rs = if strict {
            crate::endomorphism::mul_g2_public(&sig.to_projective(), r)
        } else {
            sig.to_projective().mul_limbs(&r.to_repr())
        };
        sig_sum = sig_sum + rs;
    }
    let term_refs: Vec<(&G1Affine, &G2Prepared)> = terms.iter().map(|(a, b)| (a, b)).collect();
    (multi_miller_loop(&term_refs), sig_sum)
}

/// Fold the accumulated signature sum into the pairing product with one
/// extra serial Miller term and the single final exponentiation.
fn finish_batch_pairing(
    f: crate::fp12::Fp12,
    sig_sum: &G2Projective,
) -> bool {
    let neg_g1 = G1Affine::generator().neg();
    let sig_prepared = G2Prepared::from(sig_sum.to_affine());
    let f = f * multi_miller_loop(&[(&neg_g1, &sig_prepared)]);
    match final_exponentiation(&f) {
        Some(v) => v.is_one(),
        None => false,
    }
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
///
/// v1.4: the per-item work — batched hash-to-curve (SSWU square roots
/// via the complex method with Montgomery batch inversion), the scalar
/// multiplications, G2 preparation and the Miller-loop accumulation — is
/// chunked across worker threads; the threads' Fp12 products combine
/// with `T−1` multiplications and the `Σ rᵢσᵢ` term folds in before one
/// serial final exponentiation.
pub fn verify_batch(items: &[(&PublicKey, &[u8], &Signature)]) -> bool {
    if items.is_empty() {
        return true;
    }
    let challenges = batch_challenges(items);
    for (pk, _, sig) in items.iter() {
        if pk.infinity || sig.infinity {
            return false;
        }
        if !pk.is_on_curve() || !sig.is_on_curve() {
            return false;
        }
    }
    let t = batch_threads(items.len());
    if t <= 1 {
        let (f, sig_sum) = verify_batch_chunk(items, &challenges.0, false);
        return finish_batch_pairing(f, &sig_sum);
    }
    let chunk = items.len().div_ceil(t);
    let mut f = crate::fp12::Fp12::one();
    let mut sig_sum = G2Projective::identity();
    std::thread::scope(|scope| {
        let handles: Vec<_> = items
            .chunks(chunk)
            .zip(challenges.0.chunks(chunk))
            .map(|(ic, rc)| scope.spawn(move || verify_batch_chunk(ic, rc, false)))
            .collect();
        for h in handles {
            match h.join() {
                Ok((mf, ss)) => {
                    f = f * mf;
                    sig_sum = sig_sum + ss;
                }
                Err(_) => std::process::abort(),
            }
        }
    });
    finish_batch_pairing(f, &sig_sum)
}

/// Batch verification with full subgroup soundness: one Pippenger-based G1
/// batch subgroup check over the public keys, one combined G2 check over
/// the signatures, then the same product-pairing check as [`verify_batch`]
/// — with GLV-accelerated scalar multiplications, since every point has
/// been subgroup-verified first — parallelised across chunks exactly as
/// [`verify_batch`].
pub fn verify_batch_strict(items: &[(&PublicKey, &[u8], &Signature)]) -> bool {
    if items.is_empty() {
        return true;
    }
    let pks: Vec<PublicKey> = items.iter().map(|(pk, _, _)| **pk).collect();
    let sigs: Vec<Signature> = items.iter().map(|(_, _, sig)| **sig).collect();
    for (pk, _, sig) in items.iter() {
        if pk.infinity || sig.infinity {
            return false;
        }
        if !pk.is_on_curve() || !sig.is_on_curve() {
            return false;
        }
    }
    if !crate::endomorphism::batch_subgroup_check_g2(&sigs) {
        return false;
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
    let t = batch_threads(items.len());
    if t <= 1 {
        let (f, sig_sum) = verify_batch_chunk(items, &challenges.0, true);
        return finish_batch_pairing(f, &sig_sum);
    }
    let chunk = items.len().div_ceil(t);
    let mut f = crate::fp12::Fp12::one();
    let mut sig_sum = G2Projective::identity();
    std::thread::scope(|scope| {
        let handles: Vec<_> = items
            .chunks(chunk)
            .zip(challenges.0.chunks(chunk))
            .map(|(ic, rc)| scope.spawn(move || verify_batch_chunk(ic, rc, true)))
            .collect();
        for h in handles {
            match h.join() {
                Ok((mf, ss)) => {
                    f = f * mf;
                    sig_sum = sig_sum + ss;
                }
                Err(_) => std::process::abort(),
            }
        }
    });
    finish_batch_pairing(f, &sig_sum)
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
    fn batch_verify_parallel_path() {
        // 17 items: exercises the chunked threading (uneven chunks with
        // 2 workers) and both the plain and strict paths, plus tamper
        // rejection through the parallel combination
        let items: Vec<(PublicKey, Vec<u8>, Signature)> = (0..17usize)
            .map(|i| {
                let sk = SecretKey::from_seed(format!("par batch key {}", i).as_bytes());
                let msg = format!("par batch message {}", i).into_bytes();
                let sig = sk.sign(&msg);
                (sk.public_key(), msg, sig)
            })
            .collect();
        let refs: Vec<(&PublicKey, &[u8], &Signature)> =
            items.iter().map(|(p, m, s)| (p, m.as_slice(), s)).collect();
        assert!(verify_batch(&refs));
        assert!(verify_batch_strict(&refs));
        let mut bad = items.clone();
        bad[9].1 = b"tampered message".to_vec();
        let refs: Vec<(&PublicKey, &[u8], &Signature)> =
            bad.iter().map(|(p, m, s)| (p, m.as_slice(), s)).collect();
        assert!(!verify_batch(&refs));
        assert!(!verify_batch_strict(&refs));
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
