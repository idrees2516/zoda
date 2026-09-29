//! # zoda-kzg — KZG polynomial commitments for Ethereum (EIP-4844)
//!
//! Evaluation-form (Lagrange) KZG exactly as specified by the Deneb
//! consensus specification, plus the byte-level EIP-4844 API used by
//! Ethereum execution and consensus clients.
//!
//! * Trusted setup parsing (both the 2-section and current 3-section
//!   c-kzg-4844 `trusted_setup.txt` formats) with the spec's
//!   Lagrange-form sanity check.
//! * `blob_to_kzg_commitment`, `compute_kzg_proof`, `verify_kzg_proof`,
//!   `verify_kzg_proof_batch`, `compute_blob_kzg_proof`,
//!   `verify_blob_kzg_proof(_batch)` — spec-faithful.
//! * Pippenger multi-scalar multiplication for the MSM hot paths.

pub mod eip4844;
pub mod g1fft;
pub mod msm;
pub mod srs;
pub mod spectest;

pub use eip4844::*;
pub use srs::{Setup, BLOB_BYTES, FIELD_ELEMENTS_PER_BLOB};

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_math::{Fr, ZodaRng};

    fn test_setup() -> Setup {
        // Build a small deterministic setup for tests: [tau^i] with tau
        // derived from a seed. NOT a real ceremony — unit tests only.
        Setup::from_seed_for_testing(*b"zoda-kzg-test-setup-tau-seed-000")
    }

    #[test]
    fn blob_commit_roundtrip() {
        let setup = test_setup();
        let mut rng = ZodaRng::from_seed(*b"kzg-test-seed-000000000000000000");
        let mut blob = [0u8; BLOB_BYTES];
        rng.next_bytes(&mut blob);
        // keep field elements canonical (big-endian < r): mask top byte
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f;
        }
        let c = blob_to_kzg_commitment(&blob, &setup).unwrap();
        // commitment is a valid compressed G1 point
        assert_eq!(c.len(), 48);
        assert!(zoda_bls::G1Affine::from_compressed(&c).is_some());
    }

    #[test]
    fn proof_verify_roundtrip() {
        let setup = test_setup();
        let mut rng = ZodaRng::from_seed(*b"kzg-proof-seed-00000000000000000");
        let mut blob = [0u8; BLOB_BYTES];
        rng.next_bytes(&mut blob);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f;
        }
        let commitment = blob_to_kzg_commitment(&blob, &setup).unwrap();
        let mut z_bytes = [0u8; 32];
        rng.next_bytes(&mut z_bytes);
        z_bytes[0] &= 0x3f;
        let (proof, y) = compute_kzg_proof(&blob, &z_bytes, &setup).unwrap();
        assert!(verify_kzg_proof(&commitment, &z_bytes, &y, &proof, &setup).unwrap());

        // wrong y fails
        let mut bad_y = y;
        bad_y[31] ^= 0x01;
        assert!(!verify_kzg_proof(&commitment, &z_bytes, &bad_y, &proof, &setup).unwrap());

        // wrong z (different evaluation point) fails, and a fresh proof
        // at the new point verifies
        let mut bad_z = z_bytes;
        bad_z[31] ^= 0x01;
        assert!(!verify_kzg_proof(&commitment, &bad_z, &y, &proof, &setup).unwrap());
        let (proof2, y2) = compute_kzg_proof(&blob, &bad_z, &setup).unwrap();
        assert!(verify_kzg_proof(&commitment, &bad_z, &y2, &proof2, &setup).unwrap());
    }

    #[test]
    fn blob_proof_roundtrip() {
        let setup = test_setup();
        let mut rng = ZodaRng::from_seed(*b"kzg-blob-seed-000000000000000000");
        let mut blob = [0u8; BLOB_BYTES];
        rng.next_bytes(&mut blob);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f;
        }
        let commitment = blob_to_kzg_commitment(&blob, &setup).unwrap();
        let proof = compute_blob_kzg_proof(&blob, &commitment, &setup).unwrap();
        assert!(verify_blob_kzg_proof(&blob, &commitment, &proof, &setup).unwrap());
        // tampered blob fails
        let mut bad = blob;
        bad[0] ^= 0x01;
        assert!(!verify_blob_kzg_proof(&bad, &commitment, &proof, &setup).unwrap());
    }

    #[test]
    fn batch_verification() {
        let setup = test_setup();
        let mut rng = ZodaRng::from_seed(*b"kzg-batch-seed-00000000000000000");
        let mut blobs = Vec::new();
        for _ in 0..4 {
            let mut blob = [0u8; BLOB_BYTES];
            rng.next_bytes(&mut blob);
            for chunk in blob.chunks_exact_mut(32) {
                chunk[0] &= 0x3f;
            }
            blobs.push(blob);
        }
        let commitments: Vec<_> = blobs
            .iter()
            .map(|b| blob_to_kzg_commitment(b, &setup).unwrap())
            .collect();
        let proofs: Vec<_> = blobs
            .iter()
            .zip(&commitments)
            .map(|(b, c)| compute_blob_kzg_proof(b, c, &setup).unwrap())
            .collect();
        let blob_refs: Vec<&[u8]> = blobs.iter().map(|b| b.as_slice()).collect();
        let comm_refs: Vec<&[u8]> = commitments.iter().map(|c| c.as_slice()).collect();
        let proof_refs: Vec<&[u8]> = proofs.iter().map(|p| p.as_slice()).collect();
        assert!(verify_blob_kzg_proof_batch(&blob_refs, &comm_refs, &proof_refs, &setup).unwrap());
        // tamper one proof: the tampered point may be invalid (Err) or
        // simply fail verification — both are acceptable rejections
        let mut bad_proofs = proofs.clone();
        bad_proofs[2][5] ^= 0x01;
        let bad_refs: Vec<&[u8]> = bad_proofs.iter().map(|p| p.as_slice()).collect();
        let r = verify_blob_kzg_proof_batch(&blob_refs, &comm_refs, &bad_refs, &setup);
        assert!(r.is_err() || !r.unwrap());
    }

    #[test]
    fn invalid_blob_elements_rejected() {
        let setup = test_setup();
        let mut blob = [0u8; BLOB_BYTES];
        blob[0] = 0xff; // >= r
        assert!(blob_to_kzg_commitment(&blob, &setup).is_err());
    }
}
