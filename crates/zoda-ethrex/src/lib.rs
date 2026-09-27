//! # zoda-ethrex — ethrex-compatible execution-layer types
//!
//! Byte-for-byte compatible with ethrex's EIP-4844 types
//! (`BlobTransactionSidecar`, `KZGCommitment`, `KZGProof`,
//! `VersionedHash`) so zoda DA tooling can be wired directly into
//! ethrex's block import path. See `docs/ETHEX_INTEGRATION.md`.

pub use zoda_kzg::eip4844::{kzg_to_versioned_hash, VERSIONED_HASH_VERSION_KZG};
use zoda_kzg::srs::Setup;

pub const BLOB_VERSIONED_HASHES_INDEX: usize = 0;
pub const MAX_BLOBS_PER_BLOCK: usize = 16;

/// A KZG commitment (ethrex: 48-byte newtype).
pub type KzgCommitment = [u8; 48];
/// A KZG proof (ethrex: 48-byte newtype).
pub type KzgProof = [u8; 48];
/// A versioned hash (ethrex: H256).
pub type VersionedHash = [u8; 32];
/// A blob (ethrex: 131072 bytes).
pub type Blob = [u8; 131072];

/// ethrex `BlobTransactionSidecar`: blobs + commitments + proofs.
#[derive(Clone, Debug, PartialEq)]
pub struct BlobTransactionSidecar {
    pub blob_versioned_hashes: Vec<VersionedHash>,
    pub kzg_commitments: Vec<KzgCommitment>,
    pub blobs: Vec<Blob>,
    pub kzg_proofs: Vec<KzgProof>,
}

impl BlobTransactionSidecar {
    /// Validate every blob against its commitment via KZG (the exact
    /// check ethrex runs on blob transaction import: `verify_blob_kzg_proof`).
    pub fn validate(&self, setup: &Setup) -> bool {
        if self.blobs.len() != self.kzg_commitments.len()
            || self.blobs.len() != self.kzg_proofs.len()
            || self.blobs.len() != self.blob_versioned_hashes.len()
        {
            return false;
        }
        for i in 0..self.blobs.len() {
            // versioned hash consistency
            let vh = kzg_to_versioned_hash(&self.kzg_commitments[i]);
            if vh != self.blob_versioned_hashes[i] {
                return false;
            }
            if !zoda_kzg::eip4844::verify_blob_kzg_proof(
                &self.blobs[i],
                &self.kzg_commitments[i],
                &self.kzg_proofs[i],
                setup,
            )
            .unwrap_or(false)
            {
                return false;
            }
        }
        true
    }

    /// Batch validation (single pairing for the whole sidecar).
    pub fn validate_batch(&self, setup: &Setup) -> bool {
        let blobs: Vec<&[u8]> = self.blobs.iter().map(|b| b.as_slice()).collect();
        let comms: Vec<&[u8]> = self.kzg_commitments.iter().map(|c| c.as_slice()).collect();
        let proofs: Vec<&[u8]> = self.kzg_proofs.iter().map(|p| p.as_slice()).collect();
        zoda_kzg::eip4844::verify_blob_kzg_proof_batch(&blobs, &comms, &proofs, setup)
            .unwrap_or(false)
    }

    /// Build a sidecar from blobs (computing commitments + proofs).
    pub fn from_blobs(blobs: &[Blob], setup: &Setup) -> Result<Self, String> {
        let mut commitments = Vec::with_capacity(blobs.len());
        let mut proofs = Vec::with_capacity(blobs.len());
        let mut hashes = Vec::with_capacity(blobs.len());
        for blob in blobs {
            let c = zoda_kzg::eip4844::blob_to_kzg_commitment(blob, setup)?;
            let p = zoda_kzg::eip4844::compute_blob_kzg_proof(blob, &c, setup)?;
            hashes.push(kzg_to_versioned_hash(&c));
            commitments.push(c);
            proofs.push(p);
        }
        Ok(BlobTransactionSidecar {
            blob_versioned_hashes: hashes,
            kzg_commitments: commitments,
            blobs: blobs.to_vec(),
            kzg_proofs: proofs,
        })
    }

    /// Extend the sidecar into EIP-7594 cells + proofs for the DA layer.
    pub fn to_data_column_sidecars(
        &self,
        setup: &Setup,
    ) -> Result<Vec<Vec<zoda_edas::Cell>>, String> {
        let mut out = Vec::with_capacity(self.blobs.len());
        for blob in &self.blobs {
            let (cells, _proofs) = zoda_edas::compute_cells_and_kzg_proofs(blob, setup)?;
            out.push(cells);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_blob(seed: [u8; 32]) -> Blob {
        let mut rng = zoda_math::ZodaRng::from_seed(seed);
        let mut blob = [0u8; 131072];
        rng.next_bytes(&mut blob);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f;
        }
        blob
    }

    #[test]
    fn sidecar_roundtrip() {
        let setup = Setup::from_seed_for_testing(*b"ethrex-test-setup-00000000000000");
        let blobs = vec![
            random_blob(*b"ethrex-blob-1-seed-0000000000000"),
            random_blob(*b"ethrex-blob-2-seed-0000000000000"),
        ];
        let sidecar = BlobTransactionSidecar::from_blobs(&blobs, &setup).unwrap();
        assert!(sidecar.validate(&setup));
        assert!(sidecar.validate_batch(&setup));

        // tampered blob fails
        let mut bad = sidecar.clone();
        bad.blobs[0][5] ^= 0x01;
        assert!(!bad.validate(&setup));

        // mismatched versioned hash fails
        let mut bad2 = sidecar.clone();
        bad2.blob_versioned_hashes[0][5] ^= 0x01;
        assert!(!bad2.validate(&setup));
    }

    #[test]
    fn cells_extension() {
        let setup = Setup::from_seed_for_testing(*b"ethrex-cells-setup-0000000000000");
        let blobs = vec![random_blob(*b"ethrex-cell-blob-seed-0000000000")];
        let sidecar = BlobTransactionSidecar::from_blobs(&blobs, &setup).unwrap();
        let columns = sidecar.to_data_column_sidecars(&setup).unwrap();
        assert_eq!(columns.len(), 1);
        assert_eq!(columns[0].len(), 128);
    }
}
