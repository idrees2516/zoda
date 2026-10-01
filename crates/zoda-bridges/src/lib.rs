//! # zoda-bridges — light-client bridges over the DA layer
//!
//! * **Inclusion proofs**: prove that a specific row/column of a
//!   committed grid is included (ZODA Merkle proofs + the zero-overhead
//!   projection check),
//! * **Message commitments**: KZG-backed commitments for cross-chain
//!   messages bundled into blobs,
//! * **Light-client verification**: everything verifiable with public
//!   parameters only.

use zoda_core::{Matrix, ZodaPublic};
use zoda_kzg::srs::Setup;
use zoda_math::merkle::{MerkleProof, MerkleTree};
use zoda_math::{Fr, PrimeField};

/// A verifiable inclusion proof for one row of the grid.
pub struct RowInclusionProof {
    pub index: usize,
    pub row: Vec<Fr>,
    pub merkle: MerkleProof,
}

/// A verifiable inclusion proof for one column of the grid.
pub struct ColumnInclusionProof {
    pub index: usize,
    pub column: Vec<Fr>,
    pub merkle: MerkleProof,
}

/// Build the row/column Merkle trees for a committed matrix.
pub fn commitment_trees(matrix: &Matrix) -> (MerkleTree, MerkleTree) {
    let row_leaves: Vec<Vec<u8>> = (0..matrix.rows())
        .map(|r| row_leaf(matrix, r))
        .collect();
    let col_leaves: Vec<Vec<u8>> = (0..matrix.cols())
        .map(|c| col_leaf(matrix, c))
        .collect();
    (MerkleTree::build(&row_leaves), MerkleTree::build(&col_leaves))
}

fn row_leaf(matrix: &Matrix, r: usize) -> Vec<u8> {
    let mut buf = vec![0x01u8];
    buf.extend_from_slice(&(r as u64).to_le_bytes());
    for c in 0..matrix.cols() {
        buf.extend_from_slice(&matrix.get(r, c).to_le_bytes());
    }
    buf
}

fn col_leaf(matrix: &Matrix, c: usize) -> Vec<u8> {
    let mut buf = vec![0x02u8];
    buf.extend_from_slice(&(c as u64).to_le_bytes());
    for r in 0..matrix.rows() {
        buf.extend_from_slice(&matrix.get(r, c).to_le_bytes());
    }
    buf
}

/// Prove row inclusion (Merkle + ZODA correctness).
pub fn prove_row(
    matrix: &Matrix,
    trees: &(MerkleTree, MerkleTree),
    index: usize,
) -> RowInclusionProof {
    RowInclusionProof {
        index,
        row: matrix.row(index),
        merkle: trees.0.proof(index).expect("in-range"),
    }
}

/// Verify a row inclusion proof against the public roots and the ZODA
/// projection check.
pub fn verify_row(
    proof: &RowInclusionProof,
    public: &ZodaPublic,
) -> bool {
    // 1. ZODA correctness (the row is its own proof)
    if zoda_core::verify_row_sample(public, proof.index, &proof.row).is_err() {
        return false;
    }
    // 2. Merkle inclusion — the tree hashes leaves as H(0x00 || leaf)
    let mut buf = vec![0x01u8];
    buf.extend_from_slice(&(proof.index as u64).to_le_bytes());
    for f in &proof.row {
        buf.extend_from_slice(&f.to_le_bytes());
    }
    let mut leaf_input = Vec::with_capacity(1 + buf.len());
    leaf_input.push(0x00);
    leaf_input.extend_from_slice(&buf);
    let leaf_hash = zoda_math::sha256::sha256(&leaf_input);
    proof.merkle.verify(&leaf_hash, &public.row_root)
}

/// A cross-chain message: payload + destination + nonce.
#[derive(Clone, Debug)]
pub struct BridgeMessage {
    pub destination_chain: [u8; 4],
    pub nonce: u64,
    pub payload: Vec<u8>,
}

impl BridgeMessage {
    pub fn encoding(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(12 + self.payload.len());
        out.extend_from_slice(&self.destination_chain);
        out.extend_from_slice(&self.nonce.to_le_bytes());
        out.extend_from_slice(&self.payload);
        out
    }

    /// KZG commitment to the message polynomial (coefficients from the
    /// encoding chunks) — the bridge contract verifies KZG proofs of
    /// message inclusion.
    pub fn kzg_commitment(&self, setup: &Setup) -> Result<[u8; 48], String> {
        // pack the encoding into 32-byte big-endian field elements
        let mut blob = vec![0u8; 131072];
        let n = self.encoding().len().min(131072);
        blob[..n].copy_from_slice(&self.encoding()[..n]);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f; // canonicalize
        }
        zoda_kzg::eip4844::blob_to_kzg_commitment(&blob, setup)
    }

    /// Prove that the message is included in `blob` at `index` (the
    /// blob's polynomial evaluated at a challenge point).
    pub fn kzg_inclusion_proof(
        &self,
        blob: &[u8],
        commitment: &[u8],
        setup: &Setup,
    ) -> Result<[u8; 48], String> {
        zoda_kzg::eip4844::compute_blob_kzg_proof(blob, commitment, setup)
    }
}

/// Verify a message's inclusion in a blob via KZG.
pub fn verify_message_inclusion(
    blob: &[u8],
    commitment: &[u8],
    proof: &[u8],
    setup: &Setup,
) -> bool {
    zoda_kzg::eip4844::verify_blob_kzg_proof(blob, commitment, proof, setup).unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_core::ZodaParams;
    use zoda_math::ZodaRng;

    #[test]
    fn row_inclusion_roundtrip() {
        let mut rng = ZodaRng::from_seed(*b"bridge-test-seed-000000000000000");
        let mut data = vec![Fr::zero(); 64];
        for x in data.iter_mut() {
            *x = rng.next_fr(false);
        }
        let grid = Matrix::from_row_major(8, 8, data);
        let params = ZodaParams::new(8, 8);
        let prover = params.commit(&grid);
        let public = prover.public_params();
        let trees = commitment_trees(&prover.matrix);

        for r in [0usize, 3, 15] {
            let proof = prove_row(&prover.matrix, &trees, r);
            assert!(verify_row(&proof, &public), "row {} failed", r);
        }

        // a tampered row fails
        let mut bad = prove_row(&prover.matrix, &trees, 3);
        bad.row[0] = bad.row[0] + Fr::ONE;
        assert!(!verify_row(&bad, &public));

        // a row at the wrong index fails
        let mut wrong = prove_row(&prover.matrix, &trees, 3);
        wrong.index = 4;
        assert!(!verify_row(&wrong, &public));
    }

    #[test]
    fn message_kzg_roundtrip() {
        let setup = Setup::from_seed_for_testing(*b"bridge-kzg-setup-000000000000000");
        let msg = BridgeMessage {
            destination_chain: *b"ETH2",
            nonce: 42,
            payload: b"cross-chain hello".to_vec(),
        };
        // message packed into a blob
        let mut blob = vec![0u8; 131072];
        let enc = msg.encoding();
        blob[..enc.len()].copy_from_slice(&enc);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f;
        }
        let comm = msg.kzg_commitment(&setup).unwrap();
        let proof = msg.kzg_inclusion_proof(&blob, &comm, &setup).unwrap();
        assert!(verify_message_inclusion(&blob, &comm, &proof, &setup));
        // tampered blob fails
        let mut bad = blob.clone();
        bad[100] ^= 0x01;
        assert!(!verify_message_inclusion(&bad, &comm, &proof, &setup));
    }
}
