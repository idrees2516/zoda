//! SHA-256 Merkle trees with compact audit paths — used for ZODA row/column
//! commitments and light-client bridge proofs.
//!
//! Domain separation: `leaf = H(0x00 ‖ data)`, `node = H(0x01 ‖ l ‖ r)`,
//! with duplicate-sibling padding for non-power-of-two leaf counts.

use crate::sha256::sha256;

fn hash_leaf(data: &[u8]) -> [u8; 32] {
    let mut buf = Vec::with_capacity(1 + data.len());
    buf.push(0x00);
    buf.extend_from_slice(data);
    sha256(&buf)
}

fn hash_node(l: &[u8; 32], r: &[u8; 32]) -> [u8; 32] {
    let mut buf = [0u8; 65];
    buf[0] = 0x01;
    buf[1..33].copy_from_slice(l);
    buf[33..65].copy_from_slice(r);
    sha256(&buf)
}

/// A Merkle tree over `n` leaves.
#[derive(Clone, Debug)]
pub struct MerkleTree {
    /// levels[0] = leaf hashes … levels.last = [root]
    levels: Vec<Vec<[u8; 32]>>,
}

impl MerkleTree {
    /// Build a Merkle tree over the given leaves.
    pub fn build(leaves: &[Vec<u8>]) -> Self {
        let leaf_hashes: Vec<[u8; 32]> = leaves.iter().map(|l| hash_leaf(l)).collect();
        Self::from_hashes(leaf_hashes)
    }

    /// Build from pre-hashed leaves.
    pub fn from_hashes(leaf_hashes: Vec<[u8; 32]>) -> Self {
        if leaf_hashes.is_empty() {
            return MerkleTree { levels: vec![] };
        }
        let mut levels = vec![leaf_hashes];
        while levels.last().unwrap().len() > 1 {
            let prev = levels.last().unwrap();
            let mut cur = Vec::with_capacity((prev.len() + 1) / 2);
            let mut i = 0;
            while i < prev.len() {
                let l = prev[i];
                let r = if i + 1 < prev.len() { prev[i + 1] } else { prev[i] };
                cur.push(hash_node(&l, &r));
                i += 2;
            }
            levels.push(cur);
        }
        MerkleTree { levels }
    }

    /// The Merkle root (`None` for an empty tree).
    pub fn root(&self) -> Option<[u8; 32]> {
        self.levels.last().and_then(|l| l.first().copied())
    }

    /// Generate the audit path for leaf `index`.
    pub fn proof(&self, index: usize) -> Option<MerkleProof> {
        if self.levels.is_empty() {
            return None;
        }
        let n = self.levels[0].len();
        if index >= n {
            return None;
        }
        let mut path = Vec::with_capacity(self.levels.len());
        let mut idx = index;
        for level in &self.levels[..self.levels.len() - 1] {
            let sibling = if idx % 2 == 0 {
                // right sibling (or self if duplicated)
                if idx + 1 < level.len() {
                    (level[idx + 1], false)
                } else {
                    (level[idx], true) // duplicate of self
                }
            } else {
                (level[idx - 1], true) // left sibling
            };
            path.push(sibling);
            idx /= 2;
        }
        Some(MerkleProof { path })
    }

    /// Number of leaves.
    pub fn leaf_count(&self) -> usize {
        self.levels.first().map(|l| l.len()).unwrap_or(0)
    }
}

/// An audit path: sequence of (sibling hash, sibling_is_left).
#[derive(Clone, Debug, PartialEq)]
pub struct MerkleProof {
    pub path: Vec<([u8; 32], bool)>,
}

impl MerkleProof {
    /// Verify this proof for `leaf_hash` against `root`.
    pub fn verify(&self, leaf_hash: &[u8; 32], root: &[u8; 32]) -> bool {
        let mut cur = *leaf_hash;
        for (sibling, is_left) in &self.path {
            cur = if *is_left {
                hash_node(sibling, &cur)
            } else {
                hash_node(&cur, sibling)
            };
        }
        cur == *root
    }
}

/// Convenience: verify with raw leaf data.
pub fn verify_leaf(root: &[u8; 32], data: &[u8], proof: &MerkleProof) -> bool {
    proof.verify(&hash_leaf(data), root)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merkle_roundtrip() {
        let leaves: Vec<Vec<u8>> = (0..13u32)
            .map(|i| i.to_le_bytes().to_vec())
            .collect();
        let t = MerkleTree::build(&leaves);
        let root = t.root().unwrap();
        for i in 0..leaves.len() {
            let p = t.proof(i).unwrap();
            assert!(verify_leaf(&root, &leaves[i], &p), "proof failed i={}", i);
        }
        // wrong leaf fails
        let p = t.proof(0).unwrap();
        assert!(!verify_leaf(&root, b"imposter", &p));
    }
}
