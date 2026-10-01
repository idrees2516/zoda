//! EIP-7594 column-cell custody: per-blob cells and KZG proofs organized
//! by column, with verification on ingest and a compact disk format.

use zoda_edas::{verify_cell_kzg_proof_batch, Cell, NUMBER_OF_COLUMNS};
use zoda_kzg::srs::Setup;
use zoda_math::sha256::sha256;

const CELL_MAGIC: &[u8; 8] = b"ZODACEL1";
const CELL_VERSION: u32 = 1;
const PROOF_BYTES: usize = 48;

/// One column's custody: for every blob in the slot, the cell at this
/// column and its KZG proof.
pub struct CellCustody {
    pub slot: u64,
    /// per-blob KZG commitments (48 bytes each)
    pub commitments: Vec<Vec<u8>>,
    /// `[column]` → one `(cell, proof)` per blob, or None if not held
    columns: Vec<Option<Vec<(Cell, [u8; PROOF_BYTES])>>>,
}

impl CellCustody {
    pub fn new(slot: u64, blobs: usize) -> CellCustody {
        CellCustody {
            slot,
            commitments: vec![Vec::new(); blobs],
            columns: vec![None; NUMBER_OF_COLUMNS],
        }
    }

    pub fn blobs(&self) -> usize {
        self.commitments.len()
    }

    /// Set the blob commitments (needed for verification).
    pub fn set_commitments(&mut self, commitments: Vec<Vec<u8>>) {
        assert_eq!(commitments.len(), self.commitments.len());
        self.commitments = commitments;
    }

    /// Store one column after verifying every cell against the blob
    /// commitments with a single pairing check.
    pub fn put_column_verified(
        &mut self,
        column: usize,
        items: &[(Cell, [u8; PROOF_BYTES])],
        setup: &Setup,
    ) -> Result<bool, String> {
        if column >= NUMBER_OF_COLUMNS {
            return Err("column index out of range".to_string());
        }
        if items.len() != self.blobs() {
            return Err("one (cell, proof) per blob required".to_string());
        }
        let commitments: Vec<&[u8]> = self
            .commitments
            .iter()
            .map(|c| c.as_slice())
            .collect();
        if commitments.iter().any(|c| c.len() != 48) {
            return Err("blob commitments not set".to_string());
        }
        let indices = vec![column as u64; items.len()];
        let cells: Vec<Cell> = items.iter().map(|(c, _)| *c).collect();
        let proofs: Vec<&[u8]> = items.iter().map(|(_, p)| p.as_slice()).collect();
        if !verify_cell_kzg_proof_batch(&commitments, &indices, &cells, &proofs, setup)? {
            return Err("cell proof batch failed".to_string());
        }
        let fresh = self.columns[column].is_none();
        self.columns[column] = Some(items.to_vec());
        Ok(fresh)
    }

    /// Store without verification (trusted local data only).
    pub fn put_column_unchecked(
        &mut self,
        column: usize,
        items: Vec<(Cell, [u8; PROOF_BYTES])>,
    ) -> Result<bool, String> {
        if column >= NUMBER_OF_COLUMNS {
            return Err("column index out of range".to_string());
        }
        if items.len() != self.blobs() {
            return Err("one (cell, proof) per blob required".to_string());
        }
        let fresh = self.columns[column].is_none();
        self.columns[column] = Some(items);
        Ok(fresh)
    }

    pub fn get_column(&self, column: usize) -> Option<&[(Cell, [u8; PROOF_BYTES])]> {
        self.columns[column].as_deref()
    }

    pub fn available_columns(&self) -> Vec<usize> {
        (0..NUMBER_OF_COLUMNS)
            .filter(|&c| self.columns[c].is_some())
            .collect()
    }

    pub fn column_count(&self) -> usize {
        self.columns.iter().filter(|c| c.is_some()).count()
    }

    /// Custody compliance for an assigned column set.
    pub fn custody_compliance(&self, assigned: &[usize]) -> (usize, usize) {
        (
            assigned.len(),
            assigned
                .iter()
                .filter(|&&c| c < NUMBER_OF_COLUMNS && self.columns[c].is_some())
                .count(),
        )
    }

    /// Exact serialized size.
    pub fn byte_size(&self) -> usize {
        let blobs = self.blobs();
        8 + 4 + 8 + 4
            + blobs * 48
            + NUMBER_OF_COLUMNS / 8
            + self.column_count() * blobs * (BYTES_PER_CELL + PROOF_BYTES)
            + 32
    }

    fn bitfield(&self) -> [u8; NUMBER_OF_COLUMNS / 8] {
        let mut bits = [0u8; NUMBER_OF_COLUMNS / 8];
        for (c, col) in self.columns.iter().enumerate() {
            if col.is_some() {
                bits[c / 8] |= 1 << (c % 8);
            }
        }
        bits
    }

    /// Serialize to the compact `ZODACEL1` format.
    pub fn encode(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(self.byte_size());
        buf.extend_from_slice(CELL_MAGIC);
        buf.extend_from_slice(&CELL_VERSION.to_le_bytes());
        buf.extend_from_slice(&self.slot.to_le_bytes());
        buf.extend_from_slice(&(self.blobs() as u32).to_le_bytes());
        for c in &self.commitments {
            buf.extend_from_slice(c);
        }
        buf.extend_from_slice(&self.bitfield());
        for col in self.columns.iter().flatten() {
            for (cell, proof) in col {
                buf.extend_from_slice(cell);
                buf.extend_from_slice(proof);
            }
        }
        let checksum = sha256(&buf);
        buf.extend_from_slice(&checksum);
        buf
    }

    /// Parse + checksum-verify a `ZODACEL1` payload.
    pub fn decode(bytes: &[u8]) -> Result<CellCustody, String> {
        if bytes.len() < 8 + 4 + 8 + 4 + 32 {
            return Err("file too small".into());
        }
        let split = bytes.len() - 32;
        if sha256(&bytes[..split])[..] != bytes[split..] {
            return Err("checksum mismatch".into());
        }
        let mut cur: &[u8] = &bytes[..split];
        fn take<'a>(cur: &mut &'a [u8], n: usize) -> Result<&'a [u8], String> {
            if cur.len() < n {
                return Err("truncated".into());
            }
            let out = &cur[..n];
            *cur = &cur[n..];
            Ok(out)
        }
        if take(&mut cur, 8)? != CELL_MAGIC {
            return Err("not a cell custody file".into());
        }
        let version = u32::from_le_bytes(take(&mut cur, 4)?.try_into().unwrap());
        if version != CELL_VERSION {
            return Err("unsupported version".into());
        }
        let slot = u64::from_le_bytes(take(&mut cur, 8)?.try_into().unwrap());
        let blobs =
            u32::from_le_bytes(take(&mut cur, 4)?.try_into().unwrap()) as usize;
        if blobs == 0 || blobs > 4096 {
            return Err("implausible blob count".into());
        }
        let mut commitments = Vec::with_capacity(blobs);
        for _ in 0..blobs {
            commitments.push(take(&mut cur, 48)?.to_vec());
        }
        let bits = take(&mut cur, NUMBER_OF_COLUMNS / 8)?.to_vec();
        let present: usize = (0..NUMBER_OF_COLUMNS)
            .filter(|&c| bits[c / 8] & (1 << (c % 8)) != 0)
            .count();
        let item = BYTES_PER_CELL + PROOF_BYTES;
        if cur.len() != present * blobs * item {
            return Err("cell section size mismatch".into());
        }
        let mut columns: Vec<Option<Vec<(Cell, [u8; PROOF_BYTES])>>> =
            vec![None; NUMBER_OF_COLUMNS];
        for c in 0..NUMBER_OF_COLUMNS {
            if bits[c / 8] & (1 << (c % 8)) != 0 {
                let mut items = Vec::with_capacity(blobs);
                for _ in 0..blobs {
                    let cell_bytes = take(&mut cur, BYTES_PER_CELL)?;
                    let mut cell: Cell = [0u8; BYTES_PER_CELL];
                    cell.copy_from_slice(cell_bytes);
                    let mut proof = [0u8; PROOF_BYTES];
                    proof.copy_from_slice(take(&mut cur, PROOF_BYTES)?);
                    items.push((cell, proof));
                }
                columns[c] = Some(items);
            }
        }
        Ok(CellCustody {
            slot,
            commitments,
            columns,
        })
    }
}

const BYTES_PER_CELL: usize = 2048;

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_edas::compute_cells_and_kzg_proofs;
    use zoda_kzg::srs::Setup;
    use zoda_math::ZodaRng;

    #[test]
    fn format_roundtrip_and_integrity() {
        let blobs = 2;
        let mut custody = CellCustody::new(7, blobs);
        // fake cells/proofs (format-level test; verification tested below)
        let mk = |seed: u8| -> (Cell, [u8; 48]) {
            let mut cell = [0u8; BYTES_PER_CELL];
            for (i, b) in cell.iter_mut().enumerate() {
                *b = seed.wrapping_add(i as u8);
            }
            let mut proof = [seed; 48];
            proof[0] = 1;
            (cell, proof)
        };
        custody.set_commitments(vec![vec![1u8; 48], vec![2u8; 48]]);
        custody.put_column_unchecked(0, vec![mk(3), mk(4)]).unwrap();
        custody.put_column_unchecked(129.min(NUMBER_OF_COLUMNS - 1), vec![mk(5), mk(6)])
            .unwrap();
        assert_eq!(custody.byte_size(), {
            let mut b = custody.encode();
            let len = b.len();
            b.truncate(len);
            len
        });
        let bytes = custody.encode();
        let back = CellCustody::decode(&bytes).expect("decode");
        assert_eq!(back.slot, 7);
        assert_eq!(back.available_columns(), custody.available_columns());
        assert_eq!(back.get_column(0), custody.get_column(0));
        // corruption caught
        let mut bad = bytes.clone();
        let n = bad.len();
        bad[n - 100] ^= 0xa5;
        assert!(CellCustody::decode(&bad).is_err());
    }

    #[test]
    fn verified_ingest() {
        let setup = Setup::from_seed_for_testing(*b"arch-cells-toy-tau-0000000000000");
        let mut blob = [0u8; 131072];
        let mut rng = ZodaRng::from_seed(*b"arch-cells-blob-seed-00000000000");
        rng.next_bytes(&mut blob);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f; // keep field elements canonical
        }
        let (cells, proofs) = compute_cells_and_kzg_proofs(&blob, &setup).unwrap();
        let poly = zoda_edas::blob_to_monomial(&blob, &setup).unwrap();
        let commitment = setup.commit_coeffs(&poly).to_compressed();
        let mut custody = CellCustody::new(1, 1);
        custody.set_commitments(vec![commitment]);
        // honest column verifies and stores
        let col = 5usize;
        assert!(
            custody
                .put_column_verified(col, &[(cells[col], proofs[col])], &setup)
                .is_ok()
        );
        // tampered cell is rejected
        let mut tampered = cells[col];
        tampered[10] ^= 0xff;
        assert!(
            custody
                .put_column_verified(6, &[(tampered, proofs[6])], &setup)
                .is_err()
        );
        assert!(custody.get_column(6).is_none());
    }
}
