//! # zoda-archival — custody storage, reconstruction and pruning
//!
//! An archival / custodial node stores a subset of a grid's columns (its
//! **custody duty**), serves them to samplers, triggers reconstruction
//! once ≥ k of 2k columns are verified, and manages its storage budget.
//!
//! * [`GridCustody`] — one grid's column set: verify-on-insert, exact byte
//!   accounting, reconstruction that PERSISTS the recovered columns (each
//!   re-verified before acceptance), and custody-duty compliance tracking.
//! * [`CustodyStore`] — many grids under a byte budget with slot-ordered
//!   (LRU) pruning.
//! * [`DiskCustody`] — a directory of self-contained, checksummed custody
//!   files written atomically (temp + rename) and re-verified on load.
//! * [`CellCustody`] — the EIP-7594 equivalent: per-blob cell columns
//!   (cell + KZG proof) with the same compact-format discipline.
//!
//! ## Wire/disk format
//!
//! Custody files are compact and self-verifying: a fixed header (grid id,
//! shape, the full public verification parameters), an availability
//! bitfield, the raw 32-byte field elements of exactly the present
//! columns, and a SHA-256 checksum over everything preceding. No
//! framing, no hex, no redundant per-column metadata.

pub mod cells;
pub mod disk;
pub mod format;
pub mod store;

pub use cells::CellCustody;
pub use disk::DiskCustody;
pub use store::{CustodyStore, GridCustody, ReconOutcome};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::store::CustodyCompliance;
    use zoda_core::{Matrix, ZodaParams};
    use zoda_math::{Fr, PrimeField, ZodaRng};

    fn grid(m: usize, k: usize, seed: &[u8]) -> (ZodaParams, zoda_core::ZodaProver, zoda_core::ZodaPublic) {
        let mut s32 = [0u8; 32];
        s32[..seed.len().min(32)].copy_from_slice(&seed[..seed.len().min(32)]);
        let mut rng = ZodaRng::from_seed(s32);
        let mut data = vec![Fr::zero(); m * k];
        for x in data.iter_mut() {
            *x = rng.next_fr(false);
        }
        let g = Matrix::from_row_major(m, k, data);
        let params = ZodaParams::new(m, k);
        let prover = params.commit(&g);
        let public = prover.public_params();
        (params, prover, public)
    }

    #[test]
    fn store_reconstruct_and_persist() {
        let (params, prover, public) = grid(8, 8, b"arch-store-1");
        let mut custody = GridCustody::new([7u8; 32], 42, params, public);
        // store only k = 8 of 16 columns (verified)
        for c in 0..8 {
            assert!(custody.put_column_verified(c, prover.matrix.col(c)).unwrap());
        }
        assert_eq!(custody.available_columns().len(), 8);
        assert!(!custody.is_complete());
        // reconstruct — MUST persist the recovered columns
        let outcome = custody.try_reconstruct().expect("reconstruct");
        assert_eq!(outcome.recovered_columns, 8);
        assert!(custody.is_complete());
        assert_eq!(custody.available_columns().len(), 16);
        // recovered columns verify and match the original
        for c in 8..16 {
            let col = custody.get_column(c).unwrap();
            assert_eq!(col, &prover.matrix.col(c));
        }
        // second reconstruct is a no-op
        let again = custody.try_reconstruct().expect("reconstruct 2");
        assert_eq!(again.recovered_columns, 0);
    }

    #[test]
    fn tampered_column_rejected_on_insert() {
        let (params, prover, public) = grid(8, 8, b"arch-store-2");
        let mut custody = GridCustody::new([1u8; 32], 1, params, public);
        let mut bad = prover.matrix.col(3);
        bad[0] = bad[0] + Fr::ONE;
        assert!(custody.put_column_verified(3, bad).is_err());
        assert!(custody.get_column(3).is_none());
    }

    #[test]
    fn compliance_report() {
        let (params, prover, public) = grid(8, 8, b"arch-store-3");
        let mut custody = GridCustody::new([2u8; 32], 1, params, public);
        for c in [0usize, 2, 4, 6] {
            custody.put_column_verified(c, prover.matrix.col(c)).unwrap();
        }
        let report: CustodyCompliance = custody.custody_compliance(&[0, 2, 4, 6, 8]);
        assert_eq!(report.assigned, 5);
        assert_eq!(report.present, 4);
        assert_eq!(report.missing, vec![8]);
        assert!((report.compliance() - 0.8).abs() < 1e-12);
    }

    #[test]
    fn budget_pruning_evicts_oldest() {
        let (params_a, prover_a, public_a) = grid(8, 8, b"arch-store-4");
        let (params_b, prover_b, public_b) = grid(8, 8, b"arch-store-5");
        let mut store = CustodyStore::new(usize::MAX);
        let mut a = GridCustody::new([10u8; 32], 100, params_a, public_a);
        let mut b = GridCustody::new([20u8; 32], 200, params_b, public_b);
        for c in 0..8 {
            a.put_column_verified(c, prover_a.matrix.col(c)).unwrap();
            b.put_column_verified(c, prover_b.matrix.col(c)).unwrap();
        }
        let b_bytes = b.byte_size();
        store.put_grid(a);
        store.put_grid(b);
        let total = store.total_bytes();
        // budget only fits the newest grid
        store.set_byte_budget(total - b_bytes + 1);
        let evicted = store.prune_to_budget();
        assert_eq!(evicted, vec![[10u8; 32]]);
        assert_eq!(store.grid_count(), 1);
        assert!(store.total_bytes() <= store.byte_budget());
    }

    #[test]
    fn compact_format_roundtrip() {
        let (params, prover, public) = grid(8, 8, b"arch-fmt-1");
        let mut custody = GridCustody::new([9u8; 32], 7, params, public);
        for c in 0..9 {
            custody.put_column_verified(c, prover.matrix.col(c)).unwrap();
        }
        let bytes = format::encode_grid(&custody);
        // exact-size accounting must match the serialized length
        assert_eq!(custody.byte_size(), bytes.len());
        let back = format::decode_grid(&bytes).expect("decode");
        assert_eq!(back.grid_id, custody.grid_id);
        assert_eq!(back.slot, custody.slot);
        assert_eq!(back.available_columns(), custody.available_columns());
        for c in custody.available_columns() {
            assert_eq!(back.get_column(c), custody.get_column(c));
        }
        // corrupted bytes are caught by the checksum
        let mut corrupt = bytes.clone();
        let idx = corrupt.len() - 40;
        corrupt[idx] ^= 0xff;
        assert!(format::decode_grid(&corrupt).is_err());
    }

    #[test]
    fn disk_roundtrip_and_delete() {
        let dir = std::env::temp_dir().join(format!("zoda-arch-test-{}", std::process::id()));
        let disk = DiskCustody::open(&dir).expect("open");
        let (params, prover, public) = grid(8, 8, b"arch-disk-1");
        let mut custody = GridCustody::new([33u8; 32], 11, params, public);
        for c in 0..8 {
            custody.put_column_verified(c, prover.matrix.col(c)).unwrap();
        }
        disk.save_grid(&custody).expect("save");
        let loaded = disk.load_grid(&[33u8; 32]).expect("load").expect("exists");
        assert_eq!(loaded.available_columns(), vec![0, 1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(loaded.get_column(7).map(|v| v.clone()), Some(prover.matrix.col(7)));
        assert_eq!(disk.list_grids(), vec![[33u8; 32]]);
        let _ = disk.delete_grid(&[33u8; 32]);
        assert!(disk.load_grid(&[33u8; 32]).expect("load").is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
