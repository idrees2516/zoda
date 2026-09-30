//! In-memory custody: one grid's verified columns, multi-grid stores with
//! byte budgets, and reconstruction that persists its result.

use std::collections::HashMap;
use zoda_core::reconstruct2::{PartialGrid, ReconError};
use zoda_core::{verify_column_sample, Matrix, ReconStats, ZodaParams, ZodaPublic};
use zoda_math::Fr;

/// Custody of one grid: the verified columns this node holds.
pub struct GridCustody {
    pub grid_id: [u8; 32],
    /// slot / ordering key for pruning (higher = newer)
    pub slot: u64,
    pub params: ZodaParams,
    pub public: ZodaPublic,
    columns: Vec<Option<Vec<Fr>>>,
}

/// Result of a successful reconstruction pass.
#[derive(Clone, Copy, Debug)]
pub struct ReconOutcome {
    /// columns that were absent before and are now stored
    pub recovered_columns: usize,
    pub stats: ReconStats,
}

impl GridCustody {
    pub fn new(grid_id: [u8; 32], slot: u64, params: ZodaParams, public: ZodaPublic) -> GridCustody {
        let ext_cols = params.extended_cols();
        GridCustody {
            grid_id,
            slot,
            params,
            public,
            columns: vec![None; ext_cols],
        }
    }

    /// Store a column AFTER verifying it against the public projections.
    /// Returns `true` if the column is new, `false` if it replaced an
    /// identical existing entry. A tampered column is rejected with Err.
    pub fn put_column_verified(&mut self, index: usize, col: Vec<Fr>) -> Result<bool, String> {
        if index >= self.columns.len() {
            return Err("column index out of range".to_string());
        }
        if col.len() != self.params.extended_rows() {
            return Err("column has wrong height".to_string());
        }
        verify_column_sample(&self.public, index, &col)?;
        let fresh = self.columns[index].is_none();
        self.columns[index] = Some(col);
        Ok(fresh)
    }

    pub fn get_column(&self, index: usize) -> Option<&Vec<Fr>> {
        self.columns[index].as_ref()
    }

    pub fn available_columns(&self) -> Vec<usize> {
        (0..self.columns.len())
            .filter(|&i| self.columns[i].is_some())
            .collect()
    }

    pub fn missing_columns(&self) -> Vec<usize> {
        (0..self.columns.len())
            .filter(|&i| self.columns[i].is_none())
            .collect()
    }

    pub fn is_complete(&self) -> bool {
        self.columns.iter().all(|c| c.is_some())
    }

    /// Availability bitfield over the extended columns (1 = present),
    /// byte-packed, column `i` → bit `i % 8` of byte `i / 8`.
    pub fn column_bitfield(&self) -> Vec<u8> {
        let n = self.columns.len();
        let mut out = vec![0u8; (n + 7) / 8];
        for (i, c) in self.columns.iter().enumerate() {
            if c.is_some() {
                out[i / 8] |= 1 << (i % 8);
            }
        }
        out
    }

    pub fn set_column_bitfield(&mut self, bits: &[u8]) {
        for (i, c) in self.columns.iter_mut().enumerate() {
            let present = bits.get(i / 8).map(|b| b & (1 << (i % 8)) != 0).unwrap_or(false);
            if !present {
                *c = None;
            }
        }
    }

    pub fn column_count(&self) -> usize {
        self.columns.iter().filter(|c| c.is_some()).count()
    }

    /// The custody-duty compliance of this node for an assigned column set.
    pub fn custody_compliance(&self, assigned: &[usize]) -> CustodyCompliance {
        let mut missing = Vec::new();
        let mut present = 0usize;
        for &c in assigned {
            if c < self.columns.len() && self.columns[c].is_some() {
                present += 1;
            } else {
                missing.push(c);
            }
        }
        CustodyCompliance {
            assigned: assigned.len(),
            present,
            missing,
        }
    }

    /// Exact serialized size of this custody (the compact disk format).
    pub fn byte_size(&self) -> usize {
        crate::format::grid_file_size(
            self.params.m,
            self.params.k,
            self.column_count(),
        )
    }

    /// Reconstruct the grid from the stored columns once ≥ k are present,
    /// **persisting every recovered column back into custody** (each is
    /// re-verified against the public projections before acceptance).
    pub fn try_reconstruct(&mut self) -> Result<ReconOutcome, String> {
        let k = self.params.k;
        let avail = self.available_columns();
        if avail.len() < k {
            return Err(format!(
                "need at least {} columns to reconstruct, have {}",
                k,
                avail.len()
            ));
        }
        let mut partial = PartialGrid::new(self.params.extended_rows(), self.params.extended_cols());
        for &c in &avail {
            partial.insert_column(c, self.columns[c].as_ref().expect("listed above"));
        }
        match zoda_core::reconstruct_2d(&mut partial, &self.params, &self.public) {
            Ok(stats) => {
                let full: Matrix = partial.to_matrix();
                let mut recovered = 0usize;
                for c in 0..self.columns.len() {
                    if self.columns[c].is_none() {
                        let col = full.col(c);
                        // defense in depth: re-verify before storing
                        verify_column_sample(&self.public, c, &col)?;
                        self.columns[c] = Some(col);
                        recovered += 1;
                    }
                }
                Ok(ReconOutcome {
                    recovered_columns: recovered,
                    stats,
                })
            }
            Err(ReconError::Incomplete { stats, fetch_next }) => Err(format!(
                "reconstruction incomplete after {} rounds ({} cells recovered); fetch {} more cells",
                stats.rounds,
                stats.cells_recovered,
                fetch_next.len()
            )),
            Err(ReconError::CorruptRow(r)) => Err(format!("inconsistent data in row {}", r)),
            Err(ReconError::CorruptColumn(c)) => Err(format!("inconsistent data in column {}", c)),
            Err(ReconError::ProjectionMismatch(m)) => {
                Err(format!("recovered grid failed verification: {}", m))
            }
        }
    }
}

/// Custody-duty compliance report.
#[derive(Clone, Debug)]
pub struct CustodyCompliance {
    pub assigned: usize,
    pub present: usize,
    pub missing: Vec<usize>,
}

impl CustodyCompliance {
    /// Fraction of the assigned duty currently satisfied.
    pub fn compliance(&self) -> f64 {
        if self.assigned == 0 {
            1.0
        } else {
            self.present as f64 / self.assigned as f64
        }
    }
}

/// A multi-grid custody store with a byte budget.
pub struct CustodyStore {
    grids: HashMap<[u8; 32], GridCustody>,
    byte_budget: usize,
}

impl CustodyStore {
    /// `byte_budget = 0` means unlimited.
    pub fn new(byte_budget: usize) -> CustodyStore {
        CustodyStore {
            grids: HashMap::new(),
            byte_budget,
        }
    }

    pub fn byte_budget(&self) -> usize {
        self.byte_budget
    }

    pub fn set_byte_budget(&mut self, budget: usize) {
        self.byte_budget = budget;
    }

    pub fn put_grid(&mut self, grid: GridCustody) {
        self.grids.insert(grid.grid_id, grid);
    }

    pub fn grid(&self, grid_id: &[u8; 32]) -> Option<&GridCustody> {
        self.grids.get(grid_id)
    }

    pub fn grid_mut(&mut self, grid_id: &[u8; 32]) -> Option<&mut GridCustody> {
        self.grids.get_mut(grid_id)
    }

    pub fn remove_grid(&mut self, grid_id: &[u8; 32]) -> Option<GridCustody> {
        self.grids.remove(grid_id)
    }

    pub fn grid_count(&self) -> usize {
        self.grids.len()
    }

    pub fn total_bytes(&self) -> usize {
        self.grids.values().map(|g| g.byte_size()).sum()
    }

    /// Ingest a verified column into a grid, creating the custody entry
    /// on first sight (the caller supplies the grid's params/public).
    pub fn ingest_verified_column(
        &mut self,
        grid_id: &[u8; 32],
        slot: u64,
        params: &ZodaParams,
        public: &ZodaPublic,
        index: usize,
        col: Vec<Fr>,
    ) -> Result<bool, String> {
        let grid = self
            .grids
            .entry(*grid_id)
            .or_insert_with(|| GridCustody::new(*grid_id, slot, ZodaParams::new(params.m, params.k), public.clone()));
        grid.put_column_verified(index, col)
    }

    /// Evict oldest-slot grids until the store fits its byte budget.
    /// Returns the evicted grid ids (oldest first). Unlimited budgets
    /// evict nothing.
    pub fn prune_to_budget(&mut self) -> Vec<[u8; 32]> {
        let mut evicted = Vec::new();
        if self.byte_budget == 0 {
            return evicted;
        }
        while self.total_bytes() > self.byte_budget {
            let victim = self
                .grids
                .values()
                .map(|g| (g.slot, g.grid_id))
                .min()
                .map(|(_, id)| id);
            match victim {
                Some(id) => {
                    self.grids.remove(&id);
                    evicted.push(id);
                }
                None => break,
            }
        }
        evicted
    }
}
