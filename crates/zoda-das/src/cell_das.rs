//! Cell-level 2D DAS over EIP-7594 blob grids (Fulu PeerDAS flavour).
//!
//! A slot's blob set forms a 2D grid: **rows** are blobs (each extended
//! to 128 cells by the FK20-proven Reed–Solomon code), **columns** are
//! the 128 cell indices across all blobs. Availability is per blob: a
//! blob is ambiguous only if ≥ 65 of its 128 cells are hidden
//! ([`crate::availability::BLOB_MIN_HIDDEN_CELLS`]).
//!
//! The sampler draws column indices (custody-first), fetches every blob's
//! cell in each drawn column with its KZG cell proof, and batch-verifies
//! the whole round with a single pairing ([`verify_cell_kzg_proof_batch`]).
//! It also accumulates per-column cell custody — the columns this node
//! must store and serve per its custody-group assignment.

use crate::availability::{blob_column_miss_probability, BLOB_MIN_HIDDEN_CELLS};
use std::collections::HashSet;
use zoda_edas::{
    compute_cells_and_kzg_proofs, verify_cell_kzg_proof_batch, Cell, CELLS_PER_EXT_BLOB,
    NUMBER_OF_COLUMNS,
};
use zoda_kzg::srs::{Setup, FIELD_ELEMENTS_PER_BLOB};
use zoda_math::ZodaRng;

/// Blob payload bytes (4096 field elements x 32).
pub const BLOB_BYTES: usize = FIELD_ELEMENTS_PER_BLOB * 32;

/// Map `f` over `0..n` (in order), split across cores when worthwhile.
pub(crate) fn parallel_map<T, F>(n: usize, heavy: bool, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Sync + Send,
{
    let cpus = std::thread::available_parallelism()
        .map(|c| c.get())
        .unwrap_or(1);
    if !heavy || n < 2 || cpus < 2 {
        return (0..n).map(|i| f(i)).collect();
    }
    let chunks = cpus.min(n);
    let per = (n + chunks - 1) / chunks;
    let mut results: Vec<Option<Vec<T>>> = (0..chunks).map(|_| None).collect();
    let f = &f;
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(chunks);
        let mut rest: &mut [Option<Vec<T>>] = &mut results;
        for t in 0..chunks {
            let start = t * per;
            if start >= n {
                break;
            }
            let end = (start + per).min(n);
            let (slot, tail) = rest.split_at_mut(1);
            rest = tail;
            handles.push(scope.spawn(move || {
                slot[0] = Some((start..end).map(|i| f(i)).collect());
            }));
        }
        for h in handles {
            h.join().expect("parallel_map thread panicked");
        }
    });
    let mut out = Vec::with_capacity(n);
    for part in results.into_iter().flatten() {
        out.extend(part);
    }
    out
}

/// A committed 2D blob grid: `blobs.len()` rows × 128 columns.
pub struct BlobGrid {
    pub commitments: Vec<Vec<u8>>,
    /// `[blob][column]` — 2048-byte cells
    pub cells: Vec<Vec<Cell>>,
    /// `[blob][column]` — 48-byte compressed KZG proofs
    pub proofs: Vec<Vec<[u8; 48]>>,
}

impl BlobGrid {
    /// Compute every blob's cells, commitments and FK20 cell proofs
    /// (rows are independent; large slots are encoded in parallel).
    pub fn new(blobs: &[Vec<u8>], setup: &Setup) -> Result<BlobGrid, String> {
        if blobs.is_empty() {
            return Err("empty blob grid".into());
        }
        for b in blobs {
            if b.len() != BLOB_BYTES {
                return Err(format!("blob length must be {}, got {}", BLOB_BYTES, b.len()));
            }
        }
        let setup = &setup;
        let rows: Vec<(Vec<Cell>, Vec<[u8; 48]>, Vec<u8>)> =
            parallel_map(blobs.len(), blobs.len() >= 2, |i| {
                let (cells, proofs) = compute_cells_and_kzg_proofs(&blobs[i], setup)
                    .expect("blob length checked above");
                // the blob commitment is the commitment of the (padded)
                // monomial form — c-kzg's blob_to_kzg_commitment semantics
                let commitment = {
                    let poly = zoda_edas::blob_to_monomial(&blobs[i], setup)
                        .expect("blob length checked above");
                    setup.commit_coeffs(&poly).to_compressed().to_vec()
                };
                (cells, proofs, commitment)
            });
        let mut commitments = Vec::with_capacity(rows.len());
        let mut cells = Vec::with_capacity(rows.len());
        let mut proofs = Vec::with_capacity(rows.len());
        for (c, p, com) in rows {
            debug_assert_eq!(c.len(), CELLS_PER_EXT_BLOB);
            cells.push(c);
            proofs.push(p);
            commitments.push(com);
        }
        Ok(BlobGrid {
            commitments,
            cells,
            proofs,
        })
    }

    pub fn blobs(&self) -> usize {
        self.cells.len()
    }

    pub fn columns(&self) -> usize {
        NUMBER_OF_COLUMNS
    }

    pub fn cell(&self, blob: usize, column: usize) -> &Cell {
        &self.cells[blob][column]
    }

    pub fn proof(&self, blob: usize, column: usize) -> &[u8; 48] {
        &self.proofs[blob][column]
    }
}

/// The oracle a cell sampler consults: one (blob, column) cell.
pub trait CellOracle {
    /// `None` = the peer refuses / does not have the cell.
    fn fetch_cell(&mut self, blob: usize, column: usize) -> Option<(Cell, [u8; 48])>;
    /// Identifies the serving peer (for scoring).
    fn peer_id(&self) -> usize {
        0
    }
}

/// Honest full-grid oracle with optional withholding and corruption
/// (tests and adversarial simulations).
pub struct GridOracle<'a> {
    pub grid: &'a BlobGrid,
    /// (blob, column) pairs the adversary refuses to serve
    pub hidden: HashSet<(usize, usize)>,
}

impl<'a> GridOracle<'a> {
    pub fn honest(grid: &'a BlobGrid) -> GridOracle<'a> {
        GridOracle {
            grid,
            hidden: HashSet::new(),
        }
    }

    /// Withhold a whole blob's extension (blob fully unavailable).
    pub fn withholding_blob(grid: &'a BlobGrid, blob: usize) -> GridOracle<'a> {
        let mut hidden = HashSet::new();
        for c in 0..NUMBER_OF_COLUMNS {
            hidden.insert((blob, c));
        }
        GridOracle { grid, hidden }
    }

    /// The minimal ambiguity: hide `BLOB_MIN_HIDDEN_CELLS` cells of one blob.
    pub fn withholding_minimal(grid: &'a BlobGrid, blob: usize) -> GridOracle<'a> {
        let mut hidden = HashSet::new();
        for c in 0..BLOB_MIN_HIDDEN_CELLS {
            hidden.insert((blob, c));
        }
        GridOracle { grid, hidden }
    }
}

impl<'a> CellOracle for GridOracle<'a> {
    fn fetch_cell(&mut self, blob: usize, column: usize) -> Option<(Cell, [u8; 48])> {
        if self.hidden.contains(&(blob, column)) {
            None
        } else {
            Some((*self.grid.cell(blob, column), *self.grid.proof(blob, column)))
        }
    }
}

/// Per-round / per-session bookkeeping of a cell sampling session.
#[derive(Clone, Debug, Default)]
pub struct CellSessionReport {
    pub columns_drawn: Vec<usize>,
    pub cells_fetched: usize,
    pub cells_verified: usize,
    pub cells_rejected: usize,
    pub cells_missing: usize,
    /// per-blob count of verified cells (reconstruction readiness)
    pub blob_verified: Vec<usize>,
    /// per-column count of verified cells (custody readiness)
    pub column_verified: Vec<usize>,
    /// worst-case escape probability of a withheld blob given the drawn
    /// columns (exact hypergeometric bound)
    pub miss_probability: f64,
}

impl CellSessionReport {
    pub fn available(&self) -> bool {
        self.cells_rejected == 0
    }

    pub fn confidence(&self) -> f64 {
        1.0 - self.miss_probability
    }
}

/// Session configuration.
#[derive(Clone, Debug)]
pub struct CellSessionConfig {
    /// column draws per round (spec default: `SAMPLES_PER_SLOT = 8`)
    pub columns_per_round: usize,
    /// target worst-case miss probability (default 2^-20)
    pub target_miss: f64,
    /// maximum rounds (budget)
    pub max_rounds: usize,
    /// own custody columns — always drawn first
    pub custody_columns: Vec<usize>,
}

impl Default for CellSessionConfig {
    fn default() -> Self {
        CellSessionConfig {
            columns_per_round: zoda_edas::SAMPLES_PER_SLOT as usize,
            target_miss: 2f64.powi(-20),
            max_rounds: 16,
            custody_columns: Vec::new(),
        }
    }
}

/// Custody columns for a node, derived from its custody-group assignment
/// (EIP-7594: groups → column indices).
pub fn custody_columns_for(node_id: [u8; 32], custody_group_count: u64) -> Vec<usize> {
    zoda_edas::get_custody_groups(node_id, custody_group_count)
        .into_iter()
        .flat_map(|g| zoda_edas::compute_columns_for_custody_group(g))
        .map(|c| c as usize)
        .filter(|&c| c < NUMBER_OF_COLUMNS)
        .collect()
}

/// Run a full cell-level sampling session over a blob grid:
/// rounds of column draws (custody columns first), batch verification of
/// every fetched cell, and exact confidence accounting.
pub fn sample_cell_session(
    grid_meta: &BlobGridMeta,
    oracle: &mut dyn CellOracle,
    config: &CellSessionConfig,
    rng: &mut ZodaRng,
    setup: &Setup,
) -> CellSessionReport {
    let blobs = grid_meta.blobs;
    let mut report = CellSessionReport {
        blob_verified: vec![0; blobs],
        column_verified: vec![0; NUMBER_OF_COLUMNS],
        ..Default::default()
    };
    let mut drawn: HashSet<usize> = HashSet::new();
    // custody columns first (in deterministic order)
    let mut queue: Vec<usize> = config
        .custody_columns
        .iter()
        .copied()
        .filter(|&c| c < NUMBER_OF_COLUMNS)
        .collect();
    queue.sort_unstable();
    queue.dedup();

    let mut rounds = 0usize;
    loop {
        // top up the queue with fresh uniform columns
        let mut fresh = 0usize;
        while queue.len() < config.columns_per_round && fresh < 4 * NUMBER_OF_COLUMNS {
            // Fisher–Yates-style draw of an undrawn column
            let c = (rng.next_u64() as usize) % NUMBER_OF_COLUMNS;
            if !drawn.contains(&c) && !queue.contains(&c) {
                queue.push(c);
            }
            fresh += 1;
        }
        if queue.is_empty() {
            break;
        }
        let take = queue.len().min(config.columns_per_round);
        let round_cols: Vec<usize> = queue.drain(..take).collect();
        // fetch every blob's cell in each drawn column
        let mut commitments: Vec<&[u8]> = Vec::new();
        let mut cell_indices: Vec<u64> = Vec::new();
        let mut cells: Vec<Cell> = Vec::new();
        let mut proofs: Vec<[u8; 48]> = Vec::new();
        let mut owners: Vec<(usize, usize)> = Vec::new();
        for &c in &round_cols {
            for b in 0..blobs {
                report.cells_fetched += 1;
                match oracle.fetch_cell(b, c) {
                    Some((cell, proof)) => {
                        commitments.push(&grid_meta.commitments[b]);
                        cell_indices.push(c as u64);
                        cells.push(cell);
                        proofs.push(proof);
                        owners.push((b, c));
                    }
                    None => report.cells_missing += 1,
                }
            }
        }
        // one batch verification for the whole round (single pairing)
        let proof_refs: Vec<&[u8]> = proofs.iter().map(|p| p.as_slice()).collect();
        let ok = commitments.is_empty()
            || verify_cell_kzg_proof_batch(&commitments, &cell_indices, &cells, &proof_refs, setup)
                .unwrap_or(false);
        if ok {
            for &(b, c) in &owners {
                report.blob_verified[b] += 1;
                report.column_verified[c] += 1;
            }
            report.cells_verified += owners.len();
        } else {
            // identify the offending cells: bisect the batch
            for (i, &(b, c)) in owners.iter().enumerate() {
                let one = verify_cell_kzg_proof_batch(
                    &[commitments[i]],
                    &[cell_indices[i]],
                    &[cells[i]],
                    &[proofs[i].as_slice()],
                    setup,
                )
                .unwrap_or(false);
                if one {
                    report.blob_verified[b] += 1;
                    report.column_verified[c] += 1;
                    report.cells_verified += 1;
                } else {
                    report.cells_rejected += 1;
                }
            }
        }
        for &c in &round_cols {
            drawn.insert(c);
            report.columns_drawn.push(c);
        }
        rounds += 1;
        // exact confidence for the union of drawn columns
        report.miss_probability =
            blob_column_miss_probability(drawn.len(), BLOB_MIN_HIDDEN_CELLS, NUMBER_OF_COLUMNS);
        if report.cells_rejected > 0 {
            break; // fail closed
        }
        if report.miss_probability <= config.target_miss {
            break;
        }
        if rounds >= config.max_rounds {
            break;
        }
    }
    report
}

/// The public metadata of a blob grid a sampler needs (commitments only —
/// the sampler never needs the cells themselves).
#[derive(Clone)]
pub struct BlobGridMeta {
    pub commitments: Vec<Vec<u8>>,
    pub blobs: usize,
}

impl BlobGridMeta {
    pub fn from_grid(grid: &BlobGrid) -> BlobGridMeta {
        BlobGridMeta {
            commitments: grid.commitments.clone(),
            blobs: grid.blobs(),
        }
    }
}

/// Reconstruct a blob's cells from any 64 of 128 (EIP-7594
/// `recover_cells_and_kzg_proofs`) and return the full 128 with proofs.
pub fn recover_blob(
    cell_indices: &[u64],
    cells: &[Cell],
    setup: &Setup,
) -> Result<(Vec<Cell>, Vec<[u8; 48]>), String> {
    zoda_edas::recover_cells_and_kzg_proofs(cell_indices, cells, setup)
}

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_kzg::srs::Setup;

    fn toy_setup() -> Setup {
        Setup::from_seed_for_testing(*b"cell-das-toy-setup-tau-000000000")
    }

    fn blobs(n: usize) -> Vec<Vec<u8>> {
        (0..n)
            .map(|i| {
                let mut seed = [0u8; 32];
                seed[0] = i as u8;
                seed[1] = 42;
                let mut rng = ZodaRng::from_seed(seed);
                let mut b = vec![0u8; BLOB_BYTES];
                for chunk in b.chunks_mut(32) {
                    // random bytes only in the LOW half: BE value < 2^128 < r
                    let r = rng.next_u64().to_le_bytes();
                    let r2 = rng.next_u64().to_le_bytes();
                    chunk[16..24].copy_from_slice(&r);
                    chunk[24..32].copy_from_slice(&r2);
                }
                b
            })
            .collect()
    }

    #[test]
    fn honest_grid_session_accepts() {
        let setup = toy_setup();
        let bl = blobs(2);
        let grid = BlobGrid::new(&bl, &setup).expect("grid");
        let meta = BlobGridMeta::from_grid(&grid);
        let mut oracle = GridOracle::honest(&grid);
        let mut rng = ZodaRng::from_seed(*b"cell-das-sess-seed-0000000000000");
        let config = CellSessionConfig::default();
        let report = sample_cell_session(&meta, &mut oracle, &config, &mut rng, &setup);
        assert!(report.available());
        assert_eq!(report.cells_rejected, 0);
        assert!(report.cells_verified >= 2 * 8); // 2 blobs x >=8 columns
        // 8 columns give ~2^-12; the session must keep going toward 2^-20
        assert!(report.miss_probability <= config.target_miss);
        assert!(report.columns_drawn.len() >= 8);
    }

    #[test]
    fn withheld_blob_detected() {
        let setup = toy_setup();
        let bl = blobs(2);
        let grid = BlobGrid::new(&bl, &setup).expect("grid");
        let meta = BlobGridMeta::from_grid(&grid);
        // hide a whole blob (every column serves nothing for blob 0)
        let mut oracle = GridOracle::withholding_blob(&grid, 0);
        let mut rng = ZodaRng::from_seed(*b"cell-das-hide-seed-0000000000000");
        let report = sample_cell_session(
            &meta,
            &mut oracle,
            &CellSessionConfig::default(),
            &mut rng,
            &setup,
        );
        // a fully hidden blob shows up as missing fetches (not rejections)
        assert!(report.cells_missing > 0);
        // the session still verifies the honest blob's cells
        assert!(report.blob_verified[1] > 0);
    }

    #[test]
    fn minimal_withholding_drives_more_rounds() {
        let setup = toy_setup();
        let bl = blobs(1);
        let grid = BlobGrid::new(&bl, &setup).expect("grid");
        let meta = BlobGridMeta::from_grid(&grid);
        let mut oracle = GridOracle::withholding_minimal(&grid, 0);
        let mut rng = ZodaRng::from_seed(*b"cell-das-min-seed-00000000000000");
        let config = CellSessionConfig {
            max_rounds: 16,
            ..Default::default()
        };
        let report = sample_cell_session(&meta, &mut oracle, &config, &mut rng, &setup);
        // 65 of 128 columns serve nothing; any draw touching them misses.
        // The session runs to its confidence target or round budget.
        assert!(report.cells_missing > 0 || report.columns_drawn.len() > 8);
    }

    #[test]
    fn custody_columns_first() {
        // a node with the minimum custody duty (4 groups, 1 column each)
        let cols = custody_columns_for([7u8; 32], 4);
        assert_eq!(cols.len(), 4, "cols = {:?}", cols);
        assert!(cols.iter().all(|&c| c < 128));
        // deterministic
        assert_eq!(cols, custody_columns_for([7u8; 32], 4));
        // different nodes see different columns
        let other = custody_columns_for([9u8; 32], 4);
        assert_ne!(cols, other);
    }
}
