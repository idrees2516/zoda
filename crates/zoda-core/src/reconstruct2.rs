//! 2D erasure decoding for the ZODA tensor code — arbitrary cell-loss
//! patterns, not just whole missing columns.
//!
//! The extended grid is a `2m × 2k` tensor codeword: every row is a
//! length-`2k` Reed–Solomon codeword (degree `< k`) and every column is a
//! length-`2m` codeword (degree `< m`). The decoder runs the classical
//! row/column fixpoint of product codes:
//!
//! ```text
//! loop
//!   every row    with ≥ k of 2k cells present  → interpolate, fill the rest
//!   every column with ≥ m of 2m cells present  → interpolate, fill the rest
//! until a full pass recovers nothing
//! ```
//!
//! ## Correctness and security
//!
//! * Each line is interpolated from the first `k` (resp. `m`) present
//!   cells and then **checked against every remaining present cell** — a
//!   line containing a tampered cell fails closed with its index.
//! * After the fixpoint, every row and column of the result is verified
//!   against the public Fiat–Shamir projections
//!   ([`verify_row_sample`] / [`verify_column_sample`]) — this catches
//!   consistent-but-wrong decodings (exactly-`k` adversarial cells), which
//!   no amount of re-encoding can detect, because the projections bind
//!   the codeword to the original commitment.
//! * If the fixpoint stalls before completion the caller learns which
//!   cells to fetch next ([`PartialGrid::next_cells_to_fetch`]).
//!
//! [`verify_row_sample`]: crate::sample::verify_row_sample

use crate::encode::Matrix;
use crate::params::{ZodaParams, ZodaPublic};
use crate::sample::{verify_column_sample, verify_row_sample};
use std::collections::HashMap;
use std::rc::Rc;
use zoda_math::{batch_invert, FftDomain, Fr, PrimeField};

/// A partially received `rows × cols` grid: values plus a presence bitmap.
#[derive(Clone, Debug)]
pub struct PartialGrid {
    rows: usize,
    cols: usize,
    values: Vec<Fr>,
    /// presence bitmap, row-major, `words_per_row()` u64 words per row
    present: Vec<u64>,
}

/// 64-bit words needed for a `cols`-wide presence row.
#[inline]
fn words_for(cols: usize) -> usize {
    (cols + 63) / 64
}

impl PartialGrid {
    pub fn new(rows: usize, cols: usize) -> PartialGrid {
        PartialGrid {
            rows,
            cols,
            values: vec![Fr::zero(); rows * cols],
            present: vec![0; rows * words_for(cols)],
        }
    }

    /// Everything present (e.g. the prover's own view).
    pub fn from_full(m: &Matrix) -> PartialGrid {
        let mut p = PartialGrid::new(m.rows(), m.cols());
        for r in 0..m.rows() {
            let row = m.row_slice(r);
            for (c, v) in row.iter().enumerate() {
                p.insert_cell(r, c, *v);
            }
        }
        p
    }

    /// A full view with a set of cells erased (indices must be distinct).
    pub fn with_erasures(m: &Matrix, erased: &[(usize, usize)]) -> PartialGrid {
        let mut p = PartialGrid::from_full(m);
        for &(r, c) in erased {
            p.erase_cell(r, c);
        }
        p
    }

    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    #[inline]
    pub fn is_present(&self, r: usize, c: usize) -> bool {
        let w = c / 64;
        let b = c % 64;
        self.present[r * words_for(self.cols) + w] & (1u64 << b) != 0
    }

    #[inline]
    pub fn insert_cell(&mut self, r: usize, c: usize, v: Fr) {
        assert!(r < self.rows && c < self.cols, "cell index out of range");
        let w = words_for(self.cols);
        self.values[r * self.cols + c] = v;
        self.present[r * w + c / 64] |= 1u64 << (c % 64);
    }

    #[inline]
    pub fn erase_cell(&mut self, r: usize, c: usize) {
        let w = words_for(self.cols);
        self.present[r * w + c / 64] &= !(1u64 << (c % 64));
    }

    #[inline]
    pub fn get(&self, r: usize, c: usize) -> Option<Fr> {
        if self.is_present(r, c) {
            Some(self.values[r * self.cols + c])
        } else {
            None
        }
    }

    /// Insert a full received column (height must match).
    pub fn insert_column(&mut self, c: usize, col: &[Fr]) {
        assert_eq!(col.len(), self.rows, "column height mismatch");
        for (r, v) in col.iter().enumerate() {
            self.insert_cell(r, c, *v);
        }
    }

    /// Insert a full received row (width must match).
    pub fn insert_row(&mut self, r: usize, row: &[Fr]) {
        assert_eq!(row.len(), self.cols, "row width mismatch");
        for (c, v) in row.iter().enumerate() {
            self.insert_cell(r, c, *v);
        }
    }

    /// The `(column, value)` pairs present in row `r`.
    pub fn row_cells(&self, r: usize) -> (Vec<usize>, Vec<Fr>) {
        let mut idx = Vec::new();
        let mut vals = Vec::new();
        for c in 0..self.cols {
            if self.is_present(r, c) {
                idx.push(c);
                vals.push(self.values[r * self.cols + c]);
            }
        }
        (idx, vals)
    }

    /// The `(row, value)` pairs present in column `c`.
    pub fn col_cells(&self, c: usize) -> (Vec<usize>, Vec<Fr>) {
        let mut idx = Vec::new();
        let mut vals = Vec::new();
        for r in 0..self.rows {
            if self.is_present(r, c) {
                idx.push(r);
                vals.push(self.values[r * self.cols + c]);
            }
        }
        (idx, vals)
    }

    pub fn present_count(&self) -> usize {
        self.present.iter().map(|w| w.count_ones() as usize).sum()
    }

    pub fn is_complete(&self) -> bool {
        self.present_count() == self.rows * self.cols
    }

    /// Number of present cells in row `r` / column `c`.
    pub fn row_present(&self, r: usize) -> usize {
        let w = words_for(self.cols);
        self.present[r * w..(r + 1) * w]
            .iter()
            .map(|x| x.count_ones() as usize)
            .sum()
    }

    pub fn col_present(&self, c: usize) -> usize {
        (0..self.rows).filter(|&r| self.is_present(r, c)).count()
    }

    /// The recovered grid (only meaningful after a successful decode).
    pub fn to_matrix(&self) -> Matrix {
        Matrix::from_row_major(self.rows, self.cols, self.values.clone())
    }

    /// Greedy repair guidance: the missing cells of the lines closest to
    /// decodability first (a line `t` cells short of the threshold yields
    /// `t` fetches that unlock it). Returns at most `budget` cells.
    pub fn next_cells_to_fetch(&self, k: usize, m: usize, budget: usize) -> Vec<(usize, usize)> {
        // (deficit, is_row, index) — lines with the smallest positive deficit
        let mut lines: Vec<(usize, bool, usize)> = Vec::new();
        for r in 0..self.rows {
            let p = self.row_present(r);
            if p < k {
                lines.push((k - p, true, r));
            }
        }
        for c in 0..self.cols {
            let p = self.col_present(c);
            if p < m {
                lines.push((m - p, false, c));
            }
        }
        lines.sort_unstable();
        let mut out = Vec::new();
        'outer: for &(deficit, is_row, idx) in &lines {
            let mut taken = 0;
            if is_row {
                for c in 0..self.cols {
                    if !self.is_present(idx, c) {
                        out.push((idx, c));
                        taken += 1;
                        if taken == deficit || out.len() == budget {
                            break;
                        }
                    }
                }
            } else {
                for r in 0..self.rows {
                    if !self.is_present(r, idx) {
                        out.push((r, idx));
                        taken += 1;
                        if taken == deficit || out.len() == budget {
                            break;
                        }
                    }
                }
            }
            if out.len() == budget {
                break 'outer;
            }
        }
        out
    }
}

/// Statistics of one 2D reconstruction.
#[derive(Clone, Copy, Debug, Default)]
pub struct ReconStats {
    pub rounds: usize,
    pub rows_decoded: usize,
    pub cols_decoded: usize,
    pub cells_recovered: usize,
}

/// Why a 2D reconstruction failed.
#[derive(Clone, Debug)]
pub enum ReconError {
    /// A present cell is inconsistent with the line's codeword — data
    /// corruption (or the wrong grid parameters). Index of the offending
    /// row/column included.
    CorruptRow(usize),
    CorruptColumn(usize),
    /// The fixpoint stalled before recovering every cell. The partial
    /// grid (still valid) plus fetch guidance is in the payload.
    Incomplete {
        stats: ReconStats,
        fetch_next: Vec<(usize, usize)>,
    },
    /// The recovered grid failed the public projection checks — an
    /// adversary served a consistent but wrong decoding.
    ProjectionMismatch(String),
}

/// Run the row/column fixpoint over `partial`, then verify the result
/// against the public projections. On success every cell is present and
/// cryptographically bound to the committed codeword.
pub fn reconstruct_2d(
    partial: &mut PartialGrid,
    params: &ZodaParams,
    public: &ZodaPublic,
) -> Result<ReconStats, ReconError> {
    let rows = partial.rows();
    let cols = partial.cols();
    if rows != 2 * params.m || cols != 2 * params.k {
        return Err(ReconError::ProjectionMismatch(
            "partial grid shape does not match params".into(),
        ));
    }
    let mut stats = ReconStats::default();
    let mut row_decoder = LineDecoder::new(&params.row_domain, params.k);
    let mut col_decoder = LineDecoder::new(&params.col_domain, params.m);

    loop {
        let mut changed = false;
        stats.rounds += 1;
        // ---- row pass: rows with ≥ k present cells and at least one gap
        for r in 0..rows {
            let present = partial.row_present(r);
            if present < params.k || present == cols {
                continue;
            }
            let (idx, vals) = partial.row_cells(r);
            match row_decoder.decode(&idx, &vals) {
                Ok(full) => {
                    for c in 0..cols {
                        if !partial.is_present(r, c) {
                            partial.insert_cell(r, c, full[c]);
                            stats.cells_recovered += 1;
                            changed = true;
                        }
                    }
                    stats.rows_decoded += 1;
                }
                Err(_) => return Err(ReconError::CorruptRow(r)),
            }
        }
        // ---- column pass
        for c in 0..cols {
            let present = partial.col_present(c);
            if present < params.m || present == rows {
                continue;
            }
            let (idx, vals) = partial.col_cells(c);
            match col_decoder.decode(&idx, &vals) {
                Ok(full) => {
                    for r in 0..rows {
                        if !partial.is_present(r, c) {
                            partial.insert_cell(r, c, full[r]);
                            stats.cells_recovered += 1;
                            changed = true;
                        }
                    }
                    stats.cols_decoded += 1;
                }
                Err(_) => return Err(ReconError::CorruptColumn(c)),
            }
        }
        if !changed {
            break;
        }
        if stats.rounds > 2 * (rows + cols) {
            // defensive: cannot happen for genuine tensor codewords
            break;
        }
    }

    if !partial.is_complete() {
        let fetch_next =
            partial.next_cells_to_fetch(params.k, params.m, 4 * params.k.max(params.m));
        return Err(ReconError::Incomplete { stats, fetch_next });
    }

    // ---- final binding: every line must satisfy the public projections.
    let matrix = partial.to_matrix();
    for r in 0..rows {
        if let Err(e) = verify_row_sample(public, r, matrix.row_slice(r)) {
            return Err(ReconError::ProjectionMismatch(format!("row {}: {}", r, e)));
        }
    }
    for c in 0..cols {
        if let Err(e) = verify_column_sample(public, c, &matrix.col(c)) {
            return Err(ReconError::ProjectionMismatch(format!("column {}: {}", c, e)));
        }
    }
    Ok(stats)
}

/// Verify that a complete matrix is a tensor codeword consistent with the
/// public projections (re-encoding check + projection check per line).
pub fn verify_codeword(matrix: &Matrix, params: &ZodaParams, public: &ZodaPublic) -> Result<(), String> {
    for r in 0..matrix.rows() {
        let row = matrix.row_slice(r);
        let re = crate::encode::rs_encode_vector(&row[..params.k], &params.row_domain);
        if &re[..] != row {
            return Err(format!("row {} is not a row-codeword", r));
        }
    }
    for c in 0..matrix.cols() {
        let col = matrix.col(c);
        let re = crate::encode::rs_encode_vector(&col[..params.m], &params.col_domain);
        if re != col {
            return Err(format!("column {} is not a column-codeword", c));
        }
    }
    for r in 0..matrix.rows() {
        verify_row_sample(public, r, matrix.row_slice(r))?;
    }
    for c in 0..matrix.cols() {
        verify_column_sample(public, c, &matrix.col(c))?;
    }
    Ok(())
}

/// One line (row or column) interpolator with a basis cache keyed by the
/// presence pattern: rows sharing the same available columns reuse the
/// O(k²) Lagrange basis.
struct LineDecoder {
    domain: FftDomain<Fr>,
    k: usize,
    cache: HashMap<Vec<u64>, Rc<Vec<Vec<Fr>>>>,
}

impl LineDecoder {
    fn new(domain: &FftDomain<Fr>, k: usize) -> LineDecoder {
        LineDecoder {
            domain: FftDomain::new(domain.size()),
            k,
            cache: HashMap::new(),
        }
    }

    /// Interpolate the line from its FIRST k present cells, evaluate on
    /// the full domain, and check consistency with every remaining
    /// present cell. `idx` must be sorted ascending; `idx.len() >= k`.
    fn decode(&mut self, idx: &[usize], vals: &[Fr]) -> Result<Vec<Fr>, ()> {
        debug_assert!(idx.len() >= self.k);
        let k = self.k;
        let width = self.domain.size();
        let roots = self.domain.roots_of_unity();
        let use_idx = &idx[..k];
        let use_vals = &vals[..k];
        // basis cache key: presence word of the first k indices
        let mut key = vec![0u64; words_for_key(width)];
        for &i in use_idx {
            key[i / 64] |= 1u64 << (i % 64);
        }
        let basis = match self.cache.get(&key) {
            Some(b) => b.clone(),
            None => {
                let xs: Vec<Fr> = use_idx.iter().map(|&i| roots[i]).collect();
                let b = Rc::new(lagrange_basis(&xs, &self.domain));
                self.cache.insert(key, b.clone());
                b
            }
        };
        // coeffs = Σ_j y_j · L_j  (O(k²) weighted sums, zero-skipping)
        let mut coeffs = vec![Fr::zero(); k];
        for j in 0..k {
            let y = use_vals[j];
            if y.is_zero() {
                continue;
            }
            let lj = &basis[j];
            for t in 0..k {
                coeffs[t] = coeffs[t] + y * lj[t];
            }
        }
        let mut full = vec![Fr::zero(); width];
        self.domain.fft_padded(&coeffs, &mut full);
        // consistency: every extra present cell must agree
        for (&i, &v) in idx.iter().zip(vals.iter()).skip(k) {
            if full[i] != v {
                return Err(());
            }
        }
        Ok(full)
    }
}

fn words_for_key(width: usize) -> usize {
    (width + 63) / 64
}

/// The monomial-coefficient Lagrange basis for the nodes `xs`:
/// `L_j(t)` coefficients of degree < xs.len(), normalized to L_j(x_i) = δ_ij.
fn lagrange_basis(xs: &[Fr], domain: &FftDomain<Fr>) -> Vec<Vec<Fr>> {
    let k = xs.len();
    // Z = vanishing polynomial of the nodes; L_j = Z / (X − x_j) · 1/Z'(x_j)
    let z = zoda_math::poly::vanishing_poly(xs);
    let zprime = zoda_math::poly::derivative(&z);
    let mut denominators: Vec<Fr> = xs
        .iter()
        .map(|x| zoda_math::poly::eval(&zprime, *x))
        .collect();
    batch_invert(&mut denominators);
    let mut basis = Vec::with_capacity(k);
    for j in 0..k {
        // exact synthetic division of Z (degree k) by (X − x_j)
        let mut q = vec![Fr::zero(); k];
        let mut carry = Fr::zero();
        for t in (0..k).rev() {
            carry = z[t + 1] + carry * xs[j];
            q[t] = carry;
        }
        for t in 0..k {
            q[t] = q[t] * denominators[j];
        }
        basis.push(q);
    }
    let _ = domain; // (kept for symmetry / future coset variants)
    basis
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::ZodaParams;
    use zoda_math::ZodaRng;

    fn grid(m: usize, k: usize, seed: &[u8]) -> (ZodaParams, Matrix, ZodaPublic) {
        let mut seed32 = [0u8; 32];
        seed32[..seed.len()].copy_from_slice(seed);
        let mut rng = ZodaRng::from_seed(seed32);
        let mut data = vec![Fr::zero(); m * k];
        for x in data.iter_mut() {
            *x = rng.next_fr(false);
        }
        let g = Matrix::from_row_major(m, k, data);
        let params = ZodaParams::new(m, k);
        let prover = params.commit(&g);
        (params, prover.matrix.clone(), prover.public_params())
    }

    #[test]
    fn scattered_cells_recover() {
        // 8x8 data → 16x16 grid; erase 60% of cells in a scattered
        // pattern (no full column survives, rows damaged too)
        let (params, z, public) = grid(8, 8, b"recon2-scattered-000000000000000");
        let mut rng = ZodaRng::from_seed(*b"recon2-scatter-rng-0000000000000");
        let mut erased = Vec::new();
        for r in 0..16 {
            for c in 0..16 {
                if rng.next_u64() % 10 < 6 {
                    erased.push((r, c));
                }
            }
        }
        let mut partial = PartialGrid::with_erasures(&z, &erased);
        let stats = reconstruct_2d(&mut partial, &params, &public).expect("decode");
        assert!(stats.cells_recovered > 0);
        let recovered = partial.to_matrix();
        for r in 0..16 {
            for c in 0..16 {
                assert_eq!(recovered.get(r, c), z.get(r, c), "({},{})", r, c);
            }
        }
    }

    #[test]
    fn rectangle_withholding_stalls_and_guides() {
        // Hide an (m+1)x(k+1) rectangle — the minimal ambiguous set.
        // The fixpoint must stall with fetch guidance, never guess wrong.
        let (params, z, public) = grid(8, 8, b"recon2-rectangle-000000000000000");
        let mut erased = Vec::new();
        for r in 0..9 {
            for c in 0..9 {
                erased.push((r, c));
            }
        }
        let mut partial = PartialGrid::with_erasures(&z, &erased);
        match reconstruct_2d(&mut partial, &params, &public) {
            Err(ReconError::Incomplete { stats, fetch_next }) => {
                // the ambiguous rectangle: no line reaches its threshold,
                // so the fixpoint must recover NOTHING and report guidance
                assert_eq!(stats.cells_recovered, 0);
                assert_eq!(stats.rounds, 1);
                assert_eq!(partial.present_count(), 16 * 16 - 81);
                assert!(!fetch_next.is_empty());
                // every suggested cell is genuinely missing
                for &(r, c) in &fetch_next {
                    assert!(!partial.is_present(r, c));
                }
            }
            other => panic!("expected Incomplete, got {:?}", other.is_ok()),
        }
    }

    #[test]
    fn tampered_cell_detected() {
        let (params, z, public) = grid(8, 8, b"recon2-tamper-000000000000000000");
        // keep 12 of 16 cells in every row, tamper one
        let mut partial = PartialGrid::new(16, 16);
        for r in 0..16 {
            for c in 0..16 {
                if c < 12 {
                    partial.insert_cell(r, c, z.get(r, c));
                }
            }
        }
        partial.insert_cell(4, 11, z.get(4, 11) + Fr::ONE);
        match reconstruct_2d(&mut partial, &params, &public) {
            Err(ReconError::CorruptRow(4)) => {}
            other => panic!("expected CorruptRow(4), got ok={:?}", other.is_ok()),
        }
    }

    #[test]
    fn consistent_but_wrong_decoding_caught_by_projections() {
        // Give every row EXACTLY k honest cells of a DIFFERENT codeword:
        // rows decode "consistently" but the projection check must fail.
        let (params, z, public) = grid(8, 8, b"recon2-forged-000000000000000000");
        let mut forged = z.clone();
        // swap two full rows — still a valid tensor codeword of the WRONG
        // grid (row permutation breaks row projections at swapped indices)
        let r0: Vec<Fr> = (0..16).map(|c| z.get(3, c)).collect();
        let r1: Vec<Fr> = (0..16).map(|c| z.get(5, c)).collect();
        for c in 0..16 {
            forged.set(3, c, r1[c]);
            forged.set(5, c, r0[c]);
        }
        let mut partial = PartialGrid::new(16, 16);
        for r in 0..16 {
            for c in 0..8 {
                partial.insert_cell(r, c, forged.get(r, c)); // only data half
            }
        }
        match reconstruct_2d(&mut partial, &params, &public) {
            Err(ReconError::ProjectionMismatch(_)) => {}
            other => panic!("expected ProjectionMismatch, got ok={:?}", other.is_ok()),
        }
    }

    #[test]
    fn codeword_verification() {
        let (params, z, public) = grid(8, 8, b"recon2-verify-000000000000000000");
        assert!(verify_codeword(&z, &params, &public).is_ok());
        let mut bad = z.clone();
        bad.set(7, 7, bad.get(7, 7) + Fr::ONE);
        assert!(verify_codeword(&bad, &params, &public).is_err());
    }

    #[test]
    fn full_columns_path_equivalence() {
        // The classic scenario (k full columns) through the generic engine
        let (params, z, public) = grid(16, 16, b"recon2-cols-00000000000000000000");
        let mut partial = PartialGrid::new(32, 32);
        for c in 0..16 {
            partial.insert_column(c, &z.col(c));
        }
        let stats = reconstruct_2d(&mut partial, &params, &public).expect("decode");
        assert!(stats.rows_decoded >= 32);
        assert_eq!(partial.to_matrix(), z);
    }
}
