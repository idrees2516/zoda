//! Exact data-availability theory for 2D grids.
//!
//! ## The minimal withholding set of a tensor code
//!
//! The `2m × 2k` extended grid is a product code: rows are length-`2k`
//! Reed–Solomon codewords of degree `< k` (minimum distance `k+1`), columns
//! are length-`2m` codewords of degree `< m` (minimum distance `m+1`).
//! The tensor code's minimum distance is the product `(m+1)(k+1)`.
//!
//! **Uniqueness.** The available cells determine the codeword uniquely iff
//! they do not contain — equivalently, the hidden set does not *include*
//! the support of — a nonzero codeword. Every hidden set smaller than
//! `(m+1)(k+1)` cells is therefore unambiguous.
//!
//! **Achievability.** Minimum-weight codewords of a systematic RS code on
//! the roots-of-unity domain vanish at *any* chosen `k` of the `2k`
//! evaluation points (take the degree-`k` polynomial with those roots —
//! wait, degree `k` exceeds the code's degree bound, so use the degree-`k−1`
//! polynomial vanishing on a chosen `k−1`-subset: weight exactly `k+1`).
//! Their outer products are tensor codewords of weight `(m+1)(k+1)` whose
//! supports are combinatorial rectangles. Hence the adversary's optimal
//! strategy is to hide an **(m+1) × (k+1) rectangle** — nothing smaller
//! can ever be ambiguous, and this always is.
//!
//! ## Sampling bounds
//!
//! * Cell sampling (`s` uniform distinct cells of the `4mk` grid):
//!   the rectangle escapes iff every sample lands outside it —
//!   `P = C(4mk − (m+1)(k+1), s) / C(4mk, s)` (hypergeometric, computed
//!   in log space).
//! * Line sampling (ZODA: whole rows/columns): hiding the rectangle means
//!   refusing `m+1` rows AND `k+1` columns, so the escape probability of
//!   `u_r` distinct row draws and `u_c` distinct column draws is the
//!   product of the two hypergeometric tails.
//! * EIP-7594 blob grids: each blob is an independent length-128 cell
//!   row of a degree-`<64` code (in 64-eval cell units); ambiguity needs
//!   ≥ 65 of 128 cells hidden, and column sampling detects it unless
//!   every sampled column avoids the hidden columns.

use zoda_math::stats::log_binomial;

/// Cells an adversary must hide to make a `2m × 2k` tensor grid
/// ambiguous — an `(m+1) × (k+1)` rectangle.
#[inline]
pub fn min_withholding_cells(m: usize, k: usize) -> usize {
    (m + 1) * (k + 1)
}

/// Escape probability of the minimal rectangle against `s` uniform
/// DISTINCT cell samples (exact, log-space evaluation).
pub fn cell_miss_probability(s: usize, m: usize, k: usize) -> f64 {
    let n = 4 * m * k;
    let h = min_withholding_cells(m, k);
    if s == 0 {
        return 1.0;
    }
    if h >= n || s > n - h {
        return 0.0;
    }
    (log_binomial((n - h) as u64, s as u64) - log_binomial(n as u64, s as u64)).exp()
}

/// Escape probability against `s` uniform cell samples WITH replacement:
/// `(1 − H/N)^s`. Slightly more conservative than the distinct variant.
pub fn cell_miss_probability_replacement(s: usize, m: usize, k: usize) -> f64 {
    let n = 4 * m * k;
    if n == 0 {
        return 1.0;
    }
    let h = min_withholding_cells(m, k) as f64;
    (1.0 - h / n as f64).powi(s as i32)
}

/// Smallest number of distinct cell samples achieving
/// `miss probability ≤ target` (binary search over the exact formula).
pub fn cell_samples_for_target(target: f64, m: usize, k: usize) -> usize {
    let n = 4 * m * k;
    let mut lo = 0usize;
    let mut hi = n;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if cell_miss_probability(mid, m, k) <= target {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// Line-level escape probability: the adversary refuses `m+1` rows and
/// `k+1` columns; `u_r` distinct row draws and `u_c` distinct column
/// draws all miss their hidden sets.
pub fn line_miss_probability(u_r: usize, u_c: usize, m: usize, k: usize) -> f64 {
    // P(all u_r draws avoid a fixed (m+1)-set of 2m rows)
    let rows = if u_r == 0 {
        1.0
    } else {
        let n = 2 * m as u64;
        let hidden = (m + 1) as u64;
        if u_r as u64 > n - hidden {
            0.0
        } else {
            (log_binomial(n - hidden, u_r as u64) - log_binomial(n, u_r as u64)).exp()
        }
    };
    let cols = if u_c == 0 {
        1.0
    } else {
        let n = 2 * k as u64;
        let hidden = (k + 1) as u64;
        if u_c as u64 > n - hidden {
            0.0
        } else {
            (log_binomial(n - hidden, u_c as u64) - log_binomial(n, u_c as u64)).exp()
        }
    };
    rows * cols
}

// ---------------------------------------------------------------------------
// EIP-7594 blob grids
// ---------------------------------------------------------------------------

/// Cells of one extended blob (128) that an adversary must hide to make
/// the blob's polynomial ambiguous: 65 cells = 65·64 = 4160 > 4095
/// evaluations, strictly more than the 4096-evaluation code can pin down.
pub const BLOB_MIN_HIDDEN_CELLS: usize = 65;
/// Extended cells per blob (`CELLS_PER_EXT_BLOB` from the spec).
pub const BLOB_CELLS: usize = 128;

/// Escape probability of one withheld blob against `s` distinct sampled
/// COLUMNS (Fulu-style sampling: a sampled column serves every blob's cell
/// in that column). The blob is detected unless all `s` sampled columns
/// avoid the `hidden` hidden ones.
pub fn blob_column_miss_probability(s: usize, hidden: usize, columns: usize) -> f64 {
    if s == 0 {
        return 1.0;
    }
    if hidden >= columns || s > columns - hidden {
        return 0.0;
    }
    (
        log_binomial((columns - hidden) as u64, s as u64)
            - log_binomial(columns as u64, s as u64)
    )
        .exp()
}

/// Smallest column-sample count with per-blob miss probability ≤ target.
pub fn blob_column_samples_for_target(target: f64, hidden: usize, columns: usize) -> usize {
    let mut lo = 0usize;
    let mut hi = columns;
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if blob_column_miss_probability(mid, hidden, columns) <= target {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangle_size() {
        assert_eq!(min_withholding_cells(128, 128), 129 * 129);
        assert_eq!(min_withholding_cells(8, 8), 81);
    }

    #[test]
    fn cell_miss_decays_geometrically() {
        // square grid: rectangle covers ~1/4 of the cells, so the miss
        // probability decays like (3/4)^s
        let m = 128;
        let k = 128;
        let p75 = cell_miss_probability(1, m, k);
        let exact = (4.0 * (m * k) as f64 - 129.0 * 129.0) / (4.0 * (m * k) as f64);
        assert!((p75 - exact).abs() < 1e-12);
        let p10 = cell_miss_probability(10, m, k);
        let p20 = cell_miss_probability(20, m, k);
        assert!(p10 > 0.02 && p10 < 0.08, "p10 = {}", p10);
        assert!(p20 < p10 * p10 + 1e-9); // roughly squaring each doubling
        // 2^-40 target on a 128x128 grid: (1 - 16641/65536)^s ≤ 2^-40
        // ⇒ s = 95
        let s = cell_samples_for_target(2f64.powi(-40), m, k);
        assert!(s > 80 && s < 110, "s = {}", s);
        assert!(cell_miss_probability(s, m, k) <= 2f64.powi(-40));
        assert!(cell_miss_probability(s - 1, m, k) > 2f64.powi(-40));
    }

    #[test]
    fn replacement_variant_is_conservative() {
        for s in [1usize, 5, 20, 60] {
            assert!(
                cell_miss_probability_replacement(s, 64, 64)
                    >= cell_miss_probability(s, 64, 64) - 1e-15
            );
        }
    }

    #[test]
    fn line_miss_matches_closed_form() {
        // u_r draws avoiding a fixed (m+1)-set of 2m: closed form
        // Π_{i<u_r} (m−1−i)/(2m−i)
        let m = 8usize;
        let k = 8usize;
        let u = 5usize;
        let mut closed = 1.0f64;
        for i in 0..u {
            closed *= (2.0 * m as f64 - (m + 1) as f64 - i as f64) / (2.0 * m as f64 - i as f64);
        }
        let cols_closed = {
            let mut c = 1.0f64;
            for i in 0..u {
                c *= (2.0 * k as f64 - (k + 1) as f64 - i as f64) / (2.0 * k as f64 - i as f64);
            }
            c
        };
        let got = line_miss_probability(u, u, m, k);
        assert!((got - closed * cols_closed).abs() < 1e-12 * got.max(1e-300));
        // certainty: sampling m rows covers every m+1 hidden set possible?
        // (u_r = m draws from 2m must hit any fixed (m+1)-set)
        assert_eq!(line_miss_probability(m, 0, m, k), 0.0);
    }

    #[test]
    fn blob_grid_bounds() {
        // 8 column samples vs a 65-column withholding: escape is
        // C(63,8)/C(128,8) = 0.0027088812... ≈ 2^-8.5 per blob
        let p = blob_column_miss_probability(8, BLOB_MIN_HIDDEN_CELLS, BLOB_CELLS);
        assert!((p - 0.0027088812424).abs() < 1e-9, "p = {}", p);
        assert!(p > 2f64.powi(-9) && p < 2f64.powi(-8), "p = {}", p);
        // more samples shrink it fast
        assert!(blob_column_miss_probability(16, 65, 128) < p * p);
        // target machinery
        let s = blob_column_samples_for_target(2f64.powi(-20), 65, 128);
        assert!(blob_column_miss_probability(s, 65, 128) <= 2f64.powi(-20));
        assert!(blob_column_miss_probability(s - 1, 65, 128) > 2f64.powi(-20));
    }
}
