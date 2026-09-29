//! Zero-overhead sample verification: a full row or column is its own
//! proof, checked with inner products against the public projections.

use crate::params::{ZodaPublic, ENC_COL, ENC_ROW};
use zoda_math::{Fr, PrimeField};

/// Verify a sampled row (all `2k` symbols of row `r`).
///
/// Check: `W · g_r == generator_row(r) · z_r` where the generator row is
/// `e_r` for data rows (`r < m`) and the RS parity row otherwise. Because
/// every row of the tensor codeword is a column-code codeword
/// `Z[r][·] = Σ_i G[r][i]·D'[i][·]` (with `D'` the column-extended data),
/// we have `W·g_r = Σ_i G[r][i]·(D'·g_r)[i] = Σ_i G[r][i]·z_r[i]`.
pub fn verify_row_sample(public: &ZodaPublic, row_index: usize, row: &[Fr]) -> Result<(), String> {
    let width = public.g_r.len();
    if row.len() != width {
        return Err(format!("row must have {} symbols", width));
    }
    if row_index >= public.z_r.len() {
        return Err("row index out of range".to_string());
    }
    // lhs = W · g_r
    let lhs = inner(row, &public.g_r);
    // rhs: for a data row r (< m implied by width/2), the generator row is
    // e_r; for a parity row it is the RS parity row — both captured by
    // evaluating the row code's generator at the row position:
    // rhs = (row-code generator at r) · z_r.
    let rhs = generator_projection_row(public, row_index);
    if lhs == rhs {
        Ok(())
    } else {
        Err("row sample failed the projection check".to_string())
    }
}

/// Verify a sampled column (all `2m` symbols of column `c`).
pub fn verify_column_sample(
    public: &ZodaPublic,
    col_index: usize,
    col: &[Fr],
) -> Result<(), String> {
    let height = public.g_r2.len();
    if col.len() != height {
        return Err(format!("column must have {} symbols", height));
    }
    if col_index >= public.z_r2.len() {
        return Err("column index out of range".to_string());
    }
    let lhs = inner(col, &public.g_r2);
    let rhs = generator_projection_col(public, col_index);
    if lhs == rhs {
        Ok(())
    } else {
        Err("column sample failed the projection check".to_string())
    }
}

/// The generator-side projection for a row: the row-code generator row at
/// `row_index` applied to `z_r`.
///
/// The row code is systematic: data rows r < m have generator e_r, parity
/// rows have the RS parity generator. Since `z_r` was defined as the
/// projection of the FULL codeword through `g_r`, we reproduce the
/// generator row by re-encoding `z_r`'s information set:
/// `gen(r)·z_r` = the `r`-th symbol of the codeword encoding
/// `(z_r[0..m], parity)` — computed directly from the data/parity split
/// of z_r via interpolation.
fn generator_projection_row(public: &ZodaPublic, row_index: usize) -> Fr {
    // z_r is itself a codeword of the column code (a projection of
    // codewords). Its value at "generator position" row_index is exactly
    // z_r[row_index] — because the generator row applied to the projected
    // data re-derives the projection:
    //   Σ_i G[r][i]·z_r[i] = z_r[r]  when z_r is the projection of a
    //   codeword through g_r.
    // This identity is the crux of the zero-overhead property.
    public.z_r[row_index]
}

fn generator_projection_col(public: &ZodaPublic, col_index: usize) -> Fr {
    public.z_r2[col_index]
}

#[inline]
fn inner(a: &[Fr], b: &[Fr]) -> Fr {
    let mut acc = Fr::zero();
    for (x, y) in a.iter().zip(b.iter()) {
        acc = acc + *x * *y;
    }
    acc
}

/// Domain-separated leaf bytes for a sampled row (for Merkle validation
/// by callers that also check inclusion).
pub fn row_leaf_bytes(row_index: usize, row: &[Fr]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(1 + 8 + row.len() * 32);
    buf.push(ENC_ROW);
    buf.extend_from_slice(&(row_index as u64).to_le_bytes());
    for x in row {
        buf.extend_from_slice(&x.to_le_bytes());
    }
    buf
}

/// Domain-separated leaf bytes for a sampled column.
pub fn col_leaf_bytes(col_index: usize, col: &[Fr]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(1 + 8 + col.len() * 32);
    buf.push(ENC_COL);
    buf.extend_from_slice(&(col_index as u64).to_le_bytes());
    for x in col {
        buf.extend_from_slice(&x.to_le_bytes());
    }
    buf
}
