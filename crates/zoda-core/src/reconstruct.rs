//! Erasure reconstruction of the full tensor codeword from ≥ k columns
//! (or ≥ m rows) — the piece the reference prototype left
//! `!unimplemented!()`.

use crate::encode::Matrix;
use crate::params::{ZodaParams, ZodaPublic};
use zoda_math::{FftDomain, Fr, PrimeField};

/// Reconstruct the full `2m × 2k` matrix from a set of available columns
/// (identified by index). At least `k` distinct columns are required.
///
/// Every row of the tensor codeword is a length-`2k` Reed–Solomon codeword
/// of degree `< k`, so any `k` of its symbols determine it: we interpolate
/// each row from the available positions and re-evaluate on the full
/// domain.
pub fn reconstruct(
    public: &ZodaPublic,
    available_cols: &[usize],
    partial: &Matrix,
) -> Result<Matrix, String> {
    // Infer shape from the public projections.
    let height = public.g_r2.len(); // 2m
    let width = public.g_r.len(); // 2k
    let k = width / 2;
    if partial.rows() != height || partial.cols() != width {
        return Err("partial matrix has the wrong shape".to_string());
    }
    if available_cols.len() < k {
        return Err(format!(
            "need at least {} columns to reconstruct, got {}",
            k,
            available_cols.len()
        ));
    }
    let mut seen = vec![false; width];
    for &c in available_cols {
        if c >= width || seen[c] {
            return Err("column indices must be distinct and in range".to_string());
        }
        seen[c] = true;
    }

    // Domain positions of the available columns.
    let domain = FftDomain::<Fr>::new(width);
    let roots = &domain.roots_of_unity()[..width];
    let xs: Vec<Fr> = available_cols.iter().map(|&c| roots[c]).collect();

    // All rows share the same interpolation nodes, so the Lagrange basis
    // is built once (O(k²) via synthetic division of the shared vanishing
    // polynomial + a single batch inversion) and every row reduces to an
    // O(k²) weighted sum — no per-row subproduct tree.
    let k = xs.len();
    let z = zoda_math::poly::vanishing_poly(&xs);
    let zprime = zoda_math::poly::derivative(&z);
    let mut denominators: Vec<Fr> = xs
        .iter()
        .map(|x| zoda_math::poly::eval(&zprime, *x))
        .collect();
    zoda_math::batch_invert(&mut denominators);
    // L_j = Z(X) / (X − x_j), scaled by 1/Z'(x_j)
    let mut basis: Vec<Vec<Fr>> = Vec::with_capacity(k);
    for j in 0..k {
        // synthetic division of Z (degree k) by (X − x_j): exact, since
        // Z(x_j) = 0. Top-down recurrence q_{t−1} = z_t + x_j·q_t.
        let mut qq = vec![Fr::zero(); k];
        let mut carry = Fr::zero();
        for t in (0..k).rev() {
            carry = z[t + 1] + carry * xs[j];
            qq[t] = carry;
        }
        for t in 0..k {
            qq[t] = qq[t] * denominators[j];
        }
        basis.push(qq);
    }

    let mut out = Matrix::zeros(height, width);
    // per row: coeffs = Σ_j y_j · L_j, then evaluate over the full domain
    let rows: Vec<Vec<Fr>> = (0..height)
        .map(|r| {
            let ys: Vec<Fr> = available_cols.iter().map(|&c| partial.get(r, c)).collect();
            let mut coeffs = vec![Fr::zero(); k];
            for j in 0..k {
                let y = ys[j];
                if y.is_zero() {
                    continue;
                }
                for t in 0..k {
                    coeffs[t] = coeffs[t] + y * basis[j][t];
                }
            }
            let mut evals = vec![Fr::zero(); width];
            domain.fft_padded(&coeffs, &mut evals);
            evals
        })
        .collect();
    for (r, evals) in rows.into_iter().enumerate() {
        for (c, v) in evals.iter().enumerate() {
            out.set(r, c, *v);
        }
    }
    Ok(out)
}

/// Reconstruct the original `m × k` data matrix from the recovered full
/// matrix (the top-left quadrant).
pub fn data_from_full(full: &Matrix, params: &ZodaParams) -> Matrix {
    let mut out = Matrix::zeros(params.m, params.k);
    for r in 0..params.m {
        for c in 0..params.k {
            out.set(r, c, full.get(r, c));
        }
    }
    out
}

/// Column-index sets valid for reconstruction (any k distinct columns).
pub fn default_available(width: usize, k: usize) -> Vec<usize> {
    (0..k.min(width)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolates_correctly() {
        // the underlying primitive: k points determine the row
        let mut rng = zoda_math::ZodaRng::from_seed(*b"recon-test-seed-0000000000000000");
        let k = 8;
        let width = 2 * k;
        let domain = FftDomain::<Fr>::new(width);
        let roots = &domain.roots_of_unity()[..width];
        let coeffs: Vec<Fr> = (0..k).map(|_| rng.next_fr(false)).collect();
        let mut evals = vec![Fr::zero(); width];
        domain.fft_padded(&coeffs, &mut evals);
        // keep k of the 2k positions (the "hard" half)
        let keep: Vec<usize> = (k..width).collect();
        let xs: Vec<Fr> = keep.iter().map(|&i| roots[i]).collect();
        let ys: Vec<Fr> = keep.iter().map(|&i| evals[i]).collect();
        let rec = zoda_math::poly::interpolate(&xs, &ys);
        assert_eq!(rec, coeffs);
    }
}
