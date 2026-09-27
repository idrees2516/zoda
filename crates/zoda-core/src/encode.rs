//! Matrix container and the NTT-based systematic tensor encoder.

use zoda_math::{FftDomain, Fr, PrimeField};

/// A row-major matrix over `Fr`.
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<Fr>,
}

impl Matrix {
    pub fn from_row_major(rows: usize, cols: usize, data: Vec<Fr>) -> Matrix {
        assert_eq!(data.len(), rows * cols);
        Matrix { rows, cols, data }
    }

    pub fn zeros(rows: usize, cols: usize) -> Matrix {
        Matrix {
            rows,
            cols,
            data: vec![Fr::zero(); rows * cols],
        }
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
    pub fn get(&self, r: usize, c: usize) -> Fr {
        self.data[r * self.cols + c]
    }

    #[inline]
    pub fn set(&mut self, r: usize, c: usize, v: Fr) {
        self.data[r * self.cols + c] = v;
    }

    pub fn row(&self, r: usize) -> Vec<Fr> {
        self.data[r * self.cols..(r + 1) * self.cols].to_vec()
    }

    pub fn col(&self, c: usize) -> Vec<Fr> {
        (0..self.rows).map(|r| self.data[r * self.cols + c]).collect()
    }

    pub fn data(&self) -> &[Fr] {
        &self.data
    }

    /// Serialise (for Merkle commitments): each element as 32 LE bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.data.len() * 32);
        for x in &self.data {
            out.extend_from_slice(&x.to_le_bytes());
        }
        out
    }
}

/// Encode one length-n vector into its length-2n systematic RS codeword on
/// the roots-of-unity domain: interpolate the n data points
/// `(ω^0..ω^(n-1), data)` and evaluate on the full 2n-domain.
///
/// O(n log² n) via the subproduct tree — replacing the reference
/// prototype's O(n²) Vandermonde multiplication.
pub fn rs_encode_vector(data: &[Fr], domain: &FftDomain<Fr>) -> Vec<Fr> {
    let n = data.len();
    debug_assert!(domain.size() == 2 * n);
    let roots = &domain.roots_of_unity()[..n];
    let coeffs = zoda_math::poly::interpolate(&roots.to_vec(), &data.to_vec());
    let mut out = vec![Fr::zero(); 2 * n];
    domain.fft_padded(&coeffs, &mut out);
    out
}

/// Extend all rows: `m × k → m × 2k`.
pub fn extend_rows(m: &Matrix, domain: &FftDomain<Fr>) -> Matrix {
    assert_eq!(domain.size(), 2 * m.cols());
    let mut out = Matrix::zeros(m.rows(), 2 * m.cols());
    for r in 0..m.rows() {
        let coded = rs_encode_vector(&m.row(r), domain);
        for (c, v) in coded.iter().enumerate() {
            out.set(r, c, *v);
        }
    }
    out
}

/// Extend all columns: `m × k → 2m × k` (transpose → row-extend →
/// transpose).
pub fn extend_columns(m: &Matrix, domain: &FftDomain<Fr>) -> Matrix {
    assert_eq!(domain.size(), 2 * m.rows());
    let mut out = Matrix::zeros(2 * m.rows(), m.cols());
    for c in 0..m.cols() {
        let coded = rs_encode_vector(&m.col(c), domain);
        for (r, v) in coded.iter().enumerate() {
            out.set(r, c, *v);
        }
    }
    out
}

/// Full tensor encoding: `m × k → 2m × 2k`.
///
/// Row-extend then column-extend: the tensor-code property (all rows are
/// row-codewords AND all columns are column-codewords) follows from the
/// linearity and shared code family of the two passes.
pub fn tensor_encode(data: &Matrix, row_domain: &FftDomain<Fr>, col_domain: &FftDomain<Fr>) -> Matrix {
    let row_ext = extend_rows(data, row_domain);
    extend_columns(&row_ext, col_domain)
}

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_math::ZodaRng;

    #[test]
    fn rs_encoding_is_systematic() {
        let mut rng = ZodaRng::from_seed(*b"rs-enc-test-seed-000000000000000");
        for log in 1..=4u32 {
            let n = 1usize << log;
            let domain = FftDomain::<Fr>::new(2 * n);
            let data: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
            let coded = rs_encode_vector(&data, &domain);
            assert_eq!(coded.len(), 2 * n);
            // systematic prefix
            for i in 0..n {
                assert_eq!(coded[i], data[i], "not systematic n={} i={}", n, i);
            }
        }
    }

    #[test]
    fn tensor_code_property() {
        // After tensor encoding, EVERY row lies in the column-code and
        // every column lies in the row-code. Check by re-encoding:
        // re-encoding a codeword row (as data) reproduces it.
        let mut rng = ZodaRng::from_seed(*b"tensor-enc-test-seed-00000000000");
        let n = 8;
        let data: Vec<Fr> = (0..n * n).map(|_| rng.next_fr(false)).collect();
        let m = Matrix::from_row_major(n, n, data);
        let d_row = FftDomain::<Fr>::new(2 * n);
        let d_col = FftDomain::<Fr>::new(2 * n);
        let z = tensor_encode(&m, &d_row, &d_col);
        assert_eq!(z.rows(), 2 * n);
        assert_eq!(z.cols(), 2 * n);
        // every row re-encodes to itself (col-code membership)
        for r in 0..2 * n {
            let row = z.row(r);
            let re = rs_encode_vector(&row[..n].to_vec(), &d_row);
            assert_eq!(&re[..], &row[..], "row {} not a codeword", r);
        }
        // every column re-encodes to itself (row-code membership)
        for c in 0..2 * n {
            let col = z.col(c);
            let re = rs_encode_vector(&col[..n].to_vec(), &d_col);
            assert_eq!(&re[..], &col[..], "col {} not a codeword", c);
        }
    }
}
