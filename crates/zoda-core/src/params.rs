//! ZODA protocol parameters, commitments and the Fiat–Shamir derivation
//! of the random projections.

use crate::encode::{tensor_encode, Matrix};
use zoda_math::merkle::MerkleTree;
use zoda_math::sha256::sha256;
use zoda_math::{FftDomain, Fr, PrimeField, ZodaRng};

/// Encoding marker for row (length-2k) vectors.
pub const ENC_ROW: u8 = 0x01;
/// Encoding marker for column (length-2m) vectors.
pub const ENC_COL: u8 = 0x02;

const ZODA_FS_DOMAIN: &[u8] = b"ZODA-TENSOR-FIAT-SHAMIR-V1";

/// Protocol parameters for one grid shape.
pub struct ZodaParams {
    /// data rows (m)
    pub m: usize,
    /// data columns (k)
    pub k: usize,
    /// FFT domain for rows (size 2k)
    pub row_domain: FftDomain<Fr>,
    /// FFT domain for columns (size 2m)
    pub col_domain: FftDomain<Fr>,
}

impl ZodaParams {
    /// Create parameters for an `m × k` data grid. Both `m` and `k` must
    /// be powers of two.
    pub fn new(m: usize, k: usize) -> ZodaParams {
        assert!(m.is_power_of_two() && k.is_power_of_two());
        ZodaParams {
            m,
            k,
            row_domain: FftDomain::<Fr>::new(2 * k),
            col_domain: FftDomain::<Fr>::new(2 * m),
        }
    }

    pub fn extended_rows(&self) -> usize {
        2 * self.m
    }

    pub fn extended_cols(&self) -> usize {
        2 * self.k
    }

    /// Encode and commit: builds the tensor codeword, the row/column
    /// Merkle trees, and the Fiat–Shamir-bound projections.
    pub fn commit(&self, data: &Matrix) -> ZodaProver {
        assert_eq!(data.rows(), self.m);
        assert_eq!(data.cols(), self.k);
        let matrix = tensor_encode(data, &self.row_domain, &self.col_domain);

        // Row/column Merkle commitments
        let row_leaves: Vec<Vec<u8>> = (0..self.extended_rows())
            .map(|r| row_leaf(&matrix, r))
            .collect();
        let col_leaves: Vec<Vec<u8>> = (0..self.extended_cols())
            .map(|c| col_leaf(&matrix, c))
            .collect();
        let row_tree = MerkleTree::build(&row_leaves);
        let col_tree = MerkleTree::build(&col_leaves);
        let row_root = row_tree.root().expect("nonempty");
        let col_root = col_tree.root().expect("nonempty");

        // Random projections, bound to the commitments (fixes the
        // reference prototype's unbound transcript).
        let (g_r, g_r2, z_r, z_r2) = derive_public(self, &row_root, &col_root);

        ZodaProver {
            matrix,
            row_tree,
            col_tree,
            public: None,
        }
        .with_public(self, g_r, g_r2, z_r, z_r2, row_root, col_root)
    }
}

fn row_leaf(matrix: &Matrix, r: usize) -> Vec<u8> {
    let mut buf = Vec::with_capacity(1 + matrix.cols() * 32);
    buf.push(ENC_ROW);
    buf.extend_from_slice(&(r as u64).to_le_bytes());
    for c in 0..matrix.cols() {
        buf.extend_from_slice(&matrix.get(r, c).to_le_bytes());
    }
    buf
}

fn col_leaf(matrix: &Matrix, c: usize) -> Vec<u8> {
    let mut buf = Vec::with_capacity(1 + matrix.rows() * 32);
    buf.push(ENC_COL);
    buf.extend_from_slice(&(c as u64).to_le_bytes());
    for r in 0..matrix.rows() {
        buf.extend_from_slice(&matrix.get(r, c).to_le_bytes());
    }
    buf
}

/// Derive `g_r`, `g_r2` (challenge vectors) and `z_r`, `z_r2` (their
/// projections through the encoder) from the commitment roots.
fn derive_public(
    params: &ZodaParams,
    row_root: &[u8; 32],
    col_root: &[u8; 32],
) -> (Vec<Fr>, Vec<Fr>, Vec<Fr>, Vec<Fr>) {
    // seed = H(DOMAIN || row_root || col_root || m || k)
    let mut buf = Vec::with_capacity(ZODA_FS_DOMAIN.len() + 64 + 16);
    buf.extend_from_slice(ZODA_FS_DOMAIN);
    buf.extend_from_slice(row_root);
    buf.extend_from_slice(col_root);
    buf.extend_from_slice(&(params.m as u64).to_le_bytes());
    buf.extend_from_slice(&(params.k as u64).to_le_bytes());
    let seed = sha256(&buf);
    let mut rng = ZodaRng::from_seed(seed);

    // g_r: length 2k over the row-code (column extension of data rows)
    let g_r: Vec<Fr> = (0..params.extended_cols()).map(|_| rng.next_fr(false)).collect();
    // g_r2: length 2m over the column-code
    let g_r2: Vec<Fr> = (0..params.extended_rows()).map(|_| rng.next_fr(false)).collect();

    // z_r = the g_r-combination through the COLUMN code applied to the
    // data rows: z_r[r] = Σ_j data[r][j]·(g_r as seen through the
    // systematic generator). Concretely, encode g_r with the row code
    // and read its first m... — following the reference construction:
    // z_r = col_encoding · g_r where col_encoding is the m→2m column
    // extension of the data matrix. We compute it lazily from the
    // committed matrix in the prover instead; here we return the raw
    // challenges and let the prover (who knows the data) compute z_r.
    // For the PUBLIC parameters the prover publishes z_r and z_r2.
    (g_r, g_r2, vec![], vec![])
}

/// Everything the prover keeps (the encoded matrix + trees).
pub struct ZodaProver {
    pub matrix: Matrix,
    pub row_tree: MerkleTree,
    pub col_tree: MerkleTree,
    /// public parameters (present when built via `ZodaParams::commit`)
    pub public: Option<ZodaPublic>,
}

impl ZodaProver {
    /// Compute the projections and assemble the public parameters.
    fn with_public(
        self,
        params: &ZodaParams,
        g_r: Vec<Fr>,
        g_r2: Vec<Fr>,
        _z_r: Vec<Fr>,
        _z_r2: Vec<Fr>,
        row_root: [u8; 32],
        col_root: [u8; 32],
    ) -> ZodaProver {
        // z_r[r] = Σ_j Z[r][j]·g_r[j] — a random projection of row r
        // through the column-code generator (all rows of Z are
        // column-code codewords, so this equals the generator-row
        // projection used by the protocol).
        let z_r: Vec<Fr> = (0..params.extended_rows())
            .map(|r| {
                let mut acc = Fr::zero();
                for (j, &g) in g_r.iter().enumerate() {
                    acc = acc + self.matrix.get(r, j) * g;
                }
                acc
            })
            .collect();
        let z_r2: Vec<Fr> = (0..params.extended_cols())
            .map(|c| {
                let mut acc = Fr::zero();
                for (r, &g) in g_r2.iter().enumerate() {
                    acc = acc + self.matrix.get(r, c) * g;
                }
                acc
            })
            .collect();

        ZodaProver {
            matrix: self.matrix,
            row_tree: self.row_tree,
            col_tree: self.col_tree,
            public: Some(ZodaPublic {
                g_r,
                g_r2,
                z_r,
                z_r2,
                row_root,
                col_root,
            }),
        }
    }

    pub fn public_params(&self) -> ZodaPublic {
        self.public.clone().expect("prover built via ZodaParams::commit")
    }

    /// Merkle proof for row `r`.
    pub fn row_proof(&self, r: usize) -> Option<zoda_math::merkle::MerkleProof> {
        self.row_tree.proof(r)
    }

    /// Merkle proof for column `c`.
    pub fn col_proof(&self, c: usize) -> Option<zoda_math::merkle::MerkleProof> {
        self.col_tree.proof(c)
    }
}

/// The public verification parameters.
#[derive(Clone)]
pub struct ZodaPublic {
    /// Random vector over the column code (length 2k), Fiat–Shamir bound.
    pub g_r: Vec<Fr>,
    /// Random vector over the row code (length 2m).
    pub g_r2: Vec<Fr>,
    /// Row-side projection: `z_r[r] = Σ_j Z[r][j]·g_r[j]` (length 2m).
    pub z_r: Vec<Fr>,
    /// Column-side projection: `z_r2[c] = Σ_r Z[r][c]·g_r2[r]` (length 2k).
    pub z_r2: Vec<Fr>,
    /// Merkle root over all rows.
    pub row_root: [u8; 32],
    /// Merkle root over all columns.
    pub col_root: [u8; 32],
}
