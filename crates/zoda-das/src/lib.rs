//! # zoda-das — data availability sampling over the ZODA tensor code
//!
//! * [`SamplingPlan`] — how many row/column samples achieve a target
//!   detection probability against a withholding adversary,
//! * [`sample_session`] — run a sampling session against a peer oracle,
//!   verifying every sample with the zero-overhead checks,
//! * reconstruction triggers once enough columns are available.

pub use zoda_core::{Matrix, ZodaParams, ZodaPublic};
use zoda_math::stats;
use zoda_math::{Fr, ZodaRng};

/// A sampling plan for one grid.
#[derive(Clone, Debug)]
pub struct SamplingPlan {
    pub rows: usize,
    pub columns: usize,
    /// probability that a fully-withheld dimension escapes detection
    pub miss_probability: f64,
}

impl SamplingPlan {
    /// Plan `r` row samples and `c` column samples for a grid with
    /// `ext_rows` rows and `ext_cols` columns.
    pub fn new(rows: usize, columns: usize, ext_rows: usize, ext_cols: usize) -> SamplingPlan {
        let p = stats::two_dimensional_detection_prob(rows as u64, ext_rows as u64, columns as u64, ext_cols as u64);
        SamplingPlan {
            rows,
            columns,
            miss_probability: 1.0 - p,
        }
    }

    /// The smallest plan achieving miss probability ≤ target for both
    /// dimensions independently.
    pub fn for_target(target: f64, ext_rows: usize, ext_cols: usize) -> SamplingPlan {
        let rows = stats::samples_for_detection(1.0 / ext_rows as f64, 1, target);
        let columns = stats::samples_for_detection(1.0 / ext_cols as f64, 1, target);
        SamplingPlan::new(rows as usize, columns as usize, ext_rows, ext_cols)
    }
}

/// The oracle a sampler consults: given a row or column index, return the
/// full vector (or None if the peer refuses / does not have it).
pub trait SampleOracle {
    fn fetch_row(&mut self, index: usize) -> Option<Vec<Fr>>;
    fn fetch_column(&mut self, index: usize) -> Option<Vec<Fr>>;
}

/// The deterministic full-data oracle (tests, honest nodes).
pub struct FullOracle<'a> {
    pub matrix: &'a Matrix,
}

impl<'a> SampleOracle for FullOracle<'a> {
    fn fetch_row(&mut self, index: usize) -> Option<Vec<Fr>> {
        Some(self.matrix.row(index))
    }
    fn fetch_column(&mut self, index: usize) -> Option<Vec<Fr>> {
        Some(self.matrix.col(index))
    }
}

/// An adversarial oracle that withholds a set of rows and columns.
pub struct WithholdingOracle<'a> {
    pub matrix: &'a Matrix,
    pub hidden_rows: Vec<usize>,
    pub hidden_cols: Vec<usize>,
}

impl<'a> SampleOracle for WithholdingOracle<'a> {
    fn fetch_row(&mut self, index: usize) -> Option<Vec<Fr>> {
        if self.hidden_rows.contains(&index) {
            None
        } else {
            Some(self.matrix.row(index))
        }
    }
    fn fetch_column(&mut self, index: usize) -> Option<Vec<Fr>> {
        if self.hidden_cols.contains(&index) {
            None
        } else {
            Some(self.matrix.col(index))
        }
    }
}

/// Outcome of a sampling session.
#[derive(Clone, Debug)]
pub struct SessionReport {
    pub rows_sampled: usize,
    pub columns_sampled: usize,
    pub rows_verified: usize,
    pub columns_verified: usize,
    pub rows_rejected: usize,
    pub columns_rejected: usize,
    pub missing: usize,
    /// empirical upper bound on undetected withholding probability
    pub confidence: f64,
}

impl SessionReport {
    pub fn available(&self) -> bool {
        self.rows_rejected == 0 && self.columns_rejected == 0
    }
}

/// Run a sampling session: draw the plan's indices (seeded), fetch, and
/// verify each sample with the ZODA zero-overhead checks.
pub fn sample_session(
    public: &ZodaPublic,
    oracle: &mut dyn SampleOracle,
    plan: &SamplingPlan,
    rng: &mut ZodaRng,
) -> SessionReport {
    let mut report = SessionReport {
        rows_sampled: 0,
        columns_sampled: 0,
        rows_verified: 0,
        columns_verified: 0,
        rows_rejected: 0,
        columns_rejected: 0,
        missing: 0,
        confidence: 1.0,
    };
    let ext_rows = public.z_r.len();
    let ext_cols = public.z_r2.len();
    // Distinct draws within a session (a session never asks the same
    // peer for the same row/column twice).
    let row_draws: Vec<usize> =
        rng.shuffle_indices(ext_rows)[..plan.rows.min(ext_rows)].to_vec();
    let col_draws: Vec<usize> =
        rng.shuffle_indices(ext_cols)[..plan.columns.min(ext_cols)].to_vec();
    for r in row_draws {
        report.rows_sampled += 1;
        match oracle.fetch_row(r) {
            Some(row) => match zoda_core::verify_row_sample(public, r, &row) {
                Ok(()) => report.rows_verified += 1,
                Err(_) => report.rows_rejected += 1,
            },
            None => report.missing += 1,
        }
    }
    for c in col_draws {
        report.columns_sampled += 1;
        match oracle.fetch_column(c) {
            Some(col) => match zoda_core::verify_column_sample(public, c, &col) {
                Ok(()) => report.columns_verified += 1,
                Err(_) => report.columns_rejected += 1,
            },
            None => report.missing += 1,
        }
    }
    // Miss bound for distinct (without-replacement) draws: a specific
    // hidden row escapes detection with probability
    // P(miss) = prod_{i=0}^{k-1} (n-1-i)/(n-i) = C(n-1,k)/C(n,k).
    let miss_prob = |n: usize, k: usize| -> f64 {
        if k >= n {
            return 0.0;
        }
        let mut p = 1.0f64;
        for i in 0..k {
            p *= (n as f64 - 1.0 - i as f64) / (n as f64 - i as f64);
        }
        p
    };
    let p_row_miss = miss_prob(ext_rows, report.rows_sampled);
    let p_col_miss = miss_prob(ext_cols, report.columns_sampled);
    report.confidence = 1.0 - p_row_miss * p_col_miss;
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_core::{ZodaParams, ZodaProver};
    use zoda_math::{Fr, PrimeField};

    fn setup_grid(m: usize, k: usize) -> (ZodaParams, ZodaProver, ZodaPublic) {
        let mut rng = ZodaRng::from_seed(*b"das-test-seed-000000000000000000");
        let mut data = vec![Fr::zero(); m * k];
        for x in data.iter_mut() {
            *x = rng.next_fr(false);
        }
        let grid = Matrix::from_row_major(m, k, data);
        let params = ZodaParams::new(m, k);
        let prover = params.commit(&grid);
        let public = prover.public_params();
        (params, prover, public)
    }

    #[test]
    fn honest_session_accepts() {
        let (_, prover, public) = setup_grid(8, 8);
        let mut oracle = FullOracle {
            matrix: &prover.matrix,
        };
        let mut rng = ZodaRng::from_seed(*b"das-session-seed-000000000000000");
        let plan = SamplingPlan::new(6, 6, 16, 16);
        let report = sample_session(&public, &mut oracle, &plan, &mut rng);
        assert!(report.available());
        assert_eq!(report.rows_verified, 6);
        assert_eq!(report.columns_verified, 6);
    }

    #[test]
    fn corrupted_session_rejects() {
        let (_, prover, public) = setup_grid(8, 8);
        // corrupt one row of the matrix the oracle serves
        let mut bad = prover.matrix.clone();
        bad.set(3, 3, bad.get(3, 3) + Fr::ONE);
        let mut oracle = FullOracle { matrix: &bad };
        let mut rng = ZodaRng::from_seed(*b"das-bad-seed-0000000000000000000");
        let plan = SamplingPlan::new(16, 0, 16, 16);
        let report = sample_session(&public, &mut oracle, &plan, &mut rng);
        // sampling all 16 rows must hit the bad row
        assert!(report.rows_rejected >= 1, "corruption not detected");
    }

    #[test]
    fn withholding_detected() {
        let (_, prover, public) = setup_grid(8, 8);
        let mut oracle = WithholdingOracle {
            matrix: &prover.matrix,
            hidden_rows: vec![5],
            hidden_cols: vec![],
        };
        let mut rng = ZodaRng::from_seed(*b"das-hide-seed-000000000000000000");
        let plan = SamplingPlan::new(16, 16, 16, 16);
        let report = sample_session(&public, &mut oracle, &plan, &mut rng);
        assert_eq!(report.missing, 1);
        // with all rows sampled the hidden row is definitely noticed
        assert!(report.confidence > 0.99);
    }

    #[test]
    fn plan_targets() {
        let plan = SamplingPlan::for_target(1e-9, 128, 128);
        assert!(plan.miss_probability <= 1e-9);
    }
}
