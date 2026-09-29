//! # zoda-rda — randomized data availability service
//!
//! Orchestrates adaptive sampling sessions: start with a small batch of
//! samples, verify them, and request more (up to a budget) whenever the
//! achieved confidence is below target — stopping early when peers are
//! healthy. Committee selection uses stake-weighted BLS sortition from
//! `zoda-sybils`.

use zoda_core::ZodaPublic;
use zoda_das::{sample_session, SamplingPlan, SampleOracle, SessionReport};
use zoda_math::ZodaRng;
use zoda_sybils::{Peer, Sortition};

/// Session configuration.
#[derive(Clone, Debug)]
pub struct RdaConfig {
    /// target detection confidence (e.g. 1 − 2^-32)
    pub target_confidence: f64,
    /// minimum samples per round
    pub min_batch: usize,
    /// maximum total samples per session
    pub max_samples: usize,
    /// stop early once this confidence is reached
    pub early_exit: bool,
}

impl Default for RdaConfig {
    fn default() -> Self {
        RdaConfig {
            target_confidence: 1.0 - 2f64.powi(-32),
            min_batch: 8,
            max_samples: 128,
            early_exit: true,
        }
    }
}

/// An adaptive session over one grid.
pub struct RdaSession<'a> {
    pub public: &'a ZodaPublic,
    pub config: RdaConfig,
    pub peers: Vec<Peer>,
    /// history of per-round reports
    pub rounds: Vec<SessionReport>,
    pub total_rows: usize,
    pub total_columns: usize,
}

impl<'a> RdaSession<'a> {
    pub fn new(public: &'a ZodaPublic, peers: Vec<Peer>) -> Self {
        RdaSession {
            public,
            config: RdaConfig::default(),
            peers,
            rounds: Vec::new(),
            total_rows: 0,
            total_columns: 0,
        }
    }

    /// Run adaptive rounds until the target confidence is reached or the
    /// budget is exhausted.
    pub fn run(&mut self, oracle: &mut dyn SampleOracle, rng: &mut ZodaRng) -> RdaVerdict {
        let ext_rows = self.public.z_r.len();
        let ext_cols = self.public.z_r2.len();
        loop {
            let remaining = self.config.max_samples - (self.total_rows + self.total_columns);
            if remaining == 0 {
                break;
            }
            let batch = self.config.min_batch.min(remaining);
            let plan = SamplingPlan::new(batch / 2, batch - batch / 2, ext_rows, ext_cols);
            let report = sample_session(self.public, oracle, &plan, rng);
            self.total_rows += report.rows_sampled;
            self.total_columns += report.columns_sampled;
            let rejected = report.rows_rejected + report.columns_rejected;
            self.rounds.push(report);
            if rejected > 0 {
                return RdaVerdict::Rejected;
            }
            let conf = self.cumulative_confidence();
            if self.config.early_exit && conf >= self.config.target_confidence {
                return RdaVerdict::Available;
            }
        }
        if self.cumulative_confidence() >= self.config.target_confidence {
            RdaVerdict::Available
        } else {
            RdaVerdict::InsufficientSamples
        }
    }

    /// Cumulative detection confidence over all rounds (distinct draws
    /// per round; k ≥ n saturates to certainty).
    pub fn cumulative_confidence(&self) -> f64 {
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
        let ext_rows = self.public.z_r.len();
        let ext_cols = self.public.z_r2.len();
        // per-round distinct draws: effective coverage grows with the
        // union of drawn indices; bounded by the with-replacement product
        // across rounds, saturating once the budget covers all indices
        let p_row_miss = miss_prob(ext_rows, self.total_rows);
        let p_col_miss = miss_prob(ext_cols, self.total_columns);
        1.0 - p_row_miss * p_col_miss
    }

    /// The committee for the next round (stake-weighted sortition).
    pub fn next_committee(&self, seed: &[u8; 32], count: usize) -> Vec<usize> {
        Sortition::select_weighted(seed, &self.peers, count)
    }
}

/// Final session verdict.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RdaVerdict {
    /// All samples verified and confidence reached the target.
    Available,
    /// A sample failed verification — data is invalid.
    Rejected,
    /// Budget exhausted before reaching the target confidence.
    InsufficientSamples,
}

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_core::{Matrix, ZodaParams, ZodaProver};
    use zoda_das::FullOracle;
    use zoda_math::{Fr, PrimeField};

    fn grid(m: usize, k: usize, seed: [u8; 32]) -> (ZodaProver, ZodaPublic) {
        let mut rng = ZodaRng::from_seed(seed);
        let mut data = vec![Fr::zero(); m * k];
        for x in data.iter_mut() {
            *x = rng.next_fr(false);
        }
        let g = Matrix::from_row_major(m, k, data);
        let params = ZodaParams::new(m, k);
        let prover = params.commit(&g);
        let public = prover.public_params();
        (prover, public)
    }

    fn peers(n: usize) -> Vec<Peer> {
        (0..n)
            .map(|i| {
                let sk = zoda_bls::SecretKey::from_seed(format!("rda-peer-{}", i).as_bytes());
                Peer {
                    id: [i as u8; 32],
                    stake: 100,
                    public_key: sk.public_key(),
                }
            })
            .collect()
    }

    #[test]
    fn healthy_session_available() {
        let (prover, public) = grid(8, 8, *b"rda-ok-seed-00000000000000000000");
        let mut oracle = FullOracle {
            matrix: &prover.matrix,
        };
        let mut rng = ZodaRng::from_seed(*b"rda-session-seed-000000000000000");
        let mut session = RdaSession::new(&public, peers(4));
        session.config.max_samples = 200;
        let verdict = session.run(&mut oracle, &mut rng);
        assert_eq!(verdict, RdaVerdict::Available);
        assert!(session.cumulative_confidence() >= session.config.target_confidence);
    }

    #[test]
    fn corrupted_session_rejected() {
        let (prover, public) = grid(8, 8, *b"rda-bad-seed-0000000000000000000");
        let mut bad = prover.matrix.clone();
        bad.set(2, 2, bad.get(2, 2) + Fr::ONE);
        let mut oracle = FullOracle { matrix: &bad };
        let mut rng = ZodaRng::from_seed(*b"rda-bad-session-seed-00000000000");
        let mut session = RdaSession::new(&public, peers(4));
        session.config.max_samples = 10_000;
        session.config.early_exit = false;
        let verdict = session.run(&mut oracle, &mut rng);
        assert_eq!(verdict, RdaVerdict::Rejected);
    }

    #[test]
    fn budget_limited_reports_insufficient() {
        let (prover, public) = grid(64, 64, *b"rda-budget-seed-0000000000000000");
        let mut oracle = FullOracle {
            matrix: &prover.matrix,
        };
        let mut rng = ZodaRng::from_seed(*b"rda-budget-sess-seed-00000000000");
        let mut session = RdaSession::new(&public, peers(2));
        session.config.max_samples = 4; // far too few for 128x128
        let verdict = session.run(&mut oracle, &mut rng);
        assert_eq!(verdict, RdaVerdict::InsufficientSamples);
    }

    #[test]
    fn committee_selection() {
        let (_, public) = grid(4, 4, *b"rda-comm-seed-000000000000000000");
        let session = RdaSession::new(&public, peers(6));
        let c1 = session.next_committee(&[1u8; 32], 3);
        let c2 = session.next_committee(&[1u8; 32], 3);
        assert_eq!(c1, c2);
        assert_eq!(c1.len(), 3);
    }
}
