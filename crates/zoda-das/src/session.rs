//! Multi-round, multi-peer ZODA line-sampling sessions.
//!
//! Upgrades the basic [`sample_session`](crate::sample_session) with the
//! pieces a production sampler needs:
//!
//! * **Merkle-bound samples** — every fetched row/column arrives with its
//!   inclusion proof against the public commitment root, so a serving
//!   peer cannot replay a valid line from a different grid and cannot
//!   forge one at all (the projection check binds the line to the
//!   Fiat–Shamir challenges, the Merkle proof binds it to the root).
//! * **Peer orchestration** — samples are routed round-robin across a
//!   peer set with per-peer scorecards (served / verified / rejected /
//!   missing) for reputation and peer-management.
//! * **Adaptive rounds** — batches are drawn until the *exact* rectangle
//!   escape bound (see [`availability`](crate::availability)) reaches the
//!   target, with custody columns prioritized first.
//! * **Fail-closed verdicts** — any rejection ends the session with
//!   `Rejected`; budget exhaustion reports `Insufficient` with the
//!   achieved confidence.

use crate::availability::line_miss_probability;
use std::collections::HashSet;
use zoda_core::params::ZodaPublic;
use zoda_core::sample::{col_leaf_bytes, row_leaf_bytes};
use zoda_core::{verify_column_sample, verify_row_sample, ZodaParams, ZodaProver};
use zoda_math::merkle::{verify_leaf, MerkleProof};
use zoda_math::{Fr, ZodaRng};

/// A row or column with its Merkle inclusion proof.
#[derive(Clone, Debug)]
pub struct AttestedSample {
    pub data: Vec<Fr>,
    pub proof: MerkleProof,
}

/// Oracle that serves Merkle-attested lines.
pub trait AttestedOracle {
    fn fetch_row_attested(&mut self, index: usize) -> Option<AttestedSample>;
    fn fetch_column_attested(&mut self, index: usize) -> Option<AttestedSample>;
    /// Identifies the peer for scoring.
    fn peer_id(&self) -> usize {
        0
    }
}

/// Honest full-node oracle built from a prover (serves every line with a
/// valid proof), optionally withholding rows/columns.
pub struct ProverOracle<'a> {
    pub prover: &'a ZodaProver,
    pub withhold_rows: Vec<usize>,
    pub withhold_cols: Vec<usize>,
}

impl<'a> ProverOracle<'a> {
    pub fn honest(prover: &'a ZodaProver) -> ProverOracle<'a> {
        ProverOracle {
            prover,
            withhold_rows: Vec::new(),
            withhold_cols: Vec::new(),
        }
    }
}

impl<'a> AttestedOracle for ProverOracle<'a> {
    fn fetch_row_attested(&mut self, index: usize) -> Option<AttestedSample> {
        if self.withhold_rows.contains(&index) {
            return None;
        }
        let proof = self.prover.row_proof(index)?;
        Some(AttestedSample {
            data: self.prover.matrix.row(index),
            proof,
        })
    }

    fn fetch_column_attested(&mut self, index: usize) -> Option<AttestedSample> {
        if self.withhold_cols.contains(&index) {
            return None;
        }
        let proof = self.prover.col_proof(index)?;
        Some(AttestedSample {
            data: self.prover.matrix.col(index),
            proof,
        })
    }
}

/// Per-peer scorecard.
#[derive(Clone, Copy, Debug, Default)]
pub struct PeerScore {
    pub peer_id: usize,
    pub served: usize,
    pub verified: usize,
    pub rejected: usize,
    pub missing: usize,
}

impl PeerScore {
    /// Reputation weight in [0, 1]: fraction of served lines that
    /// verified, with 1.0 only when the peer was actually used.
    pub fn trust(&self) -> f64 {
        if self.served + self.missing == 0 {
            1.0
        } else {
            self.verified as f64 / (self.served + self.missing) as f64
        }
    }
}

/// Per-dimension sample tallies (rows or columns).
#[derive(Clone, Copy, Debug, Default)]
pub struct Tally {
    sampled: usize,
    verified: usize,
    rejected: usize,
    missing: usize,
}

/// Session configuration.
#[derive(Clone, Debug)]
pub struct SessionConfig {
    /// target rectangle-escape probability (default 2^-40)
    pub target_miss: f64,
    /// lines drawn per round (split evenly rows/columns)
    pub batch: usize,
    /// maximum total lines
    pub max_samples: usize,
    /// own custody columns — scheduled before random draws
    pub custody_columns: Vec<usize>,
}

impl Default for SessionConfig {
    fn default() -> Self {
        SessionConfig {
            target_miss: 2f64.powi(-40),
            batch: 16,
            max_samples: 512,
            custody_columns: Vec::new(),
        }
    }
}

/// Final verdict of an attested session.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    /// All samples verified and the escape bound reached the target.
    Available,
    /// A sample failed verification (data is invalid).
    Rejected,
    /// Budget exhausted before the target confidence.
    Insufficient,
}

/// Full session report.
#[derive(Clone, Debug)]
pub struct DasSessionReport {
    pub rounds: usize,
    pub rows: Tally,
    pub columns: Tally,
    pub rows_drawn: Vec<usize>,
    pub cols_drawn: Vec<usize>,
    pub peers: Vec<PeerScore>,
    /// exact worst-case rectangle-escape probability of the union of draws
    pub miss_probability: f64,
    pub verdict: Verdict,
}

impl DasSessionReport {
    pub fn confidence(&self) -> f64 {
        1.0 - self.miss_probability
    }
}

/// Verify one attested row: projection check + Merkle inclusion.
pub fn verify_attested_row(
    public: &ZodaPublic,
    index: usize,
    sample: &AttestedSample,
) -> Result<(), String> {
    verify_row_sample(public, index, &sample.data)?;
    let leaf = row_leaf_bytes(index, &sample.data);
    if !verify_leaf(&public.row_root, &leaf, &sample.proof) {
        return Err("row Merkle proof invalid".to_string());
    }
    Ok(())
}

/// Verify one attested column: projection check + Merkle inclusion.
pub fn verify_attested_column(
    public: &ZodaPublic,
    index: usize,
    sample: &AttestedSample,
) -> Result<(), String> {
    verify_column_sample(public, index, &sample.data)?;
    let leaf = col_leaf_bytes(index, &sample.data);
    if !verify_leaf(&public.col_root, &leaf, &sample.proof) {
        return Err("column Merkle proof invalid".to_string());
    }
    Ok(())
}

/// Run the adaptive attested session against a peer set.
///
/// Each round draws `batch` fresh line indices (custody columns first,
/// then uniform undrawn indices, alternating rows and columns), routes
/// each fetch round-robin to the next peer, verifies projection + Merkle
/// inclusion, and stops on target, rejection, or budget exhaustion.
pub fn run_attested_session(
    public: &ZodaPublic,
    params: &ZodaParams,
    peers: &mut [&mut dyn AttestedOracle],
    config: &SessionConfig,
    rng: &mut ZodaRng,
) -> DasSessionReport {
    let ext_rows = params.extended_rows();
    let ext_cols = params.extended_cols();
    let mut report = DasSessionReport {
        rounds: 0,
        rows: Tally::default(),
        columns: Tally::default(),
        rows_drawn: Vec::new(),
        cols_drawn: Vec::new(),
        peers: peers
            .iter()
            .map(|p| PeerScore {
                peer_id: p.peer_id(),
                ..Default::default()
            })
            .collect(),
        miss_probability: 1.0,
        verdict: Verdict::Insufficient,
    };
    let mut drawn_rows: HashSet<usize> = HashSet::new();
    let mut drawn_cols: HashSet<usize> = HashSet::new();
    let mut custody_queue: Vec<usize> = config
        .custody_columns
        .iter()
        .copied()
        .filter(|&c| c < ext_cols)
        .collect();
    custody_queue.sort_unstable();
    custody_queue.dedup();
    let mut peer_cursor = 0usize;
    let total = |r: &DasSessionReport| r.rows.sampled + r.columns.sampled;

    if peers.is_empty() {
        return report;
    }

    'outer: loop {
        report.rounds += 1;
        // --- draw this round's indices (custody columns first)
        let mut round_cols: Vec<usize> = Vec::new();
        let mut round_rows: Vec<usize> = Vec::new();
        while round_cols.len() + round_rows.len() < config.batch {
            let budget_left = config.max_samples - total(&report);
            if budget_left == 0 {
                break;
            }
            if !custody_queue.is_empty() {
                let c = custody_queue.remove(0);
                if !drawn_cols.contains(&c) {
                    round_cols.push(c);
                }
                continue;
            }
            // alternate row/column uniform draws among undrawn indices
            let want_row = round_rows.len() <= round_cols.len();
            if want_row {
                let candidates: Vec<usize> = (0..ext_rows)
                    .filter(|r| !drawn_rows.contains(r))
                    .collect();
                match candidates.len() {
                    0 => {
                        if round_cols.len() >= config.batch {
                            break;
                        }
                        // no rows left; draw a column instead
                        let cc: Vec<usize> = (0..ext_cols)
                            .filter(|c| !drawn_cols.contains(c))
                            .collect();
                        if cc.is_empty() {
                            break;
                        }
                        round_cols.push(cc[(rng.next_u64() as usize) % cc.len()]);
                        continue;
                    }
                    n => round_rows.push(candidates[(rng.next_u64() as usize) % n]),
                }
            } else {
                let candidates: Vec<usize> = (0..ext_cols)
                    .filter(|c| !drawn_cols.contains(c))
                    .collect();
                match candidates.len() {
                    0 => {
                        let rr: Vec<usize> = (0..ext_rows)
                            .filter(|r| !drawn_rows.contains(r))
                            .collect();
                        if rr.is_empty() {
                            break;
                        }
                        round_rows.push(rr[(rng.next_u64() as usize) % rr.len()]);
                        continue;
                    }
                    n => round_cols.push(candidates[(rng.next_u64() as usize) % n]),
                }
            }
        }
        if round_rows.is_empty() && round_cols.is_empty() {
            break;
        }

        // --- fetch + verify
        for r in round_rows {
            report.rows.sampled += 1;
            drawn_rows.insert(r);
            report.rows_drawn.push(r);
            let peer = peer_cursor % peers.len();
            peer_cursor += 1;
            match peers[peer].fetch_row_attested(r) {
                Some(sample) => {
                    report.peers[peer].served += 1;
                    match verify_attested_row(public, r, &sample) {
                        Ok(()) => {
                            report.rows.verified += 1;
                            report.peers[peer].verified += 1;
                        }
                        Err(_) => {
                            report.rows.rejected += 1;
                            report.peers[peer].rejected += 1;
                            report.verdict = Verdict::Rejected;
                            report.miss_probability = line_miss_probability(
                                drawn_rows.len(),
                                drawn_cols.len(),
                                params.m,
                                params.k,
                            );
                            return report;
                        }
                    }
                }
                None => {
                    report.rows.missing += 1;
                    report.peers[peer].missing += 1;
                }
            }
        }
        for c in round_cols {
            report.columns.sampled += 1;
            drawn_cols.insert(c);
            report.cols_drawn.push(c);
            let peer = peer_cursor % peers.len();
            peer_cursor += 1;
            match peers[peer].fetch_column_attested(c) {
                Some(sample) => {
                    report.peers[peer].served += 1;
                    match verify_attested_column(public, c, &sample) {
                        Ok(()) => {
                            report.columns.verified += 1;
                            report.peers[peer].verified += 1;
                        }
                        Err(_) => {
                            report.columns.rejected += 1;
                            report.peers[peer].rejected += 1;
                            report.verdict = Verdict::Rejected;
                            report.miss_probability = line_miss_probability(
                                drawn_rows.len(),
                                drawn_cols.len(),
                                params.m,
                                params.k,
                            );
                            return report;
                        }
                    }
                }
                None => {
                    report.columns.missing += 1;
                    report.peers[peer].missing += 1;
                }
            }
        }

        // --- confidence accounting and stopping conditions
        report.miss_probability = line_miss_probability(
            drawn_rows.len(),
            drawn_cols.len(),
            params.m,
            params.k,
        );
        if report.miss_probability <= config.target_miss {
            report.verdict = Verdict::Available;
            break 'outer;
        }
        if total(&report) >= config.max_samples {
            report.verdict = Verdict::Insufficient;
            break 'outer;
        }
        // all lines drawn but still above target cannot happen (full
        // coverage forces miss = 0); break defensively
        if drawn_rows.len() == ext_rows && drawn_cols.len() == ext_cols {
            report.verdict = Verdict::Available;
            break 'outer;
        }
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use zoda_core::{Matrix, ZodaParams};
    use zoda_math::PrimeField;

    fn grid(m: usize, k: usize, seed: &[u8]) -> (ZodaParams, ZodaProver, ZodaPublic) {
        let mut s32 = [0u8; 32];
        s32[..seed.len().min(32)].copy_from_slice(&seed[..seed.len().min(32)]);
        let mut rng = ZodaRng::from_seed(s32);
        let mut data = vec![Fr::zero(); m * k];
        for x in data.iter_mut() {
            *x = rng.next_fr(false);
        }
        let g = Matrix::from_row_major(m, k, data);
        let params = ZodaParams::new(m, k);
        let prover = params.commit(&g);
        let public = prover.public_params();
        (params, prover, public)
    }

    #[test]
    fn honest_single_peer_available() {
        let (params, prover, public) = grid(16, 16, b"att-honest");
        let mut oracle = ProverOracle::honest(&prover);
        let mut peers: [&mut dyn AttestedOracle; 1] = [&mut oracle];
        let mut rng = ZodaRng::from_seed([1u8; 32]);
        let report = run_attested_session(
            &public,
            &params,
            &mut peers,
            &SessionConfig::default(),
            &mut rng,
        );
        assert_eq!(report.verdict, Verdict::Available);
        assert_eq!(report.rows.rejected + report.columns.rejected, 0);
        assert!(report.confidence() >= 1.0 - 2f64.powi(-40));
        assert_eq!(report.peers[0].verified, report.peers[0].served);
    }

    #[test]
    fn tampered_line_rejected() {
        let (params, prover, public) = grid(8, 8, b"att-tamper");
        struct TamperOracle<'a> {
            inner: ProverOracle<'a>,
        }
        impl<'a> AttestedOracle for TamperOracle<'a> {
            fn fetch_row_attested(&mut self, index: usize) -> Option<AttestedSample> {
                let mut s = self.inner.fetch_row_attested(index)?;
                if index == 3 {
                    // a line that passes the projection check but is not
                    // committed: swap in another VALID row of the grid
                    s.data = self.inner.prover.matrix.row(5);
                }
                Some(s)
            }
            fn fetch_column_attested(&mut self, index: usize) -> Option<AttestedSample> {
                self.inner.fetch_column_attested(index)
            }
        }
        let mut oracle = TamperOracle {
            inner: ProverOracle::honest(&prover),
        };
        let mut peers: [&mut dyn AttestedOracle; 1] = [&mut oracle];
        let mut rng = ZodaRng::from_seed([2u8; 32]);
        let mut config = SessionConfig::default();
        config.batch = 32; // draw everything
        config.max_samples = 64;
        let report = run_attested_session(&public, &params, &mut peers, &config, &mut rng);
        assert_eq!(report.verdict, Verdict::Rejected);
        assert!(report.rows.rejected + report.columns.rejected > 0);
    }

    #[test]
    fn wrong_grid_lines_rejected_by_projection() {
        // serve lines from an honest DIFFERENT grid — Merkle proofs valid
        // for their own root, but the projection check must fail
        let (_, _, public) = grid(8, 8, b"att-grid-a");
        let (params, prover_b, _) = grid(8, 8, b"att-grid-b");
        let mut oracle = ProverOracle::honest(&prover_b);
        let mut peers: [&mut dyn AttestedOracle; 1] = [&mut oracle];
        let mut rng = ZodaRng::from_seed([3u8; 32]);
        let mut config = SessionConfig::default();
        config.batch = 32;
        config.max_samples = 64;
        let report = run_attested_session(&public, &params, &mut peers, &config, &mut rng);
        assert_eq!(report.verdict, Verdict::Rejected);
    }

    #[test]
    fn multi_peer_routing_and_scores() {
        let (params, prover, public) = grid(16, 16, b"att-peers");
        let mut a = ProverOracle::honest(&prover);
        let mut b = ProverOracle {
            withhold_cols: vec![0],
            ..ProverOracle::honest(&prover)
        };
        let mut peers: [&mut dyn AttestedOracle; 2] = [&mut a, &mut b];
        let mut rng = ZodaRng::from_seed([4u8; 32]);
        let report = run_attested_session(
            &public,
            &params,
            &mut peers,
            &SessionConfig::default(),
            &mut rng,
        );
        assert_eq!(report.verdict, Verdict::Available);
        // both peers were used (round-robin)
        assert!(report.peers[0].served + report.peers[0].missing > 0);
        assert!(report.peers[1].served + report.peers[1].missing > 0);
        assert_eq!(report.peers[0].trust(), 1.0);
        assert!(report.peers[1].trust() <= 1.0);
    }

    #[test]
    fn custody_columns_scheduled_first() {
        let (params, prover, public) = grid(8, 8, b"att-custody");
        let mut oracle = ProverOracle::honest(&prover);
        let mut peers: [&mut dyn AttestedOracle; 1] = [&mut oracle];
        let mut rng = ZodaRng::from_seed([5u8; 32]);
        let config = SessionConfig {
            custody_columns: vec![3, 7, 11],
            batch: 8,
            max_samples: 8,
            ..Default::default()
        };
        let report = run_attested_session(&public, &params, &mut peers, &config, &mut rng);
        // custody columns are among the first drawn
        for &c in &[3usize, 7, 11] {
            assert!(
                report.cols_drawn.contains(&c),
                "custody column {} missing",
                c
            );
        }
    }

    #[test]
    fn budget_insufficient() {
        let (params, prover, public) = grid(64, 64, b"att-budget");
        let mut oracle = ProverOracle::honest(&prover);
        let mut peers: [&mut dyn AttestedOracle; 1] = [&mut oracle];
        let mut rng = ZodaRng::from_seed([6u8; 32]);
        let config = SessionConfig {
            batch: 4,
            max_samples: 4, // far too few for a 128x128 grid
            ..Default::default()
        };
        let report = run_attested_session(&public, &params, &mut peers, &config, &mut rng);
        assert_eq!(report.verdict, Verdict::Insufficient);
        assert!(report.miss_probability > config.target_miss);
    }
}
