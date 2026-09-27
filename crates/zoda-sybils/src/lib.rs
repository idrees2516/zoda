//! # zoda-sybils — sybil resistance for DA participants
//!
//! * **Sortition**: BLS signatures over a beacon seed select sampling
//!   committees and custody assignments (a signature over a secret seed
//!   is unpredictable and verifiable by anyone),
//! * **Stake-weighted sampling**: peers are selected with probability
//!   proportional to stake (alias-reservoir method),
//! * **Peer scoring**: windowed reputation with slashing-style decay for
//!   peers that serve invalid samples.

use zoda_bls::{PublicKey, SecretKey, Signature};
use zoda_math::sha256::sha256;
use zoda_math::{Fr, PrimeField, ZodaRng};

/// A registered validator/peer.
#[derive(Clone, Debug)]
pub struct Peer {
    pub id: [u8; 32],
    pub stake: u64,
    pub public_key: PublicKey,
}

/// Deterministic sortition: which peers are selected for a duty, derived
/// from a BLS signature over the (slot, duty) seed.
pub struct Sortition;

impl Sortition {
    /// The beacon seed: sign `domain || slot || duty` with the committee
    /// key. Any party can verify the signature and derive the same seed.
    pub fn derive_seed(
        domain: &[u8],
        slot: u64,
        duty: &[u8],
        sk: &SecretKey,
    ) -> ([u8; 32], Signature) {
        let mut msg = Vec::with_capacity(domain.len() + 8 + duty.len());
        msg.extend_from_slice(domain);
        msg.extend_from_slice(&slot.to_le_bytes());
        msg.extend_from_slice(duty);
        let sig = sk.sign(&msg);
        (sha256(&msg), sig)
    }

    /// Select `count` peers weighted by stake using the seed (VRF-style
    /// lottery without replacement).
    pub fn select_weighted(seed: &[u8; 32], peers: &[Peer], count: usize) -> Vec<usize> {
        // alias-reservoir weighted sampling without replacement
        let mut rng = ZodaRng::from_seed(*seed);
        let total: u128 = peers.iter().map(|p| p.stake as u128 + 1).sum();
        let mut chosen: Vec<usize> = Vec::with_capacity(count);
        let mut remaining = total;
        let mut alive: Vec<usize> = (0..peers.len()).collect();
        while chosen.len() < count && !alive.is_empty() {
            // draw a uniform point in [0, remaining)
            let x = rng.next_below(remaining.min(u64::MAX as u128) as u64) as u128;
            let mut acc: u128 = 0;
            let mut picked = alive.len() - 1;
            for (pos, &i) in alive.iter().enumerate() {
                acc += peers[i].stake as u128 + 1;
                if x < acc {
                    picked = pos;
                    break;
                }
            }
            let peer = alive.swap_remove(picked);
            remaining -= peers[peer].stake as u128 + 1;
            chosen.push(peer);
        }
        chosen
    }

    /// Uniform (unweighted) selection — used for sampling committees
    /// before stake accounting is available.
    pub fn select_uniform(seed: &[u8; 32], n: usize, count: usize) -> Vec<usize> {
        let mut rng = ZodaRng::from_seed(*seed);
        rng.shuffle_indices(n)[..count.min(n)].to_vec()
    }
}

/// Windowed peer scoring (EigenTrust-style exponential decay).
#[derive(Clone, Debug)]
pub struct PeerScorer {
    scores: Vec<f64>,
    decay: f64,
}

impl PeerScorer {
    pub fn new(peers: usize, decay: f64) -> PeerScorer {
        PeerScorer {
            scores: vec![1.0; peers],
            decay,
        }
    }

    pub fn score(&self, peer: usize) -> f64 {
        self.scores[peer]
    }

    /// Record a good interaction (a sample that verified).
    pub fn reward(&mut self, peer: usize, amount: f64) {
        self.scores[peer] += amount;
    }

    /// Record a bad interaction (an invalid or missing sample).
    pub fn penalize(&mut self, peer: usize, amount: f64) {
        self.scores[peer] = (self.scores[peer] - amount).max(0.0);
    }

    /// Age all scores one epoch.
    pub fn decay_epoch(&mut self) {
        for s in self.scores.iter_mut() {
            *s = 1.0 + (*s - 1.0) * self.decay;
        }
    }

    /// Peers sorted by score (best first).
    pub fn ranking(&self) -> Vec<usize> {
        let mut idx: Vec<usize> = (0..self.scores.len()).collect();
        idx.sort_by(|a, b| self.scores[*b].total_cmp(&self.scores[*a]));
        idx
    }
}

/// A BLS-signed sampling assignment: proves the assigner is who they say
/// and binds the assigned indices.
#[derive(Clone)]
pub struct SamplingAssignment {
    pub slot: u64,
    pub indices: Vec<usize>,
    pub signature: Signature,
}

pub fn sign_assignment(sk: &SecretKey, slot: u64, indices: &[usize]) -> SamplingAssignment {
    let mut msg = Vec::new();
    msg.extend_from_slice(b"zoda-sampling-assignment-v1");
    msg.extend_from_slice(&slot.to_le_bytes());
    for &i in indices {
        msg.extend_from_slice(&(i as u64).to_le_bytes());
    }
    SamplingAssignment {
        slot,
        indices: indices.to_vec(),
        signature: sk.sign(&msg),
    }
}

pub fn verify_assignment(pk: &PublicKey, a: &SamplingAssignment) -> bool {
    let mut msg = Vec::new();
    msg.extend_from_slice(b"zoda-sampling-assignment-v1");
    msg.extend_from_slice(&a.slot.to_le_bytes());
    for &i in &a.indices {
        msg.extend_from_slice(&(i as u64).to_le_bytes());
    }
    zoda_bls::verify(pk, &msg, &a.signature)
}

/// Map a hash into Fr (for challenge derivations).
pub fn hash_to_fr(data: &[u8]) -> Fr {
    Fr::from_be_bytes_mod_order(&sha256(data))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peers(n: usize) -> Vec<Peer> {
        (0..n)
            .map(|i| {
                let sk = SecretKey::from_seed(format!("peer-key-{}", i).as_bytes());
                Peer {
                    id: sha256(format!("peer-{}", i).as_bytes()),
                    stake: 10 * (i as u64 + 1),
                    public_key: sk.public_key(),
                }
            })
            .collect()
    }

    #[test]
    fn sortition_is_deterministic() {
        let ps = peers(8);
        let a = Sortition::select_weighted(&[1u8; 32], &ps, 3);
        let b = Sortition::select_weighted(&[1u8; 32], &ps, 3);
        assert_eq!(a, b);
        let c = Sortition::select_weighted(&[2u8; 32], &ps, 3);
        assert_ne!(a, c);
    }

    #[test]
    fn stake_weighting_counts() {
        // a peer with the majority of stake should dominate selections
        let mut ps = peers(3);
        ps[0].stake = 1_000_000;
        ps[1].stake = 1;
        ps[2].stake = 1;
        let mut wins = [0usize; 3];
        for s in 0u64..40 {
            let seed = sha256(&s.to_le_bytes());
            for i in Sortition::select_weighted(&seed, &ps, 1) {
                wins[i] += 1;
            }
        }
        assert!(wins[0] > wins[1] + wins[2], "stake weighting failed: {:?}", wins);
    }

    #[test]
    fn assignment_signatures() {
        let sk = SecretKey::from_seed(b"assignment-key-000000000000000000000000");
        let pk = sk.public_key();
        let a = sign_assignment(&sk, 123, &[1, 5, 9]);
        assert!(verify_assignment(&pk, &a));
        // wrong slot fails
        let mut bad = a.clone();
        bad.slot = 124;
        assert!(!verify_assignment(&pk, &bad));
        // wrong indices fail
        let mut bad2 = a;
        bad2.indices[0] = 2;
        assert!(!verify_assignment(&pk, &bad2));
    }

    #[test]
    fn scorer_ranks_honesty() {
        let mut s = PeerScorer::new(3, 0.9);
        s.penalize(1, 0.9);
        s.reward(2, 5.0);
        let ranking = s.ranking();
        assert_eq!(ranking[0], 2);
        assert_eq!(ranking[2], 1);
        s.decay_epoch();
        // ordering survives decay
        assert_eq!(s.ranking()[0], 2);
    }
}
