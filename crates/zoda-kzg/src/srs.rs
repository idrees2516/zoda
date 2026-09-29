//! Trusted setup (structured reference string) handling.
//!
//! Supports the c-kzg-4844 `trusted_setup.txt` formats:
//! * 3-section (current): `n1 n2` then G1-Lagrange, G2-monomial,
//!   G1-monomial (the EIP-7594 extension).
//! * 2-section (legacy): `n1 n2` then G1-Lagrange, G2-monomial — the
//!   monomial G1 points are then derived via an inverse NTT.
//!
//! The spec's sanity check (`e(G1[1], G2[0]) != e(G1[0], G2[1])` ⇒ the G1
//! points really are in Lagrange form) is always enforced.

use crate::msm;
use std::path::Path;
use zoda_bls::g1::G1Affine;
use zoda_bls::g2::G2Affine;
use zoda_bls::pairing::G2Prepared;
use zoda_math::u256::hex_to_bytes;
use zoda_math::{FftDomain, Fr, PrimeField};

pub const FIELD_ELEMENTS_PER_BLOB: usize = 4096;
pub const BYTES_PER_FIELD_ELEMENT: usize = 32;
pub const BLOB_BYTES: usize = FIELD_ELEMENTS_PER_BLOB * BYTES_PER_FIELD_ELEMENT;
pub const KZG_SETUP_G2_LENGTH: usize = 65;
/// The size of the extended (EIP-7594) evaluation domain: 8192.
pub const FIELD_ELEMENTS_PER_EXT_BLOB: usize = 2 * FIELD_ELEMENTS_PER_BLOB;

/// The precomputed KZG context.
pub struct Setup {
    /// Roots of unity of the 8192-domain, natural order (`roots[i] = ω^i`).
    pub roots_of_unity: Vec<Fr>,
    /// Bit-reversed roots of unity (4096-domain order for blob evals).
    pub brp_roots_of_unity: Vec<Fr>,
    /// G1 trusted setup in Lagrange form, bit-reversed.
    pub g1_lagrange_brp: Vec<G1Affine>,
    /// G1 trusted setup in monomial form (derived or parsed).
    pub g1_monomial: Vec<G1Affine>,
    /// G2 trusted setup in monomial form (65 points).
    pub g2_monomial: Vec<G2Affine>,
    /// Cached G2 generator (affine).
    pub g2_gen: G2Affine,
    /// FFT domain of the extended blob size (8192).
    pub ext_domain: FftDomain<Fr>,
    /// FFT domain of the blob size (4096).
    pub blob_domain: FftDomain<Fr>,
}

impl Setup {
    /// Load from a `trusted_setup.txt` file (2- or 3-section format).
    pub fn load_file(path: &Path) -> Result<Setup, String> {
        let data = std::fs::read_to_string(path)
            .map_err(|e| format!("failed to read trusted setup: {}", e))?;
        Setup::parse(&data)
    }

    /// Parse trusted setup text.
    pub fn parse(text: &str) -> Result<Setup, String> {
        let mut tok = text.split_whitespace();
        let n1: usize = tok
            .next()
            .and_then(|t| t.parse().ok())
            .ok_or("missing G1 count")?;
        let n2: usize = tok
            .next()
            .and_then(|t| t.parse().ok())
            .ok_or("missing G2 count")?;
        if n1 != FIELD_ELEMENTS_PER_BLOB || n2 != KZG_SETUP_G2_LENGTH {
            return Err(format!(
                "unsupported trusted setup dimensions {} x {}",
                n1, n2
            ));
        }
        let mut g1_lagrange = Vec::with_capacity(n1);
        for _ in 0..n1 {
            let hex = tok.next().ok_or("truncated G1 section")?;
            let bytes = hex_to_bytes(hex).ok_or("bad G1 hex")?;
            let p = G1Affine::from_compressed(&bytes)
                .ok_or("invalid G1 point in trusted setup")?;
            g1_lagrange.push(p);
        }
        let mut g2_monomial = Vec::with_capacity(n2);
        for _ in 0..n2 {
            let hex = tok.next().ok_or("truncated G2 section")?;
            let bytes = hex_to_bytes(hex).ok_or("bad G2 hex")?;
            let p = G2Affine::from_compressed(&bytes)
                .ok_or("invalid G2 point in trusted setup")?;
            g2_monomial.push(p);
        }
        // optional third section: G1 monomial
        let g1_monomial = match tok.next() {
            Some(hex) => {
                let mut g1m = Vec::with_capacity(n1 + 1);
                let bytes = hex_to_bytes(hex).ok_or("bad G1 hex")?;
                g1m.push(
                    G1Affine::from_compressed(&bytes).ok_or("invalid G1 monomial point")?,
                );
                for _ in 1..n1 {
                    let hex = tok.next().ok_or("truncated G1 monomial section")?;
                    let bytes = hex_to_bytes(hex).ok_or("bad G1 hex")?;
                    g1m.push(
                        G1Affine::from_compressed(&bytes)
                            .ok_or("invalid G1 monomial point")?,
                    );
                }
                g1m
            }
            None => {
                // derive the monomial form from the Lagrange form
                derive_g1_monomial(&g1_lagrange)?
            }
        };
        Setup::build(g1_lagrange, g1_monomial, g2_monomial)
    }

    fn build(
        g1_lagrange: Vec<G1Affine>,
        g1_monomial: Vec<G1Affine>,
        g2_monomial: Vec<G2Affine>,
    ) -> Result<Setup, String> {
        if g1_lagrange.len() != FIELD_ELEMENTS_PER_BLOB
            || g1_monomial.len() != FIELD_ELEMENTS_PER_BLOB
            || g2_monomial.len() != KZG_SETUP_G2_LENGTH
        {
            return Err("setup length mismatch".to_string());
        }
        // Spec sanity check: if e(G1[1], G2[0]) == e(G1[0], G2[1]) the G1
        // points are in monomial form — reject them.
        {
            let g2_0_prepared = G2Prepared::from(g2_monomial[0]);
            let g2_1_prepared = G2Prepared::from(g2_monomial[1]);
            let is_monomial = zoda_bls::pairing::pairing_check(&[
                (&g1_lagrange[1], &g2_0_prepared),
                (&(-g1_lagrange[0]).to_projective().to_affine(), &g2_1_prepared),
            ]);
            if is_monomial {
                return Err(
                    "trusted setup G1 appears to be in monomial form; \
                     Lagrange form is required"
                        .to_string(),
                );
            }
        }

        let ext_domain = FftDomain::<Fr>::new(FIELD_ELEMENTS_PER_EXT_BLOB);
        let blob_domain = FftDomain::<Fr>::new(FIELD_ELEMENTS_PER_BLOB);

        // natural-order roots of unity of the 8192 domain (plus wrap entry)
        let roots_of_unity = ext_domain.roots_of_unity().to_vec();
        // BRP'd 4096-domain roots (for blob evaluation form)
        let mut brp = blob_domain.roots_of_unity()[..FIELD_ELEMENTS_PER_BLOB].to_vec();
        bit_reverse(&mut brp);

        // BRP the Lagrange points
        let mut g1_lagrange_brp = g1_lagrange;
        bit_reverse_points(&mut g1_lagrange_brp);

        let g2_gen = g2_monomial[0];
        Ok(Setup {
            roots_of_unity,
            brp_roots_of_unity: brp,
            g1_lagrange_brp,
            g1_monomial,
            g2_monomial,
            g2_gen,
            ext_domain,
            blob_domain,
        })
    }

    /// A deterministic toy setup for tests and benchmarks: powers of a
    /// random-but-fixed τ (WITH the "toxic waste" known — obviously NOT for
    /// production use; the real ceremony setup must be loaded from file).
    pub fn from_seed_for_testing(seed: [u8; 32]) -> Setup {
        let mut rng = zoda_math::ZodaRng::from_seed(seed);
        let tau = rng.next_fr(true);
        // monomial G1: [tau^i]
        let mut g1_monomial = Vec::with_capacity(FIELD_ELEMENTS_PER_BLOB);
        let mut p = Fr::ONE;
        for _ in 0..FIELD_ELEMENTS_PER_BLOB {
            g1_monomial.push(
                zoda_bls::g1::G1Projective::generator()
                    .mul_fr(&p)
                    .to_affine(),
            );
            p = p * tau;
        }
        // G2: [tau^0, tau^1, ..., tau^64]
        let mut g2_monomial = Vec::with_capacity(KZG_SETUP_G2_LENGTH);
        let mut p = Fr::ONE;
        for _ in 0..KZG_SETUP_G2_LENGTH {
            g2_monomial.push(
                zoda_bls::g2::G2Projective::generator()
                    .mul_fr(&p)
                    .to_affine(),
            );
            p = p * tau;
        }
        // Lagrange form: inverse NTT of the monomial points over the
        // 4096-domain.
        let g1_lagrange = lagrange_from_monomial(&g1_monomial);
        Setup::build(g1_lagrange, g1_monomial, g2_monomial)
            .expect("test setup construction failed")
    }

    /// Commit to an evaluation-form polynomial over the blob domain via MSM.
    pub fn commit_evaluations(&self, evals: &[Fr]) -> G1Affine {
        debug_assert_eq!(evals.len(), FIELD_ELEMENTS_PER_BLOB);
        msm::g1_lincomb(&self.g1_lagrange_brp, evals)
    }

    /// Commit to a coefficient-form polynomial (MSM over monomial points).
    pub fn commit_coeffs(&self, coeffs: &[Fr]) -> G1Affine {
        debug_assert!(coeffs.len() <= self.g1_monomial.len());
        msm::g1_lincomb(&self.g1_monomial[..coeffs.len()], coeffs)
    }
}

/// Decode a "0x…" hex string into bytes (public helper for the spec
/// vector harness).
pub fn hex_bytes(s: &str) -> Option<Vec<u8>> {
    zoda_math::u256::hex_to_bytes(s)
}

/// Bit-reverse a slice.
pub fn bit_reverse<T>(v: &mut [T]) {
    let n = v.len();
    if !n.is_power_of_two() {
        return;
    }
    for i in 0..n {
        let j = reverse_bits(n, i);
        if j > i {
            v.swap(i, j);
        }
    }
}

fn reverse_bits(width: usize, mut x: usize) -> usize {
    let mut r = 0;
    let bits = width.trailing_zeros();
    for _ in 0..bits {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    r
}

fn bit_reverse_points(v: &mut [G1Affine]) {
    bit_reverse(v);
}

/// Derive the monomial-form G1 points from Lagrange form via a "G1 IFFT":
/// the i-th monomial point = Σⱼ L_j(τ^i)·G1_lag[j] computed through the
/// FFT relationship (equivalent to c-kzg's lagrange→monomial conversion).
fn derive_g1_monomial(g1_lagrange: &[G1Affine]) -> Result<Vec<G1Affine>, String> {
    let n = FIELD_ELEMENTS_PER_BLOB;
    if g1_lagrange.len() != n {
        return Err("bad lagrange length".into());
    }
    // We need the BRP'd Lagrange form and then an inverse FFT over G1:
    // monomial[i] = (1/n) Σ_k lag_brp_fft... The direct construction:
    // treat the Lagrange points as evaluations over the BRP domain and
    // run the inverse FFT over the group (projective).
    let mut brp = g1_lagrange.to_vec();
    bit_reverse_points(&mut brp);
    // inverse FFT over G1 (recursive, matching the field-domain convention)
    let mut coeffs = g1_ifft(&brp);
    bit_reverse_points(&mut coeffs);
    let _ = n;
    Ok(coeffs)
}

/// Inverse (unscaled) FFT over G1 points using the 4096-domain roots.
fn g1_ifft(input: &[G1Affine]) -> Vec<G1Affine> {
    let n = input.len();
    debug_assert!(n.is_power_of_two());
    let domain = FftDomain::<Fr>::new(n);
    // inverse = forward with inverse roots, then scale by 1/n
    let mut out = g1_fft_with(input, &inverse_roots(&domain));
    let inv_n = Fr::from_u64(n as u64).invert().unwrap();
    for p in out.iter_mut() {
        *p = p.to_projective().mul_fr(&inv_n).to_affine();
    }
    out
}

fn inverse_roots(domain: &FftDomain<Fr>) -> Vec<Fr> {
    domain.roots_of_unity()[..domain.size()]
        .iter()
        .map(|r| r.invert().unwrap())
        .collect()
}

/// Forward FFT over G1 with given twiddle order — thin wrapper around
/// the Jacobian engine in `crate::g1fft` (unit-twiddle and identity
/// skips, single batch normalization).
pub fn g1_fft_with(input: &[G1Affine], roots: &[Fr]) -> Vec<G1Affine> {
    let jacs: Vec<crate::msm::G1J> =
        input.iter().map(|p| crate::msm::G1J::from_affine(p)).collect();
    let out = crate::g1fft::g1_fft_jac(&jacs, roots);
    crate::msm::G1J::batch_normalize(&out)
}

/// Lagrange-form G1 points from monomial form: L_i(τ) = IFFT of τ-powers;
/// equivalently g1_lagrange = FFT^{-1}(monomial) — computed as a forward
/// FFT with inverse twiddles scaled by 1/n.
pub fn lagrange_from_monomial(g1_monomial: &[G1Affine]) -> Vec<G1Affine> {
    let n = g1_monomial.len();
    let domain = FftDomain::<Fr>::new(n);
    let mut out = g1_fft_with(g1_monomial, &inverse_roots(&domain));
    // scale by 1/n
    let inv_n = Fr::from_u64(n as u64).invert().unwrap();
    for p in out.iter_mut() {
        *p = p.to_projective().mul_fr(&inv_n).to_affine();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_setup_consistency() {
        let s = Setup::from_seed_for_testing(*b"srs-test-seed-000000000000000000");
        // The sanity check must accept our (Lagrange-form) setup.
        // Committing to the polynomial "1" (first Lagrange basis) must give
        // the first Lagrange point.
        let mut evals = vec![Fr::zero(); FIELD_ELEMENTS_PER_BLOB];
        evals[0] = Fr::ONE;
        // NOTE: g1_lagrange_brp is BRP'd; evals must be interpreted in BRP
        // order too — the identity basis element is index 0 either way.
        let c = s.commit_evaluations(&evals);
        assert_eq!(c, s.g1_lagrange_brp[0]);
    }

    #[test]
    fn monomial_lagrange_roundtrip() {
        let s = Setup::from_seed_for_testing(*b"srs-test-seed-000000000000000000");
        // lagrange_from_monomial(g1_monomial) == g1_lagrange (pre-BRP)
        // We can't access the pre-BRP lagrange directly, so verify via a
        // commitment: commit_evaluations(evals) == commit_coeffs(coeffs)
        // where coeffs = IFFT(evals).
        let mut rng = zoda_math::ZodaRng::from_seed(*b"srs-roundtrip-000000000000000000");
        let mut evals = vec![Fr::zero(); FIELD_ELEMENTS_PER_BLOB];
        for e in evals.iter_mut() {
            *e = rng.next_fr(false);
        }
        let mut brp_evals = evals.clone();
        bit_reverse(&mut brp_evals);
        let mut coeffs = vec![Fr::zero(); FIELD_ELEMENTS_PER_BLOB];
        s.blob_domain.ifft(&brp_evals, &mut coeffs);
        let c_eval = s.commit_evaluations(&evals);
        let c_coeff = s.commit_coeffs(&coeffs);
        assert_eq!(c_eval, c_coeff);
    }

    #[test]
    fn parse_mainnet_setup_file() {
        // only run if the reference trusted setup is available
        let path = std::path::Path::new("/home/z/my-project/refs/spec/trusted_setup.txt");
        if !path.exists() {
            return;
        }
        let s = Setup::load_file(path).expect("parse mainnet setup");
        assert_eq!(s.g1_lagrange_brp.len(), 4096);
        assert_eq!(s.g2_monomial.len(), 65);
        assert_eq!(s.g1_monomial.len(), 4096);
        // G2[0] must be the standard generator
        let g2_hex: String = {
            let b = s.g2_monomial[0].to_compressed();
            b.iter().map(|x| format!("{:02x}", x)).collect()
        };
        assert!(g2_hex.starts_with("93e02b6052719f60"), "G2[0] != generator");
    }
}
