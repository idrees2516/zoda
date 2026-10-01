//! The EIP-4844 KZG API, exactly following the Deneb polynomial-commitments
//! specification (evaluation-form KZG with bit-reversed roots of unity).

use crate::msm;
use crate::srs::{Setup, BLOB_BYTES, FIELD_ELEMENTS_PER_BLOB};
use zoda_bls::g1::G1Projective;
use zoda_bls::pairing::{pairing_check, G2Prepared};
use zoda_math::sha256::sha256;
use zoda_math::{batch_invert, Fr, PrimeField};

pub const BYTES_PER_COMMITMENT: usize = 48;
pub const BYTES_PER_PROOF: usize = 48;
pub const BYTES_PER_FIELD_ELEMENT: usize = 32;
pub const KZG_ENDIANNESS_BIG: bool = true;
/// PRIMITIVE_ROOT_OF_UNITY = 7.
pub const PRIMITIVE_ROOT_OF_UNITY: u64 = 7;

pub const FIAT_SHAMIR_PROTOCOL_DOMAIN: &[u8] = b"FSBLOBVERIFY_V1_";
pub const RANDOM_CHALLENGE_KZG_BATCH_DOMAIN: &[u8] = b"RCKZGBATCH___V1_";
pub const VERSIONED_HASH_VERSION_KZG: u8 = 0x01;

/// A KZG commitment (compressed G1, 48 bytes).
pub type KzgCommitment = [u8; 48];
/// A KZG proof (compressed G1, 48 bytes).
pub type KzgProof = [u8; 48];
/// A blob: 4096 field elements, 32 bytes each, big-endian.
pub type Blob = [u8; BLOB_BYTES];

/// `hash_to_bls_field`: SHA-256 then reduce mod r.
pub fn hash_to_bls_field(data: &[u8]) -> Fr {
    Fr::from_be_bytes_mod_order(&sha256(data))
}

/// `bytes_to_bls_field` with validation (< r, 32 bytes).
///
/// Fast path: parse the four 64-bit limbs directly (big-endian), check
/// canonicality with three limb comparisons, and enter Montgomery form
/// with a single multiplication — an order of magnitude cheaper than
/// reduce-then-re-encode, which matters because every blob operation
/// starts by parsing 4096 field elements.
pub fn bytes_to_bls_field(b: &[u8]) -> Result<Fr, String> {
    if b.len() != BYTES_PER_FIELD_ELEMENT {
        return Err("field element must be 32 bytes".to_string());
    }
    let mut limbs = [0u64; 4]; // little-endian canonical limbs
    for i in 0..4 {
        let mut w = [0u8; 8];
        w.copy_from_slice(&b[24 - i * 8..32 - i * 8]);
        limbs[i] = u64::from_be_bytes(w);
    }
    // canonicality: value < r (top-down limb comparison; equality —
    // value == r — is also rejected)
    let m = zoda_math::fr::FR_MODULUS;
    let mut i = 4;
    let mut all_equal = true;
    while i > 0 {
        i -= 1;
        if limbs[i] != m[i] {
            all_equal = false;
            if limbs[i] > m[i] {
                return Err("field element is not canonical (>= r)".to_string());
            }
            break;
        }
    }
    if all_equal {
        return Err("field element is not canonical (>= r)".to_string());
    }
    Ok(Fr::from_repr_limbs(limbs))
}

/// `bls_field_to_bytes`: big-endian 32 bytes.
pub fn fr_to_bytes(x: Fr) -> [u8; 32] {
    let le = x.to_le_bytes();
    let mut be = [0u8; 32];
    for (i, byte) in le.iter().rev().enumerate() {
        be[i] = *byte;
    }
    be
}

/// `blob_to_polynomial`: interpret the blob as 4096 field elements.
pub fn blob_to_polynomial(blob: &[u8]) -> Result<Vec<Fr>, String> {
    if blob.len() != BLOB_BYTES {
        return Err("blob must be 131072 bytes".to_string());
    }
    let mut poly = Vec::with_capacity(FIELD_ELEMENTS_PER_BLOB);
    for i in 0..FIELD_ELEMENTS_PER_BLOB {
        poly.push(bytes_to_bls_field(
            &blob[i * 32..(i + 1) * 32],
        )?);
    }
    Ok(poly)
}

/// `compute_challenge` (Fiat-Shamir for blob proofs).
pub fn compute_challenge(blob: &[u8], commitment: &[u8; 48]) -> Fr {
    let mut data = Vec::with_capacity(16 + 16 + blob.len() + 48);
    data.extend_from_slice(FIAT_SHAMIR_PROTOCOL_DOMAIN);
    // the degree separator is 16 bytes big-endian (the batch-verifier
    // transcript uses 8 bytes — an intentional spec asymmetry)
    data.extend_from_slice(&(FIELD_ELEMENTS_PER_BLOB as u128).to_be_bytes());
    data.extend_from_slice(blob);
    data.extend_from_slice(commitment);
    hash_to_bls_field(&data)
}

/// `evaluate_polynomial_in_evaluation_form`: barycentric evaluation of the
/// blob polynomial (evaluations in bit-reversed order) at z.
pub fn evaluate_polynomial_in_evaluation_form(
    polynomial: &[Fr],
    z: Fr,
    setup: &Setup,
) -> Fr {
    // BRP roots of unity
    let roots = &setup.brp_roots_of_unity;
    // in-domain fast path
    for (i, &w) in roots.iter().enumerate() {
        if w == z {
            return polynomial[i];
        }
    }
    let n = FIELD_ELEMENTS_PER_BLOB;
    let mut denominators = Vec::with_capacity(n);
    for &w in roots.iter() {
        denominators.push(z - w);
    }
    batch_invert(&mut denominators);
    let mut acc = Fr::zero();
    for i in 0..n {
        acc = acc + polynomial[i] * roots[i] * denominators[i];
    }
    // r = z^n - 1
    let mut zn = Fr::ONE;
    let mut e = n as u64;
    let mut base = z;
    while e > 0 {
        if e & 1 == 1 {
            zn = zn * base;
        }
        base = base * base;
        e >>= 1;
    }
    let inv_n = Fr::from_u64(n as u64).invert().unwrap();
    acc * (zn - Fr::ONE) * inv_n
}

/// `compute_quotient_eval_within_domain`.
pub fn compute_quotient_eval_within_domain(
    z: Fr,
    polynomial: &[Fr],
    y: Fr,
    setup: &Setup,
) -> Fr {
    let roots = &setup.brp_roots_of_unity;
    let mut result = Fr::zero();
    let mut denominators = Vec::with_capacity(FIELD_ELEMENTS_PER_BLOB);
    for &w in roots.iter() {
        if w == z {
            continue;
        }
        denominators.push(z * (z - w));
    }
    batch_invert(&mut denominators);
    let mut di = 0;
    for (i, &w) in roots.iter().enumerate() {
        if w == z {
            continue;
        }
        let f_i = polynomial[i] - y;
        result = result + f_i * w * denominators[di];
        di += 1;
    }
    result
}

/// `blob_to_kzg_commitment` (public method).
pub fn blob_to_kzg_commitment(blob: &[u8], setup: &Setup) -> Result<[u8; 48], String> {
    if blob.len() != BLOB_BYTES {
        return Err("blob must be 131072 bytes".to_string());
    }
    let poly = blob_to_polynomial(blob)?;
    // spec: g1_lincomb(bit_reversal_permutation(KZG_SETUP_G1_LAGRANGE), poly)
    // our g1_lagrange_brp is already BRP'd; the blob evaluations are in
    // natural order, so the MSM is against the BRP'd points directly.
    let c = msm::g1_lincomb(&setup.g1_lagrange_brp, &poly);
    Ok(to_bytes48(c))
}

/// `compute_kzg_proof` (public method): proof that p(z) = y.
pub fn compute_kzg_proof(
    blob: &[u8],
    z_bytes: &[u8],
    setup: &Setup,
) -> Result<([u8; 48], [u8; 32]), String> {
    if blob.len() != BLOB_BYTES {
        return Err("blob must be 131072 bytes".to_string());
    }
    if z_bytes.len() != BYTES_PER_FIELD_ELEMENT {
        return Err("z must be 32 bytes".to_string());
    }
    let polynomial = blob_to_polynomial(blob)?;
    let z = bytes_to_bls_field(z_bytes)?;
    let (proof, y) = compute_kzg_proof_impl(&polynomial, z, setup);
    Ok((to_bytes48(proof), fr_to_bytes(y)))
}

/// `compute_kzg_proof_impl`.
pub fn compute_kzg_proof_impl(
    polynomial: &[Fr],
    z: Fr,
    setup: &Setup,
) -> (zoda_bls::g1::G1Affine, Fr) {
    let roots = &setup.brp_roots_of_unity;
    let y = evaluate_polynomial_in_evaluation_form(polynomial, z, setup);

    // quotient in evaluation form: q(x_i) = (p(x_i) - p(z)) / (x_i - z)
    let mut quotient = vec![Fr::zero(); FIELD_ELEMENTS_PER_BLOB];
    let mut denominators = Vec::with_capacity(FIELD_ELEMENTS_PER_BLOB);
    let mut special = None;
    for (i, &w) in roots.iter().enumerate() {
        let d = w - z;
        if d.is_zero() {
            special = Some(i);
            denominators.push(Fr::ONE); // placeholder
        } else {
            denominators.push(d);
        }
    }
    batch_invert(&mut denominators);
    for i in 0..FIELD_ELEMENTS_PER_BLOB {
        let a = polynomial[i] - y;
        quotient[i] = a * denominators[i];
    }
    if let Some(i) = special {
        quotient[i] = compute_quotient_eval_within_domain(z, polynomial, y, setup);
    }
    let proof = msm::g1_lincomb(&setup.g1_lagrange_brp, &quotient);
    (proof, y)
}

/// `verify_kzg_proof` (public method).
pub fn verify_kzg_proof(
    commitment_bytes: &[u8],
    z_bytes: &[u8],
    y_bytes: &[u8],
    proof_bytes: &[u8],
    setup: &Setup,
) -> Result<bool, String> {
    let commitment = bytes_to_kzg_commitment(commitment_bytes)?;
    let z = bytes_to_bls_field(z_bytes)?;
    let y = bytes_to_bls_field(y_bytes)?;
    let proof = bytes_to_kzg_proof(proof_bytes)?;
    Ok(verify_kzg_proof_impl(&commitment, z, y, &proof, setup))
}

/// `validate_kzg_g1` + parse: a canonical compressed encoding of a
/// curve point **in the order-r subgroup** ([r]P = O; the pairing cannot
/// detect h-torsion components, so the explicit check is required for
/// soundness, matching c-kzg's `validate_kzg_g1`).
pub fn bytes_to_kzg_commitment(b: &[u8]) -> Result<zoda_bls::g1::G1Affine, String> {
    if b.len() != BYTES_PER_COMMITMENT {
        return Err("commitment must be 48 bytes".to_string());
    }
    let p = zoda_bls::g1::G1Affine::from_compressed(b).ok_or("invalid G1 point")?;
    if !p.is_identity() {
        let r_mul = p
            .to_projective()
            .mul_limbs(&zoda_math::fr::FR_MODULUS);
        if !r_mul.is_identity() {
            return Err("G1 point not in subgroup".to_string());
        }
    }
    Ok(p)
}

pub fn bytes_to_kzg_proof(b: &[u8]) -> Result<zoda_bls::g1::G1Affine, String> {
    bytes_to_kzg_commitment(b)
}

/// Parse a canonical compressed G1 point WITHOUT the r-subgroup check —
/// for batch entry points that aggregate the checks via
/// [`batch_subgroup_check_g1`] (one Pippenger MSM + one full scalar
/// multiplication instead of n individual [r]P checks).
pub fn parse_g1_canonical(b: &[u8]) -> Result<zoda_bls::g1::G1Affine, String> {
    if b.len() != BYTES_PER_COMMITMENT {
        return Err("commitment must be 48 bytes".to_string());
    }
    zoda_bls::g1::G1Affine::from_compressed(b).ok_or("invalid G1 point".to_string())
}

/// Batch r-subgroup check with **small-subgroup-safe** coefficients:
/// `S = Σ cᵢ·Pᵢ` with `cᵢ = 1 + h₁·kᵢ` (kᵢ 128-bit random, from a
/// transcript of the compressed points), then a single `[r]S == O` test.
///
/// The `≡ 1 (mod h₁)` structure is essential: plain random coefficients
/// let a small-order point vanish from the combination whenever
/// `ord(P) | cᵢ` (probability ~1/ord — grindable for tiny orders), while
/// `cᵢ ≡ 1 (mod h₁)` guarantees every h₁-order component contributes
/// itself. Coefficients stay ~254 bits (h₁ ≈ 2^126, kᵢ 128-bit), so the
/// Pippenger MSM costs the same as with plain scalars, and one MSM plus
/// one full multiplication replaces n individual `[r]P` checks.
pub fn batch_subgroup_check_g1(points: &[zoda_bls::g1::G1Affine]) -> bool {
    if points.is_empty() {
        return true;
    }
    for p in points {
        if !p.infinity && !p.is_on_curve() {
            return false;
        }
    }
    let mut transcript = Vec::new();
    for p in points {
        transcript.extend_from_slice(&p.to_compressed());
    }
    transcript.extend_from_slice(&(points.len() as u64).to_be_bytes());
    let seed = zoda_math::sha256::sha256(&transcript);
    let mut rng = zoda_math::ZodaRng::from_seed(seed);

    // h₁ (3 limbs) and 2^256 mod r, both computed once
    let h1 = zoda_bls::endomorphism::g1_cofactor_h1();
    let two_pow_256 = {
        // 2^256 mod r via five squarings-style doublings: 2^256 = (2^32)^8
        // — simplest exact route: fold 1 << 256 through from_u64 doubling
        let mut v = Fr::from_u64(1);
        for _ in 0..256 {
            v = v + v;
        }
        v
    };

    // c = 1 + h₁·k with a fresh 128-bit k per point; c < 2^(126+128+1)
    // fits five limbs, reduced as c_lo (256 bits) + c[4]·2^256
    let coeffs: Vec<Fr> = (0..points.len())
        .map(|_| {
            let k0 = rng.next_u64();
            let k1 = rng.next_u64() & 0x7fff_ffff_ffff_ffff;
            // h₁·k = h₁·(k0 + k1·2^64), schoolbook into 6 limbs
            let mut c = [0u64; 6];
            for (src, kk) in [(0usize, k0), (1, k1)] {
                let mut carry = 0u128;
                for i in 0..3 {
                    let t = (h1[i] as u128) * (kk as u128) + carry + (c[i + src] as u128);
                    c[i + src] = t as u64;
                    carry = t >> 64;
                }
                let mut j = 3 + src;
                while carry > 0 && j < 6 {
                    let t = c[j] as u128 + carry;
                    c[j] = t as u64;
                    carry = t >> 64;
                    j += 1;
                }
            }
            // + 1 with carry propagation
            let mut i = 0;
            loop {
                let (v, ov) = c[i].overflowing_add(1);
                c[i] = v;
                if !ov {
                    break;
                }
                i += 1;
            }
            // reduce: c = c_lo + c[4]·2^256 (c[5] must be zero)
            debug_assert_eq!(c[5], 0);
            let lo_bytes: Vec<u8> = (0..4).rev().flat_map(|i| c[i].to_be_bytes()).collect();
            let lo = Fr::from_be_bytes_mod_order(&lo_bytes);
            let hi = Fr::from_u64(c[4]);
            lo + two_pow_256 * hi
        })
        .collect();
    let s = crate::msm::msm_jacobian(points, &coeffs);
    s.mul_limbs(&zoda_math::fr::FR_MODULUS).is_identity()
}

/// `verify_kzg_proof_impl`: e(P − y·g1, g2_neg) · e(π, [X]₂ − z·[1]₂) == 1.
pub fn verify_kzg_proof_impl(
    commitment: &zoda_bls::g1::G1Affine,
    z: Fr,
    y: Fr,
    proof: &zoda_bls::g1::G1Affine,
    setup: &Setup,
) -> bool {
    // X_minus_z = G2_SETUP[1] - z·G2gen
    let z_neg = -z;
    let x_minus_z = setup.g2_monomial[1]
        .to_projective()
        .add(&G2ProjectiveGen().mul_fr(&z_neg))
        .to_affine();
    // P_minus_y = commitment - y·g1
    let y_neg = -y;
    let p_minus_y = commitment
        .to_projective()
        .add(&G1Projective::generator().mul_fr(&y_neg))
        .to_affine();
    // pairing check: e(P−y, −g2) · e(π, X−z) == 1
    let neg_g2 = G2Prepared::from((-setup.g2_monomial[0]).to_projective().to_affine());
    pairing_check(&[
        (&p_minus_y, &neg_g2),
        (&proof, &G2Prepared::from(x_minus_z)),
    ])
}

use zoda_bls::g2::G2Projective as G2P2;
#[allow(non_snake_case)]
fn G2ProjectiveGen() -> G2P2 {
    G2P2::generator()
}

/// `verify_kzg_proof_batch`.
pub fn verify_kzg_proof_batch(
    commitments: &[zoda_bls::g1::G1Affine],
    zs: &[Fr],
    ys: &[Fr],
    proofs: &[zoda_bls::g1::G1Affine],
    setup: &Setup,
) -> bool {
    assert_eq!(commitments.len(), zs.len());
    assert_eq!(commitments.len(), ys.len());
    assert_eq!(commitments.len(), proofs.len());
    let n = commitments.len();
    if n == 0 {
        return true;
    }
    // random challenge r
    let mut data = Vec::with_capacity(16 + 8 + 8 + n * 128);
    data.extend_from_slice(RANDOM_CHALLENGE_KZG_BATCH_DOMAIN);
    data.extend_from_slice(&(FIELD_ELEMENTS_PER_BLOB as u64).to_be_bytes());
    data.extend_from_slice(&(n as u64).to_be_bytes());
    for i in 0..n {
        data.extend_from_slice(&to_bytes48(commitments[i]));
        data.extend_from_slice(&fr_to_bytes(zs[i]));
        data.extend_from_slice(&fr_to_bytes(ys[i]));
        data.extend_from_slice(&to_bytes48(proofs[i]));
    }
    let r = hash_to_bls_field(&data);
    // r_powers
    let mut r_powers = Vec::with_capacity(n);
    let mut cur = Fr::ONE;
    for _ in 0..n {
        r_powers.push(cur);
        cur = cur * r;
    }
    // proof_lincomb = Σ r^i π_i
    let proof_lincomb = msm::g1_lincomb(proofs, &r_powers);
    // proof_z_lincomb = Σ r^i z_i π_i
    let zr: Vec<Fr> = zs
        .iter()
        .zip(r_powers.iter())
        .map(|(z, rp)| *z * *rp)
        .collect();
    let proof_z_lincomb = msm::g1_lincomb(proofs, &zr);
    // C_minus_ys
    let c_minus: Vec<zoda_bls::g1::G1Affine> = commitments
        .iter()
        .zip(ys.iter())
        .map(|(c, y)| {
            c.to_projective()
                .add(&G1Projective::generator().mul_fr(&(-*y)))
                .to_affine()
        })
        .collect();
    let c_minus_lincomb = msm::g1_lincomb(&c_minus, &r_powers);
    // lhs = c_minus_lincomb + proof_z_lincomb
    let lhs = c_minus_lincomb
        .to_projective()
        .add(&proof_z_lincomb.to_projective())
        .to_affine();
    // e(Σr^i π_i, -[s]₂) · e(lhs, [1]₂) == 1
    let neg_s_g2 = G2Prepared::from((-setup.g2_monomial[1]).to_projective().to_affine());
    let g2_gen = G2Prepared::from(setup.g2_monomial[0]);
    pairing_check(&[(&proof_lincomb, &neg_s_g2), (&lhs, &g2_gen)])
}

/// `compute_blob_kzg_proof` (public method).
pub fn compute_blob_kzg_proof(
    blob: &[u8],
    commitment_bytes: &[u8],
    setup: &Setup,
) -> Result<[u8; 48], String> {
    if blob.len() != BLOB_BYTES {
        return Err("blob must be 131072 bytes".to_string());
    }
    let _ = bytes_to_kzg_commitment(commitment_bytes)?;
    let polynomial = blob_to_polynomial(blob)?;
    let evaluation_challenge = compute_challenge(blob, &to_arr48(commitment_bytes)?);
    let (proof, _) = compute_kzg_proof_impl(&polynomial, evaluation_challenge, setup);
    Ok(to_bytes48(proof))
}

/// `verify_blob_kzg_proof` (public method).
pub fn verify_blob_kzg_proof(
    blob: &[u8],
    commitment_bytes: &[u8],
    proof_bytes: &[u8],
    setup: &Setup,
) -> Result<bool, String> {
    if blob.len() != BLOB_BYTES {
        return Err("blob must be 131072 bytes".to_string());
    }
    let commitment = bytes_to_kzg_commitment(commitment_bytes)?;
    let polynomial = blob_to_polynomial(blob)?;
    let evaluation_challenge = compute_challenge(blob, &to_arr48(commitment_bytes)?);
    let y = evaluate_polynomial_in_evaluation_form(&polynomial, evaluation_challenge, setup);
    let proof = bytes_to_kzg_proof(proof_bytes)?;
    Ok(verify_kzg_proof_impl(
        &commitment,
        evaluation_challenge,
        y,
        &proof,
        setup,
    ))
}

/// `verify_blob_kzg_proof_batch` (public method).
pub fn verify_blob_kzg_proof_batch(
    blobs: &[&[u8]],
    commitments_bytes: &[&[u8]],
    proofs_bytes: &[&[u8]],
    setup: &Setup,
) -> Result<bool, String> {
    if blobs.len() != commitments_bytes.len() || blobs.len() != proofs_bytes.len() {
        return Err("input length mismatch".to_string());
    }
    let mut commitments = Vec::with_capacity(blobs.len());
    let mut zs = Vec::with_capacity(blobs.len());
    let mut ys = Vec::with_capacity(blobs.len());
    let mut proofs = Vec::with_capacity(blobs.len());
    for i in 0..blobs.len() {
        if blobs[i].len() != BLOB_BYTES {
            return Err("blob must be 131072 bytes".to_string());
        }
        // subgroup checks are aggregated at the end (one MSM + one [r]P)
        let commitment = parse_g1_canonical(commitments_bytes[i])?;
        let polynomial = blob_to_polynomial(blobs[i])?;
        let z = compute_challenge(blobs[i], &to_arr48(commitments_bytes[i])?);
        let y = evaluate_polynomial_in_evaluation_form(&polynomial, z, setup);
        commitments.push(commitment);
        zs.push(z);
        ys.push(y);
        proofs.push(parse_g1_canonical(proofs_bytes[i])?);
    }
    let mut all_points = commitments.clone();
    all_points.extend_from_slice(&proofs);
    if !batch_subgroup_check_g1(&all_points) {
        return Err("G1 point not in subgroup".to_string());
    }
    Ok(verify_kzg_proof_batch(&commitments, &zs, &ys, &proofs, setup))
}

/// The EIP-4844 versioned hash: sha256(commitment)[1..] with 0x01 version.
pub fn kzg_to_versioned_hash(commitment: &[u8; 48]) -> [u8; 32] {
    let h = sha256(commitment);
    let mut out = [0u8; 32];
    out[0] = VERSIONED_HASH_VERSION_KZG;
    out[1..].copy_from_slice(&h[1..]);
    out
}

fn to_bytes48(p: zoda_bls::g1::G1Affine) -> [u8; 48] {
    let v = p.to_compressed();
    let mut out = [0u8; 48];
    out.copy_from_slice(&v);
    out
}

fn to_arr48(b: &[u8]) -> Result<[u8; 48], String> {
    if b.len() != 48 {
        return Err("commitment must be 48 bytes".to_string());
    }
    let mut out = [0u8; 48];
    out.copy_from_slice(b);
    Ok(out)
}
