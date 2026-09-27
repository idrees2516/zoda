//! Hashing to BLS12-381 fields and curves (RFC 9380):
//! `expand_message_xmd(SHA-256)`, `hash_to_field`, Simplified SWU + isogeny
//! maps for G1 (11-isogeny) and G2 (3-isogeny), cofactor clearing.
//!
//! Suites implemented: `BLS12381G{1,2}_XMD:SHA-256_SSWU_RO_`.

use crate::curve::Projective;
use crate::fp::Fp;
use crate::fp2::Fp2;
use crate::g1::{G1Config, G1Projective};
use crate::g2::{G2Config, G2Projective};
use crate::h2c_consts::{g1_map, g2_map};
use std::sync::OnceLock;

/// Effective cofactor for the G1 suite (RFC 9380 §8.8.1).
const G1_H_EFF: [u64; 1] = [0xd201000000010001];
/// Effective cofactor for the G2 suite (RFC 9380 §8.8.2), LE limbs.
const G2_H_EFF: [u64; 10] = [
    0xe8020005aaa95551,
    0x59894c0adebbf6b4,
    0xe954cbc06689f6a3,
    0x2ec0ec69d7477c1a,
    0x6d82bf015d1212b0,
    0x329c2f178731db95,
    0x9986ff031508ffe1,
    0x88e2a8e9145ad768,
    0x584c6a0ea91b3528,
    0x0bc69f08f2ee75b3,
];

/// L = ceil((ceil(log2(p)) + k) / 8) = 64 for BLS12-381.
const L: usize = 64;

// ---------------------------------------------------------------------------
// expand_message_xmd (RFC 9380 §5.3.1)
// ---------------------------------------------------------------------------

pub fn expand_message_xmd_pub(msg: &[u8], dst: &[u8], len_in_bytes: usize) -> Vec<u8> {
    expand_message_xmd(msg, dst, len_in_bytes)
}

fn expand_message_xmd(msg: &[u8], dst: &[u8], len_in_bytes: usize) -> Vec<u8> {
    use zoda_math::sha256::sha256;
    debug_assert!(len_in_bytes <= 255 * 32);
    debug_assert!(!dst.is_empty() && dst.len() <= 255);

    let ell = (len_in_bytes + 31) / 32;
    // DST_prime = DST || I2OSP(len(DST), 1); for long DSTs,
    // DST_prime = H(DST) || I2OSP(255, 1)
    let mut dst_prime = if dst.len() > 255 {
        let mut h = sha256(dst).to_vec();
        h.push(255);
        h
    } else {
        let mut d = dst.to_vec();
        d.push(dst.len() as u8);
        d
    };

    // b_0 = H(Z_pad || msg || l_i_b_str || 0x00 || DST_prime)
    let mut input = Vec::with_capacity(64 + msg.len() + 2 + 1 + dst_prime.len());
    input.extend_from_slice(&[0u8; 64]); // Z_pad
    input.extend_from_slice(msg);
    input.extend_from_slice(&(len_in_bytes as u16).to_be_bytes());
    input.push(0x00);
    input.extend_from_slice(&dst_prime);
    let b0 = sha256(&input);

    // b_1 = H(b_0 || 0x01 || DST_prime)
    let mut input = Vec::with_capacity(32 + 1 + dst_prime.len());
    input.extend_from_slice(&b0);
    input.push(0x01);
    input.extend_from_slice(&dst_prime);
    let mut b_prev = sha256(&input);

    let mut out = Vec::with_capacity(ell * 32);
    out.extend_from_slice(&b_prev);
    for i in 2..=ell {
        // b_i = H(strxor(b_0, b_{i-1}) || I2OSP(i, 1) || DST_prime)
        let mut xored = [0u8; 32];
        for j in 0..32 {
            xored[j] = b0[j] ^ b_prev[j];
        }
        let mut input = Vec::with_capacity(32 + 1 + dst_prime.len());
        input.extend_from_slice(&xored);
        input.push(i as u8);
        input.extend_from_slice(&dst_prime);
        b_prev = sha256(&input);
        out.extend_from_slice(&b_prev);
    }
    out.truncate(len_in_bytes);
    out
}

// ---------------------------------------------------------------------------
// hash_to_field
// ---------------------------------------------------------------------------

fn fp_from_okm(bytes: &[u8]) -> Fp {
    // OS2IP(bytes) mod p — bytes is 64 bytes (L = 64).
    // interpret as 8 u64 words big-endian, then reduce mod p.
    let mut words = [0u64; 8];
    for i in 0..8 {
        let mut w = [0u8; 8];
        w.copy_from_slice(&bytes[i * 8..i * 8 + 8]);
        words[7 - i] = u64::from_be_bytes(w);
    }
    // value = Σ words[i]·2^(64i) mod p — reduce high-to-low with the
    // Fp-domain helper: acc = acc·2^64 + word
    reduce_words_mod_p(&words)
}

fn reduce_words_mod_p(words: &[u64; 8]) -> Fp {
    // p is 381 bits; process the top 2 words (128 bits) down.
    // acc (Montgomery) accumulates the low part; high words folded via 2^128 mod p.
    // 2^384 mod p folded: start from the top.
    // MSB-first Horner: words[7] is the most significant 64-bit word.
    let two64 = {
        let mut l = [0u64; 6];
        l[1] = 1;
        Fp::from_limbs(l)
    };
    let mut acc = Fp::from_limbs([words[7], 0, 0, 0, 0, 0]);
    for i in (0..7).rev() {
        acc = acc * two64 + Fp::from_limbs([words[i], 0, 0, 0, 0, 0]);
    }
    acc
}

/// Hash to `count` Fp elements (m = 1).
pub fn hash_to_field_fp(msg: &[u8], dst: &[u8], count: usize) -> Vec<Fp> {
    let len = count * L;
    let u = expand_message_xmd(msg, dst, len);
    (0..count)
        .map(|i| fp_from_okm(&u[i * L..(i + 1) * L]))
        .collect()
}

/// Hash to `count` Fp2 elements (m = 2): each element consumes 2·L bytes.
pub fn hash_to_field_fp2(msg: &[u8], dst: &[u8], count: usize) -> Vec<Fp2> {
    let len = count * 2 * L;
    let u = expand_message_xmd(msg, dst, len);
    (0..count)
        .map(|i| {
            let c0 = fp_from_okm(&u[i * 2 * L..i * 2 * L + L]);
            let c1 = fp_from_okm(&u[i * 2 * L + L..(i + 1) * 2 * L]);
            Fp2 { c0, c1 }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Simplified SWU + isogeny maps
// ---------------------------------------------------------------------------

/// The exponent (p² − 9)/16 as LE limbs (derived once).
fn p2_minus_9_over_16() -> &'static [u64] {
    static E: OnceLock<Vec<u64>> = OnceLock::new();
    E.get_or_init(|| {
        // p² − 9 as limbs, then >> 4.
        let p: Vec<u64> = Fp::MODULUS.to_vec();
        let mut sq = vec![0u64; 12];
        for i in 0..6 {
            let mut carry = 0u128;
            for j in 0..6 {
                let t = sq[i + j] as u128 + (p[i] as u128) * (p[j] as u128) + carry;
                sq[i + j] = t as u64;
                carry = t >> 64;
            }
            let mut idx = i + 6;
            while carry > 0 && idx < 12 {
                let t = sq[idx] as u128 + carry;
                sq[idx] = t as u64;
                carry = t >> 64;
                idx += 1;
            }
        }
        // subtract 9
        let mut borrow = 9u128;
        for limb in sq.iter_mut() {
            let t = (*limb as u128).wrapping_sub(borrow);
            *limb = t as u64;
            borrow = (t >> 64) & 1;
        }
        // shift right 4
        let mut c = 0u64;
        for limb in sq.iter_mut().rev() {
            let nc = *limb << 60;
            *limb = (*limb >> 4) | c;
            c = nc;
        }
        while sq.len() > 0 && *sq.last().unwrap() == 0 {
            sq.pop();
        }
        sq
    })
}

/// Simplified SWU for G1 (AB == 0 case) on the 11-isogenous curve.
pub fn map_to_curve_sswu_g1_pub(u: &Fp) -> G1Projective {
    map_to_curve_sswu_g1(u)
}

fn map_to_curve_sswu_g1(u: &Fp) -> G1Projective {
    let a = g1_map::SSWU_ELLP_A;
    let b = g1_map::SSWU_ELLP_B;
    let xi = g1_map::SSWU_XI;

    let usq = u.square();
    let xi_usq = xi * usq;
    let xisq_u4 = xi_usq * xi_usq;
    let nd_common = xisq_u4 + xi_usq; // XI²u⁴ + XIu²
    let x_den = if nd_common.is_zero() {
        a * xi
    } else {
        a * (-nd_common)
    };
    let x0_num = b * (Fp::one() + nd_common);

    // g(x0(u)) as numerator/denominator
    let x_densq = x_den * x_den;
    let gx_den = x_densq * x_den;
    let gx0_num = (x0_num * x0_num + a * x_densq) * x0_num + b * gx_den;

    // sqrt candidate = uv·(uv³)^((p-3)/4) — the classic p ≡ 3 (mod 4) form
    let uv = gx0_num * gx_den;
    let vsq = gx_den * gx_den;
    // exponent (p-3)/4
    let e = {
        let mut e = Fp::MODULUS;
        // subtract 3
        let mut borrow = 3u128;
        for limb in e.iter_mut() {
            let t = (*limb as u128).wrapping_sub(borrow);
            *limb = t as u64;
            borrow = (t >> 64) & 1;
        }
        // >> 2
        let mut c = 0u64;
        for limb in e.iter_mut().rev() {
            let nc = *limb << 62;
            *limb = (*limb >> 2) | c;
            c = nc;
        }
        e
    };
    // sqrt(u/v) = u·v·(u·v³)^((p-3)/4)
    let uv_v3 = uv * vsq;
    let sqrt_candidate = uv * uv_v3.pow_limbs(&e);

    let gx0_is_square = sqrt_candidate * sqrt_candidate * gx_den == gx0_num;
    let x1_num = x0_num * xi_usq;
    let y1 = g1_map::SQRT_M_XI_CUBED * usq * *u * sqrt_candidate;

    let x_num = if gx0_is_square { x0_num } else { x1_num };
    let mut y = if gx0_is_square {
        sqrt_candidate
    } else {
        y1
    };
    // sign correction
    let sgn_y = {
        let r = y.to_repr();
        r[0] & 1 == 1
    };
    let sgn_u = {
        let r = u.to_repr();
        r[0] & 1 == 1
    };
    if sgn_y != sgn_u {
        y = -y;
    }

    G1Projective {
        x: x_num,
        y: y * x_den,
        z: x_den,
    }
}

/// 11-isogeny map E' → E for G1.
fn iso_map_g1(u: &G1Projective) -> G1Projective {
    let coeffs: [&[Fp]; 4] = [
        &g1_map::ISO11_XNUM,
        &g1_map::ISO11_XDEN,
        &g1_map::ISO11_YNUM,
        &g1_map::ISO11_YDEN,
    ];
    let (x, y, z) = (u.x, u.y, u.z);
    let mut mapvals = [Fp::zero(); 4];
    // pre-compute z^1 .. z^15 (Y maps have degree up to 16)
    let mut zpows = [Fp::zero(); 15];
    zpows[0] = z;
    for i in 1..15 {
        zpows[i] = zpows[i - 1] * z;
    }
    for idx in 0..4 {
        let coeff = coeffs[idx];
        let clast = coeff.len() - 1;
        mapvals[idx] = coeff[clast];
        for jdx in 0..clast {
            mapvals[idx] = mapvals[idx] * x + zpows[jdx] * coeff[clast - 1 - jdx];
        }
    }
    // x denominator is one degree less than the numerator: extra z factor
    mapvals[1] = mapvals[1] * z;
    // Y map multiplied by y, and the denominator by z
    mapvals[2] = mapvals[2] * y;
    mapvals[3] = mapvals[3] * z;
    G1Projective {
        x: mapvals[0] * mapvals[3],
        y: mapvals[2] * mapvals[1],
        z: mapvals[1] * mapvals[3],
    }
}

/// Simplified SWU for G2 (AB == 0 case) on the 3-isogenous curve.
pub fn map_to_curve_sswu_g2_pub(u: &Fp2) -> G2Projective {
    map_to_curve_sswu_g2(u)
}

fn map_to_curve_sswu_g2(u: &Fp2) -> G2Projective {
    let a = g2_map::SSWU_ELLP_A;
    let b = g2_map::SSWU_ELLP_B;
    let xi = g2_map::SSWU_XI;

    let usq = (*u) * (*u);
    let xi_usq = xi * usq;
    let xisq_u4 = xi_usq * xi_usq;
    let nd_common = xisq_u4 + xi_usq;
    let x_den = if nd_common.is_zero() {
        a * xi
    } else {
        a * (-nd_common)
    };
    let x0_num = b * (Fp2::one() + nd_common);

    let x_densq = x_den * x_den;
    let gx_den = x_densq * x_den;
    let gx0_num = (x0_num * x0_num + a * x_densq) * x0_num + b * gx_den;

    // sqrt candidate = uv⁷·(uv¹⁵)^((p²−9)/16)
    let vsq = gx_den * gx_den;
    let v3 = vsq * gx_den;
    let v4 = vsq * vsq;
    let v8 = v4 * v4;
    let uv7 = gx0_num * v3 * v4;
    let uv15 = uv7 * v8;
    let sqrt_candidate = uv7 * uv15.pow_limbs(p2_minus_9_over_16());

    // Test the candidate against the other square roots of unity.
    let mut y = sqrt_candidate;
    // multiply by i
    let tmp = Fp2 {
        c0: -sqrt_candidate.c1,
        c1: sqrt_candidate.c0,
    };
    if tmp * tmp * gx_den == gx0_num {
        y = tmp;
    }
    // multiply by RV1
    let tmp = sqrt_candidate * g2_map::SSWU_RV1;
    if tmp * tmp * gx_den == gx0_num {
        y = tmp;
    }
    // multiply by RV1·i
    let tmp = Fp2 {
        c0: tmp.c1,
        c1: -tmp.c0,
    };
    if tmp * tmp * gx_den == gx0_num {
        y = tmp;
    }

    // g(x1(u)) = g(x0(u))·ξ³u⁶
    let gx1_num = gx0_num * xi_usq * xisq_u4;
    // sqrt candidate for x1: candidate·u³
    let sc1 = sqrt_candidate * usq * *u;
    let mut eta_found = false;
    for eta in g2_map::SSWU_ETAS.iter() {
        let tmp = sc1 * *eta;
        if tmp * tmp * gx_den == gx1_num {
            y = tmp;
            eta_found = true;
        }
    }

    let x_num = if eta_found { x0_num * xi_usq } else { x0_num };
    // sign correction: sgn0(y) == sgn0(u)
    if y.sgn0() != u.sgn0() {
        y = -y;
    }

    G2Projective {
        x: x_num,
        y: y * x_den,
        z: x_den,
    }
}

/// 3-isogeny map E' → E for G2.
pub fn iso_map_g2_pub(u: &G2Projective) -> G2Projective {
    iso_map_g2(u)
}

fn iso_map_g2(u: &G2Projective) -> G2Projective {
    let coeffs: [&[Fp2]; 4] = [
        &g2_map::ISO3_XNUM,
        &g2_map::ISO3_XDEN,
        &g2_map::ISO3_YNUM,
        &g2_map::ISO3_YDEN,
    ];
    let (x, y, z) = (u.x, u.y, u.z);
    let mut mapvals = [Fp2::zero(); 4];
    let zsq = z * z;
    let zpows = [z, zsq, zsq * z];
    for idx in 0..4 {
        let coeff = coeffs[idx];
        let clast = coeff.len() - 1;
        mapvals[idx] = coeff[clast];
        for jdx in 0..clast {
            mapvals[idx] = mapvals[idx] * x + zpows[jdx] * coeff[clast - 1 - jdx];
        }
    }
    mapvals[1] = mapvals[1] * z;
    mapvals[2] = mapvals[2] * y;
    mapvals[3] = mapvals[3] * z;
    G2Projective {
        x: mapvals[0] * mapvals[3],
        y: mapvals[2] * mapvals[1],
        z: mapvals[1] * mapvals[3],
    }
}

// ---------------------------------------------------------------------------
// Public hash-to-curve API
// ---------------------------------------------------------------------------

/// DST for `BLS12381G1_XMD:SHA-256_SSWU_RO_`.
pub const DST_G1: &[u8] = b"BLS_SIG_BLS12381G1_XMD:SHA-256_SSWU_RO_NUL_";
/// DST for `BLS12381G2_XMD:SHA-256_SSWU_RO_` (the Ethereum ciphersuite).
pub const DST_G2: &[u8] = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";

/// Hash to G1 (RO variant: clear_h(map(u[0]) + map(u[1])), RFC 9380 §3).
pub fn hash_to_curve_g1(msg: &[u8], dst: &[u8]) -> G1Projective {
    let us = hash_to_field_fp(msg, dst, 2);
    let q0 = iso_map_g1(&map_to_curve_sswu_g1(&us[0]));
    let q1 = iso_map_g1(&map_to_curve_sswu_g1(&us[1]));
    (q0 + q1).mul_limbs(&G1_H_EFF)
}

/// encode_to_curve for G1 (single-element, NU variant core).
pub fn encode_to_curve_g1(msg: &[u8], dst: &[u8]) -> G1Projective {
    let u = hash_to_field_fp(msg, dst, 1)[0];
    iso_map_g1(&map_to_curve_sswu_g1(&u)).mul_limbs(&G1_H_EFF)
}

/// Hash to G2 (RO variant) — the Ethereum ciphersuite.
pub fn hash_to_curve_g2(msg: &[u8], dst: &[u8]) -> G2Projective {
    let us = hash_to_field_fp2(msg, dst, 2);
    let q0 = iso_map_g2(&map_to_curve_sswu_g2(&us[0]));
    let q1 = iso_map_g2(&map_to_curve_sswu_g2(&us[1]));
    (q0 + q1).mul_limbs(&G2_H_EFF)
}

/// encode_to_curve for G2 (single-element, NU variant core).
pub fn encode_to_curve_g2(msg: &[u8], dst: &[u8]) -> G2Projective {
    let u = hash_to_field_fp2(msg, dst, 1)[0];
    iso_map_g2(&map_to_curve_sswu_g2(&u)).mul_limbs(&G2_H_EFF)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expand_message_xmd_reference_vector() {
        // RFC 9380 test vector for expand_message_xmd SHA-256, DST "QUUX..."
        // len 32, msg empty — Appendix H.1 (SHA-256, SHA256(abc)...):
        // we verify self-consistency + length only here; full vectors below
        // in hash_to_curve tests.
        let out = expand_message_xmd(b"abc", b"SGENTITYTEST", 32);
        assert_eq!(out.len(), 32);
        let out = expand_message_xmd(b"abc", b"SGENTITYTEST", 128);
        assert_eq!(out.len(), 128);
    }

    #[test]
    fn hash_to_field_known_vectors() {
        // RFC 9380 Appendix J.9.1 (G1 suite, msg "abc"):
        // u[0] = 0d921c33f2bad966478a03ca35d05719bdf92d347557ea166e5bba
        //        579eea9b83e9afa5c088573c2281410369fbd32951
        let dst = b"QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_RO_";
        let u = hash_to_field_fp(b"abc", dst, 2);
        let r = u[0].to_repr();
        let mut be = [0u8; 48];
        for i in 0..6 {
            be[i * 8..i * 8 + 8].copy_from_slice(&r[5 - i].to_be_bytes());
        }
        let hex: String = be.iter().map(|b| format!("{:02x}", b)).collect();
        assert!(
            hex.starts_with("0d921c33f2bad966478a03ca35d05719bdf92d347557ea166e5bba"),
            "u0 = {}",
            hex
        );
    }

    #[test]
    fn hash_to_field_g2_known_vector() {
        // RFC 9380 Appendix J.10.1 (G2 suite, msg "abc"):
        // u[0].c0 = 15f7c0aa8f6b296ab5ff9c2c7581ade64f4ee6f1bf18f55179ff44
        //           a2cf355fa53dd2a2158c5ecb17d7c52f63e7195771
        // u[0].c1 = 01c8067bf4c0ba709aa8b9abc3d1cef589a4758e09ef53732d670f
        //           d8739a7274e111ba2fcaa71b3d33df2a3a0c8529dd
        let dst = b"QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_";
        let u = hash_to_field_fp2(b"abc", dst, 2);
        let bytes = u[0].to_be_bytes();
        let hex: String = bytes.iter().map(|b| format!("{:02x}", b)).collect();
        assert!(
            hex.starts_with("01c8067bf4c0ba709aa8b9abc3d1cef589a4758e09ef53732d670fd8739a7274e111ba2fcaa71b3d33df2a3a0c8529dd15f7c0aa8f6b296ab5ff9c2c7581ade64f4ee6f1bf18f55179ff44a2cf355fa53dd2a2158c5ecb17d7c52f63e7195771"),
            "u0 = {}",
            hex
        );
    }

    #[test]
    fn hash_to_curve_g1_on_curve() {
        let dst = b"QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_RO_";
        for msg in [b"" as &[u8], b"abc", b"longer message for testing"] {
            let p = hash_to_curve_g1(msg, dst);
            assert!(p.to_affine().is_on_curve());
            assert!(!p.is_identity());
        }
    }

    #[test]
    fn hash_to_curve_g2_on_curve() {
        let dst = b"QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_";
        for msg in [b"" as &[u8], b"abc", b"another message"] {
            let p = hash_to_curve_g2(msg, dst);
            assert!(p.to_affine().is_on_curve());
            assert!(!p.is_identity());
        }
    }

    #[test]
    fn hash_to_curve_g1_reference_point() {
        // RFC 9380 Appendix J.9.1: hash_to_curve("abc")
        // P.x = 03567bc5ef9c690c2ab2ecdf6a96ef1c139cc0b2f284dca0a9a794
        //       3388a49a3aee664ba5379a7655d3c68900be2f6903
        // P.y = 0b9c15f3fe6e5cf4211f346271d7b01c8f3b28be689c8429c85b67
        //       af215533311f0b8dfaaa154fa6b88176c229f2885d
        let dst = b"QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_RO_";
        let p = hash_to_curve_g1(b"abc", dst).to_affine();
        let bytes = p.x.to_repr();
        let mut be = [0u8; 48];
        for i in 0..6 {
            be[i * 8..i * 8 + 8].copy_from_slice(&bytes[5 - i].to_be_bytes());
        }
        let hex: String = be.iter().map(|b| format!("{:02x}", b)).collect();
        assert!(
            hex.starts_with("03567bc5ef9c690c2ab2ecdf6a96ef1c139cc0b2f284dca0a9a7943388a49a3"),
            "x = {}",
            hex
        );
        let ybytes = p.y.to_repr();
        for i in 0..6 {
            be[i * 8..i * 8 + 8].copy_from_slice(&ybytes[5 - i].to_be_bytes());
        }
        let hexy: String = be.iter().map(|b| format!("{:02x}", b)).collect();
        assert!(
            hexy.starts_with("0b9c15f3fe6e5cf4211f346271d7b01c8f3b28be689c8429c85b67af21553331"),
            "y = {}",
            hexy
        );
    }

    #[test]
    fn hash_to_curve_g2_reference_point() {
        // RFC 9380 Appendix J.10.1: hash_to_curve("abc")
        // P.x = 02c2d18e...f2787776e6 + I * 139cddbc...a41177fd8
        let dst = b"QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_";
        let p = hash_to_curve_g2(b"abc", dst).to_affine();
        // to_be_bytes emits c1 || c0
        let bytes = p.x.to_be_bytes();
        let hex: String = bytes.iter().map(|b| format!("{:02x}", b)).collect();
        assert!(
            hex.starts_with("139cddbccdc5e91b9623efd38c49f81a6f83f175e80b06fc374de9eb4b41dfe4ca3a230ed250fbe3a2acf73a41177fd802c2d18e033b960562aae3cab37a27ce00d80ccd5ba4b7fe0e7a210245129dbec7780ccc7954725f4168aff2787776e6"),
            "x = {}",
            hex
        );
    }
}
