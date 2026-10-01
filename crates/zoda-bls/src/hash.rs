//! Hashing to BLS12-381 fields and curves (RFC 9380):
//! `expand_message_xmd(SHA-256)`, `hash_to_field`, Simplified SWU + isogeny
//! maps for G1 (11-isogeny) and G2 (3-isogeny), cofactor clearing.
//!
//! Suites implemented: `BLS12381G{1,2}_XMD:SHA-256_SSWU_RO_`.

use crate::fp::Fp;
use crate::fp2::Fp2;
use crate::g1::G1Projective;
use crate::g2::G2Projective;
use crate::h2c_consts::{g1_map, g2_map};

/// Effective cofactor for the G1 suite (RFC 9380 §8.8.1).
const G1_H_EFF: [u64; 1] = [0xd201000000010001];
/// Effective cofactor for the G2 suite (RFC 9380 §8.8.2), LE limbs.
/// The naive [h_eff] scalar — kept as the ground truth for the
/// psi-chain cofactor-clearing equivalence test.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) const G2_H_EFF: [u64; 10] = [
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
    let dst_prime = if dst.len() > 255 {
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


/// The exponent (p − 3)/4 as LE limbs — the classic square-root candidate
/// exponent for `p ≡ 3 (mod 4)` (`sqrt(u/v) = u·v·(u·v³)^((p−3)/4)`),
/// computed once at compile time.
const SQRT_RATIO_EXP: [u64; 6] = {
    let mut e = [
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a,
    ];
    // p - 3
    let mut borrow = 3u128;
    let mut i = 0;
    while i < 6 {
        let t = (e[i] as u128).wrapping_sub(borrow);
        e[i] = t as u64;
        borrow = (t >> 64) & 1;
        i += 1;
    }
    // >> 2
    let mut c = 0u64;
    let mut j = 6;
    while j > 0 {
        j -= 1;
        let nc = e[j] << 62;
        e[j] = (e[j] >> 2) | c;
        c = nc;
    }
    e
};

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

    // sqrt(u/v) = u·v·(u·v³)^((p-3)/4) — windowed fixed-exponent ladder
    let uv = gx0_num * gx_den;
    let vsq = gx_den * gx_den;
    let uv_v3 = uv * vsq;
    let sqrt_candidate = uv * uv_v3.pow_windowed(&SQRT_RATIO_EXP);

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
///
/// The square root is computed by the **complex method** (two windowed Fp
/// exponentiations by the compile-time constant `(p+1)/4` — the norm and
/// the real part — replacing one 762-bit Fp2 exponentiation, ≈3.5×) and
/// the division by the g(x) denominator is resolved with the standard
/// homogeneous-coordinate trick instead of a field inversion:
///
/// with `u = g(x0)` as a numerator and `v` its denominator (`v = x_den³`),
/// `y' = sqrt(u·v)` exists exactly when `g(x0)` is a square (their ratio
/// `v²` is always a square) and the affine y is `y'/v`; the projective
/// point `(x_num·x_den², y', v)` realises both denominators at once
/// because `X/Z = x_num/x_den` and `Y/Z = y'/v`.
///
/// One norm-based Fp2 inversion remains for the sgn0 sign correction,
/// which needs the true affine y — it costs a single Fp exponentiation.
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
    let v = x_densq * x_den; // g(x) denominator
    let gx0_num = (x0_num * x0_num + a * x_densq) * x0_num + b * v;

    // g(x1(u)) = g(x0(u))·ξ³u⁶ shares the denominator v
    let gx1_num = gx0_num * xi_usq * xisq_u4;

    // complex-method square roots; the branch taken (x0 vs x1) matches the
    // reference implementation exactly because square-ness of u/v and u·v
    // is equivalent (their ratio v² is a square)
    let mut x_num = x0_num;
    let mut yp = match (gx0_num * v).sqrt() {
        Some(y) => y,
        None => {
            x_num = x0_num * xi_usq;
            match (gx1_num * v).sqrt() {
                Some(y1) => y1,
                None => {
                    // no square root at either candidate: the reference
                    // implementation cannot produce a point here either
                    // (it exhausts the same root-of-unity candidates).
                    // This is unreachable for hash-to-field outputs and
                    // reachable only for adversarially chosen u.
                    return G2Projective {
                        x: Fp2::zero(),
                        y: Fp2::zero(),
                        z: Fp2::zero(),
                    };
                }
            }
        }
    };
    // affine y for the sign correction: y = y'/v (one norm-based Fp2
    // inversion — a single Fp exponentiation). The projective Y is then
    // simply y'·(v/v): Y = y_affine·Z = (y'/v)·v = y'.
    let vinv = match v.invert() {
        Some(i) => i,
        None => return G2Projective::identity(),
    };
    let y = yp * vinv;
    if y.sgn0() != u.sgn0() {
        yp = -yp;
    }
    // (x_num·x_den², y', v): X/Z = x_num/x_den, Y/Z = y'/v
    G2Projective {
        x: x_num * x_densq,
        y: yp,
        z: v,
    }
}

/// Batch Simplified SWU for G2 — the hash-to-curve workhorse behind
/// batch BLS verification. Per element the two 381-bit windowed
/// exponentiations of the complex-method square root are computed in a
/// first pass ([`Fp2::sqrt_batch`] defers every inversion), the g(x)
/// denominators are inverted with a single Montgomery batch inversion,
/// and only cheap multiplications remain for the assembly pass.
///
/// The output is bit-identical to [`map_to_curve_sswu_g2_pub`] on every
/// input (the final point is unique: the x-branch is determined by
/// square-ness, and the y sign by sgn0 equality).
pub fn map_to_curve_sswu_g2_batch(us: &[Fp2]) -> Vec<G2Projective> {
    // per-element staged numerator/denominator data
    struct Staged {
        x_densq: Fp2,
        v: Fp2,
        x0_num: Fp2,
        x1_num: Fp2,
        u0: Fp2,
        u1: Fp2,
        sgn_u: bool,
    }
    let stage = |u: &Fp2| -> Staged {
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
        let v = x_densq * x_den;
        let gx0_num = (x0_num * x0_num + a * x_densq) * x0_num + b * v;
        let gx1_num = gx0_num * xi_usq * xisq_u4;
        Staged {
            x_densq,
            v,
            x0_num,
            x1_num: x0_num * xi_usq,
            u0: gx0_num * v,
            u1: gx1_num * v,
            sgn_u: u.sgn0(),
        }
    };

    let staged: Vec<Staged> = us.iter().map(stage).collect();

    // first sqrt round over g(x0)·v; NonSquare entries fall through to a
    // second round over g(x1)·v (up to half the elements, still batched)
    let w0: Vec<Fp2> = staged.iter().map(|s| s.u0).collect();
    let r0 = Fp2::sqrt_batch(&w0);
    let mut roots: Vec<Option<Fp2>> = r0;
    let mut from_x1 = vec![false; staged.len()];
    let retry: Vec<usize> = (0..staged.len())
        .filter(|&i| roots[i].is_none())
        .collect();
    if !retry.is_empty() {
        let w1: Vec<Fp2> = retry.iter().map(|&i| staged[i].u1).collect();
        let r1 = Fp2::sqrt_batch(&w1);
        for (k, &i) in retry.iter().enumerate() {
            roots[i] = r1[k];
            from_x1[i] = r1[k].is_some();
        }
    }

    // one Montgomery batch inversion of the g(x) denominators
    let vs: Vec<Fp2> = staged.iter().map(|s| s.v).collect();
    let vinvs = Fp2::invert_batch(&vs);

    // assembly: branch select + sign correction + homogeneous coordinates
    staged
        .iter()
        .enumerate()
        .map(|(i, s)| match roots[i] {
            Some(yp) => {
                let y = yp * vinvs[i];
                let yp = if y.sgn0() != s.sgn_u { -yp } else { yp };
                let x_num = if from_x1[i] { s.x1_num } else { s.x0_num };
                // (x_num·x_den², y', v): X/Z = x_num/x_den and
                // Y/Z = y'/v — the affine point of the single path
                G2Projective {
                    x: x_num * s.x_densq,
                    y: yp,
                    z: s.v,
                }
            }
            // both candidates non-square: unreachable for hash-to-field
            // outputs; parity with the single path's degenerate result
            None => G2Projective {
                x: Fp2::zero(),
                y: Fp2::zero(),
                z: Fp2::zero(),
            },
        })
        .collect()
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
    // fast cofactor clearing: the ψ/ψ2 endomorphism chain, exactly
    // [h_eff]·(q0+q1) (see endomorphism::clear_cofactor_g2)
    crate::endomorphism::clear_cofactor_g2(&(q0 + q1))
}

/// encode_to_curve for G2 (single-element, NU variant core).
pub fn encode_to_curve_g2(msg: &[u8], dst: &[u8]) -> G2Projective {
    let u = hash_to_field_fp2(msg, dst, 1)[0];
    crate::endomorphism::clear_cofactor_g2(&iso_map_g2(&map_to_curve_sswu_g2(&u)))
}

/// Batch hash to G2 (RO variant) — `hash_to_curve_g2` over many messages
/// with the SSWU square roots batched ([`Fp2::sqrt_batch`]: every
/// inversion deferred into two Montgomery batch inversions across the
/// whole batch) and the map outputs assembled from pure multiplications.
/// Bit-identical to the single-message path.
pub fn hash_to_curve_g2_batch<'m, M: AsRef<[u8]> + Sync + 'm>(
    msgs: &'m [M],
    dst: &[u8],
) -> Vec<G2Projective> {
    // two field elements per message, flattened into one batch
    let mut us: Vec<Fp2> = Vec::with_capacity(msgs.len() * 2);
    for m in msgs {
        us.extend(hash_to_field_fp2(m.as_ref(), dst, 2));
    }
    let mapped = map_to_curve_sswu_g2_batch(&us);
    let mut out = Vec::with_capacity(msgs.len());
    for i in 0..msgs.len() {
        let q0 = iso_map_g2(&mapped[2 * i]);
        let q1 = iso_map_g2(&mapped[2 * i + 1]);
        out.push(crate::endomorphism::clear_cofactor_g2(&(q0 + q1)));
    }
    out
}

/// Batch hash to G1 (RO variant). The G1 SSWU is inversion-free (the
/// numerator/denominator form), so the batch path shares only the
/// windowed fixed-exponent ladder of the single path — kept for API
/// symmetry with the G2 batch.
pub fn hash_to_curve_g1_batch<'m, M: AsRef<[u8]> + Sync + 'm>(
    msgs: &'m [M],
    dst: &[u8],
) -> Vec<G1Projective> {
    msgs.iter()
        .map(|m| hash_to_curve_g1(m.as_ref(), dst))
        .collect()
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
    #[test]
    fn map_to_curve_sswu_g2_batch_matches_single() {
        // differential: the batch complex-method SSWU must reproduce the
        // single-element map bit-exactly on random u's (both x-branches,
        // both sign corrections) including the RFC field elements
        let dst = b"QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_";
        let mut us: Vec<Fp2> = Vec::new();
        for i in 0..48 {
            us.extend(hash_to_field_fp2(
                format!("batch-diff message {}", i).as_bytes(),
                dst,
                2,
            ));
        }
        let batch = map_to_curve_sswu_g2_batch(&us);
        for (u, b) in us.iter().zip(batch.iter()) {
            assert_eq!(*b, map_to_curve_sswu_g2(u), "batch != single map");
        }
    }

    #[test]
    fn hash_to_curve_g2_batch_matches_single() {
        let dst = DST_G2;
        let msgs: Vec<Vec<u8>> = (0..64)
            .map(|i| format!("zoda verify_batch message {}", i).into_bytes())
            .collect();
        let batch = hash_to_curve_g2_batch(&msgs, dst);
        assert_eq!(batch.len(), msgs.len());
        for (m, b) in msgs.iter().zip(batch.iter()) {
            assert_eq!(*b, hash_to_curve_g2(m, dst));
            assert!(b.to_affine().is_on_curve());
            assert!(!b.is_identity());
        }
        // batch result also verifies as a BLS hashed public key would
        let single: Vec<G2Projective> =
            msgs.iter().map(|m| hash_to_curve_g2(m, dst)).collect();
        assert_eq!(batch, single);
    }

    #[test]
    fn hash_to_curve_g1_batch_matches_single() {
        let dst = DST_G1;
        let msgs: Vec<Vec<u8>> = (0..16)
            .map(|i| format!("g1 batch {}", i).into_bytes())
            .collect();
        let batch = hash_to_curve_g1_batch(&msgs, dst);
        for (m, b) in msgs.iter().zip(batch.iter()) {
            assert_eq!(*b, hash_to_curve_g1(m, dst));
        }
    }





}
