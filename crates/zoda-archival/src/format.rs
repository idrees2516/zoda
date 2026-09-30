//! Compact, checksummed custody file formats.
//!
//! Everything is little-endian, fixed-width, and self-describing: the
//! public verification parameters are embedded so a loaded file can be
//! re-verified standalone, and a trailing SHA-256 covers the whole
//! payload.

use crate::store::GridCustody;
use zoda_core::{ZodaParams, ZodaPublic};
use zoda_math::sha256::sha256;
use zoda_math::{Fr, PrimeField};

const GRID_MAGIC: &[u8; 8] = b"ZODACST1";
const FORMAT_VERSION: u32 = 1;

/// Exact size of a grid custody file for an `m × k` grid with `present`
/// of `2k` columns stored.
pub fn grid_file_size(m: usize, k: usize, present: usize) -> usize {
    let ext_rows = 2 * m;
    let ext_cols = 2 * k;
    let frs = 2 * ext_cols + 2 * ext_rows; // g_r, g_r2, z_r, z_r2
    8 // magic
        + 4 // version
        + 32 // grid id
        + 8 // slot
        + 8 // m
        + 8 // k
        + frs * 32
        + 32 // row root
        + 32 // col root
        + (ext_cols + 7) / 8 // bitfield
        + present * ext_rows * 32 // cells
        + 32 // checksum
}

fn put_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_fr(buf: &mut Vec<u8>, f: Fr) {
    buf.extend_from_slice(&f.to_le_bytes());
}

fn take_u32<'a>(cur: &mut &'a [u8]) -> Result<u32, String> {
    if cur.len() < 4 {
        return Err("truncated u32".into());
    }
    let v = u32::from_le_bytes(cur[..4].try_into().unwrap());
    *cur = &cur[4..];
    Ok(v)
}

fn take_u64<'a>(cur: &mut &'a [u8]) -> Result<u64, String> {
    if cur.len() < 8 {
        return Err("truncated u64".into());
    }
    let v = u64::from_le_bytes(cur[..8].try_into().unwrap());
    *cur = &cur[8..];
    Ok(v)
}

fn take_bytes<'a>(cur: &mut &'a [u8], n: usize) -> Result<&'a [u8], String> {
    if cur.len() < n {
        return Err(format!("truncated {}-byte field", n));
    }
    let out = &cur[..n];
    *cur = &cur[n..];
    Ok(out)
}

fn take_fr(cur: &mut &[u8]) -> Result<Fr, String> {
    let b = take_bytes(cur, 32)?;
    let le: [u8; 32] = b.try_into().unwrap();
    let f = Fr::from_le_bytes_mod_order(&le);
    if f.to_le_bytes() != le {
        return Err("non-canonical field element".into());
    }
    Ok(f)
}

fn take_frs(cur: &mut &[u8], n: usize) -> Result<Vec<Fr>, String> {
    (0..n).map(|_| take_fr(cur)).collect()
}

/// Serialize a grid custody into the compact format.
pub fn encode_grid(custody: &GridCustody) -> Vec<u8> {
    let m = custody.params.m;
    let k = custody.params.k;
    let ext_cols = 2 * k;
    let mut buf = Vec::with_capacity(grid_file_size(m, k, custody.column_count()));
    buf.extend_from_slice(GRID_MAGIC);
    put_u32(&mut buf, FORMAT_VERSION);
    buf.extend_from_slice(&custody.grid_id);
    put_u64(&mut buf, custody.slot);
    put_u64(&mut buf, m as u64);
    put_u64(&mut buf, k as u64);
    for f in &custody.public.g_r {
        put_fr(&mut buf, *f);
    }
    for f in &custody.public.g_r2 {
        put_fr(&mut buf, *f);
    }
    for f in &custody.public.z_r {
        put_fr(&mut buf, *f);
    }
    for f in &custody.public.z_r2 {
        put_fr(&mut buf, *f);
    }
    buf.extend_from_slice(&custody.public.row_root);
    buf.extend_from_slice(&custody.public.col_root);
    buf.extend_from_slice(&custody.column_bitfield());
    for c in 0..ext_cols {
        if let Some(col) = custody.get_column(c) {
            for f in col {
                put_fr(&mut buf, *f);
            }
        }
    }
    let checksum = sha256(&buf);
    buf.extend_from_slice(&checksum);
    buf
}

/// Parse and checksum-verify a grid custody file. The returned custody is
/// re-verifiable standalone (public parameters embedded).
pub fn decode_grid(bytes: &[u8]) -> Result<GridCustody, String> {
    if bytes.len() < 8 + 4 + 32 + 8 + 8 + 8 + 32 {
        return Err("file too small".into());
    }
    // checksum first
    let split = bytes.len() - 32;
    let expected = sha256(&bytes[..split]);
    if expected[..] != bytes[split..] {
        return Err("checksum mismatch".into());
    }
    let mut cur: &[u8] = &bytes[..split];
    let magic = take_bytes(&mut cur, 8)?;
    if magic != GRID_MAGIC {
        return Err("not a ZODA custody file".into());
    }
    let version = take_u32(&mut cur)?;
    if version != FORMAT_VERSION {
        return Err(format!("unsupported format version {}", version));
    }
    let grid_id: [u8; 32] = take_bytes(&mut cur, 32)?.try_into().unwrap();
    let slot = take_u64(&mut cur)?;
    let m = take_u64(&mut cur)? as usize;
    let k = take_u64(&mut cur)? as usize;
    if m == 0 || k == 0 || !m.is_power_of_two() || !k.is_power_of_two() {
        return Err("invalid grid shape".into());
    }
    let ext_rows = 2 * m;
    let ext_cols = 2 * k;
    let _ = ext_rows;
    let g_r = take_frs(&mut cur, ext_cols)?;
    let g_r2 = take_frs(&mut cur, ext_rows)?;
    let z_r = take_frs(&mut cur, ext_rows)?;
    let z_r2 = take_frs(&mut cur, ext_cols)?;
    let row_root: [u8; 32] = take_bytes(&mut cur, 32)?.try_into().unwrap();
    let col_root: [u8; 32] = take_bytes(&mut cur, 32)?.try_into().unwrap();
    let bitfield = take_bytes(&mut cur, (ext_cols + 7) / 8)?.to_vec();
    let present: usize = (0..ext_cols)
        .filter(|&c| bitfield[c / 8] & (1 << (c % 8)) != 0)
        .count();
    if cur.len() != present * ext_rows * 32 {
        return Err("cell section size mismatch".into());
    }
    let mut columns: Vec<Option<Vec<Fr>>> = vec![None; ext_cols];
    for c in 0..ext_cols {
        if bitfield[c / 8] & (1 << (c % 8)) != 0 {
            columns[c] = Some(take_frs(&mut cur, ext_rows)?);
        }
    }
    let params = ZodaParams::new(m, k);
    let public = ZodaPublic {
        g_r,
        g_r2,
        z_r,
        z_r2,
        row_root,
        col_root,
    };
    let mut custody = GridCustody::new(grid_id, slot, params, public);
    for (c, col) in columns.into_iter().enumerate() {
        if let Some(col) = col {
            custody.put_column_verified(c, col)?;
        }
    }
    Ok(custody)
}
