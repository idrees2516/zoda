//! # zoda-edas — EIP-7594 / PeerDAS erasure-coded data availability
//!
//! The cell layer of EIP-7594 (Fulu): extend blobs to 8192 evaluations,
//! slice into 128 cells of 64 field elements, produce a KZG proof per
//! cell via the **FK20 algorithm** (O(n log n) instead of 128 separate
//! quotient MSMs), batch-verify sampled cells with a single pairing, and
//! reconstruct missing cells from any 64 of 128.
//!
//! Constants follow the consensus spec:
//! `FIELD_ELEMENTS_PER_BLOB = 4096`, `FIELD_ELEMENTS_PER_CELL = 64`,
//! `CELLS_PER_EXT_BLOB = 128`, `NUMBER_OF_COLUMNS = 128`.

use zoda_kzg::msm::G1J;

pub mod pipeline;
use zoda_kzg::srs::{Setup, FIELD_ELEMENTS_PER_BLOB, FIELD_ELEMENTS_PER_EXT_BLOB};
use zoda_math::{batch_invert, FftDomain, Fr, PrimeField};
use zoda_math::ntt::bit_reversal_permutation_typed;

pub const FIELD_ELEMENTS_PER_CELL: usize = 64;
pub const CELLS_PER_EXT_BLOB: usize = FIELD_ELEMENTS_PER_EXT_BLOB / FIELD_ELEMENTS_PER_CELL;
pub const CELLS_PER_BLOB: usize = FIELD_ELEMENTS_PER_BLOB / FIELD_ELEMENTS_PER_CELL;
pub const NUMBER_OF_COLUMNS: usize = CELLS_PER_EXT_BLOB;
pub const SAMPLES_PER_SLOT: u64 = 8;
pub const NUMBER_OF_CUSTODY_GROUPS: u64 = 128;
pub const CUSTODY_REQUIREMENT: u64 = 4;
pub const BYTES_PER_CELL: usize = FIELD_ELEMENTS_PER_CELL * 32;

/// A cell: 64 field elements, big-endian 32 bytes each.
pub type Cell = [u8; BYTES_PER_CELL];

fn fr_to_be32(x: Fr) -> [u8; 32] {
    let le = x.to_le_bytes();
    let mut be = [0u8; 32];
    for (i, b) in le.iter().rev().enumerate() {
        be[i] = *b;
    }
    be
}

fn fr_from_be32(b: &[u8]) -> Option<Fr> {
    let mut le = [0u8; 32];
    for (i, x) in b.iter().rev().enumerate() {
        le[i] = *x;
    }
    let fr = Fr::from_le_bytes_mod_order(&le);
    // canonical check
    if fr_to_be32(fr) == *b {
        Some(fr)
    } else {
        None
    }
}

fn cells_to_fr(cells: &[Cell]) -> Result<Vec<Vec<Fr>>, String> {
    cells
        .iter()
        .map(|c| {
            (0..FIELD_ELEMENTS_PER_CELL)
                .map(|i| {
                    fr_from_be32(&c[i * 32..(i + 1) * 32])
                        .ok_or_else(|| "non-canonical cell element".to_string())
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect()
}

fn frs_to_cell(frs: &[Fr]) -> Cell {
    let mut cell = [0u8; BYTES_PER_CELL];
    for (i, f) in frs.iter().enumerate() {
        cell[i * 32..(i + 1) * 32].copy_from_slice(&fr_to_be32(*f));
    }
    cell
}

/// Convert a blob (evaluation form over the BRP'd 4096-domain) to
/// monomial coefficients — the spec's `poly_lagrange_to_monomial`.
pub fn blob_to_monomial(blob: &[u8], setup: &Setup) -> Result<Vec<Fr>, String> {
    let evals = zoda_kzg::eip4844::blob_to_polynomial(blob)?;
    let mut brp = evals;
    bit_reversal_permutation_typed(&mut brp);
    let mut coeffs = vec![Fr::zero(); FIELD_ELEMENTS_PER_BLOB];
    setup.blob_domain.ifft(&brp, &mut coeffs);
    Ok(coeffs)
}

/// The DAS extension: evaluate the degree-<4096 blob polynomial on the
/// full 8192-domain and bit-reverse. Returns the 128 cells.
pub fn compute_cells(blob: &[u8], setup: &Setup) -> Result<Vec<Cell>, String> {
    let poly = blob_to_monomial(blob, setup)?;
    // pad to 8192 and forward FFT over the extended domain
    let mut coeffs = vec![Fr::zero(); FIELD_ELEMENTS_PER_EXT_BLOB];
    coeffs[..poly.len()].copy_from_slice(&poly);
    let mut evals = vec![Fr::zero(); FIELD_ELEMENTS_PER_EXT_BLOB];
    setup.ext_domain.fft_in_place(&mut coeffs);
    evals.copy_from_slice(&coeffs);
    bit_reversal_permutation_typed(&mut evals);
    // slice into cells
    Ok((0..CELLS_PER_EXT_BLOB)
        .map(|i| frs_to_cell(&evals[i * FIELD_ELEMENTS_PER_CELL..(i + 1) * FIELD_ELEMENTS_PER_CELL]))
        .collect())
}

/// Compute all cells and their FK20 cell proofs.
pub fn compute_cells_and_kzg_proofs(
    blob: &[u8],
    setup: &Setup,
) -> Result<(Vec<Cell>, Vec<[u8; 48]>), String> {
    let cells = compute_cells(blob, setup)?;
    let poly = blob_to_monomial(blob, setup)?;
    let proofs = fk20_cell_proofs(&poly, setup);
    let proofs = bit_reverse_proofs(proofs);
    Ok((
        cells,
        proofs
            .iter()
            .map(|p| {
                let v = p.to_compressed();
                let mut out = [0u8; 48];
                out.copy_from_slice(&v);
                out
            })
            .collect(),
    ))
}

fn bit_reverse_proofs(mut v: Vec<zoda_bls::g1::G1Affine>) -> Vec<zoda_bls::g1::G1Affine> {
    bit_reversal_permutation_typed(&mut v);
    v
}

/// FK20 cell proofs for a degree-<4096 polynomial (in monomial form):
/// returns 128 G1 proofs in natural order (caller bit-reverses).
///
/// Phase 1: 64 circulant Toeplitz columns → FFT → MSM against the
/// precomputed setup columns. Phase 2: inverse G1-FFT (unscaled) and a
/// final forward G1-FFT of size 128.
///
/// Performance notes (vs a from-scratch implementation):
/// * the 128×64 phase-1 MSMs run against **precomputed Straus window
///   tables** (6-bit windows over the fixed setup columns, built once
///   per setup and cached) — each row costs 43 windowed rounds of
///   mixed additions + 255 shared doublings, with no bucket
///   decomposition per call;
/// * the rows are spread across the available cores;
/// * both G1-FFTs run through the Jacobian engine in `zoda-kzg::g1fft`
///   (unit-twiddle and identity skips, batch normalization).
pub fn fk20_cell_proofs(poly: &[Fr], setup: &Setup) -> Vec<zoda_bls::g1::G1Affine> {
    let n = FIELD_ELEMENTS_PER_BLOB;
    debug_assert!(poly.len() <= n);
    let r = CELLS_PER_BLOB; // 64
    let l = FIELD_ELEMENTS_PER_CELL; // 64
    let circ = 2 * r; // 128
    let domain_size = circ;

    // circulant coefficients (length 128) for each of the 64 offsets
    let mut w = vec![vec![Fr::zero(); domain_size]; l];
    for (i, wi) in w.iter_mut().enumerate().take(l) {
        // out[0] = p[d - i]; out[2r - j] = p[d - j*l - i] for j = 1..r-1
        let d = n - 1;
        wi[0] = poly.get(d - i).copied().unwrap_or(Fr::zero());
        for j in 1..r - 1 {
            let idx = d - j * l - i;
            if idx < poly.len() {
                wi[2 * r - j] = poly[idx];
            }
        }
    }
    // forward FFTs of the circulant rows
    let dom = FftDomain::<Fr>::new(domain_size);
    let mut w_fft = vec![vec![Fr::zero(); domain_size]; l];
    for (i, wi) in w.iter().enumerate().take(l) {
        dom.fft(&wi, &mut w_fft[i]);
    }
    // transpose into per-row coefficient lists, prescaled by 1/128
    let inv_n = Fr::from_u64(domain_size as u64).invert().unwrap();
    let mut coeffs = vec![vec![Fr::zero(); l]; domain_size];
    for j in 0..domain_size {
        for i in 0..l {
            coeffs[j][i] = w_fft[i][j] * inv_n;
        }
    }
    // setup columns and their Straus tables (cached per setup)
    let columns = setup_columns(setup);
    let tables = fk20_tables(setup);

    // phase-1 row MSMs, spread over the available cores
    let mut u: Vec<G1J> = vec![G1J::identity(); domain_size];
    {
        let nthreads = std::thread::available_parallelism()
            .map(|c| c.get())
            .unwrap_or(1)
            .max(1);
        let per = (domain_size + nthreads - 1) / nthreads;
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for t in 0..nthreads {
                let start = t * per;
                if start >= domain_size {
                    break;
                }
                let end = (start + per).min(domain_size);
                let cols = &columns[start..end];
                let tabs = &tables[start..end];
                let cfs = &coeffs[start..end];
                handles.push(scope.spawn(move || {
                    (0..end - start)
                        .map(|k| strauss_row_msm(&cols[k], &tabs[k], &cfs[k]))
                        .collect::<Vec<G1J>>()
                }));
            }
            let mut off = 0;
            for h in handles {
                let part = h.join().expect("fk20 row msm panicked");
                u[off..off + part.len()].copy_from_slice(&part);
                off += part.len();
            }
        });
    }

    // inverse G1-FFT (unscaled — the 1/128 prescaling is absorbed above)
    let inv_roots: Vec<Fr> = dom.roots_of_unity()[..domain_size]
        .iter()
        .map(|x| x.invert().unwrap())
        .collect();
    let mut u2 = zoda_kzg::g1fft::g1_fft_jac(&u, &inv_roots);
    // zero the second half (v has degree r - 1)
    for x in u2.iter_mut().take(domain_size).skip(r) {
        *x = G1J::identity();
    }
    // final forward G1-FFT of size 128
    let roots: Vec<Fr> = dom.roots_of_unity()[..domain_size].to_vec();
    let out = zoda_kzg::g1fft::g1_fft_jac(&u2, &roots);
    G1J::batch_normalize(&out)
}

/// Table-driven Straus MSM for one FK20 row: Σ coeffs[i]·points[i]
/// against precomputed 6-bit multiple tables
/// (`tables[i][d-1] = d·points[i]` for d in 1..=63). All window digits
/// are extracted once per call.
fn strauss_row_msm(
    points: &[zoda_bls::g1::G1Affine],
    tables: &[Vec<zoda_bls::g1::G1Affine>],
    coeffs: &[Fr],
) -> G1J {
    const W: usize = 6;
    const NWIN: usize = (256 + W - 1) / W; // 43
    let n = points.len();
    // window digits, window-major
    let mut digits = vec![0u8; NWIN * n];
    for (i, sc) in coeffs.iter().enumerate() {
        let rr = sc.to_repr();
        for w in 0..NWIN {
            let lo = w * W;
            let mut d = 0usize;
            for b in 0..W {
                let bit = lo + b;
                if bit >= 256 {
                    break;
                }
                if (rr[bit / 64] >> (bit % 64)) & 1 == 1 {
                    d |= 1 << b;
                }
            }
            digits[w * n + i] = d as u8;
        }
    }
    let mut acc = G1J::identity();
    for w in (0..NWIN).rev() {
        if !acc.is_identity() {
            for _ in 0..W {
                acc = acc.double();
            }
        }
        let dm = &digits[w * n..(w + 1) * n];
        for (i, &d) in dm.iter().enumerate() {
            if d != 0 && !points[i].infinity {
                acc = acc.add_mixed(&tables[i][(d - 1) as usize]);
            }
        }
    }
    acc
}

/// Process-wide FK20 caches, keyed by the Setup's G1 monomial first
/// point (cheap identity for the toy setups; real setups are loaded
/// once per process anyway). **Shared across threads** (Arc), so worker
/// threads — batched recovery, threaded pipelines — never rebuild the
/// ~49 MB Straus table that a thread_local cache would force on every
/// fresh thread. The build runs under the cache Mutex: concurrent
/// first-callers serialise on it once and then all observe the entry
/// (the build itself is single-threaded, so nothing is lost).
fn fk20_column_cache() -> &'static std::sync::Mutex<
    std::collections::HashMap<u64, std::sync::Arc<Vec<Vec<zoda_bls::g1::G1Affine>>>>,
> {
    static C: std::sync::OnceLock<
        std::sync::Mutex<
            std::collections::HashMap<
                u64,
                std::sync::Arc<Vec<Vec<zoda_bls::g1::G1Affine>>>,
            >,
        >,
    > = std::sync::OnceLock::new();
    C.get_or_init(|| std::sync::Mutex::new(std::collections::HashMap::new()))
}

/// Straus window tables for the FK20 column points:
/// `tables[row][offset][d-1] = d·columns[row][offset]` for d in 1..=63.
/// 63 affine multiples per point — ~49 MB for the 128×64 matrix,
/// built once per setup (a one-time cost mirroring c-kzg's
/// fixed-base precompute option).
fn fk20_table_cache() -> &'static std::sync::Mutex<
    std::collections::HashMap<u64, std::sync::Arc<Vec<Vec<Vec<zoda_bls::g1::G1Affine>>>>>,
> {
    static C: std::sync::OnceLock<
        std::sync::Mutex<
            std::collections::HashMap<
                u64,
                std::sync::Arc<Vec<Vec<Vec<zoda_bls::g1::G1Affine>>>>,
            >,
        >,
    > = std::sync::OnceLock::new();
    C.get_or_init(|| std::sync::Mutex::new(std::collections::HashMap::new()))
}

/// Cache key: the setup's process-unique construction id (NOT a content
/// probe — the monomial points of distinct setups share `[τ⁰] = g1`).
fn setup_key(setup: &Setup) -> u64 {
    setup.id
}

/// The FK20 setup columns: x_ext_fft_columns[row][offset] for row in
/// 0..128, offset in 0..64 — computed from the monomial G1 setup
/// (cached across calls).
fn setup_columns(setup: &Setup) -> Vec<Vec<zoda_bls::g1::G1Affine>> {
    setup_columns_cached(setup).to_vec()
}

fn setup_columns_cached(setup: &Setup) -> std::sync::Arc<Vec<Vec<zoda_bls::g1::G1Affine>>> {
    let key = setup_key(setup);
    let cache = fk20_column_cache();
    let mut c = cache
        .lock()
        .expect("fk20 column cache poisoned");
    if let Some(cols) = c.get(&key) {
        return cols.clone();
    }
    let cols = std::sync::Arc::new(setup_columns_uncached(setup));
    c.insert(key, cols.clone());
    cols
}

fn setup_columns_uncached(setup: &Setup) -> Vec<Vec<zoda_bls::g1::G1Affine>> {
    let n = FIELD_ELEMENTS_PER_BLOB;
    let r = CELLS_PER_BLOB;
    let l = FIELD_ELEMENTS_PER_CELL;
    let circ = 2 * r;
    let mut out = vec![vec![zoda_bls::g1::G1Affine::identity(); l]; circ];
    // x[i] = g1_monomial[4096 - 64 - 1 - offset - i*64] for i < 63
    let mut x = vec![zoda_bls::g1::G1Affine::identity(); r];
    for offset in 0..l {
        for i in 0..(r - 1) {
            let start = n - l - 1 - offset;
            let j = start - i * l;
            x[i] = setup.g1_monomial[j];
        }
        // points = FFT_128(x padded to 128)
        let mut padded: Vec<G1J> = vec![G1J::identity(); circ];
        for (i, xi) in x.iter().enumerate() {
            padded[i] = G1J::from_affine(xi);
        }
        let dom = FftDomain::<Fr>::new(circ);
        let roots: Vec<Fr> = dom.roots_of_unity()[..circ].to_vec();
        let pts = zoda_kzg::g1fft::g1_fft_jac(&padded, &roots);
        let affine = G1J::batch_normalize(&pts);
        for (row, p) in affine.iter().enumerate() {
            out[row][offset] = *p;
        }
    }
    out
}

/// Build (and cache) the Straus 6-bit window tables for the FK20
/// columns. All 8,192 × 63 multiples are produced with one inversion
/// via a single batch normalization.
fn fk20_tables(setup: &Setup) -> std::sync::Arc<Vec<Vec<Vec<zoda_bls::g1::G1Affine>>>> {
    let key = setup_key(setup);
    let cache = fk20_table_cache();
    let mut c = cache.lock().expect("fk20 table cache poisoned");
    if let Some(t) = c.get(&key) {
        return t.clone();
    }
    let columns = setup_columns(setup);
    const PER: usize = 63; // multiples 1..=63 per point
    let total: usize = columns.len() * FIELD_ELEMENTS_PER_CELL * PER;
    let mut flat: Vec<G1J> = Vec::with_capacity(total);
    for row in &columns {
        for p in row {
            let base = G1J::from_affine(p);
            let mut cur = base;
            flat.push(base); // 1·P
            for _ in 1..PER {
                cur = cur.add(&base);
                flat.push(cur);
            }
        }
    }
    let norm = G1J::batch_normalize(&flat);
    let mut tables: Vec<Vec<Vec<zoda_bls::g1::G1Affine>>> =
        Vec::with_capacity(columns.len());
    let mut it = norm.into_iter();
    for row in &columns {
        let mut row_t = Vec::with_capacity(row.len());
        for _ in row {
            let tab: Vec<zoda_bls::g1::G1Affine> = (&mut it).take(PER).collect();
            row_t.push(tab);
        }
        tables.push(row_t);
    }
    let t = std::sync::Arc::new(tables);
    c.insert(key, t.clone());
    t
}

/// Forward FFT over G1 (homogeneous in/out) — compatibility wrapper
/// around the Jacobian engine in `zoda-kzg::g1fft`.
pub fn g1_fft(
    input: &[zoda_bls::g1::G1Projective],
    roots: &[Fr],
) -> Vec<zoda_bls::g1::G1Projective> {
    let jacs: Vec<G1J> = input.iter().map(G1J::from_homogeneous).collect();
    zoda_kzg::g1fft::g1_fft_jac(&jacs, roots)
        .into_iter()
        .map(|p| p.to_homogeneous())
        .collect()
}

/// `verify_cell_kzg_proof_batch` — the EIP-7594 batch verification.
pub fn verify_cell_kzg_proof_batch(
    commitments: &[&[u8]],
    cell_indices: &[u64],
    cells: &[Cell],
    proofs: &[&[u8]],
    setup: &Setup,
) -> Result<bool, String> {
    use zoda_bls::g1::G1Affine;
    use zoda_bls::pairing::{multi_miller_loop, final_exponentiation, G2Prepared};

    let n = cells.len();
    if commitments.len() != n || cell_indices.len() != n || proofs.len() != n {
        return Err("input length mismatch".to_string());
    }
    for &ci in cell_indices {
        if ci as usize >= CELLS_PER_EXT_BLOB {
            return Err("cell index out of range".to_string());
        }
    }
    if n == 0 {
        return Ok(true);
    }

    // Deduplicate commitments
    let mut unique: Vec<&[u8]> = Vec::with_capacity(n);
    let mut commitment_index = vec![0usize; n];
    for (i, c) in commitments.iter().enumerate() {
        match unique.iter().position(|u| u == c) {
            Some(pos) => commitment_index[i] = pos,
            None => {
                unique.push(c);
                commitment_index[i] = unique.len() - 1;
            }
        }
    }

    // Fiat–Shamir challenge (RCKZGCBATCH__V1_)
    let mut data = Vec::new();
    data.extend_from_slice(b"RCKZGCBATCH__V1_");
    data.extend_from_slice(&(FIELD_ELEMENTS_PER_BLOB as u64).to_be_bytes());
    data.extend_from_slice(&(FIELD_ELEMENTS_PER_CELL as u64).to_be_bytes());
    data.extend_from_slice(&(unique.len() as u64).to_be_bytes());
    data.extend_from_slice(&(n as u64).to_be_bytes());
    for c in &unique {
        data.extend_from_slice(c);
    }
    for i in 0..n {
        data.extend_from_slice(&(commitment_index[i] as u64).to_be_bytes());
        data.extend_from_slice(&cell_indices[i].to_be_bytes());
        data.extend_from_slice(&cells[i]);
        data.extend_from_slice(proofs[i]);
    }
    let r = zoda_kzg::eip4844::hash_to_bls_field(&data);
    let r_powers = {
        let mut v = Vec::with_capacity(n);
        let mut cur = Fr::ONE;
        for _ in 0..n {
            v.push(cur);
            cur = cur * r;
        }
        v
    };

    // Parse points (subgroup checks are aggregated below: one MSM and a
    // single [r]S multiplication replace n individual [r]P checks)
    let proofs_g1: Vec<G1Affine> = proofs
        .iter()
        .map(|p| zoda_kzg::eip4844::parse_g1_canonical(p))
        .collect::<Result<_, _>>()?;
    let commitments_g1: Vec<G1Affine> = unique
        .iter()
        .map(|c| zoda_kzg::eip4844::parse_g1_canonical(c))
        .collect::<Result<_, _>>()?;
    {
        let mut all = commitments_g1.clone();
        all.extend_from_slice(&proofs_g1);
        if !zoda_kzg::eip4844::batch_subgroup_check_g1(&all) {
            return Err("G1 point not in subgroup".to_string());
        }
    }

    // proof_lincomb = Σ r^i π_i
    let proof_lincomb = zoda_kzg::msm::g1_lincomb(&proofs_g1, &r_powers);

    // Σ r^i z_i π_i  (the coset-weighted proof sum)
    let weighted: Vec<Fr> = cell_indices
        .iter()
        .zip(r_powers.iter())
        .map(|(&ci, rp)| {
            // h_k^n = (ζ^{brp(ci)}·n-th root...) — computed via the
            // roots-of-unity table
            let brp = reverse_bits_limited(CELLS_PER_EXT_BLOB, ci as usize);
            let idx = (brp * FIELD_ELEMENTS_PER_CELL) % FIELD_ELEMENTS_PER_EXT_BLOB;
            let h = setup.roots_of_unity[idx];
            h * *rp
        })
        .collect();
    let weighted_proof_sum = zoda_kzg::msm::g1_lincomb(&proofs_g1, &weighted);

    // commitment weights per unique commitment
    let mut commit_weights = vec![Fr::zero(); unique.len()];
    for i in 0..n {
        commit_weights[commitment_index[i]] = commit_weights[commitment_index[i]] + r_powers[i];
    }
    let commit_sum = zoda_kzg::msm::g1_lincomb(&commitments_g1, &commit_weights);

    // aggregated interpolation polynomial commitment:
    // Σ over used columns: r-weighted cells → interpolate per column →
    // commit via monomial MSM
    let mut column_cells = vec![vec![Fr::zero(); FIELD_ELEMENTS_PER_CELL]; CELLS_PER_EXT_BLOB];
    let mut column_used = vec![false; CELLS_PER_EXT_BLOB];
    let cell_frs = cells_to_fr(cells)?;
    for i in 0..n {
        let c = cell_indices[i] as usize;
        column_used[c] = true;
        for j in 0..FIELD_ELEMENTS_PER_CELL {
            column_cells[c][j] = column_cells[c][j] + cell_frs[i][j] * r_powers[i];
        }
    }
    let mut interp = vec![Fr::zero(); FIELD_ELEMENTS_PER_CELL];
    let dom = FftDomain::<Fr>::new(FIELD_ELEMENTS_PER_CELL);
    for c in 0..CELLS_PER_EXT_BLOB {
        if !column_used[c] {
            continue;
        }
        // bit-reverse the column, inverse FFT over the 64-domain
        let mut col = column_cells[c].clone();
        bit_reversal_permutation_typed(&mut col);
        let mut coeffs = vec![Fr::zero(); FIELD_ELEMENTS_PER_CELL];
        dom.ifft(&col, &mut coeffs);
        // shift by h_k^{-1}
        let brp = reverse_bits_limited(CELLS_PER_EXT_BLOB, c);
        let idx = (FIELD_ELEMENTS_PER_EXT_BLOB - brp) % FIELD_ELEMENTS_PER_EXT_BLOB;
        let inv_h = setup.roots_of_unity[idx];
        let mut sh = Fr::ONE;
        for k in 0..FIELD_ELEMENTS_PER_CELL {
            coeffs[k] = coeffs[k] * sh;
            sh = sh * inv_h;
        }
        for k in 0..FIELD_ELEMENTS_PER_CELL {
            interp[k] = interp[k] + coeffs[k];
        }
    }
    let interp_commit = setup.commit_coeffs(&interp);

    // final_g1 = commit_sum - interp_commit + weighted_proof_sum
    let final_g1 = commit_sum
        .to_projective()
        .add(&interp_commit.to_projective().neg())
        .add(&weighted_proof_sum.to_projective())
        .to_affine();

    // pairing check: e(final_g1, g2_gen) * e(proof_lincomb, -[s^64]_2) == 1
    let power_of_s = setup.g2_monomial[FIELD_ELEMENTS_PER_CELL];
    let neg_s = G2Prepared::from((-power_of_s).to_projective().to_affine());
    let g2_gen = G2Prepared::from(setup.g2_monomial[0]);
    let ml = multi_miller_loop(&[(&final_g1, &g2_gen), (&proof_lincomb, &neg_s)]);
    match final_exponentiation(&ml) {
        Some(f) => Ok(f.is_one()),
        None => Ok(false),
    }
}

fn reverse_bits_limited(order: usize, n: usize) -> usize {
    zoda_math::ntt::reverse_bits_limited(order, n)
}

/// `recover_cells_and_kzg_proofs`: given at least 64 of the 128 cells
/// (identified by index), reconstruct everything.
pub fn recover_cells_and_kzg_proofs(
    cell_indices: &[u64],
    cells: &[Cell],
    setup: &Setup,
) -> Result<(Vec<Cell>, Vec<[u8; 48]>), String> {
    if cell_indices.len() != cells.len() {
        return Err("index/cell length mismatch".to_string());
    }
    let n = cell_indices.len();
    if n < CELLS_PER_BLOB {
        return Err(format!("need at least {} cells, got {}", CELLS_PER_BLOB, n));
    }
    if n > CELLS_PER_EXT_BLOB {
        return Err("too many cells".to_string());
    }
    for (i, &ci) in cell_indices.iter().enumerate() {
        if ci as usize >= CELLS_PER_EXT_BLOB {
            return Err("cell index out of range".to_string());
        }
        if i > 0 && cell_indices[i] <= cell_indices[i - 1] {
            return Err("cell indices must be strictly ascending".to_string());
        }
    }
    if n == CELLS_PER_EXT_BLOB {
        let poly = cells_to_fr(cells)?;
        // treat as evaluations, produce proofs
        let coeffs = evals_to_coeffs(&poly.concat(), setup);
        let proofs = fk20_cell_proofs(&coeffs, setup);
        let proofs = bit_reverse_proofs(proofs);
        return Ok((
            cells.to_vec(),
            proofs
                .iter()
                .map(|p| {
                    let v = p.to_compressed();
                    let mut o = [0u8; 48];
                    o.copy_from_slice(&v);
                    o
                })
                .collect(),
        ));
    }

    // full-codeword erasure recovery over the 8192-domain
    let total = FIELD_ELEMENTS_PER_EXT_BLOB;
    let mut known = vec![Fr::zero(); total];
    let mut have = vec![false; total];
    for (k, &ci) in cell_indices.iter().enumerate() {
        let frs = {
            let mut v = vec![Fr::zero(); FIELD_ELEMENTS_PER_CELL];
            for j in 0..FIELD_ELEMENTS_PER_CELL {
                v[j] = fr_from_be32(&cells[k][j * 32..(j + 1) * 32])
                    .ok_or("non-canonical cell element")?;
            }
            v
        };
        for (j, f) in frs.iter().enumerate() {
            // cell ci covers BRP-ordered positions; work in BRP space
            // then un-BRP at the end
            let pos = ci as usize * FIELD_ELEMENTS_PER_CELL + j;
            known[pos] = *f;
            have[pos] = true;
        }
    }

    // Vanishing polynomial over the missing cells' COSETS:
    // Z(X) = short_Z(X^64) where short_Z vanishes at {h_c^64} for the
    // missing cells c (h_c = ω^{brp7(c)}). The X^64 spread makes Z vanish
    // on every point of each missing coset.
    let mut missing_cells: Vec<usize> = Vec::new();
    for c in 0..CELLS_PER_EXT_BLOB {
        if !have[c * FIELD_ELEMENTS_PER_CELL] {
            missing_cells.push(c);
        }
    }
    let roots: Vec<Fr> = missing_cells
        .iter()
        .map(|&c| {
            let brp = reverse_bits_limited(CELLS_PER_EXT_BLOB, c);
            setup.roots_of_unity[(brp * FIELD_ELEMENTS_PER_CELL) % total]
        })
        .collect();
    let short_vanish = zoda_math::poly::vanishing_poly(&roots);
    let mut vanish = vec![Fr::zero(); total];
    for (i, &v) in short_vanish.iter().enumerate() {
        vanish[i * FIELD_ELEMENTS_PER_CELL] = v;
    }

    // E·Z in evaluation form (BRP space), convert to coefficients
    let mut ez = vec![Fr::zero(); total];
    let mut z_evals = vec![Fr::zero(); total];
    setup.ext_domain.fft(&vanish, &mut z_evals);
    bit_reversal_permutation_typed(&mut z_evals);
    for i in 0..total {
        ez[i] = known[i] * z_evals[i];
    }
    bit_reversal_permutation_typed(&mut ez);
    let mut ez_coeffs = vec![Fr::zero(); total];
    setup.ext_domain.ifft(&ez, &mut ez_coeffs);

    // divide by Z on a coset
    let shift = Fr::from_u64(7);
    let mut ez_coset = vec![Fr::zero(); total];
    let mut z_coset = vec![Fr::zero(); total];
    setup.ext_domain.coset_fft(&ez_coeffs, &mut ez_coset, shift);
    setup.ext_domain.coset_fft(&vanish, &mut z_coset, shift);
    let mut q_coset = ez_coset.clone();
    batch_invert(&mut z_coset);
    for i in 0..total {
        q_coset[i] = q_coset[i] * z_coset[i];
    }
    let mut q_coeffs = vec![Fr::zero(); total];
    setup.ext_domain.coset_ifft(&q_coset, &mut q_coeffs, shift);

    // evaluate to get all 8192 points (BRP order)
    let mut evals = vec![Fr::zero(); total];
    setup.ext_domain.fft(&q_coeffs, &mut evals);
    bit_reversal_permutation_typed(&mut evals);

    // build cells
    let out_cells: Vec<Cell> = (0..CELLS_PER_EXT_BLOB)
        .map(|i| frs_to_cell(&evals[i * FIELD_ELEMENTS_PER_CELL..(i + 1) * FIELD_ELEMENTS_PER_CELL]))
        .collect();

    // proofs: monomial form then FK20
    // the polynomial in monomial form: coeffs from the 4096-low part of q
    let poly_coeffs: Vec<Fr> = q_coeffs[..FIELD_ELEMENTS_PER_BLOB].to_vec();
    let mut proofs = fk20_cell_proofs(&poly_coeffs, setup);
    bit_reversal_permutation_typed(&mut proofs);
    let proof_bytes: Vec<[u8; 48]> = proofs
        .iter()
        .map(|p| {
            let v = p.to_compressed();
            let mut o = [0u8; 48];
            o.copy_from_slice(&v);
            o
        })
        .collect();
    Ok((out_cells, proof_bytes))
}

/// Convert evaluation form (4096 evals, natural order) to coefficients.
fn evals_to_coeffs(evals: &[Fr], setup: &Setup) -> Vec<Fr> {
    let mut brp = evals.to_vec();
    bit_reversal_permutation_typed(&mut brp);
    let mut coeffs = vec![Fr::zero(); FIELD_ELEMENTS_PER_BLOB];
    setup.blob_domain.ifft(&brp, &mut coeffs);
    coeffs
}

/// `get_custody_groups` from the Fulu spec: deterministic custody group
/// selection from a node id.
pub fn get_custody_groups(node_id: [u8; 32], custody_group_count: u64) -> Vec<u64> {
    assert!(custody_group_count <= NUMBER_OF_CUSTODY_GROUPS);
    if custody_group_count == NUMBER_OF_CUSTODY_GROUPS {
        return (0..NUMBER_OF_CUSTODY_GROUPS).collect();
    }
    let mut current = node_id;
    let mut groups: Vec<u64> = Vec::new();
    while groups.len() < custody_group_count as usize {
        let h = zoda_math::sha256::sha256(&current);
        let v = u64::from_le_bytes(h[..8].try_into().unwrap());
        let g = v % NUMBER_OF_CUSTODY_GROUPS;
        if !groups.contains(&g) {
            groups.push(g);
        }
        // increment the 32-byte counter
        for b in current.iter_mut().rev() {
            if *b == 255 {
                *b = 0;
            } else {
                *b += 1;
                break;
            }
        }
    }
    groups.sort_unstable();
    groups
}

/// `compute_columns_for_custody_group`: columns per custody group.
pub fn compute_columns_for_custody_group(custody_group: u64) -> Vec<u64> {
    assert!(custody_group < NUMBER_OF_CUSTODY_GROUPS);
    let per_group = NUMBER_OF_COLUMNS as u64 / NUMBER_OF_CUSTODY_GROUPS;
    (0..per_group)
        .map(|i| NUMBER_OF_CUSTODY_GROUPS * i + custody_group)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_blob(seed: [u8; 32]) -> Vec<u8> {
        let mut rng = zoda_math::ZodaRng::from_seed(seed);
        let mut blob = vec![0u8; 131072];
        rng.next_bytes(&mut blob);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f;
        }
        blob
    }

    #[test]
    fn cells_roundtrip_and_recover() {
        let setup = Setup::from_seed_for_testing(*b"edas-test-setup-tau-000000000000");
        let blob = test_blob(*b"edas-cell-test-seed-000000000000");
        let (cells, proofs) = compute_cells_and_kzg_proofs(&blob, &setup).unwrap();
        assert_eq!(cells.len(), 128);
        assert_eq!(proofs.len(), 128);

        // commitment to the blob
        let comm = zoda_kzg::eip4844::blob_to_kzg_commitment(&blob, &setup).unwrap();

        // verify all 128 cells in one batch
        let idx: Vec<u64> = (0..128).collect();
        let comm_refs: Vec<&[u8]> = idx.iter().map(|_| comm.as_slice()).collect();
        let proof_refs: Vec<&[u8]> = proofs.iter().map(|p| p.as_slice()).collect();
        assert!(verify_cell_kzg_proof_batch(&comm_refs, &idx, &cells, &proof_refs, &setup).unwrap());

        // verify a 16-cell sample
        let sample_idx: Vec<u64> = (0..16).map(|i| i as u64 * 7 % 128).collect();
        let sample_cells: Vec<Cell> = sample_idx.iter().map(|&i| cells[i as usize]).collect();
        let sample_proofs: Vec<&[u8]> =
            sample_idx.iter().map(|&i| proofs[i as usize].as_slice()).collect();
        let sample_comms: Vec<&[u8]> = sample_idx.iter().map(|_| comm.as_slice()).collect();
        assert!(verify_cell_kzg_proof_batch(
            &sample_comms,
            &sample_idx,
            &sample_cells,
            &sample_proofs,
            &setup
        )
        .unwrap());

        // tampered cell fails
        let mut bad = sample_cells.clone();
        bad[0][0] ^= 0x01;
        assert!(!verify_cell_kzg_proof_batch(
            &sample_comms,
            &sample_idx,
            &bad,
            &sample_proofs,
            &setup
        )
        .unwrap());

        // recover from 64 cells
        let keep_idx: Vec<u64> = (0..64).map(|i| i as u64).collect();
        let keep_cells: Vec<Cell> = keep_idx.iter().map(|&i| cells[i as usize]).collect();
        let (rec_cells, rec_proofs) =
            recover_cells_and_kzg_proofs(&keep_idx, &keep_cells, &setup).unwrap();
        assert_eq!(rec_cells.len(), 128);
        assert_eq!(rec_cells, cells, "recovered cells differ");
        let rec_proof_refs: Vec<&[u8]> = rec_proofs.iter().map(|p| p.as_slice()).collect();
        let all_idx: Vec<u64> = (0..128).map(|i| i as u64).collect();
        assert!(verify_cell_kzg_proof_batch(&comm_refs, &all_idx, &rec_cells, &rec_proof_refs, &setup).unwrap());
    }

    #[test]
    fn custody_groups_deterministic() {
        let a = get_custody_groups([7u8; 32], 4);
        let b = get_custody_groups([7u8; 32], 4);
        assert_eq!(a, b);
        assert_eq!(a.len(), 4);
        let c = get_custody_groups([8u8; 32], 4);
        assert_ne!(a, c);
        // groups extend monotonically with count
        let d = get_custody_groups([7u8; 32], 5);
        assert!(a.iter().all(|g| d.contains(g)));
    }
}

/// Official consensus-spec-tests harness for the EIP-7594 (Fulu) cell
/// suites: `compute_cells`, `compute_cells_and_kzg_proofs`,
/// `verify_cell_kzg_proof_batch`, `recover_cells_and_kzg_proofs`,
/// replayed against the mainnet trusted setup. Skipped unless the
/// vector archive is present (see `zoda-kzg::spectest`).
#[cfg(test)]
mod spec_vectors {
    use super::*;
    use zoda_kzg::spectest::{for_each_case, load_mainnet_setup, vectors_dir, Sv};

    fn setup() -> Option<zoda_kzg::srs::Setup> {
        let dir = vectors_dir()?;
        Some(load_mainnet_setup(&dir).expect("mainnet trusted setup"))
    }

    fn parse_cells(v: &Sv) -> Vec<Cell> {
        v.as_list()
            .expect("cell list")
            .iter()
            .map(|c| c.hex_array::<BYTES_PER_CELL>().expect("cell hex"))
            .collect()
    }

    /// Parse a cell list; `None` marks structurally malformed entries
    /// (wrong byte length) — such inputs must make the API error.
    fn parse_cells_lenient(v: &Sv) -> Option<Vec<Cell>> {
        let list = v.as_list().expect("cell list");
        let mut out = Vec::with_capacity(list.len());
        for c in list {
            match c.hex_array::<BYTES_PER_CELL>() {
                Some(cell) => out.push(cell),
                None => return None,
            }
        }
        Some(out)
    }

    fn parse_proofs(v: &Sv) -> Vec<[u8; 48]> {
        v.as_list()
            .expect("proof list")
            .iter()
            .map(|c| c.hex_array::<48>().expect("proof hex"))
            .collect()
    }

    fn blob_of(case: &Sv) -> Vec<u8> {
        case.get("input")
            .and_then(|i| i.get("blob"))
            .and_then(|b| b.hex())
            .expect("blob hex")
    }

    #[test]
    fn spec_vectors_compute_cells() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(&dir, "fulu", "compute_cells", |name, case| {
            let blob = blob_of(case);
            let out = case.get("output").expect("output");
            if out.is_null() {
                assert!(
                    compute_cells(&blob, &setup).is_err(),
                    "{}: expected an error",
                    name
                );
            } else {
                let expect = parse_cells(out);
                let got = compute_cells(&blob, &setup)
                    .unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                assert_eq!(got.len(), expect.len(), "{}: cell count", name);
                for (i, (g, e)) in got.iter().zip(expect.iter()).enumerate() {
                    assert_eq!(g, e, "{}: cell {}", name, i);
                }
            }
            n += 1;
        });
        assert_eq!(count, n);
        eprintln!("compute_cells: {} cases passed", n);
    }

    #[test]
    fn spec_vectors_compute_cells_and_kzg_proofs() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(
            &dir,
            "fulu",
            "compute_cells_and_kzg_proofs",
            |name, case| {
                let blob = blob_of(case);
                let out = case.get("output").expect("output");
                if out.is_null() {
                    assert!(
                        compute_cells_and_kzg_proofs(&blob, &setup).is_err(),
                        "{}: expected an error",
                        name
                    );
                } else {
                    let lst = out.as_list().expect("[cells, proofs]");
                    let expect_cells = parse_cells(&lst[0]);
                    let expect_proofs = parse_proofs(&lst[1]);
                    let (cells, proofs) = compute_cells_and_kzg_proofs(&blob, &setup)
                        .unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                    assert_eq!(cells.len(), expect_cells.len(), "{}: cells", name);
                    for (i, (g, e)) in cells.iter().zip(expect_cells.iter()).enumerate() {
                        assert_eq!(g, e, "{}: cell {}", name, i);
                    }
                    assert_eq!(proofs.len(), expect_proofs.len(), "{}: proofs", name);
                    for (i, (g, e)) in proofs.iter().zip(expect_proofs.iter()).enumerate() {
                        assert_eq!(g, e, "{}: proof {}", name, i);
                    }
                }
                n += 1;
            },
        );
        assert_eq!(count, n);
        eprintln!("compute_cells_and_kzg_proofs: {} cases passed", n);
    }

    #[test]
    fn spec_vectors_verify_cell_kzg_proof_batch() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(
            &dir,
            "fulu",
            "verify_cell_kzg_proof_batch",
            |name, case| {
                let input = case.get("input").expect("input");
                let commitments: Vec<Vec<u8>> = input
                    .get("commitments")
                    .and_then(|v| v.as_list())
                    .expect("commitments")
                    .iter()
                    .map(|c| c.hex().expect("hex"))
                    .collect();
                let cell_indices: Vec<u64> = input
                    .get("cell_indices")
                    .and_then(|v| v.as_list())
                    .expect("cell_indices")
                    .iter()
                    .map(|c| match c {
                        Sv::Int(i) => *i,
                        _ => panic!("bad cell index"),
                    })
                    .collect();
                let cells = parse_cells_lenient(input.get("cells").expect("cells"));
                let out = case.get("output").expect("output");
                match cells {
                    None => {
                        assert!(out.is_null(), "{}: malformed cells must error", name);
                    }
                    Some(cells) => {
                        let proofs: Vec<Vec<u8>> = input
                            .get("proofs")
                            .and_then(|v| v.as_list())
                            .expect("proofs")
                            .iter()
                            .map(|c| c.hex().expect("hex"))
                            .collect();
                        let comm_refs: Vec<&[u8]> =
                            commitments.iter().map(|c| c.as_slice()).collect();
                        let proof_refs: Vec<&[u8]> =
                            proofs.iter().map(|c| c.as_slice()).collect();
                        let got = verify_cell_kzg_proof_batch(
                            &comm_refs,
                            &cell_indices,
                            &cells,
                            &proof_refs,
                            &setup,
                        );
                        match out {
                            Sv::Null => assert!(got.is_err(), "{}: expected an error", name),
                            Sv::Bool(expect) => {
                                let v = got
                                    .unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                                assert_eq!(v, *expect, "{}", name);
                            }
                            other => panic!("{}: unexpected output {:?}", name, other),
                        }
                    }
                }
                n += 1;
            },
        );
        assert_eq!(count, n);
        eprintln!("verify_cell_kzg_proof_batch: {} cases passed", n);
    }

    #[test]
    fn spec_vectors_recover_cells_and_kzg_proofs() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(
            &dir,
            "fulu",
            "recover_cells_and_kzg_proofs",
            |name, case| {
                let input = case.get("input").expect("input");
                let cell_indices: Vec<u64> = input
                    .get("cell_indices")
                    .and_then(|v| v.as_list())
                    .expect("cell_indices")
                    .iter()
                    .map(|c| match c {
                        Sv::Int(i) => *i,
                        _ => panic!("bad cell index"),
                    })
                    .collect();
                let cells = parse_cells_lenient(input.get("cells").expect("cells"));
                let out = case.get("output").expect("output");
                match cells {
                    None => {
                        assert!(out.is_null(), "{}: malformed cells must error", name);
                    }
                    Some(cells) => {
                        if out.is_null() {
                            assert!(
                                recover_cells_and_kzg_proofs(&cell_indices, &cells, &setup)
                                    .is_err(),
                                "{}: expected an error",
                                name
                            );
                        } else {
                            let lst = out.as_list().expect("[cells, proofs]");
                            let expect_cells = parse_cells(&lst[0]);
                            let expect_proofs = parse_proofs(&lst[1]);
                            let (rcells, rproofs) =
                                recover_cells_and_kzg_proofs(&cell_indices, &cells, &setup)
                                    .unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                            assert_eq!(rcells.len(), expect_cells.len(), "{}: cells", name);
                            for (i, (g, e)) in rcells.iter().zip(expect_cells.iter()).enumerate() {
                                assert_eq!(g, e, "{}: cell {}", name, i);
                            }
                            assert_eq!(rproofs.len(), expect_proofs.len(), "{}: proofs", name);
                            for (i, (g, e)) in rproofs.iter().zip(expect_proofs.iter()).enumerate() {
                                assert_eq!(g, e, "{}: proof {}", name, i);
                            }
                        }
                    }
                }
                n += 1;
            },
        );
        assert_eq!(count, n);
        eprintln!("recover_cells_and_kzg_proofs: {} cases passed", n);
    }
}
