//! zoda-bench: micro + real-world benchmarks for the whole DA stack.
//!
//! Run with `cargo run --release -p zoda-bench [-- <filter>]`.

use std::time::Instant;
use zoda_bls::pairing::{pairing, G2Prepared};
use zoda_core::{Matrix, ZodaParams};
use zoda_math::{FftDomain, Fr, PrimeField, ZodaRng};

struct Timed {
    name: String,
    ns: u128,
    note: String,
}

impl Timed {
    fn print(&self) {
        println!(
            "{:<58} {:>12.3} {:>10}   {}",
            self.name,
            self.ns as f64 / 1000.0,
            "µs",
            self.note
        );
    }
}

fn bench<T>(name: &str, note: &str, iters: usize, mut f: impl FnMut() -> T) -> Timed {
    // warmup
    let _ = f();
    let t0 = Instant::now();
    let mut last = None;
    for _ in 0..iters {
        last = Some(f());
    }
    std::hint::black_box(&last);
    let dt = t0.elapsed().as_nanos() / iters as u128;
    Timed {
        name: name.to_string(),
        ns: dt,
        note: format!("{} [{} iters]", note, iters),
    }
}

fn fmt_bytes(n: f64) -> String {
    if n > 1e9 {
        format!("{:.2} GB/s", n / 1e9)
    } else if n > 1e6 {
        format!("{:.2} MB/s", n / 1e6)
    } else {
        format!("{:.2} KB/s", n / 1e3)
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let filter = args.get(1).map(|s| s.as_str()).unwrap_or("");
    let mut results: Vec<Timed> = Vec::new();
    let mut rng = ZodaRng::from_seed(*b"zoda-bench-seed-0000000000000000");

    // ---------------- field + NTT ----------------
    if filter.is_empty() || filter == "math" {
        let mut x = rng.next_fr(false);
        let y = rng.next_fr(false);
        results.push(bench("Fr mul (Montgomery, 4 limbs)", "field core", 1_000_000, || {
            x = x * y;
            x
        }));
        let xs = rng.next_fr(false);
        results.push(bench("Fr invert (Fermat)", "batch-invert backbone", 500, || {
            xs.invert()
        }));
        for log in [12u32, 13] {
            let n = 1usize << log;
            let dom = FftDomain::<Fr>::new(n);
            let mut v: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
            results.push(bench(
                &format!("NTT size {}", n),
                &format!("{}/iter", fmt_bytes(n as f64 * 32.0 * 2.0)),
                200,
                || {
                    dom.fft_in_place(&mut v);
                    v[0]
                },
            ));
        }
    }

    // ---------------- BLS + pairing ----------------
    if filter.is_empty() || filter == "bls" {
        let g1 = zoda_bls::g1::G1Affine::generator();
        let g2 = zoda_bls::g2::G2Affine::generator();
        let sk = zoda_bls::SecretKey::from_seed(b"bench-key-000000000000000000000000");
        let msg = b"benchmark message";
        results.push(bench("G1 scalar mul (windowed, 255-bit)", "sign/pk path", 50, || {
            zoda_bls::g1::G1Projective::generator().mul_fr(&rng.next_fr(false))
        }));
        results.push(bench("BLS sign (hash_to_curve G2 + mul)", "ETH2 ciphersuite", 20, || {
            sk.sign(msg)
        }));
        let pk = sk.public_key();
        let sig = sk.sign(msg);
        results.push(bench("BLS verify (2 pairings)", "ETH2 ciphersuite", 10, || {
            zoda_bls::verify(&pk, msg, &sig)
        }));
        let g2p = G2Prepared::from(g2);
        results.push(bench("pairing e(g1, g2) (miller + final exp)", "fountain chain", 10, || {
            pairing(&g1, &g2)
        }));
        let _ = g2p;
    }

    // ---------------- GLV / endomorphism / batch verification ----------------
    if filter.is_empty() || filter == "glv" {
        use zoda_bls::g1::G1Projective;
        use zoda_bls::g2::G2Projective;
        let k = rng.next_fr(false);
        let p1 = G1Projective::generator().mul_fr(&k).to_affine();
        let q2 = G2Projective::generator().mul_fr(&k).to_affine();
        results.push(bench("G1 scalar mul (plain, 255-bit)", "reference path", 50, || {
            p1.to_projective().mul_fr(&k)
        }));
        results.push(bench("G1 scalar mul (GLV 2-dim)", "public points", 50, || {
            zoda_bls::endomorphism::mul_g1_public(&p1.to_projective(), &k)
        }));
        results.push(bench("G2 scalar mul (plain, 255-bit)", "reference path", 20, || {
            q2.to_projective().mul_fr(&k)
        }));
        results.push(bench("G2 scalar mul (GLV 2-dim)", "public points", 20, || {
            zoda_bls::endomorphism::mul_g2_public(&q2.to_projective(), &k)
        }));
        results.push(bench("hash_to_curve G2 (psi-chain cofactor)", "RFC 9380 RO suite", 50, || {
            zoda_bls::hash::hash_to_curve_g2(b"glv bench message", zoda_bls::sig::DST)
        }));
        // batched h2c: complex-method sqrt + Montgomery batch inversion
        let hmsgs: Vec<Vec<u8>> = (0..128usize)
            .map(|i| format!("glv bench message {}", i).into_bytes())
            .collect();
        results.push(bench("hash_to_curve G2 x128 (batched)", "per message; batch sqrt+inv", 2, || {
            zoda_bls::hash::hash_to_curve_g2_batch(&hmsgs, zoda_bls::sig::DST)
        }));
        // G1 subgroup checks at the verify_cell_kzg_proof_batch scale
        let pts: Vec<zoda_bls::g1::G1Affine> = (0..129usize)
            .map(|i| G1Projective::generator().mul_limbs(&[(i + 1) as u64, 0x9e37, 0, 0]).to_affine())
            .collect();
        results.push(bench("G1 subgroup check x129 (individual)", "pre-v1.3 path", 3, || {
            pts.iter().all(|p| zoda_bls::endomorphism::subgroup_check_g1(p))
        }));
        results.push(bench("G1 subgroup check x129 (batched MSM)", "edas verify path", 20, || {
            zoda_kzg::eip4844::batch_subgroup_check_g1(&pts)
        }));
        // batch BLS verification
        let items: Vec<(zoda_bls::PublicKey, Vec<u8>, zoda_bls::Signature)> = (0..16usize)
            .map(|i| {
                let sk = zoda_bls::SecretKey::from_seed(format!("glv bench key {}", i).as_bytes());
                let pk = sk.public_key();
                let msg = format!("glv bench message {}", i).into_bytes();
                let sig = sk.sign(&msg);
                (pk, msg, sig)
            })
            .collect();
        let refs: Vec<(&zoda_bls::PublicKey, &[u8], &zoda_bls::Signature)> = items
            .iter()
            .map(|(p, m, s)| (p, m.as_slice(), s))
            .collect();
        results.push(bench("BLS verify x16 (individual)", "reference", 2, || {
            refs.iter().all(|(p, m, s)| zoda_bls::verify(p, m, s))
        }));
        results.push(bench("BLS verify_batch x16", "1 final exp amortised", 3, || {
            zoda_bls::verify_batch(&refs)
        }));
        results.push(bench("BLS verify_batch_strict x16", "batch subgroup + GLV", 3, || {
            zoda_bls::verify_batch_strict(&refs)
        }));
    }

    // ---------------- EigenDA-style batch pipeline throughput ----------------
    if filter.is_empty() || filter == "eigenda" {
        println!("building KZG test setup (4096-point, deterministic tau)…");
        let t0 = Instant::now();
        let setup = zoda_kzg::srs::Setup::from_seed_for_testing(*b"bench-eigenda-tau-00000000000000");
        println!("  setup built in {:.1} ms", t0.elapsed().as_secs_f64() * 1000.0);
        const NBLOBS: usize = 8; // 1 MiB batch (EigenDA-style operator duty)
        let mut blobs: Vec<Vec<u8>> = Vec::with_capacity(NBLOBS);
        for _ in 0..NBLOBS {
            let mut b = vec![0u8; 131072];
            rng.next_bytes(&mut b);
            for chunk in b.chunks_exact_mut(32) {
                chunk[0] &= 0x3f;
            }
            blobs.push(b);
        }
        // commit phase
        let t0 = Instant::now();
        let comms: Vec<[u8; 48]> = blobs
            .iter()
            .map(|b| zoda_kzg::eip4844::blob_to_kzg_commitment(b, &setup).unwrap())
            .collect();
        let commit_s = t0.elapsed().as_secs_f64();
        // extend + FK20 prove phase
        let t0 = Instant::now();
        let mut all_cells: Vec<zoda_edas::Cell> = Vec::with_capacity(NBLOBS * 128);
        let mut all_proofs: Vec<[u8; 48]> = Vec::with_capacity(NBLOBS * 128);
        let mut all_indices: Vec<u64> = Vec::with_capacity(NBLOBS * 128);
        let mut all_comms: Vec<&[u8]> = Vec::with_capacity(NBLOBS * 128);
        for (i, b) in blobs.iter().enumerate() {
            let (cells, proofs) = zoda_edas::compute_cells_and_kzg_proofs(b, &setup).unwrap();
            for (ci, (cell, proof)) in cells.into_iter().zip(proofs.into_iter()).enumerate() {
                all_cells.push(cell);
                all_proofs.push(proof);
                all_indices.push(ci as u64);
                all_comms.push(&comms[i]);
            }
        }
        let prove_s = t0.elapsed().as_secs_f64();
        // batch verify phase (attester duty)
        let t0 = Instant::now();
        let proof_refs: Vec<&[u8]> = all_proofs.iter().map(|p| p.as_slice()).collect();
        let ok = zoda_edas::verify_cell_kzg_proof_batch(
            &all_comms,
            &all_indices,
            &all_cells,
            &proof_refs,
            &setup,
        )
        .unwrap();
        let verify_s = t0.elapsed().as_secs_f64();
        assert!(ok, "eigenda pipeline verification failed");
        let total_s = commit_s + prove_s + verify_s;
        let mib = (NBLOBS as f64) * 131072.0 / 1_048_576.0;
        println!();
        println!("  EigenDA-style batch pipeline ({} blobs, {:.0} MiB):", NBLOBS, mib);
        println!("    commit (blob->KZG)        {:>8.1} ms   {:.2} MB/s", commit_s * 1e3, mib * 1_048_576.0 / commit_s / 1e6);
        println!("    extend + FK20 prove       {:>8.1} ms   {:.2} MB/s", prove_s * 1e3, mib * 1_048_576.0 / prove_s / 1e6);
        println!("    batch verify ({} cells) {:>8.1} ms   {:.2} MB/s", all_cells.len(), verify_s * 1e3, mib * 1_048_576.0 / verify_s / 1e6);
        println!("    end-to-end                {:>8.1} ms   {:.2} MB/s", total_s * 1e3, mib * 1_048_576.0 / total_s / 1e6);
        results.push(Timed {
            name: format!("eigenda pipeline ({} blobs e2e)", NBLOBS),
            ns: (total_s * 1e9) as u128,
            note: format!("{:.2} MB/s end-to-end", mib * 1_048_576.0 / total_s / 1e6),
        });
    }

    // ---------------- KZG / EIP-4844 ----------------
    if filter.is_empty() || filter == "kzg" {
        println!("building KZG test setup (4096-point, deterministic tau)…");
        let t0 = Instant::now();
        let setup = zoda_kzg::srs::Setup::from_seed_for_testing(*b"bench-kzg-tau-000000000000000000");
        println!("  setup built in {:.1} ms", t0.elapsed().as_secs_f64() * 1000.0);
        let mut blob = [0u8; 131072];
        rng.next_bytes(&mut blob);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f;
        }
        let mut x = rng.next_fr(false);
        results.push(bench("blob_to_kzg_commitment (4096 MSM)", "EIP-4844 hot path", 20, || {
            zoda_kzg::eip4844::blob_to_kzg_commitment(&blob, &setup)
        }));
        let comm = zoda_kzg::eip4844::blob_to_kzg_commitment(&blob, &setup).unwrap();
        results.push(bench("compute_kzg_proof", "quotient + MSM", 20, || {
            let le = x.to_le_bytes();
            let mut z = [0u8; 32]; // z is I2OSP'd big-endian
            for i in 0..32 {
                z[i] = le[31 - i];
            }
            x = x + Fr::ONE;
            zoda_kzg::eip4844::compute_kzg_proof(&blob, &z, &setup)
        }));
        let z = {
            let mut t = [0u8; 32];
            t[31] = 7;
            t
        };
        let (proof, y) = zoda_kzg::eip4844::compute_kzg_proof(&blob, &z, &setup).unwrap();
        results.push(bench("verify_kzg_proof (single)", "2 pairings", 20, || {
            zoda_kzg::eip4844::verify_kzg_proof(&comm, &z, &y, &proof, &setup)
        }));
        // batch
        let mut blobs = Vec::new();
        for i in 0..6 {
            let mut b = [0u8; 131072];
            rng.next_bytes(&mut b);
            for chunk in b.chunks_exact_mut(32) {
                chunk[0] &= 0x3f;
            }
            blobs.push(b);
            let _ = i;
        }
        let comms: Vec<_> = blobs
            .iter()
            .map(|b| zoda_kzg::eip4844::blob_to_kzg_commitment(b, &setup).unwrap())
            .collect();
        let proofs: Vec<_> = blobs
            .iter()
            .zip(comms.iter())
            .map(|(b, c)| zoda_kzg::eip4844::compute_blob_kzg_proof(b, c, &setup).unwrap())
            .collect();
        let bl: Vec<&[u8]> = blobs.iter().map(|b| b.as_slice()).collect();
        let cl: Vec<&[u8]> = comms.iter().map(|c| c.as_slice()).collect();
        let pl: Vec<&[u8]> = proofs.iter().map(|p| p.as_slice()).collect();
        results.push(bench("verify_blob_kzg_proof_batch (6 blobs)", "3 MSM + 1 pairing", 10, || {
            zoda_kzg::eip4844::verify_blob_kzg_proof_batch(&bl, &cl, &pl, &setup)
        }));
    }

    // ---------------- EIP-7594 cells + FK20 ----------------
    if filter.is_empty() || filter == "edas" {
        println!("building EDAS setup…");
        let setup = zoda_kzg::srs::Setup::from_seed_for_testing(*b"bench-edas-tau-00000000000000000");
        let mut blob = [0u8; 131072];
        rng.next_bytes(&mut blob);
        for chunk in blob.chunks_exact_mut(32) {
            chunk[0] &= 0x3f;
        }
        results.push(bench("compute_cells (FFT 8192 + BRP)", "EIP-7594 extension", 20, || {
            zoda_edas::compute_cells(&blob, &setup)
        }));
        results.push(bench("compute_cells_and_kzg_proofs (FK20)", "128 cell proofs, O(n log n)", 5, || {
            zoda_edas::compute_cells_and_kzg_proofs(&blob, &setup)
        }));
        let (cells, proofs) = zoda_edas::compute_cells_and_kzg_proofs(&blob, &setup).unwrap();
        let comm = zoda_kzg::eip4844::blob_to_kzg_commitment(&blob, &setup).unwrap();
        let idx: Vec<u64> = (0..128).collect();
        let cr: Vec<&[u8]> = idx.iter().map(|_| comm.as_slice()).collect();
        let pr: Vec<&[u8]> = proofs.iter().map(|p| p.as_slice()).collect();
        results.push(bench("verify_cell_kzg_proof_batch (128 cells)", "1 pairing amortised", 5, || {
            zoda_edas::verify_cell_kzg_proof_batch(&cr, &idx, &cells, &pr, &setup)
        }));
        let keep_idx: Vec<u64> = (0..64).collect();
        let keep: Vec<[u8; 2048]> = keep_idx.iter().map(|&i| cells[i as usize]).collect();
        results.push(bench("recover_cells_and_kzg_proofs (from 64)", "erasure decode + FK20", 3, || {
            zoda_edas::recover_cells_and_kzg_proofs(&keep_idx, &keep, &setup)
        }));
    }

    // ---------------- ZODA tensor core ----------------
    if filter.is_empty() || filter == "core" {
        for (m, k) in [(16usize, 16usize), (64, 64)] {
            let mut data = vec![Fr::zero(); m * k];
            for x in data.iter_mut() {
                *x = rng.next_fr(false);
            }
            let grid = Matrix::from_row_major(m, k, data);
            let params = ZodaParams::new(m, k);
            results.push(bench(
                &format!("ZODA commit ({}x{} grid, full encode)", m, k),
                "2x NTT passes + Merkle + FS",
                5,
                || params.commit(&grid),
            ));
            let prover = params.commit(&grid);
            let public = prover.public_params();
            let row = prover.matrix.row(3);
            results.push(bench(
                &format!("ZODA verify row sample ({}x{})", m, k),
                "O(2k) inner product",
                100,
                || zoda_core::verify_row_sample(&public, 3, &row),
            ));
        }
        // reconstruction from half
        let (m, k) = (32usize, 32usize);
        let mut data = vec![Fr::zero(); m * k];
        for x in data.iter_mut() {
            *x = rng.next_fr(false);
        }
        let grid = Matrix::from_row_major(m, k, data);
        let params = ZodaParams::new(m, k);
        let prover = params.commit(&grid);
        let public = prover.public_params();
        let kept: Vec<usize> = (0..k).collect();
        results.push(bench(
            &format!("ZODA reconstruct from {} cols ({}x{})", k, m, k),
            "per-row interpolation",
            3,
            || zoda_core::reconstruct(&public, &kept, &prover.matrix),
        ));
    }

    // ---------------- 2D erasure coding + DAS deep dive ----------------
    if filter.is_empty() || filter == "das2d" {
        // RS vector encode across the FFT-path sizes
        for log in [7u32, 9, 11] {
            let n = 1usize << log;
            let domain = FftDomain::<Fr>::new(2 * n);
            let data: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
            let mb = (n * 32 * 2) as f64;
            results.push(bench(
                &format!("rs_encode_vector n={}", n),
                "O(n log n) FFT interpolation",
                20,
                || zoda_core::encode::rs_encode_vector(&data, &domain),
            ));
            let last = results.last().unwrap();
            println!(
                "    -> rs_encode_vector n={} throughput: {}",
                n,
                fmt_bytes(mb / (last.ns as f64 / 1e9))
            );
        }
        // schoolbook vs FFT at n = 1024 (the optimization this session)
        {
            let n = 1024usize;
            let data: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
            let enc = zoda_core::encode::geometric_encoder_public(n);
            results.push(bench(
                "interp n=1024 [schoolbook O(n^2)]",
                "power sums + correlation loops",
                3,
                || enc.interpolate_fast(&data),
            ));
            results.push(bench(
                "interp n=1024 [FFT O(n log n)]",
                "3 size-2n transforms",
                20,
                || enc.interpolate_fft(&data),
            ));
        }
        // transpose
        {
            let n = 256usize;
            let data: Vec<Fr> = (0..n * n).map(|_| rng.next_fr(false)).collect();
            let m = Matrix::from_row_major(n, n, data);
            results.push(bench(
                &format!("Matrix transpose {}x{}", n, n),
                "blocked, parallel",
                20,
                || m.transpose(),
            ));
        }
        // 2D scattered-erasure reconstruction
        {
            let (m, k) = (32usize, 32usize);
            let mut data = vec![Fr::zero(); m * k];
            for x in data.iter_mut() {
                *x = rng.next_fr(false);
            }
            let grid = Matrix::from_row_major(m, k, data);
            let params = ZodaParams::new(m, k);
            let prover = params.commit(&grid);
            let public = prover.public_params();
            let mut rng2 = ZodaRng::from_seed(*b"das2d-bench-seed-000000000000000");
            let mut erased = Vec::new();
            for r in 0..64 {
                for c in 0..64 {
                    if rng2.next_u64() % 10 < 6 {
                        erased.push((r, c));
                    }
                }
            }
            results.push(bench(
                "reconstruct_2d 64x64 (60% cells erased)",
                "row/col fixpoint + full verify",
                3,
                || {
                    let mut partial = zoda_core::PartialGrid::with_erasures(&prover.matrix, &erased);
                    zoda_core::reconstruct_2d(&mut partial, &params, &public)
                },
            ));
        }
        // attested multi-peer session (16x16 grid)
        {
            let (m, k) = (16usize, 16usize);
            let mut data = vec![Fr::zero(); m * k];
            for x in data.iter_mut() {
                *x = rng.next_fr(false);
            }
            let grid = Matrix::from_row_major(m, k, data);
            let params = ZodaParams::new(m, k);
            let prover = params.commit(&grid);
            let public = prover.public_params();
            results.push(bench(
                "attested DAS session 16x16 (1 peer)",
                "projection + Merkle, adaptive rounds",
                5,
                || {
                    let mut oracle = zoda_das::ProverOracle::honest(&prover);
                    let mut peers: [&mut dyn zoda_das::AttestedOracle; 1] = [&mut oracle];
                    let mut r = ZodaRng::from_seed(*b"das2d-sess-seed-0000000000000000");
                    zoda_das::run_attested_session(
                        &public,
                        &params,
                        &mut peers,
                        &zoda_das::SessionConfig::default(),
                        &mut r,
                    )
                },
            ));
        }
        // archival: custody ingest + reconstruct + serialize
        {
            let (m, k) = (32usize, 32usize);
            let mut data = vec![Fr::zero(); m * k];
            for x in data.iter_mut() {
                *x = rng.next_fr(false);
            }
            let grid = Matrix::from_row_major(m, k, data);
            let params = ZodaParams::new(m, k);
            let prover = params.commit(&grid);
            let public = prover.public_params();
            results.push(bench(
                "custody put_column_verified 32x32",
                "projection check + store",
                200,
                || {
                    let mut c = zoda_archival::GridCustody::new(
                        [1u8; 32],
                        1,
                        ZodaParams::new(m, k),
                        public.clone(),
                    );
                    c.put_column_verified(0, prover.matrix.col(0)).unwrap();
                    c
                },
            ));
            results.push(bench(
                "custody try_reconstruct 32x32 (k cols)",
                "fixpoint decode + persist + verify",
                3,
                || {
                    let mut c = zoda_archival::GridCustody::new(
                        [1u8; 32],
                        1,
                        ZodaParams::new(m, k),
                        public.clone(),
                    );
                    for i in 0..k {
                        c.put_column_verified(i, prover.matrix.col(i)).unwrap();
                    }
                    c.try_reconstruct().unwrap()
                },
            ));
        }
    }

    // ---------------- post-quantum lattice ----------------
    if filter.is_empty() || filter == "pq" {
        for params in [
            zoda_pq::PqParams::fast(),
            zoda_pq::PqParams::nist_l1(),
        ] {
            let deg = 4 * params.n;
            let poly: Vec<zoda_math::Fq> =
                (0..deg).map(|_| zoda_math::Fq(rng.next_below(zoda_pq::rq::Q as u64) as u32)).collect();
            let note = format!("n={}, {} chunks", params.n, (deg + params.n - 1) / params.n);
            results.push(bench(
                &format!("PQ commit [{}]", params.level),
                &format!("MLWE commit; {}", note),
                10,
                || {
                    let mut r = ZodaRng::from_seed(*b"pq-bench-seed-000000000000000000");
                    zoda_pq::LatticePcs::commit(&params, &poly, &mut r).0
                },
            ));
            let (com, key) = {
                let mut r = ZodaRng::from_seed(*b"pq-bench-seed-000000000000000000");
                zoda_pq::LatticePcs::commit(&params, &poly, &mut r)
            };
            let zeta = zoda_pq::ring_point_from_seed(b"bench-point", params.n);
            let v = zoda_pq::LatticePcs::eval_at_ring(&poly, &zeta);
            results.push(bench(
                &format!("PQ open [{}]", params.level),
                &format!("sigma-protocol + FS; {}", note),
                5,
                || {
                    let mut r = ZodaRng::from_seed(*b"pq-open-seed-0000000000000000000");
                    zoda_pq::LatticePcs::open(&params, &key, &zeta, &v, &mut r)
                },
            ));
            let mut r = ZodaRng::from_seed(*b"pq-open-seed-0000000000000000000");
            let proof = zoda_pq::LatticePcs::open(&params, &key, &zeta, &v, &mut r);
            results.push(bench(
                &format!("PQ verify [{}]", params.level),
                &format!("2 equations; {}", note),
                10,
                || zoda_pq::LatticePcs::verify(&params, &com, &zeta, &v, &proof),
            ));
        }
    }

    // ---------------- v1.4: mainnet EigenDA pipeline throughput ----------------
    if filter.is_empty() || filter == "throughput" {
        use zoda_edas::pipeline::{DisperserPipeline, EigenDaConfig, RetrieverPipeline};
        let cores = std::thread::available_parallelism().map(|c| c.get()).unwrap_or(1);
        println!("EigenDA mainnet pipeline throughput ({} cores)", cores);
        // real mainnet configuration; the real mainnet ceremony SRS when
        // the fetch script has been run, the deterministic test tau
        // otherwise (benchmarked identically — the SRS content does not
        // change MSM costs)
        let cfg = EigenDaConfig::from_toml_file(std::path::Path::new(
            "config/eigenda/mainnet.toml",
        ))
        .unwrap_or_else(|_| EigenDaConfig::default());
        let setup_path = std::path::Path::new("spec-vectors/trusted_setup.txt");
        let (setup, srs_note) = if setup_path.exists() {
            println!("  loading REAL mainnet ceremony trusted setup…");
            let t0 = Instant::now();
            let s = zoda_kzg::srs::Setup::load_file(setup_path).expect("mainnet SRS loads");
            println!(
                "  setup loaded + verified in {:.0} ms",
                t0.elapsed().as_secs_f64() * 1e3
            );
            (s, "mainnet ceremony SRS")
        } else {
            println!("  (spec-vectors/trusted_setup.txt absent — deterministic test tau)");
            (
                zoda_kzg::srs::Setup::from_seed_for_testing(
                    *b"bench-throughput-tau-00000000000",
                ),
                "deterministic test tau",
            )
        };
        let disperser = DisperserPipeline::new(cfg.clone(), &setup);
        let retriever = RetrieverPipeline::new(cfg.clone(), &setup);

        // one EigenDA v1 mainnet blob = 2 MiB = 16 zoda cell-blobs;
        // benchmark 16 and 64 payload blobs (2 and 8 MiB)
        for nblobs in [16usize, 64] {
            let mut blobs: Vec<Vec<u8>> = Vec::with_capacity(nblobs);
            for _ in 0..nblobs {
                let mut b = vec![0u8; 131072];
                rng.next_bytes(&mut b);
                for chunk in b.chunks_exact_mut(32) {
                    chunk[0] &= 0x3f;
                }
                blobs.push(b);
            }
            let mib = nblobs as f64 / 8.0;
            let bytes = nblobs as f64 * 131072.0;
            let refs: Vec<&[u8]> = blobs.iter().map(|b| b.as_slice()).collect();

            // disperser: commit + extend + FK20 prove (internal core-level
            // parallelism; blobs processed sequentially)
            let (dispersals, timings) = disperser.disperse_batch(&refs).unwrap();
            let commit = &timings[0];
            let prove = &timings[1];
            let total = &timings[2];

            // operator: batch verify every cell of the batch
            let mut cells: Vec<zoda_edas::Cell> = Vec::with_capacity(nblobs * 128);
            let mut proofs: Vec<[u8; 48]> = Vec::with_capacity(nblobs * 128);
            let mut idx: Vec<u64> = Vec::with_capacity(nblobs * 128);
            let mut comms: Vec<[u8; 48]> = Vec::with_capacity(nblobs * 128);
            for d in &dispersals {
                comms.extend(std::iter::repeat(d.commitment).take(128));
                cells.extend_from_slice(&d.cells);
                proofs.extend_from_slice(&d.proofs);
                idx.extend((0..128u64).map(|i| i));
            }
            let t0 = Instant::now();
            let ok = retriever.verify_cells(&comms, &cells, &proofs, &idx).unwrap();
            let verify = t0.elapsed().as_secs_f64();
            assert!(ok, "throughput pipeline verification failed");

            // custody: reconstruction of every blob from its custody half,
            // threaded across blobs (per-blob recovery is serial)
            let t0 = Instant::now();
            let recovered: Vec<usize> = {
                let halves: Vec<Vec<zoda_edas::Cell>> = cells
                    .chunks(128)
                    .map(|c| c.iter().step_by(2).copied().collect())
                    .collect();
                let half_idx: Vec<u64> = (0..64u64).map(|i| i * 2).collect();
                let retriever_ref = &retriever;
                let half_idx_ref = &half_idx;
                std::thread::scope(|scope| {
                    let handles: Vec<_> = halves
                        .chunks(cores.max(1))
                        .map(|h| {
                            scope.spawn(move || {
                                h.iter()
                                    .map(|half| {
                                        retriever_ref
                                            .recover(half, half_idx_ref)
                                            .map(|(c, _)| c.len())
                                            .unwrap_or(0)
                                    })
                                    .collect::<Vec<usize>>()
                            })
                        })
                        .collect();
                    handles
                        .into_iter()
                        .flat_map(|h| h.join().unwrap())
                        .collect()
                })
            };
            let rec_secs = t0.elapsed().as_secs_f64();
            assert!(
                recovered.iter().all(|n| *n == 128),
                "recovery returned wrong cell counts"
            );

            let e2e_secs = total.secs + verify;
            let node_secs = verify + rec_secs;
            println!(
                "  {} x 128 KiB payloads ({:.0} MiB, {}):",
                nblobs, mib, srs_note
            );
            println!(
                "    commit (blob -> KZG)      {:>8.1} ms   {:>7.2} MB/s",
                commit.secs * 1e3,
                commit.mb_per_s()
            );
            println!(
                "    extend + FK20 prove       {:>8.1} ms   {:>7.2} MB/s",
                prove.secs * 1e3,
                prove.mb_per_s()
            );
            println!(
                "    batch verify ({} cells) {:>8.1} ms   {:>7.2} MB/s",
                cells.len(),
                verify * 1e3,
                bytes / verify / 1e6
            );
            println!(
                "    custody reconstruct (50%) {:>7.1} ms   {:>7.2} MB/s",
                rec_secs * 1e3,
                bytes / rec_secs / 1e6
            );
            println!(
                "    disperser e2e             {:>8.1} ms   {:>7.2} MB/s",
                total.secs * 1e3,
                total.mb_per_s()
            );
            println!(
                "    node e2e (verify+recover) {:>7.1} ms   {:>7.2} MB/s",
                node_secs * 1e3,
                bytes / node_secs / 1e6
            );
            println!(
                "    full e2e (commit+prove+verify) {:.1} ms  {:.2} MB/s",
                e2e_secs * 1e3,
                bytes / e2e_secs / 1e6
            );
            results.push(Timed {
                name: format!("mainnet disperser e2e ({}x128KiB)", nblobs),
                ns: (total.secs * 1e9) as u128,
                note: format!("{:.2} MB/s", total.mb_per_s()),
            });
            results.push(Timed {
                name: format!("mainnet node e2e ({}x128KiB)", nblobs),
                ns: (node_secs * 1e9) as u128,
                note: format!("{:.2} MB/s", bytes / node_secs / 1e6),
            });
        }
    }

    println!();
    println!("{:=<100}", "=");
    println!("{:<58} {:>12} {:>10}   {}", "BENCHMARK", "TIME", "UNIT", "NOTES");
    println!("{:=<100}", "=");
    for r in &results {
        r.print();
    }
    println!("{:=<100}", "=");
    println!("host: {} cores", std::thread::available_parallelism().map(|c| c.get()).unwrap_or(1));
}
