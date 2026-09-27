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
        let mut y = rng.next_fr(false);
        results.push(bench("Fr mul (Montgomery, 4 limbs)", "field core", 1_000_000, || {
            x = x * y;
            x
        }));
        let mut xs = rng.next_fr(false);
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
            let mut z = [0u8; 32];
            z.copy_from_slice(&x.to_le_bytes());
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
