# zoda

**Zero-Overhead Data Availability** — a production-grade Rust implementation of the ZODA tensor-code data availability protocol ([eprint 2025/034](https://eprint.iacr.org/2025/034)), with full **Ethereum EIP-4844 / EIP-7594** integration, a **post-quantum lattice-based polynomial commitment** variant, sampling, custody, archival and bridge services.

```
tensor-code DA core  ·  KZG (EIP-4844)  ·  cells + FK20 (EIP-7594)  ·  lattice PQ commitments
sampling & custody   ·  sybil resistance  ·  archival reconstruction  ·  light-client bridges
```

[![CI](https://github.com/idrees2516/zoda/actions/workflows/ci.yml/badge.svg)](./.github/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](./LICENSE)
![Rust](https://img.shields.io/badge/rust-1.75+-orange)

## Why

Data availability sampling (DAS) lets light nodes verify that *block data was published* without downloading it. Today's designs pay for that confidence in either bandwidth (sample data) or trust (committee signatures). **ZODA** makes the samples *prove their own correctness*: a sampled row or column of a slightly-modified tensor encoding is simultaneously data **and** proof — zero incremental communication, a few extra field operations per sample.

This repository turns the research prototype into a production system:

| The paper's prototype | this implementation |
|---|---|
| O(n²) Vandermonde-matrix encoding | **O(n log n) NTT-based** systematic encoding |
| random vectors unbound to the data | **Fiat–Shamir bound** to the row/column commitments |
| `decode` was literally `!unimplemented!()` | **full erasure reconstruction** from ≥ ½ columns |
| no peer layer | sampling, custody, sybil resistance, archival, bridges |
| pairing-only commitments | KZG (EIP-4844/7594) **+ post-quantum lattice variant** |

## Layout

```
crates/
├── zoda-math      fields (Fr/Goldilocks/Fq), NTT, poly trees, SHA-256/Keccak, Merkle, RNG, stats
├── zoda-bls       BLS12-381: Fp2/Fp6/Fp12 tower, G1/G2, optimal ate pairing, RFC 9380
│                  hash-to-curve, BLS signatures — zero dependencies
├── zoda-kzg       EIP-4844 evaluation-form KZG, trusted setup (mainnet format), Pippenger MSM
├── zoda-core      ZODA tensor protocol: O(n log n) encoding, commitments, sample
│                  proofs, full-column + arbitrary-cell 2D reconstruction
├── zoda-pq        lattice post-quantum polynomial commitments (Module-LWE/SIS, BDLOP + Lyubashevsky)
├── zoda-edas      EIP-7594: DAS extension, 128 cells, FK20 proofs, batch verify, recovery, custody
├── zoda-das       2D sampling: exact availability theory, attested multi-peer line
│                  sessions, EIP-7594 cell sessions with column custody
├── zoda-sybils    BLS sortition, stake-weighted selection, peer scoring
├── zoda-rda       adaptive randomized-DA sessions with confidence accounting
├── zoda-archival  custody storage: verify-on-insert, persisting reconstruction,
│                  compact checksummed files, atomic disk writes, budget pruning
├── zoda-bridges   light-client inclusion proofs, KZG-backed bridge messages
├── zoda-ethrex    ethrex-compatible blob-transaction sidecars
├── zoda           facade re-exporting the whole stack
└── zoda-bench     the benchmark suite (`cargo run --release -p zoda-bench`)
```

## Quick start

```bash
cargo test --workspace --release    # 159 tests incl. 320 official spec vectors (when fetched)
cargo run --release -p zoda-bench  # full benchmark table
```

Commit a data grid and verify a sample:

```rust
use zoda::core::{Matrix, ZodaParams};
use zoda::math::{Fr, PrimeField, ZodaRng};

let mut rng = ZodaRng::from_seed(*b"demo-seed-00000000000000000000000000");
let mut data = vec![Fr::zero(); 64 * 64];
for x in data.iter_mut() { *x = rng.next_fr(false); }
let grid = Matrix::from_row_major(64, 64, data);

let params = ZodaParams::new(64, 64);
let prover = params.commit(&grid);          // encode + commit + Fiat–Shamir
let public = prover.public_params();

let row = prover.matrix.row(17);            // the sample is its own proof
assert!(zoda::core::verify_row_sample(&public, 17, &row).is_ok());
```

EIP-4844 blobs:

```rust
use zoda::kzg::srs::Setup;

let setup = Setup::load_file(std::path::Path::new("trusted_setup.txt"))?; // mainnet format
let blob = /* 131072 bytes, each 32-byte element < r */;
let commitment = zoda::kzg::blob_to_kzg_commitment(&blob, &setup)?;
let (proof, y)   = zoda::kzg::compute_kzg_proof(&blob, &z, &setup)?;
assert!(zoda::kzg::verify_kzg_proof(&commitment, &z, &y, &proof, &setup)?);
```

EIP-7594 cells (PeerDAS):

```rust
let (cells, proofs) = zoda::edas::compute_cells_and_kzg_proofs(&blob, &setup)?; // 128 cells
assert!(zoda::edas::verify_cell_kzg_proof_batch(&[&commitment], &[0, 1, 42],
    &[cells[0], cells[1], cells[42]], &[&proofs[0], &proofs[1], &proofs[42]], &setup)?);
// reconstruct everything from any 64 of 128 cells:
let (all_cells, all_proofs) = zoda::edas::recover_cells_and_kzg_proofs(&keep_idx, &keep_cells, &setup)?;
```

Post-quantum polynomial commitments:

```rust
let params = zoda::pq::PqParams::nist_l1();
let (com, key) = zoda::pq::LatticePcs::commit(&params, &poly, &mut rng); // no trusted setup
let zeta = zoda::pq::ring_point_from_seed(b"point", params.n);
let v = zoda::pq::LatticePcs::eval_at_ring(&poly, &zeta);
let proof = zoda::pq::LatticePcs::open(&params, &key, &zeta, &v, &mut rng);
assert!(zoda::pq::LatticePcs::verify(&params, &com, &zeta, &v, &proof));
```

## Performance (2 vCPU host, release + LTO)

| operation | time | vs v1.0.0 |
|---|---|---|
| `Fr` multiplication | **25 ns** | — |
| NTT-8192 | 1.9 ms | — |
| pairing (Miller + fountain final exp) | 2.0 ms | — |
| BLS verify (2 pairings) | 4.1 ms | 1.1x |
| `blob_to_kzg_commitment` (4096 MSM) | 69 ms | **5.6x** |
| ZODA verify row sample (64×64 grid) | **2.6 µs** | — |
| ZODA commit (64×64 grid) | **15.5 ms** | **16.5x** |
| ZODA reconstruct (32 of 64 columns) | 2.4 ms | **14.7x** |
| RS interpolation kernel (n = 1024) | **1.2 ms** | **49x vs the O(n²) loops** |
| `reconstruct_2d` 64×64, 60% cells erased | 16 ms | new: arbitrary cell loss + full verify |
| attested DAS session (16×16, 1 peer) | **438 µs** | new: projections + Merkle + peer scores |
| custody ingest (verified column) | 297 µs | new |
| `compute_cells` (EIP-7594) | 2.9 ms | 1.2x |
| FK20 (128 cell proofs) | 0.31 s | **6.0x** |
| `verify_cell_kzg_proof_batch` (128 cells) | 49 ms† | 1.7x† |
| recover 64→128 cells | 0.28 s | **6.6x** |
| PQ commit / open / verify (L1) | 5.3 / 13 / 13 ms | **9.5 / 4.6 / 3.0x** |

† now includes the spec-mandated subgroup validation of all 129 parsed
points (~24 ms on this stack); excluding it the batch check runs in
~25 ms (3.4x).

Full table, before/after comparison and the optimization notes:
[`docs/BENCHMARKS.md`](docs/BENCHMARKS.md).

## Documentation

| doc | contents |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | system design, crate graph, data flows |
| [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md) | the math inside every hot path |
| [`docs/POST_QUANTUM.md`](docs/POST_QUANTUM.md) | the lattice commitment: construction, soundness, parameters, limits |
| [`docs/ETHEREUM_COMPLIANCE.md`](docs/ETHEREUM_COMPLIANCE.md) | EIP-4844 / EIP-7594 conformance, validation status |
| [`docs/ETHEX_INTEGRATION.md`](docs/ETHEX_INTEGRATION.md) | wiring into the ethrex execution client |
| [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md) | methodology + full results |
| [`docs/DEPLOYMENT.md`](docs/DEPLOYMENT.md) | building, deploying, Docker, ops notes |
| [`SECURITY.md`](SECURITY.md) | threat model and security notes |

## Status & honesty

* **Official consensus-spec-tests conformance: 320/320 vector cases pass bit-exactly** against the real mainnet ceremony setup — all six Deneb EIP-4844 suites (253 cases) and all four Fulu EIP-7594 cell suites (67 cases). Run it yourself: `scripts/fetch_kzg_vectors.sh && cargo test -p zoda-kzg -p zoda-edas --release spec_vectors`. The harness caught and fixed three real spec deviations (challenge-transcript encoding, missing subgroup validation, non-canonical cell handling).
* All Ethereum-flavoured crypto is validated against **official vectors**: RFC 9380 hash-to-curve (G1+G2), the mainnet KZG trusted setup parse + spec sanity check, and the full EIP-4844/7594 vector suites above.
* The lattice PQ commitment is a **research prototype**: parameter sets are first estimates, not audited. See `docs/POST_QUANTUM.md`.
* The group operations (Jacobian engine, Pippenger MSM, FK20 with Straus tables) are optimized pure Rust with no assembly and no unsafe code; they remain behind blst/c-kzg's ADX-assembly paths by roughly the expected asm/no-asm margin. Hot paths are documented with optimization notes in [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md).

## License

Apache-2.0. Protocol references: [ZODA (eprint 2025/034)](https://eprint.iacr.org/2025/034), EIP-4844, EIP-7594, RFC 9380, BDLOP/Lyubashevsky lattice commitments.
