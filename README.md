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

| operation | time | note |
|---|---|---|
| `Fr` multiplication (ADX/BMI2 Montgomery) | **30 ns** | runtime-dispatched `mulx`/`adcx`/`adox` |
| pairing (Miller + fountain final exp) | **1.6 ms** | lazy-reduction Fp2 |
| BLS verify (2-term pairing_check) | **2.9 ms** | 1.4x vs v1.2 |
| BLS `verify_batch` x16 | **15.6 ms** | 2.2x (batched h2c + parallel Miller accumulation) |
| BLS sign | **1.0 ms** | complex-method SSWU + Jacobian psi-chain |
| G1 scalar mul (GLV 2-dim, public points) | **133 µs** | derived beta/lambda endomorphism |
| `blob_to_kzg_commitment` (4096 MSM) | 65 ms | Pippenger, 2 threads |
| `verify_cell_kzg_proof_batch` (128 cells) | **31.8 ms** | batched subgroup checks (4.3x on the checks) |
| ZODA verify row sample (64x64 grid) | **2.6 µs** | |
| ZODA commit (64x64 grid) | **15.5 ms** | |
| `reconstruct_2d` 64x64, 60% cells erased | **10.8 ms** | parallel row/column fixpoint passes |
| EigenDA-style batch verify (1024 cells) | **8.45 MB/s** | operator attestation throughput |
| hash_to_curve G2 batched | **378 µs/msg** | complex-method sqrt + Montgomery batch inversion |
| mainnet pipeline node e2e (8 MiB) | **0.49 MB/s** | verify + custody reconstruction, real ceremony SRS |

Full tables, methodology and per-version history: [docs/BENCHMARKS.md](docs/BENCHMARKS.md).

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
* Since v1.3 the field layer closes most of the portable-Rust vs assembly gap with runtime-dispatched ADX/BMI2 Montgomery paths (`mulx`/`adcx`/`adox` codegen, no hand-written assembly, no external crates); GLV endomorphisms, lazy-reduction Fp2, batch verification and parallel 2D decode are documented in [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md); v1.4 adds the complex-method SSWU hash-to-curve, parallel Miller-loop accumulation and the real-mainnet EigenDA pipeline ([`config/eigenda/mainnet.toml`](config/eigenda/mainnet.toml), verified endpoints + Ethereum ceremony SRS).

## License

Apache-2.0. Protocol references: [ZODA (eprint 2025/034)](https://eprint.iacr.org/2025/034), EIP-4844, EIP-7594, RFC 9380, BDLOP/Lyubashevsky lattice commitments.
