# Benchmarks

## Methodology

* Host: 2 vCPU cloud instance, rustc 1.98.1, release profile with
  thin-LTO and codegen-units=1. Numbers move ±20% between runs on this
  box; the table below is a representative full-suite run.
* Harness: `cargo run --release -p zoda-bench` — warm-up pass, then the
  mean over the listed iteration count (deterministic seeded inputs; no
  OS entropy in the measurements).
* Filter with `zoda-bench math`, `bls`, `kzg`, `edas`, `core`, `das2d` or
  `pq` to run a single group.
* KZG/EDAS runs use a deterministic 4096-point toy setup (`tau` from a
  seed; the mainnet setup file yields the same shapes — point counts and
  widths are identical), and the **official consensus-spec-tests vector
  suite replays bit-exactly against the real mainnet ceremony setup**
  (`scripts/fetch_kzg_vectors.sh` + `cargo test --release --workspace`).
* Verification paths include the **spec-mandated subgroup validation**
  of every parsed G1 point (`[r]P = O`, matching c-kzg's
  `validate_kzg_g1`) — roughly 0.19 ms per commitment/proof on this
  pure-Rust stack.

## Results (v1.2.0 — 2D DAS + erasure-coding deep dive)

New this release: the O(n log n) RS interpolation kernel (both quadratic
loops of the transposed-Vandermonde identity became butterfly
transforms), the transposed cache-friendly column encode, the 2D
scattered-erasure fixpoint decoder, attested multi-peer sampling
sessions, cell-level blob-grid sessions, and the custody storage layer.

```
====================================================================================================
BENCHMARK                                                          TIME       UNIT   NOTES
====================================================================================================
Fr mul (Montgomery, 4 limbs)                                      0.025         µs   field core [1000000 iters]
Fr invert (Fermat)                                                7.324         µs   batch-invert backbone [500 iters]
NTT size 4096                                                   910.376         µs   262.14 KB/s/iter [200 iters]
NTT size 8192                                                  1951.490         µs   524.29 KB/s/iter [200 iters]
G1 scalar mul (windowed, 255-bit)                               215.784         µs   sign/pk path [50 iters]
BLS sign (hash_to_curve G2 + mul)                              2125.106         µs   ETH2 ciphersuite [20 iters]
BLS verify (2 pairings)                                        4085.512         µs   ETH2 ciphersuite [10 iters]
pairing e(g1, g2) (miller + final exp)                         2014.090         µs   fountain chain [10 iters]
blob_to_kzg_commitment (4096 MSM)                             69009.346         µs   EIP-4844 hot path [20 iters]
compute_kzg_proof                                             24505.301         µs   quotient + MSM [20 iters]
verify_kzg_proof (single)                                      3231.029         µs   2 pairings [20 iters]
verify_blob_kzg_proof_batch (6 blobs)                         15469.866         µs   3 MSM + 1 pairing [10 iters]
compute_cells (FFT 8192 + BRP)                                 2924.138         µs   EIP-7594 extension [20 iters]
compute_cells_and_kzg_proofs (FK20)                          279098.129         µs   128 cell proofs, O(n log n) [5 iters]
verify_cell_kzg_proof_batch (128 cells)                       48751.942         µs   1 pairing amortised [5 iters]
recover_cells_and_kzg_proofs (from 64)                       278123.319         µs   erasure decode + FK20 [3 iters]
ZODA commit (16x16 grid, full encode)                          2194.321         µs   2x NTT passes + Merkle + FS [5 iters]
ZODA verify row sample (16x16)                                    0.693         µs   O(2k) inner product [100 iters]
ZODA commit (64x64 grid, full encode)                         15454.741         µs   2x NTT passes + Merkle + FS [5 iters]
ZODA verify row sample (64x64)                                    2.567         µs   O(2k) inner product [100 iters]
ZODA reconstruct from 32 cols (32x32)                          2432.799         µs   per-row interpolation [3 iters]
rs_encode_vector n=128                                          105.810         µs   O(n log n) FFT interpolation [20 iters]
rs_encode_vector n=512                                          690.404         µs   O(n log n) FFT interpolation [20 iters]
rs_encode_vector n=2048                                        3380.964         µs   O(n log n) FFT interpolation [20 iters]
interp n=1024 [schoolbook O(n^2)]                             59140.728         µs   power sums + correlation loops [3 iters]
interp n=1024 [FFT O(n log n)]                                 1159.739         µs   3 size-2n transforms [20 iters]
Matrix transpose 256x256                                        718.317         µs   blocked, parallel [20 iters]
reconstruct_2d 64x64 (60% cells erased)                       15962.040         µs   row/col fixpoint + full verify [3 iters]
attested DAS session 16x16 (1 peer)                             438.307         µs   projection + Merkle, adaptive rounds [5 iters]
custody put_column_verified 32x32                               296.881         µs   projection check + store [200 iters]
custody try_reconstruct 32x32 (k cols)                         3205.422         µs   fixpoint decode + persist + verify [3 iters]
PQ commit [FAST]                                                256.510         µs   MLWE commit; n=64, 4 chunks [10 iters]
PQ open [FAST]                                                  273.073         µs   sigma-protocol + FS; n=64, 4 chunks [5 iters]
PQ verify [FAST]                                                164.340         µs   2 equations; n=64, 4 chunks [10 iters]
PQ commit [L1]                                                 5288.005         µs   MLWE commit; n=1024, 4 chunks [10 iters]
PQ open [L1]                                                  12902.398         µs   sigma-protocol + FS; n=1024, 4 chunks [5 iters]
PQ verify [L1]                                                12679.530         µs   2 equations; n=1024, 4 chunks [10 iters]
====================================================================================================
```

Highlights of the v1.1.0 → v1.2.0 pass:

* **Interpolation kernel 49x at n = 1024** (59.1 ms → 1.16 ms): the
  power-sum loop is one NTT and the coefficient correlation is an exact
  cyclic convolution with a cached NTT(Q).
* **ZODA commit 64×64: 39 ms → 15.5 ms** (2.5x) — the FFT interpolation
  compounds across all 128 row and 128 column encodes; from v1.0.0 the
  cumulative speedup is 255 ms → 15.5 ms (**16.5x**).
* **reconstruct_2d recovers a 64×64 grid with 60% of its cells erased**
  in 15.9 ms, including the full final projection verification of all
  128 lines.
* **Whole attested sampling session over a 16×16 grid in 438 µs** —
  adaptive rounds, projection + Merkle checks, peer scoring and exact
  confidence accounting included.
* **Custody ingest at 297 µs/column** (verified), full custody
  reconstruction (k columns, fixpoint + persist + re-verify) at 3.2 ms.

## Results (v1.1.0 — after the optimization pass)

```
====================================================================================================
BENCHMARK                                                          TIME       UNIT   NOTES
====================================================================================================
Fr mul (Montgomery, 4 limbs)                                      0.023         µs   field core [1000000 iters]
Fr invert (Fermat)                                                7.279         µs   batch-invert backbone [500 iters]
NTT size 4096                                                   814.197         µs   262.14 KB/s/iter [200 iters]
NTT size 8192                                                  1757.341         µs   524.29 KB/s/iter [200 iters]
G1 scalar mul (windowed, 255-bit)                               188.045         µs   Jacobian w=5, batch-normalized table [50 iters]
BLS sign (hash_to_curve G2 + mul)                              1803.817         µs   ETH2 ciphersuite [20 iters]
BLS verify (2 pairings)                                        3415.896         µs   ETH2 ciphersuite [10 iters]
pairing e(g1, g2) (miller + final exp)                         1573.511         µs   fountain chain [10 iters]
blob_to_kzg_commitment (4096 MSM)                              64705.975         µs   Pippenger, 2 threads [20 iters]
compute_kzg_proof                                              22497.368         µs   quotient + MSM [20 iters]
verify_kzg_proof (single)                                       3236.152         µs   2 pairings + subgroup checks [20 iters]
verify_blob_kzg_proof_batch (6 blobs)                          16658.203         µs   3 MSM + 1 pairing [10 iters]
compute_cells (FFT 8192 + BRP)                                  2943.795         µs   EIP-7594 extension [20 iters]
compute_cells_and_kzg_proofs (FK20)                           305412.353         µs   128 cell proofs, Straus tables [5 iters]
verify_cell_kzg_proof_batch (128 cells)                        48891.922         µs   1 pairing + 129 subgroup checks [5 iters]
recover_cells_and_kzg_proofs (from 64)                        280541.940         µs   erasure decode + FK20 [3 iters]
ZODA commit (16x16 grid, full encode)                           2264.793         µs   chirp-z encode + Merkle + FS [5 iters]
ZODA verify row sample (16x16)                                    0.691         µs   O(2k) inner product [100 iters]
ZODA commit (64x64 grid, full encode)                          38798.650         µs   chirp-z encode, 2 threads [5 iters]
ZODA verify row sample (64x64)                                    2.581         µs   O(2k) inner product [100 iters]
ZODA reconstruct from 32 cols (32x32)                           2439.084         µs   shared Lagrange basis [3 iters]
PQ commit [FAST]                                                261.256         µs   MLWE commit; n=64, 4 chunks [10 iters]
PQ open [FAST]                                                  278.348         µs   sigma-protocol + FS; n=64, 4 chunks [5 iters]
PQ verify [FAST]                                                172.785         µs   2 equations; n=64, 4 chunks [10 iters]
PQ commit [L1]                                                  5227.926         µs   MLWE commit; n=1024, 4 chunks [10 iters]
PQ open [L1]                                                   12973.498         µs   sigma-protocol + FS; n=1024, 4 chunks [5 iters]
PQ verify [L1]                                                 12517.499         µs   2 equations; n=1024, 4 chunks [10 iters]
====================================================================================================
host: 2 cores
```

## v1.0.0 → v1.1.0 (the optimization pass)

Same host, same harness, before → after:

| Path | v1.0.0 (µs) | v1.1.0 (µs) | Speedup |
|---|---:|---:|---:|
| blob_to_kzg_commitment (4096 MSM) | 361,861 | 64,706 | **5.6x** |
| compute_kzg_proof | 131,759 | 22,497 | **5.9x** |
| verify_blob_kzg_proof_batch (6) | 49,252 | 16,658 | **3.0x** |
| compute_cells_and_kzg_proofs (FK20) | 1,835,124 | 305,412 | **6.0x** |
| recover_cells_and_kzg_proofs | 1,860,086 | 280,542 | **6.6x** |
| verify_cell_kzg_proof_batch (128) | 84,679 | 48,892 | 1.7x† |
| ZODA commit (16×16) | 12,611 | 2,265 | **5.6x** |
| ZODA commit (64×64) | 254,719 | 38,799 | **6.6x** |
| ZODA reconstruct (32 cols) | 35,943 | 2,439 | **14.7x** |
| PQ commit [L1] | 49,863 | 5,228 | **9.5x** |
| PQ open [L1] | 60,039 | 12,973 | **4.6x** |
| PQ verify [L1] | 37,013 | 12,517 | **3.0x** |
| PQ commit/open/verify [FAST] | 394/442/262 | 261/278/173 | 1.5–1.6x |
| BLS sign / verify | 2,262 / 3,763 | 1,804 / 3,416 | 1.25x / 1.10x |

† `verify_cell_kzg_proof_batch` and `verify_kzg_proof` now perform the
**spec-required subgroup validation** on every parsed point (the v1.0.0
numbers skipped it — a soundness gap found by the official vector
suite). The checks cost ≈ 0.19 ms × 129 points ≈ 24 ms of the cell-batch
time; excluding them the measured batch verification is ~25 ms, a 3.4x
improvement. c-kzg performs the same checks (at blst-asm speed).

## What changed under the hood

1. **Jacobian group engine** (`zoda-bls`) — `dbl-2007-bl` (2M+5S),
   `madd-2007-bl` (7M+4S), `add-2007-bl` (11M+5S) with explicit ±P edge
   handling, Montgomery batch normalization, and windowed scalar
   multiplication on batch-normalized tables. All scalar-mul paths
   (BLS sign/verify, FFT butterflies, MSMs) route through it.
2. **Pippenger MSM rewrite** (`zoda-kzg::msm`) — precomputed window
   digits, classic all-levels bucket reduction with pristine-bucket
   mixed additions, adaptive window schedule (`c ≈ log₂n − log₂log₂n`),
   persistent thread-local bucket arena, identity-point skips, and
   point-chunked thread parallelism.
3. **FK20 with Straus tables** (`zoda-edas`) — the 128×64 phase-1
   row-MSMs run against per-point 6-bit multiple tables built once per
   setup (~49 MB, mirroring c-kzg's fixed-base precompute option) with
   rows spread over cores; both G1-FFTs run through the Jacobian engine
   with unit-twiddle and identity skips.
4. **Negacyclic NTT ring arithmetic** (`zoda-pq`) — Rq multiplication
   twists by ψⁱ and runs a cyclic size-n NTT (q−1 has 2-adicity 13);
   a density-aware dispatcher keeps sparse operands (narrow challenges,
   monomial evaluation points) on the schoolbook path, where three
   transforms would cost more than they save.
5. **Chirp-z systematic encoder** (`zoda-core`) — interpolation on the
   geometric node set `{ωⁱ}` via the transposed-Vandermonde identity
   (cached barycentric weights, cached `Q(t) = Π(1−ω^j t)`, per-vector
   power sums), replacing the per-vector subproduct tree; row/column
   encoding is thread-parallel.
6. **Shared-Lagrange reconstruction** (`zoda-core`) — the recovery path
   builds the Lagrange basis once for the shared erasure nodes
   (synthetic division + one batch inversion) and reduces every row to
   an O(k²) weighted sum.
7. **Fast canonical parsing** (`zoda-kzg`) — blob field elements parse
   with direct limb extraction + a three-limb comparison instead of
   reduce-then-re-encode.

## Remaining optimization headroom

1. **Pairing stack** — 1.57 ms per pairing vs ~0.4 ms for blst with ADX
   assembly; the Miller loop and final exponentiation are the entire
   cost of single-proof verification. A lazy-reduction Fp2/Fp6 tower is
   the identified path.
2. **GLV scalar multiplication** — the G1 endomorphism
   `φ(x, y) = (βx, y)` halves the ~256 doublings per butterfly
   multiplication in the FK20 G1-FFTs and the subgroup checks (each
   `[r]P` currently pays a full-width scalar multiplication).
3. **Field layer** — dedicated Montgomery squaring and lazy reduction
   in the 6-limb CIOS core (~10–20% on every curve formula), and ADX
   assembly paths where the target allows.
4. **Bluestein upgrade** of the chirp-z encoder (O(n log n) power sums
   and product) for grid sides beyond n = 256, where the current O(n²)
   per-vector form starts to bite.
