# Benchmarks

## Methodology

All numbers: release profile with LTO, 2 vCPU host, `zoda-bench`
(`cargo run --release -p zoda-bench [-- <group>]`; groups: `math`, `bls`,
`glv`, `kzg`, `edas`, `core`, `das2d`, `pq`, `eigenda`,
`throughput`). Timings are wall
clock averages after one warmup call.

## Results (v1.4.0 — h2c SSWU complex method, parallel Miller accumulation, EigenDA mainnet pipeline)

### Hash-to-curve (the verify_batch hotspot)

The G2 Simplified-SWU square root is now the **complex method** — two
windowed Fp exponentiations by the compile-time `(p+1)/4` constant
(norm, then real part) instead of one 762-bit Fp2 exponentiation — and
the division by the g(x) denominator is resolved by the homogeneous
coordinate trick `(x_num·x_den², y', x_den³)`, so the geometry needs no
field inversion at all. The batch path defers every inversion into two
Montgomery batch inversions across the whole batch (Fp for the sqrt
denominators, Fp2 for the sign corrections). The Budroni–Pintore
cofactor chain now runs natively in Jacobian coordinates with
Hamming-weight-aware 64-bit `c1` multiplications (BLS_X has 6 set bits)
instead of the generic W=5 windowed path whose table build + field
inversion cost more than the scalar work at that width.

| benchmark | v1.3.0 | v1.4.0 | change |
|---|---|---|---|
| hash_to_curve G2 (single) | 809 µs | **517 µs** | 1.57x |
| hash_to_curve G2 (batched, per message) | ~620 µs | **378 µs** | 2.14x cumulative |
| map_to_curve_sswu_g2 (the SSWU core) | ~270 µs | **134 µs** | ~2x |
| clear_cofactor_g2 (psi-chain) | 463 µs* | **243 µs** | 1.9x |
| BLS sign | 1,391 µs | **1,030 µs** | 1.35x |

\* component-isolated measurement; the v1.3 h2c total was
sqrt-chain-bound.

### Batch BLS verification (parallel Miller-loop accumulation)

`verify_batch` / `verify_batch_strict` now chunk the per-item work —
batched hash-to-curve, the plain (or GLV) scalar multiplications, G2
preparation and each chunk's share of the Miller-loop accumulation —
across worker threads. The per-chunk Fp12 products combine with `T−1`
multiplications, the `Σ r_i·σ_i` term folds in serially, and a single
final exponentiation closes the check. `multi_miller_loop_parallel`
exposes the same chunked accumulation for any product-pairing check.

| benchmark | v1.3.0 | v1.4.0 | change |
|---|---|---|---|
| BLS verify_batch x16 | 34.7 ms | **15.6 ms** | 2.2x |
| BLS verify_batch_strict x16 | 45.1 ms | **28.9 ms** | 1.6x |
| BLS verify x16 (individual, reference) | 47.1 ms | 41.6 ms | — |

### EigenDA mainnet pipeline throughput (real ceremony SRS)

`zoda-bench throughput` runs the full disperser → operator pipeline
(`zoda_edas::pipeline`) on the shipped `config/eigenda/mainnet.toml`
with the **real Ethereum mainnet KZG ceremony setup** — dispersal
(commitment + FK20 cell proofs), single-pairing batch verification of
every cell, and custody reconstruction from a 50% slice, threaded
across blobs:

| stage (8 MiB batch, 2 vCPU) | throughput |
|---|---|
| commit (blob → KZG, parallel MSM) | 2.03 MB/s |
| extend + FK20 prove | 0.47 MB/s |
| batch verify (8,192 cells, one pairing) | 11.0 MB/s |
| custody reconstruct (50% erasure, threaded) | 0.51 MB/s |
| disperser e2e (commit + prove) | 0.38 MB/s |
| node e2e (verify + recover) | 0.49 MB/s |
| full e2e (commit + prove + verify) | 0.37 MB/s |

**The 100 MB/s question, answered honestly.** On this 2-vCPU host the
full pipeline runs at ~0.4 MB/s and the verification-only path at
~11 MB/s; FK20 proving is the binding stage at ~0.24 MB/s per MiB-scale
batch (0.47 MB/s at 8 MiB as the Straus tables amortise). Per-core
FK20 throughput of ~0.25 MB/s is within ~2–3x of hand-tuned assembly
implementations (blst-class c-kzg FK20 runs ≈0.6–1 MB/s per modern
core); the remaining gap is the portable field layer (51 ns/Fp-mul vs
~15–20 ns for ADX assembly). Reaching 100 MB/s end-to-end therefore
requires *both* an assembly field layer *and* ~50–100 cores of this
class — or, on this host, serving the operator path (verify +
reconstruction) which is not FK20-bound. Two real defects found and
fixed while chasing this target are documented in ALGORITHMS.md: the
thread-local FK20 table caches (every worker thread rebuilt the 49 MB
Straus table; 7x on threaded recovery) and the setup-identity cache
key collision (`[τ⁰] = g1` for every setup — a latent correctness bug
since v1.1, previously masked by thread placement).

## Results (v1.3.0 — GLV, ADX field layer, batch verification, parallel 2D decode)

### Pairing-stack and scalar multiplication

| benchmark | v1.2.0 | v1.3.0 | change |
|---|---|---|---|
| pairing e(g1,g2) (Miller + fountain FE) | 2,014 µs | **1,604 µs** | 1.26x |
| BLS verify (2-term pairing_check) | 4,086 µs | **2,915 µs** | 1.40x |
| BLS sign (h2c G2 + scalar mul) | 2,133 µs | **1,391 µs** | 1.53x |
| hash_to_curve G2 (incl. cofactor) | ~2.1 ms | **809 µs** | ~2.6x (psi-chain) |
| G1 scalar mul (plain 255-bit) | 216 µs | 177 µs | 1.22x (field layer) |
| G1 scalar mul (GLV 2-dim) | — | **133 µs** | 1.63x vs v1.2 plain |
| G2 scalar mul (plain 255-bit) | — | 522 µs | reference |
| G2 scalar mul (GLV 2-dim) | — | **477 µs** | 1.09x |

The pairing stack gains come from the ADX/BMI2 Montgomery paths (runtime
`#[target_feature]` clones emitting `mulx`/`adcx`/`adox`) and lazy-reduction
Karatsuba Fp2 multiplication (three wide products + two REDCs instead of
three fully reduced products).

### Verification paths

| benchmark | v1.2.0 | v1.3.0 | change |
|---|---|---|---|
| G1 subgroup check x129 (individual) | 24.9 ms | — | pre-v1.3 path |
| G1 subgroup check x129 (batched MSM) | — | **5.7 ms** | 4.3x |
| verify_cell_kzg_proof_batch (128 cells) | 48.8 ms | **31.8 ms** | 1.53x |
| verify_blob_kzg_proof_batch (6 blobs) | 15.5 ms | 15.5 ms | flat (MSM-bound) |
| BLS verify x16 (individual) | — | 47.1 ms | reference |
| BLS verify_batch x16 | — | **34.7 ms** | 1.35x |
| BLS verify_batch_strict x16 | — | 45.1 ms | subgroup soundness included |

### 2D reconstruction and the EigenDA-style pipeline

| benchmark | v1.2.0 | v1.3.0 | change |
|---|---|---|---|
| reconstruct_2d 64x64 (60% erased) | 16.0 ms | **10.8 ms** | 1.48x (parallel passes) |
| FK20 compute_cells_and_kzg_proofs | ~310 ms | **280 ms** | 1.11x |
| recover_cells_and_kzg_proofs (64) | ~290 ms | 292 ms | flat |
| blob_to_kzg_commitment (4096 MSM) | 64.7 ms | 64.7 ms | flat (MSM-bound) |
| compute_kzg_proof | "24.5 ms"* | **65.6 ms** | *corrected, see note |

**Measurement correction**: the pre-v1.3 `compute_kzg_proof` bench fed a
little-endian encoding of a possibly-unreduced scalar as `z`; the canonical
big-endian parser rejected roughly half the draws, so the recorded average
mixed ~85 µs fast-path errors with ~65 ms real runs. v1.3 fixes the harness
(the spec vectors always exercised the real path — 320/320 remain green)
and reports the true MSM-bound cost.

### EigenDA-style operator pipeline (new)

One batch = 8 blobs (1 MiB): commit (blob -> KZG), extend into 128 cells +
FK20 proofs each, then a single batched cell verification (1024 cells):

| phase | time | throughput |
|---|---|---|
| commit (blob -> KZG) | 513 ms | 2.04 MB/s |
| extend + FK20 prove | 6,675 ms | 0.16 MB/s |
| batch verify (1024 cells) | 124 ms | **8.45 MB/s** |
| end-to-end | 7,313 ms | 0.14 MB/s |

This is the base case for EigenLayer-level operator duties on a 2-vCPU
host: attestation (batch verify) runs at 8.45 MB/s; proving throughput is
setup-dominated (see headroom below).


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

1. **Pairing stack** — 1.6 ms per pairing vs ~0.4 ms for blst with ADX
   assembly; the Miller loop and final exponentiation are the entire
   cost of single-proof verification (the *multi*-pairing accumulation
   is now parallel — v1.4; a single pairing is inherently serial).
2. **FK20 / MSM field layer** — the portable CIOS core (51 ns/Fp-mul,
   squares not yet specialised) sits ~2.5–3x from ADX assembly; FK20
   proving is the binding stage of the disperser path. An Fp squaring
   macro and a G2 Pippenger (for `Σ rᵢσᵢ` in verify_batch) are the
   identified next steps.
3. **Bluestein upgrade** of the chirp-z encoder (O(n log n) power sums
   and product) for grid sides beyond n = 256, where the current O(n²)
   per-vector form starts to bite.
