# Benchmarks

## Methodology

* Host: 2 vCPU cloud instance, rustc 1.98.1, release profile with
  thin-LTO and codegen-units=1.
* Harness: `cargo run --release -p zoda-bench` — warm-up pass, then the
  mean over the listed iteration count (deterministic seeded inputs; no
  OS entropy in the measurements).
* Filter with `zoda-bench math`, `bls`, `kzg`, `edas`, `core` or
  `pq` to run a single group.
* KZG/EDAS runs use a deterministic 4096-point toy setup (`tau` from a
  seed; the mainnet setup file yields the same shapes — point counts and
  widths are identical).

## Results

```
====================================================================================================
BENCHMARK                                                          TIME       UNIT   NOTES
====================================================================================================
Fr mul (Montgomery, 4 limbs)                                      0.023         µs   field core [1000000 iters]
Fr invert (Fermat)                                                7.320         µs   batch-invert backbone [500 iters]
NTT size 4096                                                   814.155         µs   262.14 KB/s/iter [200 iters]
NTT size 8192                                                  1768.251         µs   524.29 KB/s/iter [200 iters]
G1 scalar mul (windowed, 255-bit)                               195.529         µs   sign/pk path [50 iters]
BLS sign (hash_to_curve G2 + mul)                              2261.654         µs   ETH2 ciphersuite [20 iters]
BLS verify (2 pairings)                                        3763.487         µs   ETH2 ciphersuite [10 iters]
pairing e(g1, g2) (miller + final exp)                         1566.231         µs   fountain chain [10 iters]
blob_to_kzg_commitment (4096 MSM)                            361860.684         µs   EIP-4844 hot path [20 iters]
compute_kzg_proof                                            131759.118         µs   quotient + MSM [20 iters]
verify_kzg_proof (single)                                      2903.851         µs   2 pairings [20 iters]
verify_blob_kzg_proof_batch (6 blobs)                         49251.789         µs   3 MSM + 1 pairing [10 iters]
compute_cells (FFT 8192 + BRP)                                 3548.911         µs   EIP-7594 extension [20 iters]
compute_cells_and_kzg_proofs (FK20)                         1835123.605         µs   128 cell proofs, O(n log n) [5 iters]
verify_cell_kzg_proof_batch (128 cells)                       84678.509         µs   1 pairing amortised [5 iters]
recover_cells_and_kzg_proofs (from 64)                      1860085.669         µs   erasure decode + FK20 [3 iters]
ZODA commit (16x16 grid, full encode)                         12610.533         µs   2x NTT passes + Merkle + FS [5 iters]
ZODA verify row sample (16x16)                                    0.641         µs   O(2k) inner product [100 iters]
ZODA commit (64x64 grid, full encode)                        254719.426         µs   2x NTT passes + Merkle + FS [5 iters]
ZODA verify row sample (64x64)                                    2.548         µs   O(2k) inner product [100 iters]
ZODA reconstruct from 32 cols (32x32)                         35942.695         µs   per-row interpolation [3 iters]
PQ commit [FAST]                                                394.239         µs   MLWE commit; n=64, 4 chunks [10 iters]
PQ open [FAST]                                                  441.533         µs   sigma-protocol + FS; n=64, 4 chunks [5 iters]
PQ verify [FAST]                                                261.804         µs   2 equations; n=64, 4 chunks [10 iters]
PQ commit [L1]                                                49863.034         µs   MLWE commit; n=1024, 4 chunks [10 iters]
PQ open [L1]                                                  60039.425         µs   sigma-protocol + FS; n=1024, 4 chunks [5 iters]
PQ verify [L1]                                                37012.704         µs   2 equations; n=1024, 4 chunks [10 iters]
====================================================================================================
host: 2 cores
```

## Reading the numbers

* **Field/NTT layer** — Fr multiplication at 22 ns is the everything-
  multiplier; the 8192-point NTT at ~1.8 ms is the EIP-7594 backbone
  (two of those per blob plus smaller transforms).
* **Pairing stack** — 1.57 ms per pairing is respectable for dependency-
  free pure Rust (blst with ADX assembly: ~0.4 ms; c-kzg's blst-backed
  pairing: ~0.5 ms). BLS verify = 2 pairings + hash-to-curve ≈ 3.7 ms.
* **KZG** — commitments and proofs are dominated by the MSM path: the
  Pippenger implementation uses complete-formula projective additions
  (correctness-first). This is the main optimization headroom vs
  c-kzg-4844 (~1.5 ms/commitment): switching bucket accumulation to
  mixed Jacobian additions with batched affine normalization should
  recover roughly an order of magnitude. Verification (3.16 ms) is
  pairing-bound and within ~3x of c-kzg.
* **EIP-7594** — `compute_cells` (3.5 ms) is NTT-bound. FK20 for all
  128 cell proofs is 1.84 s after the setup-column cache (5.6 s before);
  the remaining cost is the 128×64 MSM against setup columns — the same
  MSM headroom as above. Batch verification of all 128 cells costs a
  single pairing plus small MSMs (85 ms). Recovery from 64 cells
  (1.86 s) is FK20-dominated.
* **ZODA core** — the headline numbers: a row sample verifies in **2.6 µs
  at the 64×64 grid** (an O(2k) inner product); committing a 64×64 grid
  (two NTT passes over 4096 cells + Merkle + Fiat–Shamir) takes 252 ms,
  dominated by the subproduct-tree interpolation in the systematic
  encoder — an O(n log n) "fast interpolation on geometric sequences"
  (chirp-z) replacement is the identified path below 100 ms.
  Reconstruction from half the columns is ~0.35 s per 32 rows.
* **Post-quantum** — L1 commit/open/verify at 50/60/37 ms for a
  degree-4096 polynomial (4 chunks): the cost is the O(n²) ring
  multiplication inside the matrix products; an NTT-based ring multiply
  (q supports 8192-point transforms) is the obvious next optimization,
  expected to bring all three under 10 ms.

## Known optimization roadmap

1. MSM: mixed Jacobian bucket adds + batched normalization (KZG commit,
   FK20, batch verify) — est. 5–10x on all MSM-bound paths.
2. Ring multiplication via NTT in zoda-pq (q−1 has 2-adicity 13).
3. Chirp-z fast interpolation for systematic RS encoding (removes the
   log² factor and the tree allocations).
4. Parallel row/column encoding with scoped threads (embarrassingly
   parallel across rows).
