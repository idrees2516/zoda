# Security Policy

## Threat model

zoda secures data availability statements against:
* **Withholding adversaries** — provers that publish commitments but
  withhold/bloat the data. Detection is probabilistic (sampling) with
  exact, documented miss probabilities per session.
* **Invalid-data adversaries** — provers that publish rows, columns,
  cells or proofs that do not match their commitments. All verification
  paths are cryptographic (tensor-code projection checks, KZG pairing
  checks, lattice σ-protocol equations) with per-forgery success
  probability ≤ 2^−128 (field-sized) except where noted.
* **Sybil sampling** — committees and custody assignments are selected
  by BLS-signed sortition over beacon seeds.

## Assumptions

* SHA-256 behaves as a random oracle (Fiat–Shamir transforms).
* For the pairing layer: BLS12-381 discrete logarithms hold (pre-quantum).
* For the lattice layer: Module-LWE / Module-SIS hold (post-quantum);
  parameter sets are first estimates — see docs/POST_QUANTUM.md before
  relying on them.
* Sampling assumes an honest peer will serve *some* of its assigned data;
  the adaptive session accounts for missing responses in its confidence.

## Known limitations (be honest, be safe)

1. **Not constant-time.** Rejection sampling in zoda-pq and some
   conditional paths are timing-variable. Do not use for
   side-channel-sensitive deployments without an audit.
2. **Lattice parameters are estimates.** Not core-SVP-analyzed.
3. **FK20 cell proofs are not yet byte-diffed** against the official
   consensus-spec-tests archives (end-to-end + naive cross-checks pass;
   the vector harness is on the roadmap).
4. **Toy setups are for tests.** `Setup::from_seed_for_testing` embeds
   its own "toxic waste" — never use it in production; load the mainnet
   ceremony file.

## Reporting

Open a GitHub security advisory (Security → Report a vulnerability) or a
private issue. Please include: component crate, code path, input, and
impact analysis against the threat model above.

## Disclosure

Coordinated disclosure; fixes land before public disclosure for
exploitable issues.

## v1.3 — GLV, batch verification and subgroup checking notes

* **GLV fast paths are public-point only.** `zoda_bls::endomorphism::mul_g1_public` / `mul_g2_public` are variable-time and mathematically valid only on the r-torsion (the endomorphism-eigenvalue identity `φ = [λ]` fails off G1/G2). They are used exclusively after subgroup verification (SRS elements past `validate_kzg_g1`/the batch check, hash-to-curve outputs, checked keys/signatures) and never for secret scalars — signing keeps the plain windowed path. Subgroup checks themselves never use GLV: on a non-r-torsion point the decomposed multiplication does not compute `[r]P`, so the check would be vacuous.
* **Batch subgroup checks and small subgroups.** Random-combination batching (`[r]·ΣcᵢPᵢ == O`) is unsound when the curve cofactor has small factors: a point of order d vanishes from the combination whenever `d | cᵢ` (~1/d per attempt, trivially grindable for d = 3). The G1 batch check therefore uses coefficients `cᵢ = 1 + h₁·kᵢ` (kᵢ 128-bit, transcript-derived), which guarantees every h₁-order component contributes itself; a failing point is detected except with ~2⁻¹²⁸ probability per grinding attempt. The G2 side uses per-point `[r]P` checks because h₂ ≈ 2⁵⁰⁷ would make `1 + h₂k` coefficients costlier than the direct check. This was caught by the official spec vector `verify_cell_kzg_proof_batch_case_invalid_commitment_2` (a small-order commitment), not by review — a good argument for keeping vector suites wired into CI.
* **Batch BLS verification semantics.** `verify_batch` computes `Π e(rᵢ·pkᵢ, H(mᵢ)) · e(−g1, Σ rᵢ·σᵢ)` with per-item scalar multiplications on the **plain** path: untrusted points keep exactly the h-torsion blindness semantics of single verification (h-components pair to 1; the G1-part of a key still constrains the equation). `verify_batch_strict` performs batch subgroup checks first (independent challenges for the check and the pairing product) and then uses the GLV fast paths, which are sound on the now-verified points.
* **Challenge derivation.** All batch randomness is a SHA-256 transcript over the complete inputs (points, messages, lengths), so it is unpredictable before the batch is fixed; grinding a passing bad batch costs ~2¹²⁸ per attempt.

## v1.4 — parallel verification and hash-to-curve rewrite notes

* **Parallel paths are bit-equal, not merely equivalent.** `multi_miller_loop_parallel`, the chunked `verify_batch`/`verify_batch_strict`, `map_to_curve_sswu_g2_batch` and `hash_to_curve_g2_batch` are pinned by differential tests against the serial/single paths (product-of-products associativity for the Miller accumulation; unique SSWU output for the map — the x-branch is determined by square-ness, the y sign by sgn0 equality), and the RFC 9380 J.9.1/J.10.1 vectors still pass bit-exactly through the rewritten paths. The batch-challenge transcript is unchanged; thread count does not enter verification decisions.
* **The h2c rewrite preserves the degenerate-input behaviour** of the reference implementation (both-sides non-square `u`, probability < 2⁻³⁸¹ from hash_to_field), returning the same degenerate output as the v1.3 root-of-unity search rather than inventing new behaviour.
* **FK20 caches are process-wide and keyed by setup construction identity.** The pre-v1.4 key (first monomial point) collided across distinct Setups — `[τ⁰] = g1` for every ceremony — a latent cross-setup correctness hazard in multi-setup processes, found by running the 320-vector suite single-threaded; see ALGORITHMS.md.
* **EigenDA mainnet configuration** carries published constants only; the two v1 contract addresses flagged in `config/eigenda/mainnet.toml` must be re-verified against EigenCloud docs before being trusted for value.
