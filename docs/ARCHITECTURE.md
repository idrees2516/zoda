# Architecture

## Overview

zoda is a layered system: a math foundation, a cryptography layer, the ZODA
protocol core, Ethereum compatibility shims, and service layers built on
top. Dependencies point strictly downward.

```
            ┌────────────────────────────────────────────────┐
            │                services                        │
            │  zoda-das  zoda-rda  zoda-sybils  zoda-archival│
            │                    zoda-bridges                │
            └───────┬──────────┬───────────┬─────────┬───────┘
                    │          │           │         │
            ┌───────▼──────────▼──────┐ ┌──▼─────────▼──────┐
            │   protocol cores        │ │ Ethereum layer    │
            │  zoda-core (tensor DA)  │ │ zoda-kzg (4844)   │
            │  zoda-pq   (lattice)    │ │ zoda-edas (7594)  │
            │  zoda-das (sampling)    │ │ zoda-ethrex       │
            └───────┬─────────────────┘ └──┬────────────────┘
                    │                      │
            ┌───────▼──────────────────────▼──────┐
            │            zoda-bls                 │
            │  BLS12-381 tower, G1/G2, pairing,   │
            │  RFC 9380 hash-to-curve, BLS sigs   │
            └───────────────────┬─────────────────┘
                                │
            ┌───────────────────▼─────────────────┐
            │            zoda-math                │
            │  fields, NTT, poly trees, hashes,   │
            │  Merkle, RNG, statistics            │
            └─────────────────────────────────────┘
```

## Crate responsibilities

### zoda-math
The arithmetic foundation. `mont_field!` code-generates Montgomery-form
prime fields (CIOS multiplication, Fermat inversion in the Montgomery
domain) — both the 4-limb BLS12-381 scalar field `Fr` and the 6-limb base
field `Fp` come from this single audited macro. On top: Goldilocks and the
lattice modulus `Fq` (q = 8380417), iterative per-stage-twiddle NTTs
(natural order in/out, matching c-kzg conventions), subproduct-tree
multipoint evaluation and interpolation, Montgomery batch inversion,
SHA-256/Keccak-256, SHA-256 Merkle trees with domain separation, a
deterministic CTR-mode RNG, and log-space binomial statistics used by the
sampling planners.

### zoda-bls
A dependency-free BLS12-381 library. The Fp2/Fp6/Fp12 tower uses the
canonical non-residues (ξ = 1+u, w² = v). Point arithmetic uses the
complete Renes–Costello–Batina formulas (no edge cases). The pairing is
the optimal ate: Beuchat et al. Miller loop (algorithms 26/27) with the
Fuentes–Castañeda–Rodríguez-Henríquez "fountain" final exponentiation.
**Frobenius coefficients are derived and self-verified at initialization**
from the tower definition (Γ₁[k] = ξ^((p^k−1)/3), γ₁ = w^(p−1), composed
by Fp6-Frobenius), eliminating an entire class of transcription bugs —
the derivation asserts γ_k² = Γ₁[k] for every k on first use.
Hash-to-curve implements RFC 9380 suites for G1 and G2 (expand_message_xmd,
SSWU + isogeny maps) and is validated against the RFC's own test vectors,
including the expand_message_xmd vectors of Appendix K.

### zoda-kzg
The Deneb-specification evaluation-form KZG: `blob_to_kzg_commitment`,
`compute_kzg_proof`, `verify_kzg_proof`(+batch), blob proofs, versioned
hashes. The trusted setup parser accepts both the 2-section and current
3-section c-kzg-4844 `trusted_setup.txt` formats, enforces the spec's
Lagrange-form sanity check (a pairing check), and derives the monomial
form via a G1 inverse FFT when the third section is absent. MSMs use
Pippenger with a thread-parallel path for large batches.

### zoda-core
The ZODA protocol itself. Encoding is systematic Reed–Solomon over
roots-of-unity domains: `rs_encode_vector` runs **O(n log n)** — the
transposed-Vandermonde identity's two quadratic kernels both collapse
into butterfly transforms (the power sums are a single size-2n NTT; the
coefficient correlation is an exact size-2n cyclic convolution with a
cached NTT(Q)), then one evaluation NTT. The tensor encoding is two
passes (extend rows, then columns — the column pass transposed and
blocked for cache locality on large grids). Commitments are Merkle
trees over rows and columns; the random projections `g_r`, `g_r2` (and
their encodings `z_r`, `z_r2`) are Fiat–Shamir-derived from the two
roots, binding every sample check to the committed data. Sampling
verification is a single O(width) inner product per sample.
Reconstruction handles both the fast full-column path (shared Lagrange
basis) and `reconstruct_2d` — the row/column **fixpoint decoder for
arbitrary scattered cell loss**, with per-line consistency checks
against present cells and a final projection verification of the
recovered grid, plus fetch guidance when a decode stalls.

### zoda-pq
The post-quantum polynomial commitment: BDLOP Module-LWE commitments plus
Lyubashevsky-style σ-protocol openings with narrow ring challenges and
Fiat–Shamir. See `POST_QUANTUM.md` for the construction and analysis.

### zoda-edas
EIP-7594: the DAS extension (evaluate the degree-<4096 blob polynomial on
the 8192-domain, bit-reverse, slice into 128 cells of 64 field elements),
FK20 cell proofs (two phases: circulant-FFT Toeplitz columns + MSMs, then
an unscaled inverse G1-FFT and a final G1-FFT — O(n log n) for all 128
proofs instead of 128 separate quotient MSMs), the batch verifier
(aggregated interpolation polynomial + one pairing for any number of
cells), erasure recovery from any 64 cells (vanishing polynomial on the
X^64-spread, coset division), and the custody-group scheduling helpers.

### zoda-das
The 2D sampling engines and the availability theory. `availability`
holds the exact minimal-withholding-set analysis (the (m+1)×(k+1)
rectangle for tensor grids, the 65-of-128 cell bound for EIP-7594 blobs)
with log-space hypergeometric escape bounds — no approximations.
`session` runs attested line sessions: Merkle-proof-bound row/column
samples, multi-peer round-robin routing with per-peer scorecards,
custody-first scheduling, adaptive rounds driven by the exact escape
bound, fail-closed verdicts. `cell_das` is the Fulu-flavoured engine
over blob grids: column draws, per-blob FK20 cell proofs, single-pairing
batch verification per round, per-column custody accumulation.

### zoda-rda, zoda-sybils, zoda-bridges
Service layers. `zoda-rda` adapts sample counts to a confidence target
with early exit. `zoda-sybils` provides BLS-signed sortition,
stake-weighted selection and windowed peer scoring. `zoda-bridges`
packages ZODA row-inclusion proofs and KZG-backed message commitments
for light clients.

### zoda-archival
Custody storage with a security-first discipline: verify-on-insert
(projection checks for ZODA columns, batch cell-proof verification for
EIP-7594 columns), reconstruction that persists and re-verifies its
result, exact byte accounting, and slot-ordered pruning to a byte
budget. Two compact, checksummed, self-describing file formats
(`ZODACST1` grids, `ZODACEL1` cell columns) written atomically by
`DiskCustody`; the public verification parameters are embedded in every
file so loaded custody is re-verifiable standalone.

### zoda-ethrex
Byte-compatible `BlobTransactionSidecar` with the validation ethrex runs
on the import path (`verify_blob_kzg_proof` per blob, versioned-hash
consistency), plus batch validation and the EIP-7594 cell extension used
by data-column sidecars.

## Key data flows

### Blob → columns (block building)
```
blob ─ blob_to_kzg_commitment (EIP-4844)         ──► commitment per blob
    └► blob_to_monomial (IFFT) ─ FFT_8192 ─ BRP ─► 128 cells
                     └► FK20 ─ BRP ────────────► 128 cell proofs
```

### Sampling (block validation, light node)
```
plan (exact rectangle-escape bound) ─► draws (custody columns first)
   ZODA lines: projection check + Merkle inclusion ─► per-peer scorecards
   blob cells: fetch (blob, column) + FK20 proof ─► one batch pairing/round
                                                        ─► verdict + confidence
```

### Reconstruction (custody node / archival)
```
≥ k of 2k columns ─► per-row interpolation ─► full tensor codeword
scattered cells   ─► row/col fixpoint + projection verify ─► full grid
                     (stall ─► guidance: which cells to fetch next)
≥ 64 of 128 cells ─► vanishing-poly erasure decode ─► 128 cells + proofs
recovered custody ─► re-verify ─► persist (atomic, checksummed)
```

## Design decisions worth knowing

* **Zero external dependencies** for the entire workspace. Every line of
  field arithmetic, hashing and curve logic is in-tree and testable.
* **Self-deriving constants**: Frobenius coefficients and the KZG
  Lagrange↔monomial conversions are computed from definitions and
  verified at startup, not transcribed.
* **Correctness-first point arithmetic**: complete formulas everywhere;
  the MSM path trades some speed for never having an edge case.
* **Fiat–Shamir everywhere randomness enters**: the ZODA projections, the
  KZG batch challenges, the lattice σ-protocol challenges.
