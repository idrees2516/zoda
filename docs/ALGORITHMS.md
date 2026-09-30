# Algorithms

This is the math inside every hot path, with references.

## Field arithmetic (zoda-math)

**Montgomery multiplication (CIOS).** All field elements are kept in
Montgomery form x̄ = xR mod p; multiplication is the coarsely-integrated
operand-scanning (CIOS) algorithm with a final conditional subtraction —
fully reduced results, no data-dependent branches. `Fermat inversion`
computes a^(p−2) directly in the Montgomery domain (squaring chain over
the exponent limbs), used by batch inversion (one inversion + 3(n−1)
multiplications for n elements).

**NTT.** Iterative radix-2 Cooley–Tukey with per-stage contiguous twiddle
tables (sequential access in every butterfly pass), a precomputed
bit-reversal permutation, and coset variants (`coset_fft` /
`coset_ifft`) used by EIP-7594 recovery. The convention matches c-kzg-4844
exactly: `fft(x)[i] = Σⱼ x[j]·ω^{ij}` with ω the primitive
FIELD_ELEMENTS_PER_EXT_BLOB-th root — so field data interoperate with
consensus structures bit-for-bit.

**Subproduct trees.** Multipoint evaluation and interpolation over
arbitrary (power-of-two sized) point sets in O(n log² n): build the
product tree of (X − xᵢ), evaluate down with polynomial remainders,
interpolate up with the Lagrange-merge `P_{S∪T} = P_S·Z_T + P_T·Z_S`.

## BLS12-381 (zoda-bls)

**Tower.** Fp2 = Fp[u]/(u²+1), Fp6 = Fp2[v]/(v³−ξ) with ξ = 1+u,
Fp12 = Fp6[w]/(w²−v). Karatsuba throughout; the Fp12 cyclotomic square is
the Granger–Scott fp4-square decomposition.

**Point arithmetic.** Two coordinate systems, two roles:

* **Homogeneous + RCB** (complete Renes–Costello–Batina formulas,
  eprint 2015/1060, algorithms 7/8/9): add(P,P), add(P,−P) and the
  identity need no special cases — used for on-curve validation and
  wherever inputs are untrusted multiples.
* **Jacobian** (the EFD "2007-bl" family): doubling 2M+5S
  (`dbl-2007-bl`), mixed addition 7M+4S (`madd-2007-bl`), general
  addition 11M+5S (`add-2007-bl`) — the hot-path engine behind every
  MSM, FFT butterfly and scalar multiplication. The ±P degeneracies
  (`H = 0`) are detected explicitly and routed to doubling / identity,
  and Montgomery's trick batch-normalizes whole point arrays with a
  single field inversion. Windowed (w = 5) scalar multiplication runs
  on a projectively-built, batch-normalized table.

**Optimal ate pairing.** Miller loop over the 64-bit BLS parameter
x = −0xd201000000010000 using Beuchat–López-Tehraní doubled-line
coefficients (eprint 2010/354, algorithms 26/27), line evaluations applied
through the sparse `mul_by_014`, and the Fuentes–Castañeda–Rodríguez-
Henríquez final exponentiation: the easy part (p⁶−1)(p²+1) by Frobenius
and conjugation, the hard part by the "fountain" addition chain — five
cyclotomic exponentiations by |x| plus Frobenius maps.

**Hash-to-curve (RFC 9380).** expand_message_xmd(SHA-256) → hash_to_field
→ Simplified SWU on the isogenous curve E' → 11-isogeny (G1) /
3-isogeny (G2) → cofactor clearing by the effective cofactors h_eff.

## KZG on BLS12-381 (zoda-kzg)

Evaluation-form KZG exactly as the Deneb consensus spec: the blob is the
polynomial's evaluations on the bit-reversed 4096-domain; commitments are
MSMs against the (bit-reversed) Lagrange-form trusted setup; proofs at a
challenge z are the quotient (p(x) − p(z))/(x − z) in evaluation form,
with the special in-domain case handled by the spec's
`compute_quotient_eval_within_domain`. Batch verification folds n
statements into 3 MSMs and one pairing via r-power random linear
combinations.

**MSM.** Pippenger's bucket method with the full modern treatment:
window digits are extracted once per scalar (not per window), bucket
accumulation runs in Jacobian coordinates against freshly-copied
affine buckets (Z = 1), the classic all-levels running-prefix reduction
computes Σ_j j·B_j with pristine buckets taking the cheaper mixed-add
path, the window schedule adapts as c ≈ log₂n − log₂log₂n + 1
(clamped to [4, 13]), a persistent thread-local bucket arena is reset
only at touched entries, identity points are skipped, and large
batches are chunked across threads (each chunk runs the full window
schedule).

**G1-FFT.** Cooley–Tukey over Jacobian points for point-polynomials
P(X) = Σ pᵢXⁱ: (n/2)·log₂n twiddle scalar-multiplications with
unit-twiddle and identity skips, natural-order output, batch
normalization — the FK20 transform backbone.

## ZODA tensor code (zoda-core)

**Systematic RS on the roots domain — O(n log n).** Data occupies the
first k of 2k domain points: interpolate, then evaluate (one NTT). The
interpolation nodes {ωⁱ} are a **geometric sequence**, so the
transposed-Vandermonde identity applies: with cached barycentric weights
wᵢ = 1/Z'(ωⁱ) and the cached Q(t) = Πⱼ(1 − ω^j t), each vector's
monomial coefficients are cₖ = Σ_a Q_a·β_{n−1−k−a} from the u-weighted
power sums βₘ = Σᵢ uᵢω^{im} (u = w ⊙ data). Both quadratic kernels of
that identity collapse into butterfly transforms:

* βₘ = Σᵢ uᵢω^{im} is *literally* the m-th output of the size-2n NTT of
  u — one transform replaces the n² power-sum loop;
* cₖ = R[n−1−k] where R = Q·B is an **exact** product (deg Q + deg B =
  n + (n−1) = 2n−1, so a size-2n cyclic convolution has no wraparound):
  NTT(Q) is cached per domain size, NTT(B) is the second transform, one
  INTT recovers R.

Three size-2n transforms per vector plus the final evaluation NTT —
measured **49x over the schoolbook loops at n = 1024** (59 ms → 1.2 ms),
with the crossover at n ≈ 40 (below it the schoolbook loops win and are
used). The row and column passes run in parallel across cores.

**Tensor encoding.** Extend all rows (m×k → m×2k), then all columns
(→ 2m×2k). By linearity every row of the result is a column-code
codeword and every column a row-code codeword. The column pass runs in
a **transposed layout** on large grids: a blocked 32×32 transpose makes
every column a contiguous row, the encodes stream, and a second
transpose restores the orientation — two linear cache-friendly passes
instead of one cache miss per matrix element.

**Zero-overhead sampling.** Let Z be the tensor codeword and g_r, g_r2
random vectors (Fiat–Shamir from the Merkle roots). The prover publishes
z_r[r] = Σⱼ Z[r][j]·g_r[j] and z_r2[c] = Σᵣ Z[r][c]·g_r2[r]. A sampler
holding the full row W of index r checks W·g_r = z_r[r]; a column
analogously. A forged row can only pass if it matches the projection of
the committed row — probability 1/|F| per forgery attempt after the
commitments fix g_r.

**Reconstruction (full columns).** Any k of the 2k evaluations of a
degree-<k polynomial determine it. All rows share the same erasure
nodes, so the **Lagrange basis is built once** (synthetic division of
the shared vanishing polynomial + one batch inversion) and every row
reduces to an O(k²) weighted sum of the basis before a single
re-evaluation NTT.

**2D erasure decoding (arbitrary cell loss).** `reconstruct_2d` runs the
classical product-code fixpoint: every row with ≥ k of 2k cells present
is interpolated (from its first k cells) and filled; every column with
≥ m of 2m cells likewise; repeat until a full pass recovers nothing.
Rows sharing an availability pattern reuse the cached O(k²) Lagrange
basis (keyed by the presence bitmask). Security: every decoded line is
**checked against all its remaining present cells** (a tampered cell
fails closed with its index), and the completed grid is verified line
by line against the public Fiat–Shamir projections — catching the
consistent-but-wrong decodings (exactly-k adversarial cells) that no
re-encoding check can detect. On a stall the decoder reports the cells
whose fetch would unlock the closest-to-decodable line
(`next_cells_to_fetch`).

## EIP-7594 cells + FK20 (zoda-edas)

**DAS extension.** blob (evaluation form) → monomial coefficients
(IFFT-4096) → zero-pad to 8192 → FFT-8192 → bit-reverse → 128 cells × 64
field elements. Each cell covers a coset {h_c·ω₆₄ʲ} of the
roots-of-unity domain, h_c = ω^rev₇(c).

**FK20 cell proofs** (Fypiński–Krawczyk, eprint 2023/033). For all 128
cells in O(n log n): (1) for each of 64 offsets, the circulant
coefficients of the quotient Toeplitz matrix are FFT'd over the 128-domain
and multiplied component-wise against the FFT'd setup columns
(precomputed and cached), giving the v(X) coefficients via MSMs and an
unscaled inverse G1-FFT; (2) a final forward G1-FFT evaluates v at the
128 cell positions. Proofs are bit-reversed to match cell order.
The 128×64 phase-1 MSMs run against per-point **Straus window tables**
(6-bit multiples of the fixed setup columns, built once per setup and
cached — the same trade c-kzg makes with its fixed-base precompute),
with rows spread over cores; both G1-FFTs use the Jacobian engine with
unit-twiddle and identity skips.

**Batch verification.** Deduplicate commitments, derive the
RCKZGCBATCH__V1_ challenge, form the r-power-weighted sums of proofs
(coset-weighted by h_c⁶⁴), commitments, and the aggregated interpolation
polynomial (per-column IFFT-64 shifted by h_c⁻¹), then check one pairing:
e(Σ(Cᵢ − Rᵢ + hᵢ⁶⁴πᵢ), [1]₂)·e(−Σπᵢ, [s⁶⁴]₂) = 1.

**Recovery.** Given ≥ 64 cells: the vanishing polynomial of the missing
cosets is short_Z(X⁶⁴) (spread by 64 so it vanishes on whole cosets);
multiply the known evaluations by Z, convert to coefficients, divide by Z
on the 7-coset (avoiding zeros), convert back, and re-evaluate.

## Lattice post-quantum commitments (zoda-pq)

BDLOP Module-LWE commitments with Lyubashevsky σ-protocol openings — see
`POST_QUANTUM.md` for the full treatment: the ring R_q = Z_q[X]/(Xⁿ+1),
the commitment equation C = A·r + B·m, narrow ring challenges, rejection
sampling, and the R-linearity subtlety of ring-point evaluation
(L(x·m) = x(ζ)·L(m) for admissible ζ with ζⁿ = −1).

**Ring multiplication** is **density-aware**: dense × dense operands
route through the negacyclic NTT (twist by ψⁱ, cyclic size-n transform,
pointwise product, untwist — q − 1 has 2-adicity 13, so transforms up to
length 8192 exist), while sparse operands (narrow weight-ω challenges,
monomial evaluation points, scalars) stay on the zero-skipping
schoolbook path with the sparser side driving the outer loop — three
NTT passes cost more than they save below ≈ 3.5·log₂n nonzeros.

## 2D availability theory (zoda-das::availability)

**The minimal withholding set.** The 2m×2k tensor code has minimum
distance (m+1)(k+1) (product-code distance = d_row·d_col). Nothing
smaller can be ambiguous. Conversely, minimum-weight RS codewords
vanish at *any* chosen k−1 of the 2k domain points (degree-(k−1)
polynomials with prescribed roots), so their outer products are tensor
codewords whose supports are **(m+1)×(k+1) rectangles** — hiding such a
rectangle always leaves two codewords consistent with everything
available. The adversary's optimal strategy is therefore exactly a
rectangle, and every bound below is exact, not heuristic.

* **Cell sampling**: s uniform distinct cells of the 4mk grid escape the
  rectangle with probability C(4mk−(m+1)(k+1), s)/C(4mk, s)
  (hypergeometric, evaluated in log space). A 128×128 grid needs 95
  cell samples for 2⁻⁴⁰.
* **Line sampling (ZODA)**: hiding the rectangle means refusing m+1
  rows AND k+1 columns; u_r row draws and u_c column draws all miss
  with probability the product of the two hypergeometric tails
  C(2m−(m+1), u_r)/C(2m, u_r) · C(2k−(k+1), u_c)/C(2k, u_c).
* **EIP-7594 blob grids**: a blob is ambiguous only when ≥ 65 of its
  128 cells are hidden (65·64 = 4160 > 4095 evaluations); column
  sampling detects it unless every sampled column avoids the hidden
  ones — C(63, s)/C(128, s) ≈ 2⁻⁸·⁵ at the spec's 8 samples per slot.

`log_binomial` uses an exact product (not Stirling-series differences)
for min(k, n−k) ≤ 1024: the cancellation-free form keeps hypergeometric
*ratios* accurate to ~1e-15, which the confidence machinery depends on.

## Sampling engines (zoda-das)

**Attested line sessions** (`run_attested_session`): adaptive rounds of
row/column draws (custody columns scheduled first, then uniform undrawn
indices), round-robin routed across a peer set with per-peer
scorecards. Every fetched line carries its **Merkle inclusion proof**:
the projection check binds the line to the Fiat–Shamir challenges, the
Merkle proof binds it to the committed root — together they make
cross-grid replay and forgery both fail closed. Verdicts are exact:
`Available` (escape bound ≤ target), `Rejected` (any sample failed),
`Insufficient` (budget exhausted, achieved confidence reported).

**Cell sessions over blob grids** (`sample_cell_session`): the Fulu
2D layout — each blob is a 128-cell row, columns are cell indices
across all blobs. Rounds draw columns (custody-first), fetch every
blob's cell in each drawn column with its FK20 proof, and batch-verify
the whole round with a **single pairing** (`verify_cell_kzg_proof_batch`,
with a bisection fallback that pinpoints the bad cell if the batch
fails). Per-blob and per-column custody counters accumulate for
reconstruction readiness.

## Custody storage (zoda-archival)

Verify-on-insert everywhere (projection check for ZODA columns, batch
cell-proof verification for EIP-7594 columns). Reconstruction
**persists its result**: recovered columns are re-verified then stored,
so custody only ever contains lines bound to the commitment. The
compact file format is self-describing and self-verifying: header
(grid id, shape, the full public verification parameters), an
availability bitfield, raw 32-byte field elements of exactly the
present columns, and a trailing SHA-256. `DiskCustody` writes
atomically (temp + fsync + rename); `CustodyStore` enforces byte
budgets by evicting oldest-slot grids. The EIP-7594 equivalent
(`CellCustody`) stores per-blob (cell, proof) pairs per column with
the same discipline — 2104 bytes per custodied column-blob slot, no
framing overhead.
