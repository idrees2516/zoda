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

**Systematic RS on the roots domain.** Data occupies the first k of 2k
domain points: interpolate, then evaluate (one NTT). The interpolation
nodes {ωⁱ} are a **geometric sequence**, so the transposed-Vandermonde
identity applies: with cached barycentric weights wᵢ = 1/Z'(ωⁱ) and the
cached Q(t) = Πⱼ(1 − ω^j t), each vector's monomial coefficients are
cₖ = Σ_a Q_a·β_{n−1−k−a} from the u-weighted power sums
βₘ = Σᵢ uᵢω^{im} (u = w ⊙ data) — no subproduct tree, no per-vector
allocations, ~10x faster at k = 64. The row and column passes run in
parallel across cores.

**Tensor encoding.** Extend all rows (m×k → m×2k), then all columns
(→ 2m×2k). By linearity every row of the result is a column-code
codeword and every column a row-code codeword.

**Zero-overhead sampling.** Let Z be the tensor codeword and g_r, g_r2
random vectors (Fiat–Shamir from the Merkle roots). The prover publishes
z_r[r] = Σⱼ Z[r][j]·g_r[j] and z_r2[c] = Σᵣ Z[r][c]·g_r2[r]. A sampler
holding the full row W of index r checks W·g_r = z_r[r]; a column
analogously. A forged row can only pass if it matches the projection of
the committed row — probability 1/|F| per forgery attempt after the
commitments fix g_r.

**Reconstruction.** Any k of the 2k evaluations of a degree-<k polynomial
determine it. All rows share the same erasure nodes, so the **Lagrange
basis is built once** (synthetic division of the shared vanishing
polynomial + one batch inversion) and every row reduces to an O(k²)
weighted sum of the basis before a single re-evaluation NTT.

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

## Sampling statistics (zoda-das / zoda-rda)

Distinct (without-replacement) draws per session: the probability that a
specific hidden row of n escapes k distinct samples is exactly
C(n−1,k)/C(n,k) = Πᵢ (n−1−i)/(n−i) — computed in floating point and
saturated to 0 at k = n. Sampling plans solve for the smallest k meeting
a target miss probability; the adaptive session accumulates confidence
across rounds with early exit.
