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

**Point arithmetic.** Renes–Costello–Batina complete formulas for a = 0
curves (eprint 2015/1060, algorithms 7/8/9) in homogeneous projective
coordinates: add(P,P), add(P,−P) and the identity all work with no
special cases.

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

**MSM.** Pippenger with a size-adaptive window; the parallel path splits
the batch across threads and sums the per-shard results.

## ZODA tensor code (zoda-core)

**Systematic RS on the roots domain.** Data occupies the first k of 2k
domain points: interpolate (subproduct tree) then evaluate (one NTT).
Systematic because the first k evaluations reproduce the input.

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
determine it: interpolate each row from the available positions and
re-evaluate on the full domain.

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

## Sampling statistics (zoda-das / zoda-rda)

Distinct (without-replacement) draws per session: the probability that a
specific hidden row of n escapes k distinct samples is exactly
C(n−1,k)/C(n,k) = Πᵢ (n−1−i)/(n−i) — computed in floating point and
saturated to 0 at k = n. Sampling plans solve for the smallest k meeting
a target miss probability; the adaptive session accumulates confidence
across rounds with early exit.
