# Post-Quantum Polynomial Commitments (zoda-pq)

## Why a lattice variant

KZG's binding rests on the discrete-logarithm hardness of pairings —
breakable by a sufficiently large quantum computer (Shor). ZODA's tensor
code is already post-quantum-friendly (its proofs are information-
theoretic inner-product checks), but the *blob commitment* layer (which
Ethereum binds into consensus) is not. `zoda-pq` provides a drop-in
polynomial commitment with the same commit/open/verify shape whose
assumptions — Module-LWE and Module-SIS — are believed quantum-resistant,
and which needs **no trusted setup**.

## Construction

### Ring and notation

R_q = Z_q[X]/(Xⁿ+1) with q = 8380417 (the Dilithium modulus; 2-adicity 13
supports the negacyclic structure). Vectors of ring elements are written
bold; ‖·‖_∞ is the max-coefficient norm.

### Commit (BDLOP)

Public matrices A ∈ R_q^{k̄×m̄} and B ∈ R_q^{k̄×k} are expanded from a
domain-separated seed (SHA-256 counter mode with rejection sampling for
uniformity). For a message m = (m₁..m_k) ∈ R_q^k — the coefficient chunks
of the committed polynomial — and short randomness r (coefficients from
CBD_η):

    C = A·r + B·m ∈ R_q^{k̄}

* **Hiding**: C is statistically close to uniform given r under the
  Module-LWE assumption (a fresh r per commitment).
* **Binding**: finding (r, m) ≠ (r', m') with equal C is a Module-SIS
  instance.
* No trusted setup: A, B are expandable from public randomness.

### Open (Lyubashevsky σ-protocol, non-interactive via Fiat–Shamir)

To prove that the committed polynomial evaluates to v at an admissible
ring point ζ (see below):

1. Sample masking m₁ ← R_q^k uniform, r₁ short; publish
   C₁ = A·r₁ + B·m₁ and v₁ = L(m₁) where L is the evaluation functional.
2. Derive the challenge x by Fiat–Shamir over the transcript
   (C, C₁, ζ, v): a **narrow ring element** — coefficients in {−1,0,1},
   weight ≤ 32 (challenge space ≈ 2^160, rejection-free to sample).
3. Respond with r₂ = r₁ + x·r (rejection-sampled until ‖r₂‖_∞ ≤ β) and
   m₂ = m₁ + x·m.
4. The verifier checks
   * **Equation 1**: A·r₂ + B·m₂ = C₁ + x·C  (commitment homomorphism)
   * **Equation 2**: L(m₂) = v₁ + x(ζ)·v  (evaluation linearity)
   * **Shortness**: ‖r₂‖_∞ ≤ β

### The evaluation functional and admissible ring points

L is R-linear only at ring points ζ with ζⁿ = −1: the canonical family is
ζ_t = t·X with t an n-th root of unity mod q (since (tX)ⁿ = tⁿXⁿ = −1).
For such ζ, substitution X ↦ ζ is a ring homomorphism, which is exactly
what Equation 2 relies on: the challenge acts through the *character*
x ↦ x(ζ) (L(x·m) = x(ζ)·L(m)), so the verifier evaluates the challenge
at ζ before multiplying v. A scalar evaluation point z maps to the ring
point via the same machinery; the committed polynomial's value at ζ is a
full ring element (the honest cost of lattice-based commitments: values
are n coefficients, not one field element).

### Soundness sketch

Forking yields two accepting transcripts (r₂, m₂), (r₂', m₂') for
challenges x ≠ x'. Their difference gives a Module-SIS solution for the
statement (A, B, C) together with the evaluation relation — the standard
Lyubashevsky argument. The challenge space size (≈ 2^160) sets the
soundness error; the Fiat–Shamir transform makes it non-interactive in
the random-oracle model (SHA-256).

## Parameter sets

| level | n    | m̄ | k̄ | η | β  | notes |
|-------|------|----|----|---|----|-------|
| FAST  | 64   | 3  | 2  | 2 | 128 | tests/benchmarks only |
| L1    | 1024 | 4  | 2  | 2 | 128 | ~128-bit MSIS/MLWE target |
| L3    | 1024 | 5  | 2  | 3 | 160 | ~192-bit target |
| L5    | 1024 | 6  | 3  | 4 | 192 | ~256-bit target |

These are **first-estimate** configurations in the BDLOP/Lyubashevsky
tradition (Dilithium-style modulus, CBD noise, weight-32 challenges).
A rigorous core-SVP analysis for the exact (n, m̄, k̄, ω, β) combinations
— and the corresponding proof sizes — is future work; treat them as
engineering starting points, not security claims.

## Sizes (L1, degree-4096 polynomial, n = 1024, 4 chunks)

| object | size |
|---|---|
| commitment C (k̄=2 ring elements) | 2·1024·23/8 ≈ **5.9 KB** |
| proof (m₂: 4 + r₂: 4 ring elements, C₁: 2, v₁: 1) | ≈ **32 KB** |
| reference (pairing KZG) | 48 B commitment / 48 B proof |

This is the honest lattice tax at naive parameters: the proof carries the
masked message. Production lattice PCS designs (LaBRADOR, Greyhound-style
inner-product folding) shrink proofs to a few KB at these degrees —
implementing an amortized inner-product argument on top of the BDLOP
commitment is the natural next step (the current API was structured to
allow it: `LatticePcs::open` already produces a single σ-protocol
transcript over a random linear combination of chunk relations).

## What is NOT claimed

* The parameters have not been independently analyzed against the best
  known lattice attacks (BKZ core-SVP); do not deploy at NIST levels
  without that analysis.
* The implementation has not been audited and contains timing-variable
  operations (rejection sampling restarts are visible; the API is not
  constant-time).
* Strong Fiat–Shamir binds (C, C₁, ζ, v) but not the full proof context;
  strengthening the transcript is a one-line change held off until the
  parameter analysis is done.

## Testing

The test suite validates: honest round-trips at every level, wrong-point
and wrong-value rejection, tampered-proof rejection, hiding (same
message, different randomness → different commitments) and determinism
(same randomness stream → identical commitments), ring-homomorphism of
the evaluation at admissible points, and chunking for degree ≥ n.
