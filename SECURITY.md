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
