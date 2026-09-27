# Ethereum Compliance

zoda targets the Deneb (EIP-4844) and Fulu (EIP-7594 / PeerDAS)
specifications for everything that touches consensus.

## EIP-4844 — blob KZG (zoda-kzg)

Implemented per the Deneb `polynomial-commitments.md` specification:

| spec function | status |
|---|---|
| `hash_to_bls_field` (SHA-256 → int mod r) | ✔ |
| `bytes_to_bls_field` (canonical, < r) | ✔ rejects ≥ r |
| `bls_field_to_bytes` (big-endian) | ✔ |
| `validate_kzg_g1` / `bytes_to_kzg_{commitment,proof}` | ✔ |
| `blob_to_polynomial` | ✔ |
| `compute_challenge` (`FSBLOBVERIFY_V1_`) | ✔ |
| `evaluate_polynomial_in_evaluation_form` (barycentric, BRP domain) | ✔ |
| `compute_quotient_eval_within_domain` | ✔ |
| `blob_to_kzg_commitment` | ✔ |
| `compute_kzg_proof` / `verify_kzg_proof(_impl)` | ✔ |
| `verify_kzg_proof_batch` (`RCKZGBATCH___V1_`) | ✔ |
| `compute_blob_kzg_proof` / `verify_blob_kzg_proof(_batch)` | ✔ |
| `kzg_to_versioned_hash` (0x01 version byte) | ✔ |

Constants: FIELD_ELEMENTS_PER_BLOB = 4096, BYTES_PER_FIELD_ELEMENT = 32,
KZG_ENDIANNESS = big, PRIMITIVE_ROOT_OF_UNITY = 7, KZG_SETUP_G2_LENGTH =
65.

**Bit-level conventions** match c-kzg-4844: bit-reversed
roots-of-unity/Lagrange setup, big-endian field encodings, 48-byte
compressed G1 / 96-byte compressed G2 (ZCash serialization).

### Trusted setup

The parser accepts:
* the current 3-section `trusted_setup.txt` (G1-Lagrange, G2-monomial,
  G1-monomial — the EIP-7594 extension), and
* the legacy 2-section format, deriving G1-monomial via a G1 inverse FFT.

On load it enforces the spec sanity check — if e(G1[1], G2[0]) =
e(G1[0], G2[1]) the G1 points are in monomial form and the setup is
rejected — and validates every point on its curve. The mainnet
`trusted_setup.txt` from c-kzg-4844 parses successfully in the test suite
(the test is skipped when the reference file is not present).

## EIP-7594 / Fulu — cells and PeerDAS (zoda-edas)

Per `polynomial-commitments-sampling.md` and c-kzg-4844's eip7594 sources:

| function | status |
|---|---|
| `compute_cells` (DAS extension: IFFT-4096 → FFT-8192 → BRP → 128 cells) | ✔ |
| `compute_cells_and_kzg_proofs` (FK20, O(n log n)) | ✔ |
| `compute_verify_cell_kzg_proof_batch_challenge` (`RCKZGCBATCH__V1_`) | ✔ |
| `verify_cell_kzg_proof_batch` (single pairing) | ✔ |
| `recover_cells_and_kzg_proofs` (≥ 64 of 128 cells) | ✔ |
| `get_custody_groups` / `compute_columns_for_custody_group` | ✔ |

Constants: FIELD_ELEMENTS_PER_EXT_BLOB = 8192,
FIELD_ELEMENTS_PER_CELL = 64, CELLS_PER_EXT_BLOB = 128,
NUMBER_OF_COLUMNS = 128, SAMPLES_PER_SLOT = 8,
NUMBER_OF_CUSTODY_GROUPS = 128, CUSTODY_REQUIREMENT = 4.

**Validation status.** The full pipeline is exercised end-to-end in the
test suite (128 cells + FK20 proofs, batch verification of all cells, a
spread 16-cell sample, tamper rejection, and recovery from exactly 64
cells with re-verification of the recovered proofs). A single FK20 cell
proof is additionally cross-checked against a naively-computed quotient
commitment. What is *not* yet done: byte-for-byte diffing against the
official consensus-spec-tests KZG vectors (the vectors are distributed as
large archives; the harness can consume them once downloaded —
see "Next steps").

## BLS12-381 and RFC 9380 (zoda-bls)

The hash-to-curve suites `BLS12381G{1,2}_XMD:SHA-256_SSWU_RO_` are
validated against the RFC's own test vectors (Appendix J.9/J.10,
including u[0]/u[1] field elements and full hash_to_curve points), and
`expand_message_xmd` against Appendix K.1. The BLS signature scheme is
the IETF ciphersuite `BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_`
(Ethereum's flavour: G2 signatures, G1 public keys), with on-curve
validation plus optional r-subgroup checks in `verify_strict`.

## ethrex compatibility (zoda-ethrex)

`BlobTransactionSidecar` mirrors ethrex's type of the same name: 48-byte
commitments/proofs, 32-byte versioned hashes, 131072-byte blobs. The
`validate` method performs exactly the checks ethrex runs when importing
a blob transaction (versioned-hash consistency +
`verify_blob_kzg_proof` per blob), and `validate_batch` the batched
equivalent. See `ETHEX_INTEGRATION.md` for wiring instructions.

## Next steps (compliance roadmap)

1. **Official vector harness**: consume the consensus-spec-tests
   `kzg_mainnet` archives for `blob_to_kzg_commitment`,
   `compute_(blob_)kzg_proof`, `verify_*` and the EIP-7594 cell suites;
   wire as an optional test (skipped unless the archive is present).
2. **c-kzg FFI conformance runner**: differential test against the
   reference C library for randomized inputs.
3. **SSZ encodings** for sidecar containers (currently types only).
