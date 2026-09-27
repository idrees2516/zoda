# Deployment

## Building

```bash
# Full workspace (zero external dependencies — no system libraries needed)
cargo build --release

# Tests (includes RFC 9380 + mainnet-setup validation)
cargo test --workspace --release

# Benchmarks
cargo run --release -p zoda-bench
```

Rust 1.75+ (tested on 1.98). No C toolchain, no OpenSSL, no assembly —
the whole stack is portable pure Rust (the only platform-specific call is
reading `/dev/urandom` for production key material, with a deterministic
fallback for locked-down environments).

## Artifacts

| what | how |
|---|---|
| library crate | `zoda` (facade re-exporting the full stack) |
| benchmark binary | `target/release/zoda-bench` |
| Docker image | `docker build -t zoda .` (multi-stage, distroless runtime) |

## Configuring the trusted setup

KZG operations need `trusted_setup.txt`:

```rust
let setup = zoda_kzg::srs::Setup::load_file(std::path::Path::new(
    &std::env::var("ZODA_TRUSTED_SETUP")?.))?;   // mainnet file from c-kzg-4844
```

Deployment checklist:
1. Download the mainnet `trusted_setup.txt` (c-kzg-4844 repository) and
   pin its SHA-256 in your deployment config.
2. Load once at startup (~0.2 s parse + validation incl. the spec's
   Lagrange-form pairing check).
3. Treat the setup as a public constant: distribute it with the binary
   or fetch with integrity verification — never accept it over an
   unauthenticated channel.

## Running a DA node

The service layers compose:

* **Custody**: `zoda_edas::get_custody_groups(node_id, count)` → store
  the assigned columns (`zoda_archival::ColumnStore`), serve them on
  request.
* **Sampling**: `zoda_rda::RdaSession` with a confidence target — run
  per slot against your peer set (`zoda_das::SampleOracle`
  implementations map onto your network stack).
* **Reconstruction**: when ≥ 64 cells of a blob arrive,
  `zoda_edas::recover_cells_and_kzg_proofs` and re-serve the recovered
  columns (cross-seeding per the Fulu spec).
* **Peer management**: `zoda_sybils::PeerScorer` penalties for invalid
  samples; sortition for sampling committees.

A node binary wiring these together against a real P2P stack (libp2p)
is intentionally not in scope yet — the layers are trait-based
(`SampleOracle`, custody stores) so the network binding is a thin
adapter.

## Docker

```bash
docker build -t zoda .
docker run --rm zoda zoda-bench          # quick smoke test + numbers
```

The image is a two-stage build: cargo-chef-cached builder, then a
distroless-style minimal runtime carrying only the benchmark binary and
the trusted setup slot (mount at `/data/trusted_setup.txt`).

## Production notes

* **Key material**: `zoda_bls::SecretKey::from_seed` is for tests;
  production keys should come from your KMS via
  `SecretKey::from_scalar`.
* **Constant time**: the field and curve code is branch-reduced but not
  audited constant-time; assume side-channel hardening is your
  responsibility (same posture as most research-grade Rust crypto).
* **Logging**: none by default (a library); services return structured
  results (`SessionReport`, verdicts) for your telemetry.
* **Resource ceiling**: the largest fixed allocations are the 8192-point
  transforms (~256 KB per NTT) and the trusted setup (~1 MB in memory).
