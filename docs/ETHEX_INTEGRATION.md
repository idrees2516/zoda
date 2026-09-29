# ethrex Integration

[ethrex](https://github.com/lambdaclass/ethrex) is LambdaClass's Ethereum
execution client. `zoda-ethrex` provides the bridge types; this document
shows where they plug in.

## Type mapping

| ethrex type | zoda-ethrex type | notes |
|---|---|---|
| `BlobTransactionSidecar` | `BlobTransactionSidecar` | same field layout |
| `KZGCommitment` ([u8;48]) | `KzgCommitment` | byte-identical |
| `KZGProof` ([u8;48]) | `KzgProof` | byte-identical |
| `VersionedHash` (H256) | `VersionedHash` | byte-identical |
| `Blob` (Bytes / [u8;131072]) | `Blob` | byte-identical |

## Block import path

ethrex validates blob transactions when applying a block body:

```rust
// ethrex: crate/ethrex/core/types/blobs_bundle.rs (validation) —
// replace the c-kzg/kzg-rs call with:
fn validate_blob_transactions(sidecar: &BlobTransactionSidecar) -> Result<(), Error> {
    let setup = zoda_kzg::srs::Setup::load_file(trusted_setup_path())?;
    if sidecar.validate(&setup) { Ok(()) } else { Err(Error::InvalidBlobProof) }
}
```

`validate` checks, per blob: the versioned hash matches
`kzg_to_versioned_hash(commitment)`, and the blob's KZG proof verifies
against the commitment. `validate_batch` performs the whole sidecar with
3 MSMs + a single pairing.

## DA-side integration (PeerDAS)

When ethrex grows data-column support, the block builder path is:

```rust
let sidecar = BlobTransactionSidecar::from_blobs(&blobs, &setup)?;   // commitments + proofs
let columns = sidecar.to_data_column_sidecars(&setup)?;              // 128 cells per blob
```

and the node-side path (custody + sampling):

```rust
use zoda_edas::{verify_cell_kzg_proof_batch, recover_cells_and_kzg_proofs};

// verify gossip column sidecars (any number of cells, one pairing)
verify_cell_kzg_proof_batch(&commitments, &cell_indices, &cells, &proofs, &setup)?;

// once 64+ cells of a blob are present, reconstruct everything
let (cells, proofs) = recover_cells_and_kzg_proofs(&known_idx, &known_cells, &setup)?;
```

Custody scheduling follows the spec helpers:

```rust
let groups = zoda_edas::get_custody_groups(node_id, custody_group_count);
let columns = groups.iter().flat_map(|g| zoda_edas::compute_columns_for_custody_group(*g));
```

## Trusted setup

Point both clients at the same mainnet `trusted_setup.txt`
(c-kzg-4844 `src/trusted_setup.txt`). zoda accepts the 3-section format
directly and derives what it needs.

## Performance notes

Measured on 2 vCPU (see BENCHMARKS.md): a single blob commitment is
~360 ms with the pure-Rust MSM (vs ~1.5 ms for c-kzg's assembly-backed
path). If ethrex needs mainnet-grade blob throughput today, use zoda for
the protocol/sampling layers and c-kzg/blst for the raw KZG hot path —
the types convert trivially since everything is byte-oriented. The
MSM/Pippenger path is the documented optimization target for closing the
gap in pure Rust (windowed affine additions with batch normalization).
