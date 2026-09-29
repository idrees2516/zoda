//! # zoda — Zero-Overhead Data Availability
//!
//! A production-grade Rust implementation of the ZODA tensor-code data
//! availability protocol ([eprint 2025/034](https://eprint.iacr.org/2025/034))
//! with full Ethereum EIP-4844/EIP-7594 integration and a
//! post-quantum lattice-based polynomial commitment variant.
//!
//! ## Quick tour
//!
//! ```
//! use zoda::core::{Matrix, ZodaParams};
//! use zoda::math::{Fr, ZodaRng, PrimeField};
//!
//! // Commit an 8x8 data grid
//! let mut rng = ZodaRng::from_seed(*b"zoda-doc-example-seed-0000000000");
//! let mut data = vec![Fr::zero(); 64];
//! for x in data.iter_mut() { *x = rng.next_fr(false); }
//! let grid = Matrix::from_row_major(8, 8, data);
//! let params = ZodaParams::new(8, 8);
//! let prover = params.commit(&grid);
//! let public = prover.public_params();
//!
//! // A sampled row is its own proof — verify with the public parameters
//! let row = prover.matrix.row(5);
//! assert!(zoda::core::verify_row_sample(&public, 5, &row).is_ok());
//!
//! // A forged row fails
//! let mut bad = row.clone();
//! bad[0] = bad[0] + Fr::ONE;
//! assert!(zoda::core::verify_row_sample(&public, 5, &bad).is_err());
//! ```
//!
//! See the crate docs of each module for the component APIs and
//! `docs/` for architecture deep-dives.

pub mod math {
    pub use zoda_math::*;
}

pub mod bls {
    pub use zoda_bls::*;
}

pub mod kzg {
    pub use zoda_kzg::*;
}

pub mod core {
    pub use zoda_core::*;
}

pub mod pq {
    pub use zoda_pq::*;
}

pub mod das {
    pub use zoda_das::*;
}

pub mod edas {
    pub use zoda_edas::*;
}

pub mod sybils {
    pub use zoda_sybils::*;
}

pub mod rda {
    pub use zoda_rda::*;
}

pub mod archival {
    pub use zoda_archival::*;
}

pub mod bridges {
    pub use zoda_bridges::*;
}

pub mod ethrex {
    pub use zoda_ethrex::*;
}

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
