//! Parameter sets for the lattice polynomial commitment.

/// A parameter configuration.
#[derive(Clone, Copy, Debug)]
pub struct PqParams {
    /// Ring dimension (R_q = Z_q[X]/(X^n+1)).
    pub n: usize,
    /// Width of the MLWE randomness vector r (hiding).
    pub m_bar: usize,
    /// Number of ring elements in the commitment (binding amplification).
    pub k_bar: usize,
    /// CBD parameter for the short randomness.
    pub eta: u32,
    /// Rejection-sampling bound β for the response norm.
    pub beta: u32,
    /// Security level label.
    pub level: &'static str,
}

impl PqParams {
    /// NIST level 1 (~128-bit post-quantum target).
    pub fn nist_l1() -> PqParams {
        PqParams {
            n: 1024,
            m_bar: 4,
            k_bar: 2,
            eta: 2,
            beta: 240,
            level: "L1",
        }
    }

    /// NIST level 3 (~192-bit).
    pub fn nist_l3() -> PqParams {
        PqParams {
            n: 1024,
            m_bar: 5,
            k_bar: 2,
            eta: 3,
            beta: 320,
            level: "L3",
        }
    }

    /// NIST level 5 (~256-bit).
    pub fn nist_l5() -> PqParams {
        PqParams {
            n: 1024,
            m_bar: 6,
            k_bar: 3,
            eta: 4,
            beta: 400,
            level: "L5",
        }
    }

    /// A fast configuration for tests and benchmarks.
    pub fn fast() -> PqParams {
        PqParams {
            n: 64,
            m_bar: 3,
            k_bar: 2,
            eta: 2,
            beta: 80,
            level: "FAST",
        }
    }

    /// Sample a random coefficient-form polynomial of the given degree
    /// (coefficients uniform mod q).
    pub fn random_poly(&self, rng: &mut zoda_math::ZodaRng, degree: usize) -> Vec<zoda_math::Fq> {
        (0..=degree)
            .map(|_| zoda_math::Fq(rng.next_below(crate::rq::Q as u64) as u32))
            .collect()
    }
}
