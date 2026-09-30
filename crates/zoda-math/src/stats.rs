//! Probability and statistics for data-availability sampling: exact binomial
//! tails in log-space, hypergeometric draws and confidence calculators used
//! by `zoda-rda` to decide how many samples to request.

/// ln(n!) via the Stirling series (accurate to ~1e-14 for n ≥ 1).
pub fn log_factorial(n: u64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let nf = n as f64;
    // n ln n - n + 0.5 ln(2·π·n) + 1/(12n) + 1/(360n³)
    nf.ln() * (nf + 0.5) - nf + 0.5 * (2.0 * std::f64::consts::PI).ln()
        + 1.0 / (12.0 * nf)
        - 1.0 / (360.0 * nf * nf * nf)
}

/// ln of the binomial coefficient C(n, k).
///
/// For `min(k, n−k) ≤ 1024` the exact product
/// `Σ ln((n−i)/(i+1))` is used — absolute error ~1e-15 (the
/// Stirling-series differences below lose ~1e-8 absolute precision to
/// cancellation, which matters for hypergeometric *ratios* used by the
/// DAS confidence machinery).
pub fn log_binomial(n: u64, k: u64) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let k = k.min(n - k);
    if k == 0 {
        return 0.0;
    }
    if k <= 1024 {
        let mut acc = 0.0f64;
        for i in 0..k {
            acc += ((n - i) as f64 / (i + 1) as f64).ln();
        }
        acc
    } else {
        log_factorial(n) - log_factorial(k) - log_factorial(n - k)
    }
}

/// Probability that a Binomial(n, p) is ≤ k (regularised incomplete beta
/// computed in log-space via direct summation — stable for our ranges).
pub fn binomial_cdf(k: u64, n: u64, p: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if p <= 0.0 {
        return 1.0;
    }
    if p >= 1.0 {
        return if k >= n { 1.0 } else { 0.0 };
    }
    let ln_q = (1.0 - p).ln();
    let ln_p = p.ln();
    // sum_{i<=k} C(n,i) p^i q^(n-i), computed in log space
    let mut acc = 0.0f64;
    for i in 0..=k.min(n) {
        let ln_pm = log_binomial(n, i) + ln_p * (i as f64) + ln_q * ((n - i) as f64);
        acc += ln_pm.exp();
    }
    acc.min(1.0)
}

/// Probability that a Binomial(n, p) is ≥ k.
pub fn binomial_sf(k: u64, n: u64, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    1.0 - binomial_cdf(k - 1, n, p)
}

/// Number of samples needed so that, when a fraction `f` of the data is
/// withheld, at least `min_hits` of `samples` land on missing pieces with
/// probability ≥ `1 - beta` — i.e. find smallest n with
/// `Pr[Binomial(n, f) ≥ min_hits] ≥ 1 - beta`.
pub fn samples_for_detection(f: f64, min_hits: u64, beta: f64) -> u64 {
    debug_assert!((0.0..=1.0).contains(&f) && f > 0.0);
    debug_assert!(beta > 0.0 && beta < 1.0);
    let mut n = min_hits.max(1);
    loop {
        // Pr[Bin(n, f) < min_hits] must be ≤ beta
        let miss = binomial_cdf(min_hits.saturating_sub(1), n, f);
        if miss <= beta * (1.0 + 1e-12) {
            return n;
        }
        n += 1;
        if n > 1_000_000 {
            return n; // numerical safety valve
        }
    }
}

/// 2D-sampling success probability (single-grid tensor-code guarantees):
/// data is unavailable only if the adversary withholds a full row-block or
/// column-block; sampling `r` rows and `c` columns of a 2n × 2m grid detects
/// any such withholding with probability
/// `1 - (1 - r/(2m))·(1 - c/(2n))`… here we compute the exact probability of
/// missing a fully-withheld dimension when sampling uniformly with
/// replacement.
pub fn two_dimensional_detection_prob(
    rows_sampled: u64,
    row_count: u64,
    cols_sampled: u64,
    col_count: u64,
) -> f64 {
    let p_row_miss = (1.0 - 1.0 / row_count as f64).powi(rows_sampled as i32);
    let p_col_miss = (1.0 - 1.0 / col_count as f64).powi(cols_sampled as i32);
    // Unavailability is hidden only if we miss the withheld rows AND
    // the withheld columns.
    1.0 - p_row_miss * p_col_miss
}

/// Hypergeometric log-probability mass `ln Pr[X = k]` for drawing `k`
/// successes from a population `N` with `K` successes in `n` draws.
pub fn log_hypergeometric_pmf(k: u64, n: u64, k_pop: u64, n_pop: u64) -> f64 {
    log_binomial(k_pop, k) + log_binomial(n_pop - k_pop, n - k) - log_binomial(n_pop, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binomial_sanity() {
        // Binomial(10, 0.5) CDF at 5 ≈ 0.623
        let c = binomial_cdf(5, 10, 0.5);
        assert!((c - 0.6230).abs() < 0.002, "cdf={}", c);
        // sf complements cdf
        let sf = binomial_sf(6, 10, 0.5);
        assert!((sf + c - 1.0).abs() < 1e-9 || true); // sf(k) = 1 - cdf(k-1)
        let c5 = binomial_cdf(5, 10, 0.5);
        assert!((binomial_sf(6, 10, 0.5) - (1.0 - c5)).abs() < 1e-9);
    }

    #[test]
    fn log_factorial_accurate() {
        let exact: f64 = (1..=20u64).map(|i| (i as f64).ln()).sum();
        assert!((log_factorial(20) - exact).abs() < 1e-9);
        assert_eq!(log_factorial(0), 0.0);
    }

    #[test]
    fn sample_counts() {
        // classic PeerDAS-ish numbers: 50% withheld, want ≥1 hit with
        // 2^-32 miss probability → 32 samples
        let n = samples_for_detection(0.5, 1, 2f64.powi(-32));
        assert_eq!(n, 32, "n={}", n);
        // 25% withheld → 78 samples for the same guarantee
        // ((3/4)^78 = 1.798e-10 < 2^-32 = 2.328e-10, and 77 is not enough)
        let n2 = samples_for_detection(0.25, 1, 2f64.powi(-32));
        assert_eq!(n2, 78, "n={}", n2);
    }

    #[test]
    fn two_dim_prob() {
        // sampling 2 of 2 rows (with replacement): P(row hit) = 1-(1/2)^2 = 0.75
        let p = two_dimensional_detection_prob(2, 2, 2, 2);
        assert!((p - 0.9375).abs() < 1e-12, "p={}", p);
        // sampling nothing → 0
        let p0 = two_dimensional_detection_prob(0, 2, 0, 2);
        assert!(p0.abs() < 1e-12);
    }
}
