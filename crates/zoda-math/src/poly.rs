//! Polynomial utilities over any [`PrimeField`]: Horner and barycentric
//! evaluation, subproduct-tree interpolation and multi-point evaluation
//! (O(n log² n)), vanishing polynomials and synthetic division.
//!
//! These power the ZODA tensor-code encoder, erasure reconstruction and the
//! data-availability sampling probability engine.

use crate::{batch_invert, PrimeField};

/// Dense coefficient-form polynomial (index = degree).
pub type Poly<F> = Vec<F>;

/// Evaluate `p` at `x` (Horner's rule).
pub fn eval<F: PrimeField>(p: &[F], x: F) -> F {
    let mut acc = F::zero();
    for c in p.iter().rev() {
        acc = acc * x + *c;
    }
    acc
}

/// Multiply two dense polynomials.
pub fn mul<F: PrimeField>(a: &[F], b: &[F]) -> Poly<F> {
    if a.is_empty() || b.is_empty() {
        return vec![F::zero(); 0];
    }
    let mut out = vec![F::zero(); a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        if ai.is_zero() {
            continue;
        }
        for (j, &bj) in b.iter().enumerate() {
            out[i + j] = out[i + j] + ai * bj;
        }
    }
    out
}

/// Fast multi-point evaluation of polynomial `p` (deg < n) at `n` points
/// using the remainder tree. `points.len()` must be a power of two.
pub fn multi_evaluate<F: PrimeField>(p: &[F], points: &[F]) -> Vec<F> {
    let n = points.len();
    debug_assert!(n.is_power_of_two());
    if n == 1 {
        return vec![eval(p, points[0])];
    }
    // polys[l][i] = product of (X - x) over subtree i at level l
    // (l = 0: leaves, one per point; last level: single root).
    let levels = n.trailing_zeros() as usize + 1;
    let mut polys: Vec<Vec<Poly<F>>> = Vec::with_capacity(levels);
    let mut leaves = Vec::with_capacity(n);
    for &x in points {
        leaves.push(vec![-x, F::ONE]);
    }
    polys.push(leaves);
    let mut sz = n >> 1;
    while sz >= 1 {
        let prev = polys.last().unwrap();
        let mut cur = Vec::with_capacity(sz);
        for i in 0..sz {
            cur.push(mul(&prev[2 * i], &prev[2 * i + 1]));
        }
        polys.push(cur);
        if sz == 1 {
            break;
        }
        sz >>= 1;
    }

    // Descend the tree with polynomial remainders, starting at the root.
    let mut cur = vec![trim(p.to_vec())];
    let mut l = polys.len() - 1;
    while l > 0 {
        let mut next = Vec::with_capacity(cur.len() * 2);
        for (i, c) in cur.iter().enumerate() {
            next.push(rem(c, &polys[l - 1][2 * i]));
            next.push(rem(c, &polys[l - 1][2 * i + 1]));
        }
        cur = next;
        l -= 1;
    }
    debug_assert_eq!(cur.len(), n);
    cur.iter()
        .map(|c| c.first().copied().unwrap_or(F::zero()))
        .collect()
}

/// Polynomial remainder `a mod b` (synthetic-division style long division).
/// `b` must be non-empty with nonzero leading coefficient.
pub fn rem<F: PrimeField>(a: &[F], b: &[F]) -> Poly<F> {
    debug_assert!(!b.is_empty() && !b[b.len() - 1].is_zero());
    let db = b.len() - 1;
    if a.is_empty() {
        return vec![];
    }
    let mut r = a.to_vec();
    let inv_lead = b[db].invert().unwrap();
    if r.len() < b.len() {
        return r;
    }
    for i in (db..r.len()).rev() {
        let coef = r[i] * inv_lead;
        if coef.is_zero() {
            continue;
        }
        let shift = i - db;
        for j in 0..=db {
            r[shift + j] = r[shift + j] - coef * b[j];
        }
    }
    r.truncate(db);
    while r.len() > 0 && r[r.len() - 1].is_zero() {
        r.pop();
    }
    r
}

/// Fast interpolation: returns the coefficient-form polynomial of degree
/// `< points.len()` passing through `(points[i], values[i])`.
/// `points.len()` must be a power of two and points must be distinct.
pub fn interpolate<F: PrimeField>(points: &[F], values: &[F]) -> Poly<F> {
    let n = points.len();
    assert_eq!(points.len(), values.len());
    debug_assert!(n.is_power_of_two());
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![values[0]];
    }

    // Build subproduct tree (levels of vanishing polynomials).
    let levels = n.trailing_zeros() as usize + 1;
    let mut polys: Vec<Vec<Poly<F>>> = Vec::with_capacity(levels);
    let mut leaves = Vec::with_capacity(n);
    for &x in points {
        leaves.push(vec![-x, F::ONE]);
    }
    polys.push(leaves);
    let mut sz = n >> 1;
    while sz >= 1 {
        let prev = polys.last().unwrap();
        let mut cur = Vec::with_capacity(sz);
        for i in 0..sz {
            cur.push(mul(&prev[2 * i], &prev[2 * i + 1]));
        }
        polys.push(cur);
        if sz == 1 {
            break;
        }
        sz >>= 1;
    }

    // derivative of the root vanishing polynomial, evaluated at all points
    let root = polys.last().unwrap()[0].clone();
    let droot = derivative(&root);
    let denominators = multi_evaluate(&droot, points);
    // invert batch
    let mut denoms = denominators;
    batch_invert(&mut denoms);

    // Linear combination form: P(X) = Σ values[i]/droot(x_i) · Π_{j≠i}(X - x_j)
    // computed by merging remainders up the tree (Lagrange via remainders).
    // leaf-level Lagrange "basis remainders": at the leaf level each node
    // polynomial is (X - x_i); we store the coefficient values.
    let mut cur: Vec<Poly<F>> = points
        .iter()
        .zip(values.iter())
        .zip(denoms.iter())
        .map(|((&_x, &v), &d)| vec![v * d])
        .collect();

    // Linear (Lagrange) merge up the tree:
    //   P_{S1∪S2} = P_{S1}·Z_{S2} + P_{S2}·Z_{S1}
    // where Z_S is the vanishing polynomial of the subtree point set.
    let mut level = 1usize;
    while level < polys.len() {
        let prev_tree = &polys[level - 1];
        let mut next = Vec::with_capacity(cur.len() / 2);
        let mut i = 0;
        while i < cur.len() {
            let za = &prev_tree[i];
            let zb = &prev_tree[i + 1];
            let ra = &cur[i];
            let rb = &cur[i + 1];
            let left = mul(ra, zb);
            let right = mul(rb, za);
            next.push(add_mod(&left, &right, zb));
            i += 2;
        }
        cur = next;
        level += 1;
    }
    cur.into_iter().next().unwrap_or_default()
}

/// `a - b` keeping the result reduced mod `m` (polynomial sense).
#[allow(dead_code)]
fn sub_mod<F: PrimeField>(a: &[F], b: &[F], _m: &[F]) -> Poly<F> {
    let n = a.len().max(b.len());
    let mut out = vec![F::zero(); n];
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(F::zero());
        let y = b.get(i).copied().unwrap_or(F::zero());
        out[i] = x - y;
    }
    while out.len() > 0 && out[out.len() - 1].is_zero() {
        out.pop();
    }
    out
}

fn add_mod<F: PrimeField>(a: &[F], b: &[F], _m: &[F]) -> Poly<F> {
    let n = a.len().max(b.len());
    let mut out = vec![F::zero(); n];
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(F::zero());
        let y = b.get(i).copied().unwrap_or(F::zero());
        out[i] = x + y;
    }
    while out.len() > 0 && out[out.len() - 1].is_zero() {
        out.pop();
    }
    out
}

/// Polynomial division: returns `(quotient, remainder)`.
pub fn divmod<F: PrimeField>(a: &[F], b: &[F]) -> (Poly<F>, Poly<F>) {
    debug_assert!(!b.is_empty() && !b[b.len() - 1].is_zero());
    let db = b.len() - 1;
    let mut r = a.to_vec();
    if a.len() < b.len() {
        return (vec![], trim(a.to_vec()));
    }
    let inv_lead = b[db].invert().unwrap();
    let mut q = vec![F::zero(); a.len() - db];
    for i in (db..r.len()).rev() {
        let coef = r[i] * inv_lead;
        q[i - db] = coef;
        if coef.is_zero() {
            continue;
        }
        let shift = i - db;
        for j in 0..=db {
            r[shift + j] = r[shift + j] - coef * b[j];
        }
    }
    r.truncate(db);
    (trim(q), trim(r))
}

/// Derivative of a dense polynomial.
pub fn derivative<F: PrimeField>(p: &[F]) -> Poly<F> {
    if p.len() <= 1 {
        return vec![];
    }
    p[1..]
        .iter()
        .enumerate()
        .map(|(i, &c)| c * F::from_u64(i as u64 + 1))
        .collect()
}

/// Trim trailing zero coefficients.
pub fn trim<F: PrimeField>(mut p: Poly<F>) -> Poly<F> {
    while p.len() > 0 && p[p.len() - 1].is_zero() {
        p.pop();
    }
    p
}

/// Vanishing polynomial over arbitrary points (O(n²) product; used for small
/// root sets such as EIP-7594 cell recovery).
pub fn vanishing_poly<F: PrimeField>(roots: &[F]) -> Poly<F> {
    let mut p = vec![F::ONE];
    for &r in roots {
        // p *= (X - r)
        let mut next = vec![F::zero(); p.len() + 1];
        for (i, &c) in p.iter().enumerate() {
            next[i + 1] = next[i + 1] + c;
            next[i] = next[i] - c * r;
        }
        p = next;
    }
    p
}

/// Barycentric evaluation of a polynomial given in evaluation form over a
/// multiplicative subgroup domain (size n, natural order roots of unity).
///
/// Returns `P(z)` in O(n) given `evals[i] = P(ω^i)`.
pub fn barycentric_eval_from_domain<F: PrimeField>(
    evals: &[F],
    z: F,
    domain_roots: &[F],
) -> F {
    let n = evals.len();
    debug_assert_eq!(domain_roots.len(), n);
    let inv_n = F::from_u64(n as u64).invert().unwrap();
    // If z is on the domain, direct lookup
    for (i, &w) in domain_roots.iter().enumerate() {
        if w == z {
            return evals[i];
        }
    }
    // f(z) = (z^n - 1)/n · Σ f_i·ω_i / (z - ω_i)
    let mut denom = Vec::with_capacity(n);
    for &w in domain_roots {
        denom.push(z - w);
    }
    batch_invert(&mut denom);
    let mut acc = F::zero();
    for i in 0..n {
        acc = acc + evals[i] * domain_roots[i] * denom[i];
    }
    // z^n - 1 via repeated squaring on z
    let mut zn = F::ONE;
    let mut e = n as u64;
    let mut base = z;
    while e > 0 {
        if e & 1 == 1 {
            zn = zn * base;
        }
        base = base * base;
        e >>= 1;
    }
    acc * (zn - F::ONE) * inv_n
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Fr;

    fn rand_polys(seed: &[u8; 32], n: usize) -> (Vec<Fr>, Vec<Fr>) {
        let mut rng = crate::ZodaRng::from_seed(*seed);
        let p: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
        let pts: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
        (p, pts)
    }

    #[test]
    fn eval_and_divmod() {
        let (p, _) = rand_polys(b"poly-test-1-divmod-0000000000000", 17);
        let x = Fr::from_u64(99);
        let v = eval(&p, x);
        let d: Vec<Fr> = vec![Fr::ONE, Fr::from_u64(0), Fr::from_u64(0), Fr::ONE]; // X^3+1
        let (q, r) = divmod(&p, &d);
        // p == q*d + r
        let qd = mul(&q, &d);
        let mut recomp = vec![Fr::zero(); qd.len().max(r.len())];
        for i in 0..recomp.len() {
            let a = qd.get(i).copied().unwrap_or(Fr::zero());
            let b = r.get(i).copied().unwrap_or(Fr::zero());
            recomp[i] = a + b;
        }
        for i in 0..p.len() {
            assert_eq!(p[i], recomp[i]);
        }
        // and at the root of d, r must agree with p
        for trial in 0..3 {
            let _ = trial;
            let root = d[0].neg(); // -d[0] is a root of X^3+1
            assert_eq!(eval(&r, root), eval(&p, root));
        }
    }

    #[test]
    fn multipoint_and_interpolation_roundtrip() {
        for log in 1..=6u32 {
            let n = 1usize << log;
            let mut rng = crate::ZodaRng::from_seed(*b"poly-test-2-mpi-0000000000000000");
            let p: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
            let pts: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
            let vals = multi_evaluate(&p, &pts);
            for i in 0..n {
                assert_eq!(vals[i], eval(&p, pts[i]), "n={} i={}", n, i);
            }
            let rec = interpolate(&pts, &vals);
            assert_eq!(rec.len(), p.len(), "n={} length", n);
            for i in 0..n {
                assert_eq!(rec[i], p[i], "n={} coeff {}", n, i);
            }
        }
    }

    #[test]
    fn vanishing() {
        let mut rng = crate::ZodaRng::from_seed(*b"poly-test-3-van-0000000000000000");
        let roots: Vec<Fr> = (0..9).map(|_| rng.next_fr(false)).collect();
        let vp = vanishing_poly(&roots);
        for &r in &roots {
            assert!(eval(&vp, r).is_zero());
        }
    }

    #[test]
    fn barycentric() {
        let n = 64;
        let d = crate::FftDomain::<Fr>::new(n);
        let mut rng = crate::ZodaRng::from_seed(*b"poly-test-4-bar-0000000000000000");
        let coeffs: Vec<Fr> = (0..n).map(|_| rng.next_fr(false)).collect();
        let mut evals = vec![Fr::zero(); n];
        d.fft(&coeffs, &mut evals);
        let roots = &d.roots_of_unity()[..n];
        // on-domain
        assert_eq!(
            barycentric_eval_from_domain(&evals, roots[13], roots),
            evals[13]
        );
        // off-domain
        let z = rng.next_fr(false);
        assert_eq!(
            barycentric_eval_from_domain(&evals, z, roots),
            eval(&coeffs, z)
        );
    }
}
