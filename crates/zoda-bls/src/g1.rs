//! The BLS12-381 G1 group: y² = x³ + 4 over Fp.

use crate::curve::{Affine, CurveConfig, Projective};
use crate::fp::Fp;

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct G1Config;

impl CurveConfig for G1Config {
    type Base = Fp;

    fn b() -> Fp {
        Fp::from_u64(4)
    }
    fn b3() -> Fp {
        Fp::from_u64(12)
    }
    fn generator_xy() -> (Fp, Fp) {
        // Canonical BLS12-381 G1 generator (Montgomery limbs).
        let x = Fp([
            0x5cb38790fd530c16,
            0x7817fc679976fff5,
            0x154f95c7143ba1c1,
            0xf0ae6acdf3d0e747,
            0xedce6ecc21dbf440,
            0x120177419e0bfb75,
        ]);
        let y = Fp([
            0xbaac93d50ce72271,
            0x8c22631a7918fd8e,
            0xdd595f13570725ce,
            0x51ac582950405194,
            0x0e1c8c3fad0059c0,
            0x0bbc3efc5008a26a,
        ]);
        (x, y)
    }
    fn elem_bytes() -> usize {
        48
    }
    fn base_to_be_bytes(x: &Fp, out: &mut [u8]) {
        let r = x.to_repr();
        for i in 0..6 {
            out[i * 8..i * 8 + 8].copy_from_slice(&r[5 - i].to_be_bytes());
        }
    }
    fn base_from_be_bytes(bytes: &[u8]) -> Option<Fp> {
        if bytes.len() != 48 {
            return None;
        }
        let mut limbs = [0u64; 6];
        for i in 0..6 {
            let mut w = [0u8; 8];
            w.copy_from_slice(&bytes[i * 8..i * 8 + 8]);
            limbs[5 - i] = u64::from_be_bytes(w);
        }
        // validate < p
        for i in (0..6).rev() {
            if limbs[i] > Fp::MODULUS[i] {
                return None;
            }
            if limbs[i] < Fp::MODULUS[i] {
                break;
            }
        }
        Some(Fp::from_repr_limbs(limbs))
    }
}

pub type G1Projective = Projective<G1Config>;
pub type G1Affine = Affine<G1Config>;

impl G1Projective {
    /// Multiply by the Fr scalar (little-endian limbs).
    pub fn mul_fr(&self, scalar: &zoda_math::Fr) -> Self {
        self.mul_limbs(&scalar.to_repr())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generator_on_curve() {
        let g = G1Affine::generator();
        assert!(g.is_on_curve());
        assert!(!g.is_identity());
    }

    #[test]
    fn group_axioms() {
        let g = G1Projective::generator();
        let two_g = g.double();
        assert_eq!(g.add(&g).to_affine(), two_g.to_affine());
        // (a+b)+c == a+(b+c)
        let a = G1Projective::generator().mul_limbs(&[7, 0, 0, 0]);
        let b = G1Projective::generator().mul_limbs(&[11, 0, 0, 0]);
        let c = G1Projective::generator().mul_limbs(&[13, 0, 0, 0]);
        assert_eq!(
            a.add(&b).add(&c).to_affine(),
            a.add(&b.add(&c)).to_affine()
        );
        // scalar mult is a homomorphism
        let a7 = G1Projective::generator().mul_limbs(&[7, 0, 0, 0]);
        let a11 = G1Projective::generator().mul_limbs(&[11, 0, 0, 0]);
        let a18 = G1Projective::generator().mul_limbs(&[18, 0, 0, 0]);
        assert_eq!(a7.add(&a11).to_affine(), a18.to_affine());
        // negation
        assert!(a7.add(&a7.neg()).is_identity());
        // all computed points on the curve
        for pt in [&a, &b, &c, &a18, &two_g] {
            assert!(pt.to_affine().is_on_curve());
        }
    }

    #[test]
    fn compression_roundtrip() {
        let g = G1Affine::generator();
        let bytes = g.to_compressed();
        assert_eq!(bytes.len(), 48);
        let back = G1Affine::from_compressed(&bytes).unwrap();
        assert_eq!(back, g);

        // a random-ish point
        let p = G1Projective::generator()
            .mul_limbs(&[0xdeadbeef, 0x12345678, 0x9abcdef0, 0x0fedcba9])
            .to_affine();
        let bytes = p.to_compressed();
        let back = G1Affine::from_compressed(&bytes).unwrap();
        assert_eq!(back, p);

        // infinity
        let inf = G1Affine::identity();
        let bytes = inf.to_compressed();
        let back = G1Affine::from_compressed(&bytes).unwrap();
        assert_eq!(back, inf);

        // invalid: x not on curve
        let mut bad = G1Affine::generator().to_compressed();
        bad[47] ^= 0x01; // flip a low bit of x
        if G1Affine::from_compressed(&bad).is_some() {
            // not every x flip leaves the curve; only fail if it parsed AND
            // is not on the curve
            let parsed = G1Affine::from_compressed(&bad).unwrap();
            assert!(parsed.is_on_curve());
        }
    }

    #[test]
    fn known_generator_hex() {
        // The compressed G1 generator must start with 0x97 (compressed,
        // not infinity, not lex-largest).
        let bytes = G1Affine::generator().to_compressed();
        assert_eq!(bytes[0] & 0x80, 0x80);
        let hex: String = bytes.iter().map(|b| format!("{:02x}", b)).collect();
        assert!(hex.starts_with("97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb"), "got {}", hex);
    }
}
