//! The BLS12-381 G2 group: y² = x³ + 4(1 + u) over Fp2.

use crate::curve::{Affine, CurveConfig, Projective};
use crate::fp::Fp;
use crate::fp2::Fp2;

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct G2Config;

impl CurveConfig for G2Config {
    type Base = Fp2;

    fn b() -> Fp2 {
        Fp2 {
            c0: Fp::from_u64(4),
            c1: Fp::from_u64(4),
        }
    }
    fn b3() -> Fp2 {
        Fp2 {
            c0: Fp::from_u64(12),
            c1: Fp::from_u64(12),
        }
    }
    fn generator_xy() -> (Fp2, Fp2) {
        // Canonical BLS12-381 G2 generator (Montgomery limbs).
        let x = Fp2 {
            c0: Fp([
                0xf5f28fa202940a10,
                0xb3f5fb2687b4961a,
                0xa1a893b53e2ae580,
                0x9894999d1a3caee9,
                0x6f67b7631863366b,
                0x058191924350bcd7,
            ]),
            c1: Fp([
                0xa5a9c0759e23f606,
                0xaaa0c59dbccd60c3,
                0x3bb17e18e2867806,
                0x1b1ab6cc8541b367,
                0xc2b6ed0ef2158547,
                0x11922a097360edf3,
            ]),
        };
        let y = Fp2 {
            c0: Fp([
                0x4c730af860494c4a,
                0x597cfa1f5e369c5a,
                0xe7e6856caa0a635a,
                0xbbefb5e96e0d495f,
                0x07d3a975f0ef25a2,
                0x0083fd8e7e80dae5,
            ]),
            c1: Fp([
                0xadc0fc92df64b05d,
                0x18aa270a2b1461dc,
                0x86adac6a3be4eba0,
                0x79495c4ec93da33a,
                0xe7175850a43ccaed,
                0x0b2bc2a163de1bf2,
            ]),
        };
        (x, y)
    }
    fn elem_bytes() -> usize {
        96
    }
    fn base_to_be_bytes(x: &Fp2, out: &mut [u8]) {
        let c1 = x.c1.to_repr();
        let c0 = x.c0.to_repr();
        // G2 serialisation order: c1 first, then c0
        for i in 0..6 {
            out[i * 8..i * 8 + 8].copy_from_slice(&c1[5 - i].to_be_bytes());
        }
        for i in 0..6 {
            out[48 + i * 8..48 + i * 8 + 8].copy_from_slice(&c0[5 - i].to_be_bytes());
        }
    }
    fn base_from_be_bytes(bytes: &[u8]) -> Option<Fp2> {
        if bytes.len() != 96 {
            return None;
        }
        let mut c1_limbs = [0u64; 6];
        let mut c0_limbs = [0u64; 6];
        for i in 0..6 {
            let mut w = [0u8; 8];
            w.copy_from_slice(&bytes[i * 8..i * 8 + 8]);
            c1_limbs[5 - i] = u64::from_be_bytes(w);
            let mut w = [0u8; 8];
            w.copy_from_slice(&bytes[48 + i * 8..48 + i * 8 + 8]);
            c0_limbs[5 - i] = u64::from_be_bytes(w);
        }
        let check = |limbs: &[u64; 6]| -> bool {
            for i in (0..6).rev() {
                if limbs[i] > Fp::MODULUS[i] {
                    return false;
                }
                if limbs[i] < Fp::MODULUS[i] {
                    return true;
                }
            }
            true
        };
        if !check(&c1_limbs) || !check(&c0_limbs) {
            return None;
        }
        Some(Fp2 {
            c0: Fp::from_repr_limbs(c0_limbs),
            c1: Fp::from_repr_limbs(c1_limbs),
        })
    }
}

pub type G2Projective = Projective<G2Config>;
pub type G2Affine = Affine<G2Config>;

impl G2Projective {
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
        let g = G2Affine::generator();
        assert!(g.is_on_curve());
    }

    #[test]
    fn group_axioms() {
        let g = G2Projective::generator();
        assert_eq!(g.add(&g).to_affine(), g.double().to_affine());
        let a = G2Projective::generator().mul_limbs(&[7, 0, 0, 0]);
        let b = G2Projective::generator().mul_limbs(&[11, 0, 0, 0]);
        let c = G2Projective::generator().mul_limbs(&[13, 0, 0, 0]);
        assert_eq!(
            a.add(&b).add(&c).to_affine(),
            a.add(&b.add(&c)).to_affine()
        );
        assert!(a.add(&a.neg()).is_identity());
        let a18 = G2Projective::generator().mul_limbs(&[18, 0, 0, 0]);
        assert_eq!(a.add(&b).to_affine(), a18.to_affine());
        for pt in [&a, &b, &c, &a18] {
            assert!(pt.to_affine().is_on_curve());
        }
    }

    #[test]
    fn compression_roundtrip() {
        let g = G2Affine::generator();
        let bytes = g.to_compressed();
        assert_eq!(bytes.len(), 96);
        let back = G2Affine::from_compressed(&bytes).unwrap();
        assert_eq!(back, g);

        let p = G2Projective::generator()
            .mul_limbs(&[0xcafebabe, 0xfeedface, 0x0badc0de, 0x0defaced])
            .to_affine();
        let bytes = p.to_compressed();
        let back = G2Affine::from_compressed(&bytes).unwrap();
        assert_eq!(back, p);

        let inf = G2Affine::identity();
        let back = G2Affine::from_compressed(&inf.to_compressed()).unwrap();
        assert_eq!(back, inf);
    }

    #[test]
    fn known_generator_hex() {
        // The compressed G2 generator must start with 0x93 (compressed,
        // not infinity, not lex-largest).
        let bytes = G2Affine::generator().to_compressed();
        let hex: String = bytes.iter().map(|b| format!("{:02x}", b)).collect();
        assert!(hex.starts_with("93e02b6052719f607dacd3a088274f6559"), "got {}", hex);
    }
}
