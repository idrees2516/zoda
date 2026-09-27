//! Keccak-256 (the original Keccak padding, as used by Ethereum) plus
//! SHA3-256. Used for Ethereum-flavoured commitments in the sybil-resistance
//! and bridge modules.

const RC: [u64; 24] = [
    0x0000000000000001,
    0x0000000000008082,
    0x800000000000808a,
    0x8000000080008000,
    0x000000000000808b,
    0x0000000080000001,
    0x8000000080008081,
    0x8000000000008009,
    0x000000000000008a,
    0x0000000000000088,
    0x0000000080008009,
    0x000000008000000a,
    0x000000008000808b,
    0x800000000000008b,
    0x8000000000008089,
    0x8000000000008003,
    0x8000000000008002,
    0x8000000000000080,
    0x000000000000800a,
    0x800000008000000a,
    0x8000000080008081,
    0x8000000000008080,
    0x0000000080000001,
    0x8000000080008008,
];

/// Rotation offsets for the 24 rho steps (tiny-keccak layout).
const RHO: [u32; 24] = [
    1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44,
];

/// Pi-step destination order (tiny-keccak layout).
const PILN: [usize; 24] = [
    10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1,
];

#[inline(always)]
fn keccak_f(a: &mut [u64; 25]) {
    for &rc in RC.iter() {
        // theta
        let mut c = [0u64; 5];
        for x in 0..5 {
            c[x] = a[x] ^ a[x + 5] ^ a[x + 10] ^ a[x + 15] ^ a[x + 20];
        }
        let mut d = [0u64; 5];
        for x in 0..5 {
            d[x] = c[(x + 4) % 5] ^ c[(x + 1) % 5].rotate_left(1);
        }
        for y in 0..5 {
            for x in 0..5 {
                a[5 * y + x] ^= d[x];
            }
        }
        // rho + pi
        let mut last = a[1];
        for i in 0..24 {
            let j = PILN[i];
            let tmp = a[j];
            a[j] = last.rotate_left(RHO[i]);
            last = tmp;
        }
        // chi
        for y in 0..5 {
            let row = [a[5 * y], a[5 * y + 1], a[5 * y + 2], a[5 * y + 3], a[5 * y + 4]];
            for x in 0..5 {
                a[5 * y + x] = row[x] ^ ((!row[(x + 1) % 5]) & row[(x + 2) % 5]);
            }
        }
        // iota
        a[0] ^= rc;
    }
}

fn sponge(data: &[u8], rate: usize, pad_byte: u8, out_len: usize) -> Vec<u8> {
    let mut st = [0u64; 25];
    let mut offset = 0;
    // absorb full blocks
    let mut process_block = |st: &mut [u64; 25], block: &[u8]| {
        for (i, chunk) in block.chunks_exact(8).enumerate() {
            st[i] ^= u64::from_le_bytes(chunk.try_into().unwrap());
        }
        keccak_f(st);
    };
    while data.len() - offset >= rate {
        process_block(&mut st, &data[offset..offset + rate]);
        offset += rate;
    }
    // final block with padding
    let mut block = [0u8; 200];
    let rem = data.len() - offset;
    block[..rem].copy_from_slice(&data[offset..]);
    block[rem] = pad_byte;
    block[rate - 1] |= 0x80;
    process_block(&mut st, &block[..rate]);
    // squeeze
    let mut out = Vec::with_capacity(out_len);
    let mut pos = 0;
    while out.len() < out_len {
        if pos == 0 {
            for i in 0..out_len.min(rate) / 8 {
                out.extend_from_slice(&st[i].to_le_bytes());
            }
            pos = rate;
            if out.len() >= out_len {
                break;
            }
        } else {
            keccak_f(&mut st);
            for i in 0..out_len.min(rate) / 8 {
                out.extend_from_slice(&st[i].to_le_bytes());
            }
        }
    }
    out.truncate(out_len);
    out
}

/// Keccak-256 (Ethereum flavour: pad byte 0x01).
pub fn keccak256(data: &[u8]) -> [u8; 32] {
    let out = sponge(data, 136, 0x01, 32);
    let mut h = [0u8; 32];
    h.copy_from_slice(&out);
    h
}

/// SHA3-256 (FIPS 202, pad byte 0x06).
pub fn sha3_256(data: &[u8]) -> [u8; 32] {
    let out = sponge(data, 136, 0x06, 32);
    let mut h = [0u8; 32];
    h.copy_from_slice(&out);
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(b: &[u8]) -> String {
        b.iter().map(|x| format!("{:02x}", x)).collect()
    }

    #[test]
    fn keccak_vectors() {
        assert_eq!(
            hex(&keccak256(b"")),
            "c5d2460186f7233c927e7db2dcc703c0e500b653ca82273b7bfad8045d85a470"
        );
        assert_eq!(
            hex(&keccak256(b"abc")),
            "4e03657aea45a94fc7d47ba826c8d667c0d1e6e33a64a036ec44f58fa12d6c45"
        );
        // > rate boundary
        let data = vec![0x42u8; 200];
        let h = keccak256(&data);
        assert_eq!(h.len(), 32);
    }

    #[test]
    fn sha3_vectors() {
        assert_eq!(
            hex(&sha3_256(b"")),
            "a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a"
        );
        assert_eq!(
            hex(&sha3_256(b"abc")),
            "3a985da74fe225b2045c172d6bd390bd855f086e3e9d525b46bfe24511431532"
        );
    }
}
