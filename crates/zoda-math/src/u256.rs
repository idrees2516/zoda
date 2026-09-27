//! 256-bit unsigned integer helpers (comparison, shifts, hex parsing) used
//! by byte-oriented APIs across the stack.

/// Compare two 4-limb little-endian integers.
#[inline]
pub fn cmp4(a: &[u64; 4], b: &[u64; 4]) -> core::cmp::Ordering {
    for i in (0..4).rev() {
        if a[i] != b[i] {
            return a[i].cmp(&b[i]);
        }
    }
    core::cmp::Ordering::Equal
}

/// Parse a hex string (with or without `0x` prefix, any whitespace/newlines)
/// into bytes. Returns `None` on invalid characters.
pub fn hex_to_bytes(s: &str) -> Option<Vec<u8>> {
    let cleaned: String = s
        .chars()
        .filter(|c| !c.is_whitespace() && *c != ',')
        .collect();
    let s = cleaned.strip_prefix("0x").unwrap_or(&cleaned);
    if s.len() % 2 != 0 {
        return None;
    }
    let mut out = Vec::with_capacity(s.len() / 2);
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let hi = hex_val(bytes[i])?;
        let lo = hex_val(bytes[i + 1])?;
        out.push(hi << 4 | lo);
        i += 2;
    }
    Some(out)
}

#[inline]
fn hex_val(c: u8) -> Option<u8> {
    match c {
        b'0'..=b'9' => Some(c - b'0'),
        b'a'..=b'f' => Some(c - b'a' + 10),
        b'A'..=b'F' => Some(c - b'A' + 10),
        _ => None,
    }
}

/// Lower-hex encode bytes.
pub fn bytes_to_hex(b: &[u8]) -> String {
    let mut s = String::with_capacity(b.len() * 2);
    for x in b {
        s.push_str(&format!("{:02x}", x));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hex_roundtrip() {
        let h = "0x0123abcd";
        let b = hex_to_bytes(h).unwrap();
        assert_eq!(b, vec![0x01, 0x23, 0xab, 0xcd]);
        assert_eq!(bytes_to_hex(&b), "0123abcd");
        assert!(hex_to_bytes("zz").is_none());
        // whitespace + newlines (trusted setup file format)
        let messy = "a0 41\n3c\t0d, ca";
        let b2 = hex_to_bytes(messy).unwrap();
        assert_eq!(b2, vec![0xa0, 0x41, 0x3c, 0x0d, 0xca]);
    }

    #[test]
    fn ordering() {
        let a = [1u64, 0, 0, 0];
        let b = [2u64, 0, 0, 0];
        let c = [1u64, 1, 0, 0];
        assert_eq!(cmp4(&a, &b), core::cmp::Ordering::Less);
        assert_eq!(cmp4(&b, &a), core::cmp::Ordering::Greater);
        assert_eq!(cmp4(&a, &c), core::cmp::Ordering::Less);
        assert_eq!(cmp4(&a, &a), core::cmp::Ordering::Equal);
    }
}
