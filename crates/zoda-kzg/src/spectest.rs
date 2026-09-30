//! Official consensus-spec-tests KZG vector harness.
//!
//! Consumes the `kzg_mainnet` suites from the consensus-spec-tests
//! release archives (`tests/general/<fork>/kzg/...`) plus the mainnet
//! ceremony `trusted_setup.txt`, and replays every case through this
//! crate's EIP-4844 API. `zoda-edas` reuses the parser and directory
//! plumbing for the EIP-7594 (Fulu) cell suites.
//!
//! The tests are **optional**: they are skipped unless the vector
//! archive is present. Point the harness at an archive with
//! `ZODA_KZG_VECTORS=<dir>` (and optionally
//! `ZODA_TRUSTED_SETUP=<file>`), or place both under
//! `<workspace>/spec-vectors/`. `scripts/fetch_kzg_vectors.sh` does the
//! download.
//!
//! YAML subset understood (the consensus tests' "simple YAML"): two
//! levels of mappings, block lists (possibly nested, as emitted for
//! `[cells, proofs]` outputs), inline flow lists (possibly wrapped over
//! continuation lines, as emitted for `cell_indices`), and single-quoted
//! scalar strings.

use crate::srs::Setup;
use std::path::{Path, PathBuf};

/// A parsed YAML value.
#[derive(Clone, Debug, PartialEq)]
pub enum Sv {
    Str(String),
    Bool(bool),
    Null,
    Int(u64),
    List(Vec<Sv>),
    Map(Vec<(String, Sv)>),
}

impl Sv {
    pub fn get(&self, key: &str) -> Option<&Sv> {
        match self {
            Sv::Map(pairs) => pairs.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Sv::Str(s) => Some(s),
            _ => None,
        }
    }
    pub fn as_list(&self) -> Option<&[Sv]> {
        match self {
            Sv::List(v) => Some(v),
            _ => None,
        }
    }
    pub fn is_null(&self) -> bool {
        matches!(self, Sv::Null)
    }
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Sv::Bool(b) => Some(*b),
            _ => None,
        }
    }
    /// "0x…" → bytes.
    pub fn hex(&self) -> Option<Vec<u8>> {
        self.as_str().and_then(|s| crate::srs::hex_bytes(s))
    }
    /// A fixed-size byte array.
    pub fn hex_array<const N: usize>(&self) -> Option<[u8; N]> {
        let v = self.hex()?;
        if v.len() != N {
            return None;
        }
        let mut out = [0u8; N];
        out.copy_from_slice(&v);
        Some(out)
    }
}

/// Parse the simple-YAML subset used by consensus-spec-tests.
pub fn parse_simple_yaml(text: &str) -> Sv {
    // raw lines with (indent, content); empty lines dropped
    let mut lines: Vec<(usize, String)> = Vec::new();
    for raw in text.lines() {
        if raw.trim().is_empty() {
            continue;
        }
        let indent = raw.len() - raw.trim_start().len();
        lines.push((indent, raw.trim_start().to_string()));
    }
    // join wrapped flow lists: while a line has unbalanced '[' / ']', merge
    let mut joined: Vec<(usize, String)> = Vec::with_capacity(lines.len());
    let mut i = 0;
    while i < lines.len() {
        let (ind, mut content) = lines[i].clone();
        let mut depth = bracket_depth(&content);
        let mut j = i + 1;
        while depth > 0 && j < lines.len() {
            content.push(' ');
            content.push_str(&lines[j].1);
            depth = bracket_depth(&content);
            j += 1;
        }
        joined.push((ind, content));
        i = j;
    }
    let mut idx = 0;
    if joined.is_empty() {
        return Sv::Null;
    }
    let first_indent = joined[0].0;
    parse_block(&mut joined, &mut idx, first_indent)
}

fn bracket_depth(s: &str) -> i32 {
    let mut d = 0i32;
    let mut in_quote = false;
    for ch in s.chars() {
        if ch == '\'' {
            in_quote = !in_quote;
        }
        if !in_quote {
            if ch == '[' {
                d += 1;
            } else if ch == ']' {
                d -= 1;
            }
        }
    }
    d
}

/// Parse the value block starting at `lines[idx]` whose first line sits
/// at indentation `indent` (a map, a list, or — for single scalar lines
/// — a scalar).
fn parse_block(lines: &mut [(usize, String)], idx: &mut usize, indent: usize) -> Sv {
    if *idx >= lines.len() {
        return Sv::Null;
    }
    let (ind, content) = lines[*idx].clone();
    let _ = ind;
    if content == "-" || content.starts_with("- ") {
        parse_list(lines, idx, indent)
    } else if content.contains(':') && !starts_flow(&content) {
        parse_map(lines, idx, indent)
    } else {
        // bare scalar line
        *idx += 1;
        parse_scalar(&content)
    }
}

fn starts_flow(s: &str) -> bool {
    s.starts_with('[')
}

fn parse_list(lines: &mut [(usize, String)], idx: &mut usize, indent: usize) -> Sv {
    let mut items = Vec::new();
    while *idx < lines.len() {
        let (ind, content) = lines[*idx].clone();
        if ind != indent || !(content == "-" || content.starts_with("- ")) {
            break;
        }
        if content == "-" {
            // item value is a nested block on the following lines
            *idx += 1;
            if *idx < lines.len() && lines[*idx].0 > indent {
                let child_indent = lines[*idx].0;
                items.push(parse_block(lines, idx, child_indent));
            } else {
                items.push(Sv::Null);
            }
            continue;
        }
        let rest = content[2..].to_string();
        if rest.is_empty() {
            items.push(Sv::Null);
            *idx += 1;
            continue;
        }
        if rest.starts_with("- ") || starts_flow(&rest) || rest.contains(": ") {
            // the item is itself a composite: rewrite the line with the
            // "- " replaced by deeper indentation and reparse in place
            lines[*idx] = (indent + 2, rest);
            items.push(parse_block(lines, idx, indent + 2));
        } else {
            // plain scalar item (possibly quoted)
            *idx += 1;
            items.push(parse_scalar(&rest));
        }
    }
    Sv::List(items)
}

fn parse_map(lines: &mut [(usize, String)], idx: &mut usize, indent: usize) -> Sv {
    let mut pairs: Vec<(String, Sv)> = Vec::new();
    while *idx < lines.len() {
        let (ind, content) = lines[*idx].clone();
        if ind != indent {
            break;
        }
        if content == "-" || content.starts_with("- ") {
            break;
        }
        let Some(colon) = content.find(':') else {
            break;
        };
        let key = content[..colon].trim().to_string();
        let value_part = content[colon + 1..].trim().to_string();
        if value_part.is_empty() {
            // value is a nested block: a map on deeper indentation, or a
            // list at indentation >= this key's (YAML allows same-indent
            // block lists under a key)
            *idx += 1;
            if *idx < lines.len() {
                let (nind, ncontent) = lines[*idx].clone();
                let is_list = ncontent == "-" || ncontent.starts_with("- ");
                if is_list && nind >= indent {
                    pairs.push((key, parse_list(lines, idx, nind)));
                    continue;
                }
                if nind > indent {
                    pairs.push((key, parse_block(lines, idx, nind)));
                    continue;
                }
            }
            pairs.push((key, Sv::Null));
        } else {
            *idx += 1;
            pairs.push((key, parse_scalar(&value_part)));
        }
    }
    Sv::Map(pairs)
}

fn parse_scalar(s: &str) -> Sv {
    let s = s.trim();
    if s.len() >= 2 && ((s.starts_with('\'') && s.ends_with('\'')) || (s.starts_with('"') && s.ends_with('"'))) {
        return Sv::Str(s[1..s.len() - 1].to_string());
    }
    if s.starts_with('[') && s.ends_with(']') {
        let inner = &s[1..s.len() - 1];
        let mut items = Vec::new();
        for tok in split_flow(inner) {
            let t = tok.trim();
            if !t.is_empty() {
                items.push(parse_scalar(t));
            }
        }
        return Sv::List(items);
    }
    match s {
        "null" | "~" | "" => Sv::Null,
        "true" => Sv::Bool(true),
        "false" => Sv::Bool(false),
        _ => {
            if let Ok(i) = s.parse::<u64>() {
                Sv::Int(i)
            } else {
                Sv::Str(s.to_string())
            }
        }
    }
}

fn split_flow(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut depth = 0i32;
    let mut in_quote = false;
    for ch in s.chars() {
        if ch == '\'' {
            in_quote = !in_quote;
        }
        if !in_quote {
            if ch == '[' {
                depth += 1;
            } else if ch == ']' {
                depth -= 1;
            } else if ch == ',' && depth == 0 {
                out.push(core::mem::take(&mut cur));
                continue;
            }
        }
        cur.push(ch);
    }
    if !cur.trim().is_empty() {
        out.push(cur);
    }
    out
}

// ---------------------------------------------------------------------------
// Directory plumbing
// ---------------------------------------------------------------------------

/// The vector archive root: `ZODA_KZG_VECTORS` or `<workspace>/spec-vectors`.
pub fn vectors_dir() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("ZODA_KZG_VECTORS") {
        let p = PathBuf::from(p);
        if p.is_dir() {
            return Some(p);
        }
        return None;
    }
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../spec-vectors");
    if root.is_dir() {
        Some(root)
    } else {
        None
    }
}

/// Load the mainnet ceremony setup used by the vector suites.
pub fn load_mainnet_setup(dir: &Path) -> Result<Setup, String> {
    let path = match std::env::var("ZODA_TRUSTED_SETUP") {
        Ok(p) => PathBuf::from(p),
        Err(_) => {
            let local = dir.join("trusted_setup.txt");
            if local.is_file() {
                local
            } else {
                return Err(format!(
                    "trusted_setup.txt not found in {} (set ZODA_TRUSTED_SETUP)",
                    dir.display()
                ));
            }
        }
    };
    Setup::load_file(&path)
}

/// Run `f` over every case of a suite directory
/// (`<fork>/kzg/<suite>/kzg-mainnet/<case>/data.yaml`).
pub fn for_each_case<F>(dir: &Path, fork: &str, suite: &str, mut f: F) -> usize
where
    F: FnMut(&str, &Sv),
{
    let suite_dir = dir.join("tests/general").join(fork).join("kzg").join(suite);
    // descend into the configuration directory (kzg-mainnet)
    let Ok(cfg_entries) = std::fs::read_dir(&suite_dir) else {
        return 0;
    };
    let config_dir: Option<PathBuf> = cfg_entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| p.is_dir() && !p.file_name().unwrap_or_default().to_string_lossy().starts_with("._"));
    let Some(config_dir) = config_dir else {
        return 0;
    };
    let Ok(entries) = std::fs::read_dir(&config_dir) else {
        return 0;
    };
    let mut case_dirs: Vec<PathBuf> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.is_dir()
                && !p
                    .file_name()
                    .and_then(|n| n.to_str())
                    .map(|n| n.starts_with("._"))
                    .unwrap_or(false)
        })
        .collect();
    case_dirs.sort();
    let mut count = 0;
    for case_dir in case_dirs {
        let data = case_dir.join("data.yaml");
        let Ok(text) = std::fs::read_to_string(&data) else {
            continue;
        };
        let name = case_dir
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("?")
            .to_string();
        let parsed = parse_simple_yaml(&text);
        f(&name, &parsed);
        count += 1;
    }
    count
}

// ---------------------------------------------------------------------------
// Deneb (EIP-4844) suites — run against this crate's API
// ---------------------------------------------------------------------------

#[cfg(test)]
mod deneb_tests {
    use super::*;
    use crate::eip4844::*;

    fn setup() -> Option<Setup> {
        let dir = vectors_dir()?;
        Some(load_mainnet_setup(&dir).expect("mainnet trusted setup"))
    }

    /// input.blob as bytes (length-checked by the API).
    fn blob_of(case: &Sv) -> Vec<u8> {
        case.get("input")
            .and_then(|i| i.get("blob"))
            .and_then(|b| b.hex())
            .expect("blob hex")
    }

    #[test]
    fn spec_vectors_blob_to_kzg_commitment() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(&dir, "deneb", "blob_to_kzg_commitment", |name, case| {
            let blob = blob_of(case);
            let out = case.get("output").expect("output");
            if out.is_null() {
                assert!(
                    blob_to_kzg_commitment(&blob, &setup).is_err(),
                    "{}: expected an error",
                    name
                );
            } else {
                let expect: [u8; 48] = out.hex_array().expect("commitment hex");
                let got = blob_to_kzg_commitment(&blob, &setup)
                    .unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                assert_eq!(got, expect, "{}", name);
            }
            n += 1;
        });
        assert_eq!(count, n);
        eprintln!("blob_to_kzg_commitment: {} cases passed", n);
    }

    #[test]
    fn spec_vectors_compute_kzg_proof() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(&dir, "deneb", "compute_kzg_proof", |name, case| {
            let blob = blob_of(case);
            let z: Vec<u8> = case
                .get("input")
                .and_then(|i| i.get("z"))
                .and_then(|z| z.hex())
                .expect("z hex");
            let out = case.get("output").expect("output");
            if out.is_null() {
                assert!(
                    compute_kzg_proof(&blob, &z, &setup).is_err(),
                    "{}: expected an error",
                    name
                );
            } else {
                let lst = out.as_list().expect("output list");
                let expect_proof: [u8; 48] = lst[0].hex_array().expect("proof hex");
                let expect_y: [u8; 32] = lst[1].hex_array().expect("y hex");
                let (proof, y) = compute_kzg_proof(&blob, &z, &setup)
                    .unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                assert_eq!(proof, expect_proof, "{}: proof", name);
                assert_eq!(y, expect_y, "{}: y", name);
            }
            n += 1;
        });
        assert_eq!(count, n);
        eprintln!("compute_kzg_proof: {} cases passed", n);
    }

    #[test]
    fn spec_vectors_verify_kzg_proof() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(&dir, "deneb", "verify_kzg_proof", |name, case| {
            let input = case.get("input").expect("input");
            let commitment = input.get("commitment").and_then(|v| v.hex()).expect("c");
            let z = input.get("z").and_then(|v| v.hex()).expect("z");
            let y = input.get("y").and_then(|v| v.hex()).expect("y");
            let proof = input.get("proof").and_then(|v| v.hex()).expect("proof");
            let out = case.get("output").expect("output");
            let got = verify_kzg_proof(&commitment, &z, &y, &proof, &setup);
            match out {
                Sv::Null => assert!(got.is_err(), "{}: expected an error", name),
                Sv::Bool(expect) => {
                    assert!(
                        got.is_ok(),
                        "{}: unexpected error {:?}",
                        name,
                        got.err()
                    );
                    assert_eq!(got.unwrap(), *expect, "{}", name);
                }
                other => panic!("{}: unexpected output {:?}", name, other),
            }
            n += 1;
        });
        assert_eq!(count, n);
        eprintln!("verify_kzg_proof: {} cases passed", n);
    }

    #[test]
    fn spec_vectors_compute_blob_kzg_proof() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(&dir, "deneb", "compute_blob_kzg_proof", |name, case| {
            let blob = blob_of(case);
            let commitment: Vec<u8> = case
                .get("input")
                .and_then(|i| i.get("commitment"))
                .and_then(|v| v.hex())
                .expect("commitment hex");
            let out = case.get("output").expect("output");
            if out.is_null() {
                assert!(
                    compute_blob_kzg_proof(&blob, &commitment, &setup).is_err(),
                    "{}: expected an error",
                    name
                );
            } else {
                let expect: [u8; 48] = out.hex_array().expect("proof hex");
                let got = compute_blob_kzg_proof(&blob, &commitment, &setup)
                    .unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                assert_eq!(got, expect, "{}", name);
            }
            n += 1;
        });
        assert_eq!(count, n);
        eprintln!("compute_blob_kzg_proof: {} cases passed", n);
    }

    #[test]
    fn spec_vectors_verify_blob_kzg_proof() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(&dir, "deneb", "verify_blob_kzg_proof", |name, case| {
            let input = case.get("input").expect("input");
            let blob = input.get("blob").and_then(|v| v.hex()).expect("blob");
            let commitment = input.get("commitment").and_then(|v| v.hex()).expect("c");
            let proof = input.get("proof").and_then(|v| v.hex()).expect("proof");
            let out = case.get("output").expect("output");
            let got = verify_blob_kzg_proof(&blob, &commitment, &proof, &setup);
            match out {
                Sv::Null => assert!(got.is_err(), "{}: expected an error", name),
                Sv::Bool(expect) => {
                    let v = got.unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                    assert_eq!(v, *expect, "{}", name);
                }
                other => panic!("{}: unexpected output {:?}", name, other),
            }
            n += 1;
        });
        assert_eq!(count, n);
        eprintln!("verify_blob_kzg_proof: {} cases passed", n);
    }

    #[test]
    fn spec_vectors_verify_blob_kzg_proof_batch() {
        let Some(setup) = setup() else {
            eprintln!("skipping: consensus-spec-tests archive not present");
            return;
        };
        let dir = vectors_dir().unwrap();
        let mut n = 0;
        let count = for_each_case(
            &dir,
            "deneb",
            "verify_blob_kzg_proof_batch",
            |name, case| {
                let input = case.get("input").expect("input");
                let blobs: Vec<Vec<u8>> = input
                    .get("blobs")
                    .and_then(|v| v.as_list())
                    .expect("blobs list")
                    .iter()
                    .map(|b| b.hex().expect("blob hex"))
                    .collect();
                let commitments: Vec<Vec<u8>> = input
                    .get("commitments")
                    .and_then(|v| v.as_list())
                    .expect("commitments list")
                    .iter()
                    .map(|b| b.hex().expect("commitment hex"))
                    .collect();
                let proofs: Vec<Vec<u8>> = input
                    .get("proofs")
                    .and_then(|v| v.as_list())
                    .expect("proofs list")
                    .iter()
                    .map(|b| b.hex().expect("proof hex"))
                    .collect();
                let out = case.get("output").expect("output");
                let blob_refs: Vec<&[u8]> = blobs.iter().map(|b| b.as_slice()).collect();
                let comm_refs: Vec<&[u8]> = commitments.iter().map(|b| b.as_slice()).collect();
                let proof_refs: Vec<&[u8]> = proofs.iter().map(|b| b.as_slice()).collect();
                let got = verify_blob_kzg_proof_batch(&blob_refs, &comm_refs, &proof_refs, &setup);
                match out {
                    Sv::Null => assert!(got.is_err(), "{}: expected an error", name),
                    Sv::Bool(expect) => {
                        let v = got.unwrap_or_else(|e| panic!("{}: unexpected error {}", name, e));
                        assert_eq!(v, *expect, "{}", name);
                    }
                    other => panic!("{}: unexpected output {:?}", name, other),
                }
                n += 1;
            },
        );
        assert_eq!(count, n);
        eprintln!("verify_blob_kzg_proof_batch: {} cases passed", n);
    }
}
