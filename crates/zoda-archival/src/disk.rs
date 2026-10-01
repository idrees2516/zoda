//! Directory-backed custody storage with atomic writes.

use crate::format;
use crate::store::GridCustody;
use std::path::{Path, PathBuf};

/// A directory of self-contained custody files
/// (`<hex(grid_id)>.zoda`), written atomically (temp + rename).
pub struct DiskCustody {
    dir: PathBuf,
}

fn hex_encode(id: &[u8; 32]) -> String {
    let mut s = String::with_capacity(64);
    for b in id {
        s.push_str(&format!("{:02x}", b));
    }
    s
}

fn hex_decode(s: &str) -> Option<[u8; 32]> {
    if s.len() != 64 {
        return None;
    }
    let mut out = [0u8; 32];
    for i in 0..32 {
        out[i] = u8::from_str_radix(&s[2 * i..2 * i + 2], 16).ok()?;
    }
    Some(out)
}

impl DiskCustody {
    /// Open (creating if needed) a custody directory.
    pub fn open(dir: &Path) -> Result<DiskCustody, String> {
        std::fs::create_dir_all(dir).map_err(|e| format!("mkdir {}: {}", dir.display(), e))?;
        Ok(DiskCustody {
            dir: dir.to_path_buf(),
        })
    }

    fn grid_path(&self, id: &[u8; 32]) -> PathBuf {
        self.dir.join(format!("{}.zoda", hex_encode(id)))
    }

    /// Persist a grid custody atomically: write to a unique temp file in
    /// the same directory, fsync, then rename over the destination.
    pub fn save_grid(&self, custody: &GridCustody) -> Result<(), String> {
        let bytes = format::encode_grid(custody);
        let dest = self.grid_path(&custody.grid_id);
        let tmp = dest.with_extension(format!(
            "tmp.{}",
            std::process::id()
        ));
        std::fs::write(&tmp, &bytes).map_err(|e| format!("write {}: {}", tmp.display(), e))?;
        if let Ok(f) = std::fs::File::open(&tmp) {
            let _ = f.sync_all();
        }
        std::fs::rename(&tmp, &dest).map_err(|e| format!("rename: {}", e))?;
        Ok(())
    }

    /// Load and checksum-verify a grid custody file
    /// (`Ok(None)` if not present).
    pub fn load_grid(&self, grid_id: &[u8; 32]) -> Result<Option<GridCustody>, String> {
        let path = self.grid_path(grid_id);
        match std::fs::read(&path) {
            Ok(bytes) => format::decode_grid(&bytes).map(Some),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(format!("read {}: {}", path.display(), e)),
        }
    }

    /// All stored grid ids.
    pub fn list_grids(&self) -> Vec<[u8; 32]> {
        let mut out = Vec::new();
        if let Ok(entries) = std::fs::read_dir(&self.dir) {
            for e in entries.flatten() {
                let name = e.file_name();
                let name = name.to_string_lossy();
                if let Some(id) = name.strip_suffix(".zoda") {
                    if let Some(id) = hex_decode(id) {
                        out.push(id);
                    }
                }
            }
        }
        out.sort_unstable();
        out
    }

    /// Total bytes used by custody files.
    pub fn total_bytes(&self) -> usize {
        self.list_grids()
            .into_iter()
            .filter_map(|id| std::fs::metadata(self.grid_path(&id)).ok())
            .map(|m| m.len() as usize)
            .sum()
    }

    /// Delete one grid's file.
    pub fn delete_grid(&self, grid_id: &[u8; 32]) -> Result<(), String> {
        let path = self.grid_path(grid_id);
        match std::fs::remove_file(&path) {
            Ok(()) => Ok(()),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(e) => Err(format!("delete {}: {}", path.display(), e)),
        }
    }

    /// Enforce a byte budget: evict oldest-slot grids (file modification
    /// is not a proxy for slot order, so the caller passes slot ids via
    /// the loaded files' headers). Returns evicted grid ids.
    pub fn prune_to_budget(&self, budget: usize) -> Result<Vec<[u8; 32]>, String> {
        if budget == 0 {
            return Ok(Vec::new());
        }
        let mut ids = self.list_grids();
        let mut total = self.total_bytes();
        if total <= budget {
            return Ok(Vec::new());
        }
        // oldest slot first
        let mut with_slot: Vec<(u64, [u8; 32])> = Vec::with_capacity(ids.len());
        for id in ids.drain(..) {
            let slot = self
                .load_grid(&id)?
                .map(|g| g.slot)
                .unwrap_or(u64::MAX);
            with_slot.push((slot, id));
        }
        with_slot.sort_unstable();
        let mut evicted = Vec::new();
        for (_, id) in with_slot {
            if total <= budget {
                break;
            }
            if let Ok(meta) = std::fs::metadata(self.grid_path(&id)) {
                total -= meta.len() as usize;
            }
            self.delete_grid(&id)?;
            evicted.push(id);
        }
        Ok(evicted)
    }
}
