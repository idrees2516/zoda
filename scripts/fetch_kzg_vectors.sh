#!/usr/bin/env bash
# Fetch the official consensus-spec-tests KZG vector archive and the
# mainnet ceremony trusted setup into ./spec-vectors/, enabling the
# optional spec-vector test suites in zoda-kzg (Deneb / EIP-4844) and
# zoda-edas (Fulu / EIP-7594).
#
# Usage:  scripts/fetch_kzg_vectors.sh [release-tag]
#         (default release: v1.5.0)
#
# The vector tests are skipped unless ./spec-vectors/ exists, so CI and
# local builds stay hermetic by default. Total download ≈ 22 MB.

set -euo pipefail

RELEASE="${1:-v1.5.0}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/spec-vectors"

mkdir -p "$OUT"

echo ">> downloading consensus-spec-tests general package ($RELEASE, ~22 MB)…"
curl -fL --retry 3 -o "$OUT/general.tar.gz" \
    "https://github.com/ethereum/consensus-spec-tests/releases/download/$RELEASE/general.tar.gz"

echo ">> extracting the deneb (EIP-4844) and fulu (EIP-7594) kzg suites…"
tar -xzf "$OUT/general.tar.gz" -C "$OUT" \
    --wildcards \
    "tests/general/deneb/kzg/*" \
    "tests/general/fulu/kzg/*"
rm -f "$OUT/general.tar.gz"

if [ ! -f "$OUT/trusted_setup.txt" ]; then
    echo ">> downloading the mainnet ceremony trusted setup (~800 KB)…"
    curl -fL --retry 3 -o "$OUT/trusted_setup.txt" \
        "https://raw.githubusercontent.com/ethereum/c-kzg-4844/main/src/trusted_setup.txt"
fi

echo ">> done. vectors in $OUT/tests/general/{deneb,fulu}/kzg"
echo ">> run: cargo test -p zoda-kzg -p zoda-edas --release spec_vectors -- --nocapture"
