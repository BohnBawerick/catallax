#!/usr/bin/env bash
# Runs inside cx-offload: builds the census observer and keeps it with a source snapshot.
set -euo pipefail
K=$HOME/offload/keep/cx-p5-speed-owner
mkdir -p "$K/bin" "$K/logs" "$K/census"
cargo build --release -p catallax-import --example opened_import
cp "$CARGO_TARGET_DIR/release/examples/opened_import" "$K/bin/opened_import"
sha256sum "$K/bin/opened_import" > "$K/bin/opened_import.sha256"
git rev-parse HEAD > "$K/bin/opened_import.rev"
rm -rf "${K:?}/src"
mkdir -p "$K/src"
git archive HEAD tools docs/research | tar -x -C "$K/src"
echo "census build done at $(git rev-parse HEAD)"
