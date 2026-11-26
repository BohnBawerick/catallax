#!/usr/bin/env bash
# Runs inside cx-offload: builds the step 5.1 binaries and keeps them with a source snapshot.
# usage: build-all.sh [census-only]
set -euo pipefail
K=$HOME/offload/keep/cx-p5-speed-owner
mkdir -p "$K/bin" "$K/logs" "$K/census"
export CARGO_PROFILE_RELEASE_DEBUG=line-tables-only
names="examples/opened_import examples/predictor_corpus"
cargo build --release -p catallax-import --example opened_import
if [ "${1:-}" = census-only ]; then
  cargo build --release -p catallax-cli --example predictor_corpus
else
  cargo build --release -p catallax-cli --bins --example predictor_corpus --example work_profile --example app_work_profile --example work_span
  names="$names examples/work_profile examples/app_work_profile examples/work_span catallax-verify"
fi
for name in $names; do
  cp "$CARGO_TARGET_DIR/release/$name" "$K/bin/$(basename "$name")"
done
( cd "$K/bin" && for name in $names; do shasum -a 256 "$(basename "$name")"; done > SHA256SUMS )
git rev-parse HEAD > "$K/bin/REVISION"
rm -rf "${K:?}/src"
mkdir -p "$K/src"
git archive HEAD | tar -x -C "$K/src"
echo "build-all done at $(git rev-parse HEAD)"
