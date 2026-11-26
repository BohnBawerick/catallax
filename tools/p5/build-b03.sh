#!/usr/bin/env bash
# Runs on bohn-03 inside its checkout: builds the corpus runner and the verifier.
set -euo pipefail
export CARGO_PROFILE_RELEASE_DEBUG=line-tables-only CARGO_INCREMENTAL=0
cargo build --release -j 3 -p catallax-cli --bins --example predictor_corpus
echo "bohn-03 build done at $(git rev-parse HEAD)"
