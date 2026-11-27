#!/usr/bin/env bash
set -euo pipefail

# Run on the build station with a lane-specific CARGO_TARGET_DIR. The source
# directory and selection are provided on that station, not copied into Git.
exec cargo run --release -p catallax-cli --example predictor_corpus -- census-preview "$@"
