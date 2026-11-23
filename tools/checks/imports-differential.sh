#!/usr/bin/env bash
# Run with the same toolchain, profile and CATALLAX_THREADS on both revisions.
set -euo pipefail

case ${1-} in
    record)
        output=${2:?usage: imports-differential.sh record OUTPUT_DIR}
        mkdir -p "$output"
        export CATALLAX_IMPORTS_DIFFERENTIAL=1 CATALLAX_THREADS=4
        cargo test --release -j "${CARGO_BUILD_JOBS:-2}" -p catallax-compare --test import_spline -- \
            --nocapture --test-threads=1 > "$output/spline.log" 2>&1 || {
            cat "$output/spline.log" >&2
            exit 1
        }
        cargo test --release -j "${CARGO_BUILD_JOBS:-2}" -p catallax-field --test scan_fit \
            --test scan_p3_round_trip -- --nocapture --test-threads=1 > "$output/scan.log" 2>&1 || {
            cat "$output/scan.log" >&2
            exit 1
        }
        for suite in spline scan; do
            sed -n 's/^.*\(imports-differential .*\)$/\1/p' "$output/$suite.log" |
                LC_ALL=C sort > "$output/$suite.txt"
            test -s "$output/$suite.txt"
        done
        ;;
    compare)
        main=${2:?usage: imports-differential.sh compare MAIN_DIR CANDIDATE_DIR}
        candidate=${3:?usage: imports-differential.sh compare MAIN_DIR CANDIDATE_DIR}
        for suite in spline scan; do
            diff -u "$main/$suite.txt" "$candidate/$suite.txt"
        done
        ;;
    *)
        printf 'usage: %s record OUTPUT_DIR | compare MAIN_DIR CANDIDATE_DIR\n' "$0" >&2
        exit 2
        ;;
esac
