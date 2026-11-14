> Historical B02 implementation evidence at `32e460b40e6918b0c42d1666446217cd7ffb67b6`.
> For the parameter-identity repair and current focused results, see [P2/VALIDATION.md](P2/VALIDATION.md).
> Original capture bytes are preserved; earlier results do not qualify the new guard by themselves.

# B02 executed evidence

Baseline: `3a2f4e0c65feae548f4980d5c520a0c330acf6fb`.
Branch: `phase1/B02-seam-start`.
Date: 2026-11-14 UTC.

Commands are portable substitutions of the executed saved-environment commands. Set `CARGO_HOME_DIR` and `RUSTUP_HOME_DIR` to existing authorized tool roots; run from the repository root. No toolchain probe was needed.

```sh
CARGO_HOME="$CARGO_HOME_DIR" RUSTUP_HOME="$RUSTUP_HOME_DIR" "$CARGO_HOME_DIR/bin/cargo" test --offline -p catallax-trace --test b02_seam -- --nocapture
CARGO_HOME="$CARGO_HOME_DIR" RUSTUP_HOME="$RUSTUP_HOME_DIR" "$CARGO_HOME_DIR/bin/cargo" test --offline -p catallax-trace --test intersector --test drop_curve --test curve_surface
```

Final focused run: exit 0; 3 tests passed, 0 failed, in 10.31 seconds. Observed worker summaries:

```text
B02 workers=1 cases=128 partitions=[128] complete_outcomes=identical; application_work=unmetered
B02 workers=2 cases=128 partitions=[64, 64] complete_outcomes=identical; application_work=unmetered
B02 workers=4 cases=128 partitions=[32, 32, 32, 32] complete_outcomes=identical; application_work=unmetered
```

Existing integration run: exit 0. `curve_surface`: 6 passed in 11.89 seconds; `drop_curve`: 3 passed in 5.51 seconds; `intersector`: 6 passed in 75.43 seconds. Those suites ran after the production repair; subsequent production edits were rustfmt only. Saved fixtures and comparison allowances are unchanged. Times are observations on a shared cloud executor, not benchmarks.

Baseline regression: the translated/clipped regression test alone was run on the frozen baseline `surface_trace.rs`, exit 101 with `NotPorted(Intersect)` on the first full-period case. Candidate bytes were restored unconditionally, and the full new suite passed afterward. The baseline log is a publication-sanitized copy; original bytes are retained privately. The final and existing-suite summaries here are transcribed from observed tool output, not represented as raw stdout files.

Static validation: changed Rust files formatted; local relative links resolve; public outputs contain no machine paths; changed-artifact SHA-256 digests verified; `git diff --check` passes. Only the B02 source, regression, documentation and evidence are included. Disposable restore-copy scratch was removed; baseline raw evidence and reusable build cache are retained.

Limits: no measured application work, request-level allowance validation, production parallelism, donor-oracle regeneration, application Boolean/STEP replay, all-thread qualification or geometric acceptance certificate. No security, account or permission changes; no denied operations, hosted gates, CI launch, PR creation, merge or deployment.
