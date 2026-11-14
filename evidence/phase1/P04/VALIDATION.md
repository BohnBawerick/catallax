# P04 executed validation

Baseline: `f4db345f2e7c51835d1ebeba9f6446d117144117`, 2026-11-14 UTC.
Worktree: `/workspace/p04-catallax`; branch `phase1/P04-port-consumers`.
Compiler: `rustc 1.100.0 (0ef1bc386 2026-11-09)`.

## Rust saved-gate subset

Executed successfully, exit 0:

```sh
CARGO_HOME=/workspace/.cargo RUSTUP_HOME=/workspace/.rustup /workspace/.cargo/bin/cargo test --offline -p catallax-trace --test drop_curve --test intersector --test curve_surface --test surface_drop --test curve_on_surface
```

| Integration target | Passed | Failed | Test harness elapsed seconds |
| --- | ---: | ---: | ---: |
| curve_on_surface | 1 | 0 | 0.39 |
| curve_surface | 6 | 0 | 12.78 |
| drop_curve | 3 | 0 | 5.64 |
| intersector | 6 | 0 | 75.31 |
| surface_drop | 1 | 0 | 0.84 |

This is a transcription of observed tool output, not a retained raw stdout log. Total 17 tests passed. The harnesses include their existing 1/2/4-worker checks where coded; no blanket all-thread qualification is inferred. The drop gate still expects two `drefuse` rows. No saved expectations or binary64 allowances were edited. No CLI/Boolean request, STEP, external import, local oracle, sanitizer or memory-fault probe was run. Times are harness observations on a shared cloud machine, not performance acceptance.

Initial launch diagnostics: plain `cargo` was absent from PATH (exit 127). Calling `/workspace/.cargo/bin/cargo` without the saved-environment roots failed to create `/home/agent/.rustup` (exit 1, read-only filesystem). Setting the existing authorized `/workspace/.cargo` and `/workspace/.rustup` roots allowed the offline command above. No permissions changed or escalation requested. Cargo emitted the pre-existing CLI dual-bin manifest warning; it compiled `libm` and `catallax-trace` locally.

## Static and baseline checks

`git ls-remote origin refs/heads/main` returned the baseline SHA before branch creation and again after inspection. `git fetch origin main` updated origin/main to that SHA. First-parent history and tracked source/fixture files support the landed-step table. Static row counts and hashes are recorded in `source-manifest.json`; hashes establish inspected input identity, not donor correctness. `git diff --check` passed before committing. Only P04 documentation/evidence paths are included in the deliverable commit.

No hosted gate state was read or retried. No manual CI was launched. No PR was created or merged. Publication head verification, if successful, is reported in the handoff, not represented as qualification here.
