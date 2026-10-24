---
name: verify-catallax
description: Drive Catallax's deterministic CLI to verify questions, recipe builds and two-kernel comparisons, replay regressions, compare strict STEP output and measure interactive latency.
---

# Verify Catallax

Read [the feature map](references/features/README.md), then choose the relevant entry. The protocol is [VERIFY-INTERFACE.md](../../../docs/VERIFY-INTERFACE.md).

## Launch

From the repository root, run `cargo build --release -p catallax-cli --bins -j 2`. Use `nice -n 10` on a shared machine. The binary is `target/release/catallax-verify`. The executable driver is `.claude/skills/verify-catallax/scripts/control-catallax.py`; run it with `--help` and a subcommand's `--help` for the complete command reference. Python 3.11 or newer and Linux `taskset` are required for the documented timing workflow.

There is no server. `run` starts its own stream workers and closes them when done. Every worker handles several requests without process startup between them. It owns only the processes it starts. Ctrl-C stops those workers; incomplete requests can be resumed. Never share a run output directory between concurrent commands. Cache directories may be shared.

## Doctor

Run `.claude/skills/verify-catallax/scripts/control-catallax.py doctor --binary target/release/catallax-verify --json`. Require `about.request` to equal `catallax-request/1`. Save the binary hash and select a core from `allowed_cores`. Doctor reads identity and affinity only. It does not create evidence or change the kernel.

## Drive

Choose an external directory for `PROOF`, then run:

```sh
DRIVER=.claude/skills/verify-catallax/scripts/control-catallax.py
EXAMPLES=.claude/skills/verify-catallax/references/examples
"$DRIVER" run "$EXAMPLES/block.json" --binary target/release/catallax-verify --out "$PROOF/block" --json
"$DRIVER" run "$EXAMPLES/recipe.json" --source recipe --binary target/release/catallax-verify --out "$PROOF/recipe" --json
"$DRIVER" run "$EXAMPLES/recipe.json" --source compare --binary target/release/catallax-verify --out "$PROOF/compare" --json
```

For private corpora pass paths at runtime. JSONL corpus rows may wrap the question under `question`. Directories recursively contribute question and recipe files; saved manifests and pack indexes are ignored. Use `--workers 1`, `2` and `4` with separate `--out` directories and `--fresh` for independent determinism proofs. Use `--cache "$PROOF/cache"` for repeated regression checks. Reusing the same `--out` resumes completed verified artifacts.

Compare main and candidate builds using `compare "$PROOF/main/run.json" "$PROOF/candidate/run.json" --json`. It checks ordered inputs, manifests excluding stats and build identity, and strict STEP bytes. Compare same-build runs to check deterministic work too. An internal error fails comparison even when it occurs in both runs.

For timing use `pair INPUT --baseline BASE_BINARY --candidate CANDIDATE_BINARY --cores CORE --out "$PROOF/pairs" --json`. This does not read caches. It alternates order on one core, with a warmup and three measured pairs by default. Every pair records the one-, five- and fifteen-minute load averages at its start and end. Use `legacy` to time the old process-per-question CLI. Do not describe a shared machine as idle without evidence.

Each completed request includes `operation_log` and writes `operations.json` plus a portable `question.json`. Replay the latter with `catallax-verify replay QUESTION --out DIR`, or pass it to the driver's `run`. Both preserve its original work allowance. Compare manifests, STEP bytes and work; exclude operation `wall_ms` when comparing logs.

Run `"$DRIVER" parts parts --binary target/release/catallax-verify --out "$PROOF/parts" --max-work 1000000 --json` for the target scoreboard. It runs B-rep, field and sampled comparison independently, including unsupported recipes. `parts/edit-checks.json` adds complete P0/seam and P1/P2/P4 parameter rebuilds and face/edge name resolution; inspect `edit_check.status` and its feature waits. Field registry lookup is not a boundary certificate. For real STEP imports pass `--occt-python` pointing at a scratch environment with `cadquery-ocp`. Inspect every refusal and `step_check`; missing bindings and absent exports are not passing checks. A passing import does not prove volume deviation, SolidDraft interaction or complete target-part acceptance. See `parts/README.md` and the protocol for the recorded dimensions and remaining gates.

Configured edits checkpoint in top-level `edit_checks` before display work; finished geometry rows also carry `edit_check`. Inspect report `status` before claiming the full scoreboard completed. Budget-limited checks retain their first `blocked_operation`.

For the released-tool comparison, use the [head-to-head rerun guide](scripts/head-to-head/README.md). Build the opt-in `head_to_head` example, download the checksum-pinned release outside the repository, and run `head-to-head parts --binary target/release/examples/head_to_head --reference RELEASE_BINARY --core CORE --out DIR`. This benchmark separates construction, edit, display and requested measurement without instrumenting application builds. Read raw correctness coverage before making any speed claim.

## Evidence

Keep proof outside the repository. Capture the command, doctor output, exit code, `config.json`, `inputs.json`, `run.json`, per-case request/response/manifest/STEP files and matched-pair samples. Check the files themselves, including STEP presence or removal. A response alone does not prove the artifact side effect.

Exercise the real CLI through the driver. Do not substitute internal setters or mock geometry. `report RUN_JSON --json` separates correctness outcomes from slow cases and groups slow cases by operation, shapes and refusal cause. Cached timing is labeled. A cache hit proves artifact reuse, not fresh latency or independent determinism. A cancelled process, missing response or panic cannot prove a geometric refusal. Two-kernel agreement is sampled evidence, not a global distance certificate.

## Cleanup

Run `"$DRIVER" cleanup --out "$PROOF" --dry-run --json`, then `"$DRIVER" cleanup --out "$PROOF" --json`. Cleanup reports the policy and preserves evidence; normal completion and interruption already close owned worker processes. It never kills by process name or deletes proof files. Confirm `"$PROOF/block/run.json"` and its STEP file still exist after cleanup. A failed run leaves complete artifacts available for resume.

## Helpers

`control-catallax.py --help` is the canonical interface. Built-in source adapters cover question/2, recipe/3 and recipe/4, opt-in certified-recipe/1, explicit-accuracy `catallax-print/1` field STL/OBJ output, independent field compilation, sampled comparison, opt-in `name-edit` wrappers and portable question/1 replay. Print requests retain the field timing target; strict comparison includes `print.stl` bytes. A trusted custom Python module passed with `--adapter` exports `load(Path)` yielding source/problem objects; it cannot change the kernel's exchange semantics.

`catallax-stand-in/1` wrappers use source `recipe`. They request field moment enclosures and an explicit massless-envelope/unit-density STEP stand-in, independently of print meshes. Check `stand_in.field_properties`, `surrogate_properties`, `exported_property_bounds` and `within_field_enclosures`. On independent import, sum only solid mass properties; a closed envelope shell must remain surface geometry. Density/material metadata and density-weighted option B are not implemented.

Run `python3 -B -m unittest discover -s .claude/skills/verify-catallax/scripts -p 'test_*.py' -v` after building to verify the driver against real commands. Set `CATALLAX_VERIFY_BINARY` to test a preserved executable.

Use `/maintain-verification-skill` when commands or feature coverage change.
