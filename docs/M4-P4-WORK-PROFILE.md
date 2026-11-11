# P4 work profile and allowance verdict

P4 still refuses at the unchanged 1M scoreboard allowance. At source
`0013f75c4078e0fb0b06b02946d456e863858b7a`, the complete verification request
uses **1,733,492 work units**. Exact construction alone uses 1,554,975. The
eager application build with display uses 1,819,042 and fits its separate 20M
allowance. M4 remains open.

This is the count and verdict landing for [PLAN item 17](PLAN.md#114-the-22-october-next-items-on-11-november).
The input is the unchanged [canonical P4 recipe](../parts/P4-shelled-tray.json).
[Portable evidence](evidence/m4-p4-work-profile.json) records the source and
binary hashes, request outcomes, operation counts, stage attribution and fresh
thread comparisons. Raw evidence retained privately.

## Allowances and operation counts

| Entry | Existing allowance | Work used | Result |
| --- | ---: | ---: | --- |
| Verification, recipe through STEP | 1,000,000 | 999,995 | `numbers.operation-budget` in the rounding boolean |
| Same verification request | 20,000,000 | 1,733,492 | Answered, exact, 15 faces |
| Same verification request | 100,000,000 | 1,733,492 | Same manifest and STEP bytes as at 20M |
| Application, exact body and topology check | 20,000,000 | 1,554,975 | Built, 15 faces |
| Application, eager build with 0.01 mm display | 20,000,000 | 1,819,042 | Built, 194 display vertices and 384 triangles |

The 100M request measures completed work under the existing diagnostic
allowance. It does not turn the 1M refusal into a scoreboard pass. The eager
display result does not qualify the separate unproved-preview lifecycle.

| Verification operation | At 1M | At 20M and 100M |
| --- | ---: | ---: |
| `blank`, block | 0 | 0 |
| `tray`, shell with top open | 63,289 | 63,289 |
| `rounding`, cylinder | 0 | 0 |
| `rounded`, intersection | 936,706, refused | 1,491,686 |
| Result summary | Unreached | 32 |
| STEP export | Unreached | 178,485 |
| Total | 999,995 | 1,733,492 |

The shell completes before the failing operation. The failed request leaves no
STEP artifact. The answered requests produce the same STEP SHA-256:
`f9a8b9a15a576a4989690f75e1bd9f697e052299b80b33c1bfb4fdc4634c7381`.

## Where the work goes

The stage rows below include every boolean in the recipe, including the shell's
cavity construction. Each column partitions the same request total as the
operation table above. Other construction is work outside those boolean stages.
The classify and sew work at 1M belongs to shell construction; the final
rounding boolean reaches imprinting and then refuses.

| Stage | At 1M | Completed request at 100M |
| --- | ---: | ---: |
| Intersect | 462,478 | 462,478 |
| Imprint | 527,653 | 879,723 |
| Classify | 5,721 | 208,631 |
| Sew | 309 | 309 |
| Other construction | 3,834 | 3,834 |
| Summary | 0 | 32 |
| STEP export | 0 | 178,485 |
| Total | 999,995 | 1,733,492 |

Imprinting accounts for 50.7% of the completed request. Its largest exclusive
regions are probe depths, 493,114 units, and Bernstein barrier signs, 244,995.
Those counts are included in the imprint row. The most expensive numeric caller
is `catallax-boolean/src/probe_path.rs:31`, with 388,728 units in separated
enclosures. Circular-chart evaluation at `catallax-intersect/src/chart.rs:205`
spends another 180,471 units in enclosures.

| Numeric decision class, completed request | Work | Share |
| --- | ---: | ---: |
| Separated and other enclosures | 944,665 | 54.5% |
| Enclosures still containing zero | 87,863 | 5.1% |
| Nonzero signs | 291,465 | 16.8% |
| Exact-zero signs | 224,988 | 13.0% |
| Root isolation | 110,435 | 6.4% |
| Cancellation checkpoints | 74,076 | 4.3% |
| Exact reductions | 0 | 0.0% |

Both observers attribute all request work. Their manifests, work counts and
STEP bytes match the verification binary with the observers inactive. The
application profile attributes the extra 264,067 units above exact construction
to display work outside the recipe operations.

## Verdict and next work

The 1M capability gate fails. Even eliminating all 63,289 units of shell work
would leave 1,491,686 units of exact construction in the final boolean, before
summary or export. A shell-local saving cannot close the gap. The next cost
reduction must reach that boolean's intersection and imprint proofs, especially
probe-depth and barrier work, while preserving their witnesses and decisions.
More threads preserve the answers but do not reduce the deterministic charge.

The analytic curved-face offset and revolved-shell gaps of item 17 remain
separate landings. This report adds no curved-face capability or spline offset
route. Kernel code, tolerances, allowances and saved replay records are unchanged.

## Reproduce and check

Build the existing development observers and verifier:

```sh
cargo build --release -p catallax-cli --bin catallax-verify \
  --features work-span --example work_profile --example work_span \
  --example app_work_profile -j 2
```

Use an external evidence directory for `PROOF`. Run doctor, then repeat the
verification and both profiles at 1M and 100M. Also run verification at 20M.

```sh
DRIVER=.claude/skills/verify-catallax/scripts/control-catallax.py
INPUT=parts/P4-shelled-tray.json
python3 -B "$DRIVER" doctor --binary target/release/catallax-verify --json
CATALLAX_THREADS=1 python3 -B "$DRIVER" run "$INPUT" \
  --binary target/release/catallax-verify --max-work 1000000 \
  --workers 1 --fresh --out "$PROOF/p4-1m" --json
CATALLAX_THREADS=1 target/release/examples/work_profile "$INPUT" \
  --out "$PROOF/profile-1m" --max-work 1000000
CATALLAX_THREADS=1 target/release/examples/work_span "$INPUT" \
  --out "$PROOF/span-1m" --max-work 1000000
CATALLAX_THREADS=1 target/release/examples/app_work_profile "$INPUT" \
  --max-work 20000000 --exact-only
CATALLAX_THREADS=1 target/release/examples/app_work_profile "$INPUT" \
  --max-work 20000000
```

The recorded fresh matrix has four requests per run: P4 at its three existing
allowances, plus the unchanged `parts/P0/near-miss.json` control at 1M. Exported
`catallax-question/1` wrappers retain each request's own allowance. The control
lets all four persistent workers start. It answers at 7,076 units in every run.

| Persistent workers | Kernel threads per worker | Result |
| ---: | ---: | --- |
| 1 | 1 | Baseline |
| 2 | 2 | Equal |
| 4 | 4 | Equal |
| 1 | 12, all available | Equal |

All workers exited zero, with no cache hits. The driver's `compare` checks
ordered inputs, normalized manifests, deterministic work and strict STEP bytes.
Operation logs also agree after removing clocks. Refused requests have no STEP
file at every thread count. The two profilers separately reproduce each P4
verification outcome, work count, normalized manifest and STEP artifact.

This is a work-count observation, not a wall-clock qualification. SolidDraft,
the unproved-preview lifecycle, independent OpenCASCADE import, field comparison
and the protected 477-case replay were not run for this documentation-only
landing. The preserved replay records are untouched. These remaining gates
prevent a claim that P4 or M4 is complete.
