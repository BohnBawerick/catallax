# Edge hunt

Generate independent regression questions, run only the kernel, and group its refusals, timeouts and crashes. This does not check whether answered geometry is correct. It does not use an external tester's code or inputs.

Python 3.11 or newer, standard library only, on Linux. Python keeps process scheduling, JSONL recovery and reports out of the kernel workspace without adding Cargo dependencies. Each kernel invocation remains single-threaded. The tool uses a bounded thread pool only to wait for child processes.

## Run locally

Build the release binary from the revision being investigated. Keep generated material outside the checkout.

```sh
cargo build --release -p catallax-cli
RUN=/tmp/edge-hunt-baseline
python3 -B tools/edge-hunt/edge_hunt.py generate --seed 20261020 --count 5004 --out "$RUN/questions"
nice -n 15 python3 -B tools/edge-hunt/edge_hunt.py run \
  "$RUN/questions/catalog.jsonl" --workers 4 --cpu-budget 13500 \
  --results "$RUN/results.jsonl"
python3 -B tools/edge-hunt/edge_hunt.py report "$RUN/results.jsonl" --out "$RUN/report"
```

The default binary is `target/release/catallax`. Pass `--binary /absolute/path/to/catallax` for another build. Each answer runs through `nice -n 15`, including when the runner itself was not niced. The default outer wall limit is 60 seconds. Questions keep the interface's default 50-second internal deadline and request default measurement and STEP export. A smaller `--wall-limit` is available for smoke tests, never above 60 seconds.

`--cpu-budget` limits child CPU seconds for this invocation. Before starting work, the scheduler reserves the full wall limit for each active single-threaded kernel process, plus CPU already charged for finished children by `getrusage`. It stops admitting questions rather than changing their deadlines. Exit 75 means the corpus is incomplete at that budget. Resume with a new explicitly chosen invocation budget; budgets do not silently renew inside one invocation. The final console record includes measured child CPU seconds. If installed, `/usr/bin/time -v` also accounts for the Python process. An interrupted invocation can spend unrecorded CPU on unfinished questions, so retain these logs when enforcing an aggregate budget.

Four workers and 13,500 child CPU seconds leave room below four core-hours for a small smoke run and tool overhead. This is not a guarantee that any given 5,004-case corpus finishes within that budget. Slow cases are the subject of the run, not grounds to shorten their judge-compatible deadlines. No cloud commands are used.

## Coverage and reproducibility

The generator follows `docs/JUDGE-INTERFACE.md`. These strata remain useful regression families, not product acceptance grid cells:

| Dimension | Values |
| --- | --- |
| Pose family | general, tangent, coaxial, equal-radius, along-seam, near-miss |
| Ordered shape pair | all nine pairs of block, cylinder and sphere |
| First operation | unite, subtract, intersect |
| Boolean count | 1, 3, 4, 5 |
| Scale exponent | 0, -30, 30, -15, 15, -6, 6 |

Every 648 questions covers each family, ordered pair, operation and chain length once. Scale changes every 648 questions. Every 4,536 questions covers their complete product. Counts divisible by six give equal family populations. Chains use new operands and cycle operations and shapes at each step, rather than repeatedly aliasing the same bodies.

Each case has an independent RNG derived from SHA-256 of the seed and index. Changing scheduling or generating a longer corpus leaves the existing prefix unchanged. Integer radii vary between 1 and 3 before scaling. All geometry is serialized as exact rational strings, never binary floats.

The centered shapes have support radius `r` in x. Tangent cases place their centers exactly `r1+r2` apart in x. Near misses add a positive exact gap of `10^-6`, `10^-15` or `10^-30` before scaling. Coaxial cases share the z axis and have an axial offset. For blocks this means aligned center axes, not a periodic surface axis. Equal-radius cases give cylinders and spheres the same radius and blocks the same half-size. Along-seam cases put the second shape's y-minimum on the first shape's local `y=0` seam plane, with positive-x contact/intersection. Blocks use that same reference plane even though they have no periodic seam.

General cases overlap at rational center offsets and can have relative rotations. All families also receive shared rigid rotations, preserving exact tangency and near misses. Quaternions include identity, `[1,1,1,1]` and `[1/2,1/3,1/5,1/7]`, which gives a non-axis-aligned exact rational matrix. Only general cases rotate the operands relative to each other. Chain operands need not preserve the first pair's named relation; their actual poses remain in the question.

Question files contain only `exact-brep-question/2` fields. The format explicitly rejects family labels. `catalog.jsonl` is the sidecar with one object per question:

- `seed`, `index`, `family` and `parameters`, including initial shape pair, first operation, chain length, scale, radii, gap and quaternions;
- `file`, a relative filename for direct kernel use;
- the complete `question` and its canonical JSON `question_sha256`.

The catalog's embedded question is the runner's authority. Editing the standalone file alone does not change a run. To investigate a hand-edited reproducer, run it directly with `catallax answer` or generate a new catalog with matching hashes. Use a new output directory when generating another corpus. Existing catalogs are never overwritten.

## Resume and results

Repeat the same `run` command to resume. Completed question IDs are skipped only if their question hash, binary SHA-256, diagnostic mode and wall limit match. The result file is exclusively locked, preventing competing writers. Resume repairs a torn final JSONL append, but rejects complete malformed records, duplicate IDs and hash mismatches. SIGINT or SIGTERM kills only this runner's active process groups, leaves unfinished questions absent from results, and exits 130. A hard kill may leave children alive until their internal deadline; use graceful interruption.

Every completed question produces one flushed, fsynced JSONL object with format `catallax-edge-hunt-result/1`. It embeds the catalog entry and records:

- `status`: answered, refused, invalid-question, timeout or crash;
- manifest status, refusal code, stage, step and free-text detail;
- the entire `stats` block, including kernel `wall_ms`, `decisions`, `max_bits`, `work` and `refused`, or null when unavailable;
- measured subprocess `wall_ms`, exit code, up to 16 KiB of stderr and any manifest-reading error;
- `cpu_ms` from this exact child's `wait4` user plus system time, and `cpu_wall_ratio`, CPU time divided by subprocess wall time;
- the full manifest, binary SHA-256, outer wall limit and diagnostic mode.

Nonzero exits are crashes even with a manifest. A missing, malformed or wrongly identified manifest is also a crash. Reaching the outer limit is a timeout and records the killed process's negative signal return code. A kernel `numbers.deadline` refusal stays a refusal, not a runner timeout. Invalid generated questions are grouped too, so a generator defect cannot disappear from the report.

CPU attribution uses `wait4` on each worker's own PID, not differences of process-wide counters while other workers finish. A low CPU/wall ratio on a deadline indicates that the process spent much of its wall budget off-CPU. A ratio near one with an imprint deadline points toward computation rather than scheduler contention, but is not itself a correctness diagnosis. CPU samples are measurement only and never enter outcome comparison.

`--diagnostics` enables Rust panic backtraces and extracts the first kernel source location printed by the panic hook. Ordinary runs set `RUST_BACKTRACE=0`; no profiler, stack sampling or kernel instrumentation runs. Refusals generally have no source location because the public manifest exposes no such field. The tool reports null instead of inferring a location from the refusal text. Diagnostic and nondiagnostic result sets cannot be resumed or compared as equivalent runs.

## Groups and reproducers

`report` writes `report.md`, `summary.json` and up to three exact question files per cause under `examples/<group-id>/`. It starts with code/stage rollups, each linked to its finer groups. These rollups retain failure count, share of failures, total wall time, observed CPU time and aggregate CPU/wall ratio. Fine keys are operation, ordered shape pair, family, refusal code, stage and optional source location. The JSON ranks fine groups by count and lists the independent total-time ordering; Markdown gives both ranks.

For a named failing boolean step, the key uses that step's operation and operand types. An earlier result is labeled `result(block,cylinder,...)`, preserving its shape lineage without pretending it is an untouched primitive. A refusal after the boolean steps or an outer timeout has no known current step. Its key explicitly says `unknown-step:<operation sequence>` and uses the initial shape pair. Reports never guess the failing operation from partial timing data.

Deadline groups, both `numbers.deadline` and outer wall limits, also split by `identity` versus `rotated` quaternion poses and scale band. A case is rotated if any operand has a nonidentity rotation, including later chain operands. Bands use the catalog's uniform scale exponent: tiny at or below -15, small from -14 to -1, unit at 0, large from 1 to 14, and huge at or above 15. The full exponent and quaternion remain in each case's parameters.

Older result lines without `cpu_ms` remain readable. Their CPU is unknown, not zero. Reports state the number of CPU observations and compute aggregate ratios only over those observations; they never estimate per-question CPU from the old invocation-wide total.

Each group retains the smallest observed parameter sizes, breaking ties by question hash. Size is the sum of exact geometric numerator and denominator bit lengths plus the step count. This is a cheap selection rule, not delta debugging. Every example is an unmodified format-2 question, ready for:

```sh
target/release/catallax answer /path/to/report/examples/<group>/<question>.json --out /tmp/edge-hunt-repro
```

Results embed all geometry, so reports remain reproducible after the original generated question directory is removed. Generate each report into its own directory to avoid retaining stale examples from a previous report.

## Compare a fix

Run the unchanged catalog with the new binary into a new result file. Never resume old results with a different binary.

```sh
python3 -B tools/edge-hunt/edge_hunt.py run "$RUN/questions/catalog.jsonl" \
  --binary target/release/catallax --workers 4 --cpu-budget 13500 \
  --results /tmp/edge-hunt-after.jsonl
python3 -B tools/edge-hunt/edge_hunt.py compare "$RUN/results.jsonl" \
  /tmp/edge-hunt-after.jsonl --out /tmp/edge-hunt-comparison
```

Comparison requires identical complete question sets, hashes, metadata, wall limits and diagnostic modes. It writes `comparison.md` and `comparison.json`, with every cause marked vanished, shrank, appeared, grew or unchanged. Every question with a changed outcome is listed with its full before/after stable outcome. Kernel stats and build identity are excluded; answered topology, properties, marks, export bounds, exactness and refusals are included. Thus an answered-to-answered output change is visible too. STEP file bytes are not retained or compared, only the manifest's export description. This establishes which reported outcomes changed, not whether new geometry is mathematically correct. It cannot replace independent correctness verification.

## Small tests

```sh
python3 -B -m unittest discover -s tools/edge-hunt -v
```

Tests cover geometric contact/gap invariants, strata, chain references, torn writes, resume identity, crash/timeout classification, step attribution, smallest reproducers and comparison regressions. They do not change kernel source or kernel tests. Runtime is a few seconds, including one deliberate one-second timeout.
