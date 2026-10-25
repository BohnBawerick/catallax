# Interactive latency

Measure latency separately from geometric correctness and retain slow cases for the kernel lanes.

## Sub-features

- `latency.pairs`: pinned, alternating AB/BA measurements.
- `latency.legacy`: original process-per-question elapsed replay.
- `latency.report`: target violations grouped by operation, shapes and outcome.

## How to get to it (user POV)

Driver commands `pair`, `legacy` and `report`.

## Driving it with control-catallax

Preconditions: doctor identified an allowed `CORE`, both stream binaries exist, external `PROOF` directory.

- Pair: `control-catallax.py pair "$INPUT" --baseline "$BASE" --candidate "$CANDIDATE" --cores "$CORE" --out "$PROOF/pairs" --json`. Require one warmup and three measured pairs per case in `pairs.json`, with `load_average_before` and `load_average_after` on each pair.
- Legacy: `control-catallax.py legacy "$INPUT" --binary "$LEGACY" --cores "$CORE" --out "$PROOF/legacy" --json`. Capture total elapsed time and each exit status.
- Report: `control-catallax.py report "$PROOF/one/run.json" --json`. Inspect `outcomes` separately from `slow`. Targets are 100 ms stock, 500 ms boolean, 1,000 ms recipe and 2,000 ms full comparison.

### B-rep edit display

The application benchmark keeps timing outside the code that builds a part:

```sh
cargo build --release -p catallax-app --bench display -j 2 --message-format=json
```

Use the `executable` from Cargo's `compiler-artifact` record for the `display` bench and record its hash. A filename glob can select another build configuration. Pass a Catallax or SolidDraft recipe JSON path, an external output directory, and a positive repetition count. Optional fourth and fifth arguments request a chord tolerance in millimetres and angular deviation in radians. Without them it uses the unchanged `MeshOptions::default()`. Pin the executable to one allowed core with `taskset -c "$CORE"` and use `nice -n 10` on a shared machine. Default mode builds and checks the exact body, then times `mesh_body` separately. `samples.json` records build/display milliseconds, vertex and triangle counts, boundary segment counts, requested controls, per-face bounds and vertex allowance. `display.obj` and `mesh.txt` retain the final mesh, including exact reported bounds and edge polylines.

The chord argument `drawing` derives `1e-5 * max(bounding_diagonal_mm, 1)` from a separately measured native bounding box. Bounding-box work is outside the mesh timer. This mode is not available with `--incremental`, `--application` or `--edit-from`. Public fine-drawing reproductions live in `crates/catallax-app/benches/fine_mesh`; the chamfered collar uses `parts/P2-collar.json`. Use the handoff's explicit chord values when reproducing an existing application request rather than silently substituting the bbox-derived mode.

For edits use the application regression's bracket flange heights 60 and 72, and P3's selected top-edge fillet radii 1 and 3/2. Preserve baseline and candidate executables built with identical release settings and Cargo feature sets. Alternate AB/BA order with one warmup and at least three measured pairs per input; record affinity, executable hashes and load averages. Compare both OBJ and complete mesh bytes, not only counts. Build timings remain separate and must not be described as display latency.

Put `--incremental` before the input to time every body-producing prefix through
`catallax_app::build_incremental`. Put `--application` there instead to time one
final-result `catallax_app::build`. Both use the unchanged 20M application work
allowance and report combined `edit_ms` or `app_ms`. Comparing these modes changes
both intermediate display requests and numeric request scopes. It measures a
scheduling counterfactual, not an implemented SolidDraft optimization or a
complete click. Compare the final OBJ bytes with identical explicit display controls.

For a real suffix edit, put `--edit-from "$BEFORE_NATIVE"` before the new input.
The benchmark constructs the old checkpoint, then times one bounded
`build_incremental` call for the new recipe, including its retained-checkpoint
clone and display. Each iteration separately verifies that its OBJ equals a
fresh new build and that replaying the retained old checkpoint reproduces the
old OBJ. Those checks and the old build are outside `edit_ms`; process-wide CPU
time includes them and must not be labelled edit CPU time. This is a kernel API
measurement, not a complete UI click.


### Recipe Boolean stages

Build `catallax-features` with `--bench rebuild --release --message-format=json`
and select that bench's executable from Cargo's artifact record. For example:

```sh
"$REBUILD" --input "$INPUT" --repeat 3 --stages --step-output "$PROOF/answer.step"
```

The input is a native Catallax recipe. Whole-build rows include executed
operations and `decision_work`; the optional stage pass independently rebuilds
each prefix and times intersection, imprint, classification, selection and
sewing. Do not add these alternate passes together. This raw profile has no
application aggregate work scope; use the bounded application modes above for
that path.


## Gotchas

Use identical inputs and allowances. Record ambient load and affinity. Cached wall times cannot establish a fresh latency result. A faster refusal is not evidence that geometry became faster. Group labels suggest profiling targets; they do not identify measured hot functions.
