# Target-part head-to-head

This is scoreboard tooling. The Rust example and Python driver do not change a kernel or application build path. SolidCraft scripts were authored from the fixed dimensions in `parts/README.md`. No SolidCraft source is copied or vendored. Only the released executable is downloaded into an external directory.

## Rerun

Run from the repository root. Keep executables, downloaded archives and output outside the repository. Use an allowed core from the verification doctor's `allowed_cores` or `os.sched_getaffinity(0)`.

```sh
export CARGO_BUILD_JOBS=2
nice -n 10 cargo build --release -p catallax-cli --example head_to_head
DRIVER=.claude/skills/verify-catallax/scripts/control-catallax.py
python3 -B "$DRIVER" head-to-head-reference --out /tmp/h2h/reference --json
uv venv /tmp/h2h/occt
uv pip install --python /tmp/h2h/occt/bin/python cadquery-ocp==8.0.1.1.0
nice -n 10 python3 -B "$DRIVER" head-to-head parts \
  --binary target/release/examples/head_to_head \
  --reference /tmp/h2h/reference/solidcraft-cli \
  --occt-python /tmp/h2h/occt/bin/python \
  --core 0 --max-work 1000000 --out /tmp/h2h/run --json
```

`reference.json` pins the versioned Linux x86-64 release URL, archive SHA-256 and executable SHA-256. Download verifies the archive before extracting only its named executable, then verifies the executable. Each benchmark checks the executable again before starting. Other platforms require a reviewed pin, not an unchecked executable or a source build substituted for the release.

For a light smoke run, add `--part P0/near-miss.json --repeats 1`. The default covers all thirteen P0–P5 cases, with one warmup and three measured pairs per case. Every case has a committed scalar edit, including P1/P2/P4's existing scoreboard edits. `cases/` pins each original recipe's SHA-256 and carries the authored reference commands. A changed target requires reviewing its script and independent definition before updating the hash. P5 expands the four countersink cutters into explicitly named copies at the same positions. It does not omit the pattern's geometry.

Use a fresh output directory for each run. `head-to-head.json` checkpoints finished pairs, `report.md` is the table, and `cases/` retains commands, standard streams, phase reports, exported meshes, STEP and saved reference designs. A timeout is incomplete, never a typed geometric refusal. `--timeout` guards each external command; it never changes a geometric tolerance or work allowance. Use a dedicated benchmark host for the full run rather than a crowded development machine.

## What the columns mean

Catallax construction and edit each evaluate the complete recipe. Every phase has a separate fixed work allowance. The example measures the public feature builder, B-rep display mesh, requested volume/area/centroid and STEP export directly. Names are resolved on the two built results. Each phase retains its deterministic work and operation log; construction/edit count the executed recipe operations. There is no incremental edit cache claimed for Catallax. A refused construction does not prevent the independent edited recipe from being attempted. No successful prefix stands in for a whole target.

SolidCraft uses the pinned release's `run --results --timings --profile`. Construction saves a native design. A native `parameters.change` edits that design and reports `recomputed` and `restored`. Display OBJ export, requested measurement, name listing and STEP export run as separate commands against the saved original or edited design. Internal command timers exclude reopening and saving; process times including those costs are retained separately. Build and edit can already perform tessellation. Their `embedded_display_ms` records that cost without subtracting it from the time the command actually takes. The display column is an additional OBJ export, including serialization, not a claim of identical mesh algorithms or a display-only speedup. The release exposes no deterministic work counter, so its work is null with a reason.

Failures retain stderr and the command's whole-process refusal upper bound when the release does not emit per-command timings. They do not invent a refusal-only timer or call a timeout a refusal. Failed display or measurement does not erase a completed build. Warmup rows remain in the data but never enter the medians. Order alternates on one pinned core, with `RAYON_NUM_THREADS=1`, and each pair records load averages before and after. Correctness checks run after that timing window. Shared-machine results are not quiet-core latency claims. Timings from failed or incomplete phases are not replaced by zero.

Successful-phase medians require every measured sample to answer that phase. A fast display or measurement refusal appears in the separate refusal column, never as a successful mesh or query time. Raw rows keep each phase's outcome and refusal latency even when no successful median exists.

## Correctness and coverage

`check.py` constructs the targets independently in OpenCASCADE from their dimensions. It does not import either tool's result as its reference. For P1, two round slot ends plus the straight bridge give the same capsule as the overlapping cutters in the recipe. P4's planar inward shell is constructed by subtracting its exact rectangular cavity before applying the cylindrical corner trim. P5 uses the same four through holes, countersinks and junction round. The independent construction must pass topology validation before its properties are used.

Requested volume, area and all centroid components are compared with that independent construction. Catallax's exact rational interval endpoints are retained, rather than reduced to their midpoints. The independent calculation uses doubles; its diagnostic allowance is `max(1e-6, abs(value)*1e-8)`. This allowance belongs only to this external check. It is not a modelling tolerance, a certificate, a retry rule or permission to accept a kernel answer. Unavailable bindings or properties remain unmeasured. A passing numerical comparison is not exact-answer proof.

Each exported STEP is opened without repair in OpenCASCADE. The checker records root transfer, topology validity and the producing tool's expected unique face count. Independently constructed and exported models can split faces differently, so the producer's own topology count is checked rather than requiring both tools to have the same partition. A fixed 216-point grid compares imported and independent solid membership. Boundary classifications remain unresolved and are counted. Passing samples do not prove global set equality. Name checks require every original result face and edge name to resolve after the edit; they record missing and unnamed entities. String-name survival is not a proof of semantic face correspondence.

The existing `parts` command still runs independent B-rep, field and sampled comparison. This command compares B-rep target construction and its consumers with the released reference tool. It neither replaces field certification nor pretends to exercise SolidDraft's UI. No speed ratio is valid for a target that did not complete and pass the relevant correctness checks in both tools.

## Prior art

[OpenCASCADE's property API](https://dev.opencascade.org/doc/refman/html/class_b_rep_g_prop.html) distinguishes exact supporting geometry from triangulation and warns that volume integration does not itself validate a closed solid. The checker therefore validates the independent body before using its values and checks STEP separately. [SMLib's measurement documentation](https://smlibpy.classcad.ch/docs/analysis) separates integration accuracy from the accuracy of the constructed shape. The report keeps numeric checks distinct from exact-answer claims.

[truck](https://github.com/ricosjp/truck) and [monstertruck](https://github.com/ens1/monstertruck) separate modelling, meshing and exchange. The benchmark keeps those costs separate rather than calling a kernel build a displayed part. [Fidget](https://github.com/mkeeter/fidget) and [libfive](https://libfive.com) measure implicit evaluation and meshing, not the same B-rep edit and STEP workload. Their rendering numbers are not imported as CAD speed comparisons. Plantinga and Vegter's [isotopic approximation paper](https://www.math.rug.nl/~gert/documents/2004/pv-sgp-2004.pdf) requires interval bounds for its guarantee; a finite membership grid here is explicitly sampled evidence. No implementation code was taken from these projects.
