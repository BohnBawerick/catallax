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

## Five-tool comparison

`head-to-head` also runs CadQuery 2.6.1 on OpenCASCADE 7.8.1, smlibpy 0.4.0 on its bundled SMLib, and Pali Ijo 0.1.2. `tools.json` records release URLs, licenses and SHA-256 values. CadQuery has released headless Python wheels and the solid operations needed here. The public MIT [smlibpy](https://github.com/awv-informatik/smlibpy) wheel bundles SMLib under Apache-2.0; no closed tool is substituted. [Pali Ijo](https://github.com/braydennilsen/PaliCAD/tree/main/packages/pali_ijo) has no GitHub release or published pub.dev package. Its official source snapshot is pinned by commit and archive checksum, then compiled with the pinned Dart SDK without modifying the kernel.

Prepare outside the checkout, with Python 3.11, `uv`, Linux x86-64 and a C library compatible with the pinned wheels:

```sh
TOOLS=/tmp/catallax-comparison-tools
python3 -B .claude/skills/verify-catallax/scripts/head-to-head/prepare.py --out "$TOOLS"
```

Preparation checks the direct artifact hashes, installs the hash-locked Python dependencies, verifies Pali's upstream dependency lock, and writes `runtime.json` with the native runner checksum and installed Python/checker file hashes. The benchmark verifies that receipt before executing any additional tools. Do not edit an installation to change a result.

On the designated quiet benchmark host, after building the opt-in example and fetching the unchanged SolidCraft pin:

```sh
nice -n 10 python3 -B "$DRIVER" head-to-head parts \
  --binary target/release/examples/head_to_head \
  --reference /tmp/h2h/reference/solidcraft-cli \
  --tool-runtime "$TOOLS/runtime.json" \
  --machine-class 'quiet work PC' --core 0 --max-work 100000000 \
  --repeats 3 --out /tmp/h2h/five-tools --json
```

Select an allowed core rather than assuming core 0 is available. Tool order rotates and reverses on odd samples. Construction and edit are complete target-script rebuilds for the new tools. Pali uses the same authored operation graph through its public Dart API, not a modified CLI or a mesh substitute. `targets.py` contains the thirteen authored recipes beside the existing SolidCraft scripts; `python_tools.py` and `pali_parts.dart` implement those recipes for each tool. The P1 capsule and P5 expanded pattern have the exact same sets as the original recipe. P4 invokes each tool's native shell operation. Circular edge picks exclude artificial seam lines and include all ring pieces meeting the pick.

Successful display is a separate STL mesh export for the new tools. Meshing algorithms and densities differ; display is not used for build speed ratios. Every phase has its own outcome. Unsupported whole parts are `not-expressible`; timeouts are incomplete; other build errors remain failed rebuilds. No exported successful prefix stands in for a target. New tools' deterministic work is unreported, not zero.

The additional tools have no persistent face/edge name contract across these complete rebuilds. Pali's ancestry IDs are process-local, not stable document-independent rebuild names. Their name checks remain not run rather than comparing traversal indices or inventing names.

Publication imports only a completed thirteen-target, five-tool run labelled with the designated benchmark machine class:

```sh
python3 -B site/render/comparison.py --import-run /tmp/h2h/five-tools/head-to-head.json
python3 -B site/render/comparison.py
```

The first command projects public results into `site/comparison-results.json`, excluding installation paths and subprocess streams. The second regenerates `site/comparison.html` using only committed result and pin files. The page records source, tool, recipe and script hashes, machine class, CPU model, load range and measurement scope. Never publish local smoke timings as quiet-host evidence.


## What the columns mean

Catallax construction and edit each evaluate the complete recipe. Every phase has a separate fixed work allowance. The example measures the public feature builder, B-rep display mesh, requested volume/area/centroid and STEP export directly. Names are resolved on the two built results. Each phase retains its deterministic work and operation log; construction/edit count the executed recipe operations. There is no incremental edit cache claimed for Catallax. A refused construction does not prevent the independent edited recipe from being attempted. No successful prefix stands in for a whole target.

SolidCraft uses the pinned release's `run --results --timings --profile`. Construction saves a native design. A native `parameters.change` edits that design and reports `recomputed` and `restored`. Display OBJ export, requested measurement, name listing and STEP export run as separate commands against the saved original or edited design. Internal command timers exclude reopening and saving; process times including those costs are retained separately. Build and edit can already perform tessellation. Their `embedded_display_ms` records that cost without subtracting it from the time the command actually takes. The display column is an additional OBJ export, including serialization, not a claim of identical mesh algorithms or a display-only speedup. The release exposes no deterministic work counter, so its work is null with a reason.

Failures retain stderr and the command's whole-process refusal upper bound when the release does not emit per-command timings. They do not invent a refusal-only timer or call a timeout a refusal. Failed display or measurement does not erase a completed build. Warmup rows remain in the data but never enter the medians. Order alternates on one pinned core, with `RAYON_NUM_THREADS=1`, and each pair records load averages before and after. Correctness checks run after that timing window. Shared-machine results are not quiet-core latency claims. Timings from failed or incomplete phases are not replaced by zero.

Successful-phase medians require every measured sample to answer that phase. A fast display or measurement refusal appears in the separate refusal column, never as a successful mesh or query time. Raw rows keep each phase's outcome and refusal latency even when no successful median exists.

## Correctness and coverage

`check.py` uses independent scalar definitions in `analytic.py`, not an OpenCASCADE-built reference. Eight target families have elementary property formulas. Circular blends use smooth one-dimensional section quadrature, checked at two resolutions. Four intersecting curved P0 families have analytic membership but no implemented property formula: block-sphere, crossing-pipes, general-cone and general-torus. P5's junction blend reaches the countersink openings at its flange contact, so an untrimmed full torus is not an independent definition of the complete fillet. P5 properties and membership remain unmeasured rather than checking an approximate blend or the unfilleted prefix. Historical reports used an independently built OpenCASCADE reference; do not relabel their results as analytic checks.

Requested volume, area and all centroid components are compared with those analytic values. Catallax's exact rational interval endpoints are retained rather than reduced to midpoints. The diagnostic allowance is unchanged at `max(1e-6, abs(value)*1e-8)`. It belongs only to the external check, not a modeling tolerance, certificate, retry rule or permission to accept a kernel answer. Unavailable properties remain unmeasured.

Each exported STEP is opened without repair in OpenCASCADE 8.0.1.1.0. OpenCASCADE is both the CadQuery competitor's kernel and the STEP checker, so STEP validity is not independent of that implementation family. Property and expected membership values do not use OpenCASCADE. The checker records root transfer, topology validity and the producing tool's unique face count. A fixed 216-point grid uses independent analytic bounds and solid predicates. Boundary classifications remain unresolved and counted; passing samples do not prove set equality. Name checks require every original result face and edge string name to resolve after the edit, not semantic face correspondence.

The table separates STEP files opened, checks completed and full passes. A valid file can still change the producer's expected face count, and an unavailable checker has a zero checked denominator.

The existing `parts` command still runs independent B-rep, field and sampled comparison. This command compares B-rep target construction and its consumers with the released reference tool. It neither replaces field certification nor pretends to exercise SolidDraft's UI. No speed ratio is valid for a target that did not complete and pass the relevant correctness checks in both tools.

## Prior art

[OpenCASCADE's property API](https://dev.opencascade.org/doc/refman/html/class_b_rep_g_prop.html) distinguishes exact supporting geometry from triangulation and warns that volume integration does not itself validate a closed solid. The checker therefore validates the independent body before using its values and checks STEP separately. [SMLib's measurement documentation](https://smlibpy.classcad.ch/docs/analysis) separates integration accuracy from the accuracy of the constructed shape. The report keeps numeric checks distinct from exact-answer claims.

[truck](https://github.com/ricosjp/truck) and [monstertruck](https://github.com/ens1/monstertruck) separate modelling, meshing and exchange. The benchmark keeps those costs separate rather than calling a kernel build a displayed part. [Fidget](https://github.com/mkeeter/fidget) and [libfive](https://libfive.com) measure implicit evaluation and meshing, not the same B-rep edit and STEP workload. Their rendering numbers are not imported as CAD speed comparisons. Plantinga and Vegter's [isotopic approximation paper](https://www.math.rug.nl/~gert/documents/2004/pv-sgp-2004.pdf) requires interval bounds for its guarantee; a finite membership grid here is explicitly sampled evidence. No implementation code was taken from these projects.

The multi-tool pass also checked [smlibpy's kernel and accuracy notes](https://smlibpy.classcad.ch/docs/kernel) before retaining its released chamfer discrepancies and curved-Boolean refusals. Its measurement API distinguishes numeric integration error from fitted-geometry gaps. The benchmark therefore keeps the existing diagnostic allowance, records discrepancies, and never changes dimensions or repairs an exported answer to make a tool pass. The analytic predicates follow the functional-representation approach described by Fidget and libfive, but no meshing rate from those projects is presented as a B-rep build result. No upstream implementation code, example scripts or assets were copied.
