# R8 spline certification after the blanket gate

Date: 12 November 2026. Scope: exactly the 532 `EDGE_CURVE` and 196
`ADVANCED_FACE` first refusals selected from the placement-corrected census
in [the source comparison](https://github.com/BohnBawerick/catallax/pull/430),
revision `0951013dab04f79ad0ca4b3779c3112e8be711f5`.

**New real-model admissions: zero.** All 728 models return typed refusals
at their next proof steps. The import route now admits the exact straight
spline edge and complete embedded patch controls. Those controls establish
capability in the kernel, without changing the real-model admission count.

The [per-model table](R8-spline-certification.csv) records the next step,
work and source coverage. The [manifest](R8-spline-certification.json) pins
executables, production source files, observer, archive and worker checks.
Raw evidence retained privately. No model files or geometry are published.

## Where the original sample stopped

Baseline source: `cdc0a3884c59b2bea538c230e37639bd74cdb9a2`.
`random.Random(20261111)` selects 12 sorted public ids from the edge bucket,
then 12 from the face bucket. All 24 reproduce the original refusal:
12 `EDGE_CURVE (spline solid certification)` and 12
`ADVANCED_FACE (spline solid certification)`.

Recognition had already read the exact source definitions. `build_solid`
then rejected any spline unconditionally, before incidence,
self-intersection or topology construction. No failed geometric certificate
explained those original refusals. Three synthetic production-entry tests
also failed at that gate before the change.

## Admission checks now attempted

- Spline curves need unique finite endpoint preimages, consistent parameter
  direction and a strictly monotone projection over every knot span. The
  projection proof checks the exact homogeneous derivative numerator.
  Collinear certified arcs can use an exact line carrier.
- Curve incidence uses exact homogeneous polynomial substitution on the
  spans crossed by the physical edge. Unused source spans do not invalidate
  a finite edge. Cone nappe membership still needs its own proof.
- Spline surfaces pass the existing `RationalPatch` embedding certificate.
  `RationalPatchFace` checks complete physical boundary rows, endpoint
  identity and material-relative orientation. It supplies genuine incidence
  facts to the existing M10 self-intersection machinery.
- Finite spline boundaries supply intersection breakpoints through their
  endpoints. They no longer ask an unused whole-carrier equality question
  in that path. Ring handling and the existing analytic query order remain.
- Broad-phase spline face bounds include the full positive-weight surface
  control hull, including interior controls. Boundary vertices alone do not
  bound a bulging patch.

Every returned body still passes self-intersection and topology checks.
A support certificate alone does not admit a solid. The
[interface fragment](../interface-changes/2026-11-12-step-spline-certification.md)
states the construction scope. Spline admission adds no fit, repair or
nonzero local tolerance. Existing tolerances and work allowances remain.

## Next first refusals

| First refused step | Models |
|---|---:|
| Analytic edge-on-face incidence | 133 |
| Spline finite endpoint incidence | 126 |
| Analytic vertex-on-curve incidence | 120 |
| Spline whole-patch embedding | 108 |
| Cone bound is not one complete rim circle | 62 |
| Closed edge has no supported closed carrier | 58 |
| Spline edge-on-face incidence | 33 |
| Cone axis length is not rational | 28 |
| Import operation reaches its 20M work cap | 26 |
| Spline curve embedding projection unproved | 12 |
| Ellipse carrier unsupported | 12 |
| Curve/surface pair unsupported | 4 |
| Spline complete boundary certificate unproved | 3 |
| Cone vertex loop is off its axis | 2 |
| Cone rim is not a complete circle about its axis | 1 |
| **Total** | **728** |

The 412 incidence refusals identify exact source mismatches under this
import contract. They are not independent judgements that the original CAD
parts are unusable. Curved repair would need its own measured, recorded
bound under M24. No distance threshold was used to turn these refusals
into admissions.

The 108 surface refusals carry `shapes.rational-patch-unproven`; the 12
curve projection refusals carry `import.spline-embedding-unproven`; the
three boundary refusals carry `numbers.geometric-gap-unproven`. These
mean the available checker did not prove its requirement. They do not
establish that a different checker cannot certify the source.

The initial 60-second observation cap left six processes incomplete.
Longer observation returned four work-limit refusals and two unproved
patch-support refusals. Final qualification uses a 300-second external
observation cap, with the same 20,000,000-unit import allowance. All 728
return typed results. Wall time is observation information, not a geometry
or speed gate. This shared-machine run makes no latency claim.

## Validation

- All 137 import tests pass, including the interior-control bound and finite
  trim regressions. Strict import Clippy and workspace formatting pass.
- Seven collector tests pass. Explicit model-id filtering preserves the
  seeded sample order and rejects unknown or duplicate ids. The timeout
  checkpoint test injects a timeout after its checkpoint, so interpreter
  startup scheduling cannot decide its result.
- The lowest 24 selected ids repeat with 1, 2 and 4 workers. Ordered records,
  work and raw observer stdout/stderr match the 12-worker qualification.
  This is subset worker evidence, not a full 728-model worker comparison.
- The three synthetic controls agree at 1, 2 and 4 kernel threads and test
  workers: straight spline edge admission uses 7,355 work units; complete
  patch admission uses 33,676; off-face spline refusal uses 2,428.

The 477 saved replay answers and STEP regressions remain required PR CI
checks. The import-only measurements here do not substitute for those gates.

## Reproduce the subset

Use the archive acquisition and licence scope in the
[original census method](R8-public-step-census.md). Obtain the pinned chunk
in a fresh private scratch directory outside the checkout, with at least
30 GB free. `CENSUS_ROOT`, `CENSUS_SCRATCH` and `CENSUS_BINARY` below are
caller-supplied paths. Build the `r8_census` release example at the reviewed
source revision, then derive the public id filter from this report:

```sh
python3 - "$CENSUS_ROOT/docs/milestones/R8-spline-certification.csv" \
  "$CENSUS_SCRATCH/spline-model-ids.txt" <<'PY'
import csv, pathlib, sys
rows = list(csv.DictReader(open(sys.argv[1])))
pathlib.Path(sys.argv[2]).write_text(''.join(r['model_id'] + '\n' for r in rows))
PY
python3 -B "$CENSUS_ROOT/tools/research/r8_step_census.py" \
  --scratch "$CENSUS_SCRATCH" --binary "$CENSUS_BINARY" \
  --count 3000 --seed 20261111 \
  --model-ids "$CENSUS_SCRATCH/spline-model-ids.txt" \
  --workers 12 --timeout 300 --validate-workers 24
```

The filter is applied after the original seeded 3,000-model selection, so
an id outside that sample refuses. Only the 728 requested members are
extracted. The collector retains raw observations privately and removes
the archive and extracted model files in its cleanup block.
