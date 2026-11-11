# R8 planar self-intersection audit

Date: 12 November 2026. Cohort: the 745 planar-polygon first refusals identified
by the [saved-pair triage](https://github.com/BohnBawerick/catallax/pull/430).

**No new models import.** All 745 now return `not_solid` because a polygon
boundary is not exactly incident to its declared plane. The old planar pass
applied point-location without checking that precondition and reported
`self_intersecting` instead. This corrects the refusal; it does not repair or
admit the source geometry.

The [per-model table](R8-planar-self-intersection.csv) contains every model's
classification and new first refusal. Its `sampled` column identifies the
seeded 30-model sample. The [manifest](R8-planar-self-intersection.json) records
the sample ids, source and executable hashes, allowances and worker comparisons.

## What the sample establishes

Select `random.Random(20261111).sample(sorted(planar_ids), 30)` from the exact
745-model triage cohort. Fresh baseline imports reproduce all 30 saved refusals,
including the same solid and face-pair references.

Every sampled pair shares exactly one declared edge and two declared vertices.
Every pair also has at least one boundary vertex for which
`normal dot (vertex - plane_origin) != 0`, in exact rational arithmetic.
The classification of every sample is therefore `nonincident_boundary`.

| Sample classification | Models |
|---|---:|
| Verified true interior crossing | 0 |
| Verified legitimate touching contact | 0 |
| Inconsistent plane/boundary data; geometric relation unverified | 30 |

These are not 30 proven touching contacts. Nor do their existing diagnostics
prove 30 true crossings. The file declares a shared boundary whose endpoints
do not lie on both declared carriers. Inferring the intended CAD geometry
would require a separate, bounded source-repair operation. No tolerance,
projection, carrier replacement or topology exemption is introduced here.

The source audit was extended to all 745 saved pairs. All have nonincident
boundaries. Of these, 741 share declared edges; the remaining four share a
declared vertex. All 745 audits complete. An independent rational parser agrees
with the production decoder on the seeded sample after exact METRE conversion.
Raw evidence retained privately.

## Cause and change

The planar pass intersects two carrier planes, splits their intersection line
at boundary events, then classifies interval midpoints. Boundary membership
uses three-dimensional segments, while interior membership uses a planar
projection. A boundary off its carrier violates the assumptions of that test.
Projection can then report a midpoint as interior even beside a declared shared
edge.

Before constructing a planar polygon, check both endpoints of every straight
boundary edge against the declared plane. Any nonzero residual returns the
existing typed `NotSolid` refusal with the vertex, edge and face named. Curved
boundary dispatch and the exact crossing predicate are unchanged.

The reduced regression starts with a clean unit box and changes one component
of its top-plane normal to exactly `1e-80`. The baseline returns
`SelfIntersecting`; the candidate returns `NotSolid`. The nonzero discrepancy
is never rounded away.

Exact shared edges, shared vertices, coplanar adjacency and partial boundary
contacts remain allowed. True interior crossings still refuse. A separate
regression checks two faces that share an edge and also cross elsewhere, so
declared adjacency cannot hide a crossing.

## Import rerun

Each model uses the public `import_step_bodies` entry, a fresh model, the actual
source hash and the unchanged 20,000,000-unit allowance. No check is bypassed.

| Outcome | Saved census cohort | Candidate rerun |
|---|---:|---:|
| Admitted models | 0 | 0 |
| First refusal `self_intersecting` | 745 | 0 |
| First refusal `not_solid` | 0 | 745 |
| Resource limits | 0 | 0 |
| External observation timeouts | 0 | 0 |

The new first refusal is the planar endpoint/plane incidence check in every
model. Later solids, incidence checks and geometric relations remain
unobserved after that refusal. This is zero admission gain, not a completed
source-repair capability.

The main rerun uses 12 processes. The same seeded 30 models repeat with 1, 2
and 4 workers and `CATALLAX_THREADS` set to 1, 2 and 4 respectively. Public
records, work counts and private stdout/stderr bytes match the main run in all
three repeats. This proves sample determinism, not a full-cohort thread sweep.
No STEP output is claimed for refused imports.

## Controls and evidence

The self-intersection controls pass all 14 tests. They include oblique touching
cases, a true crossing, a shared edge with a crossing elsewhere and the tiny
nonincidence regression. The collector controls check cohort selection, unsafe
archive members and exclusion of private paths and diagnostics from public rows.

The 477 saved input/answer records and STEP regression fixtures are unmodified.
Fresh replay and STEP regression results are provided by the PR check workflow.
Source data, archive members and raw geometry are not published. The dataset
attribution and unresolved per-model licence scope remain those of the
[original census](R8-public-step-census.md#source-and-licence-scope).

## Reproduce

Use the retained census selection, raw refusal diagnostics and pair-details
mapping. Extract only the 745 selected members from the checksum-pinned ABC
chunk into private scratch, following the
[census acquisition checks](R8-public-step-census.md#reproduce). Check free disk
space first. Set the following paths at runtime; evidence stays outside the
checkout.

```sh
: "${PLANAR_SCRATCH:?Set private scratch outside the checkout}"
: "${CENSUS_SELECTION:?Set the retained selection JSON path}"
: "${CENSUS_REFUSALS:?Set the retained raw diagnostics directory}"
: "${CENSUS_PAIRS:?Set the retained pair-details JSON path}"
PLANAR_ROOT="$(git rev-parse --show-toplevel)"
export CARGO_TARGET_DIR="$PLANAR_SCRATCH/target"
cargo build --manifest-path "$PLANAR_ROOT/Cargo.toml" --release \
  -p catallax-import --example r8_census --example r8_planar_contact
python3 -B "$PLANAR_ROOT/tools/research/planar_self_intersection.py" \
  --models "$PLANAR_SCRATCH/selected-step" --selection "$CENSUS_SELECTION" \
  --saved-refusals "$CENSUS_REFUSALS" --pair-details "$CENSUS_PAIRS" \
  --audit-binary "$CARGO_TARGET_DIR/release/examples/r8_planar_contact" \
  --import-binary "$CARGO_TARGET_DIR/release/examples/r8_census" \
  --out "$PLANAR_SCRATCH/results" --workers 12 --seed 20261111 --sample 30
cargo test --manifest-path "$PLANAR_ROOT/Cargo.toml" --release \
  -p catallax-import --test self_intersection --example r8_planar_contact
python3 -B -m unittest discover -s "$PLANAR_ROOT/tools/research" \
  -p test_planar_self_intersection.py -v
cargo clean --manifest-path "$PLANAR_ROOT/Cargo.toml"
```

Delete the downloaded archive and extracted models after evaluation. Audit
timeouts and process errors are recorded as unverified observations, never as
crossings or successful imports. The public CSV uses LF line endings.
