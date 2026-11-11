# R8 public STEP census after the METRE fix

Date: 11 November 2026. Sample: the same 3,000 ABC model ids, seed `20261111`.

The import entry accepts **189 models**, of which **186 have complete source
coverage**, or **6.20% of the sample**. Three accepted models still have
unverified representation items and are excluded from that complete-coverage
admission figure. The entry returns 2,797 typed refusals; 14 processes reach
the separate observation timeout. The earlier `METRE` unit gate is gone:
`unknown_unit` falls from 2,992 models to zero on the identical sample.

Whole-file section checks establish 2,883 complete source traversals. Of all
3,000 models, 2,002, or 66.73%, have only the five analytic face-carrier
families. Spline faces are observed in 868 models. Fillet and chamfer feature
counts remain unknown for all 3,000; R8's feature-provenance gap remains open.

The [per-model table](R8-public-step-census.csv) records counts, section scope,
work, typed outcomes and public model ids. The [run manifest](R8-public-step-census.json)
pins the source, tools, executable and archive. The
[earlier preflight](R8-public-step-preflight.md) is a separate acquisition attempt.

## Import outcomes

The production source is main snapshot
`d6ea026f00eb269467c282c9cae657cfa85e3beb`, containing the
[merged unit fix](https://github.com/BohnBawerick/catallax/pull/418).
The manifest's source/tooling snapshot is
`8f029060e6da8a94d2df6598ec3e7a177d2592b0`. Its production crates match that
main snapshot; the additional code is the research observer and collector.

| Outcome | Models | Share of all 3,000 |
|---|---:|---:|
| Entry accepted, complete source coverage | 186 | 6.20% |
| Entry accepted, source coverage unverified | 3 | 0.10% |
| Returned typed refusal | 2,797 | 93.23% |
| External observation timeout | 14 | 0.47% |

The CSV's `admitted` outcome combines the first two rows: the actual import
entry returned bodies. `coverage_complete` separates the three whose whole
source representation remains unverified. Every evaluation uses a fresh
model; partial builder mutation is never admission.

These are the first returned refusal causes, sorted by count. Removing one
cause would not prove that a subsequent import succeeds.

| First refusal cause | Models |
|---|---:|
| `malformed` during syntax/geometry recognition | 993 |
| `self_intersecting` | 587 |
| `resource_limit` | 255 |
| Unsupported spline edge solid certification | 240 |
| `not_solid`, exact incidence checks | 178 |
| Unsupported `SURFACE_OF_LINEAR_EXTRUSION` | 141 |
| Unsupported spline face solid certification | 103 |
| Unsupported `SHELL_BASED_SURFACE_MODEL` | 91 |
| `self_intersection_not_verified` | 72 |
| Unsupported `BREP_WITH_VOIDS` | 33 |
| Unsupported cone variants | 30 |
| Unsupported `SURFACE_OF_REVOLUTION` | 29 |
| Unsupported `OFFSET_SURFACE` | 22 |
| Unsupported torus variants | 21 |
| Unsupported vertex loop beside other bounds | 2 |

The manifest/CSV name the offset bucket `unsupported:other_step_type`;
all 22 underlying diagnostics identify `OFFSET_SURFACE`. Cone refusals
comprise 13 non-rational-axis-length cases, nine incomplete rim bounds,
six off-axis vertex loops and two incomplete circles about the axis.

Of the 993 `malformed` returns, 982 report the placement axis/reference
perpendicularity check, nine report non-increasing distinct knots, one
reports no `MANIFOLD_SOLID_BREP`, and one is the syntax failure at public
model id `00003087`. These are importer diagnostics, not an independent
judgement that each source file is invalid. The 178 `not_solid` returns
comprise 150 edge vertices not exactly on their curves, 16 edges not
exactly on their faces and 12 vertex-loop points not exactly on their faces.

The 255 resource-limit returns, 72 unverified-intersection returns and 14 external
timeouts are incomplete observations or proofs, not impossible-geometry
claims. Eleven timeouts retain complete source counts; three have no source
checkpoint. Work allowances and tolerances were not raised.

## Whole-model and face coverage

The sample contains 3,000 distinct model ids and byte hashes, with no exact
byte duplicates. It samples one 10,000-model chunk, not the million-model
population. Geometric duplicates are not detected. A model can contain
several solids: the observed graphs contain one `MANIFOLD_SOLID_BREP` in
1,948 models, several in 1,012, none in 36, and no completed graph in four.
The sample unit is models, not individual solids or manufactured parts.

| Source observation | Models | Share of all 3,000 |
|---|---:|---:|
| Every face has an analytic carrier, complete traversal | 2,002 | 66.73% |
| Carrier outside the five analytic families, complete traversal | 881 | 29.37% |
| Incomplete source coverage | 117 | 3.90% |
| Spline face observed, including partial traversals | 868 | 28.93% |
| Spline absence established | 2,002 | 66.73% |
| Spline need unknown | 130 | 4.33% |

The first three rows partition coverage; the last three partition spline
need. Among the 2,883 complete traversals, 2,002, or **69.44%**, have only
plane, cylinder, cone, sphere and torus carriers. Carrier-family reach does
not establish admissible trims, supported boundary curves or solid topology.

The full-file scan records exactly one supported DATA section for 2,997
models; three have no source checkpoint. No additional or unsupported DATA
section is observed among those 2,997. Source syntax parses in 2,996 models,
fails in one, and is unobserved in three. The 117 incomplete traversals are
113 models with unknown representation items, one syntax failure and three
missing checkpoints. Completed observations contain no missing references,
unreferenced faces, unresolved face layouts or assembly/mapping entities.

Earlier whole-file percentages were withdrawn because the retained first-run
records did not establish section scope. This rerun reads the files again
with the corrected guard. The restored percentages above come from fresh
section-qualified observations, not from the earlier reassessment.

The observer counts **704,829 referenced face entities**, deduplicating
shared entity references within each model. Counts include observed faces
of incompletely covered models and do not estimate their unseen geometry.

| Observed face carrier | Faces | Share of observed faces | Models containing it |
|---|---:|---:|---:|
| Plane | 354,968 | 50.36% | 2,956 |
| Cylinder | 166,637 | 23.64% | 2,609 |
| B-spline | 54,799 | 7.77% | 339 |
| Rational B-spline | 52,530 | 7.45% | 749 |
| Torus | 24,770 | 3.51% | 925 |
| Cone | 22,329 | 3.17% | 943 |
| Other | 21,324 | 3.03% | 341 |
| Sphere | 7,472 | 1.06% | 489 |

`other` is a direct carrier outside the listed families. Known spline
presence can survive incomplete traversal; unobserved geometry or an unknown
carrier cannot establish spline absence. Fillet and chamfer fields use
`unknown_no_feature_history_in_step`: surface curvature alone does not
establish feature provenance. Those unknown counts are not zeros.

M24's measured priorities now include placement recognition, exact incidence,
unsupported representations/carriers and the returned intersection outcomes.
M21 receives a section-qualified analytic-carrier candidate set of 2,002
models and 868 models with observed spline faces. This measures the B-rep
import entry; field construction and independent two-kernel checks remain
separate work. No independent source-validity oracle was used.

## Method and validation

The collector sorts public ABC model ids and selects
`random.Random(20261111).sample(ids, 3000)`, then sorts results by id.
The downloaded archive matches the earlier run in size, MD5 and SHA-256;
all selected member identities match. Where both runs completed syntax
parsing, every observed face-kind count matches.

The observer scans the entire file for DATA declarations outside comments
and strings. Complete coverage requires exactly one supported literal
`DATA;`, a complete lexical scan and no hidden parser marker in a string.
Additional or unsupported sections make completeness and absence unknown.
The production Part 21 parser itself is unchanged and reads its first literal
DATA section. The observer traverses references from B-rep, shell and shape
representation roots; unknown items and transformations prevent a complete
coverage claim.

Import calls the public recipe entry `import_step_bodies` with the actual
input hash and a fresh `Model`. Numeric decider defaults and the existing
20,000,000-unit application build allowance remain unchanged. Exhausted
work discards the returned result as `resource_limit`. No repairs are applied.
A separate 60-second process observation cap is not a speed gate or a
geometric refusal. Source counts are flushed before import.

The main run uses 12 processes. The smallest 24 selected ids are repeated
with 1, 2 and 4 workers. Public results and complete private stdout/stderr
match the 12-worker baseline in all three repeats. This is subset determinism
evidence, not a full-sample worker comparison.

Six Rust controls and six Python tests pass, including the two-DATA-section
regression and handling legacy records without section metadata. Strict
example Clippy passes on the same production source. An independent audit
checks the sample, CSV/manifest arithmetic, section qualifications, actual
refusal codes, admitted body counts, work limits and source/executable
hashes. The saved 477 replay inputs/answers are unmodified and were not rerun
by this observation-only change.

The complete archive and extracted models were deleted after evaluation;
absence was checked directly. Build output was removed with `cargo clean`.
The retry started with 63 GB free and had 61 GB before extraction. Raw
evidence retained privately.

## Source and licence scope

The source is [ABC chunk 0000](https://archive.nyu.edu/handle/2451/44309),
STEP v00, released in 2019. The [dataset authors](https://deep-geometry.github.io/abc-dataset/)
describe 10,000-model chunks. Attribution: Koch et al., *ABC: A Big CAD Model
Dataset for Geometric Deep Learning*, CVPR 2019.

NYU labels the dataset MIT. The
[author repository licence](https://github.com/deep-geometry/abc-dataset/blob/392789113158b874bfa630098f5bbab3ed38309e/LICENSE)
is MIT, copyright 2019 Deep Geometry Processing. Use here is aggregate
analysis under that dataset grant. Per-model source exceptions were **not
individually verified**; this is not a per-model licence clearance.

The authors retain creators' CAD copyrights and refer to
[Onshape terms section 1.g.ii](https://www.onshape.com/en/legal/terms-of-use),
which preserve some documents' greater reserved rights in LICENSE tabs.
Their effective date, 15 July 2020, follows the dataset release. The earlier
metadata audit found 10,000 public documents predating 7 August 2018, with
533 dated 2014 and 9,467 dated 2015. Metadata lacked licence/tab listings;
an unauthenticated tab-list request returned HTTP 403. Those individual
exceptions remain unresolved.

Only counts, typed outcomes and public model ids are published. No model
file, geometry, render, excerpt or archive-member path is included. Downloaded
metadata from the earlier audit was also deleted. Raw evidence retained privately.

## Reproduce

Requirements: Linux, Python 3.11 or later, `bsdtar`, `curl` and Rust. Use the
manifest's source/tooling revision. Set `CENSUS_SCRATCH` to a fresh private
scratch directory outside the checkout on a volume with at least 30 GB free.

```sh
CENSUS_ROOT="$(git rev-parse --show-toplevel)"
: "${CENSUS_SCRATCH:?Set a private scratch directory outside the checkout}"
mkdir -p "$CENSUS_SCRATCH"
python3 -c 'import shutil,sys; assert shutil.disk_usage(sys.argv[1]).free >= 30_000_000_000' "$CENSUS_SCRATCH"
export CARGO_TARGET_DIR="$CENSUS_SCRATCH/target"
cargo test --manifest-path "$CENSUS_ROOT/Cargo.toml" -p catallax-import --example r8_census --release
cargo build --manifest-path "$CENSUS_ROOT/Cargo.toml" -p catallax-import --example r8_census --release
python3 -B -m unittest discover -s "$CENSUS_ROOT/tools/research" -p test_r8_step_census.py -v
curl --fail --location --http1.1 --connect-timeout 30 --max-time 1800 \
  --speed-time 60 --speed-limit 1024 \
  https://archive.nyu.edu/bitstream/2451/44309/3/abc_0000_step_v00.7z \
  -o "$CENSUS_SCRATCH/abc_0000_step_v00.7z"
python3 -B "$CENSUS_ROOT/tools/research/r8_step_census.py" \
  --scratch "$CENSUS_SCRATCH" \
  --binary "$CARGO_TARGET_DIR/release/examples/r8_census" \
  --count 3000 --seed 20261111 --workers 12 --timeout 60 --validate-workers 24
cargo clean --manifest-path "$CENSUS_ROOT/Cargo.toml"
```

The collector checks free space, the expected 1,594,129,754-byte size and
published MD5 `695388be7a278c7798c8c8ae239772ac` before extraction. Actual
archive SHA-256 is in the manifest. On transfer failure, discard the partial
file and retry from zero; the server refused range resumes during this run.
Only selected members are extracted. A `finally` block deletes the archive
and extracted models on completion or a collection exception; a forced
process kill can bypass it, so inspect cleanup and remove remaining copies.
The public CSV uses LF line endings.

The [observer](../../crates/catallax-import/examples/r8_census.rs),
[collector](../../tools/research/r8_step_census.py) and
[tests](../../tools/research/test_r8_step_census.py) are checked in. The
collector's `--reassess` mode conservatively renders retained legacy records
without section metadata; it does not replace a fresh whole-file scan.
