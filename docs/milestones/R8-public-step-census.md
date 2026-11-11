# R8 public STEP census

Date: 11 November 2026.

Whole-file carrier coverage is **unknown for all 3,000 sampled ABC models**.
A review found that the original observer treated the syntax parser's first
DATA section as the whole file. The retained records do not include section
metadata and the source files were deleted, so the old 66.73% analytic-only
sample share and 69.44% complete-traversal share are withdrawn. The corrected
observer verifies section scope before making completeness or absence claims.

The observed graph prefixes contain spline faces in 868 models. None imports successfully:
2,992 stop at `unknown_unit`, one is malformed and seven exceed the external
observation timeout. Fillet and chamfer provenance remains unknown for all
3,000 models. This run supplies the import census; R8's feature-count gap
remains open.

The [per-model table](R8-public-step-census.csv) contains counts, coverage
flags, typed outcomes and public model ids. The
[run manifest](R8-public-step-census.json) records aggregate counts and
reproduction identities. The [earlier preflight](R8-public-step-preflight.md)
describes a separate executor whose archive access failed.

## Measured coverage

The sample has 3,000 distinct model ids and 3,000 distinct byte hashes, with
no exact duplicates. An ABC model can contain several solids: 1,948 sampled
model prefixes contain one `MANIFOLD_SOLID_BREP`, 1,012 contain several, 36
contain none and four have no completed syntax census. These are prefix
entity counts, not verified whole-file solid counts. The denominator is
models, not individual solids or manufactured parts.

| Source observation | Models | Share of all 3,000 |
|---|---:|---:|
| Whole-file completeness established | 0 | 0% |
| Whole-file completeness unknown | 3,000 | 100% |
| Spline face observed in a parsed prefix | 868 | 28.93% |
| Spline absence established | 0 | 0% |
| Spline need unknown | 2,132 | 71.07% |

The first two rows partition file scope; the last three partition spline
need. The old observer marked 2,883 graph prefixes complete under its local
traversal checks, including 2,002 with only analytic carriers. Those counts
do not establish whole-file coverage. The manifest's corrected complete
denominator is zero and its whole-file analytic-only share is null.

Other prefix limitations remain: 113 models have unknown representation
items, one has a syntax failure and three timed out before the source checkpoint.
There were no observed missing references, unreferenced faces, unresolved
face layouts or assembly/mapping entities. Unknown representation items
remain an explicit coverage limit.

The observer counted **704,829 referenced face entities in parsed graph
prefixes**. These counts do not estimate unobserved sections or geometry.
Faces with complex rational B-spline representations
are counted separately. `other` means a direct face carrier outside the
listed families, not a successful recognition or construction.

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

## Import outcomes and milestone routing

| First outcome | Models | Observation |
|---|---:|---|
| `unknown_unit` | 2,992 | Every diagnostic reports `METRE` as the found unit |
| `observation_timeout` | 7 | External 60-second limit; four retain complete face counts |
| `malformed` | 1 | Invalid `EDGE_CURVE` parameters, model id `00003087` |
| Admitted | 0 | No complete successful import |
| Other typed refusal, including `resource_limit` | 0 | No such returned result |

M24's immediate measured obstacle is unit recognition. These first refusals
mask later geometry checks, so this run cannot rank seam, trim or incidence
failures behind them. The 2,002 analytic-only prefixes are candidates for a
new whole-file census, not an established M21-compatible part population.
Spline faces were observed in 868 prefixes. The next run must record
whole-file section scope before publishing an analytic-only or spline-absence share.
Fillet/chamfer feature provenance still needs an independently justified
source before it can guide blend counts.

## Validation and cleanup

The corrected observer passes six Rust tests and the collector passes six
Python tests. The two-section CLI regression was reproduced before the fix
and passes afterwards. Additional tests cover named/unsupported sections,
comment/string markers and legacy records without section evidence. Strict
Clippy for the example passed. During the original run, a 24-model subset, the smallest ids in the
seeded selection, was rerun with 1, 2 and 4 workers. Every public result and
full private stdout/stderr diagnostic matched the original 12-worker run.
This is subset determinism evidence, not a full-sample worker comparison.

An independent CSV audit recomputed the seeded selection, every face total,
outcome/refusal totals, feature unknowns and work limits. It checked returned
refusal codes against private stderr and verified source/executable hashes
against the manifest. The 477 saved replay files are unmodified; this
observation-only change did not rerun them. The section correction regenerates
the table and aggregates from retained observations. It preserves observed
face counts and import outcomes, while marking whole-file completeness and
unproven absence unknown. No model files were reread. The manifest preserves
the original measurement identities and separately pins this reassessment.

Archive and extraction cleanup succeeded, confirmed both by `cleanup.json`
and by checking that no downloaded archive or STEP file remained in the
scratch data. The initial free disk was 56 GB; 54 GB remained at an observed
extraction checkpoint. Build output was removed after validation. Raw
evidence retained privately.

## Source and licence scope

The source is [ABC chunk 0000](https://archive.nyu.edu/handle/2451/44309),
STEP v00, released in 2019. The [dataset authors](https://deep-geometry.github.io/abc-dataset/)
describe chunks of 10,000 models. Attribution: Koch et al., *ABC: A Big CAD
Model Dataset for Geometric Deep Learning*, CVPR 2019.

NYU labels the dataset MIT. The
[author repository licence](https://github.com/deep-geometry/abc-dataset/blob/392789113158b874bfa630098f5bbab3ed38309e/LICENSE)
is MIT, copyright 2019 Deep Geometry Processing. This run uses the dataset
MIT grant for aggregate analysis only. Per-model source exceptions were
**not individually verified**. This is not a per-model licence clearance.

The authors retain creators' CAD copyrights and refer to
[Onshape terms section 1.g.ii](https://www.onshape.com/en/legal/terms-of-use).
Those terms preserve some documents' greater reserved rights in `LICENSE`
tabs. Their effective date, 15 July 2020, follows the dataset release. A
metadata-only audit found that all 10,000 chunk-0000 records were public
documents created before 7 August 2018, with 533 dated 2014 and 9,467 dated
2015. None contained a licence or document-element-list field. One
unauthenticated element-list request returned HTTP 403. The applicable
individual exceptions remain unresolved.

Only counts, typed outcomes and public ABC model ids are published. No model
file, geometry, render, source excerpt or archive-member path is included.
Raw evidence retained privately. The metadata archive and extracted metadata
were removed after the audit.

## Selection and import boundary

The downloaded archive is exactly one STEP chunk, `abc_0000_step_v00.7z`:
1,594,129,754 bytes, with published MD5
`695388be7a278c7798c8c8ae239772ac`. Its actual SHA-256 and the observer,
collector, syntax parser and executable hashes are in the
[run manifest](R8-public-step-census.json).

The collector sorts the chunk's public model ids, then calls
`random.Random(20261111).sample(ids, 3000)`. Results are sorted by model id.
The published CSV uses LF line endings; compare parsed rows when replaying
the pinned measurement revision, whose writer used CRLF.
The sampling unit is a distinct ABC model id. Exact byte duplicates remain
in the sample and their aggregate counts are reported. This estimates this
chunk's model-id population, not the entire million-model corpus. It does
not detect geometric duplicates or establish a distribution of users' work.

Production geometry is the local-main snapshot
`5d3bad1c423cca27bf25d8ece119e9a927e5a66d`. The manifest pins the tooling
commit used for this run. The change adds an example and a collector; it
does not modify production import code.

The corrected observer first scans the entire file for DATA declarations,
outside comments and quoted strings. Only one supported literal `DATA;`
section can qualify for complete coverage. Additional sections, named or
otherwise unsupported declarations, incomplete lexical scans and literal
parser markers inside strings make completeness and absence unknown.
The production parser itself remains unchanged and reads only its first
literal DATA section. This guard does not extend its format support.

For each input, the observer uses that production Part 21
syntax parser to traverse references from B-rep, shell and shape
representation roots. It counts each referenced face entity once, including
faces whose geometry the production recognizer would refuse. Unreferenced
faces, missing references, unresolved face layouts and assembly/mapping
entities are counted separately. Unknown representation items make coverage
incomplete. Instances are not expanded. Completely traversed source faces
can still have unsupported carriers, trims or topology.

The original records lack the new `data_sections`, `unsupported_data_sections`
and `data_section_scan_complete` fields. Their blank CSV values mean
unrecorded, not zero additional sections. The collector treats missing
section metadata as incomplete. Observed spline presence can be retained;
missing section coverage can never establish spline absence.

The analytic carrier families are plane, cylinder, cone, sphere and torus.
`all_face_kinds_buildable` means every counted face has one of these carriers
and source coverage is complete. This measures carrier-family reach. It does
not certify admissible bounds, cone variants, seam curves or solid topology.
Spline counts distinguish ordinary and rational B-spline surface entities.
A known spline face establishes spline need even if another part of the
traversal is incomplete. Otherwise incomplete coverage or an unknown carrier
leaves spline need unknown.

The observer then calls the public recipe import entry,
`import_step_bodies`, with the actual input hash and a fresh `Model`.
It keeps the numeric decider defaults and uses the existing eager application
build allowance of 20,000,000 work units. An exhausted allowance discards the
result and reports `resource_limit`. There is no repair or tolerance change.
A successful result means the entire entry returned bodies; partial builder
mutation is not admission.

The main run uses 12 processes. Each process handles one model at a time.
A separate 60-second observation timeout can interrupt a process; it is not
a geometric refusal or a speed target. Source counts are flushed before
import, so completed source traversal survives an import timeout. Missing
source counts remain unknown. Only the first import refusal per model is
tabulated. Removing it would not prove that a later import will succeed.

STEP B-rep does not provide feature history sufficient to identify fillet
and chamfer operations reliably. Both counts use the typed observation
`unknown_no_feature_history_in_step`. A cylindrical or toroidal face does
not establish a fillet, and a planar face does not establish a chamfer.
These unknowns are not zero counts. R8's feature-provenance measurement
remains open.

## Reproduce

Requirements: Linux, Python 3.11 or later, `bsdtar`, `curl` and Rust. Set
`CENSUS_SCRATCH` to a private scratch directory outside the checkout, on a
volume with at least 30 GB free. Use the corrected tooling revision in
`coverage_reassessment.source_revision`; the top-level identity records
the original measurement.

```sh
CENSUS_ROOT="$(git rev-parse --show-toplevel)"
: "${CENSUS_SCRATCH:?Set a private scratch directory outside the checkout}"
mkdir -p "$CENSUS_SCRATCH"
python3 -c 'import shutil,sys; assert shutil.disk_usage(sys.argv[1]).free >= 30_000_000_000' "$CENSUS_SCRATCH"
export CARGO_TARGET_DIR="$CENSUS_SCRATCH/target"
cargo test --manifest-path "$CENSUS_ROOT/Cargo.toml" -p catallax-import --example r8_census --release
cargo build --manifest-path "$CENSUS_ROOT/Cargo.toml" -p catallax-import --example r8_census --release
python3 -B -m unittest discover -s "$CENSUS_ROOT/tools/research" -p test_r8_step_census.py -v
curl --fail --location --max-time 1800 \
  https://archive.nyu.edu/rest/bitstreams/88598/retrieve \
  -o "$CENSUS_SCRATCH/abc_0000_step_v00.7z"
python3 -B "$CENSUS_ROOT/tools/research/r8_step_census.py" \
  --scratch "$CENSUS_SCRATCH" \
  --binary "$CARGO_TARGET_DIR/release/examples/r8_census" \
  --count 3000 --seed 20261111 --workers 12 --timeout 60 --validate-workers 24
cargo clean --manifest-path "$CENSUS_ROOT/Cargo.toml"
```

The collector verifies archive size and MD5 before extraction and checks
free disk again. It extracts only selected members. Its `finally` block
deletes the archive and extracted models, including after a collection
exception. A forced process kill can bypass cleanup, so inspect the private
`cleanup.json` and remove any remaining copies before reporting completion.
Use a fresh scratch directory for each run. The public CSV and JSON remain
there for review; publish only those two files.

To reproduce the conservative reassessment from retained private records,
without downloading or reopening models, run:

```sh
python3 -B "$CENSUS_ROOT/tools/research/r8_step_census.py" \
  --scratch "$CENSUS_SCRATCH" --reassess
```

The [Rust observer](../../crates/catallax-import/examples/r8_census.rs) has
controls for a clean box, counting through a refusal, ignoring unused
surfaces, multiple or unsupported DATA sections and incomplete coverage claims. The
[collector tests](../../tools/research/test_r8_step_census.py) check
archive-order-independent selection, unsafe/ambiguous member rejection,
denominators, legacy section metadata and retaining face counts through an observation timeout.
