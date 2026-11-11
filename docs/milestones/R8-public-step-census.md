# R8 public STEP census after the placement fix

Date: 11 November 2026. Same 3,000 ABC model ids, seed `20261111`.

**Admissions are unchanged.** The import entry accepts the same 189 model
ids as the [merged baseline](https://github.com/BohnBawerick/catallax/pull/415).
Of these, 186 have complete source coverage, or 6.20% of the sample; three
retain unverified representation items. All admitted ids and accepted body
counts are unchanged.

The 982 placement-perpendicularity diagnostics disappear after the
[placement fix](https://github.com/BohnBawerick/catallax/pull/425). Those models
now reach 964 later refusals and 18 external observation timeouts, with no
new admission. The first-refusal distribution changes substantially even
though the admitted set does not.

The [per-model table](R8-public-step-census.csv) and
[run manifest](R8-public-step-census.json) are the fresh run. Baseline source
and data remain available in the merged baseline linked above. Fillet and
chamfer provenance remains unknown for all 3,000 models.

## Side-by-side import outcomes

The baseline uses main `d6ea026f00eb269467c282c9cae657cfa85e3beb` with
source/tooling snapshot `8f029060e6da8a94d2df6598ec3e7a177d2592b0`.
This run uses main `c6ba91fda6bebee6b0bdb98c4e254af6b7cab13c` with
source/tooling snapshot `dd2b3e1411010d455d94266a2867500b8275ec57`.
Production crates match the respective main snapshots. This follow-up
changes evidence files, not import code or allowances.

| Outcome | Merged baseline | After placement fix | Change |
|---|---:|---:|---:|
| Entry accepted, complete source coverage | 186 | 186 | 0 |
| Entry accepted, source coverage unverified | 3 | 3 | 0 |
| Returned typed refusal | 2,797 | 2,782 | -15 |
| External observation timeout | 14 | 29 | +15 |
| Total evaluated | 3,000 | 3,000 | 0 |

The CSV's `admitted` outcome combines the first two rows, giving 189 in
each run. The actual import entry returned bodies; `coverage_complete`
separates the three whose source representation remains unverified.
Every evaluation uses a fresh model. Partial builder mutation is not admission.

All 189 baseline admissions remain admitted with the same body counts.
Of the baseline refusals, 2,779 still return a refusal and 18 now time out.
Of the baseline timeouts, 11 still time out and three now return a refusal.
These are bounded observations on a shared machine, not a timing qualification
or a claim that the placement fix caused a speed regression.

## First refusal causes

Counts below are first returned causes, sorted by the new count. Removing
one cause exposes later checks; it does not establish admission.

| First refusal cause | Merged baseline | After placement fix | Change |
|---|---:|---:|---:|
| `self_intersecting` | 587 | 822 | +235 |
| Unsupported spline edge solid certification | 240 | 532 | +292 |
| `resource_limit` | 255 | 334 | +79 |
| `not_solid`, exact incidence checks | 178 | 203 | +25 |
| Unsupported spline face solid certification | 103 | 196 | +93 |
| Unsupported `SURFACE_OF_LINEAR_EXTRUSION` | 141 | 193 | +52 |
| Unsupported `SHELL_BASED_SURFACE_MODEL` | 91 | 116 | +25 |
| `self_intersection_not_verified` | 72 | 85 | +13 |
| Unsupported cone variants | 30 | 81 | +51 |
| Unsupported `SURFACE_OF_REVOLUTION` | 29 | 61 | +32 |
| Unsupported torus variants | 21 | 58 | +37 |
| Unsupported `BREP_WITH_VOIDS` | 33 | 43 | +10 |
| Unsupported `OFFSET_SURFACE` | 22 | 31 | +9 |
| `malformed` during syntax/geometry recognition | 993 | 24 | -969 |
| Unsupported vertex loop beside other bounds | 2 | 3 | +1 |
| Total returned refusals | 2,797 | 2,782 | -15 |

Unsupported geometry/representation refusals total 712 before and 1,314
after. The manifest/CSV call the offset bucket `unsupported:other_step_type`;
all 31 new-run diagnostics in that bucket identify `OFFSET_SURFACE`.

Of the former 982 placement cases, the new first causes include 292 spline
edge certification refusals, 234 reported self-intersections, 93 spline face
certification refusals and 78 resource limits. All 982 are accounted for by
964 returned refusals and 18 timeouts. The earlier perpendicularity message
is absent; no parallel-placement diagnostic is observed either.

The 24 `malformed` returns comprise 22 non-increasing-distinct-knot
diagnostics, one missing `MANIFOLD_SOLID_BREP`, and the same syntax failure
at model id `00003087`. Thirteen knot diagnostics were previously hidden
behind the placement check. The 203 `not_solid` returns comprise 167 edge
vertices not exactly on their curves, 22 edges not exactly on their faces,
and 14 vertex-loop points not exactly on their faces.

Cone refusals comprise 35 non-rational-axis-length cases, 29 incomplete rim
bounds, seven off-axis vertex loops and ten incomplete circles about the
axis. These labels describe Catallax's diagnostics, not an independent
judgement that every refused source model is invalid.

The 334 resource-limit returns, 85 unverified-intersection returns and 29
external timeouts remain incomplete observations or proofs. They are not
impossible-geometry claims. Twenty-seven timeouts retain complete source
counts; two have no source checkpoint. Numeric and wall observation caps
are the same as the merged baseline.

## Source coverage and face kinds

| Observation | Merged baseline | After placement fix |
|---|---:|---:|
| Complete source traversals | 2,883 | 2,884 |
| Every face has an analytic carrier, complete traversal | 2,002 | 2,003 |
| Analytic-only share of all 3,000 models | 66.73% | 66.77% |
| Analytic-only share of complete traversals | 69.44% | 69.45% |
| Non-analytic carrier, complete traversal | 881 | 881 |
| Incomplete source coverage | 117 | 116 |
| Spline face observed, including partial traversals | 868 | 868 |
| Spline absence established | 2,002 | 2,003 |
| Spline need unknown | 130 | 129 |
| Observed referenced face entities | 704,829 | 851,349 |

This is the identical input sample, with 3,000 distinct model ids and byte
hashes. One previously missing source checkpoint now completes: model id
`00009272` contributes 146,520 planar faces and still reaches the import
observation timeout. That recovered observation accounts for the extra
analytic-only model and every added observed face. It is not a new input,
a new supported family or an admission gain. Counts match exactly wherever
both runs complete syntax parsing.

The full-file scan records one supported DATA section in 2,998 models; two
have no source checkpoint. No additional or unsupported DATA section is
observed among those 2,998. Syntax parses in 2,997 models, fails in one and
is unobserved in two. Incomplete source coverage comprises 113 models with
unknown representation items, one syntax failure and two missing checkpoints.
Completed observations have no missing references, unreferenced faces,
unresolved face layouts or assembly/mapping entities.

The observed graphs contain one `MANIFOLD_SOLID_BREP` in 1,948 models,
several in 1,013, none in 36, and no completed graph in three. The sample
unit is models, not individual solids or manufactured parts. It samples
one 10,000-model chunk, not the million-model corpus. Exact byte duplicates
are absent; geometric duplicates are not detected.

| Observed face carrier | Merged baseline faces | New faces | Share of new observed faces |
|---|---:|---:|---:|
| Plane | 354,968 | 501,488 | 58.91% |
| Cylinder | 166,637 | 166,637 | 19.57% |
| B-spline | 54,799 | 54,799 | 6.44% |
| Rational B-spline | 52,530 | 52,530 | 6.17% |
| Torus | 24,770 | 24,770 | 2.91% |
| Cone | 22,329 | 22,329 | 2.62% |
| Other | 21,324 | 21,324 | 2.50% |
| Sphere | 7,472 | 7,472 | 0.88% |

Face counts deduplicate shared entity references within a model and include
observed faces of incompletely covered models. They do not estimate unseen
geometry. The changed face-frequency denominator follows the recovered
source checkpoint, not changed model geometry. Part-level carrier reach is
therefore reported separately.

The five analytic carriers are plane, cylinder, cone, sphere and torus.
Their presence does not certify trims, supported boundary curves or closed
solid topology. `other` is a carrier outside the listed families. Known
spline presence survives incomplete traversal; unobserved geometry or an
unknown carrier cannot establish spline absence. Fillet and chamfer fields
remain `unknown_no_feature_history_in_step` for every model. Surface
curvature alone does not establish feature provenance, so these are not zeros.

M24 now has measured later refusals for the formerly blocked placement
cohort. Spline solid certification, reported intersections, exact incidence
and incomplete proofs remain the main import boundaries. M21 receives a
carrier-family candidate set, not field-construction or independent
kernel-agreement evidence. No source-validity oracle was used.

## Method, validation and cleanup

The collector sorts public ABC model ids and selects
`random.Random(20261111).sample(ids, 3000)`, then sorts results by id.
Archive size, MD5 and SHA-256 match the merged baseline, as do all selected
member identities. There are no substituted or omitted selected models.

The observer scans the entire file before making completeness or absence
claims. Qualification requires one supported literal `DATA;`, a complete
lexical scan and no hidden parser marker in a string. Additional or
unsupported sections make completeness and absence unknown. The production
Part 21 parser still reads its first literal DATA section. The observer
traverses B-rep, shell and representation references independently of the
first import refusal; unknown items and transformations prevent a complete
coverage claim.

Import uses `import_step_bodies`, the actual input hash and a fresh `Model`.
Numeric decider defaults and the existing 20,000,000-unit application build
allowance are unchanged. Exhausted work discards a returned result as
`resource_limit`. No repair or tolerance change is applied. Each process
has the same separate 60-second observation cap as the baseline. Source
counts are flushed before import. These clocks are not M26 speed gates.

The main run uses 12 processes. Repeating the smallest 24 selected ids with
1, 2 and 4 workers gives identical public results and complete private
stdout/stderr diagnostics, including work counts. This is subset determinism
evidence, not a full-sample worker comparison.

Six Rust controls, six Python tests and strict example Clippy pass on this
main snapshot. The controls include the multiple-DATA-section regression.
An independent audit checks source and executable hashes, sample identity,
CSV/manifest totals, section qualifications, actual refusal codes, admitted
ids/body counts and work limits. The saved 477 replay inputs/answers are
unmodified and were not rerun by this evidence-only follow-up.

The archive and extracted models were deleted after evaluation, and absence
was checked directly. Build output was removed with `cargo clean`. The run
started with 65 GB free and had 63 GB before extraction. Raw evidence retained privately.

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
