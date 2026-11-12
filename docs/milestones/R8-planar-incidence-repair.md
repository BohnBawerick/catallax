# R8 recorded planar source repair

Date: 12 November 2026. This reruns the same 745 models whose saved planar
refusal pairs have inconsistent plane/boundary data in the
[source audit](R8-planar-self-intersection.md).

**89 models now import with complete coverage through recorded repair.**
The remaining 656 return typed refusals. The largest certified source-to-result
upper bound is 0.000063796781 mm, or 0.063796781 micrometres, below the unchanged
10 micrometre import cap. No source data is rounded away.

This measures the public `apply` entry with no supplied approvals or declines,
a fresh model and the existing 20,000,000-unit allowance. The unrepaired
`import_step_bodies` entry retains its exact-incidence requirement. These are
recorded repairs, not newly proved exact interpretations of inconsistent data.

The [per-model table](R8-planar-incidence-repair.csv) contains all 745 source
hashes, outcomes, first-refusal stages, work counts and attained bounds.
The [manifest](R8-planar-incidence-repair.json) records the measured code
revision, source/build hashes, seeded sample and worker comparisons.

## Outcomes and remaining proof gaps

| Outcome or first refusal | Models |
|---|---:|
| Admitted with complete source coverage | 89 |
| `planar_repair_not_verified`: a shared vertex depends on a curved boundary | 468 |
| `planar_repair_not_verified`: a neighbouring feature needs a complete cover | 120 |
| `planar_repair_not_verified`: incident plane constraints have no exact solution | 33 |
| `planar_repair_not_verified`: the projected source lacks a convex single-bound cover | 31 |
| `planar_repair_not_verified`: shared edge geometry would disagree between faces | 3 |
| `self_intersection_not_verified`: a face pair remains undecided | 1 |
| Resource limits, process errors or external observation timeouts | 0 |

All 120 nearest-feature refusals lack a complete feature cover. None is a
measured violation of the tenth-of-feature limit. No model in this cohort
returns the over-10-micrometre refusal; synthetic controls exercise that limit.
The largest observed work count is 598,284 units.

Refusal measurements retain their scope: 504 are vertex/plane residuals,
31 are source-boundary correspondence bounds, and 120 are whole-face source
bounds with missing feature evidence. The remaining intersection refusal has
no such measurement. Point residuals do not establish whole-face deviation or
permission to import. These numbers do not certify the unobserved later stages
of refused models.

## Construction and source definition

Try refitting the affected planar carriers through unchanged source vertices.
If that cannot be certified, try projecting each shared vertex onto the exact
intersection of its incident planes. Refit straight boundary carriers where
needed. Do not move vertices used by curved boundaries, collapse edges or join
distinct source coordinates. Shared edge occurrences must still agree.

An inconsistent source does not define one exact trimmed face. The `/3` repair
record therefore names its interpretation, `projected_vertex_boundary/1`:
the original boundary vertices projected orthogonally onto the declared plane.
The raw vertex fan and finite original line pieces remain separate retained
constraints. Original bytes and their hash stay immutable.

Complete convex source/result covers prove two-sided geometric Hausdorff
deviation. Corresponding affine maps bound the raw fan and finite line pieces
against the result. Take the maximum of these directly measured bounds.
Do not add intermediate errors or clamp them to a limit. Recompute final
bounds against the immutable source after all repairs. Every changed face
must meet 0.01 mm and a tenth of the source/result neighbouring feature bound.
Topology, closure, orientation, positive volume and exact self-intersection
checks precede publication.

The [interface contract](../interface-changes/2026-11-12-planar-source-incidence-repair.md)
defines the construction order, bounds, typed refusals and record migration.
The record separates measured source-data discrepancy from geometric
deviation and vertex movement. It does not invent file-supplied uncertainty or
an intended design behind the file.

## Attained bounds on admitted faces

There are 6,706 changed-face records: 6,687 use plane refitting and 19 use
shared-vertex projection. The table gives nearest-rank percentiles of certified
**upper endpoints in millimetres**, not estimates of actual errors. The existing
outward enclosure grid is 0.000000000001 mm. Source discrepancy, geometric
deviation and movement are distinct measurements.

| Upper endpoint | Minimum | Median | 95th percentile | 99th percentile | Maximum |
|---|---:|---:|---:|---:|---:|
| Measured source discrepancy | 0 | 1e-12 | 0.00002726459 | 0.000038233552 | 0.000064263898 |
| Two-sided geometric deviation | 0 | 1e-12 | 0.000021059956 | 0.000038150409 | 0.000063796781 |
| Direct bound over retained source constraints | 1e-12 | 1e-12 | 0.00003033301 | 0.000038150409 | 0.000063796781 |
| Vertex movement | 0 | 0 | 0 | 0 | 2e-12 |

The direct source-to-result upper endpoint becomes the local face tolerance.
These are admitted-face bounds only. They do not bound refused faces, later
modified geometry, output rounding or a whole solid against an intended design.

## Regression and determinism evidence

The original 30-model sample uses seed 20261111. Add the first admitted model
for each construction, giving 32 checked models. Repeat them with 1, 2 and 4
workers and matching `CATALLAX_THREADS`. All output bytes, refusals, recorded
bounds and work counts match the main run. This is a sample comparison, not a
full-cohort thread sweep.

Both fresh CLI replays finish all 477 pinned requests at their existing
100,000,000-unit allowance. The ordered comparison reports `equal: true`, with
no manifest or strict STEP-byte differences. The
[replay summary](R8-planar-incidence-replay.json) records input, executable and
revision hashes. The earlier unequal progress counts came from different
start times; neither run omitted requests.

Local checks pass 144 importer tests, nine application tests, eleven observer
tests and five collector tests. Strict Clippy passes for importer/application
targets, and formatting passes. Controls cover the inclusive 10 micrometre
limit, a just-over-limit refusal, the feature limit, all projection ranks,
shared-corner repair, exact fractional coordinates, record reopening, declines,
forged bounds and budget exhaustion without partial publication. The fourteen
existing self-intersection tests still pass, including touching contacts and
true crossings. A reproduced rebuild rollback bug is also covered: a forged
record leaves the destination model unchanged.

The exact all-pairs scan remains quadratic. This slice qualifies capability
and bounded repair; it does not qualify interactive latency.

## Reproduce and data handling

Use the pinned ABC chunk and the retained 745-model selection and pair-details
mapping from the [original census](R8-public-step-census.md#reproduce).
Check free disk space before extraction. Source models and raw geometry stay
outside the checkout. Their attribution and unresolved per-model licence scope
remain those of the original census.

Build `r8_planar_repair`, then run
`tools/research/planar_incidence_repair.py` with runtime paths for the selected
models, selection, pair details, previous cohort records, executable and private
output directory. The collector checks every source hash against the previous
census, uses all available worker cores by default, and retains all first
refusals. It has no external wall-clock admission gate. The public CSV uses LF
line endings and excludes coordinates, private paths and raw diagnostics.

Raw evidence retained privately. The downloaded archive and selected model
files are removed after evaluation. No dataset geometry is committed.
