# Coupled planar/curved numerical import

The fixed census parts `00001818` and `00009734` previously stopped at
`planar_repair_not_verified`, stage `curved_dependency`. Both now return closed
numerical B-reps, take a pose rebuild and write approximate STEP. Direct raw
OCCT checks validate the source-pose and edited files without healing,
SameParameter or sewing. This qualifies these actual files, not a predicted
number of admissions from the 791-file diagnostic group.

## Method and limits

Keep shared 3D source vertices and analytic carriers; rebuild and numerically
check finite face parameters under the fixed 0.01 mm tolerance. Do not force a
shared vertex independently onto one plane and then require exact incidence on
its curved neighbour. Successful certified imports retain their existing path.

This follows OCCT's tolerance treatment in
[`ShapeFix_Edge::FixVertexTolerance` and `FixSameParameter`](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V8_0_0/src/ModelingAlgorithms/TKShHealing/ShapeFix/ShapeFix_Edge.cxx):
measure discrepancies, retain shared topology, compare parameterized candidates
and carry the attained tolerance to the affected edge and vertices. The
[`ShapeFix_Wire::FixEdgeCurves`](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V8_0_0/src/ModelingAlgorithms/TKShHealing/ShapeFix/ShapeFix_Wire.cxx)
sequence also checks/rebuilds face curves. The
[P04 SMLib consumer map](../P04/PORT-CONSUMERS.md) identifies no required import
caller of its missing projected-drop path for this analytic family. No new
spline tracer, surface intersection or Boolean algorithm is added.

The first family uses planes/cylinders and finite lines/circular arcs. It checks
shared topology, oriented finite parameter loops, endpoint distances and nine
curve/face locations per use. Planar circular boundaries use sampled projected
walks; cylinder charts must cover less than one turn. Full periodic bands and
singular vertex-loop faces remain explicit refusals. A measured feature-limit
violation is preserved; a missing exact certificate may use numerical checks.

The body and its STEP stay numerical even when every sampled residual is zero.
Pointwise check maxima are not complete source-to-result or whole-body bounds.
Global self-intersection freedom and continuous between-sample separation remain
unestablished. The existing exact source check still refuses a crossing it finds;
its unavailable exact-incidence prerequisite is replaced by the stated numerical
checks. The finite planar STEP writer retains its existing path and controls.

The analytic STEP helper writes explicit cylindrical pcurves and checks actual
decimal geometry, including STEP frame normalization and projected references.
Actual decimal pcurve mappings have separate distance samples. Length rounding,
3D direction-ratio rounding and parameter-coordinate rounding retain their own
scopes. No STEP global uncertainty or complete-body bound is invented.

## Actual real controls

[real-controls.json](real-controls.json) pins source bytes, the actual producer,
executable and resulting files. The original source observation is reproduced
with the frozen census binary; both input sources translate raw in OCCT, but its
analyzer considers those raw source shapes invalid. The positive result here is
the newly generated, independently checked output, not a source-validity claim.

Both accepted numerical bodies and both pose rebuilds contain one closed shell,
six faces, twelve edges and eight vertices. Pose changes preserve volume and
translate the measured centroid by 3, 5 and 7 mm. Each edited body retains one
history stage. The original source bytes stay unchanged.

At producer `2f2e5b8b02481f80d3286eddae16ccea9eeef78b`, work is:

| Part | Import | Source STEP | Pose rebuild | Edited STEP | Whole qualification |
|---|---:|---:|---:|---:|---:|
| 00001818 | 143,428 | 190,858 | 99,899 | 190,585 | 624,770 |
| 00009734 | 144,895 | 207,699 | 98,955 | 207,677 | 659,226 |

All stages share the unchanged 20M request allowance. The third fixed source,
`00002345`, still returns `import.numerical-unsupported-family`; no success is
claimed for its complete periodic boundary.

Authored quarter-cylinder controls have volume `pi/4`. Their actual source/pose
files pass raw OCCT validity, topology and material probes. A subsequent OCCT cut
also validates, with volume `sqrt(3)/8 + pi/12`; this is reference-kernel
continuation evidence, not a Catallax Boolean claim. Negative controls cover
over-tolerance residuals, smaller requested tolerances, winding, periodic cover,
distant periodic coordinates, exhausted requests and changed output bytes.
Existing near-feature, source-cap, crossing, rollback, planar and oblique controls
remain in the scoped suite.

## Reproduction and census

Build `catallax-import`'s `numerical_import` example in release mode and run each
immutable private `MODEL_ID.step` with a new external output directory. The
response records source admission, actual checks, stage work, pose history and
both disk STEP hashes. Validate the actual artifacts with the raw OCCT builder
used by [N03](../N03/REPORT.md) and
[the earlier oblique qualification](../import-oblique-step/README.md).

The restored `r8_census --numerical` observer and
`tools/research/import_census.py` rerun exactly the original 3,000 model IDs,
seed, individual hashes, 20M allowance and external 60-second cap. The collector
also pins the prior complete 9.20% table, retains its outcomes/first causes and
records actual transitions. It keeps full/partial/refused/timeout distinct;
missing results or process errors prevent a complete publication. Its output is
admission evidence and does not assert an edit or valid STEP for every admission.
The full rerun and resulting before/after table are a separate publication step.

Raw source geometry stays private under the
[existing licence scope](../../milestones/R8-public-step-census.md#source-and-licence-scope).
The SolidDraft replacement/reference and numerical-group work remains open in
the app slice.
