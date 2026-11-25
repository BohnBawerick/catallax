Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Hollow oblique draft bodies

The B-rep draft operation now builds a separated hollow coaxial source about
an oblique neutral plane. It consumes the existing finite separation bound,
shares four physical spline rings between two annular caps and two walls,
and retains each source face's origins. The inner wall has reversed material
sense. No disk fills either hole opening.

`DraftedCylinder::check` accepts this four-face body as well as its existing
three-face solid family. It checks the common frame and heights, strict radial
separation throughout the finite height, exact boundary rows, cap domains and
opposite edge uses. Membership excludes the hole, including its axis on the
cap planes. Volume subtracts the inner finite law from the outer law.

The B-rep display uses a common tensor grid for both walls and annular cap
strips. It checks the requested chord bound, rounding, facet orientation and
closed edge incidence. Its combined displacement must remain below a positive
three-dimensional wall-separation bound. It retains the existing typed refusals for
unsupported angular and edge-length controls.

For inner residual `F = rho - A(z) - lean.dot(ray)`, let `r` be the minimum
inner radius and `g` the minimum signed radial gap over the finite height.
A segment shorter than `r/2` stays at least `r/2` from the axis.
There the full gradient norm is bounded by `1 + 2*|lean|/r + |tangent|`.
Thus `min(r/2, g/(1 + 2*|lean|/r + |tangent|))` bounds the Euclidean wall
distance even when the two points have different axial heights.
The display and STEP displacement checks both use this spatial bound.

STEP writes each wall as two rational charts. Both annular caps share the wall
rows. The file has six faces, twelve edges and eight vertices, with a complete
file-face map back to the four model faces. Independent rounding of both walls
must satisfy their individual bounds and a shared wall-separation bound before
the writer emits any bytes. This preserves an exact model with a separate
bounded file displacement. It introduces no numerical modelling result or
larger tolerance.

This follows OCCT's shared replacement surfaces, edge curves and pcurves in
[`Draft_Modification`](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKOffset/Draft/Draft_Modification.cxx).
The [P04 SMLib consumer map](../phase1/P04/PORT-CONSUMERS.md) identifies no
missing tracer call needed by this family. Its four exact cap intersections
already exist, so assembly reuses them without a new intersection algorithm.
For display, OCCT's
[`BRepMesh_DelaunayDeflectionControlMeshAlgo`](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKMesh/BRepMesh/BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx)
checks three-dimensional facet deviation and refines failed triangles.
This family keeps its existing whole-domain derivative grid bound and adds the
analytic spatial margin before rounding, so no new tracer or sampled
distance search is needed from the P04 SMLib map.

Field compilation retains its previous separation diagnostic and does not yet
publish the oblique annulus. Touching, crossing, reversed walls and a zero
neutral-axis rate retain their existing refusals. Zero-angle and axial draft
keep their existing paths. No recipe version or application API changes.

The public input is
`crates/catallax-features/tests/fixtures/draft-annular/oblique-bound.json`.
The native regression is in `draft_oblique.rs`. Independent unhealed STEP,
polar-volume and material witnesses are in
`tools/step-regressions/test_draft_oblique.py`.
`draft_annular_clearance.rs` adds a steep thin annulus with radial slope 100,
minimum radial gap 0.0495 mm and an exact cross-height wall pair only
0.000495 mm apart, including a translated, tilted frame.

## Verification

The independent OCCT suite passed all five annular controls without healing.
Each file passed `BRepCheck_Analyzer`, matched an independent polar-volume
integral within `1e-7 mm^3`, and preserved the material and hole witnesses.
Each B-rep display had closed, consistently oriented triangle incidence and
volume within `0.5 mm^3` of the integral.

The fresh CLI runs used the unchanged 20M allowance. These work counts include
the complete recipe and STEP export. Draft and STEP are separate recorded
stages; the total also includes the sketch and request bookkeeping.

| Control | Draft work | STEP work | Total work | File displacement bound |
| --- | ---: | ---: | ---: | ---: |
| Both walls | 4,139 | 857 | 5,023 | 0 mm |
| Outer only | 3,951 | 1,073 | 5,051 | 0 mm |
| Inner only | 3,975 | 857 | 4,859 | 0 mm |
| Negative angle | 4,211 | 857 | 5,095 | 0 mm |
| Outer radius edited to 9 mm | 4,139 | 857 | 5,023 | 0 mm |

All five models are exact. These inputs also have terminating rational file
coordinates, so their file displacement bound is zero. Other admitted inputs
retain the writer's separate displacement check. These are fresh kernel
construction and export runs, not SolidDraft interaction or cold/live export
latency qualification.
