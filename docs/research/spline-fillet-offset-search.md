# Spline fillets: the offset search, and the route that replaced it

Outcome: route changed. The direct search did not finish. A different route rounds one loft edge within 59 nm.
When: 11 to 12 November 2026.

## What we were trying to achieve

A fillet between two spline faces, on a loft and on an imported part, with the whole fillet face proved within 100 nm of its rolling-ball definition (hard limit 1 µm).

The centre of the rolling ball runs along the curve where the two faces, each moved inward by the radius, meet. Everything else in the fillet hangs off that curve, so the question was how to find and prove it.

## Why it was proposed

The spline-fillet milestone, now 2.1, was about to commit to fillets on arbitrary spline faces. The port scope study had said to scope free-form fillets separately from intersections. So a scope study was added on 11 November as a bounded experiment: does the cell search for spline cuts, which works for cuts, extend to the fillet's equations?

## What was tried

1. The direct extension: five equations in six unknowns, on one edge of a ruled loft.
2. A follow-up after the owner chose the 2026 engine: enclose the two normal lengths instead of searching them, which cuts the search to four unknowns. Tried on the same loft edge, the next edge and a curved control.

## Result

1. The six-unknown search refuses at proof depth 24 after 2,910,656 of 20,000,000 units. It does not run out of work. It runs out of depth.
2. The four-unknown search also refuses at depth 24 on all three edges. The loft edge uses 5,431,053 units.

No fillet face came out of either.

## Why the route was changed

The verdict in the note: this is evidence against that search, not proof that every exact offset search is too dear. It supports changing the proof route. It does not support raising the allowance, narrowing a tolerance or dropping small pieces.

The route that worked, from 12 November: give the spline-pair engine the exact offset faces as input, prove the centre curve there, then build the fillet face on the proven curve.

- The loft edge rounds at radius 0.1 mm in 10,349,665 units. The face is proved within 59 nm.
- It exports to STEP and opens as a valid solid in OpenCASCADE without repair.
- The same edge on a loft written by OpenCASCADE and imported rounds in 9,158,476 units, within 25 nm.
- Face and edge names survive a radius change.

## What is still open

The builder rounds one convex edge between two single-patch faces whose ends meet planes. Fillets across more than one patch, extra trims and junctions still refuse. That is milestone 2.1.

## Sources

- [M27-SPLINE-FILLET-SCOPE.md](../M27-SPLINE-FILLET-SCOPE.md): https://github.com/BohnBawerick/catallax/pull/400
- [M11-OFFSET-PROOF.md](../M11-OFFSET-PROOF.md): https://github.com/BohnBawerick/catallax/pull/445
- [M11-COVERAGE.md](../M11-COVERAGE.md): https://github.com/BohnBawerick/catallax/pull/521, https://github.com/BohnBawerick/catallax/pull/523, https://github.com/BohnBawerick/catallax/pull/526

## 16 November: numerical multi-patch body and raw STEP

The phase-2.1 continuation uses metered Newton offset contacts across atlas
charts, C1 cubic centre/section fits and paired world/native-UV end cuts.
The sampled tolerance is 1e-7 mm, not a whole-face or whole-body bound.
The perpendicular-plane seam and circular-section controls pass. An earlier
curved-loft carrier/cap-only probe passed at 1e-4 mm in 29,335,696 work; that
is separate from the current, stricter body result.

Declared numerical carriers use sampled chart-box coverage and nine UV-walk
locations per coedge. Only the new fillet carrier is declared numerical in
this probe. Source-carrier checks are unchanged. A shared edge incident to
a declared numerical face uses nine native parameters for its numerical
world-distance evidence, with the same 1e-7 mm tolerance. Exact-exact edges
retain their original clipped 3D/UV knot-span sample grid. The mixed-body
control rejects an exact-exact pcurve spike between the coarse locations;
the exact chart-gap and default native continuation/repair/support-edit/
handover controls also pass.

The smooth three-section loft's fillet body admits at 96,732,996 modeled
work within the 100M request allowance. The earlier 96,751,874 count also
included writing an unfilleted source STEP fixture. Independent circle-radius
and end-plane samples pass. Its measured final-admission work is:

| Check authority | Charged work |
|---|---:|
| Numerical face-local | 1,271,040 |
| Unchanged source-face-local | 6,540,517 |
| Cross-face checks and setup | 168,954 |
| Final-admission total | 7,980,511 |

These counts cover final admission only; the same body-build allowance also
covers the preceding source and fillet work. The STEP continuity fix below
lets the admitted body write and reread numerical STEP in a separate
80,265,942-work request. That is not a combined 100M end-to-end app result.
The kernel entry point is `catallax_features::fillet::numerical::build`;
ordinary recipe evaluation has not been migrated to this numerical result.

OpenCASCADE/OCP 8.0.1.1 reads the actual 644,190-byte STEP through
`StepToTopoDS_Builder`, with healing, `SameParameter` and `TransferRoots`
off. It is valid, with seven faces and one closed solid. A later box cut
is valid and reduces volume from 8.051606364758838 to 6.516691845842079 mm³.
The file SHA-256 is
`6f1984918902e4ee5f9ee2a57593aa33f9330f7b1c8f70ab52cef0c04fe3edaa`.

The [bounded attempt record](../phase5-work-cases.md#incomplete-21-bounded-attempts)
lists radius-edit, import, housing and all 13 frozen later-cut outcomes.
Only the far-disjoint M8 cut builds. No completed 2.1 milestone or stable
radius-edit names are claimed. The loft regression is
`m11_spline_fillet::numerical_fillet_builds_on_the_smooth_loft_across_its_station_seam`.

The interpolation and separate UV/3D-curve representation follow the
[SMLib rail-blend/lift workflow](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmFilletSolver.cpp#L8465-L8520)
and [OpenCASCADE's paired intersection curves](https://dev.opencascade.org/doc/refman/html/class_int_tools___face_face.html).
Numerical face-local and declared numerical-pair admission use Catallax's
finite sampling convention, not either kernel's complete wire certification.
No reference source was ported.

### STEP continuity encoding

Independent decimal rounding can split an exact C1 join whose continuity
is carried by pole relations at a degree-fold knot. The numerical STEP
helper now removes that one redundant pole strip before rounding, only
when the homogeneous inverse-knot-insertion equality is rational and exact.
The remaining knot multiplicity encodes C1 continuity structurally.
Non-removable joins are unchanged; no rounding tolerance is enlarged.

This follows the knot-multiplicity reduction used by
[OpenCASCADE `RemoveUKnot`/`RemoveVKnot`](https://occt3d.com/dev/doc/refman/html/class_geom___b_spline_surface.html),
restricted here to zero-error removal rather than tolerance-based smoothing.
[SMLib's `N_SrfRemoveKnotConditional` and boundary-constrained removal](https://docs.nvidia.com/smlib/manual/nlib/functions/index.html)
were checked as the alternative. No reference implementation was copied.
The decimal C1 regression and the existing `catallax-step` Rust suite pass.
Direct STEP from the admitted loft also passes its numerical reread and
the independent raw OpenCASCADE read/edit above. Optional diagnostic body
handover was removed after its separate 100M refusal; no transfer-engine
performance change or repeated-proof claim is made. Numerical body/history
layout 03 retains the new carrier declarations. Tagged v02 numerical body
and history records still read with empty declarations; unversioned
positional layouts continue to refuse.

### Review corrections

Output admission and endpoint projection use the stricter of the source
and requested fillet tolerances. A loose source therefore cannot admit a
new seam beyond the requested tolerance. This conservatively rechecks the
whole output at that stricter tolerance.

Native bindings are carried by face/boundary/use occurrence. Only replaced
rails, shortened border endpoints and inserted cap uses change; untouched
UV branches are retained, including two uses of a periodic seam edge.
Inserted cap bridges take their endpoint UVs from the two adjacent edited
occurrences, which were evaluated at their retained native edge parameters;
a separate projection of the new world vertex cannot replace those UVs.
Current-main surface kinds, continuous planar-circle clearance and
certified handover migration are preserved.


## 17 November: radius-edit end cuts

The radius change from 0.1 mm to 0.15 mm exhausted its 100M allowance in
the first end-plane cut. The diagnostic retained the original 1e-7 mm
tolerance and reported the typed geometric stage, not just the work limit.

[OCCT's `math_BracketedRoot`](https://occt3d.com/dev/doc/refman/html/classmath___bracketed_root.html)
uses derivative-free Brent interpolation inside a bracket.
[SMLib's fillet solver](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmFilletSolver.cpp#L2997-L3038)
instead sets bounded Newton equations for its coupled contact problem.
The end-cut candidate followed the OCCT bounded-interpolation approach,
using a secant proposal alternating with mandatory bisection rather than
the full Brent scheme. It still exhausted 100M in the first cut, so the
candidate was parked and the original bisection retained. No reference
implementation was copied.

In that experiment binary64 interpolation only proposed the next parameter.
Returned points still passed the same 128-bit plane-residual check, with
unchanged cap fitting tolerance and paired world/native-UV checks. The
measured request split is in the phase-5 record. No successful strict radius
edit or exact intersection authority was established.


## 17 November: finite clearance before general continuation

The everyday-part scout found false overflow refusals, not exhausted work:
the guard compared complete carriers instead of the retained patches.
An infinite cylinder could meet a rounded side below its actual trim. A
boss-base torus of rolling radius `r` reached `R + 2r` in that test, although
its retained strip stopped at its contact radius `R + r`.

The ordinary recipe fillet builder now encloses complete straight-meridian
ring regions by their edge and contact circles. Existing cut hulls still
bound finite straight strips. A spherical corner gets a pyramid containing
its contact cone's sphere sector. Trimmed neighbour boxes come from
`catallax_measure::face_bounding_box`, not the neighbour's complete support.
Disjoint finite bounds prove clearance. Overlap and unavailable bounds stay
unresolved and retain the intersector or continuation refusal.

This follows the face-aware bounding approach of
[OCCT `BRepBndLib`](https://occt3d.com/dev/doc/refman/html/class_b_rep_bnd_lib.html)
and [SMLib `SmFace::CalculateBoundingBox`](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmFace.cpp#L819-L925).
SMLib uses finite edges for linear surface directions and the trimmed UV
domain otherwise. The new corner bound uses its exact positive contact cone
and Cauchy-Schwarz. No sampling or tolerance change enters these clearance
proofs, and no reference implementation was copied.

The new `trimmed_fillet_clearance` regression exercises ordinary feature
evaluation, topology, independent mass or material controls, and STEP:

| Case | Build work | STEP work |
|---|---:|---:|
| Box corner R10 | 1,889 | 89,631 |
| Box corner R15 | 1,874 | 94,614 |
| Box corner R19 | 1,874 | 93,588 |
| Boss base R4 | 13,844 | 1,973 |
| Boss base R5 | 13,588 | 1,973 |
| Boss base R6 | 13,588 | 1,973 |
| D12 hole, axis 12 mm from wall, R3 | 13,258 | 1,903 |
| Chamfered cylinder rim R2 | 30,641 | 21,155 |

All eight emitted files open as valid closed single solids through raw
OCP 8.0.1.1 `StepToTopoDS_Builder`, with healing, SameParameter and
TransferRoots off. Each takes a valid material-removing OpenCASCADE cut.
The owner's pocket rim was not among them; "Native strips between planes
and cylinders" below records it. Demo edge coverage is not inferred from
these eight cases.

The separate native rolling-ball continuation implementation follows
[OCCT `BRepBlend_Walking`](https://occt3d.com/dev/doc/refman/html/class_b_rep_blend___walking.html)
and [SMLib's constant-radius contact and section-plane equations](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmFilletSolver.cpp#L2665-L2715).
Accepted contacts and guide sections receive interval residual checks at
the declared 1e-7 mm tolerance, or 0.1 nm. It is not yet the application's fallback.
The cone rim and three-support ball checks pass; the cylinder rim still
exhausts 100M inside continuation. Further continuation work follows the
ordered clearance, intersection and arc-end repairs.

## 17 November: native angular charts and the ordinary numerical outcome

The native numerical boundary now admits regular spherical charts and reuses
the ring-torus implementation of the analytic numerical families. A full turn
requires a rectangular paired seam cell; repeated edges without a complete
period translation refuse. Sphere latitude stays strictly away from its
poles. Fitted native trims additionally need positive homogeneous weights
and a complete cropped Bézier control hull inside the supported finite
chart. These are local chart checks, not a global body-embedding proof.

The parameter convention follows [OCCT `Geom_SphericalSurface`](https://occt3d.com/dev/doc/refman/html/class_geom___spherical_surface.html):
longitude and latitude in radians. SMLib's
[`SmSphere`](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmSphere.cpp)
and [`SmTorus`](https://github.com/NVIDIA-Omniverse/usd-brep/blob/main/source/SMLib/src/SmTorus.cpp)
also retain the analytic generating circles, but express their angular
intervals in degrees. No degree convention or reference implementation was
copied into Catallax.

After rebasing onto the analytic-family fixes, independent 100M requests at
1e-7 mm produce:

| Native body | Admission work | STEP work |
|---|---:|---:|
| Ring torus, each of four cyclic seam-walk starts | 16,811 | 37,991–37,997 |
| Spherical wedge with four planar neighbours | 87,125 | 110,379 |
| Spherical wedge with fitted native UV trim | 91,129 | 4,470,183 |
| Same wedge with shifted circle parameter frame | 91,129 | 4,470,189 |
| Toroidal wedge with fitted native UV trim | 87,683 | 3,433,161 |

The five iso-curve STEP files open as valid closed single solids through raw
OCP 8.0.1.1, with healing, SameParameter and TransferRoots off. Each takes a
valid material-removing cut. The parent rebase preserves their bytes.
The three fitted-trim files also pass raw validity and subsequent cuts.
The two spherical sources emit identical bytes despite their different
circle parameter frames; the toroidal wedge also retains its concave
cylindrical neighbour. Generic native-model embedding and native fillet
assembly remain outside these carrier proofs.

Both STEP paths share one native circle phase/frame preparation and one
native surface codec. Fitted UV knot ranges move with the written circle
range, following [OCCT `BSplCLib::Reparametrize`](https://occt3d.com/dev/doc/refman/html/class_b_spl_c_lib.html#a3d67bbc71b74521a3648ee3d47d64f11)
and SMLib's [coordinated 3D/pcurve edge reparameterization](https://github.com/NVIDIA-Omniverse/usd-brep/blob/3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5/source/SMLib/src/SmProtoTopology.cpp#L1195-L1224).
The affine knot transformation is rational; decimal carriers and pcurves
still require numerical admission and measured distance checks. Closed
circles with fitted pcurves retain a typed full-period knot-domain refusal.
The shared analytic codec's rational native frame/radius requirement also
remains explicit; numerical inputs are not labelled exact.

Native source-fit coordinates come from 128-bit enclosures, not the cached
display-only `Real::to_f64` estimate. Circle endpoint phases use the same
qualified vertex coordinates as the face UV endpoints. A coarse display
estimate can otherwise make those two descriptions disagree before fitting.
Native UV fit refusals retain the source edge, face, range and residual;
the requested tolerance and full-body admission remain unchanged.

The ordinary `build_incremental_without_display` application API now
returns `EvaluatedApplicationOutcome`: either its exact application body
or an opaque numerical body with its source recipe and face origins.
`build_incremental_with_files_without_display` carries the same distinction
with resolved files. Numerical results have no exact model or reusable
exact checkpoint. The multichart loft fillet consumer passes through this
ordinary route at 96,733,000 work with one unchanged 100M allowance; its
fitted blend is geometrically displaced from the old sharp edge. This does
not prove native-chain fallback, later edits, imported-part fillets, or the
owner's pocket.

The outcome and numerical holder are application-owned types. The producer
moves admitted geometry, the authored recipe and origins into that holder;
it does not export a feature-owned alias, rebuild the body, or convert it
to an exact checkpoint.

## General rolling-ball construction, 17 November

One numerical rolling-ball path replaces the pair/end-recipe construction
gates. The contact stage now groups selected edges into unique
tangent chains using 128-bit enclosures of endpoint tangent and support-normal
gaps, scaled by radius and limited to one eighth of the requested tolerance.
Sharp vertices remain strip ends. An ambiguous tangent continuation returns
`fillet.numerical-chain`, not an arbitrarily chosen branch.

The construction follows [OCCT BRepBlend_Walking's guide-normal sections](https://occt3d.com/dev/doc/refman/html/class_b_rep_blend___walking.html)
and [SMLib's constant-radius four-UV equations](https://github.com/NVIDIA-Omniverse/usd-brep/blob/3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5/source/SMLib/src/SmFilletSolver.cpp#L2665-L2745):
three offset-equality equations plus one section-plane equation. It uses the
same regular-surface adapters for native and spline supports, rather than
selecting a surface-pair recipe. Convex/concave material side and alignment
between parameter normals and source-face normals determine signed offsets.
Retained curve traces supply sampled, fitted guides even when Boolean
intersection carriers have lost explicit line/circle identity. Endpoint
tangent lines continue the section guide outside its retained interval;
they do not assert an extension of the actual trimmed support faces.

The release regression independently checks material side, support contact
and circular-section radius on all 12 edge orientations of a box at
R0.5/R1/R2/R5. The unchanged tolerance is 1e-7 mm. Request work and the
over-20M cases are recorded in [the phase-5 work record](../phase5-work-cases.md#general-rolling-ball-contact-stage-17-november).
That regression proves untrimmed candidate strips and rails. Ordinary
multi-edge numerical continuation now calls the common constructor.
Real-domain G1/two-strip joins, owner bodies and raw STEP acceptance remain open.

The generic neighbour-cut stage follows [OCCT's two-surface walking](https://occt3d.com/dev/doc/refman/html/class_int_walk___p_walking.html)
and [SMLib's traced-curve fitting](https://github.com/NVIDIA-Omniverse/usd-brep/blob/3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5/source/SMLib/src/SmSurfaceIntersector.cpp#L4900-L5085).
It traces the fitted strip against a regular neighbouring carrier with zero
offsets, retaining one world curve and both parameter curves. Refinement
checks their physical agreement with 128-bit surface evaluations at retained
stations and eighth-span samples. These are sampled numerical checks, not a
complete branch census or a proof against the neighbouring face's trim domain.
An independent cylinder/quadratic-height test checks both UV mappings at
1e-7 mm; its request cost is recorded in the phase-5 work record. Real-domain
clipping and closed-body publication remain the next construction stage.

The capped-strip publication primitive uses those two-surface cuts against
the actual end neighbours, not an end-plane recipe. It shortens retained
border edges, replaces the picked support occurrences in their actual loops,
and inserts the cut with both fillet and neighbour parameter curves. Its
topology order follows [SMLib's trim-original-surfaces stage](https://github.com/NVIDIA-Omniverse/usd-brep/blob/3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5/source/SMLib/src/SmFilletExecutive.cpp#L3381-L3400):
intersect, trim the retained faces, assemble shared boundaries, then admit
the complete body. It retains holes and nonzero boundary indices. Numerical
endpoint constraints remain subject to the unchanged physical admission
tolerance; they are not incidence certificates. The one-edge regression
checks all six original planar face domains, including the rounded cap's
area. Network/junction publication and owner acceptance remain separate gates.

Three-support network publication reuses each common ball's canonical
support contacts as shared vertices. Each extracted section is one world
edge shared by a strip and its sphere, with the paired native parameter
curves and recorded orientation kept together. Source loop occurrences are
replaced before the complete network enters numerical solid admission.
This topology stage adds no exactness claim to the fitted geometry.

The R1 regression now admits all 12 strips and eight corners of a
20×30×40 box as one numerical solid. Independent area controls check the
six retained original plane domains. Its whole-request cost is recorded
in the [phase-5 ledger](../phase5-work-cases.md#shared-rolling-strip-and-corner-network-18-november).
This is not yet owner-rim, ordinary application or STEP acceptance.

Retargeting only the end pole of a clipped UV rail can reverse a tiny
retained knot span. The constructor now applies one orientation-preserving
affine map to the complete UV curve, preserving its knots, degree and
rational weights. We follow [OCCT's whole-curve `GeomLib::GTransform`](https://occt3d.com/dev/doc/refman/html/class_geom_lib.html#a9ee76c30e846d5a3ab60198485738bd0);
[SMLib likewise supports nonuniform transforms of base B-spline geometry](https://github.com/tangency-labs/smlib/blob/main/.agents/operations/transforms.md).
Coincident ends with incompatible targets and a reversed target chord
return typed causes. The affine arithmetic is rational; physical
same-parameter checks and the constructed fillet remain numerical at
1e-7 mm. No loop or physical admission guard was relaxed.

### Early native contact limits

A plane parallel to a circular cylinder's axis has straight intersection
generators even when the selected edge is stored as `Curve::Intersection`.
The early check therefore uses the native support geometry, not a Line
carrier test. It reuses the strict affine inequalities for the signed
offset plane and cylinder over every selected chain before strip walking
or fitting, retaining the first proved contact limit.

This follows [OCCT's analytic plane/cylinder offset construction](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/ChFiKPart/ChFiKPart_ComputeData_FilPlnCyl.cxx)
and [SMLib's signed constant-radius offset equality](https://github.com/NVIDIA-Omniverse/usd-brep/blob/3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5/source/SMLib/src/SmFilletStandardSolver.cpp#L168-L230).
It does not invoke an exact fillet producer or gate the general constructor
when the cheap proof is unavailable. A reached face is named only after
the critical original-cylinder contact lies on an actual finite boundary
line coedge shared with that input face, not an untrimmed carrier guess.

The real six-edge owner R5 request now refuses at the exact R4 contact
limit and names its original outer-side face. The ordinary application
retains `fillet.radius-too-large`, `; reaches face #57` and the terminal
`; radius limit 4 mm` hint. Contact proof, source context and the warning
that a contact bound does not prove a fitting radius precede that hint.
The existing one-third-bound test also checks the SDK's terminal parser
without losing the exact proof or the rounded legacy accessor.
Request scopes and costs are in the [phase-5 ledger](../phase5-work-cases.md#owner-radius-overflow-eligibility).

For an open selection on a cylinder and an axis-normal plane, the
shrinking native offset circle has radius R-r. Its plane-contact rail
must cross the actual sharp line border terminating the selected arc.
The strict crossing disappears at R-d, where d is the in-plane distance
from the cylinder axis to that line. The critical foot must lie on the
finite original coedge; a smoothly continuing neighbour, another selected
edge at the endpoint, a periodic edge, or an expanding offset does not
qualify for this proof.

This uses the circular-rim branch of the same OCCT construction
(its inner offset radius `cylrad - Radius`) and SMLib's signed offset
orientation. It is a necessary regular end-contact constraint, not a
claim that the two support offsets have no centres, nor a global largest
fitting radius. The typed refusal states that continuing past this contact
needs different topology. Naming the reached face uses the actual
original plane border and its input neighbour, not an untrimmed carrier.
On the five-edge owner selection the geometry gives 25-17=8 mm; none of
these dimensions or recipe roles enter the algorithm.

### Which construction rounds an edge

The typed application evaluator tries the exact fillet builders first. A body
they finish is the exact result, charged as it always was, with its
continuation checkpoint, so later features read it. A refusal of theirs is no
verdict on the fillet: `catallax_numbers::refund_on_refusal` returns what the
attempt charged and clears an exhaustion it latched, and the rolling-ball
construction starts with the allowance the attempt found. An earlier revision
sent every fillet to the rolling ball and never tried the exact builders. It
turned box fillets into fitted splines and refused every operation after a
fillet, because a numerical body has no continuation.

The rolling ball has three constructions, tried in this order:

1. One edge between two retained spline faces keeps the fitted single-edge
   construction of the first slice.
2. Supports that are planes, or a plane and a cylinder square or parallel to
   it, roll native strips (`fillet_native.rs`, below).
3. Anything else walks the general contact strip of this document.

The native R4 and R8 contact limits are proved from the supports before any of
the three builds a strip, so a radius no ball can roll at refuses in under a
million work units with the face it reaches and the largest radius that fits.

### Native strips between planes and cylinders

A ball against two planes sweeps a cylinder. Against a plane and a cylinder
whose axis is square to the plane it sweeps a torus, and where the axis is
parallel to the plane a cylinder again. These are the elementary cases
[OCCT's `ChFiKPart`](https://github.com/Open-Cascade-SAS/OCCT/tree/V7_9_0/src/ChFiKPart)
builds in closed form before it resorts to walking. Taking the carrier as it
is costs the same at every radius: the owner's rim strip that needed 324
walked stations and 42M of fitting at R3.99 is one torus.

Strips are trimmed where they meet, as the reference kernels trim after
walking. Every vertex in scope has three edges.

- **Tangent continuation.** Two selected edges whose own walls are tangent at
  the vertex (the flat top running into the rounded side) share one
  cross-section of the ball. The unselected seam between the walls is cut back
  to that section's contact.
- **Sharp corner with a concave third edge.** Two convex strips with a common
  support run on past the vertex and meet in the curve their surfaces share.
  It starts where both centre paths cross, the ball touching all three
  supports. There both strips are tangent to the common support, so the two
  surface gradients agree and say nothing of the curve's direction. It leaves
  along a direction in which the two strips bend alike: each bends as the ball
  does across its section and as the common support does along its rail. Of
  the two such lines the curve takes the one that runs into both strips. The
  rest of the curve is the first exit of one strip's section from the other's
  canal surface, a bracketed scalar root in the distance to the other centre
  path, never a Newton step on the two nearly parallel gradients.
- **Unequal reach.** Where the curve leaves one strip through that strip's
  own rail before it has crossed the other, the other strip runs on to the
  first strip's wall. That wall grows by the piece between the two, bounded by
  the wall's cut across the second strip. This is the
  [SMLib `AdjustConcaveCornerTopology`](https://github.com/NVIDIA-Omniverse/usd-brep/blob/c9c5979f04cea8c38213004d4f04df0dcfb40d05/source/SMLib/src/SmFilletCorner.cpp)
  case; no third strip or sphere is invented. The owner's two corners on the
  rounded side are of this kind, its two on the flat top are of equal reach.
- **Open end.** A strip that ends alone is cut by the face its two border
  edges share. The borders are shortened or run on along their own curves to
  the rail ends, as in OCCT's `PerformOneCorner`.

A support that is neither a plane nor such a cylinder, a vertex of another
valence, a convex third edge under two strips, or three strips at one vertex
leave the selection to the general construction.

Every position is a binary64 proposal. Shared curves that are not circles are
C2 quintics through sampled points, first and second derivatives, with spans
halved until the quintic is within 2e-10 mm of the sampled curve at the middle
of every span. The body is assembled beside the exact model without exact
incidence facts, embedded by `numerical_embedding` and admitted by
`build_with_native_trims` at 1e-7 mm; the construction is not authority for
any of it. OpenCASCADE's own fillet of the owner's part agrees with the
admitted bodies to 8e-4 mm3 of 318,000 in volume and 1e-4 mm2 in area at
R0.5, R1, R2 and R3.99, and reads their STEP files raw as one valid closed
solid of 19 faces.

The general walk's first proposal extrapolates both native UV seeds to the
actual first guide parameter, rather than passing the source-edge start to
an already extended section. Spline-support parameter padding is twice the
guide padding, leaving correction room at that section; equal limits put
the Newton solve on the edge of its retained cover. These are numerical
initialisation/extension estimates, not physical radius or completeness bounds.
Accepted stations retain their interval offset and guide-plane checks, and
the existing extension rechecks support regularity and G1 agreement.

This follows [OCCT's boundary-seeded walking](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/BRepBlend/BRepBlend_Walking.cxx#L214-L275),
which sets the first guide parameter before solving the surface coordinates.
[SMLib's corner extension](https://github.com/NVIDIA-Omniverse/usd-brep/blob/3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5/source/SMLib/src/SmFilletCorner.cpp#L658-L667)
also uses an extension estimate. Its implementation multiplies the radius
by 7.5, despite the earlier comment saying 10. Catallax does not copy that
factor or relax regularity when an extension cannot be admitted.
