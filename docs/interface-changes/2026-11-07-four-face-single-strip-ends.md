## 2026-11-07: a single planar strip ends at a four-face vertex

Pieces: features, recipe CLI

One selected straight edge of a convex polyhedron can now end at a four-face
vertex. It previously refused with `fillet.unsupported-end`. The selected
strip remains its exact cylinder. Two conic arcs on the unselected end walls
close it, meeting at the unique hit of those walls' shared input edge in the
strip's shorter normal sector. The other input edges shorten to these exact
contacts; connection alone never starts another blend run.

The input skin must have planar faces and straight edges, with every vertex
in every outward halfspace. This proves the convex material sector. Both end
walls must be transverse to the cylinder axis and face the remaining run.
Their conics project one-to-one onto the cylinder's normal circle. Excluding
the other rail selects the correct sector interval for each arc. Every
plane crossing of each complete conic is enumerated. No event may lie in the
retained arc's interior, and an interior point fixes its material side.
Thus the proof covers each complete finite arc, not sampled endpoints alone.
The retained contact rails must remain positively ordered between their ends.

The cylinder and its cuts are exact. No fitted model face or approximation
allowance is introduced. Each shared boundary is one edge, strip naming
stays `blend.edge`, and the surviving input walls keep their origins. A
single strip end creates no finite-area corner face. There is no public API,
recipe schema or overflow-policy change.

`vertex-4-one` exports a raw-valid six-face solid at an explicit 20M allowance.
It uses 4,167,270 work units, including 2,658,585 for the final fillet. At the
ordinary 1M scoreboard allowance it still exhausts construction work. The
snapshot therefore moves from `fillet.unsupported-end` to
`numbers.operation-budget`, not to a checked export at 1M. Work optimization
is deliberately deferred.

Seven material witnesses and an independent clipped-section volume integral
check the actual CLI export. Its reported file bound is
4.564144149380011770722525218291367494769e-9 mm. Original, translated,
radius-edited and shallower-pyramid controls check closed shared topology and
material membership.

## 2026-11-07: two selected strips meet at a stationary four-face corner

Two equal-radius strips sharing one planar sheet can now meet at a convex
four-face vertex. The selected sheets alone generate the spherical corner.
Its exposed normal-cone side closes with an exact planar cliff; the fourth
input sheet remains an end wall. Each selected cylinder has a retained
lower component as well as its upper run. The native solid has four cylinder
faces, six plane faces and one spherical face. Both finite-area corner
faces retain the original vertex's `corner.vertex` origin.

The constructor proves that the cliff stays through the original vertex
throughout the radius homotopy. Moving cliffs do not use this reduction.
Complete conic/plane event isolation proves every retained arc. The common
centre's already-proved cliff-plane equality is expressed directly in its
coordinates, so STEP does not refine an algebraic exact zero indefinitely.
Replacement face gaps now admit an oriented multi-edge path without adding
an allocation for existing single-edge gaps.

`vertex-4-two` checks its actual raw STEP against the independent section
and cone integral and four target-only material witnesses. It exports a
valid eleven-face solid of volume 1332.8688454588148 mm3, against the analytic
target 1332.8688454588128 mm3. The whole-file bound is
3.515278019967841804683348889091929316915e-9 mm. Total work is 20,079,415:
12,611,809 for the fillet and 7,430,755 for STEP, with 36,851 preceding work.
The ordinary 1M snapshot now records construction budget exhaustion rather
than `fillet.unsupported-corner`; only the private proof allowance is raised.

OCCT still builds. Its smooth closure has a different target set and volume,
recorded as `definition-difference`, not agreement or a waiver. Discriminant
witnesses apply only to the imported Catallax STEP. Shared checks and
tolerances are unchanged. The full derivation is in
[`fillet-coverage.md`](../research/fillet-coverage.md#partial-selected-normal-cliff-corners).

Curved or nonconvex end fans are not covered by this path; unequal earlier
mitres are unchanged. There is no recipe schema or overflow-policy change.
The independent field evaluator is unchanged.

## 2026-11-07: explicit partial-corner closure policy

The fillet data model adds optional `corner: "cliff" | "smooth"` in recipe
versions 4 and 5. Version 3 rejects an explicit policy. Absence preserves
existing behavior; authoring `"smooth"` must never execute a cliff fallback.
`FilletCorner` is exported and every direct Rust constructor is migrated.

The smooth definition preserves the full ball radius and continues an
exposed normal-cone side's two-support contact sheet. On the planar
four-face pyramid this adds the exact A/C cylindrical junction sheet,
tangent to the common sphere and A/C, and clipped by D. It is not a
radius-tapering cone or a fit of OCCT's corner spline. Both the new sheet
and the sphere use the original `corner.vertex` origin.

This data-model change is separate from construction. B-rep and field
dispatch explicitly refuse a requested smooth closure when they lack its
certified continuation or field region. Default and explicit cliff recipes
retain their existing paths. The independent full-radius target and its
section integral are recorded in the coverage research note.

## 2026-11-07: exact full-radius smooth junction construction

The B-rep builds the convex planar two-strip four-face fan with explicit
`corner: "smooth"`. The missing A/C sheet is a cylinder of the requested
radius, joined to the common A/B/C sphere by its great circle and clipped
transversely by D. All spherical contacts are shared topological vertices.
Both the sphere and the continued sheet keep the original four-face
`corner.vertex` origin; no input A/C edge is selected.

The ball is inside every input halfspace. Complete conic events certify
the D cap's angular and input-face bounds. Transverse axis projection and
affine rulings extend those bounds over the whole cylinder face. Native
regressions prove constant radii and agreeing normals along all three
sphere/cylinder great circles, including mirrored, translated, smaller-radius
and different-height controls.

The original smooth variant imports without healing as one valid
nine-face solid. Its volume and target-only discriminants match the
independent section integral in the coverage note. The whole-face STEP
rounding bound is below 100 nm. Default and explicit cliff paths are
unchanged; the separate recipe/schema contract landed before construction.

The independent field compiles the same angular cells with exact cylinder
and sphere interiors. It reuses certified outward convex-prism
halfspaces, not the carriers' natural normals, and no longer requires a
rational profile frame for this junction. Point-query regressions cover
the original, mirrored/translated and moving-centre fans.

The native continuation's separation bound includes the full affected
normal-cone cell, enclosed by the ball whose diameter runs from the common
centre to the original vertex. Other strips and generated corners must be
proved clear of that region. Unproved smooth terminations and mixed
higher-contact field fans report capability refusals; neither evaluator
silently executes the cliff policy.


