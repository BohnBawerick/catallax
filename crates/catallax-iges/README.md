# Exact-surface IGES writer

`write_iges(&Model, &[BodyId], &IgesOptions, &mut Decider, &mut Write)`
prepares a deterministic IGES 5.3 manifold-solid B-rep, proves its output bound,
and only then writes it. Kernel coordinates, file units and report distances are
millimetres. There is no importer, application command or change to STEP.

The inspected public baseline was main commit
`4b7ce8f0632feff5e327fbd1aeda164657dc640a`. It exposed retained analytic
surface/curve definitions and shared topology, but no IGES package or writer.
Its STEP tolerance is a mark threshold; this API instead refuses an output whose
proven bound exceeds the requested inclusive tolerance. PLAN section 12 item 22
agrees with this separation. No existing crate API needed modification.

## Supported cases and explicit limits

| Input | Result |
| --- | --- |
| Valid, finite, unmarked closed bodies; multiple lumps and shells | Manifold solids, one per lump |
| Planar faces with finite line/circular-arc wires, including holes | Exact type 190 supports and shared trims |
| Cylindrical/conical bands between two full vertexless circular edges | Exact type 192/194 supports, one added seam per band |
| Finite axial cylinder/cone sectors with two circular arcs and two generators | One unwrapped four-edge wire; cardinal endpoints in the retained surface frame, including 90/180/270-degree spans and seam crossings |
| Box, bored plate, public P2 collar | Independently tested with 6, 7, 6 faces |
| Rational placements and dimensions whose circle frames have terminating decimal unit coordinates | Exact definitions, followed by the certified common scale |
| A cone half-angle with a proved finite decimal degree encoding | Exact type 194 angle; rational slope 1 has angle 45 degrees |
| Other surface/curve kinds | Typed `UnsupportedSurface` / `UnsupportedCurve` |
| General algebraic coordinates, nonterminating circle frames, uncovered cone angles, other periodic trims, apex caps, nonaligned ring seams | Typed `UnrepresentableDefinition`, or a shared decision refusal |
| Marked input without a retained whole-definition certificate | `UncertifiedModel` |

These are coverage limits, not claims that the refused geometry is impossible.
In particular this first writer does **not** cover every body made of the three
surface kinds. It fails closed on uncovered definitions and periodic trims.
Non-millimetre output is not exposed. The caller supplies geometrically valid
closed bodies; `Model::check` additionally checks the public structural contract.
This is not a new geometric body validator.

Input edges remain shared. The format requires a vertex on each vertexless ring;
the report counts those additions. A periodic seam is an exact straight
rational degree-one type 126 curve, shared by its two opposite coedge uses.
No input edge or face is merged, and the source model is immutable. Faces use
190/192/194, original curves use 110/100 (sector generators use bounded degree-one
126 curves), and topology uses 502/504/508/510/514/186.
There is no tessellated or fitted surface fallback.

Finite sectors require one alternating four-edge loop. Exact rational projection
proves coaxial, constant-height circular rails, nonzero generator heights and
cardinal rays. Cone vertices must lie strictly on the positive nappe, away from
the apex. Directed circle normals, edge senses, coedge senses and face senses
select the signed angular span; accumulating that span keeps one continuous
degree chart through the seam. Opposing spans must close one rectangle. No
extra seam, face, edge or vertex is introduced for these sectors. Finite type-100
arcs use an exact orthonormal placement whose start angle is zero, so their
native ranges match bounded pcurves even when a retained circle has a different
reference ray. This reparameterizes the same exact circle without fitting it.

## Whole-domain certificate

The target is the retained exact model boundary, not an imported numeric shape.
All accepted rational recoveries from `Real` expressions require exact equality
from the shared decider; interval midpoints and continued fractions only propose
candidates. Non-rational candidates are never accepted as exact.

Let `D` be the least common multiple of the denominators of **all** metric file
parameters, including rounded conical pcurve coordinates. Select the nearest
finite decimal `q` to `1/D` at `scale_digits` significant digits (lower on a tie),
and set `s = D*q > 0`. Every emitted metric parameter is exactly `s` times its
prepared value. Because each denominator divides `D`, each emitted value is a
terminating decimal. Direction vectors become exact integer homogeneous vectors;
circle placement matrices remain exactly orthonormal. Thus every 3D support,
line, circle, vertex, seam and finite trim is the same positive similarity of
its retained definition. This coordinated rounding preserves incidence and
orientation, unlike unrelated coordinate rounding. `scale_digits` controls the
scale choice; it is not a per-coordinate significant-digit truncation.

The coordinate bound `B` is an L1 bound on the **whole** finite boundary:

- A line segment is bounded by its endpoints.
- Every point on a circle is bounded by `|center|_1 + 3*radius`; this also bounds
  any finite arc.
- A bounded planar region is contained in the convex hull of its boundary.
- Every point of a cylindrical or conical band is a convex combination of the
  corresponding points on its two boundary circles. The endpoint circle bounds
  therefore cover the entire band, including its seam.

Consequently `|s*x - x|_2 <= |s-1|*B` everywhere, in both directions under the
same correspondence. This is not a sampled residual test.

Periodic pcurves use exact degree endpoints 0 and 360. Conical slant levels can
be irrational: outward 160-bit balls bound each exact endpoint, and the maximum
absolute difference from the emitted 25-digit rational bounds its rounding.
Interpolation between endpoints bounds the entire straight pcurve. A unit-speed
cone generator converts that parameter error directly to a positional bound.
The shared seam's 3D and 2D representations have the same [0,1] domain.
Circular pcurves use a bounded [0,T] domain, where `T` is a 17-digit decimal of
the arc's radian span (or `2*pi` for a full ring). An outward ball bounds the
absolute endpoint error; multiplication by the whole-circle coordinate bound conservatively
bounds the whole-domain angular correspondence/compatibility error. These
angular and slant bounds are **added**, not replaced by their maximum.
Linear interpolation of the two rounded slant levels bounds the whole finite
sector; the cone generator is unit speed in slant. Angular interpolation covers
the whole directed arc, including a span longer than 180 degrees. Common scaling
then applies to all metric pcurve levels and all shared 3D geometry together.

Writing `E` for that sum, the report declares the rational bound

```
max_deviation = abs(s - 1) * B + s * E
```

It covers face/trim displacement and the redundant parameter/3D boundary
compatibility error. It does not bound area or volume error, and does not cover
additional arithmetic inside a downstream reader. The formal file certificate
and native import checks are separate evidence. No modelling tolerance is changed.
If this bound exceeds `tolerance`, the writer returns `ToleranceExceeded` and
never calls the sink. The inclusive equality case succeeds. I/O errors follow
`Write::write_all`: a failing sink can have accepted a prefix.

## Reproduction

Run ordinary contracts with `cargo test --release -p catallax-iges`.
The explicitly ignored test only generates files for the separate native checker:

```sh
CATALLAX_IGES_FIXTURES=/tmp/iges-fixtures cargo test --release -p catallax-iges \
  --test contract write_readback_fixtures -- --ignored
python3 crates/catallax-iges/tests/readback.py /tmp/iges-fixtures
```

Install `cadquery-ocp==8.0.1.1.0` in a scratch environment for the second command.
The checker requires the bindings and fails if they are missing. It replaces
IGESControl's healing boundary tool with the base `IGESToBRep` tool, transfers
only the type 186 root directly, and never invokes the control actor,
`XSAlgo::ProcessShape`, `ShapeFix`, sewing or application repair. It checks solid
validity, face and exact surface-kind counts, shared incidence, seam counts,
units and signed-volume diagnostics. This also exercises a rational pose and a
nonzero rounded-scale case. Ordinary workspace CI does not install native
bindings for Rust test shards; this separate check must be run explicitly.
The Rust tests check strict rejection and unchanged sinks, rational emitted
vertices, composed trim bounds, determinism, unsupported inputs, work refusal,
and unchanged retained STEP bytes before and after IGES export.

The finite-sector regression was reproduced on main
`10d57eb79110941f930c11358f663968440b1b46`: the public sketch/revolve recipe
builds, then `write_iges` returns `UnrepresentableDefinition` for its periodic
face. `partial_axial_sectors_export_through_the_recipe_path` preserves that
failing-before/passing-after boundary. The native fixtures add cylinder/cone
90/180/270-degree sectors, opposite cone axes, recharted seam crossings with
reversed edges and circle reference rays, and a nonzero common-scale cone.
The checker also checks finite domains, matching curve/pcurve ranges and sampled
3D/pcurve compatibility. Those numerical checks do not establish the bound above.
General angle encoding, apex caps and other periodic trim shapes remain refused.

## Primary references

- [Open CASCADE IGES translator documentation](https://occt3d.com/dev/doc/overview/html/occt_user_guides__iges.html): IGES 5.3, B-rep entity mapping, units and shape processing.
- [OCCT IGESSolid tools](https://github.com/Open-Cascade-SAS/OCCT/tree/V7_8_1/src/IGESSolid): authoritative parameter ordering and forms for 186/190/192/194/502/504/508/510/514.
- [OCCT base B-rep transfer](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_8_1/src/IGESToBRep/IGESToBRep_BRepEntity.cxx): shared edge lists, face/shell sense and direct transfer.
- [OCCT base boundary transfer](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_8_1/src/IGESToBRep/IGESToBRep_IGESBoundary.cxx) and [control boundary tool](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_8_1/src/IGESControl/IGESControl_IGESBoundary.cxx): the distinction necessary for a no-healing read-back.
- [NIST IGES recommended practices](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir4600.pdf): IGES record/placement interoperability.
