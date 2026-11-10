# Cuts after a fitted fillet (M8)

PLAN milestone M8 asks for three things on the P5 housing with its fillets built:

1. The later cuts build, or refuse with a reason that names the missing intersection. Each fitted fillet face is measured again after the cut and is still inside its recorded tolerance.
2. A written list of the intersection kinds these cuts need. M25 scopes the SMLib intersection port from it, and M9 ports in its order.
3. A written verdict on whether Catallax's own fitters cover the blends met so far, or general surface fitting is needed.

This page is parts 2 and 3, plus the evidence behind them. Part 1 lands in a follow-up change from the same lane, and this page will record its outcome.

## The part and the cuts

Canonical P5 (`parts/P5-flanged-housing.json`) finishes with a radius-2 junction fillet between the flange top (plane y=5) and the body wall (cylinder r=12 about y). The countersink mouths cut into the fillet zone, so the fillet is four exact torus pieces between the holes and four fitted rolling-envelope faces over the mouths. The result has 22 faces: 3 planes, 7 cylinders, 4 cones, 4 tori and 4 fitted faces. Each fitted face is a rational B-spline carrier of degree 8 by 2 with 257 by 3 controls, fitted to the rolling-ball definition. Its checked two-way bound is its recorded tolerance: 0.00000417 mm (4.2 nm) for each of the four, inside the 100 nm fit target.

Each case below is the canonical recipe plus one tool body and one more subtraction from the filleted result. Lengths are millimetres. The holes sit on radius 17 at 0, 90, 180 and 270 degrees about y, with the first at +x.

| Case | Tool | Placement | Faces of the part the tool really meets |
| --- | --- | --- | --- |
| k01 wall hole | cylinder r=2 | along x at y=20, through both walls | body wall and bore (cylinders) |
| k02 top slot | block 30 x 6 x 4 | across the body top, y from 27 | body wall, bore, body top |
| k03 rim notch | cylinder r=3 | along y at (16, 16) in x and z, between two holes | flange rim, flange top, flange bottom |
| k04 torus pocket | cylinder r=1/2 | radial at y=6, 36.87 degrees from the first hole | one torus fillet face |
| k05 plane notch | block 4 x 2.5 x 1 | x 11..15, y 5.5..8, z -0.5..0.5 | one fitted face (three block planes), body wall |
| k06 cylinder drill | cylinder r=1/2 | along x at y=6, z=0 | one fitted face |
| k07 cone drill | cone r 1 to 1/2 | along x at y=6, z=0 | one fitted face, body wall |
| k08 sphere dimple | sphere r=1 | centre (13, 6, 0) | one fitted face, one countersink cone |
| k09 ring groove | torus major 12.6, minor 1/2 | coaxial, centred on the fillet surface at y=5.6 | all four fitted faces, all four torus faces |
| k10 axis section | half space z<0 | plane z=0 through the axis | two fitted faces on their symmetry plane, and analytic faces |
| k11 offset section | half space z<1 | plane z=1, parallel to the axis | two fitted faces off their symmetry plane, and analytic faces |
| k12 disjoint tool | sphere r=3 | centre (100, 0, 0) | nothing |
| k13 fillet first | the canonical countersunk pattern | P5 reordered: fillet the plain housing, then subtract the pattern | see below |

The "really meets" column comes from the probe described under [Method](#method). It is evidence for this list only and plays no part in any kernel decision.

## Where each cut stopped before M8

At `885db576`, with the 20M application allowance, every later cut refused:

| Case | Refusal | Named pair | Is the named pair among the faces the cut really meets? |
| --- | --- | --- | --- |
| k01 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | no |
| k02 | `intersect.traced-pair-required` | oblique or offset axial torus section | no |
| k03 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | no |
| k04 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | no; the cut needs cylinder/torus |
| k05 | `intersect.traced-pair-required` | oblique or offset axial torus section | no; the cut needs spline/plane |
| k06 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | no; the cut needs spline/cylinder |
| k07 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | yes |
| k08 | `intersect.traced-pair-required` | non-coaxial cone/sphere | yes |
| k09 | `intersect.unsupported-pair` | none, "surface pair is outside piece 4 coverage" | the cut needs spline/torus |
| k10 | `intersect.unsupported-pair` | none | the cut needs spline/plane |
| k11 | `intersect.traced-pair-required` | oblique or offset axial torus section | no; the cut needs spline/plane |
| k12 | `intersect.unsupported-pair` | none | the cut meets nothing |
| k13 | `numbers.operation-budget` at step `drilled` | none | not reached |

Seven of the twelve later cuts refused on a pair whose two faces never meet. Three more named no pair at all, including k12, whose tool is 100 mm from the part. Only k07 and k08 named a kind they really need, and both need a spline kind as well. The work used was between 1,081,488 and 1,088,986 units for every later cut, so each refused early, inside the section pass.

The section pass takes each face pair and first asks for the meeting of the two whole supports. Only when that query refuses as needing tracing does it try to prove the two trimmed faces apart with a finite box. Those boxes covered cylinder and cone faces bounded by two complete coaxial circles, spheres, whole tori and Bezier extrusion bands. After the fillet, the countersink cones are no longer bounded by two complete circles, the torus box is a sphere of radius 16 around the fillet, and a fitted face has no box. A refused spline pair was never offered the box test at all. So the first refusal reported was usually a pair of faces far from the cut.

## The intersection list

These are the kinds the later cuts need, by what the cut faces really meet. Each row says what Catallax answers today and who builds it.

### Spline faces (M25 scope, M9 port, M10 proof)

The fitted face in every row is a rolling-envelope carrier: rational, degree 8 along the ball path and 2 across it, positive weights, 257 by 3 controls, a few millimetres across. The analytic face is always a recipe stock or sketch surface with exact rational definition.

| Kind | Needed by | Today |
| --- | --- | --- |
| fitted spline with plane, plane crossing the patch obliquely | k05 | `intersect.unsupported-pair` |
| fitted spline with plane, the patch's own symmetry plane | k10 | `intersect.unsupported-pair` |
| fitted spline with plane, parallel to that symmetry plane | k11 | `intersect.unsupported-pair` |
| fitted spline with cylinder, axis crossing the patch | k06 | `intersect.unsupported-pair` |
| fitted spline with cone, axis crossing the patch | k07 | `intersect.unsupported-pair` |
| fitted spline with sphere | k08 | `intersect.unsupported-pair` |
| fitted spline with torus, coaxial with the fillet's own rolling axis | k09 | `intersect.unsupported-pair` |
| curve with fitted spline: circles and ellipses of the cut tool, and traced section curves of the tool with the body | k05 to k11 | `intersect.parametric-surface`; only lines have a solver, the one the inside-or-outside test uses |

Two notes for the port order:

- Curve with surface comes first. Every row above also needs it: the boolean meets each section curve with every support of both bodies, so a tool edge or a section curve that passes through a fitted patch asks for curve-with-spline before the face pair is split.
- The k10 row is a symmetric configuration: the plane is the mirror plane of the patch, so the section is a single planar curve through the middle of the patch. A port that gets the general oblique plane right but fails here is the tangent-and-singular risk M9 names.

### Analytic pairs that refuse today

These are not spline intersections. Catallax refuses them with `intersect.traced-pair-required`, the same way it refuses other analytic pairs whose section is a curve of degree above two that it does not yet trace. Non-coaxial cylinder pairs already trace with Catallax's own branch proof.

| Kind | Needed by | Configuration in the case |
| --- | --- | --- |
| cylinder with torus, skew | k04 | cylinder axis meets the torus axis at a right angle, off the torus centre |
| cone with cylinder, not coaxial | k07 | perpendicular axes, the tool cone crossing the body wall |
| cone with sphere, not coaxial | k08 | sphere centre off the cone axis, crossing the cone near its rim |

They could come either from Catallax's own tracing, extended pair by pair the way the cylinder pair works, or from the same ported curve tracer as the spline rows, with Catallax's proof deciding. M25 should state which.

### Spline edges inside the boolean

A fitted face brings spline edges into the part: its rail on the body wall is a rational B-spline curve, and its corner rows meet the torus pieces tangentially, because the fillet is G1 across them. The boolean meets those edges on every face it rebuilds, even when the cut is nowhere near the fillet. Tracing the cuts past the first refusal found three queries on spline edges, each refused today:

| Query | Where the boolean asks it | Today |
| --- | --- | --- |
| order of the edges leaving a vertex, when two of them are tangent there | imprinting, `germ::outgoing`, at the G1 corners where a rail meets a torus row | `boolean.parametric-query`, "curve equation unavailable" |
| whether a spline edge and another curve are the same curve | imprinting's edge reuse, sewing's edge merge | `intersect.parametric-coincidence` |
| where a section curve meets a fitted carrier, for curves that are not lines | the section pass, each section curve against every support | `intersect.parametric-surface` |

The first is the one that stops the cuts that never touch a fitted face. It needs the normal acceleration of a rational B-spline curve at its end, or an implicit form of the curve near the vertex. That is a curve-with-curve germ, not a surface intersection, but it belongs to the same port: SMLib's curve evaluation answers it, and Catallax's proof decides. The second comes up only where two edges share both ends.

### Not intersections, but needed by the same cuts

A cut that crosses a fitted face leaves a trimmed fitted face. Two things then refuse that are not intersections:

- Inside or outside on a trimmed fitted face. The point test on a rolling-envelope face checks that the face still has its four chart rows (`envelope_retains_chart`) and answers nothing otherwise.
- The fitted-body certificate. The junction certificate is written for the body the fillet made, through its recorded unfilleted base. A later boolean makes a new body with no recorded base, so the certificate and the volume enclosure (`catallax-measure`, junction path) have nothing to start from.

Both belong with M10's builds and M12's tolerance rules.

### Work, not geometry

Two cases stop on the 20M application allowance instead of a missing kind. k13 spends 18.0M units in exact sign decisions of the analytic intersector while subtracting the countersinks from the exact torus fillet. The same signature costs k02 18.9M units once its false refusals are gone. Both are M26 work, after the capability is in.

## Fitter verdict

Catallax's own fitters cover every blend met so far, and general surface fitting is not needed for M8.

- The blends on P5 are the rolling-envelope bands over the countersink mouths. `RollingEnvelope::fit` makes them inside the 100 nm target; on P5 each one's bound is 0.00000417 mm. The rolling saddle between nonparallel cylinders and the rolling-sphere band are fitted the same way, by Catallax's own code. No P5 cut meets them.
- A later cut never asks for a new blend. It trims the fitted faces or leaves them alone. A trimmed fitted face keeps its carrier, and the carrier's bound holds on any part of it. Nothing in k01 to k13 refits a surface.
- k13 builds the junction fillet before the countersinks. The fillet on the plain housing is one exact torus, with no fitted face at all. Fitting enters only because the fillet is asked to roll over an existing mouth, and the fitter covers that case.

The open question moves to M11. Each fitted blend family so far, the envelope, the saddle and the sphere band, has its own Catallax fitter, written for one contact law. Blends between two spline faces, on lofts and imported parts, have no fixed contact law to write a fitter for. M11 should record whether they still fit with Catallax's own code, or whether that is where a general fitter pays for itself.

## Method

The probe built each case up to the filleted body and the tool body, then took every pair of one face from each. For each pair it recorded the exact answer of `intersect_surfaces` on the two supports, and a binary64 test of whether the two faces meet: display meshes at a 0.005 mm proven chord bound, then a segment-triangle crossing test between the two face meshes. A pair that the mesh test calls meeting and the kernel refuses is a real need. A pair that the mesh test calls apart is not, whatever the kernel says about the whole supports.

The recipes are generated by a short script from the canonical P5 recipe. The script, every probe log and the per-case refusals are retained privately; the recipes themselves land with the builds in the follow-up change.
