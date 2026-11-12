# Cuts after a fitted fillet (M8)

PLAN milestone M8 asks for three things on the P5 housing with its fillets built:

1. The later cuts build, or refuse with a reason that names the missing intersection. Each fitted fillet face is measured again after the cut and is still inside its recorded tolerance.
2. A written list of the intersection kinds these cuts need. M25 scopes the SMLib intersection port from it, and M9 ports in its order.
3. A written verdict on whether Catallax's own fitters cover the blends met so far, or general surface fitting is needed.

This page is parts 2 and 3, plus the evidence behind them. Part 1's outcome is under [Where each cut stops now](#where-each-cut-stops-now).

## The part and the cuts

Canonical P5 (`parts/P5-flanged-housing.json`) finishes with a radius-2 junction fillet between the flange top (plane y=5) and the body wall (cylinder r=12 about y). The countersink mouths cut into the fillet zone, so the fillet is four exact torus pieces between the holes and four fitted rolling-envelope faces over the mouths. The result has 22 faces: 3 planes, 7 cylinders, 4 cones, 4 tori and 4 fitted faces. Each fitted face is a rational B-spline carrier of degree 8 by 2 with 257 by 3 controls, fitted to the rolling-ball definition. Its checked two-way bound is its recorded tolerance: 0.00000417 mm (4.2 nm) for each of the four, inside the 100 nm fit target.

Each case below is the canonical recipe plus one tool body and one more subtraction from the filleted result. Lengths are millimetres. The holes sit on radius 17 at 0, 90, 180 and 270 degrees about y, with the first at +x.

| Case | Tool | Placement | Faces of the part whose meshes cross the tool's |
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
| k12 disjoint tool | sphere r=3 | centre (100, 0, 0) | none |
| k13 fillet first | the canonical countersunk pattern | P5 reordered: fillet the plain housing, then subtract the pattern | see below |

The last column comes from the mesh probe described under [Method](#method). It is evidence for this list only and plays no part in any kernel decision. A crossing between two meshes is strong evidence that the faces meet. No crossing is not evidence that they are apart.

## Where each cut stopped before M8

At `885db576`, with the 20M application allowance, every later cut refused:

| Case | Refusal | Named pair | Mesh crossings for the named kind in this case |
| --- | --- | --- | --- |
| k01 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | none observed |
| k02 | `intersect.traced-pair-required` | oblique or offset axial torus section | none observed |
| k03 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | none observed |
| k04 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | none observed; crossings are cylinder/torus |
| k05 | `intersect.traced-pair-required` | oblique or offset axial torus section | none observed; crossings include spline/plane |
| k06 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | none observed; crossings are spline/cylinder |
| k07 | `intersect.traced-pair-required` | non-coaxial cone/cylinder | observed |
| k08 | `intersect.traced-pair-required` | non-coaxial cone/sphere | observed |
| k09 | `intersect.unsupported-pair` | none, "surface pair is outside piece 4 coverage" | not probed; the follow-up names spline/torus |
| k10 | `intersect.unsupported-pair` | none | crossings include spline/plane |
| k11 | `intersect.traced-pair-required` | oblique or offset axial torus section | none observed; crossings include spline/plane |
| k12 | `intersect.unsupported-pair` | none | none observed |
| k13 | `numbers.operation-budget` at step `drilled` | none | not reached |

Seven of the twelve later cuts refused on a pair whose meshes show no crossing. That does not prove those faces apart, but it means the refusal did not name any kind the case showed crossing evidence for. Three more named no pair at all, including k12, whose tool is 100 mm from the part. Only k07 and k08 named a kind with an observed crossing, and both also have spline crossings. The work used was between 1,081,488 and 1,088,986 units for every later cut, so each refused early, inside the section pass.

The section pass takes each face pair and first asks for the meeting of the two whole supports. Only when that query refuses as needing tracing does it try to prove the two trimmed faces apart with a finite box. Those boxes covered cylinder and cone faces bounded by two complete coaxial circles, spheres, whole tori and Bezier extrusion bands. After the fillet, the countersink cones are no longer bounded by two complete circles, the torus box is a sphere of radius 16 around the fillet, and a fitted face has no box. A refused spline pair was never offered the box test at all. So the first refusal reported was usually a pair with no observed crossing, often far from the cut.

## Where each cut stops now

The boolean now tries exact finite boxes on every face pair or curve whose exact query refused, before it lets the refusal stand. The boxes cover trimmed cylinder and cone faces through their boundary edges, circular arcs shorter than a half turn, torus pieces that keep to one side of the axis, whole tori and spheres, and fitted spline carriers through their control nets. A pair the boxes separate is certified apart. A refused spline query now names its two kinds. Two solids whose lumps the boxes separate answer as disjoint once the staged boolean has refused. Imprinting no longer asks for a spline coincidence proof between two edges with different ends. All of this runs only when one of the two operands has a fitted rolling band. Every other boolean takes the same path and spends the same work as before, including a refused one whose caller catches the refusal and tries another route.

A refused curve query also tests the curve's necessary equations over the bounded support's box. One strict sign proves the curve absent, tangencies included, even when the boxes overlap. For point incidence, an envelope's refused side query can use the existing finite-carrier miss proof. This establishes non-incidence; it leaves the global side undefined. These two exclusions take k03 through the intersection stage. Its next refusal names the outgoing germ of a spline edge at imprinting.

Fresh verification-driver runs of k01 to k12 use the unchanged 20M application allowance. The answers, work, operation records and STEP bytes agree with 1, 2 and 4 kernel threads and the same number of stream workers, and with all 12 available kernel threads and one stream worker. All 477 frozen replay answers are unchanged.

| Case | Outcome | Work | Named kind | Crossing observed for that kind |
| --- | --- | --- | --- | --- |
| k01 | `numbers.operation-budget` | 19,999,984 | none | not reached |
| k02 | `numbers.operation-budget` | 19,994,736 | none | not reached |
| k03 | `boolean.parametric-query` | 14,628,693 | outgoing germ of a spline edge | not a surface-pair query |
| k04 | `intersect.traced-pair-required` | 3,282,020 | skew cylinder/torus | observed |
| k05 | `intersect.unsupported-pair` | 3,423,665 | fitted rolling-envelope spline with plane | observed |
| k06 | `intersect.unsupported-pair` | 3,424,526 | fitted rolling-envelope spline with cylinder | observed |
| k07 | `intersect.traced-pair-required` | 3,926,275 | non-coaxial cone/cylinder | observed |
| k08 | `intersect.traced-pair-required` | 2,827,514 | non-coaxial cone/sphere | observed |
| k09 | `intersect.unsupported-pair` | 2,829,876 | fitted rolling-envelope spline with torus | not probed |
| k10 | `intersect.unsupported-pair` | 2,824,607 | fitted rolling-envelope spline with plane | observed |
| k11 | `intersect.unsupported-pair` | 3,418,207 | fitted rolling-envelope spline with plane | observed |
| k12 | answered, certified manifest/4 | 12,915,677 | none | none observed |
| k13 (earlier run) | `numbers.operation-budget` at `drilled` | 19,997,104 | none | not reached |

k13 is outside the current twelve-case CLI matrix. Its refusal is also covered by the later-cut integration test; its work count above is from the earlier run.

- **k12 builds.** The tool is 100 mm away. The staged boolean still refuses on the fillet's spline edges, then the boxes certify the two solids apart, and the cut answers with the filleted body's own lumps. All four fitted faces keep their carriers and their four chart rows. Their fit bound, checked again with the fit's own check after the cut, is 0.00000417 mm for each, equal to the recorded one. The answer carries the junction certificate, STEP and the same volume enclosure as canonical P5. Its certified manifest, answer digest and STEP bytes are unchanged from the baseline before this step.
- **Seven cuts name a kind with an observed crossing.** k04 to k08, k10 and k11 stop on a pair with crossing evidence. k09 names fitted spline/torus, but its corrected groove has not had the mesh probe repeated.
- **k03 reaches imprinting.** The curve equations exclude the old circle query. Finite-carrier point exclusions then complete the intersection graph. Imprinting refuses because it cannot construct the outgoing germ of a spline edge at a fitted-fillet boundary. The missing query is the spline-edge germ listed below, rather than the previously reported circle with a fitted carrier.
- **k01, k02 and k13 run out of work.** k01 and k02 no longer refuse on a pair before the allowance ends, and the exact analytic path spends it: k01 mostly in traced-curve chart enclosures, k02 and k13 in exact sign decisions of the analytic intersector. Past that point the staged boolean meets the spline-edge queries below, as k12 shows.

### Retained carrier at M10's input checks

M10's field pair engine takes rational coefficients or its exact-offset input. P5's fitted carrier has 844 algebraic expressions among its Cartesian controls and weights. The analytic section routes take checked `RationalPatch` or `RationalAtlas` inputs, so both existing constructors were tried on that unchanged carrier.

Each attempt fits one canonical P5 band and checks the resulting input under one 20M allowance. Both stop before admitting the input, so neither reaches a section search. The [opt-in input test](../crates/catallax-shapes/tests/m8_carrier_input.rs) records the same refusal and work at 1, 2, 4 and the default thread count.

| Input route | Outcome | Work |
| --- | --- | --- |
| Whole checked patch | `numbers.operation-budget` in `RationalPatch::new` | 19,999,970 |
| Joined checked charts | `numbers.operation-budget` in `RationalAtlas::partition` | 19,999,995 |

The exact-offset input proves a support plus a signed normal offset. It supplies no input proof for this retained rolling-envelope fit. The envelope's line solver already carries the transition tangent and end weight as auxiliary variables fixed by their squared values. A section route needs those exact constants too. The test preserves the recorded carrier and tolerance. Raw evidence retained privately.

## The intersection list

Every refused face pair falls in one of three classes:

- **Crossing observed.** The two face meshes cross. Each mesh is within 0.005 mm of its face, so the faces come within 0.01 mm of each other. This is evidence of a meeting, not a proof, and it puts the kind on the list.
- **No crossing observed.** The meshes do not cross. The faces may still touch or cross within the approximation error, so the status is unknown. These kinds stay on the list, after the observed ones, until an exact exclusion removes them.
- **Certified apart.** An exact proof says the two faces do not meet: the kernel's own `Disjoint` answer for the supports, or an exact box separation of the trimmed faces. Only this class leaves M25's scope.

The tables below give the first two classes. Each row says what Catallax answers today and who builds it.

### Spline faces (M25 scope, M9 port, M10 proof)

The fitted face in every row is a rolling-envelope carrier: rational, degree 8 along the ball path and 2 across it, positive weights, 257 by 3 controls, a few millimetres across. The analytic face is always a recipe stock or sketch surface with exact rational definition.

All rows but the last two are crossing observed.

| Kind | Crossing observed in | Today |
| --- | --- | --- |
| fitted spline with plane, plane crossing the patch obliquely | k05 | `intersect.unsupported-pair` |
| fitted spline with plane, the patch's own symmetry plane | k10 | `intersect.unsupported-pair` |
| fitted spline with plane, parallel to that symmetry plane | k11 | `intersect.unsupported-pair` |
| fitted spline with cylinder, axis crossing the patch | k06 | `intersect.unsupported-pair` |
| fitted spline with cone, axis crossing the patch | k07 | `intersect.unsupported-pair` |
| fitted spline with sphere | k08 | `intersect.unsupported-pair` |
| fitted spline with torus, coaxial with the fillet's own rolling axis | k09, not probed: the groove was corrected after the probe, and the kernel names this pair | `intersect.unsupported-pair` |
| curve with fitted spline: circles and ellipses of the cut tool, and traced section curves of the tool with the body | every case with a spline crossing, by how the boolean works (below) | `intersect.parametric-surface`; only lines have a solver, the one the inside-or-outside test uses |

Two notes for the port order:

- Curve with surface comes first. Every row above also needs it: the boolean meets each section curve with every support of both bodies, so a tool edge or a section curve that passes through a fitted patch asks for curve-with-spline before the face pair is split.
- The k10 row is a symmetric configuration: the plane is the mirror plane of the patch, so the section is a single planar curve through the middle of the patch. A port that gets the general oblique plane right but fails here is the tangent-and-singular risk M9 names.

### Analytic pairs that refuse today

These are not spline intersections. Catallax refuses them with `intersect.traced-pair-required`, the same way it refuses other analytic pairs whose section is a curve of degree above two that it does not yet trace. Non-coaxial cylinder pairs already trace with Catallax's own branch proof.

Crossing observed:

| Kind | Crossing observed in | Configuration in the case |
| --- | --- | --- |
| cylinder with torus, skew | k04 | cylinder axis meets the torus axis at a right angle, off the torus centre |
| cone with cylinder, not coaxial | k07 | perpendicular axes, the tool cone crossing the body wall |
| cone with sphere, not coaxial | k08 | sphere centre off the cone axis, crossing the cone near its rim |

No crossing observed in any case, status unknown. Each refused at least once, on faces whose meshes did not cross:

| Kind | Refused in |
| --- | --- |
| plane with torus, oblique or offset from the axis | k01, k02, k05, k06, k07, k11 |
| cone with torus, skew | k07 |
| cone with cone, not coaxial | k07 |
| sphere with torus, not coaxial | k08 |

These could come either from Catallax's own tracing, extended pair by pair the way the cylinder pair works, or from the same ported curve tracer as the spline rows, with Catallax's proof deciding. M25 should state which.

### Spline edges inside the boolean

A fitted face brings spline edges into the part: its rail on the body wall is a rational B-spline curve, and its corner rows meet the torus pieces tangentially, because the fillet is G1 across them. The boolean meets those edges on every face it rebuilds, even when the cut is nowhere near the fillet. With the follow-up change, an exact box proof certifies many refused pairs apart and the cuts get further. Tracing them there found three queries on spline edges, each refused today:

| Query | Where the boolean asks it | Today |
| --- | --- | --- |
| order of the edges leaving a vertex, when two of them are tangent there | imprinting, `germ::outgoing`, at the G1 corners where a rail meets a torus row | `boolean.parametric-query`, "curve equation unavailable" |
| whether a spline edge and another curve are the same curve | imprinting's edge reuse, sewing's edge merge | `intersect.parametric-coincidence` |
| where a section curve meets a fitted carrier, for curves that are not lines | the section pass, each section curve against every support | `intersect.parametric-surface` |

The first is the one that stops the cuts with no observed crossing of a fitted face: k01 to k03, and k12, whose tool is 100 mm away. It needs the normal acceleration of a rational B-spline curve at its end, or an implicit form of the curve near the vertex. That is a curve-with-curve germ, not a surface intersection, but it belongs to the same port: SMLib's curve evaluation answers it, and Catallax's proof decides. The second comes up only where two edges share both ends.

### Not intersections, but needed by the same cuts

A cut that crosses a fitted face leaves a trimmed fitted face. Two things then refuse that are not intersections:

- Inside or outside on a trimmed fitted face. The point test on a rolling-envelope face checks that the face still has its four chart rows (`envelope_retains_chart`) and answers nothing otherwise.
- The fitted-body certificate. The junction certificate is written for the body the fillet made, through its recorded unfilleted base. A later boolean makes a new body with no recorded base, so the certificate and the volume enclosure (`catallax-measure`, junction path) have nothing to start from.

Both belong with M10's builds and M12's tolerance rules.

### Work allowance

`numbers.operation-budget` is a capability refusal under the 20M application contract. k01 and k02 still stop there; the earlier k13 trial did too. Earlier profiling located k01's work in traced-curve chart enclosures, and k02 and k13's work in exact sign decisions of the analytic intersector. They need proofs inside 20M before they can build. M26's wall-clock qualification is separate.

## Fitter verdict

Catallax's own fitters cover every blend met so far, and general surface fitting is not needed for M8.

- The blends on P5 are the rolling-envelope bands over the countersink mouths. `RollingEnvelope::fit` makes them inside the 100 nm target; on P5 each one's bound is 0.00000417 mm. The rolling saddle between nonparallel cylinders and the rolling-sphere band are fitted the same way, by Catallax's own code. No P5 cut meets them.
- A later cut never asks for a new blend. It trims the fitted faces or leaves them alone. A trimmed fitted face keeps its carrier, and the carrier's bound holds on any part of it. Nothing in k01 to k13 refits a surface.
- k13 builds the junction fillet before the countersinks. The fillet on the plain housing is one exact torus, with no fitted face at all. Fitting enters only because the fillet is asked to roll over an existing mouth, and the fitter covers that case.

The open question moves to M11. Each fitted blend family so far, the envelope, the saddle and the sphere band, has its own Catallax fitter, written for one contact law. Blends between two spline faces, on lofts and imported parts, have no fixed contact law to write a fitter for. M11 should record whether they still fit with Catallax's own code, or whether that is where a general fitter pays for itself.

## Method

The probe built each case up to the filleted body and the tool body, then took every pair of one face from each. For each pair it recorded the exact answer of `intersect_surfaces` on the two supports, and a binary64 segment-triangle crossing test between the two face meshes, made at a 0.005 mm proven chord bound. A refused pair whose meshes cross is on the list as observed. A refused pair whose meshes do not cross is on the list as unknown: the test can miss a tangency, or a crossing smaller than the approximation error. The probe proves no pair apart. The kernel's exact `Disjoint` answers do, and so would an exact box separation of the trimmed faces.

The recipes are generated by a short script from the canonical P5 recipe. The script, every probe log and the per-case refusals are retained privately; the recipes themselves land with the builds in the follow-up change.
