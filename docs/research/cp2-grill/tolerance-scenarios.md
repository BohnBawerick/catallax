# Tolerance scenarios for Catallax (units mm)

Builds on `tolerance-numbers-research.md` (same folder). Labels:
[read] = primary text read in a session; [snippet] = search-engine summary only, page not read;
[unread/memory] = memory; **Inference** = my own reasoning or arithmetic.

## 0. New evidence gathered this session

| Fact | Label |
|---|---|
| Parasolid size box 1000 units centred at origin, 1e-8 linear resolution; "large transformations" let a body inside the box be placed outside it (factory example, a few km) | [snippet] (Parasolid V35 Model Structure chapter, fetch refused; Altair CADFEKO doc says the same) |
| Ansys DesignModeler (Parasolid): normal 1 km box at 1e-8 m; "large model support" gives a 10 km box with precision loosened to 1e-5 m, "could result in merging points inadvertently". Ratios: 1e11 normal, 1e9 large mode | [read] ansyshelp.ansys.com dmBooleanOpsModel |
| OCCT forum: large structures at 0.5 mm tolerance, Fuse grew tolerances to 1.2; advice "for big objects ... use meter units ... scale down 1000 times"; move parts near the origin before Booleans | [snippet] dev.opencascade.org "Tolerance issues", "opencascade and Boolean operations" |
| OCCT forum: an edge shorter than one of its vertex tolerances makes a Boolean through it fail; fuzzy value must be "significantly smaller than the minimum length of any geometry you want to keep" | [snippet] |
| OCCT GitHub issues (late 2026): fillet with radius equal to an adjacent fillet gives a 4.83 mm vertex tolerance from 1.5e-7 inputs and passes BRepCheck; other fillet cases give 3.9e4 mm and 271 mm vertex tolerances | [snippet] issues #1701, #1695, #1639; not verified |
| SolidWorks/Parasolid practical range about 1e-5 mm to 1 km; users report a minimum feature near 1 um (one report says 100 nm); MEMS users scale models up | [snippet] comp.cad.solidworks, forum.solidworks.com, CADmunity |
| ACIS users quote dynamic range 1e-5 to 1e4 units; "smallest feature will limit the largest" | [snippet] |
| Creo relative accuracy = ratio of smallest model dimension to part size; absolute = relative x part size; default relative value is the same for all sizes; accuracy decides whether geometry in assemblies merges | [read] PTC Creo 11 help |
| Creo practice: absolute accuracy recommended for merging parts of different sizes and in manufacturing mode; relative accuracy makes small details fail on large parts | [snippet] PTC community, creo3d.com, Creo Tool Design guide |
| Convert-to-sheet-metal on STEP fails when thickness varies by more than a few hundredths of a mm (from general solid modelling tolerances) | [snippet] cadshift.com, Autodesk forums |

## 1. Does "1e-7 mm, size limit 1e4 mm" apply to Catallax?

The Parasolid (1e11) and ACIS (1e10) resolution-to-size ratios exist because f64 has about 16 digits and
iterative geometry (Newton, marching, near-tangent intersection) loses several of them. Ansys's large mode
shows the trade directly: 10x the box costs 1000x the precision.

**Inference:**
- Exact rational entities have no such ratio. A plane at 1e9 mm is exact; only bit length (speed) grows.
  The size limit does not apply to them.
- Fitted f64 entities do have it. f64 spacing at 1e5 mm is 1.5e-11 mm, so 1e-7 is still 4 digits above
  rounding, but ill-conditioned steps (near-tangent intersections, which lose about half the digits) need
  more headroom. Parasolid's 1e11 ratio is the best-tested headroom figure.
- The error depends on the coordinate magnitude in the frame where the computation runs, not on part size.
  A 2 mm part placed at 2 m from the origin computes like a 2 m part unless the operation runs in a local
  frame.

So replace the fixed limit with a position-dependent floor for f64 entities:

`eps(R) = max(1e-7, 1e-11 * R)`, R = largest |coordinate| of the entity in the frame the operation runs in.

| R (mm) | eps | fitted floor 10 eps | fitted cap 100 eps (min 1e-5) |
|---|---|---|---|
| <= 1e4 (anything up to 10 m) | 1e-7 | 1e-6 | 1e-5 |
| 3e4 (wing) | 3e-7 | 3e-6 | 3e-5 |
| 1e5 (100 m ship) | 1e-6 | 1e-5 | 1e-4 |
| 1e6 (1 km site) | 1e-5 (= Parasolid) | 1e-4 | 1e-3 |

At the Parasolid box edge this equals Parasolid's 1e-5 mm. Inside 10 m it equals OCCT. Operations run in a
frame centred on the operands' bounding box, reached by an exact rational translation, so R stays small
(a floating origin; Parasolid's "large transformations" is the same idea for assemblies).

## 2. Question 1: how the class cap and the local cap combine

Write the local rule pairwise, since coincidence is "within the sum of tolerances" (Jackson):

For any two **non-adjacent** entities a, b (no bounding relation, no shared vertex):
- error if `t_a + t_b > d(a,b) / 2`
- warn if `t_a + t_b > d(a,b) / 5`

With equal tolerances this is the proposed L/4 error and L/10 warning. Per entity:

`cap_eff(x) = min( C_class(x), L(x) / 4 )`, warning above `L(x) / 10`,

`L(x) = min distance from x to any non-adjacent entity`. This one definition covers:
- short edge: its two end vertices are non-adjacent, so length >= 2(t1 + t2) = 4t (Jackson's 4x);
- narrow or sliver face: distance between its non-adjacent boundary edges;
- thin wall: distance between opposite faces (non-adjacent);
- a feature crowding a distant feature (a hole close to a wall).

**No small feature (the owner's question).** L is then large (a 1000 mm block: L is hundreds of mm), L/4 is
far above any class cap, and `cap_eff = C_class`. The local term only ever tightens; it never loosens the
class cap. It begins to bite only when L < 4 x C_class: L < 0.04 mm for 0.01 imports, L < 0.4 mm for 0.1
opt-in imports, L < 4e-5 mm for own fitted geometry. Tolerances are measured, never set to the cap, so a
loose cap on a featureless part costs nothing.

Pairs that share only a vertex (two edges leaving a corner, a face tapering to a point) are tested only
outside the shared vertex ball plus 2(t_a + t_b). Tangent-meeting edges at a cusp stay an open question.

**Cheap computation (inference):** you never need L itself, only "is anything non-adjacent within
5(t_a + t_b)". Keep a BVH of entity boxes, each inflated by 5 t. Broad phase returns overlapping box
pairs; drop adjacent pairs by topology lookup; run an exact distance only on the rest. With t <= 0.01 the
inflation is <= 0.05 mm, so on a part with no small features the broad phase returns only adjacent pairs
and no distance is computed. Edge length and vertex-vertex distance are free. Update incrementally for the
entities an operation changed. Cost O(n log n) build, near O(k) per operation.

## 3. Reference numbers per size (inference from sourced rules)

| Scenario | size | Catallax own fitted cap | Parasolid | OCCT | SMLib ZoneTol | SMLib gap cap bbox/800 | Creo relative 0.0012 x size |
|---|---|---|---|---|---|---|---|
| MEMS | 0.1 | 1e-5 | 1e-5 | 1e-7 | 1e-6 | 2e-4 | 1.2e-4 |
| watch gear | 2 | 1e-5 | 1e-5 | 1e-7 | 1e-5 | 4e-3 | 2.4e-3 |
| bracket | 60 | 1e-5 | 1e-5 | 1e-7 | 1e-5 | 0.13 | 0.07 |
| block | 1000 | 1e-5 | 1e-5 | 1e-7 | 1e-4 | 2.1 | 1.2 |
| car body | 5000 | 1e-5 | 1e-5 | 1e-7 | 1e-4 | 10.6 | 6 |
| wing | 3e4 | 3e-5 | 1e-5 | 1e-7 | 1e-4 | 64 | 36 |
| ship | 1e5 | 1e-4 | 1e-5 | 1e-7 | 1e-3 | 212 | 120 |

(bbox diagonal taken as 1.7 x size. Spline segment count goes as tol^(-1/4): fitting at 1e-6 instead of
1e-5 costs 1.8x segments, instead of 0.01 costs 10x. [unread/memory, approximation theory])

## 4. Scenarios

Format: what the rules do / goes right / goes wrong / fix. "vs" compares with OCCT, Parasolid, SMLib.

**4.1 Bracket 60 mm, own model.** Planes and cylinders exact; fillets fitted at 1e-6, stored <= 1e-5; L
around 1 mm, so local rule silent. Export uncertainty 1e-5. Right: clean in every receiver, prismatic
Booleans exact with no growth. Wrong: nothing material; fillets 1.8x more segments than Parasolid. vs:
same as Parasolid and SMLib, better than OCCT (no fillet tolerance blow-ups; OCCT writes the average
tolerance to STEP by default). Fix: none.

**4.2 Car body 5 m.** Mostly imported class-A surfaces and 0.7 mm panels. R <= 1e4, so eps 1e-7 and fitted
cap 1e-5. Panels: L = 0.7, import tol up to 0.0175 passes silently, 0.07 warns, 0.175 errors. Right:
f64 headroom fine at 5e3 mm (spacing 9e-13). Wrong: own fitted trims on 5 m edges at 1e-6 carry about
10x the segments of an NX export at 0.0254; large STEP files and slow fitting. vs: same as Parasolid,
better than SMLib (its gap cap would accept 10 mm gaps). Fix: optional export-only refit of fitted
geometry to receiver tol / 10, re-measured before writing (open: needs edges re-fitted too).

**4.3 Aircraft wing 30 m.** R = 3e4 gives eps 3e-7, fitted cap 3e-5. Right: still inside Parasolid's box,
so receivers accept. Wrong: with a fixed 1e-7 floor and 1e-5 cap (the original rule), near-tangent
intersections at 3e4 mm would miss the cap and ops would be refused. vs: same as Parasolid, better than
OCCT (fixed 1e-7 at large coordinates; forum advice is to rescale to metres). Fix: eps(R) (section 1).

**4.4 Ship or building 100 m+.** eps 1e-6, fitted cap 1e-4. Hull blocks modelled in local frames, placed
with exact rational transforms. Right: no hard size limit; georeferenced BIM coordinates (1e9 mm) cost
nothing if operations run in local frames. Wrong: a single body over 1 km cannot go to Parasolid-based
receivers (box +-500 m); a fitted cap of 1e-4 at 100 m is above Parasolid's 1e-5, so some edges arrive
tolerant there. vs: better than OCCT and SMLib (SMLib zone 1e-3 at this size), same as Parasolid with
large transforms. Fix: export warning for bodies over 1 km; ship as assembly instances.

**4.5 Watch gear 2 mm.** Involute flanks fitted at 1e-6; tip radii around 0.01 to 0.02. L/10 = 1e-3,
far above 1e-5. Right: kernel-side fine. Wrong: in a 0.01 receiver (NX, Creo, Rhino) features below
about 0.04 (4 x receiver tol) get merged or break joins; Rhino wants features >= 10 x tol [read earlier].
vs: same as all four kernels; the failure is in the receiver. Fix: export preflight that lists features
below 10 x and 4 x the receiver tolerance.

**4.6 MEMS 0.1 mm, 1 um features, 0.1 um gaps.** Exact geometry 1e-7 gives 1000x margin on a 1e-4 gap.
Fitted cap 1e-5 hits the L/10 warning across a 1e-4 gap. Right: works unscaled. Wrong: fitted entities
near the smallest gaps warn; every mm-based receiver at 0.01 destroys the part. vs: better than Parasolid
(1e-5 resolution; features at 10x precision break its 100x separation guideline; users scale up 1000x),
same as OCCT and SMLib (1e-6 zone here). Fix: for own fitted entities make the local rule a fit target,
not only a check (refit to min(1e-6, L/100), not below 10 eps); offer a scaled export (x1000) with a
clear label, as CADFEKO and MicroStation do for Parasolid.

**4.7 Thin sheet metal 0.5 mm.** Own model exact: thickness equal exactly, so unfold recognition never
sees "almost uniform". Imports: L = 0.5 (opposite faces). 0.01 import passes; 0.05 gaps warn (t > L/10);
0.1 opt-in passes only because 0.11 + 0.11 < 0.25; at 0.3 mm walls 0.1 opt-in is refused. Right: the
rule catches a tolerance zone that spans the wall, which other kernels let through (OCCT healing max
1 mm by default). Wrong: refusals on poor sheet-metal imports where others "work" with a broken solid.
vs: better than all on own data and on honesty; worse on convenience. Fix: offer local refit (rebuild the
gappy edge from one face) before refusing.

**4.8 Big simple block, no small feature.** See section 2: class cap governs, local rule silent. Poor
import with 0.05 gaps on a 1000 mm block: 0.055 exceeds the 0.01 auto cap, so refit, then opt-in. SMLib
accepts up to 2.1 mm automatically, OCCT up to 1 mm. Right: no size-scaled loosening, so exports stay
under receivers' 0.01. Wrong: friction on parts where a big tolerance is locally harmless. vs: worse on
convenience, better on export. Fix: none to the cap (the 0.01 cap exists for receivers, not for local
safety); make the opt-in a one-click per-file choice that names the receivers it will break.

**4.9 Tight fillets near tangency.** A normal error d at curvature radius r shifts the intersection point
along the surface by about sqrt(2 r d): r = 5, d = 1e-6 gives 3e-3 mm (inference, arithmetic). If two
edges meeting at a vertex are fitted separately, their ends disagree by that much and the vertex
tolerance jumps 1000x. Equal-radius adjacent fillets are a documented OCCT failure (4.83 mm vertex,
passes BRepCheck) [snippet]. Right: the own-fitted cap 1e-5 turns that into a refused operation, not a
silent bad solid. Wrong: refusals; slow marching near tangency. vs: better than OCCT (fails loudly),
same as Parasolid (SP-curves, many segments). Fix: vertex-first construction (solve the vertex on all
its surfaces, then fit edges clamped to it); test radius equality exactly (rationals) and route
equal-radius blends to a special case.

**4.10 Long chains (50 Booleans and fillets).** Exact operands give exact results: zero growth.
Fitted: each op re-measures and stores 1.1 x (gap + delta) from measurement. Wrong if the margin were
applied to stored tolerances: 1.1^50 = 117x (inference). SMLib's 2x per re-validation compounds faster.
Deviation from design intent (a separate quantity) grows at most additively: 50 x 1e-6 = 5e-5.
vs: better than OCCT (documented growth after fillets and Fuse) and SMLib; same as Parasolid if the
application calls PK_EDGE_optimise. Fix: state that the margin applies only to measured values.

**4.11 Good STEP import, CATIA (declares 0.005) or NX (declares 0.01).** Per-edge measured. CATIA models at
0.001, so most edges land <= 1e-3 [inference from CATIA 0.001 merging distance]. NX tolerant edges reach
0.01 and tolerant vertices 0.0141 [snippet]. Wrong in the original rules: the 1e-3 "good import" cap
flags nearly every NX model as poor. Right: measuring ignores optimistic declared values (Rhino's known
join failure). vs: better than OCCT (trusts the declared value, heals to 1 mm), same as Parasolid
tolerant modelling. Fix: drop the good/poor class split; import auto cap 0.01 for edges and faces, 0.02
for vertices; report bands (<= 1e-3 tight, <= 0.01 standard); refit only edges above 1e-3.

**4.12 Poor import, 0.05 gaps.** Over the auto cap: refit first (re-intersect adjacent surfaces, extend
if needed). True surface mismatches survive refit; opt-in to 0.1 stores about 0.055. Wrong: on export to a
0.01 receiver those edges open; to CATIA 0.001, more so. vs: same end state as others (OCCT and SMLib
accept silently with larger margins), more friction, smaller margins (1.1x vs 2x). Fix: surface-
modifying heal (Parasolid-style simplify) as the recommended path.

**4.13 Mixed assembly, 2 mm part and 2 m part.** Per-entity tolerances, no global value; ops in local
frames. Merging both into one body: R about 2e3, eps 1e-7. Right: nothing changes for the small part.
vs: same as OCCT and Parasolid; better than Creo relative (0.0012 x 2000 = 2.4 mm, larger than the small
part) and SMLib (gap cap 4 mm on the merged body). Fix: none.

**4.14 Imports with slivers or micro edges below 0.01.** Validity is tested against the entity's own
tolerance (R5), not against 0.01: a 0.005 edge with 0.001 vertices is valid. Invalid ones are collapsed
(merge vertices, replace sliver by edge) when the result passes the caps and the local rule; otherwise
reported. Wrong: a fixed 0.01 "tiny" threshold would delete real micro-part features. vs: same as
Parasolid (short edge = inside vertex bubble), better than OCCT (Booleans fail through such edges).
Fix: PDQ thresholds (0.1 and 0.01) become export-profile warnings, not modelling errors.

**4.15 Export to 0.01 receivers (NX, Creo).** Own geometry max 1e-5: clean. Imported 0.01 data: same as
the receiver. Small features get the preflight warning. vs: better than OCCT (average uncertainty).
**4.16 Export to CATIA 0.001.** Own clean. Imported NX-quality edges (0.01) exceed 0.001, so CATIA join
fails. Fix: preflight refits edges above the receiver tolerance, lists residuals; uncertainty = max tol.
vs: better than all three, none of which targets a receiver.

## 5. Revised rule set

R0 Exact rational entities: no stored tolerance, no size limit; behave as eps(R) only when compared with
f64 entities. (Evidence: Parasolid exact entities = half precision [read-excerpt]; inference.)

R1 f64 resolution floor eps(R) = max(1e-7, 1e-11 R); operations run in a frame centred on the operands
(exact translation). Drop the 1e4 mm limit. Warn when a single body exceeds 1 km at export.
(Parasolid 1e11 ratio and box; Ansys large mode; OCCT forum rescaling advice.)

R2 Own fitted: fit target max(1e-6, 10 eps), tightened locally to L/100 (not below 10 eps); stored
= max(floor, 1.1(gap + delta)); cap max(1e-5, 100 eps). Over the cap: operation refused.

R3 Imports: measure every entity; declared value is a hint only. Auto cap 0.01 (edges, faces), 0.02
(vertices). Above: refit; then opt-in up to 0.1 (vertices 0.2). Report bands. (NX 0.01 and 0.0141
tolerant vertices; receivers at 0.01; CATIA heal max 0.1; Rhino declared-tol failures.)

R4 Margin 1.1 x (measured gap + error bound), applied to measurements only. Vertex >= edge >= face.
Vertex-first construction where three or more surfaces meet.

R5 Local rule, pairwise over non-adjacent entities: error t_a + t_b > d/2, warn > d/5.
cap_eff = min(C_class, L/4). BVH broad phase inflated by 5t.

R6 Short edges and slivers: invalid only by R5; collapse when the result is valid, else report.

R7 After each op: re-measure changed entities, lower by refit (pcurve first), never raise above the class
cap; warn when a tolerance rises over 10x the op's input maximum.

R8 Export with a receiver profile (0.01 NX/Creo/Rhino, 0.001 CATIA): refit entities above it, list
residuals, warn on features < 10x and < 4x receiver tol, write uncertainty = max entity tol, offer
scaled export for micro parts.

Relative tolerances: rejected for caps (receivers are absolute; Creo moved to absolute and recommends it
for assemblies; relative caps break mixed-size merges). Accepted only for the f64 floor (R1), where the
error really does scale with coordinate magnitude. SMLib's sqrt(size) growth has no stated derivation.

## 6. Where Catallax lands

Better: prismatic and sheet-metal own models (exact, zero growth, exact thickness); long operation chains;
mixed-size assemblies (vs Creo relative, SMLib); MEMS (vs Parasolid); near-tangent fillets (fails loudly
vs OCCT); exports to tighter receivers (preflight, honest uncertainty); large sites (local frames, exact
transforms, vs OCCT).
Same: brackets, watch gears, good imports, cars and wings (vs Parasolid), slivers.
Worse: convenience on poor imports and thin-walled poor imports (refusals where others give a broken but
displayable solid); file size and fitting time for large freeform models (1.8x Parasolid's segments,
about 10x an NX 0.0254 export); bodies over 1 km to Parasolid receivers (inherent to receivers).

## 7. Confidence

- eps(R) form: medium (ratio borrowed from Parasolid; Catallax's own conditioning untested).
- Pairwise local rule and BVH: medium-high (consistent with Jackson, Parasolid, OCCT failure reports).
- Import caps 0.01 / 0.02: medium-high (receiver defaults; NX tolerant vertex figure is a snippet).
- Near-tangency sqrt(2 r d) effect: high as arithmetic, untested in Catallax.
- Segment-count costs: low-medium (theory, not measured).
- OCCT GitHub issue details: low (snippets of recent issues, not opened).
- 1 um SolidWorks minimum feature: low (old forum posts).
