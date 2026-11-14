# Kernel implementation reference

Moved from the README at the 0.2.0 source freeze, 14 November 2026. This reference preserves the detailed operation domains and display contracts. The [0.2.0 release scope](https://github.com/BohnBawerick/catallax/blob/release/0.2.0/docs/milestones/M2.md) records later qualification results and parked work. Historical measurements below retain their original sources.

## Two kernels, one recipe

The recipe, the feature tree, is the master. An exact B-rep kernel and a signed distance field kernel already evaluate a shared subset of its operations. Each operation has one mathematical definition. Recipe builds name faces and edges by the operations that created them, not their position in a result list.

In design: the B-rep carries parametric parts, STEP and drawings; the field carries lattices, field-driven design and optimisation, and checks the B-rep. SolidDraft grows on both kernels as their operations land.

The kept rules are the same answer on every machine and at every thread count, no unsafe code, derived geometry stored by definition, and no wrong proven answer. Only impossible geometry is an acceptable final refusal. Other refusals are unfinished work. The destination remains a kernel comparable to OpenCASCADE, SMLib and the Rust kernels; target parts rebuilding in SolidDraft, and then rebuilding at click speed, are milestones, not the finish line.

In design: one modelling tolerance defines sameness, separate from display and export accuracy. Approximation is reserved for decisions proven too close to call and shapes with no exact form the B-rep can carry. Decisions are recorded once, objects judged the same are merged, approximated pieces carry capped local tolerances measured against their definitions, and topology is checked before publication. The modelling tolerance value and the other answered choices are in [the roadmap](archive/PLAN-2026-11.md#15-decided-on-22-october-2026); its [proposed PRD changes](archive/PLAN-2026-11.md#14-proposed-prd-changes) remain proposals except for the tolerance and stand-in decisions.

[The PRD](PRD.md) defines the product, and [the roadmap](../ROADMAP.md) orders the work. [RECIPE.md](RECIPE.md) defines the shared language; the reader accepts versions 3, 4 and 5. Version 5 adds [all-faces references and explicit migration](interface-changes/2026-11-02-all-face-references.md). The version 4 additions that evaluate are one-sided through-all cuts, bounded face/next extrusions ([PR47](https://github.com/BohnBawerick/catallax/pull/47), [PR50](https://github.com/BohnBawerick/catallax/pull/50)) and sharp multi-face thicken with `round` roles ([PR86](https://github.com/BohnBawerick/catallax/pull/86)). Rational spline paths have checked schema support without sweep construction, and loft correspondence has an explicit versioned migration guard. The remaining version 4 definitions keep their implementation limits in RECIPE.md. [FIELD.md](FIELD.md) states the field's query guarantees and unresolved cases. [COMPARE.md](COMPARE.md) documents the two-way sampled display check, not a whole-boundary or 1 nm certificate. In design: the tolerance and consistency core and the finer pre-export check. The inherited B-rep mark APIs do not yet enforce that geometric policy. The deterministic verification exchange uses work budgets; the legacy commands still have their 50-second deadline. Three work budgets are distinct contracts: the 1M allowance of the historical and `target-parts` scoreboards, the 20M application default and the 100M verification-request default. None stands in for another, and no result here is relabelled from one to another.

## What works today

The code is early, and supported operations can still refuse. Measure it by its results, not the length of its roadmap.

The [five-tool target-part comparison](https://kernel.zigrawbohnbawerick.org/comparison.html) records pinned versions, construction and edit timings, and per-part correctness checks.

- `catallax-numbers` supplies exact arithmetic, certified enclosures and decisions, polynomial roots, construction facts and marks. The shared algebraic foundation adds `RealPoly` and certified sine/cosine construction for rational-degree angles.
- `catallax-shapes` stores exact algebraic placements and analytic dimensions, including elliptic-cylinder definitions, with proven parameter enclosures. Derived points retain their definitions. These definitions do not imply complete recipe coverage for every kind.
- `catallax-topology` builds blocks, cylinders, spheres, cones and tori in append-only curved topology and checks structural faults.
- `catallax-intersect` handles plane, cylinder and sphere pairs, closed-form cone/torus configurations, and cylinder/torus and cone/torus pairs with parallel, distinct axes. Branch-aware queries and certified finite traces retain disconnected components, shared nodes, cusp limits, seams and trim endpoints. Skew cylinder/torus and cone/torus tracing remain unfinished.
- `catallax-boolean` implements unite, subtract and intersect through intersect, imprint, classify, keep and sew. Torus regions bounded by finite parallel and meridian arcs retain their seams and endpoints; P3 center sections work through the junction-radius edit. Ambiguous and unsupported cases remain unfinished work.
- `catallax-measure` computes proven body mass properties, face centroids and areas, edge lengths and half-arc-length points, and body/face axis-aligned bounds. Construction does not measure: properties run when asked, and the legacy commands still measure every answer. Mass properties on finite offset and rounded faces are proven ([PR72](https://github.com/BohnBawerick/catallax/pull/72)); the field measures lattices separately (below).
- `catallax-step` exports deterministic STEP AP242 with analytic geometry, certified spline fits for intersection branches and a deviation report. The legacy contract has interoperability limits, including no seam curves.
- `catallax-import` supports clean analytic STEP bodies through proposal, application, persisted-record rebuild, display and checked re-export ([application proof](https://github.com/BohnBawerick/catallax/pull/255)). Convex planar repairs retain source-bound whole-face local tolerances ([repair proof](https://github.com/BohnBawerick/catallax/pull/266)). SolidDraft can [import these bodies](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/141) and [build features on them](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/161). General curved repairs, periodic seam exchange, spline and multiple-solid recipe imports remain unfinished.
- `catallax-mesh` builds shared-edge display meshes with picking IDs, proven chord, angular and edge-length controls, aspect-ratio reports and OBJ export. Finite contractible torus patches use bounded minor-angle charts without raising the vertex cap; complementary and winding regions keep their periodic domains.
- STEP spline source geometry is read as exact rational knots, Euclidean controls and separate positive weights. `catallax_mesh::preview_step` displays the retained source through an explicitly unproved `PreviewMesh` with shared edges and STEP entity picking keys. Certified spline imports use a separate [admission route](interface-changes/2026-11-12-step-spline-certification.md) and still need complete solid checks. Curved repair and printing keep their own proof requirements. See the [source-preview contract](interface-changes/2026-11-07-step-spline-source-preview.md) and [third-party notices](THIRD-PARTY.md).
- `catallax-recipe` reads, checks and writes canonical `catallax-recipe/3`, `catallax-recipe/4` and `catallax-recipe/5` JSON, including optional `materials` region bindings to body steps. Version 3 rejects version 4 end-condition fields, operations, path corner axes and spline paths. Versions 3 and 4 require face-origin patterns; version 5 also permits all-faces references. Rational spline-path validation and versioned loft migration have a [case inventory](interface-changes/2026-11-10-rational-path-schema.md). The crate evaluates no geometry.
- `catallax-field` compiles recipe stock solids, exact sketch profiles, signed extrusions, two-sided and one-sided through-all cuts, rational-frame full/partial revolutions, named holes, linear and circular body patterns, mirrors and booleans into a fixed-order tape. Bounded face/next extrusions use definition-native planar first-hit cells, separated rectangular caps on stock spheres, cylinders and cones, and coaxial annular caps on stock tori. Names stay outside the tape. Sweep, cell and Boolean residuals carry conservative distance bounds, not a true-distance claim. Offset, shell, thicken and draft cells are built independently from their definitions in the domains listed below ([PR68](https://github.com/BohnBawerick/catallax/pull/68), [PR57](https://github.com/BohnBawerick/catallax/pull/57)). Native field display is a separate path from the application's B-rep display and is not accepted for interactive use: the matched P4 field display took 41.9 s in the [M1 note](milestones/M1.md).
- `catallax-features::evaluate` builds named B-rep stock solids, exact line/circle/arc profiles with holes, signed extrusions, rationally represented full and partial revolutions, named holes, exact linear and circular body copies, mirrors, booleans and target-bounded through-all cuts. Bounded face/next extrusions support planar limiting boundaries, trimmed references and discontinuity curtains, a uniformly separated curved stock-quadric cap over one rectangle, and a coaxial annular cap on a stock torus. `Built` retains an immutable recipe definition and returns the model, result and intermediate body handles; `Model::face_origins`, `edge_name` and `vertex_name` expose the names. Split faces keep their origins and merged faces take their union. `catallax_app::pattern::expand_feature` saves explicitly bound feature instances as ordinary recipe steps.
- Recipe sweeps build exact line chains in fixed orientation and meridional line/cardinal-arc chains in rotation-minimising orientation. Perpendicular joins include the entire stationary round turn; reversals use a saved `corner_axis`. Regularised union removes overlaps and zero-volume remnants, and maximal faces keep all generating origins. Display, requested properties and bounded STEP use existing analytic supports. See the [round-turn fixture](../crates/catallax-features/tests/fixtures/round-turn-sweep.json) and [supported domain](RECIPE.md#411-sweep). General curved motion remains unfinished.
- Version 4/5 rational lofts join planar line, circular-arc and exact spline profiles through the B-rep application entry, with ruled or uniform-station natural interpolation. Exact coefficient and control-hull reductions preserve analytic geometry where possible; numerical chart cuts are not physical picking edges. [The case inventory](LOFT-CASES.md) lists successful operation cases, whole-body admission proofs and typed coverage gaps. Full-circle sections use exact rational charts, so the natural three-circle lofts build with display inside the unchanged 20M application budget. General intersection solving, legacy arc-length evaluation and field loft compilation remain unfinished.
- Exact elementary fillets build named cylinder strips between planes, at straight generators between parallel cylinders, and at straight generators between a cylinder and a plane parallel to its axis, and torus strips between a plane and a perpendicular cylinder. A circle where any other two planes, cylinders, cones, spheres or tori about one axis meet, such as a frustum or chamfer rim, a countersink, a dome's rim or a revolved groove, builds an exact torus strip too, which classifies, displays, measures and writes to STEP. A closed ellipse where a plane cuts a cylinder obliquely builds one named exact rational band, which classifies and displays, answers volume, mass properties, whole-face area, centroid and bounds, measures zero against its centre ellipse in `measure_face`, and writes to STEP. A closed ellipse where a plane cuts a cone obliquely builds the same kind of band, exact and rational of degree (6, 2), whose cone rail is one of its rows. It classifies, displays, answers volume and mass properties, and writes to STEP. A closed ring where a narrow cylinder passes through a wide one at any angle builds one named fitted band, within 100 nm of the rolling ball, which classifies, displays, answers volume, area, centroid and body box, carries a manifest/4 certificate and writes to STEP. Mass properties and measures of the band's own faces still refuse. Where the narrow cylinder runs past the wide one's side, the two cross in one closed loop that winds round neither cylinder. That loop builds the same kind of band, which classifies, displays, answers volume, area, centroid and body box, carries a manifest/4 certificate and writes to STEP. Three equal strips between independent planes, all convex or all concave, meet in a named sphere sector, such as a pocket's floor corner. `Built::fillets` retains their generating definitions. P3 construction, named display faces and radius-edit names are covered. Where three mutually perpendicular convex planes meet and all three straight edges are selected, the strips close with an exact `corner.vertex` sphere sector sewn to their shared quarter-circle boundaries. A straight strip between planes is checked against the body's other faces only between its end faces, so the edges of pockets, steps and bosses round. Two straight strips of one convexity meeting at a vertex where their outer faces lean differently meet on an exact step: the deeper strip runs past the other's outer rail and is cut by that rail's face. A straight strip whose edge ends on a cylindrical or spherical wall, such as a slot's straight top edge, runs past the edge and is cut by that wall along their exact intersection. An arc strip ends at a meridian circle on a planar end face through its axis, as on a half-disc. Other junctions, such as convex and concave strips at one vertex or curved strips meeting at a sharp vertex, refuse with `fillet.unsupported-corner`, and strip ends not built yet, such as an arc strip on an end plane off its axis, with `fillet.unsupported-end`. `catallax_features::fillet_corner::SphereSector` builds the sphere sector of any ball touching three supports for other junctions to reuse. Strips whose edges and outer supports are tangent at a vertex, such as the straight and arc edges of a slot's top boundary, join on their shared meridian circle with no corner face; in `catallax-recipe/4` a selected edge rolls over such tangent boundaries along its smooth chain, keeping its one blend name. A chain rolls over an earlier fillet's cylinder strip the same way, over its end as a torus strip, or as a sphere strip whose plane contact is one point where the two radii are equal, as on a block's top after its uprights were rounded or a pocket's floor after its corners were. Rounding every edge leaves out the smooth rails of earlier fillets, which their retained definitions identify. A strip ending where two earlier strips of its radius meet at a mitre, such as an upright under a rounded top loop, closes with the corner the three would have closed together, and the mitre goes, also when every edge is selected. A ball wider than the earlier round refuses with `fillet.radius-too-large`, and a mitre of another radius with `fillet.unsupported-corner`. A plane/cylinder contact away from the axis whose ring torus would be horn or spindle uses the regular outer sheet of that torus, a `TorusSheet`, written to STEP as a trimmed toroidal surface. A radius proved too large for its edge refuses with `fillet.radius-too-large` and, where the proof gives one, the bound past which every radius is too large. Whole-part sampled comparison remains a separate acceptance gate.
- Exact elementary chamfers build named planar strips and conical rings, with ordered equal or unequal setbacks. `Built::chamfers` retains their input supports, edge, distances and exact contacts; the field constructs the affected regions independently. Roll-over, general corner patches and trimmed circular arcs remain unfinished.
- Offset, shell and thicken have exact analytic reductions in both evaluators ([PR68](https://github.com/BohnBawerick/catallax/pull/68)). Convex planar offsets use plane meetings and explicit cylindrical and spherical normal fans. Native inward offsets also accept nonconvex polygonal planar solids, with finite edge distances, rounded reentrant corners, holes and Boolean-trimmed faces. Shell measures distance to kept trims, including open rims and consumed cavities. Smooth cylindrical and radial sweeps retain their circular rims and regularise folds. Later additions are analytic conical distance domains and frustum joins ([PR74](https://github.com/BohnBawerick/catallax/pull/74)), outward joins of convex coaxial plane, cylinder and sphere cells with complete circular trims ([PR78](https://github.com/BohnBawerick/catallax/pull/78)) and version 4 sharp multi-face thicken of proved analytic sheets ([PR86](https://github.com/BohnBawerick/catallax/pull/86)). Outward nonconvex planar offsets, general trimmed curved or spline normal envelopes and uncovered conical joins return typed coverage refusals. [The coverage table](RECIPE.md#715-distance-feature-coverage) lists the implemented families and the legal geometry still refused. The P4 shelled tray builds at the 100M request default with valid 15-face STEP; it still refuses at the 1M scoreboard allowance and misses the interactive targets.
- Draft builds exact planar and axial cylindrical cases in both evaluators ([PR57](https://github.com/BohnBawerick/catallax/pull/57)): neutral-line plane replacement, single-nappe cylinder-to-cone reduction when the neutral normal is parallel to the pull, swallowed-face removal, neutral-plane retention with `keep_regions`, and planar caps against cylindrical neighbours. A convex planar body with mutually separated axial through bores between the same two caps drafts its planar walls and keeps every bore, its cap rings and their names unchanged, in any rigid placement. One axial cylindrical boss strictly inside a stationary cap stays the same way, with its rings, end and names, alone, beside one separate through bore in either history order, or with one newer through bore strictly inside its footprint that opens on its end. One countersunk through hole stays the same way: its wall, its cone, the ring they share and both cap openings keep their exact trims and names. Other oblique cylindrical neutrals, a selected bore, cone or boss, a bore or boss meeting an outer wall or another bore or boss, a bore that stops inside the part, several bosses, several bores beside a boss and a countersunk hole beside another bore or a boss are unfinished coverage. Complete P1 with its drafted variant is not accepted.
- Field lattices compile periodic beam cells with constant and linear radius grades ([lattice construction](https://github.com/BohnBawerick/catallax/pull/71)). On explicit request the field writes a watertight STL and OBJ at a stated two-sided distance bound, measures the original field definition rather than the print mesh, and exports a mass-equivalent STEP stand-in of proxy solids plus a massless envelope. Ordinary `target-parts` adapters do not print or stand in: P6 there refuses or requires the explicit print request. A negative lattice operand in subtraction still refuses the stand-in envelope.
- The field compiles one [network gyroid](https://github.com/BohnBawerick/catallax/pull/207) with its frame fixed to the part frame and a constant dimensionless level, or a [linear level graded along one part axis](https://github.com/BohnBawerick/catallax/pull/222); the same linear grade drives a cubic lattice's radius. The terminal field stage shows either as an unproved preview and exports a labelled STL certified at a two-sided distance bound. Field moments and the explicit mass stand-in cover the gyroid and the axis-graded cubic lattice. Sheets and B-rep gyroid boundaries are not implemented.
- `catallax-compare` checks every selected display mesh vertex against certified field distance bounds, projects independent field boundary samples onto exact trimmed B-rep faces, and compares origin sets and outward directions. Reports distinguish agreement, proven disagreement and unresolved work. Edit regions and cooperative cancellation do not require app integration.
- `catallax-app` is the application entry. `build`, `build_incremental` and `build_soliddraft` return the exact body, a separately certified display mesh with edge polylines, and face, edge and vertex names; `build_certified` is the distinct certified path that adds geometry-bound inward offset prisms and original-definition properties on request. Construction and display share the 20M default allowance (`application.work-limit` on exhaustion); the display certifies chord error, not the whole boundary. SolidDraft's [Catallax integration](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/3), [SMLib removal](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/6) and [hosted end-to-end application loop](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/5) have landed. Complete target-part and strict latency acceptance remain separate.
- `catallax-cli` builds `catallax` and the legacy `exact-brep` program from the same main. They read a question, build and operate, measure, and write a manifest and optional STEP file.
- `catallax-verify` runs persistent question, recipe and sampled comparison requests with whole-request deterministic work budgets, separate wall-time reporting and distinct geometry tags. Its replay driver supports caching, strict STEP comparison and worker-count checks; the legacy commands are unchanged.
- Historical release: `v0.1.0` was published as an [alpha prerelease](https://github.com/BohnBawerick/catallax/releases/tag/v0.1.0) ([PR79](https://github.com/BohnBawerick/catallax/pull/79), [PR81](https://github.com/BohnBawerick/catallax/pull/81)), unsigned, with no crates.io publication. Its figures are historical records in [the M1 note](milestones/M1.md); they describe the sources named there, not later `main`.

Not yet: complete algebraic feature frames, oblique circular sweeps and some generated radii, listed with refusal codes in [RECIPE.md section 7.14](RECIPE.md#714-b-rep-feature-geometry-bugs); general cone/torus branch tracing; and sound answers for every coplanar, equal-diameter crossing or near-tangent case. Supported operations can still refuse or take seconds. The kernel does not implement the whole PRD, and ignored regressions remain visible.

Open gates, unqualified on current `main`: complete canonical P1 acceptance, P4 at the 1M scoreboard allowance, P0 general cone and torus branch tracing, and complete target-part runs meeting the separate selectable-preview and background-proof budgets. The historical targets are 100 ms for stock, 500 ms for Booleans, 1000 ms for feature clicks and 2000 ms for background proof. Display-only or construction-only timings and focused edit checks do not close them. The [historical 1M report](milestones/2026-10-26-head-to-head-1m.md) retains its named source and refusals, not measurements of later `main`. The [milestone ledger](archive/PLAN-2026-11.md#13-in-flight-ledger) records the open gates; timing requirements moved to M26 on 11 November 2026.

Still unfinished: general tolerance/consistency-policy composition beyond the [first certified offset consumer](https://github.com/BohnBawerick/catallax/pull/53), whole-result certificates for every supported operation, general non-planar first-hit extrusion limits, remaining general-fillet and corner continuation, chamfer continuation, the remaining distance-feature families, oblique cylindrical draft, general curved sweeps, general loft intersection coverage, wider STEP import/repair domains and other field-only families. Analytic STEP import, network gyroids and SolidDraft's application loop are implemented in the scoped domains above. Settled feature definitions in RECIPE.md are specifications, not blanket implementation claims.

## Display meshes

`catallax_mesh::mesh_body(&model, body, &options, &mut decider)` returns an
`Answered<Mesh>`. `MeshOptions::chord_tolerance` is a positive exact length in
millimetres. `max_vertices` limits deterministic refinement, including temporary
periodic chart copies. The default is 0.01 mm and one million vertices.
Exact planar boundary sampling also has a per-mesh allowance of 4096 retained
samples, checked before constructing each sample's expression graph. Exhaustion
returns `mesh.resource-limit` without a partial mesh or a coarser tolerance.
This bounds the expensive retained boundary data separately from binary64
vertices. It is not a coordinate-size limit or a process-wide heap limit.

Optional `angular_deviation` bounds the angle in radians between each facet normal
and every analytic normal over its parameter triangle. `max_edge_length` bounds
every triangle edge and boundary-polyline segment in millimetres. Both bounds
include the returned binary64 coordinates. Nonpositive optional controls return
`mesh.invalid-options`.

Each face reports `max_aspect_ratio`, an optional upper bound on longest edge
divided by shortest altitude. `None` means rounding collapsed a facet to zero
area. `max_aspect_ratio` in the options requests a ratio but does not force
quality refinement across fixed boundaries. Unmet ratios and undefined facet
normals appear in `Mesh::unavailable_controls`, with the control and reason. A requested control
must either be certified or appear there; a resource or precision failure still
refuses the mesh. No option changes exact geometry or STEP output.

`Mesh::positions` is shared across faces. Each `Triangle` carries its input
`FaceId`; each `FaceMesh` identifies a contiguous triangle range and its proven
maximum chord error, aspect-ratio bound and optional angular bound. `EdgeMesh`
records the one polyline used by both incident faces, including its proven
curve-to-polyline error. Face IDs belong to the
input model, not to rebuild-stable recipe names.

Planes, circular and elliptic cylinders, cones, spheres and tori use their exact
definitions. Intersection edges, including finite parabola and hyperbola arcs,
use piece 4's certified traces. Surface bounds use Taylor's theorem and analytic
second-derivative bounds over parameter triangles. An elliptic cylinder uses its
semi-axis-scaled angular chart; its major semi-axis bounds positional curvature,
and the semi-axis ratio bounds normal rotation.
For an elliptic-cylinder strip bounded by two angularly monotone closed loops,
or a finite cylindrical fillet between two monotone end curves, the mesher joins
the existing edge samples in angular order. Strip vertex certificates use the
shared surface sampler; every triangle passes the same orientation, chord and
requested-control checks. Other trims or controls needing
interior vertices use the periodic triangulator. Toroidal bands use a bounded
minor-angle chart and the shared certified surface sampler. With an angular
control, interior columns follow the shared rails' major-angle parameters and
minor-angle subdivisions resolve normal rotation, even when the tube radius is
below the chord tolerance. This avoids staggered major-circle chords that
centroid refinement cannot repair on a thin band. Shared boundary samples are
not moved or split. A face-local diagonal cannot reuse another face's interior
chord between boundary vertices; the mesher refines that chord on the analytic
face instead.
Bounds include replacement of boundary vertices and binary64 coordinate
rounding. They describe chord error against the analytic face and edge
definitions, not a proof that polygonal trims equal the exact curved region.
All covered definitions have proven bounds, so no sampled estimate is presented
as a proof. Input and trace marks are preserved.

Outward-rounded double intervals filter surface evaluation, trace Taylor
coefficients, chord checks, aspect ratios and angular certificates. Angular
filters retain the whole-parameter normal-variation bound and include the
returned coordinates; uncertain cases retain the ball or exact-rational path.
Closed cone and bounded-torus bands with two horizontal sampled rails use
angularly ordered layers between their shared edge samples. Interior columns
follow the rails' longitude parameters; chord and normal-rotation bounds select
the layer count. Each facet passes the same complete rounded-coordinate
certificates. An uncertified proposal is discarded before the general periodic
triangulator runs. Other trims and apices retain that triangulator.
A triangle whose three vertices belong to one shared rail receives a certified
interior Steiner point instead of becoming a zero-volume fin between adjacent
faces. Every resulting triangle still passes the same chord certificate and
work limits.

The implementation is serial with fixed traversal, subdivision and output
orders. Pure-Rust angle projection and dyadic periodic coordinates avoid
platform trigonometric libraries and inconsistent seam copies. Concurrent
callers do not share state. Tests compare OBJ bytes with 1, 2 and 4 workers,
count every mesh edge twice with opposite directions, and sample triangle
interiors against independent analytic distance formulas.

Write an OBJ for Blender without changing the judge CLI:

```sh
cargo run --release -p catallax-mesh --example obj -- torus torus.obj 0.05
cargo run --release -p catallax-mesh --example obj -- cross intersection.obj 0.05
```

The example also accepts `block`, `cylinder`, `cone`, `frustum`, `sphere`,
`drill` and `lens`. Library callers write any built body with `Mesh::write_obj`.
OBJ groups are named `face_<id>`, and their triangles share vertex indices.
Coordinate-precision and vertex-budget failures return `mesh.precision-or-budget`;
the exact boundary sample allowance returns `mesh.resource-limit`, and
invalid controls return `mesh.invalid-options`. No over-tolerance or
non-watertight mesh is returned as a successful answer. A caller-supplied
wall-clock deadline can still cancel meshing through `Decider`; the default
library decider has no deadline.

Application entries `catallax_app::build`, `build_soliddraft` and `build_certified`
share a 20,000,000-unit aggregate arithmetic allowance for construction and
display unless the caller already owns an operation budget. Default exhaustion
returns `application.work-limit` without a body. An explicit caller allowance is
never reset or enlarged by a nested build. This work limit is separate from the
display sample/vertex limits and from all geometric tolerances.
Planar orientation predicates charge this allowance even when interval filtering
settles their signs without exact arithmetic. Bridge and ear searches therefore
cannot bypass cancellation or the work limit through their filtered fast path.

`catallax_app::lifecycle::Controller` can deliver a revision-bound construction
display before exact B-rep evaluation, using the independent field compiler.
Unsupported definitions or incomplete displays refuse rather than showing a
substitute shape. The caller reserves one preview/proof allocation and schedules
the returned proof job. Confirmation retains an incremental exact checkpoint and
certified display. Failure restores the last proved display, but does not allow
that older body to be exported as the failed input. Cancellation invalidates queued
work. Selection origins resolve only against the same revision's confirmed result.
See the [construction and snapback contract](interface-changes/2026-10-27-construction-preview-and-snapback.md)
for representation types, output gates, budget accounting and exercised limits.
