# Catallax

Geometry should not require a vendor's permission. Catallax is a geometry kernel in safe Rust, free for everyone, for ever.

Catallax is built for SolidDraft, the CAD application. They are one project in two repositories. [SolidDraft is public under MIT](https://github.com/BohnBawerick/soliddraft). It is not released for download yet; a release with the features is coming soon. For now it is research grade only. In design: SolidDraft runs only on Catallax, with no other kernel as a fallback. [The roadmap](docs/PLAN.md) orders the kernel work and its application integration.

Kernel site: [kernel.zigrawbohnbawerick.org](https://kernel.zigrawbohnbawerick.org/).

## Why this is free

Nearly every commercial CAD program runs on one of a handful of closed geometry kernels. Their owners licence them, meter them and can withdraw them. The maths inside is a century old and the algorithms are published. What is closed is the code, and the code is kept closed by law.

That law is not property. Property exists because things are scarce. If you copy this kernel, we still have it, and now you have one too. Nobody is poorer. Programs are not scarce, so the "intellectual property" over them has to be manufactured by the state: a grant that forbids you from arranging your own bits in a pattern someone else arranged first. Murray Rothbard wrote, "Patents, therefore, are grants of exclusive monopoly privilege by the State" in [Man, Economy, and State, chapter 10, section 7](https://mises.org/online-book/man-economy-and-state-power-and-market/10-monopoly-and-competition/7-patents-and-copyrights). Stephan Kinsella extended the point to copyright: "Ideas are not naturally scarce" in Against Intellectual Property, Journal of Libertarian Studies 15:2, 2001. We stand with that reading.

Geometry is mathematics and belongs to everyone. A public kernel lets a one-person CAD product enter the market without buying permission from a kernel vendor. Copy this code. Change it. Sell a product built on it. Compete with us.

Catallax puts AI coding agents to work on exact geometry, directed and checked by a human. Most of the inherited implementation was written that way. Whether this approach works at the scale a useful CAD kernel demands is an open question. We will find out in public.

## About the licence

Copyright attaches to code by default. We use the legal instrument available to let you use the code without asking us. MIT grants permission to use, copy, modify, publish and sell it, with the notice kept on copies and no warranty.

Licensed under MIT. See [LICENSE](LICENSE).

Catallax was started by Pedro De Paiva.

## Two kernels, one recipe

The recipe, the feature tree, is the master. An exact B-rep kernel and a signed distance field kernel already evaluate a shared subset of its operations. Each operation has one mathematical definition. Recipe builds name faces and edges by the operations that created them, not their position in a result list.

The direction is **exact by definition, approximate on purpose**. In design: the B-rep carries parametric parts, STEP and drawings; the field carries lattices, field-driven design and optimisation, and checks the B-rep. SolidDraft grows on both kernels as their operations land.

The kept rules are the same answer on every machine and at every thread count, no unsafe code, derived geometry stored by definition, and no wrong proven answer. Only impossible geometry is an acceptable final refusal. Other refusals are unfinished work. The destination remains a kernel comparable to OpenCASCADE, SMLib and the Rust kernels; target parts rebuilding in SolidDraft at click speed are milestones, not the finish line.

In design: one modelling tolerance defines sameness, separate from display and export accuracy. Approximation is reserved for decisions proven too close to call and shapes with no exact form the B-rep can carry. Decisions are recorded once, objects judged the same are merged, approximated pieces carry capped local tolerances measured against their definitions, and topology is checked before publication. The modelling tolerance value and the other answered choices are in [the roadmap](docs/PLAN.md#15-decided-on-22-october-2026); its [proposed PRD changes](docs/PLAN.md#14-proposed-prd-changes) remain proposals except for the tolerance and stand-in decisions.

[The PRD](docs/PRD.md) defines the product, and [the roadmap](docs/PLAN.md) orders the work. [RECIPE.md](docs/RECIPE.md) defines the shared language; the reader accepts versions 3 and 4. The version 4 additions that evaluate are one-sided through-all cuts, bounded face/next extrusions ([PR47](https://github.com/BohnBawerick/catallax/pull/47), [PR50](https://github.com/BohnBawerick/catallax/pull/50)) and sharp multi-face thicken with `round` roles ([PR86](https://github.com/BohnBawerick/catallax/pull/86)). Solved sketch roots, helix and spline paths, corner axes, scan values and the other version 4 definitions are specifications, not implementation. [FIELD.md](docs/FIELD.md) states the field's query guarantees and unresolved cases. [COMPARE.md](docs/COMPARE.md) documents the two-way sampled display check, not a whole-boundary or 1 nm certificate. In design: the tolerance and consistency core and the finer pre-export check. The inherited B-rep mark APIs do not yet enforce that geometric policy. The deterministic verification exchange uses work budgets; the legacy commands still have their 50-second deadline. Three work budgets are distinct contracts: the 1M allowance of the historical and `target-parts` scoreboards, the 20M application default and the 100M verification-request default. None stands in for another, and no result here is relabelled from one to another.

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
- `catallax-mesh` builds shared-edge display meshes with picking IDs, proven chord, angular and edge-length controls, aspect-ratio reports and OBJ export. Finite contractible torus patches use bounded minor-angle charts without raising the vertex cap; complementary and winding regions keep their periodic domains.
- `catallax-recipe` reads, checks and writes recipes in the canonical `catallax-recipe/3` and `catallax-recipe/4` JSON forms. Version 3 rejects version 4 end-condition fields and operations. It evaluates no geometry.
- `catallax-field` compiles recipe stock solids, exact sketch profiles, signed extrusions, two-sided and one-sided through-all cuts, rational-frame full/partial revolutions, named holes, linear and circular body patterns, mirrors and booleans into a fixed-order tape. Bounded face/next extrusions use definition-native planar first-hit cells, separated rectangular caps on stock spheres, cylinders and cones, and coaxial annular caps on stock tori. Names stay outside the tape. Sweep, cell and Boolean residuals carry conservative distance bounds, not a true-distance claim. Offset, shell, thicken and draft cells are built independently from their definitions in the domains listed below ([PR68](https://github.com/BohnBawerick/catallax/pull/68), [PR57](https://github.com/BohnBawerick/catallax/pull/57)). Native field display is a separate path from the application's B-rep display and is not accepted for interactive use: the matched P4 field display took 41.9 s in the [M1 note](docs/milestones/M1.md).
- `catallax-features::evaluate` builds named B-rep stock solids, exact line/circle/arc profiles with holes, signed extrusions, rationally represented full and partial revolutions, named holes, exact linear and circular body copies, mirrors, booleans and target-bounded through-all cuts. Bounded face/next extrusions support planar limiting boundaries, trimmed references and discontinuity curtains, a uniformly separated curved stock-quadric cap over one rectangle, and a coaxial annular cap on a stock torus. `Built` retains an immutable recipe definition and returns the model, result and intermediate body handles; `Model::face_origins`, `edge_name` and `vertex_name` expose the names. Split faces keep their origins and merged faces take their union. `catallax_app::pattern::expand_feature` saves explicitly bound feature instances as ordinary recipe steps.
- Exact elementary fillets build named cylinder strips between planes and torus strips between a plane and a perpendicular cylinder. Three equal convex strips between independent planes meet in a named sphere sector. `Built::fillets` retains their generating definitions. P3 construction, named display faces and radius-edit names are covered. Where three mutually perpendicular convex planes meet and all three straight edges are selected, the strips close with an exact `corner.vertex` sphere sector sewn to their shared quarter-circle boundaries; other three-strip junctions are refused with `fillet.unsupported-corner`. A straight strip whose edge ends on a cylindrical wall, such as a slot's straight top edge, runs past the edge and is cut by that wall along their exact intersection. Strips whose edges and outer supports are tangent at a vertex, such as the straight and arc edges of a slot's top boundary, join on their shared meridian circle with no corner face; in `catallax-recipe/4` a selected edge rolls over such tangent boundaries along its smooth chain, keeping its one blend name. A plane/cylinder contact away from the axis whose ring torus would be horn or spindle uses the regular outer sheet of that torus, a `TorusSheet`, written to STEP as a trimmed toroidal surface. Whole-part sampled comparison remains a separate acceptance gate.
- Exact elementary chamfers build named planar strips and conical rings, with ordered equal or unequal setbacks. `Built::chamfers` retains their input supports, edge, distances and exact contacts; the field constructs the affected regions independently. Roll-over, general corner patches and trimmed circular arcs remain unfinished.
- Offset, shell and thicken have exact analytic reductions in both evaluators ([PR68](https://github.com/BohnBawerick/catallax/pull/68)). Convex planar offsets use plane meetings and explicit cylindrical and spherical normal fans. Shell measures distance to kept trims, including open rims and consumed cavities. Smooth cylindrical and radial sweeps retain their circular rims and regularise folds. Later additions are analytic conical distance domains and frustum joins ([PR74](https://github.com/BohnBawerick/catallax/pull/74)), outward joins of convex coaxial plane, cylinder and sphere cells with complete circular trims ([PR78](https://github.com/BohnBawerick/catallax/pull/78)) and version 4 sharp multi-face thicken of proved analytic sheets ([PR86](https://github.com/BohnBawerick/catallax/pull/86)). General trimmed or spline envelopes and uncovered conical joins return typed coverage refusals. [The coverage table](docs/RECIPE.md#715-distance-feature-coverage) lists the implemented families and the legal geometry still refused. The P4 shelled tray builds at the 100M request default with valid 15-face STEP; it still refuses at the 1M scoreboard allowance and misses the interactive targets.
- Draft builds exact planar and axial cylindrical cases in both evaluators ([PR57](https://github.com/BohnBawerick/catallax/pull/57)): neutral-line plane replacement, single-nappe cylinder-to-cone reduction when the neutral normal is parallel to the pull, swallowed-face removal, neutral-plane retention with `keep_regions`, and planar caps against cylindrical neighbours. A convex planar body with mutually separated axial through bores between the same two caps drafts its planar walls and keeps every bore, its cap rings and their names unchanged, in any rigid placement. One axial cylindrical boss strictly inside a stationary cap stays the same way, with its rings, end and names. Other oblique cylindrical neutrals, a selected bore or boss, a bore or boss meeting an outer wall or another bore, a bore that stops inside the part, several bosses and a boss combined with a bore are unfinished coverage. Complete P1 with its drafted variant is not accepted.
- Field lattices compile periodic beam cells with constant and linear radius grades ([PR71](https://github.com/BohnBawerick/catallax/pull/71)). On explicit request the field writes a watertight STL and OBJ at a stated two-sided distance bound, measures the original field definition rather than the print mesh, and exports a mass-equivalent STEP stand-in of proxy solids plus a massless envelope. Ordinary `target-parts` adapters do not print or stand in: P6 there refuses or requires the explicit print request. Gyroids are not implemented, and a negative lattice operand in subtraction refuses the stand-in envelope.
- The field compiles one network gyroid with its frame fixed to the part frame and a constant dimensionless level, or a linear level graded along one part axis; the same linear grade drives a cubic lattice's radius. The terminal field stage shows either as an unproved preview and exports a labelled STL certified at a two-sided distance bound. Field moments and the explicit mass stand-in cover the gyroid and the axis-graded cubic lattice. Sheets and B-rep gyroid boundaries are not implemented.
- `catallax-compare` checks every selected display mesh vertex against certified field distance bounds, projects independent field boundary samples onto exact trimmed B-rep faces, and compares origin sets and outward directions. Reports distinguish agreement, proven disagreement and unresolved work. Edit regions and cooperative cancellation do not require app integration.
- `catallax-app` is the application entry. `build`, `build_incremental` and `build_soliddraft` return the exact body, a separately certified display mesh with edge polylines, and face, edge and vertex names; `build_certified` is the distinct certified path that adds geometry-bound inward offset prisms and original-definition properties on request. Construction and display share the 20M default allowance (`application.work-limit` on exhaustion); the display certifies chord error, not the whole boundary. The entry exists in this repository and in tests. SolidDraft's own integration and its end-to-end loop are tracked in the SolidDraft repository, not here.
- `catallax-cli` builds `catallax` and the legacy `exact-brep` program from the same main. They read a question, build and operate, measure, and write a manifest and optional STEP file.
- `catallax-verify` runs persistent question, recipe and sampled comparison requests with whole-request deterministic work budgets, separate wall-time reporting and distinct geometry tags. Its replay driver supports caching, strict STEP comparison and worker-count checks; the legacy commands are unchanged.
- Releases: `v0.1.0` is published as an [alpha prerelease](https://github.com/BohnBawerick/catallax/releases/tag/v0.1.0) ([PR79](https://github.com/BohnBawerick/catallax/pull/79), [PR81](https://github.com/BohnBawerick/catallax/pull/81)), unsigned, with no crates.io publication. Its figures are historical records in [the M1 note](docs/milestones/M1.md); they describe the sources named there, not later `main`.

Not yet: complete algebraic feature frames, oblique circular sweeps and some generated radii, listed with refusal codes in [RECIPE.md section 7.14](docs/RECIPE.md#714-b-rep-feature-geometry-bugs); general cone/torus branch tracing; and sound answers for every coplanar, equal-diameter crossing or near-tangent case. Supported operations can still refuse or take seconds. The kernel does not implement the whole PRD, and ignored regressions remain visible.

Open gates, unqualified on current `main`: complete P1 (it still refuses at `slot-cut-4` at 100M), P4 at the 1M scoreboard allowance, P0 general cone and torus branch tracing, any complete-click result inside the 1000 ms limit for the whole requested wait, the stricter 100 ms stock and 500 ms Boolean targets, and the 2000 ms background comparison. Display-only or construction-only timings, the historical [head-to-head report](docs/milestones/2026-10-24-head-to-head-1m.md) (source `74bf31a8`) and focused edit checks do not close them. Existing exceptions remain raw misses.

In design, not implemented: the tolerance and consistency core, whole-result geometric certificates, general non-planar first-hit extrusion limits, general fillets and corner continuation, chamfer continuation, the remaining distance-feature families, oblique cylindrical draft, sweep, loft, STEP import, gyroids and other field-only features, and SolidDraft's end-to-end verification loop. Settled feature definitions in RECIPE.md are specifications, not blanket implementation claims.

## Building

```sh
cargo build --release -p catallax-cli --bins
target/release/catallax about
target/release/catallax check question.json
target/release/catallax answer question.json --out answer
```

`check` validates without evaluating geometry. `answer` writes `answer/manifest.json` and, for a nonempty solid when requested, `answer/answer.step`. Refused and invalid questions also write a manifest and exit 0. [JUDGE-INTERFACE.md](docs/JUDGE-INTERFACE.md) defines the exact command and file contract, including its unchanged legacy names.

Local checks:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release
python3 -B -m unittest discover -s tools/edge-hunt -v
```

GitHub-hosted CI runs formatting, clippy, all enabled release tests and doctests. Crate `TESTS.md` files describe extra opt-in workloads. Timing benchmarks report measurements, never correctness gates.

## Verification

`catallax-verify` adds a persistent, deterministic exchange for questions, recipes and sampled two-kernel comparisons. Whole-request work budgets choose outcomes; wall time is reported separately. The legacy `catallax` and `exact-brep` commands retain their existing deadlines and file formats.

```sh
cargo build --release -p catallax-cli --bins
.claude/skills/verify-catallax/scripts/control-catallax.py doctor --binary target/release/catallax-verify --json
.claude/skills/verify-catallax/scripts/control-catallax.py --help
```

The [verification skill](.claude/skills/verify-catallax/SKILL.md) covers corpus replay, caching, resume, strict STEP comparison and matched-pair timing. The [versioned contract](docs/VERIFY-INTERFACE.md) defines work accounting and the expanded geometry kind tags. Pass private inputs as runtime paths and keep proof artifacts outside the repository.

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
See the [construction and snapback contract](docs/interface-changes/2026-10-27-construction-preview-and-snapback.md)
for representation types, output gates, budget accounting and exercised limits.


## Where to start

- [PRD](docs/PRD.md), [plan](docs/PLAN.md) and [research](docs/research/) explain the direction.
- [Recipe language](docs/RECIPE.md) defines operations, names and the recipe file.
- [Field queries](docs/FIELD.md) defines distance semantics, proof limits, selection and batch storage.
- [Design](docs/DESIGN.md) sets the rules; [architecture](docs/ARCHITECTURE.md) connects the B-rep crates.
- [Piece contracts](docs/pieces/) describe the current interfaces and mathematical tests.
- [Interface changes](docs/INTERFACE-CHANGES.md) records migrations.
- [Provenance](docs/PROVENANCE.md) records the publication audit and dependency licence boundary.

## Contributing

Read the design rules first. A guessed topological decision, an unproven construction fact, an `unsafe` block or a silent fallback is a defect. Keep API and format changes explicit, update callers and document the migration in a new [interface-change fragment](docs/interface-changes/README.md). Bring mathematical evidence for geometry changes and report refusals honestly.

Agent-written code is welcome. Say how it was produced and stay accountable for what it does. This repository's independent tester was built without seeing the kernel's code.

## Special thanks

This kernel stands on work that people put into the world for anyone to read: papers, retrospectives, source code and standards, published for free, which is how knowledge should move. None of it entered Catallax as code. What we took is what they wrote down.

Hanno Braun worked on [Fornjot](https://github.com/hannobraun/fornjot), an early-stage B-rep CAD kernel in Rust, from July 2020 until he shut it down in June 2026. His retrospective, ["Shutting Down Fornjot"](https://archive.hannobraun.com/fornjot/blog/shutting-down-fornjot/), names the mistakes he saw: extrapolating from early success, selling a dream before there was a product, sticking to incremental improvements, half-measures in a crisis, a vision that blurred once the application was cut, and prototyping that came too late. [The roadmap](docs/PLAN.md#2-fornjots-lessons-and-what-each-one-changes) applies those lessons to our build order: measure hard cases and real parts in SolidDraft, keep a works-today list beside sponsorship asks, prototype when an unknown answer would change the build, and end each work-based milestone with work committed or parked.

To be plain about what this is and is not: Catallax shares no code with Fornjot. It is not affiliated with Hanno Braun and he has not endorsed it. The same holds for every project and author below: their published ideas and results, not their source, with thanks.

- Matt Keeter's [Fidget](https://github.com/mkeeter/fidget) and [libfive](https://libfive.com/). The field kernel's fixed-order tapes, independent point evaluation and reusable evaluator storage follow their design; see [FIELD.md](docs/FIELD.md).
- [OpenCASCADE](https://dev.opencascade.org/), which opens and checks every STEP file Catallax exports. It never runs inside the kernel.
- [monstertruck](https://github.com/virtualritz/monstertruck), [Nacre](https://github.com/elgar328/nacre) and [ogeom-rs](https://github.com/gilbertorconde/ogeom-rs), Rust B-rep kernels whose public code was a reference point when the B-rep pieces were compared against existing implementations; see [PROVENANCE.md](docs/PROVENANCE.md).
- Raghothama and Shapiro, ["Consistent updates in dual representation systems"](https://doi.org/10.1016/S0010-4485(00)00036-1) (2000), the method behind checking that two kernels model the same solid.
- Rossignac and Requicha, "Constant-radius blending in solid modeling" (1984) and ["Offsetting operations in solid modelling"](https://faculty.cc.gatech.edu/~jarek/papers/Offsets.pdf) (1986), the rolling-ball fillet as an offset construction.
- Farouki and Sverrisson, ["Approximation of rolling-ball blends for free-form parametric surfaces"](https://doi.org/10.1016/0010-4485(96)00008-5) (1996).
- Jackson, [Parasolid's tolerant-modelling paper](https://doi.org/10.1145/218013.218067) (1995), and Hedberg et al., ["Testing the digital thread"](https://doi.org/10.1115/1.4032697) (2016), the reasons behind the tolerance caps in [the roadmap](docs/PLAN.md#34-local-tolerance-caps).
- [STEP Tools](https://www.steptools.com/docs/stp_aim/html/), public documentation of the STEP AP242 entities the exporter writes.

The research notes in [docs/research](docs/research/) record the wider reading behind the design.

## Other products

[check out my other products](https://zigrawbohnbawerick.org/)

## Sponsorship

Support the free kernel through [GitHub Sponsors](https://github.com/sponsors/BohnBawerick) or [Ko-fi](https://ko-fi.com/pedropaiva56920).

Nobody should sponsor a dream. This is what the code on `main` does today, and what it does not. The [full list](#what-works-today) above carries the caveats.

Works today:

- Blocks, cylinders, spheres, cones and tori, built exactly and by definition in safe Rust.
- Unite, subtract and intersect on plane, cylinder and sphere pairs, plus closed-form cone and torus cases, not arbitrary poses.
- Holes, planar and bounded stock-quadric extrusion end conditions, patterns, mirror, exact elementary fillets and chamfers, and analytic offset, shell, thicken and planar or axial-cylindrical draft, each only in the domains of the list above.
- An application entry that returns body, display, names and typed refusals, and field lattices with explicit print and mass-stand-in output.
- Exact algebraic placements and dimensions, rational-degree rotation construction and elliptic-cylinder definitions. Complete feature and display coverage does not follow from those definitions.
- Sketch profiles of lines, arcs and circles, extruded or revolved in the supported rational subset, with faces named by the operation that made them.
- Volume, area, centroid and unit-density centroidal inertia with proven bounds.
- STEP AP242 export with a proven deviation report, and display meshes with proven chord error.
- A signed distance field evaluator for the same solids, sketches and booleans, and a sampled two-way check between the two kernels.
- A deterministic verification CLI for questions, recipes and sampled comparisons, with replay and strict STEP checks.

Not yet: complete cone/torus and algebraic-feature coverage, all coplanar and near-tangent booleans, full P1 or P4 at their scoreboard allowances, or consistently interactive rebuilds. Valid cases still refuse or take seconds. The legacy commands can still reach their 50-second deadline; `catallax-verify` uses deterministic work allowances instead.

In design: general fillets and corner continuation, chamfer continuation, the remaining offset/shell/thicken families, oblique draft, sweep, loft, STEP import, gyroids and other field-only features, the tolerance and consistency core, and SolidDraft integration. [The roadmap](docs/PLAN.md) sets the build order and work-based milestones.
