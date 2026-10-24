# Catallax roadmap

Written 22 October 2026, against `main` after the shared algebraic foundation (#18) and the deterministic verification exchange (#12) merged. [PRD.md](PRD.md) defines the product. This file replaces the previous `docs/PLAN.md` and keeps its backbone. Every item of the previous plan appears below with a verdict, or moves word for word to [COMPARE.md](COMPARE.md) when it documents the comparison library rather than ordering work.

The companion document is SolidDraft's roadmap in the SolidDraft repository. Catallax and SolidDraft are one project in two repositories. This file orders the kernel work; SolidDraft's roadmap orders the application work, and each names the other's items where they depend on each other.

The plan was rewritten after reading Hanno Braun's retrospective ["Shutting Down Fornjot"](https://archive.hannobraun.com/fornjot/blog/shutting-down-fornjot/). His lessons are added to the existing plan, not used to replace it. Where a recommendation would change the PRD, it is listed under [Proposed PRD changes](#14-proposed-prd-changes) and waits for the project owner.

A later study of SolidCraft, a parametric CAD application on a B-rep kernel, added deliverables and exit criteria to the items in section 12 and the statement in [section 3.8](#38-what-we-take-from-solidcraft-and-what-we-do-not). It moved no item and reordered nothing.

## 1. Premise

Catallax is SolidDraft's kernel, and SolidDraft runs only on Catallax. Kernel work is chosen by where SolidDraft is going, not by what the application does today. SolidDraft grows step by step toward the most usable parts of SolidWorks and nTopology: parametric solids, drawings and assemblies on one side, and implicit modelling, lattices, field-driven thickness, topology optimisation and simulation-led shapes on the other.

The destination stays where the PRD puts it, a kernel comparable to OpenCASCADE, SMLib and the Rust kernels. Target parts that rebuild in SolidDraft at click speed are the milestones on the way and the scoreboard, not the finish line. The former 1-in-10,000 programme and grid-cell acceptance goals stay retired; their useful mathematical regressions stay.

Two kernels build one recipe. The exact B-rep carries parametric parts, STEP and drawings. The signed distance field carries lattices, field-driven design and optimisation, and checks the B-rep. Each operation has one written definition that both kernels read.

The kept rules hold: the same question gives the same answer on every machine and at every thread count; no unsafe code; derived geometry is stored by its definition; a proven answer is never wrong; only impossible geometry is an acceptable final refusal.

### 1.1 Long-term goals

Catallax aims for the PRD's parity target: a kernel comparable to OpenCASCADE, SMLib and the Rust kernels. SolidDraft grows toward parametric parts, drawings and assemblies on the SolidWorks side, and implicit modelling and lattices on the nTopology side. These goals guide the build order without setting dates.

## 2. Fornjot's lessons, and what each one changes

Hanno Braun worked on Fornjot from July 2020 to June 2026 and did not reach a usable kernel. His retrospective names six mistakes.

| Lesson | What changes here |
|---|---|
| Early progress does not extrapolate. Primitives and transversal intersections come fast; tangent, coincident and tolerance cases are a cliff. | The hard cases are the work. Refusal repair stays a standing lane, and every slow or refusing case is sorted by cause and gets the fix that matches it. A fixed set of hard cases is measured at every milestone. Progress is counted as target parts that build in the application at click speed. |
| Funding before a product sells a dream. | A works-today list stands beside every sponsorship ask. Aspirations carry an "in design" label. Each milestone publishes what changed, including the numbers that went the wrong way. |
| Incremental change alone cannot fix an architecture. Prototype constantly; integrate only what proved itself. | Prototypes are a standing practice with a fixed shape and a decision rule written before any code. Section 8 says which questions get one and which do not. |
| Half-measures sink a project. In each crisis the options were full commitment or quitting. | Full commitment: SolidDraft drops its other kernel and runs only on Catallax. Milestones end every piece of work as committed or parked, never "almost there". Direction changes only at a milestone, with evidence. |
| A generic kernel is a vaguer target than an application. | Catallax is designed for SolidDraft's intended feature set. Integration moves from the last item of the previous plan to now. |
| Prototyping came too late, and the redesign never landed. The codebase sat half-migrated. | One representation-wide redesign at a time, landed whole. The algebraic foundation landed on 22 October; the tolerance and consistency core in section 3 is the next and only one. |

His final experimental architecture stored geometry as polylines and triangles built at a caller tolerance, and gave up exact geometry on purpose. It never reached booleans, fillets or STEP. Catallax shares his data model (append-only stores, handles, explicit coincidence, keeping the generating definition) and makes the opposite bet on geometry: exact by definition, approximated only on purpose.

## 3. Kernel direction: exact by definition, approximate on purpose

### 3.1 One modelling tolerance

Catallax uses one modelling tolerance of 1 nm when it creates geometry. Two things closer than it are the same thing. The value is decided (section 15).

A floating-point kernel needs a tolerance well above the rounding noise that builds up across operations, which is why Parasolid works at 10 nm. Catallax's numbers carry no accumulated noise, so its tolerance only defines sameness. Mechanical parts are made to microns, so both values sit far below anything a designer means. Proving a gap is under 1 nm costs a few more bits of precision than proving it is under 10 nm. The expensive step is proving an exact zero, and any tolerance removes that step.

Display accuracy and export accuracy are settings of those outputs, not kernel tolerances. The tolerance is chosen before the work starts and is never enlarged because work ran long, so the same question keeps the same answer.

### 3.2 When approximation happens

Approximation never happens for geometry that is clearly apart; the floating-point filter in front of the exact decisions settles those, about four decisions in five on the existing corpus. It never happens for contact known by construction either. Builders record design intent when they create geometry (coaxial, concentric, flush, tangent), and a recorded fact answers its decision with no arithmetic.

It happens in two situations only.

1. A decision is too close to call: the gap is tiny and proving it exactly zero would be expensive. Typical causes are a contact the user meant but reached through arithmetic instead of declaring it, the same point or edge reached by two routes, near-tangency, and imported geometry. The decision is settled as "the same" when the gap is proven smaller than the modelling tolerance.
2. A shape has no exact form the B-rep can carry: fillets, offsets and blends between general surfaces, and intersections of spline surfaces. The B-rep carries an approximation within a local tolerance and keeps the definition beside it.

### 3.3 Four consistency rules

Once anything is approximated, disagreements between approximations, or between an approximation and an exact answer, cannot be ruled out. Sameness within a tolerance is not transitive; an approximate answer can contradict an exact one about the same object; an approximated face may not quite meet its neighbour. Four rules keep the model whole.

1. Decide once, everyone obeys. Each too-close-to-call decision is recorded, and every later decision reads the record instead of asking again.
2. What is declared the same becomes one object. Merging points and edges judged the same stops sameness from drifting across three nearly equal objects.
3. Approximated pieces carry their own recorded local tolerance, the way Parasolid's tolerant edges do. An imported edge does the same. Later operations respect it, and the looseness is visible rather than hidden.
4. Topology is checked before any result is published. A result is a valid solid or a failed operation, never a cracked solid.

In a floating-point kernel every number is an approximation, so disagreements can start anywhere. In Catallax they can start only where the kernel chose to approximate, and each one has a proven size.

### 3.4 Local tolerance caps

| Source of the local tolerance | Target | Cap | Above the cap |
|---|---|---|---|
| An approximation the kernel makes itself: fillets, offsets, blends | 100 nm | 1 µm | The operation refines, splits and refits, or substitutes exact geometry. If none works it fails and reports the deviation it reached. |
| Imported geometry | | 10 µm | Refused, unless the user explicitly accepts a looser tolerance for that file. The acceptance is recorded in the model. |

A local tolerance is also never more than a tenth of the shortest neighbouring edge or the narrowest neighbouring face, so it cannot swallow a feature. The number is never clamped down to the cap; a clamped tolerance leaves a broken solid. Under PRD section 3.2 a failure to meet the cap is unfinished coverage, a bug to fix, not impossible geometry.

The reasons behind the numbers: Parasolid's tolerant-modelling paper (Jackson, 1995, doi:10.1145/218013.218067) says its earlier single tolerance was below what imported mechanical parts carried, ties every limit to feature size, and says "tolerance growth should be contained where possible". ACIS creates a tolerant edge in local operations only below 100 times its base tolerance. Automotive and defence data-quality rules for exchanged models allow face gaps up to 0.01 mm (Hedberg et al., JCISE 16(2), 2016). Fine machining and resin printing hold about ±0.05 mm.

### 3.5 No growth

Parasolid's tolerances grow because later operations build on the approximations and inherit their error. Catallax keeps every shape's definition, so a later operation measures each approximated piece against its definition and resets its local tolerance to the deviation actually present. A boolean only trims an approximated face; a shell or offset of a fillet face is the same kind of shape at another radius. Local tolerances therefore do not accumulate along a feature tree.

### 3.6 Field-only features

Lattices, gyroids, graded density, field-driven thickness and smooth blends belong to the field. They do not have to become B-rep. Conversion is attempted when the user wants it. When it gets too awkward, by size, time or the caps above, the feature stays field-only: it is displayed and exported from the field, as a mesh or 3MF for printing. The B-rep then carries the PRD section 7 stand-in, the outer envelope plus bodies matching mass, centre of gravity and inertia, as decided in section 15. A smooth blend is a field feature in its own right, not an approximate fillet; smooth-minimum blends match a rolling ball only between two perpendicular planes.

### 3.7 The field as the checker

The field holds the exact version of what the B-rep has to approximate. An offset by r is everything within r of the shape. A rolling-ball fillet is everything at distance r from the path of the ball's centre, and that path is where the two faces meet once each is pushed outward by r. Between two planes the path is a line, between a plane and a cylinder a circle or ellipse, and between two cylinders a degree-4 curve. Distance to a curve is an exact expression.

1. Each approximated piece is measured against its own definition: a fillet face against distance r from both faces it rolls on, an offset face against distance r from the original. Proven bounds over the whole face, not samples alone, give the worst deviation, and that value is the piece's local tolerance.
2. The field does not measure against one min/max field for the whole model. After a boolean, a min/max field stops being a true distance near corners and inside the solid, and offsetting it gives wrong shapes. Where the field measures tolerances it uses distance to each face's definition, bounded from both sides. A one-sided "at least this far" bound can prove two things apart but never prove them close.
3. For the whole-part second opinion the field builds every step it can. Where it cannot build a step yet, it takes over the B-rep's result for that step and checks every later step independently. Only that step loses its independent check, and the local measurement in point 1 still covers it. Every check reports what it covered.

When the two kernels disagree, the definition is the truth, not either kernel. A disagreement larger than the local tolerance is a bug in one of them, and the field's distance shows which.

Commercial kernels check results with a self-checker rather than a second kernel; Open CASCADE samples each edge against its faces and compares the gap with the edge tolerance. Hybrid products such as nTopology convert between forms within a stated deviation and do not cross-check afterwards. Catallax's check is stricter, and its fallback is never weaker than theirs.

### 3.8 What we take from SolidCraft, and what we do not

[SolidCraft](https://github.com/tangency-labs/solidcraft) is a parametric CAD application written in Rust on a patched truck kernel. It is after the same thing SolidDraft and Catallax are after: a parametric modeller whose edits come back at click speed. It has one B-rep kernel. Catallax has the B-rep and the field, drives a feature through the field when that is the better route, and then rebuilds the B-rep. SolidCraft's code is offered under MIT or Apache-2.0. This roadmap takes ideas from it and no code; Catallax stays MIT.

What we borrow:

- Analytic routes first. Where a case has an exact construction, such as aligned prisms combined in their sketch plane, inset walls sewn to a rim, or a plane moved by a distance, build that and keep the general machinery for the rest.
- Recovery attempts logged. When a fast route fails and another is tried, the per-operation log records which route answered, which were rejected and why. Nothing is retried in silence.
- Edit cost measured. Construction, the edit, display, requested measurement and refusal are timed separately, with the count of operations rebuilt, so a fast command cannot hide a slow display or a slow refusal.
- Shared evaluations. A point, derivative or curve evaluated once is carried to everyone who needs it instead of being evaluated again. The answer does not change; the work does.
- Case inventories. Each operation keeps a written list of its cases, such as end conditions, neighbour contacts, seams and collapsed boundaries, with a test for each. Coverage is a list, not an impression.

What we reject:

- Nudging or loosening as acceptance. SolidCraft retries a failed boolean by moving the tool, pushing a face or relaxing a tolerance, then accepts the result when its checks pass. Catallax does not. A decision is settled once at 1 nm, the definition stays exact, and no retry changes a settled decision.
- Sampled checks standing in for proof. A fit that is within tolerance at sampled points is measured, not proven. A local tolerance in Catallax is a bound over the whole face.
- Coordinate-grid identity. Rounding solved coordinates to a grid and treating equal rounded values as the same object is a different semantics from a recorded sameness decision. Identity comes from construction origins; fingerprints are checks.

One rule reconciles them: a fast attempt may answer only once it is proved. A cheap analytic route or a reused result is welcome as the first try. Its result is then proved the same way as any other: topology checked, every decision recorded, every approximated piece measured against its definition. If the proof fails, the exact path decides and the attempt is logged. Design intent that geometry cannot prove comes from construction facts, or on import from Jev choosing among options the mathematics has already measured, never from a numeric nudge. A preview is not an answer, so this rule does not slow it down; section 7.3 says how previews are marked.

## 4. The starting point, kept

Carry forward exact arithmetic and certified decisions with deterministic work budgets, definition-based analytic geometry, append-only curved topology, the five-stage boolean engine, bounded mass properties and STEP AP242 export. Preserve the public exchange contract unchanged, including its legacy names. `catallax-features` adds named B-rep recipe stock solids and rationally represented sketch features. The independent `catallax-field` evaluator covers stock solids and booleans; see [FIELD.md](FIELD.md) for its guarantees.

The existing tests include ignored cases. Their source records the reasons; a green default suite does not mean those cases were repaired. Existing work still in the source repository must arrive as a separately reviewed export delta, never by silently replacing this tree.

## 5. The feature build order, kept

The build order merged with the settled recipe definitions stays the backbone of the kernel work. The staged path in section 6 follows it, with two changes: fillets between planes and cylinders and chamfers on plane and circular edges come forward into stage 3. These fillets are exact cylinders and tori. Chamfers produce planes and cones, need no approximation, and SolidDraft loses its current chamfer when its other kernel is removed.

RECIPE.md sections 7.1 to 7.13 now record the settled definitions and product choices, with future variable-radius and nonplanar parting-line extensions kept distinct. Its version 4 specification is docs-only; the crate still reads version 3. Land the schema and migration review separately before emitting version 4 files. Do not relabel arc-length lofts automatically.

1. Finish algebraic sketches, extrusion and revolution, then named holes and the additional extrusion end conditions. Root-defined sketch values, frames, rational-degree angles, elliptic cylinders and analytic trims consume the shared `Real`/`RealPoly` foundation. First-hit caps need certified trimmed-face intersections and finite bounds. Pattern and mirror can follow as soon as their inputs and transformed carriers work; saved feature-pattern expansion needs no separate geometric kernel.
2. Add complete trimmed-face distance queries and regularised offset construction. Build distance shells, outward-shell compositions and thicken from those queries, first on analytic reductions, then on certified pipe and envelope carriers. Keep straight thicken rims and shared origin sets. A sign-only field residual cannot supply these metrics.
3. Build the defined analytic fillet/chamfer runs and planar/cylindrical draft cases. Then implement the agreed roll-over/keep-neighbours modes, corner continuation, consumed-face history and neutral-plane draft retention with kept-region overrides. Do not claim these paths from support equations alone.
4. Land generic derived carriers with their charts, predicates, intersections, integration and export together. Then add rational lofts, fixed spline sweeps and general algebraic offset/blend/draft supports. Analytically reducible line/arc sweeps and their stationary round turns can land earlier. Helices and general rotation-minimising transport need an additional certified function/ODE representation beyond algebraic `Real`.
5. Implement STEP replacement edits and persistent-reference editing alongside import support. Add reference-only open scans and reconstruction deviation reports alongside closed mesh import; never infer a solid from the scan. No edit may infer identity from a distance tolerance alone.

At each stage compare point sets, surviving origins and refusal classifications in both kernels. Require whole-result error transfer proofs under ARCHITECTURE.md before accepting bounded geometry; display agreement is separate evidence. Reassess after shared analytic operations and their comparison proofs land, before committing to general derived carriers.

## 6. The staged path

Stages are ordered by dependency. Work inside a stage runs in parallel unless an item says otherwise.

### Stage 0: land what waited for the foundation

The foundation, the verification exchange, the settled recipe definitions and the Fornjot statement are merged. The lanes told to rebase once onto the foundation now do so and land, or park: refusal fix batch 5, and any other lane in ARCHITECTURE.md's integration order whose work is not yet on `main`. Each rebases exactly once and lands before the next dependent lane starts, or parks.

### Stage 1: cheap speed wins, both verification loops, and Catallax in SolidDraft

- Record design intent when parts are built: coaxial, concentric, flush and tangent. On the existing corpus about one decision in eight proves an exact zero the slow way, while recorded facts settle fewer than one in thirty.
- Compute proven volume, area and centroid only when asked. In the slow cases measured so far, most of the time went to proven measurement, not to the boolean; one block-and-sphere union spent 0.18 s of 12 s in the boolean. The legacy command keeps measuring every answer under its contract.
- Both verification loops from section 7.
- Catallax's application entry: a recipe goes in; the body, its display mesh with edge polylines and face origins, names and typed refusals come out. SolidDraft builds its first parts through it (sketches, extrudes, cuts, revolves and booleans) and removes its other kernel, as its roadmap describes.
- Sound booleans when operands share a face or a face lies in another's plane, in both kernels. SolidDraft's first parts unite flush faces.
- The two prototypes in section 8 start.

Exit: the first target part rebuilds in SolidDraft on Catallax inside the interactive targets.

### Stage 2: the tolerance and consistency core

Section 3 built directly, in pieces, following the migration order already recorded in ARCHITECTURE.md for the composable error policy: certificate types and the modelling tolerance; scoped residual conversions; certificates bound to valid bodies; meeting loci, contact and branch selection; composition through the five boolean stages and every feature node; measurement; export. The decide-once record, merging, local tolerances with their caps, no growth and the topology check before publishing land with the pieces that need them. Each piece is measured by the feature it unlocks, delivered together with its first consumer, without changing any proven answer. If a piece unlocks no feature, work stops and the cause is found before the next piece starts. Near-miss geometry already answers. Operation-budget refusals belong to refusal repair; self-touching results stay correct refusals under RECIPE.md section 1.5 and PRD section 3.2.

Exit: unlocked features meet their recorded tolerance contracts, with the same answer on every machine and at every thread count, and a forged mark is rejected.

### Stage 3: everyday features, and distance in both kernels

- Build order step 1: finish algebraic sketches, extrusion and revolution; named holes and countersinks; the additional extrusion end conditions; pattern and mirror. Exact fillets between planes and cylinders come forward into this stage, ahead of chamfers.
- Build order step 2: trimmed-face distance queries and regularised offsets in both kernels, then shells, outward shells and thicken. In the field this is the per-face distance that section 3.7 uses as the measuring stick.
- The field-only feature track starts once the field's per-face distance exists.

Exit: a part with holes, a chamfer and a pattern, and a shelled part, build and compare in SolidDraft.

### Stage 4: fillets and draft

Build order step 3 for general fillets and draft. General fillets are built on the distance work, owned and constructed as the fillet prototype decided. The field measures each fillet face's local tolerance. Ends and corners follow the settled definitions in RECIPE.md section 4.5.1, last. Planar and cylindrical draft land alongside.

Exit: a filleted part builds in SolidDraft with every fillet face under the 100 nm target, and the field check passes.

### Stage 5: sweeps, lofts and splines

Build order step 4: generic derived carriers with their charts, predicates, intersections, integration and export together, then rational lofts and fixed spline sweeps. Helices need the certified function representation the build order names.

### Stage 6: STEP import and editing

Build order step 5, with imported geometry carrying local tolerances up to 10 µm under section 3.4. Seam curves in STEP need the versioned exchange-contract change. STEP accuracy levels stay an export setting.

### Throughout every stage

- Refusal repair. Every slow or refusing case is sorted into one cause and gets the fix for it: unrecorded intent gets the fact recorded; a near-miss gets the tolerance rule; growing complexity gets a different representation or a bounded approximation; measurement gets moved off the rebuild path; repeated evaluation of the same point, derivative or curve gets one shared evaluation that leaves the answer unchanged. A profile of the case is enough to sort it. Cases a target part hits go first.
- The field grows alongside the B-rep toward parity, and takes over where it cannot build a step yet.
- Each new operation reaches SolidDraft as soon as it lands.
- The works-today list changes in the same pull request as the capability it describes.

## 7. Verification loops

Agents write most of the code, so the loops that check it must be fast and run without a person.

### 7.1 The kernel loop

The `catallax-verify` exchange is the kernel loop: one question or recipe in, the answer with its manifest, deterministic work counts and wall time out. It finds and fixes a bug inside the kernel in seconds. One, two and four workers give identical outcomes and work counts. Its policy, kept from the previous plan:

The policy conflict is closed for the opt-in `catallax-verify` exchange. Whole-request work budgets choose outcomes; wall-clock measurements report performance only. External cancellation leaves an incomplete run. The legacy CLI retains its 50-second deadline unchanged until callers migrate. [VERIFY-INTERFACE.md](VERIFY-INTERFACE.md) defines the versions, adapters and replay driver; [ARCHITECTURE.md](ARCHITECTURE.md#the-time-budget) records the allowance and interactive targets.

Two additions serve the loops: a per-operation log in a stable machine-readable form (how each decision was settled, refusals with codes, marks, local tolerances and timings), and a question format that the application loop can export.

### 7.2 The application loop

Headless SolidDraft, running on Catallax, driven by agents through its command layer and MCP the way a user drives it, at volume. Each run saves the recipe, rendered views, the kernel's per-operation log and a pass or fail result. Any failure exports its recipe or step as a kernel-loop question, so the fix is checked in seconds and the run repeats through the application. It is the same application users run. SolidDraft's roadmap owns this loop; Catallax owns the log and the export format.

Full suites run on GitHub's runners. Timings that count as progress come from a quiet, pinned core, as ARCHITECTURE.md requires; numbers taken under load are labelled as such.

### 7.3 Useful latency, kept

Record rebuild and query time per representative edit, including refusal latency and agreement-check cost. Speed per click matters. Optimize measured costs without weakening proofs. Safe background threads and cancellation must preserve deterministic answers and must never expose a half-built result to the application.

The interactive targets in ARCHITECTURE.md are the click budget: 100 ms for stock solids, 500 ms for booleans, 1,000 ms for recipe feature builds and 2,000 ms for the complete background comparison, including requested measurement and export.

Edit cost is timed in parts: construction, the edit itself, the display mesh, requested measurement and refusal, each beside its deterministic work count and the number of operations rebuilt. A command that returns before its display mesh or its requested comparison has not met the click budget. The budget is checked at the layer where the user waits.

While the B-rep builds an edit, the field may show a preview. A preview is marked as a preview of a named edit, carries no certificate and no names, is dropped as soon as the input changes again, and is replaced by the proven body for the exact value requested. A preview is never published as a result, and a failed build never shows its last preview as success. Previews fix responsiveness; geometric latency is still measured on its own.

## 8. Prototypes

### 8.1 When a prototype is worth it

A prototype is worth it only when all three hold.

1. The answer is unknown.
2. The answer would change what is built or in what order.
3. A throwaway version answers it far more cheaply than building it properly.

Otherwise the work is built directly, with acceptance tests. The tolerance and consistency core is built directly: it is established practice with proofs added, and it is needed whatever a prototype would say. Recording intent and measurement on request are built directly because they are plainly right. The verification loops are tooling.

### 8.2 Prototypes now

1. Fillet construction and ownership. Which kernel builds the fillet the user sees, and does a B-rep fillet face built from the offset-face path meet the 100 nm target inside the click budget? Candidates: the B-rep constructs the rolling-ball surface from its definition and the field measures it; or the field displays the fillet first and the B-rep builds it only for export and drawings. The B-rep candidate builds the strip from the offset-face contacts: the spine is where the two faces meet once each is pushed outward by r, each rail is the foot of the ball on its face, and strip and rails come from one construction so they meet exactly. Measured on the filleted target part: time to first display, deviation from the definition measured by the field, whether the strip and its trimmed neighbours share one boundary curve rather than two nearly equal ones, names after a radius edit, and STEP that opens in OpenCASCADE. Decision rule: the candidate that meets the click budget and the 100 nm target owns fillet display; if both do, the B-rep, because STEP and drawings need it anyway; if neither does, RECIPE.md section 4.5 is revisited for interactive use.
2. Intersections that grow in complexity. The equal-radius cylinder class refuses today on polynomial algebra over very large numbers. Candidates: keep the meeting curve implicit as "where these surfaces meet" and decide against it directly; or carry a certified approximation within the local tolerance and keep the definition beside it. Both candidates keep branch identity across a periodic seam by construction, treat a collapsed boundary such as an apex as its own case, and evaluate each point, derivative or curve once and share it. Measured on the equal-radius and seam cases of the hard-case set, and on translated, scaled and seam-crossing copies of each: time, refusals, topology validity and proven deviation. Decision rule: the candidate that answers the most cases inside the click budget with valid topology, no changed proven answer and the same answer on every copy is promoted; a tie goes to the implicit candidate, because it approximates nothing.

### 8.3 Prototypes later

How large a lattice can be before B-rep conversion is too awkward. The answer sets the threshold at which a lattice stays field-only, and the PRD's open stand-in size (section 11 of the PRD) with it.

### 8.4 How a prototype runs

- A one-page file under `docs/prototypes/<name>.md` comes before any code: the question, the candidates, the test cases, the numbers, the decision rule, and a work bound of one measured round plus at most one repeat. The prototype is decided before the build item it informs starts.
- The code lives on a `prototype/<name>` branch and is never merged. It may break every design rule except two: it may not report a wrong answer as right, and it may not change the measuring tool.
- It is measured with the kernel loop and, once it exists, the application loop, against a baseline taken on `main` at a named revision before anything changes.
- The verdict is promote, park, or repeat once when the measurement itself was wrong. "Promising, needs more time" is park.
- A promoted idea is rebuilt properly on `main`, with the prototype's numbers as its exit check.

## 9. Target parts

Target parts are the scoreboard and a tie-breaker. Work proceeds in the build order whether or not a target part needs it yet; when two items compete for the same hands, the one an earlier unbuilt part needs goes first. P1 to P10 below are confirmed (section 15). Each part is a `catallax-recipe` file under `parts/` with fixed dimensions, so timings compare across revisions.

| Part | What it exercises | First stage that can build it |
|---|---|---|
| P1 Bracket. A plate with a flange along one edge, a boss with a through bore, four bolt holes on a rectangular pattern and a slot through the flange. | Extrusions and booleans, including a union whose operands share a face, and a chain of about a dozen subtractions. | 1 |
| P2 Turned collar. A rectangle revolved a full turn into a ring, chamfered on both outer edges. | Revolution, torus and cone surfaces, chamfer on circular edges. | 3 |
| P3 Filleted boss. A block with a cylindrical boss on top, the junction edge filleted, the block's top edges filleted. | Plane-to-plane and plane-to-cylinder fillets, concave and convex, with transversal ends. | 3 |
| P4 Shelled tray. A block shelled to a thin wall with its top face open, then intersected with a cylinder to round its corners. | Shell with an open face, then a boolean on the shelled body. | 3 |
| P5 Flanged housing. A revolved housing with a pattern of countersunk holes and fillets where the flange meets the body. | Pattern, countersink cones in booleans, fillet against a revolved face. | 4 |

The next five are named so that no near-term decision closes them off: P6 lattice-filled bracket, P7 field-driven bracket with gyroid infill, P8 helical spring and threaded boss, P9 lofted duct, P10 topology-optimised bracket. Their recipes are written when their turn comes.

P0 is the hard-case set: equal-diameter pipes crossing at right angles, cylinders tangent along a line, a cylinder and a sphere with equal radii, a block and sphere union, poses along a seam, near-miss poses, and the general cone and torus poses that refuse today. P0 is both a gauge, measured at every milestone, and a work queue for the refusal-repair lane. Every P0 case also runs as translated, scaled and seam-crossing copies of itself; a copy that answers differently is a gauge failure. Wrong successful results are caught by point-set membership as well as by volume, because a closed body can still be wrong.

A target part counts as built when:

- the B-rep builds it from the recipe with no refusal, and every local tolerance is under its target;
- the field builds the same recipe and the comparison agrees, or the field took over named steps and the report says which;
- a full rebuild finishes inside the interactive targets;
- its STEP file opens in OpenCASCADE with the same face count and a volume within the exported deviation;
- a parameter edit (plate thickness, boss radius, fillet radius) rebuilds with every recipe name resolving to the same face;
- the application loop builds and edits it in SolidDraft and passes.

The scoreboard also carries the head-to-head columns of item 5 in section 12: the same parts built and edited in SolidCraft's released command line and in Catallax.

## 10. Process

### 10.1 Milestones

Milestones are defined only by completed work, with no dates or cadence. A milestone is reached when its must-show list is shown, whenever that is. At each one a note under `docs/milestones/` records what landed, what parked and why, the target-part and P0 numbers, and what the next milestone must show. The same note goes to sponsors.

Every piece of work in flight gets one of two verdicts. Committed means merged, or rebased once and set to land before the next dependent lane starts. Parked means the branch stays, `docs/parked/<name>.md` records what works, what does not and what resuming would take, and no further work happens until a later milestone un-parks it. Work that is "almost there" at two milestones in a row is parked.

A change of kernel direction, meaning a new representation, a new default arithmetic path or a new kernel effort, happens only at a milestone, only with a prototype result or a measured failure as its reason, and only as a written decision in ARCHITECTURE.md.

| Milestone | Must show |
|---|---|
| M1 | Stage 0 lanes landed or parked. Intent recording and measurement on request merged. The kernel loop's per-operation log. `parts/` with P0 to P5 recipes and a scoreboard report. Both prototype files written and both prototypes started. |
| M2 | P1 builds in SolidDraft on Catallax inside the interactive targets, with the other kernel removed. Coplanar booleans sound in both kernels. The application loop runs on GitHub's runners. Both prototypes decided. The first tolerance-core piece merged. |
| M3 | The near-miss class of P0 answers. Holes, end conditions, exact fillets between planes and cylinders, and chamfers in the application; P2 and P3 built. STEP and mesh export from SolidDraft through Catallax. |
| M4 | Pattern and mirror. Offset, shell and thicken in both kernels; P4 built. Names survive the parameter edit on P1, P2 and P4. |
| M5 | General fillets as the prototype decided, with local tolerances measured by the field; P5 built. Planar and cylindrical draft. P0 republished. The field-only track's first feature started. |

An incomplete must-show list means the milestone has not been reached. Progress and parked work remain visible in the ledger.

### 10.2 Landing work without a half-migrated codebase

1. One representation-wide redesign at a time. The foundation was the first. The tolerance and consistency core is the next, and nothing else of that kind starts while it is open.
2. A redesign merges only whole, behind its gate. If the gate cannot be met before the next dependent lane starts, it parks and `main` keeps the previous representation.
3. Each dependent lane rebases exactly once, onto the merged base, and lands before the next dependent lane starts or parks.
4. A lane that needs a type the base does not have asks for a recorded decision. It does not add a parallel type.
5. Independent fixes that widen no interface may land at any time.
6. No decision taken to ship an early target part may close off a later one. The recipe language keeps every operation in RECIPE.md section 4, the field stays a kernel rather than a checker only, and generic swept surfaces stay the planned path for sweeps and lofts.
7. Every merge updates the in-flight ledger in section 13. A branch not in it is parked.

### 10.3 Sponsorship tied to shipped capability

- The works-today list in README.md and on the site names only what `main` does, with a "not yet" line for what it does not.
- A pull request that adds, removes or materially speeds up a capability updates both lists in the same change. Review checks for it.
- The site names an aspiration only with an "in design" label, never as a capability.
- Sponsor links stay live with the works-today list beside them.
- The milestone note is the sponsor update, published when a milestone is reached, with the numbers, including the ones that went the wrong way.

### 10.4 How the work is produced

Most Catallax code is written by AI coding agents, directed and checked by a person. The rules above are written so an orchestrator can dispatch them: a prototype is a file, a branch, a measured round plus at most one repeat and a decision rule; a milestone is a completed must-show list and a review of the ledger; a target part is a recipe and a scoreboard row. Agent throughput also makes direction churn cheap, which is why the milestone rule is strict. Prototypes are a good fit for offloaded workers because their code is throwaway and their decision rule is written down.

## 11. Verdicts

### 11.1 The previous plan's items

| Previous item | Verdict | Reason |
|---|---|---|
| Keep the B-rep starting code | Stays | Section 4, word for word. |
| Next work 1: independent testing and refusal repair | Stays | A standing lane. Each case is sorted by cause and gets the matching fix; cases a target part hits go first. The independent tester and its exchange contract are unchanged. |
| Next work 2: settle RECIPE.md section 7 and widen exact geometry | Stays, largely done | Sections 7.1 to 7.13 are settled. The algebraic foundation widened exact geometry beyond rational frames. Section 7.14's bugs land in stage 0 and stage 3. |
| Next work 3: define and enforce the composable error policy | Changes | Becomes the tolerance and consistency core, stage 2. The 1 nm decision limit stays. The amended certificate covers exact regions and listed local exceptions; decide-once records, merging, local caps, no growth and the publication check govern acceptance. |
| Next work 4: field degenerate zeros and finite closest-point bounds | Stays, moves up | The coplanar part blocks SolidDraft's first parts and lands in stage 1. Two-sided distance to each face's definition is the field's checking role and lands in stage 3. |
| Next work 5: section 7.14 bugs, then holes through mirror | Stays | In the build order's sequence, with exact fillets between planes and cylinders and chamfers brought forward into stage 3. |
| Next work 6: display meshes, STEP import, conversions, hand-off checks, export levels, stand-ins | Changes | Display meshes with edge polylines and face origins stay at the front, because SolidDraft draws them from stage 1. B-rep to field takeover arrives with the field checking role. Field to B-rep conversion follows the field-only rule in section 3.6. STEP import is stage 6. Export accuracy levels stay an export setting. The stand-in is confirmed. Seam curves still need the versioned contract change. Interoperability is measured by opening exports in OpenCASCADE. |
| Next work 7: integrate into SolidDraft once operations are covered | Changes | Integration starts now with the operations Catallax builds today, and SolidDraft runs only on Catallax. Field-only features, printing, FEA queries and interference checks keep the PRD's hand-off guarantees and come on the field track; the FEA solver design stays separate. |
| Feature build order after the shared foundation | Stays | Section 5 as the backbone of stages 3 to 6, with exact fillets between planes and cylinders and chamfers forward. |
| Background comparison prerequisites | Moves | Documentation of `catallax-compare`, moved word for word to COMPARE.md. |
| Deterministic verification policy | Stays | Section 7.1, word for word. |
| Measure useful latency | Stays | Section 7.3. The whole-check measurements of 22 October move word for word to COMPARE.md. |
| Evidence for completion | Stays | Kept at the end, with target parts and both loops added. |

### 11.2 Work in flight and recently landed

| Item | State on 22 October 2026 | Verdict | Condition |
|---|---|---|---|
| Shared algebraic foundation | Merged (#18) | Landed | Waiting lanes rebase once onto it in stage 0. |
| Deterministic verification exchange | Merged (#12) | Landed, becomes the kernel loop | Gains the per-operation log and the exported-question format. |
| Recipe open-item definitions | Merged (#16) | Landed | Each build item implements the settled definition; none is reopened without a recorded decision. |
| Fornjot statement and works-today list | Merged (#15) | Landed | The works-today list follows section 10.3 from now on. |
| Refusal fix batch 5 | Not on `main` | Stays | Rebases once onto the foundation in stage 0. Fixes that keep a P0 or target-part number at least as good land; the rest park with their state recorded. Later batches pick cases by cause. |
| Previous roadmap pull request (#17) | Open | Replaced | This file replaces it; it closes unmerged. |

### 11.3 Queued items

| Item | Verdict | Reason and condition |
|---|---|---|
| Display meshes for new kinds: elliptic cylinders, open conics | Stays | After the foundation, for every kind it added, because SolidDraft draws every face Catallax builds. Edge polylines and face roles come first. |
| Primitive and feature name equivalence | Stays, moves up | SolidDraft's persistent references depend on it, and the target-part edit check tests it. Lands with the application entry in stage 1 and completes when the edit check passes on P1, P2 and P4 at M4. |
| Comparison-found gaps | Stays | The coplanar subtraction query and the same-body boolean assertion land in stage 1. The torus meshing allowance joins refusal repair as a measurement-class case, not parked. |
| Error-policy enforcement | Changes | Becomes stage 2, the tolerance and consistency core. |
| Codex cloud offload | Stays, outside the product roadmap | Capacity, not product. Prototypes suit offloaded workers. It changes nothing in this plan's order. |

## 12. Next items

Numbered in priority order. The project owner picks from this list at each milestone.

1. **Land the stage 0 lanes.** Build. Deliverable: refusal fix batch 5, and any other lane in ARCHITECTURE.md's integration order whose work is not yet on `main`, rebased once onto the foundation and merged, or parked with a state file. Exit: every successful exchange fixture keeps its STEP and manifest content apart from elapsed time; every enabled test passes on GitHub's runners; no lane carries a second representation. Required for M1.
2. **Record design intent at construction.** Build. Deliverable: builders and intersections record coaxial, concentric, flush and tangent facts wherever the construction guarantees them. The record includes the construction behind each face and curve, so a repeated cut into the same target reads its proven contacts instead of deciding them again. Exit: on the existing corpus the share of decisions settled by facts rises and the share proved exactly zero by arithmetic falls; a second cut into a target settles its contacts from the first cut's record with no new exact-zero proof; no proven answer changes; one, two and four workers give identical outcomes. Required for M1.
3. **Measurement on request.** Build. Deliverable: library and application paths build a body without computing its proven volume, area and centroid; measurement runs when asked and returns the same enclosures as before. The legacy command keeps measuring every answer. Exit: the block-and-sphere union builds its body without its integration cost; a measurement request returns identical bounds. Required for M1.
4. **Kernel-loop log and exported questions.** Build. Deliverable: a stable per-operation log in the exchange's output and a question format the application loop writes. Each log entry times construction, edit, display mesh, requested measurement and refusal separately, beside its work count and the number of operations rebuilt. For a boolean, offset, shell or draft it names the route that answered and the routes it rejected, with the reason. Exit: a recipe exported by SolidDraft replays through `catallax-verify` with the same outcome and work count, and a P1 rebuild's log shows the phases as separate entries with every route named. Required for M1.
5. **Target-part corpus and scoreboard.** Build. Deliverable: `parts/` with P0 to P5 recipes and a report that builds each part on both kernels, compares, times and checks STEP. The report has one column per cost: construction, each representative edit, display mesh, requested measurement and refusal latency, beside deterministic work counts, provenance coverage and refusal codes, with the click budgets of section 7.3 checked at the layer where the user waits. Head-to-head: the same parts, written as scripts for SolidCraft's released command line and committed beside the recipes, are built and edited in both tools. The columns are edit time, wrong answers (volume against the analytic value, and point-set membership), failed rebuilds, edit survival (every name resolving after the edit) and STEP validity (opens in OpenCASCADE with the expected face count). The SolidCraft binary is pinned by version and checksum, and timings follow section 7.2. This is how the approach is proved by results rather than by argument. Exit: the report runs on `main`, the first one is committed under `docs/milestones/`, and the head-to-head table covers every part both tools can build. Required for M1.
6. **Application entry for SolidDraft.** Build. Deliverable: the converter from SolidDraft's recipe to `catallax-recipe` following RECIPE.md section 6.2, and a library entry returning the body, its display mesh with edge polylines and face origins, names and typed refusals. A second deliverable is the preview channel of section 7.3: while the B-rep builds an edit, the entry returns the field's view of the recipe, or the previous body, marked as a preview of that edit, so SolidDraft can show a drag at once and deliver intermediate values. The preview carries no certificate and no names, is dropped when the input changes again, and the proven body for the exact value requested replaces it. Exit: P1 built through the entry displays and edits in SolidDraft with faces picked and named; the backend lives in SolidDraft's repository. Preview exit: a drag of P1's plate thickness shows intermediate previews and ends on the proven body for the final value; a failed final build never shows its last preview as success; one undo step returns to the previous proven body. Required for M2.
7. **Coplanar booleans and the same-body rule.** Build. Deliverable: sound union, subtraction and intersection when operands share a face or a face lies in another's plane, in both kernels, with the field's degenerate zero classified. Analytic routes come first: prisms from the same sketch plane with aligned extrusions are combined in that plane and extruded; a batch of tools proven apart from each other cuts in one pass. Each route is a fast path under the one definition in RECIPE.md section 4.4, and its result passes the same topology check and certificate as the five-stage engine; when the check fails, the engine decides. The log names the route that answered and the routes rejected. A cut that lands exactly on a face is a first-class case, and moving the tool is never the answer. Exit: the joined-blocks and coplanar reproducers in COMPARE.md agree with no unresolved observation; subtracting a body from itself reaches the set rule; every P1 boolean gives the same answer with its fast route disabled. Required for M2.
8. **Prototype: fillet construction and ownership.** Prototype, section 8.2 item 1, with the B-rep candidate built from offset-face contacts so the strip and its rails share one construction. Exit: the result section filled in with P3 numbers against a named baseline, the shared-boundary measurement, and a verdict. Required for M2 and decided before the build item it informs starts.
9. **Prototype: intersections that grow in complexity.** Prototype, section 8.2 item 2, measured on the hard cases and their translated, scaled and seam-crossing copies. Exit: the result section filled in with P0 numbers, the same verdict on every copy, and a verdict. Required for M2 and decided before the build item it informs starts.
10. **Tolerance core, first piece.** Build. Deliverable: certificate and evidence types, the modelling tolerance as a recorded setting, the decide-once record, the two reserved refusal codes, and their first consumer: an approximated offset face carrying its definition record. Offsets of planes, cylinders, spheres, cones and tori are exact, a moved plane or a changed radius, and carry no local tolerance; the first consumer is a face with no exact offset. Exit: the face is built and the field's measuring stick proves its local tolerance within the 100 nm target and 1 µm cap, or the operation fails and reports the deviation reached; no proven answer changes; outcomes are identical across worker counts. Required for M2.
11. **Tolerance core, remaining pieces.** Build. Deliverable: merging, local tolerances with the caps, no growth and the topology check before publishing, through the ARCHITECTURE.md migration order, each piece delivered with the feature it unlocks. Periodic charts keep branch identity across the seam by construction, and a collapsed boundary such as an apex or a degenerate edge is its own case, not a tolerance accident. Exit: each unlocked feature meets its recorded tolerance contract without changing proven answers; outcomes are identical across worker counts; a forged mark is rejected; an approximation above its cap fails with its deviation reported; the publication check passes on the P0 cases and their translated, scaled and seam-crossing copies. Required for M3.
12. **Holes, countersinks and end conditions.** Build. Deliverable: RECIPE.md sections 4.1.1, 4.1.2 and 4.3 implemented. Exit: the hole features of P1 and P5 build in SolidDraft inside the interactive targets. Required for M3.
13. **Exact fillets between planes and cylinders.** Build in stage 3. Deliverable: the exact cylinder and torus fillet surfaces of RECIPE.md section 4.5 for plane-plane and plane-cylinder edges, with names. The strip and its two contact rails come from one construction, so the strip's boundaries and the trimmed neighbours' new edges are the same curves. Exit: P3's plane-plane and plane-cylinder fillets build and compare in SolidDraft; each shared boundary is one edge in the topology, not two coincident edges merged afterwards; a radius edit keeps their names, and STEP opens in OpenCASCADE.
14. **Chamfer.** Build. Deliverable: RECIPE.md section 4.6 for plane-plane and circular edges, with names, and a written case inventory: equal-distance, two-distance and distance-angle chamfers; straight and circular edges and chains of them; and the end conditions where a chamfer meets a free end, a neighbouring chamfer or fillet, or a face too small for the requested distance. Each case has a test. The planes and cones are constructed from exact inputs, never recognised from sampled points. When the requested distance does not fit, the refusal names the largest distance that does; the size is never changed in silence. Consumed faces follow the agreed roll-over and keep-neighbours modes with consumed-face history. Exit: P2 builds and compares; every inventory case builds or refuses with its reason; STEP opens in OpenCASCADE. Required for M3.
15. **Name equivalence and edit survival.** Build. Deliverable: primitive and feature names equivalent across both kernels and the parameter-edit check in the scoreboard. Exit: the edit check passes on P1, P2 and P4. Required for M4.
16. **Pattern and mirror.** Build. Deliverable: RECIPE.md sections 4.13 and 4.14 with saved feature-pattern expansion. Exit: P5's hole pattern builds with names resolving per copy. Required for M4.
17. **Distance, offset, shell and thicken.** Build. Deliverable: build order step 2 in both kernels, with the field's two-sided per-face distance. The analytic paths are separate and each is proven: a convex planar shell as an intersection of half-spaces; face movement, planes by distance and cylinders, spheres, cones and tori by radius, with their meetings rebuilt exactly; a shell of a revolved body by sewing exact inset walls to its rims rather than asking a boolean to find them; and the general face by the distance definition with a fitted carrier beside it. Each path is a fast route under the one definition in RECIPE.md sections 4.7 to 4.9. Its result must equal the distance-defined result, which the field's two-sided distance proves; a face moved to the wrong place is a wrong answer, not a tolerance. Changes of topology follow the definition: opposing walls that consume the cavity can be a valid solid, and corners, openings and rims are a written case list with a test each. Exit: P4 builds and compares; the field reports an offset face's deviation from its definition; every analytic path agrees with the distance definition on the field's measure, and the log names the path that answered. Required for M4.
18. **General fillets.** Build in stage 4, as the prototype decided. Deliverable: general fillets under RECIPE.md section 4.5, then ends and corners per section 4.5.1. The construction follows item 13: strip and rails from the offset-face contacts, with a fitted carrier beside its definition where no exact surface exists. A written case inventory with a test each: a rail crossing a periodic seam, an open-chain end, a curved neighbour, a round meeting an earlier round, a rim nearly swallowed by the radius, equal-radius rounds meeting at a mitre, and corner patches. The field proves the carrier's local tolerance over the whole face; a fit checked at sampled stations is a development aid, never the acceptance. Exit: P3 and P5 build with every fillet face under 100 nm, every inventory case builds within its tolerance or refuses with the deviation reached, a radius edit keeps every name, and STEP opens in OpenCASCADE. Required for M5.
19. **Draft.** Build. Deliverable: planar and cylindrical cases of RECIPE.md section 4.10: a planar wall rotates about its line in the neutral plane, a cylindrical wall becomes a cone. A neighbouring cap is checked where the wall meets it, not by refusing unrelated faces. Swallowed faces and retained regions follow section 4.10.1. A drafted boss with a hole and a hollow section are inventory cases. Exit: a drafted variant of P1 builds and compares; each inventory case builds or refuses with its reason. Required for M5.
20. **First field-only feature.** Build. Deliverable: lattice inside a shell per RECIPE.md section 4.15, displayed and exported from the field. Exit: P6 exports a printable mesh, and the B-rep carries the PRD section 7 stand-in, the outer envelope plus bodies matching mass, centre of gravity and inertia. Starts for M5; completes when this exit is shown.
21. **Milestone notes.** Build, recurring. Deliverable: a note under `docs/milestones/` when each milestone is reached. Exit: the note records the completed must-show list, and the same text went to sponsors.

Sweeps, lofts and splines, then STEP import and editing, follow these in the build order.

## 13. In-flight ledger

Updated at every merge and every milestone. A branch not listed here is parked.

| Lane | State on 22 October 2026 | Landing condition |
|---|---|---|
| Refusal fix batch 5 | Committed, rebase once onto the foundation | Lands before the next dependent lane starts, or parks; required for M1. |
| This roadmap | Owner's six answers recorded in section 15 | Merges with those answers applied; required for M1. |

## 14. Proposed PRD changes

These change the PRD's direction. Item 3 was approved on 23 October and is applied in PRD section 3.1 and the amended composable-error ruling in ARCHITECTURE.md. The modelling tolerance and stand-in decisions are also settled as recorded below. The remaining proposals still need approval.

1. **Section 1.** "It is the kernel for SolidDraft and for anyone else who wants it" becomes "It is the kernel for SolidDraft, designed for the features SolidDraft's users drive as it grows toward the most usable of SolidWorks's parametric modelling and nTopology's implicit modelling. Anyone may use it; generalising it for others comes after it carries one application well."
2. **Section 3.8.** The target stays: comparable to OpenCASCADE, SMLib and the Rust kernels. Add: "Target parts that rebuild in SolidDraft at click speed with no wrong answer are the milestones toward it and the published scoreboard."
3. **Section 3.1, decided on 23 October 2026.** Adopt the tolerance model of section 3: 1 nm for sameness and decisions; approximation only in the two stated cases; the four consistency rules; local tolerances of 100 nm target and 1 µm hard cap for kernel-made approximations and up to 10 µm for imports with an explicit recorded acceptance above that; the tenth-of-the-nearest-feature cap; no growth. The owner approved carrier A, a certified polynomial or rational B-spline surface beside its exact definition, and the case-1 widening that permits it inside the B-rep. The whole-result certificate proves exact agreement except for its listed approximated/imported faces and bounded decisions, each within its applicable local limit, with complete coverage and valid topology. A fitted face within its cap is not also required to meet a universal 1 nm whole-result bound. The amendment covers measurement, comparison and STEP hand-offs; their declared output errors remain separate from model tolerances.
4. **Section 3.2.** Keep the rule. Add: refusal latency is reported with every refusal; a refusal that arrives after the click budget is a bug twice; failing to meet a local tolerance cap is unfinished coverage, not impossible geometry.
5. **Section 3.5.** Add the field's checking role: each approximated piece is measured against its own definition with two-sided bounds; for the whole-part check the field takes over the B-rep's result for a step it cannot build yet, and every check reports its coverage.
6. **Sections 3.3 and 7.** Field-only features need not become B-rep. Conversion is attempted on request and abandoned when it is too awkward; the feature then stays field-only. Decided: the B-rep carries the PRD section 7 stand-in, the outer envelope plus bodies matching mass, centre of gravity and inertia.
7. **Section 9.4 and RECIPE.md section 6.** SolidDraft runs only on Catallax. Its kernel slot stops being a swappable slot; its interface and recipe format converge on RECIPE.md through the section 6.2 mapping.
8. **Section 10.** Integration with SolidDraft moves from step 6 to step 2. The field kernel grows in parallel with the B-rep toward parity instead of gating each widening step.
9. **Section 11.** Close the deadline item: the verification exchange's deterministic allowance decides outcomes, and the legacy command keeps its deadline for compatibility. Record the decided 1 nm tolerance and field-only stand-in.

## 15. Decided on 22 October 2026

1. **Modelling tolerance.** 1 nm.
2. **The stand-in for field-only lattices.** Yes. The B-rep carries the PRD section 7 stand-in, the outer envelope plus bodies matching mass, centre of gravity and inertia, so drawings and STEP show where a field-only lattice is.
3. **Target parts.** P1 to P5 and P6 to P10 confirmed.
4. **Milestones, not a checkpoint cadence.** No time-based checkpoints or schedule. Milestones depend only on completed work.
5. **Sponsorship during the change.** Keep sponsor links live with the works-today list beside them.
6. **Earlier fillets.** Fillets between planes and cylinders, which are exact surfaces, come forward into stage 3, ahead of chamfer. General fillets stay in stage 4 as the fillet prototype decides.

## Evidence for completion

Use mathematical cases, independent numeric oracles, scale and placement variations, repeatability checks and measured timings. Keep refusals, marked answers and wrong proofs distinct. An independent tester built without seeing the kernel's code uses the documented exchange contract. No unpublished tester design or acceptance system forms part of this plan.

Add to that evidence the target-part scoreboard, the P0 numbers at each milestone, and passing runs of both verification loops.

---

