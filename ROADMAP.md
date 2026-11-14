# Roadmap

One roadmap for Catallax, the kernel, and SolidDraft, the app. They are one project in two repositories.
`VISION.md` in each repository says why. This file says what is done, where we are and what comes next.
Dates appear only on finished work. Nothing ahead has a date.

## How we run it

- **Most critical first.** The order is set by what hurts most against OpenCASCADE and by what unlocks the most other work. Speed is the last phase.
- Milestones are numbered by phase in the order of work: 1.1 is the first milestone of phase 1. Work on different code runs side by side. When machines or people are short, the lower number wins.
- Every phase ends at a checkpoint. A checkpoint is a real step with four parts:
  1. No new work starts. Running work merges, or is parked with a note saying what works and what does not.
  2. We gather the results: the coverage count, the edit-chain scoreboard, the STEP census, the list of what still refuses and why, prototype verdicts.
  3. We update this file. Finished items move to the history with their date.
  4. The owner picks the order of the next phase.
- Direction changes at checkpoints. Between them we follow this file. One exception: a proven answer found wrong stops work on that method at once.
- A milestone is done when its "done when" line is shown on real parts, through SolidDraft where the app has the feature.
- Earlier plans named milestones M1 to M35. [The table at the end](#old-milestone-names) maps each old name to its number here, so old notes still resolve.
- Jobs marked * need a machine to themselves. They are timing and GPU jobs.

The edit-chain scoreboard is the measure `VISION.md` asks for. The target parts P1 to P10 and ten real parts from the public STEP census are built and edited step by step in Catallax, OpenCASCADE and SMLib, with one definition of done for all three.

## Where we are

In release: Catallax 0.2.0 with SolidDraft 0.5.0, cut on 14 Nov 2026. The tags are not published yet.

The kernel builds exact solids from a recipe. It has stock solids, extrude, revolve, pattern and mirror, Booleans on analytic solids, certified fillet families, chamfers, planar shells, planar and cylindrical draft, exact lofts, and sweeps along lines, arcs and one spline. It proves where spline faces meet, cuts lofts from the side, and reads and writes STEP and IGES. It imports closed meshes, turns a scan into a clean part, reports mass with proved bounds and answers clash or clearance. A second kernel, the distance field, previews the same recipe and owns lattices and gyroids.

The app has the Edit Suite frame, a full 2D sketcher, the part features above, 20 mate types, drawings, proof lights, and STEP, IGES, STL, 3MF and OBJ export. It runs on Linux, on macOS and in the browser. Assemblies are desktop only.

### Coverage against OpenCASCADE and SMLib

100% is everything the two reference kernels model, put together as 61 capability lines. A line we have in full counts 1. A line we have in part counts half. The count says what exists, not how fast it is or how many real parts it finishes.

| | 10 Nov, previous release | 12 Nov | 14 Nov, this release |
|---|---|---|---|
| Combined scope | 48% | 52% | 55% |
| Of OpenCASCADE's 58 lines | 50% | 55% | 58% |
| Of SMLib's 52 lines | 48% | 54% | 57% |
| Lines complete, in part, not started | 15, 28, 18 | 17, 30, 14 | 17, 33, 11 |

| Area | Lines | Now | Not started | Closed by |
|---|---|---|---|---|
| Curves and surfaces | 12 | 54% | Surface fitting. Gordon, Coons and N-sided surfaces. Bend and twist. | 2.1, 3.4 |
| Intersection | 6 | 58% | None. Tangent contact and some cone and torus pairs refuse. | 1.1, 1.3 |
| Solids and topology | 7 | 64% | Surface, wire and non-manifold bodies. | 2.2, 3.2 |
| Features that add or remove material | 11 | 55% | Extrude with taper. Sweeps with guides, twist or a helix. | 1.3, 2.3 |
| Edge and wall features | 9 | 39% | Curvature-continuous blends. Remove a face, fillet or hole. | 1.4, 1.5, 2.4, 3.4 |
| File exchange | 8 | 62% | Assemblies, colours, names and annotations in files. | 1.2, 3.2 |
| Meshes and queries | 8 | 56% | Mesh Booleans and decimation. Hidden-line views. | 3.1, 3.3, 3.4 |

### The three gaps that hurt most

These set the order of phase 1.

1. **Real STEP files do not open.** Of 3,000 public models, 6% import in full. Of the 728 stopped at spline faces, one imports. In the largest group the file is off by less than 0.006 mm, and our exact rules reject it. [The census](docs/research/step-census.md) has the numbers.
2. **Booleans on spline faces do not reach a file.** Side cuts on a loft build. All 14 recorded spline Boolean requests still refuse before a STEP file is written.
3. **Fillets, shells and draft work on certified families only.** Edge and wall features are the weakest area at 39%, and no line in it is complete.

A work unit counts exact arithmetic and is the same on every machine. The allowance per edit is 20 million units, about 2 to 4 seconds of proof work on one core.

Twelve things neither reference kernel does are not in the count. They include a proof behind each result, exact numbers, the same answer on every machine, the field kernel, lattices and scan to part.

## Plan

### Phase 0: close the release

**Doing:** making the 0.2.0 and 0.5.0 cut clean. **Why first:** everything after starts from a working release. **Unlocks:** a fixed point to measure from.

- Both tags published. SolidDraft pinned to the Catallax tag.
- Every download installs, builds a part and exports it. The browser app loads, opens, edits and exports.
- Release notes list the slow steps and the two export timeouts, with working links.
- Faults found at the cut: the macOS test panic is fixed on main. The lattice and gyroid export that timed out is rerun on a quiet machine.
- Docs match the cut. This file replaces the old plan and roadmap files.
- Three quick re-checks: whether the later cuts on the filleted P5 housing already build (part of 1.3), the last 72,166 units on P4's scoreboard run, and P1 with its fast route off.

**Checkpoint 0.** This roadmap is adopted. Coverage is 55%. Review the 13 pending PRD proposals with the owner (section 14 of [the archived plan](docs/archive/PLAN-2026-11.md#14-proposed-prd-changes)).

### Phase 1: close the three biggest gaps against OpenCASCADE

**1.1 Numerical modelling under stated tolerances**
- Doing: where an exact proof is not available, the kernel builds the geometry numerically within a stated tolerance. Each result records its tolerance, the checks it passed and any repair. Its STEP file is written and marked approximate.
- Why here: it is the base under all three gaps. Real files miss exact rules by micrometres. Spline Booleans stop on a gap or a coincidence they cannot prove. OpenCASCADE and SMLib work to a tolerance throughout. We will work to one only where the exact proof is missing.
- Unlocks: milestones 1.2, 1.3 and 1.5 below, and the first bounded-answer rule, for [tangent contact](docs/research/tangent-touch.md).
- Done when: one curved housing builds as a chain of loft, spline cut, fillet and one more edit. Its STEP opens in OpenCASCADE with healing off and takes an edit there.

**1.2 Real-file import**
- Doing: opening real STEP files under those tolerances, analytic files first, then spline files. What still refuses is counted by cause.
- Why here: it is the largest single gap, and it is the first thing a user meets: a bought-in part that will not open.
- Unlocks: bought-in parts in every later feature, the ten real parts on the scoreboard, direct edits (2.4) and imports in the field (4.4).
- Done when: the 3,000-file census is rerun and the share that imports is published beside the 6%. Replacing an imported file re-points every reference.

**1.3 Cuts and joins on spline faces, with what is left of the SMLib intersection port**
- Doing: Booleans on lofts and imported freeform parts that finish and export.
- Why here: cuts and joins shape every part after its first feature, and today they stop at the first freeform face. It comes after 1.1 because 4 of the 14 refusals are unproved gaps and coincidences, which the rules of 1.1 are meant to settle.
- Unlocks: every freeform part past its first cut, and the fillets of milestone 1.5.
- Done when: the 14 recorded spline Boolean requests and the loft top cut build and reach STEP. The later cuts on the filleted housing build.

**1.4 Everyday edge and wall features**
- Doing: fillets over mixed corners, shells with side openings, hollow oblique draft, chamfer corners and general holes on ordinary prismatic and turned parts.
- Why here: two thirds of the census files have only flat and round faces, so this is where most real parts are. It touches different code from 1.1 to 1.3 and runs beside them. The owner's bar from 7 Nov stands: at least match OpenCASCADE on fillets.
- Unlocks: the five everyday target parts end to end, and the first complete line in the weakest area.
- Done when: P1 to P5 each build, edit, save, reopen and export through SolidDraft.

**1.5 Fillets between spline faces**
- Doing: rounded edges on lofted and imported freeform parts, across more than one surface patch.
- Why here: it needs milestone 1.3. The centre curve of a spline fillet is found by the same engine that intersects spline faces. [One loft edge already rounds within 59 nm](docs/research/spline-fillet-offset-search.md).
- Unlocks: finished freeform parts. P9, the lofted duct.
- Done when: a fillet across several patches builds on a loft and on an imported part, and a radius edit keeps the names.

**App, alongside:** browser first. Kernel work runs off the interface thread, commands disabled in the browser work there, and a part built within tolerance shows as finished. Why here: the browser is where a new user tries a real file.

**Checkpoint 1.** First edit-chain scoreboard against OpenCASCADE and SMLib. Census rerun. Order of phase 2.

### Phase 2: the modelling lines we do not have yet

Why this phase is second: these are whole lines of the count with nothing built. They matter less than phase 1 because a user can work around a missing surface tool, but not around a file that will not open or an edge that will not round. The order inside is the owner's order of 13 Nov.

| # | Milestone | What we are doing | Why here | Unlocks | Done when |
|---|---|---|---|---|---|
| 2.1 | NURBS toolkit | Freeform surface construction: Coons patch, surface through a grid of points, Gordon surface, guided skin, N-sided fill. | First of the new milestones in the owner's order. Three "not started" lines sit here. Surface bodies and direct edits need its tools. | Freeform surfacing. Lofts with guide curves. | Each surface builds and exports. Knot insertion, degree raise, split and join are public. First step: a dated decision that widens the 10 Nov port rule. |
| 2.2 | Surface bodies | Open surfaces and wires as bodies. Sew, trim, extend, thicken into a solid. | Second in the owner's order. Direct edits (2.4) need its trim, extend and sew. | Surface modelling. Imported files that are not closed solids. | A surface trims, extends and thickens into a solid. Faces sew into a solid outside import. |
| 2.3 | General sweeps | Extrude with taper. Sweep with a guide, with twist and along a helix. Extrude up to a curved face. | Two "not started" lines. It does not depend on 2.1 or 2.2, so it runs beside them. | Springs, threads, tapered and twisted shapes. P8, the spring and threaded boss. | P8 builds. |
| 2.4 | Direct edits | Delete a face, a fillet or a hole and heal the solid. Move, offset or replace a face. | Needs 2.1 and 2.2. It pays only once imports work (1.2). | Editing bought-in parts that have no feature tree. | Each edit works on an imported solid. |

**App, alongside:** import health report, reverse-engineering workspace, direct-edit commands, reference axis and point, slot, draft check, exact sketch solver.

**Checkpoint 2.** Scoreboard and census rerun. Order of phase 3.

### Phase 3: drawings, assemblies, exchange and the last lines

Why this phase is third: it completes the count. Each item has a working stand-in today (mesh outlines in drawings, assemblies in the app), so parts get made without it.

| # | Milestone | What we are doing | Why here | Unlocks | Done when |
|---|---|---|---|---|---|
| 3.1 | Drawing views from the kernel | Exact section curves, silhouettes on curved solids, hidden-line views. | Drawings work today from meshes. This makes them exact. | Drawings that match the model to the last digit. | SolidDraft drawings use them in place of mesh outlines. |
| 3.2 | Assemblies and exchange | The kernel holds several solids as an assembly. STEP carries structure, colours, names and annotations. | Needs import (1.2) to be solid first. | Whole products in one file, both ways. | A STEP assembly round-trips. IGES exports spheres and tori. |
| 3.3 | Clash, clearance and mass | Clash, clear by a stated gap, or unresolved, for every pair. | Exact overlap already works. Distant and nested pairs are left. | A proven gap between parts. Refusals that state the largest radius or thickness that fits. | Every pair in a test assembly gets one of the three answers. |
| 3.4 | Tools only SMLib has | Bend and twist. Curvature-continuous blends. Mesh Booleans and decimation. | Last because OpenCASCADE does not have them either. | The final three lines. | Each has a first build. |

**App, alongside:** assemblies in the browser, assembly drawings, BOM and balloons, Pack and Go, the remaining mates, drawing annotations, URDF and USD export with joints, the refusal assistant.

**Checkpoint 3.** Every one of the 61 lines has something built. We publish what still refuses and why. The owner decides what to resolve before moving on.

### Phase 4: field, lattices, printing and scans

Why this phase is fourth: this is what neither reference kernel has. It is our edge, not a gap, so it waits until the gaps are closed. It never blocks a B-rep feature in the meantime.

| # | Milestone | What we are doing | Why here | Unlocks | Done when |
|---|---|---|---|---|---|
| 4.1 | Field preview and first display * | The app draws from whichever kernel is faster, with a proved, pending, bounded or refused light. | First in the phase: everything else in it is seen through this. | A feature shows at once. | A preview never stands in for an export. |
| 4.2 | The field catches up | Field versions of fillets, sweeps and lofts. | The field must hold every part before it can print or analyse it. | Every part has both descriptions. | P5 builds in the field with no fillet refusal. |
| 4.3 | Lattices and printing | Printable lattice parts. The missing piece is [a proof that a mesh keeps its cavities](docs/research/adaptive-printing.md). | Needs 4.2. | P6 and P7. Sheet gyroids. | P6 exports a printable mesh at its recorded accuracy. |
| 4.4 | Imports and scans in both kernels | Imported spline parts in the field. [Scan fitting](docs/research/scan-fitting.md) grows past planes and cylinders. | Needs import (1.2) and the spline work (1.3, 1.5). | Bought-in parts and real scans in the field. | An imported spline part builds in the field. A real scan becomes a clean part. |
| 4.5 | FEA groundwork, first result and P10 | [Parked by the owner on 13 Nov](docs/research/fea-first-beam.md). | Starts on the owner's word at a checkpoint. | Stress on the same part the B-rep exports. | |

**App, alongside:** field-driven thickness, more grades, fields on imported bodies, infill that is denser near the skin or the load.

**Checkpoint 4.** Order of the speed work, from measurements.

### Phase 5: speed and simplification *

**5.1 Speed and simplification.**

Why last: the owner's ruling of 11 Nov. Clock targets were stalling capability work, and speeding up a method that is about to be replaced is wasted.

- P1 to P5 and their edits meet the click targets in SolidDraft: 100 ms for stock, 500 ms for a Boolean, 1,000 ms for a feature, 2,000 ms for proof.
- Cheaper proofs, in the order [the P1 work](docs/research/cheaper-proof-p1.md) ranked them: proofs from construction facts, prepared bounds, reuse across edits, coarse proofs refined where they fail.
- The GPU and every core carry what they can. [Booleans gain nothing from extra cores today](docs/research/run-ahead.md).
- Simplify. Routes that the measurements show are not needed are removed.
- Slow steps recorded at the 0.5.0 cut: loft 35 s, fillet on a loft edge 26 s, spline sweep 13 s, P5 fillet edit 4 to 6 s, field mass 3.5 to 7.7 s.
- Windows packages and signed downloads.

**Final checkpoint.**

### Lanes that never stop

- The resolve list. Every refusal of possible geometry is sorted by cause and fixed in batches. Cases a target part hits go first.
- Each new kernel operation reaches SolidDraft as it lands, in the browser first.
- Generated recipes with built-in checks run nightly and publish counts by outcome.
- The SMLib port supplies the units this roadmap needs next.

## Rules that stay fixed between checkpoints

- The allowance is 20M.
- Proven, bounded and refused answers stay separate.
- Kernel changes are general. No allowance or tolerance is loosened to make one case pass.
- Clock time never fails a change before phase 5.
- The app shows Coming soon only for scope written on this roadmap.

## History

### What was added

| Date | Added |
|---|---|
| 22 Oct 2026 | Catallax starts after the Fornjot retrospective. Exact numbers, one recipe, a B-rep kernel and a field kernel, work counted in deterministic units. |
| 25 to 29 Oct | Catallax 0.1.0 with target parts P0 to P5 and a scoreboard. IGES export. SolidDraft 0.1.0 runs on Catallax only. |
| 1 to 4 Nov | Draft past holes and bosses. Fillets between planes and cylinders, ring and loop fillets. Sweeps along lines and arcs. Planar shells. Clean analytic STEP import with repairs. Sketch editing tools, autosave, the website. Catallax 0.1.2 and 0.1.3, SolidDraft 0.2.0 and 0.3.0. |
| 5 to 9 Nov | Fillet campaign: corners, linear variable radius, cone bands and sphere rings. Exact spline import. Ellipse and spline profiles. Sketch relations, trims and text. 20 mate types, interference and clearance. Scene import. Catallax 0.1.4. |
| 10 Nov | Exact rational lofts. Names checked to survive edits. Catallax 0.1.5, SolidDraft 0.4.0. Coverage 48%. |
| 11 Nov | Three-circle lofts inside 20M. Sweep along a spline. Proved cuts of spline faces by planes, cylinders, spheres, cones and tori. Closed mesh import in both kernels. First clash and mass queries. GPU viewer. Census of 3,000 public STEP files. |
| 12 Nov | Spline faces meeting spline faces, by the 2026 method with exact proof. SMLib intersection code ported as a source of candidates. IGES import. Scan to clean part. First fillets between spline faces and on a loft edge. Coverage 52%. |
| 13 Nov | Exact sections and first silhouettes. Side Booleans on lofts. Files with several solids. P1 inside 20M, down from 41.6M to 15.3M units. The bounded answer kind. In the app: Edit Suite frame, Loft, spline sweep, Offset, Thicken, feature patterns, proof lights, STEP proof reports, scan fit, exact clash. |
| 14 Nov | Several bodies and Combine. Tolerance exit run recorded. `VISION.md` in both repositories. Release cut: Catallax 0.2.0, SolidDraft 0.5.0. Coverage 55%. |

### Decisions

| Date | Decision | Reason |
|---|---|---|
| 22 Oct | Geometry is exact by definition. One modelling tolerance of 1 nm. | One rule for sameness, and no error growth across edits. |
| 23 Oct | Fitted faces are allowed with local caps, 100 nm target and 1 µm limit. | A 1 nm bound on the whole part was out of reach for fitted faces. |
| 26 Oct | Display and proof have separate budgets. Field comparison runs on request only. | What the user sees must not wait for proof. |
| 7 Nov | Fillets must at least match OpenCASCADE. | The owner's bar for the fillet campaign. |
| 8 Nov | The fillet panel asks for edges and radius only. | Users should not have to pick kernel constructions. |
| 10 Nov | NURBS first. SMLib code is adapted only for spline intersections, spline blends and fitting. | There is no other efficient way in those areas. |
| 10 Nov | Both kernels act on the same geometry, imports included. The faster one does the job. The B-rep owns export. | Use each kernel where it is strong. |
| 11 Nov | Capability first, speed last, in one final speed milestone (now 5.1). | Clock targets were stalling capability work. |
| 11 Nov | Use every core and the GPU by design, with the same answer at any thread count. | Do not gate parallel work on speed studies. |
| 11 Nov | Spline against spline uses the 2026 method with an exact proof layer. | Fast search, with Catallax's own checks deciding. |
| 13 Nov | Fewer refusals. Tolerance levels per kind of case, and a bounded answer that is never shown as proved. | A kernel that refuses possible geometry is unfinished. |
| 13 Nov | New milestones in this order: NURBS toolkit, surface bodies, then drawing views and direct edits. | Owner's order on the release board. |
| 13 Nov | Kernel work must be general. When a general method cannot meet a strict requirement, choose an exception. | Tuning to one sample part is worthless. |
| 13 Nov | FEA parked. | Owner's ruling. |
| 14 Nov | Better algorithms. No proof-job machinery. | The point of a kernel is to create the geometry. |
| 14 Nov | Numerical B-reps under stated tolerances where exact construction is unavailable. Each answer claims what its checks establish. Approximate STEP is written and marked. | Finish the part, honestly. |
| 14 Nov | One roadmap with checkpoints. Most critical first: the biggest gaps against OpenCASCADE lead. The allowance stays 20M. Speed and GPU work last. | Stop changing direction between checkpoints, and work on what hurts most. |

### Prototypes and research

Each entry links to a short file: what we were trying to achieve, why it was proposed, what was tried, the result, and why it was adopted, parked or dropped.

| Work | Outcome | Result in one line |
|---|---|---|
| [Fillet prototypes: which kernel builds it, and patch-first](docs/research/fillet-prototypes.md) | Parked | No solid came out of either round. Certified fillet families replaced them. |
| [Hard-intersection prototype](docs/research/hard-intersections.md) | Parked | The existing code answered 15 of 56 hard inputs. The two candidates answered 10 and 7. |
| [Tolerance as a speed lever](docs/research/tolerance-as-speed-lever.md) | Dropped | Exact-zero proofs are 0.4% to 12.4% of the work, so a tolerance does not make slow cases fit. |
| [Decide-once merging](docs/research/decide-once-merging.md) | Dropped | No saved answer would have improved. |
| [Threads inside one field request](docs/research/field-threads.md) | Adopted | Same bytes and work counts on 1 to 12 processors. |
| [Running ahead on other cores](docs/research/run-ahead.md) | Adopted | STEP export about 8 times quicker on 24 threads. Booleans no quicker. |
| [SMLib intersection port scope](docs/research/smlib-intersection-port-scope.md) | Adopted | 37,255 lines needed before a plane cut. SMLib's tracer missed a tiny loop, so Catallax keeps its own complete search. |
| [Spline against spline by the 2026 method](docs/research/spline-intersection-2026-method.md) | Adopted with a proof per cell | All 24 test pairs resolve inside 20M: 20 proved, 4 tangent contacts refused by name. |
| [Spline fillets: the offset search](docs/research/spline-fillet-offset-search.md) | Route changed | The direct search hit its depth limit. The replacement rounds one loft edge within 59 nm. |
| [Tangent-touch prototype](docs/research/tangent-touch.md) | Parked | Interior touches classify. A touch on a border is undecided. First candidate for a bounded answer. |
| [Census of 3,000 public STEP files](docs/research/step-census.md) | Recorded | 6% import in full. The largest group fails on gaps under 0.006 mm in the source files. |
| [Cheaper proof on P1](docs/research/cheaper-proof-p1.md) | Adopted | 41.6M to 15.3M units by proofs from construction and reuse inside one request. |
| [Proof-cost research](docs/research/proof-cost.md) | Recorded | No fixed allowance covers every input. Sampling alone is unsound. A GPU cannot rescue a work refusal. |
| [Independent body checker audit](docs/research/body-checker-audit.md) | Closed | Checking can cost up to 25 times the build. The owner chose not to build certificate systems. |
| [Adaptive printing](docs/research/adaptive-printing.md) | Parked in parts | Works on a convex block. A distance bound cannot prove that cavities survive. |
| [Scan fitting](docs/research/scan-fitting.md) | First cut | 30,000 noisy points from P3 come back as its planes and cylinders. The torus round does not. |
| [First FEA beam](docs/research/fea-first-beam.md) | Parked | The cantilever runs within 0.9% of beam theory, unproved. |

## Closed or not scheduled

- Independent body checker and the proof verifier idea. Closed by the 14 Nov choice not to build certificate systems.
- Durable proof jobs and the bounded-modelling contract spec. Rejected on 14 Nov.
- Decide-once merging, tolerance as a speed lever and the two fillet prototypes. Dropped or parked on measurement, see the table above.
- Dual-implicit experiments and the old Draftwell queue. Closed on 12 Nov.
- Not scheduled: machine-checked core arithmetic, the field choosing which shell faces survive, USD BRep slices, a pure field kernel as its own product.
- Waiting on the owner: Windows signing, crates.io publication, when FEA restarts.

## Further reading

Existing: `VISION.md` in both repositories, `docs/ARCHITECTURE.md` for dated decisions, `docs/CAPABILITIES.md` for what works today, `docs/TARGET-PARTS.md` for P0 to P10 and when a part counts as built, `docs/milestones/` for milestone notes and the STEP census, `docs/research/` for the files linked above, and `docs/archive/PLAN-2026-11.md` for the plan this file replaces.

To be written, each short:

- `docs/COVERAGE.md`: the 61 lines one by one, the scoring method and the state of each line.
- `docs/LEFTOVERS-0.2.0.md`: everything the release sprint postponed, with its recommendation.

### Old milestone names

| Old name | Now |
|---|---|
| M1 Catallax 0.1.0 alpha | Reached with Catallax 0.1.0. See the history for 25 to 29 Oct. |
| M2 to M5 The first part in SolidDraft, everyday features, distance features, general fillets and draft | 1.4 Everyday edge and wall features |
| M6 The roadmap and the decisions written down | Reached. This file replaces its plan. |
| M7 Lofts and spline sweeps finished | Reached. See the history for 10 and 11 Nov. |
| M8 Cuts after a fitted fillet | Part of 1.3 Cuts and joins on spline faces |
| M9 Spline intersections ported | Part of 1.3 Cuts and joins on spline faces |
| M10 Catallax cuts and joins spline faces | 1.3 Cuts and joins on spline faces |
| M11 Fillets between spline faces | 1.5 Fillets between spline faces |
| M12 The tolerance rules finished | Exit run recorded on 14 Nov. What is left belongs to 1.1. |
| M13 to M15 Field preview, fast evaluator, graphics-card viewer | 4.1 Field preview and first display |
| M16 FEA groundwork and M20 the first FEA result | 4.5 FEA, parked |
| M17 Clash, clearance and mass | 3.3 Clash, clearance and mass |
| M18 The field catches up | 4.2 The field catches up |
| M19 Lattices and printing | 4.3 Lattices and printing |
| M21 The field reads imported STEP and M23 Mesh and scan to a clean part | 4.4 Imports and scans in both kernels |
| M22 Closed mesh import in both kernels | Reached. See the history for 11 Nov. |
| M24 STEP and IGES import and editing finished | 1.2 Real-file import |
| M25 The spline-intersection port scoped | Recorded. See [SMLib intersection port scope](docs/research/smlib-intersection-port-scope.md). |
| M26 Speed optimisation | 5.1 Speed and simplification |
| M27 The spline-fillet proof scoped | Recorded. See [the offset search](docs/research/spline-fillet-offset-search.md). |
| M28 Numerical modelling under stated tolerances | 1.1 Numerical modelling under stated tolerances |
| M29 NURBS toolkit | 2.1 NURBS toolkit |
| M30 Surface bodies | 2.2 Surface bodies |
| M31 General sweeps | 2.3 General sweeps |
| M32 Direct edits | 2.4 Direct edits |
| M33 Drawing views from the kernel | 3.1 Drawing views from the kernel |
| M34 Assemblies and exchange | 3.2 Assemblies and exchange |
| M35 Tools only SMLib has | 3.4 Tools only SMLib has |
