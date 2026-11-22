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
- **One owner per lane.** One strong worker owns each lane until its current milestone meets the done-when line, then takes its next milestone. Each PR is one reviewable step with green checks. The [Checkpoint 2 lane assignment](#checkpoint-2-lane-assignment) replaces the phase 2 setup, retained in history. Phase 1 ran as four slices: 1.1 with 1.3, 1.2, 1.4, and the app.
- **Walk the done-when line first.** A milestone opens by walking its "done when" line end to end, through SolidDraft where the app has the feature, and listing every missing piece on the first day.
- **Delivery until phase 5.** Large PRs that cover a whole step. Tests are the fewest that prove the change works. No speed work: a speed idea becomes a one-line note in the PR's Later section and joins the phase 5 list. Work off the "done when" line becomes such a note too. The owner's rule of 16 Nov, kept until phase 5 at Checkpoint 1. 3.0A and 3.0B are design work and are not blocked by the no-speed-work rule, Q18.
- **Tolerance first.** Build fitted geometry in f64, measure the gap, store 1.1 x gap within the caps of T0 to T8, check once, and write STEP marked approximate. The predictor chooses the exact path where it is cheap. The proof button comes at 5.2. [Checkpoint 2 decisions](docs/checkpoints/checkpoint-2-decisions.md), Q1, Q13 and Q30, replace the numerical-first rule.
- Before writing a new algorithm, check how OpenCASCADE and SMLib solve it, and say in the PR which approach you follow.
- **The allowance until phase 5.** The application's default work allowance is 100M until phase 5, so slow parts still build in SolidDraft. 20M stays the target. Every edit over 20M records its measured work per stage and joins the [list of cases for phase 5](docs/phase5-work-cases.md). A request still refuses at the allowance in force, and a refused step is not a completed build.
- A SolidDraft feedback list goes to one worker as one PR, or two when the kernel changes too.
- When a lane's model quota runs short before a checkpoint, firstmate reports it and the owner decides what gives.

**Lessons, 15 Nov.** From the first ten hours of phase 1:
- Tiny cards cost more than they changed. Setup, build, a 50-minute CI queue and review rounds outweighed each change.
- Exact-first spent 13.6M work units rounding one fillet-end vertex. OpenCASCADE computes that vertex numerically.
- Nine app items as nine PRs queued on the same test runners.
- Work authored in the cloud needed a second worker to open each PR.

**Lessons, Checkpoint 1 (16 Nov).** From the whole of phase 1, 98 merged PRs in 48 hours:
- Waiting cost more than the work. 34 of the 98 PRs waited over 2 hours between opening and merge, through two or three review rounds of 20 to 60 minutes and 30 to 45 minutes of hosted CI. Several CI failures came from the runners, not the change. Hence the narrower review hold above; a runner failure gets a rerun, not an investigation, and builds and tests run on the work PC before a push.
- One worker per slice was the slowest path until the spline work split between two workers. Hence three kernel lanes in phase 2.
- Two pieces a "done when" line needed surfaced on the last day: drawing an oblique draft in the app and re-pointing references after a file replacement. Hence walking the line first.
- The measure moved during the phase. The census headline changed definition several times, the scoreboard became one command only at the checkpoint, three of its ten real parts test nothing and P6 to P10 have no definition. Hence the measures are fixed at the start of phase 2 and their definitions held to Checkpoint 2.
- Numerical first paid. Complete STEP solid imports rose from 186 to 418, 11 of the 14 spline Boolean requests reach STEP, and P1, P4, P5 and the crossing pipes build.
- Speed ideas parked as Later notes cost nothing and gave the checkpoint its phase 5 list. A side robustness test cost an hour.
- Model quota is a real limit: the Codex weekly quota ran close to empty and one Claude account hit its five-hour limit twice. Fewer, longer lanes.

The edit-chain scoreboard is the measure `VISION.md` asks for. The target parts P1 to P10 and ten real parts from the public STEP census are built and edited step by step in Catallax, OpenCASCADE and SMLib, with one definition of done for all three.

## Where we are

Phase 1 ended at Checkpoint 1 on 16 Nov 2026. [The checkpoint report](docs/checkpoints/checkpoint-1.md) has every number below. Catallax 0.3.0 and SolidDraft 0.6.0 are cut from the measured commit.

| Measure at Checkpoint 1 | Catallax | OpenCASCADE | SMLib |
|---|---:|---:|---:|
| Edit chains done, P0 to P5 (13) | 10 | 12 | 8 |
| Edit chains done, real census parts (10) | 6 | 7 | 6 |
| STEP census files reaching a valid solid (3,000) | 374 | 2,701 | |

Three of the ten census parts count for no kernel, because the test hole lands in an existing bore. The coverage count was not recounted at Checkpoint 1 and stands at 55% from 14 Nov.

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
| Curves and surfaces | 12 | 54% | Surface fitting. Gordon, Coons and N-sided surfaces. Bend and twist. | 2.2, 3.4 |
| Intersection | 6 | 58% | None. Tangent contact and some cone and torus pairs refuse. | 1.1, 1.3 |
| Solids and topology | 7 | 64% | Surface, wire and non-manifold bodies. | 2.3, 3.2 |
| Features that add or remove material | 11 | 55% | Extrude with taper. Sweeps with guides, twist or a helix. | 1.3, 2.4 |
| Edge and wall features | 9 | 39% | Curvature-continuous blends. Remove a face, fillet or hole. | 1.4, 2.1, 2.5, 3.4 |
| File exchange | 8 | 62% | Assemblies, colours, names and annotations in files. | 1.2, 3.2 |
| Meshes and queries | 8 | 56% | Mesh Booleans and decimation. Hidden-line views. | 3.1, 3.3, 3.4 |

### The three gaps that hurt most

These set the order of phase 1.

1. **Most real STEP files do not open as verified solids.** The same 3,000-file census now has 418 complete-source solid imports (13.93%), beside the original 186 (6.20%, about 6%). All expected solids export and pass OCCT builder validity in 384 files (12.80%); builder processing may heal internally. The new consumer also retains 119 mixed/partial files and 2,202 unverified source-surface sets, for 2,739 supported source openings (91.30%). These are separate from solid admission and checked exports. OpenCASCADE's reference remains 2,701 healed-valid solid files (90.03%) under its different protocol. [The census](docs/research/step-census-opened/README.md) records every input and remaining cause.
2. **The loft top cuts build above the 20M target.** All 14 recorded spline Boolean requests now write checked approximate STEP. The three original top requests build at the 100M application allowance, using 26.78M, 27.29M and 25.95M work. Their files pass raw OpenCASCADE validity with healing off and a later cut. [Stage work and checks](docs/phase1/spline-boolean-top-cuts/REPORT.md) join the phase-5 list.
3. **Fillets, shells and draft work on certified families only.** Edge and wall features are the weakest area at 39%, and no line in it is complete.

A work unit counts exact arithmetic and is the same on every machine. The target per edit is 20 million units, about 2 to 4 seconds of proof work on one core. Until phase 5 the application allows 100 million.

Twelve things neither reference kernel does are not in the count. They include exact numbers for exact pieces, the same answer on every machine, the field kernel, lattices and scan to part.

## Plan

### Phase 0: close the release

**Doing:** making the 0.2.0 and 0.5.0 cut clean. **Why first:** everything after starts from a working release. **Unlocks:** a fixed point to measure from.

- Both tags published. SolidDraft pinned to the Catallax tag.
- Every download installs, builds a part and exports it. The browser app loads, opens, edits and exports.
- Release notes list the slow steps and the two export timeouts, with working links.
- Faults found at the cut: the macOS test panic is fixed on main. The lattice and gyroid export that timed out is rerun on a quiet machine.
- Docs match the cut. This file replaces the old plan and roadmap files.
- Three quick re-checks: whether the later cuts on the filleted P5 housing already build (deferred to 2.1), the last 72,166 units on P4's scoreboard run, and P1 with its fast route off.

**Checkpoint 0.** This roadmap is adopted. Coverage is 55%. Review the 13 pending PRD proposals with the owner (section 14 of [the archived plan](docs/archive/PLAN-2026-11.md#14-proposed-prd-changes)).

### Phase 1: close the three biggest gaps against OpenCASCADE

**1.1 Numerical modelling under stated tolerances**
- Doing: where an exact proof is unavailable, because it is missing or too costly, the kernel builds the geometry numerically within a stated tolerance. Each result records its tolerance, the checks it passed and any repair. Its STEP file is written and marked approximate.
- Why here: it is the base under all three gaps. Real files miss exact rules by micrometres. Spline Booleans stop on a gap or a coincidence they cannot prove. OpenCASCADE and SMLib work to a tolerance throughout. We will work to one only where the exact proof is missing or too costly.
- Unlocks: milestones 1.2 and 1.3 below, the spline fillets of 2.1, and the first bounded-answer rule, for [tangent contact](docs/research/tangent-touch.md).
- Done when: numerical results carry their tolerance, checks and repair history through import (1.2) and spline cuts (1.3). A STEP file marked approximate opens in OpenCASCADE with healing off and takes an edit there. The curved-housing chain (loft, spline cut, fillet across patches, one more edit, STEP opened by OpenCASCADE) closes with the spline fillets in phase 2 (2.1).

**1.2 Real-file import**
- Doing: opening real STEP files under those tolerances, analytic files first, then spline files. What still refuses is counted by cause.
- Measured: the [full frozen 3,000-file rerun](docs/research/step-census-opened/README.md) gives 418 complete-source solid imports (13.93%), beside the original 186 (6.20%, about 6%). Every expected solid exports and is valid after OCCT builder processing in 384 files (12.80%), with possible internal healing. There are also 119 mixed/partial files and 2,202 retained unverified source-surface sets; 259 inputs refuse recognition/geometry reading and two opening queries time out. No solid tier or export-valid file was lost against the landed 405/374 baseline; all ten real P06 controls are retained. Kernel verified correspondence and atomic source replacement are landed; SolidDraft replacement/reference UI remains in the app slice. The distinct OCCT reference has 2,701 healed-valid solid files (90.03%).
- Why here: it is the largest single gap, and it is the first thing a user meets: a bought-in part that will not open.
- Unlocks: bought-in parts in every later feature, the ten real parts on the scoreboard, direct edits (2.5) and imports in the field (4.4).
- Done when: the 3,000-file census is rerun and the share that imports is published beside the 6%. Replacing an imported file re-points every reference.

**1.3 Cuts and joins on spline faces, with what is left of the SMLib intersection port**
- Doing: Booleans on lofts and imported freeform parts that finish and export.
- Why here: cuts and joins shape every part after its first feature, and today they stop at the first freeform face. It comes after 1.1 because 4 of the 14 refusals are unproved gaps and coincidences, which the rules of 1.1 are meant to settle.
- Unlocks: every freeform part past its first cut, and the spline fillets of milestone 2.1.
- Done when: the 14 recorded spline Boolean requests and the loft top cut build and reach STEP. The later cuts on the filleted housing are not a criterion here. They are deferred to 2.1.

<a id="everyday-edge-and-wall"></a>**1.4 Everyday edge and wall features**
- Doing: fillets over mixed corners, shells with side openings, hollow oblique draft, chamfer corners and general holes on ordinary prismatic and turned parts.
- Why here: two thirds of the census files have only flat and round faces, so this is where most real parts are. It touches different code from 1.1 to 1.3 and runs beside them. The owner's bar from 7 Nov stands: at least match OpenCASCADE on fillets.
- Unlocks: the five everyday target parts end to end, and the first complete line in the weakest area.
- Done when: P1 to P5 each build, edit, save, reopen and export through SolidDraft.

**App, alongside:** browser first. Kernel work runs off the interface thread, commands disabled in the browser work there, and a part built within tolerance shows as finished. Why here: the browser is where a new user tries a real file.

**Checkpoint 1, 16 Nov.** Reached before the spline fillets. State of each milestone against its "done when" line:

- 1.1: done. Numerical results carry tolerance, checks and repairs through import and spline cuts, and approximate STEP opens raw in OpenCASCADE and takes a later cut. The curved-housing chain moved to 2.1.
- 1.2: done by its line. The census is published at 374 valid files (12.47%) beside the old 6.20%, and replacing an imported file re-points every reference it can match and names the rest. Import is still the largest gap against OpenCASCADE.
- 1.3: open. 11 of the 14 recorded requests reach checked STEP inside 20M. The three loft top cuts carry into phase 2, lane 1.
- 1.4: done. P1 to P5 build, edit, save, reopen and export through SolidDraft. The scoreboard still shows P3's membership point and P5's 16 unnamed faces or edges.
- App: kernel work runs off the page. Cancel and the edit after it stay over 250 ms in the browser; the fix carries into the phase 2 app lane.

### Phase 2: the modelling lines we do not have yet

Why this phase is second: these are whole lines of the count with nothing built. They matter less than phase 1 because a user can work around a missing surface tool, but not around a file that will not open or an edge that will not round.

The phase 2 lane setup is retained in [history](#phase-2-lane-setup-superseded-at-checkpoint-2). The [Checkpoint 2 assignment](#checkpoint-2-lane-assignment) governs new work. Existing phase 2 open items remain open; P8 remains parked. Coverage extensions are listed after the existing milestones, without renumbering or dropping their current done-when lines.

| # | Milestone | What we are doing | Why here | Unlocks | Done when |
|---|---|---|---|---|---|
| 2.1 | <a id="spline-fillets"></a>Fillets between spline faces | Rounded edges on lofted and imported freeform parts, across more than one surface patch. | Needs 1.3. The centre curve of a spline fillet is found by the same engine that intersects spline faces. [One loft edge already rounds within 59 nm](docs/research/spline-fillet-offset-search.md). Moved from phase 1 on 16 Nov. | Finished freeform parts. P9, the lofted duct. The curved-housing chain of 1.1. | A fillet across several patches builds on a loft and on an imported part, and a radius edit keeps the names. The curved housing builds as loft, spline cut, fillet, one more edit, its STEP opens in OpenCASCADE with healing off and takes an edit there, and the later cuts on the filleted housing (old M8) build. |
| 2.2 | <a id="nurbs-toolkit"></a>NURBS toolkit | Freeform surface construction: Coons patch, surface through a grid of points, Gordon surface, guided skin, N-sided fill. SMLib's surface-construction code may be adapted, with attribution, under the 16 Nov widening of the port rule. | Lane 3, from the first day. Three "not started" lines sit here. Surface bodies and direct edits need its tools. | Freeform surfacing. Lofts with guide curves. | Each surface builds and exports. Knot insertion, degree raise, split and join are public. |
| 2.3 | <a id="surface-bodies"></a>Surface bodies | Open surfaces and wires as bodies. Sew, trim, extend, thicken into a solid. | Lane 2. Sewing comes first, as stitching imported surface bodies into checked solids; trim, extend and thicken follow 2.2. Direct edits (2.5) need its trim, extend and sew. | Surface modelling. Imported files that are not closed solids. | A surface trims, extends and thickens into a solid. Faces sew into a solid outside import. |
| 2.4 | <a id="general-sweeps"></a>General sweeps | Extrude with taper. Sweep with a guide, with twist and along a helix. Extrude up to a curved face. | Lane 3, after 2.2. Two "not started" lines. It depends on neither 2.2 nor 2.3. | Springs, threads, tapered and twisted shapes. P8, the spring and threaded boss. | P8 builds. |
| 2.5 | Direct edits | Delete a face, a fillet or a hole and heal the solid. Move, offset or replace a face. | Lane 3, after 2.3, while lane 1 keeps 2.1. Needs 2.2 and 2.3. It pays only once imports work (1.2). | Editing bought-in parts that have no feature tree. | Each edit works on an imported solid. |

**2.4 partial delivery:** taper, guide/twist sweeps, a quarter helix and curved-face termination have build/edit/save/reopen and unhealed STEP evidence. General-motion fits are explicitly numerical. P8 clears the earlier 100M carrier-comparison budget barrier but still refuses `intersect.parametric-surface` during threaded-union clipping at 55,626,285 work; the milestone remains open. Its [stage record](docs/RECIPE.md#work-recorded-for-phase-5) is retained for Phase 5.

<a id="phase-2-app"></a>**App, alongside:** import health report, reverse-engineering workspace, direct-edit commands, reference axis and point, slot, draft check, exact sketch solver.

**Checkpoint 2.** Scoreboard and census rerun, and the coverage recount. Order of phase 3.

### Phase 3: drawings, assemblies, exchange and the last lines

Why this phase is third: it completes the count. Each item has a working stand-in today (mesh outlines in drawings, assemblies in the app), so parts get made without it.

**3.0A Tolerance-first core**

[Checkpoint 2 specification, section 4.1](docs/checkpoints/checkpoint-2-decisions.md#41-item-30a-tolerance-first-core), Q1, Q5 to Q11, Q18 and Q25. Starts before any other implementation, after owner approval of this roadmap PR. This is design work, not phase-5 speed work.

Owner Fable; every PR reviewed by the dot. Starts now.

Contents, in order:

1. The 14-edge known-answer test set and its three guards (below). Built first, so the new checker is measured before it replaces the old one.
2. The A+ checker ([decisions, section 4.2](docs/checkpoints/checkpoint-2-decisions.md#42-the-a-checker)).
3. Copy tagging ([decisions, section 4.4](docs/checkpoints/checkpoint-2-decisions.md#44-copy-tagging)).
4. Tolerance rules T0 to T8 ([decisions, section 4.3](docs/checkpoints/checkpoint-2-decisions.md#43-tolerance-rules-t0-to-t8)), including the export receiver profile.
5. Check once. Export writes the numbers and adds the rounding error to the tolerance; it does not re-admit the body. Reopening trusts the tolerance recorded at export. Today export re-admits every face (21.9M to 42M work on measured cases) and reopening admits again (47 s to 7 min on the R2 rim). The kept test file for opening speed is `owner-pocket-rim-r2.step` (private acceptance fixture, not bundled here).
6. The gap rule for imports (Q9).

Done when:

- The checker passes all 14 known-answer edges: on each, the reported gap is at least the true gap and at most 1.1 times it.
- The three guards hold (below).
- Export no longer re-admits; reopening the R2 rim file no longer re-checks; the work and time of both are recorded before and after.
- Every existing exact answer is unchanged.

The 14 edges, each with an exact known maximum gap:

1. A narrow peak between the fitter's eighths, on an edge with 100 to 300 knots.
2. A peak set only by the surface's knot lines.
3. A narrow peak on a degree-25 surface.
4. A rational end weight, and a rational rail against a cubic face-curve.
5. Origin at 1e6 to 1e12.
6. A peak exactly midway between two equal samples.
7. Seam crossing, a one-period shift, and a full closed circle.
8. A 1.0000057 parameter stretch.
9. A face-curve outside the face domain.
10. A traced section between faces meeting at under 5 degrees.
11. A fit that passes at eighths but fails at sixteenths.
12. A 256-span edge, for work cost.
13. Zero-gap controls.
14. A circle against its rational-spline copy with known sag.

Each one copies a recorded failure in OpenCASCADE or SMLib (section 9.6).

The three guards against tuning to the test:

1. Fresh cases each run. Besides the 14 fixed edges, every run generates new random edges, each with a hidden bump of known size and position.
2. The checker's point count depends only on the shape (knots, seams, poles, piece complexity), never on which edge it is. A reviewer can see this in the code.
3. Every method in the checker names its source (paper or kernel). No unexplained constants.

Coverage: infrastructure for existing lines; no line becomes full from this milestone alone.

<a id="one-model-and-predictor"></a>**3.0B One model and the predictor**

[Checkpoint 2 specification, section 4.5](docs/checkpoints/checkpoint-2-decisions.md#45-item-30b-one-model-and-the-predictor), Q12 to Q14, Q18, Q25 and Q28. Follows 3.0A, before 3.1.

Owner Fable, after 3.0A; the dot reviews every PR.

One model (Q12):

- One topology type. Vertex, edge and face gain an optional tolerance; none means exact. Geometry stays stored as it is (even fitted geometry is stored as exact rationals today), so this changes topology, not geometry.
- Exact-exact decisions keep today's exact predicates. Any decision with a tolerant participant is "distance against the sum of tolerances" through a filtered predicate: f64, then interval, then exact (the CGAL and Shewchuk pattern; Catallax already uses it in mesh and field).
- Tolerant entities never reach the exact intersector. Every slow case on record is exact arithmetic on fitted or near-tangent data.
- An exact face never takes a tolerance. Only vertices and edges forced by a tolerant neighbour do (vertex >= edge >= face).
- Re-tighten tolerant edges between exact surfaces after each operation, as Parasolid's `PK_EDGE_optimise` and ACIS tweak do.
- The `features.numerical-continuation` refusal goes away: a hole after a numerical fillet works.
- The numerical checker merges into `Model::check`. `NumericalBody` becomes a view or goes.

Why not the pure CGAL route: on curved faces, exact decisions need exact roots (degree 4 for two cylinders, 324 for two bicubic splines), so the exact fallback is impractical on splines; exact decisions on rounded floats lose design intent (two coaxial cylinders are no longer coaxial); and planes and cylinders lose their exact formulas. The useful part of CGAL, filtered decisions, is inside the one-model design above.

The predictor (Q13):

1. First fix the work meter. Report refunded speculative work in its own counter instead of erasing it (`catallax-numbers/src/operation.rs:192-233`, called at `evaluate.rs:1499`). Charge work in `catallax-step/src/fit.rs` by limb product, as the decider does (`decider.rs:1659`). Today exact STEP export took 8.5 s on 807k counted work.
2. A quick look reads degree product D, limb count L of the rationals, and proven bounding-box separation. It costs at most 1% of the operation.
3. Two cap tiers by D x L^2: low-complexity inputs get the 20M target for the exact attempt; fitted splines stored as 53-bit fractions get about 2M.
4. If exact finishes under the cap, the result is exact. If it hits the cap, the fitted path takes over.
5. No cap where no fitted path exists yet; exact is then the only way to a result.
6. The normal angle is used only to tune the fitted path and to demote spline pairs, never to demote analytic pairs (tangency is cheap in exact arithmetic).
7. No memo per recipe step and no learned model.

On the measured cases, a fixed rule by shape kind gets two results wrong that the predictor gets right: the far-disjoint loft cut (k12, builds exact at 9M) and the spline top cuts (26M to 27M, built exact). The predictor also saves 62M to 79M per doomed loft cut.

Predictor done when, on a frozen corpus (all phase-5 rows, checkpoint-1 parts, the box fillet, and a held-out ABC subset chosen before tuning):

- path choice and work are identical at 1, 2, 4 and 24 threads and on two machines;
- the quick look costs at most 1% of each operation's work;
- everything that builds exact under 20M today still does;
- refunded work stays within its tier cap;
- wall time per million work units stays within a 10x band across all rows;
- thresholds are unchanged on the held-out set.

The 2M and 20M tiers are judgement and need measuring.

Open sub-question for 3.0B: whether exact entities use no tolerance or a floor like Parasolid's 5e-6 mm. T0 says none; the research suggested measuring the M8 loft cuts under both before the code fixes it.

Done when, for one model: exact and tolerant entities coexist under the rules above; a hole after a numerical fillet works without `features.numerical-continuation`; the numerical checker is part of `Model::check`. The predictor's frozen-corpus done-when lines above also hold. Coverage: infrastructure and continuation within existing lines, not a full-line promotion.

<a id="drawings-and-pmi"></a>**Drawings and PMI, alongside and just behind 3.0B**

- Doing: AP242 semantic PMI stored on the part, automatic drawings from it, feature recognition on bare imports and a rules-based interview. [Checkpoint 2 decisions, section 4.8](docs/checkpoints/checkpoint-2-decisions.md#48-drawings-and-pmi-the-p1-bracket-test) and the [P1 bracket specification with D1 to D6](docs/checkpoints/checkpoint-2-drawings-pmi.md) define the scope. One owner across kernel and app. Fable, moving to Claude Opus if Fable quota runs short; the dot reviews.
- Why here: the owner adopted the P1 bracket test, Q3 and Q20. Runs mostly in parallel wherever it does not touch 3.0A or 3.0B code. Its roadmap place is just behind 3.0B, not a renumbering of 3.1.
- Unlocks: manufacturing information travels in STEP and produces the drawing without hand layout. The PMI owner owns annotations within 3.2; the exchange lane owns structure, names and colours.
- <a id="drawings-pmi-acceptance"></a>Done when: C1 round trip of all 14 items passes the NIST analyzer, OCCT reads all 14, and re-import, save and reopen keep their attachments; C2 bore edit 8 to 10 mm keeps its tolerance; C3 one button gives STEP plus PDF and DXF, every item once, no overlaps or hand moves, with title block; C4 no drawing information is absent from the file; C5 bare-import recognition finds bore, boss, four holes and slot, and the interview asks only for missing information. These are the specification's C1 to C5, not coverage IDs. Coverage: X8 annotations, partial until 3.2's structure, colours and names also round-trip.


| # | Milestone | What we are doing | Why here | Unlocks | Done when |
|---|---|---|---|---|---|
| 3.1 | Drawing views from the kernel | Exact section curves, silhouettes on curved solids, hidden-line views. | Drawings work today from meshes. This makes them exact. | Drawings that match the model to the last digit. | SolidDraft drawings use them in place of mesh outlines. |
| 3.2 | Assemblies and exchange | The kernel holds several solids as an assembly. STEP carries structure, colours, names and annotations. | Needs import (1.2) to be solid first. | Whole products in one file, both ways. | A STEP assembly round-trips. IGES exports spheres and tori. |
| 3.3 | Clash, clearance and mass | Clash, clear by a stated gap, or unresolved, for every pair. | Exact overlap already works. Distant and nested pairs are left. | A proven gap between parts. Refusals that state the largest radius or thickness that fits. | Every pair in a test assembly gets one of the three answers. |
| 3.4 | <a id="smlib-modelling-tools"></a>Tools only SMLib has | Bend and twist. Curvature-continuous blends. Mesh Booleans and decimation. | Last because OpenCASCADE does not have them either. | The final three lines. | Each has a first build. |
| 3.5 | <a id="one-display-mesher"></a>One display mesher | Mesh every face by one path. Sample shared edges first, triangulate each face in its own parameters and set density from second-derivative bounds. Preview, bounded and proven use the same path. | SolidDraft cannot draw the filleted M11 loft. Each new spline feature would otherwise need its own grid scheme. | Display of every spline body, whole-body two-kernel checks and per-face reuse after edits. | See below. |

**3.5 steps**

1. Freeze the display corpus and record current work, point counts, bounds and refusals. Include P1 to P5, the [phase 5 work cases](docs/phase5-work-cases.md), the three M11 loft fillets, the small-warp loft, saddle and ellipse fillets, phase 2 spline bodies, and the SolidDraft three-turn coil.
2. Put spline charts on the shared path. Bound derivatives per span. Project in floating point, starting from the neighbouring point, and include the projection residual in the bound. Split charts where first derivatives vanish.
3. Plan each edge once using the stricter of the curve chord rule and the chord and normal limits of both incident charts. Cache those checks and sample one polyline on the edge curve for both faces, without retries or exchanges between faces. Record the charged work on the frozen corpus. Reuse the analytic path's trim triangulation and keep its chord checks and the spline angle and no-fold checks. Remove the junction, spline prism and oblique family grids.
4. Make preview, bounded and proven levels of one option. Fold the preview routes into it.
5. Cache by face and tolerance and make the levels nest. After an edit, mesh only the changed faces.
6. Share the bound formula and triangle check between B-rep and field display.
7. Draw the 4.2a rolling fillets and patch shells with the shared charts. Add rolling and shell trim handling, a swept-arc chart for rolling pieces and the spline chart for shell faces.

**3.5 done when**

- Every body in the frozen corpus displays at 0.01 mm and 0.005 mm chord with no refusal, including the three M11 loft fillets.
- At 0.005 mm each M11 loft body takes under 1M display work. Display work for each analytic body is at most its baseline.
- Both faces at every shared edge use the same polyline. Every B-rep mesh passes the existing closed-mesh check and retains chord, angle and no-fold checks.
- Independent dense sampling on every face finds no point farther from the mesh than the reported bound.
- The point count predicted before meshing equals the point count produced.
- Mesh bytes and work agree at 1, 2, 4 and 24 workers, and on two machines.
- `catallax-mesh` has at least 2,000 fewer lines than at the baseline.
- SolidDraft draws the M11 loft fillet with one later feature on top (R2).
- Compare's whole-body forward pass runs on the M11 rows and the [P5 flanged housing](parts/P5-flanged-housing.json) inside 100M work.
- SolidDraft's three-turn helical sweep (2 mm wire diameter, 34 mm coil diameter, 3 mm pitch) displays inside 100M work with its periodic seam represented in the display chart. Target under 20M display work for a few turns.
- Field chart display draws each 4.2a fillet row the field builds, including the filleted M11 loft, and the patch shells without `field.display-unsupported`. The same dense-sampling test checks their proven two-sided bounds.

This sets the full bar for [coverage line Q1](docs/COVERAGE.md#meshes-and-queries-8). Keep the existing display entry points, the 20M work target and the 100M allowance.

**App, alongside:** assemblies in the browser, assembly drawings, BOM and balloons, Pack and Go, the remaining mates, drawing annotations, URDF and USD export with joints, the refusal assistant.

**Checkpoint 3.** Coverage is judged against the new done-when lines in the [coverage follow-on](#coverage-follow-on-complete-the-61-lines) and [docs/COVERAGE.md](docs/COVERAGE.md), not merely "every line has something built". We publish what still refuses and why. A line with its bar still to be set cannot be claimed full. The owner decides what to resolve before moving on.

### Phase 4: field, lattices, printing and scans

Why this phase is fourth: this is what neither reference kernel has. It is our edge, not a gap, so it waits until the gaps are closed. It never blocks a B-rep feature in the meantime.

| # | Milestone | What we are doing | Why here | Unlocks | Done when |
|---|---|---|---|---|---|
| 4.1 | Field preview and first display * | The app draws field previews only for field-only features, lattices and gyroids. The B-rep previews everything else. Remove proved and bounded lights; keep a refused message and a working indicator. | First in the phase: field-only features need a display. [Checkpoint 2 decisions](docs/checkpoints/checkpoint-2-decisions.md), Q31 and Q32. | A feature shows at once. | A preview never stands in for an export. |
| 4.2 | <a id="field-modelling"></a>The field catches up | Field versions of fillets, sweeps and lofts. The field may take over a B-rep step's result only when that is the only way to build the part, and only with no downstream problem, Q15. | The field must hold every part before it can print or analyse it. | Every part has both descriptions. | P5 builds in the field with no fillet refusal. 4.2a below defines the independent fillet route. |
| 4.3 | Lattices and printing | Printable lattice parts. The missing piece is [a proof that a mesh keeps its cavities](docs/research/adaptive-printing.md). | Needs 4.2. | P6 and P7. Sheet gyroids. | P6 exports a printable mesh at its recorded accuracy. |
| 4.4 | Imports and scans in both kernels | Imported spline parts in the field. [Scan fitting](docs/research/scan-fitting.md) grows past planes and cylinders. | Needs import (1.2) and the spline work (1.3, 2.1). | Bought-in parts and real scans in the field. | An imported spline part builds in the field. A real scan becomes a clean part. |
| 4.5 | FEA groundwork, first result and P10 | [Parked by the owner on 13 Nov](docs/research/fea-first-beam.md). | Starts on the owner's word at a checkpoint. | Stress on the same part the B-rep exports. | |

**4.2a The field builds every fillet on its own**

[Checkpoint 2 specification, section 4.6](docs/checkpoints/checkpoint-2-decisions.md#46-item-42a-the-field-builds-every-fillet-on-its-own), Q15 to Q17, Q21 and Q23. A sub-milestone inside 4.2, not a change to any old number. Any extra building block field parity needs is added in the order it is needed.

A numbered sub-milestone of 4.2, "the field catches up". Constant-radius rolling-ball fillets only; variable radius, setback vertex blends and conic or G2 sections keep their typed refusals.

Rule: the field derives its own ball paths. It does not reuse the B-rep's `rolling_pipe`, `rolling_saddle`, `rolling_sphere` or `rolling_mitre`, because a shared bug would fool both kernels. The milestone document names any shared dependency (the M10 pair search is one).

Field queries return inside or outside plus a guaranteed distance range near the surface (Q23).

Steps, easiest first:

1. Turn `rolling_region` into one general routine for distance to a path piece. Today's plane-plane and plane-cylinder answers stay bit-identical.
2. Plane with cylinder in any pose, cone or sphere: the oblique, oblique-cone and generator test files.
3. Quadric pairs: saddle, coaxial, sphere-ring and equal-crossing files. Paths have closed forms (Dupont, Lazard, Lazard, Petitjean 2008).
4. Fillet ends that touch an edge or vertex: P5 at radius 2, rim-contact, hole-ring and slot files.
5. Fillet over fillet, using the nested path.
6. Torus pairs.
7. Spline faces (below).

Done when, for steps 1 to 6:

- P5 builds in the field at radius 2 with no fillet refusal (the existing 4.2 line);
- for every constant-radius rolling-ball case in the fillet test files (list frozen at the start), catallax-compare finds every B-rep display vertex within the field's bound at the part's tolerance;
- each remaining refusal is named with its reason.

Step 7, spline faces. Viability study verdict: go, with conditions.

Measured on this machine (release build, i7-9750H; logs in [docs/research/cp2-grill/viability-logs/](docs/research/cp2-grill/viability-logs/)):

- M10 spine search 2,332,939 / 2,666,401 / 6,793,048 work on the three M11 edges, identical at 1, 2 and 4 workers, 0.14 s to 0.40 s per edge;
- full B-rep fillet 5.25M / 6.97M / 16.43M work, within 78 / 5 / 80 nm of the exact fillet.

4.2a remains open. See [the step 7 evidence and remaining refusal](docs/M42A-FIELD-FILLETS.md#steps-75-to-77-near-tangency-and-a-plane-chart).

Step 7 sub-steps:

- 7.0 The query contract is already set by Q23; no new decision is needed.
- 7.1 Spike: certified point-to-patch distance on the M27 patch pair, in catallax-field, timed over 10^4 points. This retires the biggest risk.
- 7.2 The spline-face building block for whole, untrimmed patches, with ray-crossing sign. Trimmed faces go to 4.4 with a typed refusal until then.
- 7.3 The field's own spine proof: prove the two end points on the face borders, then a parametric Krawczyk proof span by span. Tube width eps <= min(100 nm, part tolerance / 2). A build that cannot reach it refuses and reports the eps reached.
- 7.4 The swept-ball distance query on that tube.
- 7.5 Near-tangent sweep on a synthetic edge: record the smallest face angle proved inside 20M, then decide the tail policy (deflation per Leykin et al. 2006, or a labelled bound).
- 7.6 Spline-with-analytic pairs: write the analytic offset (same family) as a rational patch with offset 0.
- 7.7 Typed refusals, including a singular-contact control.

Step 7 done when, on the three M11 edges plus one spline-plane edge:

- the field proves its own eps <= min(100 nm, tolerance / 2) inside 20M work, identical at 1, 2 and 4 workers;
- every B-rep display vertex of each fillet face lies inside the field range widened by the B-rep's reported bound;
- a dense-sampling third check in tests agrees (catches a bug shared through M10);
- median query time is recorded;
- a near-tangent control refuses or returns a labelled bound.

Why a proven tube and not an estimate: the field is the yardstick for B-rep bugs (Q16). An estimated path cannot tell whether the B-rep path is wrong. It runs in tests only, and it proves only down to the tolerance, never past it.

Caveat: the three M11 test edges are single-span bilinear or biquadratic patches, so they say little about real imported parts.

Coverage: independent field parity is outside the 61 reference-kernel lines. No extra coverage points.


<a id="phase-4-app"></a>**App, alongside:** field-driven thickness, more grades, fields on imported bodies, infill that is denser near the skin or the load.

**Checkpoint 4.** Order of the speed work, from measurements.

### Phase 5: speed and simplification *

<a id="speed-and-simplification"></a>**5.1 Speed and simplification.**

Why last: the owner's ruling of 11 Nov. Clock targets were stalling capability work, and speeding up a method that is about to be replaced is wasted.

- P1 to P5 and their edits meet the click targets in SolidDraft: 100 ms for stock, 500 ms for a Boolean, 1,000 ms for a feature.
- The GPU and every core carry what they can. [Booleans gain nothing from extra cores today](docs/research/run-ahead.md).
- Simplify. Routes that the measurements show are not needed are removed.
- Slow steps recorded at the 0.5.0 cut: loft 35 s, fillet on a loft edge 26 s, spline sweep 13 s, P5 fillet edit 4 to 6 s, field mass 3.5 to 7.7 s.
- Later: the 2.4 P8 thread sweep exceeds 20M; its [per-stage record](docs/RECIPE.md#work-recorded-for-phase-5) distinguishes the completed fit from the refused whole-part Boolean. The 100M application allowance is unchanged.
- Windows packages and signed downloads.


**5.2 The proof button**

[Checkpoint 2 specification, section 4.7](docs/checkpoints/checkpoint-2-decisions.md#47-item-52-the-proof-button), Q10.

A real whole-edge proof on demand for one part (Q10 B). Late and low priority. Built from parts Catallax already has: exact range arithmetic (`evaluate_box`, `jet_box`, Ball and Interval) and a one-dimensional branch and bound that splits wherever an upper bound beats the best sample. Nothing in Catallax computes a certified maximum of |C(t) - S(P(t))| today. Expect minutes per part and "could not prove" on blends and sweeps (solidcraft's experience: one edge took 153 s, some bounds came out 60x loose, blends were not provable). VISION.md's "no proof jobs or proof services" still holds: the button is one computation on request, with nothing stored and no service.

- Done when: on request for one part, the button runs a real whole-edge proof, not sampled admission, and reports proof or "could not prove". It stores no proof job and runs no proof service. Coverage: this is outside the 61 reference-kernel lines.

### Coverage follow-on: complete the 61 lines

[Checkpoint 2 decisions, section 4.10](docs/checkpoints/checkpoint-2-decisions.md#410-coverage-the-12-changes-for-100), Q2. These additions follow the current milestones and do not reorder them. Existing demonstrations remain required, including the open 2.1 and 2.5 lines and 2.4 with P8 parked. The extensions below widen their final bars. Every coverage ID named below targets **full**, 1 point, not a first example; unresolved full bars remain "bar to be set" in [docs/COVERAGE.md](docs/COVERAGE.md), not invented acceptance criteria.

Twelve changes: seven widen an existing done-when line and five add a new milestone line. Each widened done-when line moves from a named demonstration to a generated-recipe or census pass rate, with the line's Doing scope as the test. Underneath all of them: write `docs/COVERAGE.md` with the full bar for each of the 61 lines, and make every done-when line name the coverage lines it closes and at which level.

#### Curves and surfaces

1. <a id="coverage-guided-loft"></a>**Widen 2.2 (C7, C8, C9, C10, F8).** Add least-squares fitting of curves and surfaces to points within a stated tolerance, G1 N-sided and Gordon fills, curve networks whose corners miss by up to the stated tolerance, and lofts with guide curves. Done when each builds on a generated set without refusing possible geometry. *A user can fit a smooth surface to a point cloud, fill an n-sided hole with a tangent patch, and loft through sections along guide rails.*
2. **New line: offsets and projection (C6, C11, E6).** Offset any curve, surface or whole solid, spline faces included, within a stated tolerance. Project curves onto surfaces. *A user can offset a freeform part's skin by 2 mm or wrap a sketch onto a curved face.*
3. <a id="coverage-bend-twist-and-blends"></a>**Widen 3.4 from "a first build" to general input (C12, E8, Q8).** *A user can bend or twist any body, blend faces with curvature continuity, and run Booleans and decimation on meshes.*

#### Intersection

4. <a id="coverage-intersections"></a>**New line: intersections complete (I1, I3, I4, I5).** Every analytic pair, skew cone/torus and cylinder/torus included. Spline against spline, and curve against surface, with no refusal on the generated set. Tangent and near-tangent contact answered under the bounded-answer rule. *A user's cut or join never stops because two faces touch or meet at an awkward angle.*

#### Solids and topology

5. <a id="coverage-spline-booleans-and-fillets"></a>**Widen 1.3 and 2.1 (S3, E3).** Done when Booleans and fillets on spline-faced solids pass a generated set and the census parts, not 15 named requests. *A user can cut, join and round any freeform part, imported or modelled.*
6. <a id="coverage-surface-bodies"></a>**Widen 2.3 (S6, S7, X7) and add split and imprint (S5, F10).** Sew any face set that closes. Non-manifold bodies. Sheets and wires export as sewn shells. Split a solid by a plane, sheet or body, keeping every piece. Imprint a curve on a face. *A user can sew a patchwork of imported faces into a solid, slice a part in two, and imprint a logo outline on a face.*

#### Features

7. <a id="coverage-general-sweeps"></a>**Widen 2.4 from "P8 builds" (F4, F5, F6, F7, F9, F11).** Taper, guide, twist, helix and up-to-face each build on line, arc, ellipse and spline profiles along any path. *A user can sweep any profile along any path, with a guide or twist, and extrude up to any curved face.*

#### Edge and wall features

8. <a id="coverage-edge-and-wall"></a>**New line: edge and wall features in general (E1, E2, E4, E5, E7).** Constant and variable-radius fillets with any law on any edge. Chamfer roll-over and general corners. Shells and draft on any face, spline and oblique cylinder included. Done when the census parts take each feature without refusal of possible geometry. *A user can round, chamfer, shell or draft any edge or face of a bought-in part, as in OpenCASCADE.*
9. **Widen 2.5 (E9).** Each direct edit works on any face type, through SolidDraft. *A user can delete a fillet or move a face on any imported part.*

#### File exchange

10. **New line: STEP import parity (X2, X4).** A census target set against OpenCASCADE's 90.03% healed-valid figure. *Most real STEP files a user is sent open as checked solids.*
11. **Widen 3.2 (X8).** Done when colours, names and annotations round-trip, beside the assembly. *A user can send and receive a full product file with its structure, colours and annotations.*

#### Meshes and queries

12. **New line: geometry queries (Q4, Q5), and widen 3.3 (Q6).** Closest point on any curve or surface. Curvature and local surface analysis. 3.3 done when no pair answers "unresolved" outside the stated tolerance. *A user can check curvature on a freeform face and get a proved gap for every pair in an assembly.*

With these twelve changes every line reaches 1 at its done-when line: 61 points, 100% of OpenCASCADE's 58 and SMLib's 52. The "or more" can only come from the twelve extras outside the count ("Twelve things neither reference kernel does"). They never raise the count itself above 100%.

**Final checkpoint.**

### Lanes that never stop

- The resolve list. Every refusal of possible geometry is sorted by cause and fixed in batches. Cases a target part hits go first.
- Each new kernel operation reaches SolidDraft as it lands, in the browser first.
- Generated recipes with built-in checks run nightly and publish counts by outcome.
- The SMLib port supplies the units this roadmap needs next.


### Checkpoint 2 lane assignment

[Checkpoint 2 decisions, section 8](docs/checkpoints/checkpoint-2-decisions.md#8-who-builds-what). The roadmap update PR needs owner approval before any implementation worker starts, Q29. Then 3.0A, 3.0B, drawings alongside where independent, 3.1 to 3.4, phase 4 with 4.2a inside 4.2, phase 5 and 5.2. The app lane runs alongside: R2 scenario selection, remove proved and bounded lights while keeping the refused message and a working indicator, narrow 4.1 preview, then mixed exact/tolerant bodies and continuation after fillets.

| Work | Owner | Review |
|---|---|---|
| 3.0A and 3.0B (one owner, in order) | Fable | the dot (GPT Astra), every PR |
| Drawings and PMI | Fable; moves to Claude Opus if Fable quota runs short | the dot |
| App lane | Claude Opus, or GPT 6.1 or 6.2 if available | the dot |
| ROADMAP.md update | firstmate | owner approval |

## Rules that stay fixed between checkpoints

- The allowance target is 20M. Until phase 5 the application default is 100M.
- Proven, bounded and refused answers stay separate inside the kernel. The user sees results as approximate until the proof button, Q14. The app removes proved and bounded lights, retaining the refused message and working indicator, Q31.
- Kernel changes are general. No allowance or tolerance is loosened to make one case pass.
- Clock time never fails a change before phase 5.
- The app shows Coming soon only for scope written on this roadmap.

[Checkpoint 2 working rules, section 4.9](docs/checkpoints/checkpoint-2-decisions.md#49-working-rules-r1-to-r5-with-the-r2-tuning), Q4 and Q24.

- R1: one feature, one owner across kernel and app; parallel lanes only on different features.
- R2: a kernel change to an app-used feature passes that app acceptance, with one later feature on top, before merge. Tuning (Q24): each SolidDraft scenario lists the kernel parts it uses; a Catallax pin bump runs only the scenarios whose kernel parts changed, plus the bracket scenario as a smoke test; the full set runs nightly and before every release; shared core changes or scenarios with no list run everything.
- R3: a target that fails two different approaches is parked and brought to the owner.
- R4: once a release is called, only wrong-result, crash and broken-feature fixes go in.
- R5: scoreboard and census once per checkpoint, no dry runs; coverage counted from `docs/COVERAGE.md`.

### Tolerance rules T0 to T8

[Checkpoint 2 tolerance specification, section 4.3](docs/checkpoints/checkpoint-2-decisions.md#43-tolerance-rules-t0-to-t8), Q7.

| Rule | Content |
|---|---|
| T0 | Exact entities store no tolerance and have no size limit |
| T1 | Floor for f64 geometry: eps(R) = max(1e-7 mm, 1e-11 x R), R the coordinate size in the frame the operation runs in (about 1e-5 mm at 1 km). Warn when one body exceeds 1 km. Operations run in local frames reached by exact rational transforms |
| T2 | Own fitted geometry: fit to max(1e-6 mm, 10 x eps), tightened to L/100 near small features; cap max(1e-5 mm, 100 x eps). Today it is 1e-7 mm |
| T3 | Imports: measure every entity; the file's declared value is only a hint. Automatic cap 0.01 mm for edges and faces, 0.02 mm for vertices. Above that, refit first, then the user may opt in up to 0.1 mm |
| T4 | Margin: 1.1 x (measured gap + measurement error), applied only to fresh measurements, never stacked. Where three or more surfaces meet, build the vertex first |
| T5 | Pair rule: for two entities that share no topology, error if t_a + t_b > d/2, warn if t_a + t_b > d/5. Found with a bounding-box tree inflated by 5t, so a part with no small features needs no distance computation |
| T6 | Short edges and slivers are invalid only through T5 |
| T7 | After each operation: refit the 2D copy, measure again, lower where possible. Warn if a tolerance rises more than 10x above its inputs |
| T8 | Export with a receiver profile (for example NX at 0.01 mm, CATIA at 0.001 mm): refit edges above the receiver tolerance, warn on features under 10x and 4x the receiver tolerance, write the largest entity tolerance as the file uncertainty, offer scaled export for micro parts |

The cap that applies is always the smaller of the class cap and the local term. On a part with no small features, the class cap governs.

## History

### Phase 2 lane setup, superseded at Checkpoint 2

**How phase 2 runs.** The owner chose three kernel lanes from the start at Checkpoint 1, so the longest chain, 2.2 then the rest of 2.3 then 2.5, starts on the first day. Lanes run side by side because no lane's result changes another lane's current job.

- **First, the measures.** A short job before any lane needs them. The census-part hole is placed so it removes material in the reference kernels, P8, P9 and the curved-housing chain join the scoreboard as targets, and the definitions then stay fixed to Checkpoint 2.
- **Lane 1, spline.** Finish 1.3's three loft top cuts at the 100M application default, then 2.1.
- **Lane 2, import.** First rank the measured failure causes by files gained: export gaps (files that import but do not export valid) against first-step import refusals, from the Checkpoint 1 census data. Then stitch imported surface bodies into checked solids, which is 2.3's "faces sew into a solid" moved to the front. Then work the causes in ranked order, and take 2.3's trim, extend and thicken once 2.2 lands.
- **Lane 3, new tools.** 2.2, the remaining 2.3 surface tools, then 2.5. General sweeps (2.4) remain on this lane.
- **App, alongside, on Claude Opus.** The Cancel fix first, then surface bodies shown and clearly marked, then the import health report, then each kernel item as it lands.

Kernel lanes ran on Codex models and the app lane on Claude Opus. Codex cloud reviewed and advised, without authoring PRs. A review held a merge only on new kernel algorithms, numerical admission or tolerance rules, and data-loss or release paths; everything else merged on green checks and was reviewed after. Checkpoint 2 replaces this setup with the lane table above and the dot reviewing every 3.0A and 3.0B PR from code and research, not the author's summary.


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
| 15 to 16 Nov | Phase 1. Numerical bodies under stated tolerances with approximate STEP. STEP import by tier, per-entity tolerances and source topology healing: 418 complete solid imports and 384 valid exports of 3,000 in the tier census. Spline loft Booleans: 11 of 14 recorded requests reach STEP. Mixed fillet corners, partial-oblique fillet ends, hollow drafts, three-edge chamfer corners, shells with side openings, skew holes through bevelled collars: P1 to P5 through SolidDraft. Body handover between kernel instances. Browser kernel worker. Verified STEP file replacement. First edit-chain scoreboard. Checkpoint 1. Catallax 0.3.0, SolidDraft 0.6.0. |

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
| 15 Nov | Until phase 5 the application's default allowance is 100M. 20M stays the target, and every edit over it is measured per stage and listed for phase 5. | Slow parts still build in SolidDraft while capability comes first. |
| 16 Nov | The owner moved Checkpoint 1 before the spline fillets. Milestone 1.5 is now 2.1, and the curved-housing chain moves with it. | Reach the first scoreboard sooner. |
| 16 Nov | Stitch imported surface bodies into checked solids early in phase 2. | Import is the largest gap, and 2,202 census files open only as unchecked surface sets. |
| 16 Nov | Phase 2 runs as three kernel lanes from the start. The import lane ranks measured failure causes before choosing its next batch. | Start the longest chain on the first day. Work where the files are. |
| 16 Nov | The delivery rule of large PRs, minimal tests and no speed work stands until phase 5. | It finished phase 1, and its Later notes became the phase 5 list. |
| 16 Nov | The 10 Nov port rule widens for 2.2: SMLib's surface-construction code may be adapted, with attribution, and Catallax's own checks decide what ships. | The NURBS toolkit is spline-only work. |
| 16 Nov | Checkpoints stay undated, and the scoreboard and census run only at checkpoints. When model quota runs short, the owner decides at the time. | The owner's answers at Checkpoint 1. |

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
| [SMLib intersection port scope](docs/research/smlib-intersection-port-scope.md) | Adopted | 37,255 physical lines of SMLib left to port before a plane cut. SMLib's tracer missed a tiny loop, so Catallax keeps its own complete search. |
| [Spline against spline by the 2026 method](docs/research/spline-intersection-2026-method.md) | Adopted with a proof per cell | All 24 test pairs resolve inside 20M: 20 proved, 4 tangent contacts refused by name. |
| [Spline fillets: the offset search](docs/research/spline-fillet-offset-search.md) | Route changed | The direct search hit its depth limit. The replacement rounds one loft edge within 59 nm. |
| [Tangent-touch prototype](docs/research/tangent-touch.md) | Parked | Interior touches classify. A touch on a border is undecided. First candidate for a bounded answer. |
| [Census of 3,000 public STEP files](docs/research/step-census.md) | Recorded | 369 files (12.3%) valid after OCCT builder construction, including possible internal pcurve healing, beside 2,701 OCCT healed-valid (90.0%). Conditioning: 400 admitted, 369 valid after builder construction (92.25%), 14 invalid, 15 export-refused, 2 timeouts. Original 186 (6.20%) is admission-only. Same seed, source hashes and denominator. |
| [Rational support projection qualification](docs/research/step-census-opened/README.md#rational-support-projection-qualification) | Measured, no file gain | Determinant/trace bounds certify a wider patch class, but the census stays at 550 complete solid files and 516 builder-valid exports. Six initial refusals change; later incidence and resource gates remain. |
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

## Later planning: SolidDraft feature gaps

The owner put all eleven proposals in Later. They are candidates for checkpoint review. Their letters identify the research proposals, and the sizes are rough audit estimates. Starting work requires the owner's decision at a checkpoint.

| Proposal | Use case | Rough size | Existing roadmap scope |
|---|---|---|---|
| [A. Hole stacks and callouts](docs/research/feature-gap/report.md#proposal-a) | Edit a compound fastener seat as one hole and carry its stages into a drawing callout. | M · 2–4 engineer-weeks | [1.4](#everyday-edge-and-wall) delivered ordinary holes. [3.0B](#one-model-and-predictor) covers continuation. [Drawings and PMI](#drawings-and-pmi) owns the annotation pipeline. |
| [B. Grid and adaptive feature patterns](docs/research/feature-gap/report.md#proposal-b) | Lay out holes or pockets with two pitches, skip members and replay the seed where supports vary. | M · 2–4 engineer-weeks for a grid, L with replay | [C5](#drawings-pmi-acceptance) covers recognition of bracket holes on import. The authoring grid and replay are Later. |
| [C. Slot, pocket and profile cuts](docs/research/feature-gap/report.md#proposal-c) | Edit a keyway, guide slot or radiused-floor pocket through its profile and termination. | M · 2–4 engineer-weeks for straight cuts, L for general sweeps | [Phase 2 App](#phase-2-app) owns Solid Slot. [2.4 and Coverage #7](#coverage-general-sweeps) cover sweep and termination breadth. [Coverage #5](#coverage-spline-booleans-and-fillets) and [#8](#coverage-edge-and-wall) cover spline stock and floor-round dependencies. Later planning covers the expanded profile-cut workflow. |
| [D. Revolved cuts and groove presets](docs/research/feature-gap/report.md#proposal-d) | Place a retaining-ring seat or seal gland and move it from an axis or face datum. | M · 2–4 engineer-weeks | [Coverage #8](#coverage-edge-and-wall) covers edge treatments used on groove rims. |
| [E. Edge finishing that survives edits](docs/research/feature-gap/report.md#proposal-e) | Adjust a housing's rounds and chamfers after a recess or hole changes while keeping surviving references. | L · 1–2 engineer-months per major slice | [2.1](#spline-fillets), [3.0B](#one-model-and-predictor), [3.4](#smlib-modelling-tools) and Coverage [#3](#coverage-bend-twist-and-blends), [#5](#coverage-spline-booleans-and-fillets) and [#8](#coverage-edge-and-wall) own kernel blends, chamfers and continuation. [4.2/4.2a](#field-modelling) owns field fillets. Later planning covers app controls and previews beyond those commitments. |
| [F. Rib and gusset](docs/research/feature-gap/report.md#proposal-f) | Stiffen a bracket or housing with a repeatable support that follows adjacent faces. | M · 2–4 engineer-weeks for planar ribs | [Coverage #8](#coverage-edge-and-wall) covers subsequent rib finishing. |
| [G. Shell and enclosed cavity controls](docs/research/feature-gap/report.md#proposal-g) | Lighten an enclosure while inspecting its walls, openings and sealed cavities. | L · 1–2 engineer-months for first curved coverage | [Coverage #8](#coverage-edge-and-wall) owns wider shell geometry. [Phase 4 App](#phase-4-app) covers field-driven thickness. Later planning covers cavity controls and section previews. |
| [H. Sheet metal and flat patterns](docs/research/feature-gap/report.md#proposal-h) | Produce a folded bracket, its flat DXF and bend directions from one feature history. | XL · 2–4+ engineer-months, staged | [3.4](#smlib-modelling-tools) and [Coverage #3](#coverage-bend-twist-and-blends) cover general bend and twist geometry. [Coverage #8](#coverage-edge-and-wall) supports edge and wall geometry. The bend graph and flat pattern are Later. |
| [I. Threads and coils](docs/research/feature-gap/report.md#proposal-i) | Specify tapped holes for manufacture, model printable threads and check coil fit. | M · 2–4 engineer-weeks for thread intent, XL for geometry | [2.4 and P8](#general-sweeps), [Coverage #7](#coverage-general-sweeps) and [5.1](#speed-and-simplification) own helix delivery and cost follow-up. Later planning covers thread intent and thread/coil controls. |
| [J. Guided loft and connected surfaces](docs/research/feature-gap/report.md#proposal-j) | Shape a duct along editable rails and repair and close its surfaces for manufacturing. | XL · 2–4+ engineer-months, staged | [2.2](#nurbs-toolkit), [2.3](#surface-bodies), [2.1](#spline-fillets) and Coverage [#1](#coverage-guided-loft), [#5](#coverage-spline-booleans-and-fillets) and [#6](#coverage-surface-bodies) own guided construction, surface operations and spline finishing. Later planning covers rail and section controls and diagnostics. |
| [K. Dedicated spur gear generator](docs/research/feature-gap/report.md#proposal-k) | Generate an involute spur gear for a prototype transmission with editable tooth size and backlash. | L · 1–2 engineer-months for external spur gears | Dedicated gear authoring is Later. |

## Further reading

Existing: `VISION.md` in both repositories, `docs/ARCHITECTURE.md` for dated decisions, `docs/CAPABILITIES.md` for what works today, `docs/TARGET-PARTS.md` for P0 to P10 and when a part counts as built, `docs/milestones/` for milestone notes and the STEP census, `docs/research/` for the files linked above, and `docs/archive/PLAN-2026-11.md` for the plan this file replaces. [docs/COVERAGE.md](docs/COVERAGE.md) records the 61 lines and bars; [Checkpoint 2 decisions](docs/checkpoints/checkpoint-2-decisions.md) records the adopted direction and acceptance specifications.

To be written, each short:

- `docs/LEFTOVERS-0.2.0.md`: everything the release sprint postponed, with its recommendation.

### Old milestone names

| Old name | Now |
|---|---|
| M1 Catallax 0.1.0 alpha | Reached with Catallax 0.1.0. See the history for 25 to 29 Oct. |
| M2 to M5 The first part in SolidDraft, everyday features, distance features, general fillets and draft | 1.4 Everyday edge and wall features |
| M6 The roadmap and the decisions written down | Reached. This file replaces its plan. |
| M7 Lofts and spline sweeps finished | Reached. See the history for 10 and 11 Nov. |
| M8 Cuts after a fitted fillet | Deferred to 2.1 Fillets between spline faces |
| M9 Spline intersections ported | Part of 1.3 Cuts and joins on spline faces |
| M10 Catallax cuts and joins spline faces | 1.3 Cuts and joins on spline faces |
| M11 Fillets between spline faces | 2.1 Fillets between spline faces |
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
| M29 NURBS toolkit | 2.2 NURBS toolkit |
| M30 Surface bodies | 2.3 Surface bodies |
| M31 General sweeps | 2.4 General sweeps |
| M32 Direct edits | 2.5 Direct edits |
| M33 Drawing views from the kernel | 3.1 Drawing views from the kernel |
| M34 Assemblies and exchange | 3.2 Assemblies and exchange |
| M35 Tools only SMLib has | 3.4 Tools only SMLib has |
