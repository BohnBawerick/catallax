# Checkpoint 1

Status: final. The cut run measured Catallax main `ebc42164` on 16 Nov 2026, once for the scoreboard and once for the census. The dry-run numbers stay beside them in their own sections.

This is step 2 of the checkpoint in `ROADMAP.md`: gather the results. The owner then picks the order of phase 2.

## How to rerun

Both runs are one command each. Run them on the work PC through `cx-offload`, from a clean checkout of the commit being measured. The census uses two import workers and kills any one observer above 6,144 MiB, recording that file as out of memory by name, so one huge source cannot exhaust the shared machine. The largest known source, 00005642, peaks at 4,605 MiB on its own.

```sh
TOOLS=~/cx-checkpoint1/tools        # pinned CadQuery, smlibpy and OCP install, made on first use
SOURCES=~/cx-census/selected-step   # private census corpus
SELECTION=~/cx-census/selection-private.json
python3 -B tools/checkpoint/checkpoint.py scoreboard --tools $TOOLS --sources $SOURCES --selection $SELECTION --out RUN/scoreboard
python3 -B tools/checkpoint/checkpoint.py census --tools $TOOLS --sources $SOURCES --selection $SELECTION --out RUN/census
```

After this cut, the scoreboard command also runs the phase-2 targets (P8, P9 and the curved-housing chain), and it moves a census hole off-centre when the centred hole removes nothing. [tools/checkpoint/phase2](../../tools/checkpoint/phase2/README.md) defines both. The cut-run tables below used the earlier command.

Each command builds its own binaries and writes `scoreboard.md` or `census.md` beside its raw data. Wall time on the work PC at the cut: 14.2 min for the scoreboard and 13.0 min for the census, on a warm target directory. A cold release build adds about 4 min to each command. Installing the tools the first time adds about 5 min.

## Coverage count

55% of the 61 capability lines, from the 14 Nov release count. Nobody has recounted since. The count says what exists, not how many real parts finish, so the scoreboard and census below are the measures of phase 1.

## Edit-chain scoreboard

Targets: P0 (eight hard-case recipes), P1 to P5, and the ten real census parts. P6 to P10 have no comparison definition yet ([P06 input pin](../phase1/P06/README.md)).

Each target is built, then rebuilt with one committed scalar edit, in Catallax, OpenCASCADE (CadQuery 2.6.1 on OCCT 7.8.1) and SMLib (smlibpy 0.4.0). Each kernel runs the whole target script both times. Each kernel uses its own released import for the census parts.

A census part's chain is: import the STEP file, drill one hole, then change the hole radius. The hole runs along the smallest bounding-box extent through the centre of the other two. Its radius is a tenth of the smaller remaining extent, and the edit makes it 1.5 times larger. `tools/checkpoint/cutter.py` takes the box from OCCT's reading of the source and rounds every value to 0.001 mm before any kernel runs, so all three get the same numbers. CadQuery and smlibpy read the STEP file with their default settings, which include healing; Catallax reads it without healing.

**One definition of done.** The original and the edited target both build. Both STEP exports open in OCCT 8.0.1.1 as valid closed solids with the kernel's own face count: at least one solid, every shell closed, and no face outside a solid. No independent check found a wrong value. On a census part, the drill and the radius change must each lower the volume of the exported solid, measured in OCCT, by more than a millionth of the source volume. Volume, area, centroid and the 216-point membership grid are checked only where `head-to-head/analytic.py` has a formula. The census parts have none, so the STEP check decides them.

The edit chains are one edit long. No longer chain has an owner definition yet (P06 found none), so this run does not invent one.

### Cut run

Main `ebc42164` with this tooling (head `86e71794`), work PC, 16 Nov 2026 13:03 to 13:18 UTC. 100M work per Catallax phase, 120 s per engine, one warmup and one measured sample. Wall time 14.2 min, release build included.

| Target | Catallax | OpenCASCADE | SMLib | Catallax work (build, edit) | Catallax names after edit |
| --- | --- | --- | --- | ---: | --- |
| P0 block and sphere | done | done | done | 633,934, 397,086 | kept |
| P0 crossing pipes | done | done | edit fuse failed | 3,169,337, 1,997,714 | kept |
| P0 equal cylinder and sphere | done | done | done | 2,813, 76,393 | kept |
| P0 general cone | build refused | done | done | 2,341, 2,341 | not run |
| P0 general torus | build refused | done | done | 4, 4 | not run |
| P0 line tangent | done | done | build fuse failed | 35,038, 39,208 | kept |
| P0 near miss | done | done | build fuse failed | 3,625, 3,590 | kept |
| P0 seam | done | done | done | 296,359, 280,564 | kept |
| P1 bracket | done | done | done | 13,618,847, 13,479,812 | kept |
| P2 collar | done | done | wrong volume and area | 1,345, 1,382 | kept |
| P3 filleted boss | membership check failed | done | done | 20,773, 21,591 | kept |
| P4 shelled tray | done | done | done | 917,132, 936,262 | kept |
| P5 flanged housing | done | fillet failed | fillet failed | 1,058,398, 1,057,839 | 16 unnamed |
| ABC 00000131 | drill removed nothing | drill removed nothing | drill raised the volume | 90,211, 90,131 | kept |
| ABC 00000173 | build refused | done | done | 1,989,712, 1,989,922 | not run |
| ABC 00000190 | done | done | done | 31,083, 30,963 | kept |
| ABC 00000405 | done | done | done | 415,413, 410,612 | kept |
| ABC 00000511 | done | done | done | 30,403, 30,283 | kept |
| ABC 00000605 | done | done | done | 27,077, 27,250 | kept |
| ABC 00000632 | drill removed nothing | drill removed nothing | drill raised the volume | 190,365, 190,365 | kept |
| ABC 00000633 | drill removed nothing | drill removed nothing | radius change removed nothing | 104,183, 104,705 | kept |
| ABC 00000695 | done | done | done | 86,583, 83,250 | kept |
| ABC 00000821 | done | done | drill raised the volume | 218,037, 218,237 | kept |

| Kernel | Done, P0 to P5 (13) | Done, census parts (10) | Done in the published 0.2.0 comparison (13) |
| --- | ---: | ---: | ---: |
| Catallax | 10 | 6 | 6 |
| OpenCASCADE (CadQuery 2.6.1) | 12 | 7 | 12 |
| SMLib (smlibpy 0.4.0) | 8 | 6 | 8 |

The last column applies the same rule to the [published comparison](../../site/comparison-results.json) at `3852c773`. Since then P0 crossing pipes, P1, P4 and P5 finish in Catallax. Every Catallax build used at most 13.6M work, under the 20M target.

The only verdict that moved since the dry run is Catallax's crossing-pipes edit, which now finishes inside the 120 s per-engine cap. Catallax's work counts are unchanged on P0 to P5. On the census parts they rise by 110 to 278 work.

### Dry run

Dry run on main `c52fd8f4` with this tooling, work PC, 100M work per Catallax phase, 120 s per engine, one warmup and one measured sample.

| Target | Catallax | OpenCASCADE | SMLib | Catallax work (build, edit) | Catallax names after edit |
| --- | --- | --- | --- | ---: | --- |
| P0 block and sphere | done | done | done | 633,934, 397,086 | kept |
| P0 crossing pipes | edit display timed out | done | edit fuse failed | 3,169,337, 1,997,714 | not run |
| P0 equal cylinder and sphere | done | done | done | 2,813, 76,393 | kept |
| P0 general cone | build refused | done | done | 2,341, 2,341 | not run |
| P0 general torus | build refused | done | done | 4, 4 | not run |
| P0 line tangent | done | done | build fuse failed | 35,038, 39,208 | kept |
| P0 near miss | done | done | build fuse failed | 3,625, 3,590 | kept |
| P0 seam | done | done | done | 296,359, 280,564 | kept |
| P1 bracket | done | done | done | 13,618,847, 13,479,812 | kept |
| P2 collar | done | done | wrong volume and area | 1,345, 1,382 | kept |
| P3 filleted boss | membership check failed | done | done | 20,773, 21,591 | kept |
| P4 shelled tray | done | done | done | 917,132, 936,262 | kept |
| P5 flanged housing | done | fillet failed | fillet failed | 1,058,398, 1,057,839 | 16 unnamed |
| ABC 00000131 | drill removed nothing | drill removed nothing | drill raised the volume | 90,059, 89,979 | kept |
| ABC 00000173 | build refused | done | done | 1,989,546, 1,989,756 | not run |
| ABC 00000190 | done | done | done | 30,847, 30,727 | kept |
| ABC 00000405 | done | done | done | 415,261, 410,460 | kept |
| ABC 00000511 | done | done | done | 30,167, 30,047 | kept |
| ABC 00000605 | done | done | done | 26,841, 27,014 | kept |
| ABC 00000632 | drill removed nothing | drill removed nothing | drill raised the volume | 190,157, 190,157 | kept |
| ABC 00000633 | drill removed nothing | drill removed nothing | radius change removed nothing | 103,905, 104,427 | kept |
| ABC 00000695 | done | done | done | 86,347, 83,014 | kept |
| ABC 00000821 | done | done | drill raised the volume | 217,927, 218,127 | kept |

| Kernel | Done, P0 to P5 (13) | Done, census parts (10) | Done in the published 0.2.0 comparison (13) |
| --- | ---: | ---: | ---: |
| Catallax | 9 | 6 | 6 |
| OpenCASCADE (CadQuery 2.6.1) | 12 | 7 | 12 |
| SMLib (smlibpy 0.4.0) | 8 | 6 | 8 |

The last column applies the same rule to the [published comparison](../../site/comparison-results.json) at `3852c773`. Since then P1, P4 and P5 finish in Catallax. Every Catallax build used at most 13.6M work, under the 20M target.

The census rows come from a rerun of the ten census parts with the material rule. The P0 to P5 STEP files of the dry run were checked again under the closed-solid rule, which counts the distinct faces owned by all solids together, and no verdict changed. Line tangent and near miss pass with two closed solids each in Catallax and OCCT, as does SMLib's general torus.

On 00000131, 00000632 and 00000633 the hole lands inside an existing bore, so Catallax's and OCCT's volumes do not change at all. These three edits test nothing as defined. Picking a hole position that cuts material, for example by trying a few offsets from the centre, would make them count; it is not done here.

SMLib's exported solids measure 0.28% to 0.33% more volume in OCCT than the source file, after a subtraction. That is why its rows read "raised the volume". The cause is not diagnosed.

Names are Catallax's own check: every face origin and edge name of the original result resolves after the edit. The other two kernels have no persistent naming contract, so they have no names column. P5 builds and passes its STEP check, but 16 of its faces or edges carry no name, so its name check fails.

## STEP census

The same 3,000 ABC files, seed and source hashes as the [census](../research/step-census.md). Every file goes through import on main. Every full admission is exported and checked in OCCT. A file counts only when every imported body exports and passes `BRepCheck_Analyzer` after OCCT builder construction, which can heal pcurves internally.

### Cut run

Main `ebc42164` with this tooling (head `86e71794`), work PC, 16 Nov 2026 13:18 to 13:31 UTC. The 20M import and export allowances, a 60 s cap per file, two import workers, two export workers and a 6,144 MiB cap per observer. Wall time 13.0 min: 7.9 min import, 4.9 min export, a few seconds of build. The run is complete, and no observer reached the memory cap.

| Measure | Files | Share of 3,000 |
| --- | ---: | ---: |
| Catallax, valid after OCCT builder construction | 374 | 12.47% |
| OpenCASCADE, healed-valid | 2,701 | 90.03% |
| Catallax, [previous capture](../research/step-census-circular-pcurves/README.md) | 369 | 12.30% |

Five files are new (00000174, 00004152, 00004399, 00004452, 00009819). None of the 369 is lost, and each of them has a measured verdict in this run.

| Admission | Files |
| --- | ---: |
| Full | 405 |
| Partial (some source items not verified) | 3 |
| Refused | 2,590 |
| Timeout | 2 |

The 408 admitted files split into 190 certified and 218 numerical approximate imports. Of the 405 full admissions, 374 export valid, 14 export invalid bodies, 15 refuse at `step.numerical-gap`, one refuses at `resource_limit` and one times out.

Not covered, with no verdict in this run: 00000187 and 00008201 time out at import, and 00004777 times out at export. They are not passes, losses or refusals. No run ended for memory or by a process error.

The only change from the dry run: 00002204, which timed out at export in the dry run, now refuses at `resource_limit` after 19,586,631 of its 20M export work. Every count in the first-refusal table below is unchanged.

### Dry run

Dry run on main `c52fd8f4`, the 20M import and export allowances, a 60 s cap per file, four import workers and two export workers.

| Measure | Files | Share of 3,000 |
| --- | ---: | ---: |
| Catallax, valid after OCCT builder construction | 374 | 12.47% |
| OpenCASCADE, healed-valid | 2,701 | 90.03% |
| Catallax, [previous capture](../research/step-census-circular-pcurves/README.md) | 369 | 12.30% |

Five files are new (00000174, 00004152, 00004399, 00004452, 00009819) and none of the 369 is lost.

#669 changed both census observers after this dry run, so the census command ran once more on main `634a9483` to check it still works. It gave the same 374 valid files. That run was not complete: the observer for 00005642, a 465 MB source, reached 3.9 GB and the work PC's out-of-memory guard killed it while other jobs were running. Alone, the same file refuses at `resource_limit` after 18 s. The cut run has the work PC to itself and runs under the memory cap above.

| Admission | Files |
| --- | ---: |
| Full | 405 |
| Partial (some source items not verified) | 3 |
| Refused | 2,590 |
| Timeout | 2 |

The 408 admitted files split into 190 certified and 218 numerical approximate imports. Of the 405 full admissions, 374 export valid, 14 export invalid bodies, 15 refuse at `step.numerical-gap` and two time out.

| First refusal cause | Files |
| --- | ---: |
| `incidence_gap` | 728 |
| `import.numerical-unsupported-family` | 421 |
| `numerical.invalid-trim` | 330 |
| `resource_limit` | 222 |
| `unsupported:ADVANCED_FACE` | 168 |
| `planar_repair_not_verified` | 129 |
| `self_intersecting` | 119 |
| `self_intersection_not_verified` | 93 |
| `unsupported:EDGE_CURVE` | 76 |
| `not_solid` | 69 |
| `unsupported:TOROIDAL_SURFACE` | 66 |
| `unsupported:other_step_type` | 54 |
| `unsupported:BREP_WITH_VOIDS` | 40 |
| `unsupported:SURFACE_OF_LINEAR_EXTRUSION` | 31 |
| `malformed` | 28 |
| `unsupported:ELLIPSE` | 13 |
| `unsupported:VERTEX_LOOP` | 3 |

`tools/research/import_census.py` used to require every numerical import to state exactly 0.01 mm. On current main a few state slightly less, such as 0.009999999994 mm, so the collector now checks that every stated tolerance is at most the 0.01 mm cap and that the checked residual is within it.

A refusal is the first stop only. Fixing one cause can move a file to the next cause rather than admit it.

The OCCT figure is the published reference capture: normal transfer, default `ShapeFix_Shape`, then `BRepCheck_Analyzer`. Its tolerances differ from Catallax's, so the two numbers measure coverage, not a common tolerance.

The tolerant-import tiers landed in [#669](https://github.com/BohnBawerick/catallax/pull/669) during the checkpoint. This run measures the existing import entry, which reports two categories, certified and numerical approximate. The tier consumer has its own full census at producer `3437c18d2fcc59bbc74c1a4123bf443237f99f73`, measured by the import lane and not rerun here ([step-census-opened](../research/step-census-opened/README.md)):

| Tier census measure, #669 at `3437c18d` | Files | Share of 3,000 |
| --- | ---: | ---: |
| Complete-source solid imports | 418 | 13.93% |
| Every expected solid exports and passes OCCT builder validity | 384 | 12.80% |
| Supported source openings, any tier (no checked B-rep or export authority) | 2,739 | 91.30% |

The openings count includes 119 mixed or partial files and 2,202 retained unverified source-surface sets. It is not a solid import or a valid export, so it stays out of the headline. #669 reports no loss against 405 complete solids and 374 valid exports.

## App

The browser app's kernel work now runs off the page's thread. #320, #321 and #322 landed; #322 moved browser builds, edits, Undo, Redo and opening a part into a kernel worker, whose built bodies the page adopts, with a "Building the part" bar and Cancel. On the bracket in headless Chromium with CI's flags, the longest page freeze went from 7.4 s to about 0.2-0.4 s on open, from 2.8 s and then 7 s to 0.2-0.5 s on a depth edit, and from 2.7-2.9 s to about 0 on Undo and Redo. #323, waiting on CI with its web build already green, makes the page take the worker's bodies without reading their exact geometry. That brings the depth edit to 65-190 ms. The first exact command after an edit (Measure, Check, export, section, sketch on a face) then waits behind a "Reading the part" bar for about 0.2-0.7 s. #327, green and waiting for review, repins Catallax to 8bf3999a with #657 and #650. Measured on the same app revision, #657 gives no gain on the bracket: 192 ms inside the hand-over read on the old pin and 199 ms on the new one, and every page timing within run-to-run spread. #650 lets the near-coincident flange use its early preview, so two display tests were flipped. Later: the viewport reads its pick back from the GPU every frame the pointer rests on it, which software rendering holds for 0.7-1 s. That keeps Cancel and the edit after it over 250 ms in CI's browser, so the 250 ms browser check stays out of the web job.

Since that was written, SolidDraft #327 (the repin to `8bf3999a`) and the app hand-over PRs have landed. #323 (lazy reads) was still merging at the cut; it does not change Catallax main.

P1 to P5 pass through SolidDraft on Catallax main `47fb0a2e`, run on the work PC on 16 Nov with the app's scenario runner. Each part builds, takes its committed edit, undoes and redoes, saves, reopens with the same geometry, and exports STEP from a fresh process. Every export, original, edited and reopened, is valid in OCCT with healing off. This is milestone 1.4's "done when" line. Timing budgets, which are advisory, were not met on P1 and P5.

| Part | Edit | STEP export work |
| --- | --- | ---: |
| P1 | slot dimension | 15.12M |
| P2 | chamfer 1 to 0.75 | 0.54M |
| P3 | junction fillet 2 to 2.5 | 0.93M |
| P4 | shell wall 2 to 1.5 | 1.25M |
| P5 | junction fillet 2 to 2.5 | 4.50M |

## Release

The checkpoint ends with a release of both repositories and a redeploy of the website and browser app on it, right after this report. Firstmate dispatches the release; this section lists what to cut. Both stay alpha prereleases.

**Catallax 0.3.0**, from the main commit pinned for the cut run, so the release is the commit these numbers measure. `docs/RELEASING.md` advances the minor version for a completed milestone step, and the phase-1 checkpoint is one. Since v0.2.0 it carries 78 merged changes:

- Real STEP import: numerical admission of planar, oblique, cylindrical, cone and B-spline revolution sources, shell models as typed surface bodies, per-entity tolerances, source topology healing, and opening by tier (#669). The census moves from 186 admitted files to 418 complete solid imports and 384 valid exports, against OCCT's 2,701.
- STEP export: approximate STEP for numerical bodies, curved patches with holes, cones, and spline Boolean results without an exact form (#668), with the seam, circle and pcurve repairs that took post-builder validity from 322 to 374 files.
- Spline Booleans: loft Booleans and finite cylinder sections (#664), plane contacts through spline Booleans, numerical spline-line and cone-cylinder queries, and finite spline Boolean queries that bring the recorded spline requests to eleven STEP answers out of 14 (#672).
- Everyday features: mixed fillet corners, partial-oblique fillet ends, hollow drafts, three-edge chamfer corners, shells with side openings, draft and shell around mixed openings, and skew holes through bevelled collars (#667).
- Kernel hand-over between instances (#632, #657), verified STEP file replacement with face and edge correspondence, and the 100M application allowance with a 20M target (#619).
- This report's tooling (#665).

**SolidDraft 0.6.0**, pinned to that same Catallax commit. Since 0.5.0 it carries:

- Browser kernel work off the page thread: builds, edits, Undo, Redo and opening in a kernel worker (#322), mass properties in a worker (#314), lazy reads of the adopted bodies (#323, merging at the cut; include it only if it has merged before the release).
- STEP: multi-solid file replacement as one verified transaction (#319), matched face names kept across replacement (#324).
- Modelling: part colours and transparency, Draft with its own pull direction (#326), supporting planes kept for face-mounted cuts, 100M work per request.
- Sketcher and interface: box and crossing selection, relations for two-entity selections, relation glyph hints, broken-dimension recovery, the left icon strip with the feature tree, walk-over redesigns R1 to R8, and eleven small fixes from a user walk.

The redeploy takes the browser app and website to SolidDraft 0.6.0 on Catallax 0.3.0. The release checklist's replay, evidence and package dry runs still apply to both; this list does not replace them.

## What still refuses and why

From the cut run. The P0 crossing-pipes edit, which timed out in display meshing in the dry run, now finishes.

Scoreboard, Catallax:

- **P0 general cone and general torus.** `intersect.traced-pair-required`: non-coaxial cone/cone and torus/torus intersections need certified branch tracing. Both reference kernels build them.
- **P3 filleted boss.** STEP is valid with the right face count, but OCCT's classifier puts one grid point outside the part (two after the edit) inside the exported solid. The published 0.2.0 comparison failed on the same points. Not yet diagnosed.
- **P5 flanged housing.** Done, but 16 faces or edges carry no name. The published comparison did not run this check on P5, because P5 did not build then.
- **ABC 00000173.** The drilled hole refuses at `boolean.no-rational-point`: no rational probe was certified strictly inside a face. Import alone answers.
- **ABC 00000131, 00000632, 00000633.** Not refusals: the hole as defined removes no material in any kernel, so these rows are not done for anyone.

Scoreboard, reference kernels, for scale: OCCT fails the P5 junction fillet. SMLib fails P5's fillet, the fuses at line contact in P0 line tangent, near miss and crossing pipes, and reports a P2 volume and area outside the check's allowance. Its drilled census parts 00000131, 00000632 and 00000821 export with more volume than the source.

Census, cut run: 2,590 files refuse at import, with first causes in the table above, and two time out. 31 full admissions do not export as valid solids: 14 invalid bodies, 15 `step.numerical-gap` export refusals, one `resource_limit` export refusal (00002204) and one export timeout (00004777).

Not rerun here: the 14 recorded spline Boolean requests of milestone 1.3. As measured in [#672](https://github.com/BohnBawerick/catallax/pull/672), eleven now write checked approximate STEP within the 20M allowance, and the three top requests still refuse, so milestone 1.3 stays open.

## Later: optimisation notes from merged PRs

Phase 1 does no speed work. These are the costs that merged PRs recorded, kept here for the checkpoint review and for phase 5.

- **Fillet-end vertices.** Exact-first spent 13.6M work rounding one fillet-end vertex, which OpenCASCADE computes numerically (roadmap lessons, 15 Nov).
- **Partial-oblique fillet ends.** Radius 0.75 uses 30,019,099 work, over the 20M target ([#616](https://github.com/BohnBawerick/catallax/pull/616)).
- **STEP export of curved faces.** The edited radius-0.75 mixed-corner fillet spent 8.5 s in STEP export, 99.2% of it in `step.face`, mostly cylinder and torus spans. Reusing Taylor control symmetry cut that to 7.4 s. The export uses only 807,063 work, so the clock cost sits in exact rational arithmetic under `fit::taylor_controls`, not in the work count ([#646](https://github.com/BohnBawerick/catallax/pull/646)).
- **Body handover.** P1 handover export uses 81,482,582 work (182 ms) and import 34,914,680 work (92 ms); the GUI bracket export uses 23,872,737. All three exceed the 20M target and finish under the 100M guard ([#632](https://github.com/BohnBawerick/catallax/pull/632)).
- **Face-mounted extrusions.** Keeping identity placements exact cut a captured SolidDraft cut from 19.0M to 14.6M work. Other placements may carry the same algebraic zero expressions ([#604](https://github.com/BohnBawerick/catallax/pull/604)).
- **Spline-line queries.** The first Cartesian-box pose trial went over 20M; transverse Bezier subdivision removed it. Box-based pruning elsewhere may have the same cost ([#618](https://github.com/BohnBawerick/catallax/pull/618)).
- **Import work limit.** 177 census files stopped at the 20M import allowance in the last admission capture. Four of them (00000011, 00000094, 00000205, 00000207) stop at 19,999,976 work ([#638](https://github.com/BohnBawerick/catallax/pull/638)). Optional curved source checks share one 1M prepayment, so a file with less than 1M left refuses there, and legacy planar intersection work is not fully metered ([census](../research/step-census.md)).
- **Loft numerical admission.** The below-contact loft takes 10,219,174 work for construction, source conversion and numerical admission together ([#618](https://github.com/BohnBawerick/catallax/pull/618)).
- **Tolerant import healing.** Index per-body healing logs in `tolerant/healing/checks.rs` so checks stop revisiting foreign bodies, and reuse verified source recognition in opening and replay for large parts. Both keep the stated allowance ([#669](https://github.com/BohnBawerick/catallax/pull/669), Later).
- **Small spline lofts.** Moved 1/1000-scale lofts reach sewing near 9.8M work, then spend the rest of the allowance isolating source roots; a finite-domain search could skip unrelated global roots. Near-contact 1/1000 translation exhausts in border certificates and cell arithmetic, where a bounded numerical query could replace the exact certificates ([#664](https://github.com/BohnBawerick/catallax/pull/664), Later).
- **STEP geometry source.** Cache decoded nodes in `catallax-import::geometry_source` so shared geometry is not converted again for each use ([#666](https://github.com/BohnBawerick/catallax/pull/666), Later).
- **Numerical body handover.** Intern repeated checked definitions and evidence in the `NumericalBody` handover. The initial and edited encodes cost 24.18M and 48.26M work, the decodes 10.57M and 21.05M ([#667](https://github.com/BohnBawerick/catallax/pull/667), Later).
- **Per-body import charts.** Prepare derivative expressions once per chart, not once per fitted cell ([#670](https://github.com/BohnBawerick/catallax/pull/670), Later).
- **Slow steps at the 0.5.0 cut.** Loft 35 s, fillet on a loft edge 26 s, spline sweep 13 s, P5 fillet edit 4 to 6 s, field mass 3.5 to 7.7 s (roadmap 5.1).



## Proposed order of phase 2

The owner picks. This is a proposal with one reason each, using the numbering on main since [#662](https://github.com/BohnBawerick/catallax/pull/662) moved the spline fillets to 2.1.

1. **2.1 Fillets between spline faces, with the curved-housing chain.** It carries the two promises phase 1 handed on: the multi-patch fillet and the curved housing of 1.1. Nothing in these numbers outranks it, though neither the scoreboard nor the census measures it yet; P9, the lofted duct, still has no recipe.
2. **2.2 NURBS toolkit.** Three "not started" coverage lines sit here, and 2.3 and 2.5 need its tools.
3. **2.3 Surface bodies.** Under #669, 2,202 census files open only as retained, unverified source-surface sets, and OCCT itself gets only an open shell from 44. Sewing and thickening are what turn those into parts.
4. **2.4 General sweeps.** It depends on neither 2.2 nor 2.3, so it can run beside them; P8 is the only target part it unlocks.
5. **2.5 Direct edits.** It needs the trim, extend and sew of 2.2 and 2.3, and it pays only once imports open, which the census says most files still do not.

The numbers argue for one addition beside this order. Import is still the largest gap: 12.47% of files reach a valid export in this run, and 12.80% in #669's tier census, against OCCT's 90.03%. The four largest first refusal causes, `incidence_gap`, `import.numerical-unsupported-family`, `numerical.invalid-trim` and `resource_limit`, cover 1,701 files, and all four are 1.2 import work rather than a phase-2 line. I propose that import repair keeps running beside 2.1 as part of the resolve list, so it does not stop at the checkpoint.
