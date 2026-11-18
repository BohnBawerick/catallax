# Checkpoint 2

Status: reconciliation in progress. Numbers are draft until the cut run after the SolidDraft catch-up lands. This report gathers results; the owner chooses the next order. It does not close milestones whose acceptance still refuses.

## How to rerun

Run both commands from a clean checkout, through the build offload helper. Keep at least 20 GB free and remove the task's build output afterward. Pass the comparison installation, census sources and frozen selection as runtime paths. No source files belong in this repository.

```sh
python3 -B tools/checkpoint/checkpoint.py scoreboard --tools "$TOOLS" --sources "$SOURCES" --selection "$SELECTION" --out "$RUN/scoreboard"
python3 -B tools/checkpoint/checkpoint.py census --tools "$TOOLS" --sources "$SOURCES" --selection "$SELECTION" --out "$RUN/census"
```

The scoreboard retains the phase-2 definitions fixed in [the measures PR](https://github.com/BohnBawerick/catallax/pull/677). It runs P0 to P5, P8, P9, the curved housing and ten census parts against OpenCASCADE and SMLib. P6 and P7 have no comparison definition in this command. They are not silently counted as passing. The three formerly empty census drills now use the committed material-cutting offsets.

The census command now calls the existing `opened_import_census.py` build and collector. It measures the default per-body importer, including tolerant bodies and selected whole-source export, rather than rerunning the legacy all-or-nothing import entry. Two workers, 20M opening work, aggregate 20M export work, a 60-second observation cap and 6,144 MiB per observer remain unchanged. The collector binds the same 3,000 source hashes and selection. Each command writes its Markdown table and raw receipts outside the checkout.

## Coverage count

The last capability count remains 55% of 61 lines, from 14 Nov. This report does not turn merged APIs into a new percentage without the line-by-line recount. The scoreboard and census measure completed parts and files, not the capability percentage.

## Edit-chain scoreboard

One definition of done, unchanged from Checkpoint 1: the original and edited targets both build and their STEP exports open as valid closed solids in OCP 8.0.1.1, with the producer's face count and no independent wrong-value verdict. On census parts both the drill and radius edit must lower OCCT-measured volume by more than a millionth of source volume. Names are a separate Catallax check, not a gate invented for the reference kernels.

Catallax gets 100M work per phase. Each engine gets 120 seconds, with one warmup and one measured sample. The curved housing also reports the raw healing-off read separately; that extra read does not change the fixed scoreboard gate.

### Dry run, draft numbers

Main `403081f8`, 18 Nov 2026. The unchanged scoreboard completed end to end in 898.91 seconds, 15.0 minutes, including its build. These are not final cut numbers.

| Target | Catallax | OpenCASCADE | SMLib | Catallax work, build / edit | Names |
| --- | --- | --- | --- | ---: | --- |
| P0 block and sphere | done | done | done | 633,934 / 397,086 | kept |
| P0 crossing pipes | done | done | edit fuse failed | 3,169,337 / 1,997,714 | kept |
| P0 equal cylinder and sphere | done | done | done | 2,813 / 76,393 | kept |
| P0 general cone | build refused | done | done | 2,341 / 2,341 | not run |
| P0 general torus | build refused | done | done | 4 / 4 | not run |
| P0 line tangent | done | done | fuse failed | 35,038 / 39,208 | kept |
| P0 near miss | done | done | fuse failed | 3,625 / 3,590 | kept |
| P0 seam | done | done | done | 296,359 / 280,564 | kept |
| P1 bracket | done | done | done | 13,618,847 / 13,479,812 | kept |
| P2 collar | done | done | volume and area failed | 1,345 / 1,382 | kept |
| P3 filleted boss | membership check failed | done | done | 20,773 / 21,591 | kept |
| P4 shelled tray | done | done | done | 917,132 / 936,262 | kept |
| P5 flanged housing | done | fillet failed | fillet failed | 1,063,025 / 1,062,466 | failed |
| P8 spring and threaded boss | threaded union refused | done | done | 55,626,285 / 55,626,285 | not run |
| P9 lofted duct | outer loft refused | done | done | 4,793,088 / 4,793,088 | not run |
| Curved housing | spline cut refused | done | done | 1,511,266 / 1,511,274 | not run |
| ABC 00000131 | done | done | drill raised volume | 469,022 / 479,854 | kept |
| ABC 00000173 | drill refused | done | done | 1,989,712 / 1,989,922 | not run |
| ABC 00000190 | done | done | done | 31,083 / 30,963 | kept |
| ABC 00000405 | done | done | done | 415,413 / 410,612 | kept |
| ABC 00000511 | done | done | done | 30,403 / 30,283 | kept |
| ABC 00000605 | done | done | done | 27,077 / 27,250 | kept |
| ABC 00000632 | done | done | done | 697,371 / 1,128,070 | kept |
| ABC 00000633 | done | done | done | 112,187 / 113,026 | kept |
| ABC 00000695 | done | done | done | 86,583 / 83,250 | kept |
| ABC 00000821 | done | done | drill raised volume | 218,037 / 218,237 | kept |

| Kernel | Done, P0 to P5 of 13 | Done, phase-2 targets of 3 | Done, census parts of 10 |
| --- | ---: | ---: | ---: |
| Catallax | 10 | 0 | 9 |
| OpenCASCADE, CadQuery 2.6.1 | 12 | 3 | 10 |
| SMLib, smlibpy 0.4.0 | 8 | 3 | 8 |

The three corrected drill locations explain the census-chain gain from Checkpoint 1's 6/10 to 9/10 in Catallax; this is a repaired measure, not three newly implemented kernel families. P0 to P5 remain 10/13. Both reference kernels complete the three new targets, but their curved-housing exports fail the separately reported healing-off read at 1e-7 mm.

## STEP census

The headline is complete-source, occurrence-preserving, all-solids-valid export. An isolated valid prototype is not a whole-file pass. The report also keeps any observed valid solid, prototype-only validity and incomplete observation separate.

The reference remains OpenCASCADE's 2,701 healed-valid files of 3,000, 90.03%, from normal transfer and default ShapeFix. Different readers and tolerances prevent a common-tolerance or source-fidelity claim. Checkpoint 1's legacy result, 374 valid exports, 12.47%, is historical. It must not be relabelled as a new strict occurrence-preserving measurement.

Certified bodies, checked numerical bodies and sampled-native tolerant bodies are reported separately. A full file's group is its least-checked body: any sampled-native body places it in that group; otherwise any numerical body places it in numerical; only an entirely certified file is certified-only. These are disjoint file groups. Body counts include partial files and are not file counts. Retained source surfaces are not admitted solids.

The phase-2 import evidence needs particular care. [The per-body repair PR](https://github.com/BohnBawerick/catallax/pull/697) rejected its earlier native-first ordering after it lost former checked admissions. Its corrected targeted observations and 114-file selected-document checks are evidence at their recorded producers, not a full census of the final merged main. The new checkpoint run, not a sum of those cohorts, supplies the current headline.

### Legacy dry run, historical protocol

The unchanged legacy command also completed end to end on main `403081f8`, 18 Nov 2026, in 1,003.88 seconds, 16.7 minutes. It records 647 full admissions, five partials, 2,346 refusals and two import timeouts. Of the full admissions, 601 export valid after builder processing, 40 export invalid, five refuse at the export work limit and one times out. The 652 admitted files are labelled 190 certified and 462 numerical approximate by that old entry.

| Historical-protocol measure | Files | Share of 3,000 |
| --- | ---: | ---: |
| Legacy import/export dry run, main `403081f8` | 601 | 20.03% |
| Legacy import/export, Checkpoint 1 cut | 374 | 12.47% |
| OpenCASCADE healed-valid reference | 2,701 | 90.03% |

No previously valid file is measured lost against the old 369-file reference. Import times out on 00000187 and 00008201; export times out on 00004777. None reached the observer memory cap. This diagnostic run is not the new tolerant-tier headline. The default-tier collector is running separately after the approved wrapper change; no current tier total is inferred from these legacy rows.

## What phase 2 delivered against ROADMAP.md

| Roadmap item | Landed delivery | Acceptance still open or limited |
| --- | --- | --- |
| Measures first | [Fixed material-cutting census holes and added P8, P9 and curved-housing targets](https://github.com/BohnBawerick/catallax/pull/677). | P6 and P7 remain outside this comparison; no longer edit chain was defined. |
| 1.3 carry-over, loft top cuts | [All 14 recorded spline Boolean requests reach approximate STEP](https://github.com/BohnBawerick/catallax/pull/674). The three original top cuts use 26.78M, 27.29M and 25.95M work at the 100M application allowance. | The three are above the 20M target. Their recorded raw-reader and later-cut proof is not a fresh checkpoint rerun. |
| 2.1 spline fillets | [Numerical multi-patch fillet construction and raw STEP proof](https://github.com/BohnBawerick/catallax/pull/686), then [ordinary application rolling-ball fallback](https://github.com/BohnBawerick/catallax/pull/698). Exact builders stay first; numerical results stay approximate at 1e-7 mm. Native strips round plane/cylinder rims up to the neighbouring support and refuse beyond it with the limiting face and radius. | Partial. The earlier multi-patch loft builds, but its radius edit exhausts 100M. Imported multi-patch acceptance, edited-name preservation, the complete curved-housing chain and all later housing cuts were not established by that delivery. A working rim is not proof of the entire milestone. |
| 2.2 NURBS toolkit | [Coons, point-grid, Gordon, guided skin and N-sided fill; public knot insertion, degree raise, split and join](https://github.com/BohnBawerick/catallax/pull/676). All six emitted surface artifacts passed the recorded healing-off reader. | Delivered as finite surface definitions and library APIs. N-sided fills are C0, not G1; no solid or global embedding claim follows from surface export. |
| 2.3 surface bodies | [Retained-import sewing](https://github.com/BohnBawerick/catallax/pull/678), then [surface trim, extend, signed thickening, free wires and constructed-face sewing](https://github.com/BohnBawerick/catallax/pull/681). The recorded six-face sewn box is independent of STEP input. | Library acceptance is demonstrated. These are not new recipe operators or SolidDraft surface-edit commands. Signed thickening/export costs about 29.15M, above the target. |
| 2.4 general sweeps | [Taper, guide/twist laws, numerical helix sweeps and curved-face termination](https://github.com/BohnBawerick/catallax/pull/688), with [further P8 carrier repairs](https://github.com/BohnBawerick/catallax/pull/694). Build/edit/save/reopen and raw STEP evidence exists for the delivered families. | Open. P8 does not complete its threaded union, so spring construction and whole-part STEP are not accepted. |
| 2.5 direct edits | [Delete a face, fillet or hole; move, offset or replace a face](https://github.com/BohnBawerick/catallax/pull/687), retaining the original admitted source and publishing reconstructed candidates atomically. | Kernel/library API, not a recipe or SolidDraft command. Six real-import artifact pairs were qualified before review corrections; current focused controls are separate. Do not claim that historical six-pair run was repeated on final main. |
| Import resolve lane | [Multiple holes and complete circular rims](https://github.com/BohnBawerick/catallax/pull/680), [numerical export checks](https://github.com/BohnBawerick/catallax/pull/682), [cylindrical period charts](https://github.com/BohnBawerick/catallax/pull/684), [default native tolerant solids](https://github.com/BohnBawerick/catallax/pull/691), [torus, offset and extrusion families](https://github.com/BohnBawerick/catallax/pull/692), then [checked-first per-body repair and selected-document export](https://github.com/BohnBawerick/catallax/pull/697). | Checked, numerical, sampled-native and retained surfaces remain different authorities. Cap-certified incidence records conservative upper bounds, not measured residuals. Mixed geometry units across prototypes still refuse. |
| Display and edit fixes | [Reuse unchanged face triangulations](https://github.com/BohnBawerick/catallax/pull/689), [preserve authored legacy sketch endpoints](https://github.com/BohnBawerick/catallax/pull/690) and [correct analytic triangulation metrics](https://github.com/BohnBawerick/catallax/pull/696). | No general click-latency qualification follows. Mesh bounds and refused display results remain explicit. |
| App alongside | [Surface import display and import health](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/329), [immediate Cancel and import-message corrections](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/335), [numerical bodies retained as approximate](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/336), and the display repin. | The final [fillet acceptance and repin](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/339) is catch-up work at the start of reconciliation. Its merge status is checked at the cut. The two kernel commands do not prove app build/edit/Undo/save/reopen behavior. |

## Known gaps and what still refuses

The dry scoreboard pins the current failures more precisely:

- P0 general cone and torus still need `intersect.traced-pair-required` branch tracing.
- P3 has one wrong construction membership value; its STEP validity is not enough to pass the fixed rule.
- P5 finishes geometry and export but fails the names check.
- P8 stops at `intersect.parametric-surface`, before the spring, at 55,626,285 work in both variants.
- P9 stops at `loft.admissibility-unproven`, before its inner subtraction.
- The curved housing stops at `intersect.unsupported-pair` between its rational atlas and Bezier extrusion, before the fillet.
- ABC 00000173 still stops at `boolean.no-rational-point` on the drill.

- P8 spring and threaded boss is accepted as unfinished. The latest hand-off reports the threaded union stopping at 99,879,180 of 100M work. This is consumed work before refusal, not a completed union. The earlier merged carrier-repair receipt stopped at `intersect.parametric-surface` after 55,626,285 work. The fresh scoreboard must retain whichever result its pinned producer actually gives; these observations are not interchangeable.
- Assemblies whose prototypes use different length units still refuse on purpose. Distinct placement contexts for a common prototype are supported; conflicting prototype geometry units are not silently rescaled.
- Numerical fillets do not establish whole-body approximation bounds, global branch completeness or exact rolling-ball authority. `keep_neighbours` and smooth-corner constructions can still refuse. A radius limit is a geometric answer, not a timeout.
- The multi-patch radius edit, imported fillet acceptance and full housing chain remain separate from the working plane/cylinder rim acceptance. The [partial 2.1 receipt](https://github.com/BohnBawerick/catallax/pull/686) records those remaining gates.
- In import, source opening, admitted solid, valid isolated prototype and complete occurrence-preserving export are four different results. Resource exhaustion, timeout, memory termination and unobserved siblings must not be described as geometric invalidity.

## Later: optimisation notes from merged PRs

These extend [Checkpoint 1's list](checkpoint-1.md#later-optimisation-notes-from-merged-prs). They are deferred work, not changes made by this reconciliation.

- [NURBS toolkit](https://github.com/BohnBawerick/catallax/pull/676): reduce rational-product degree and redundant Bezier seams without changing the definitions.
- [Import sewing](https://github.com/BohnBawerick/catallax/pull/678): reuse source recognition across whole-file and body attempts without changing the work policy.
- [Surface thickening](https://github.com/BohnBawerick/catallax/pull/681): profile the positive and negative thickening/export chains, about 29.15M work each, against the 20M target.
- [Multi-patch fillets](https://github.com/BohnBawerick/catallax/pull/686): exact offset proof across perpendicular-plane atlas seams and optional numerical-body handover encoding each exhaust 100M. The successful loft build and separate export use 96,732,996 and 79,621,750 work. Separate requests do not establish a combined 100M app request.
- [Native tolerant import](https://github.com/BohnBawerick/catallax/pull/691): investigate native-fit work. Its historical targeted replay does not replace a uniform full census.
- [Per-body repairs](https://github.com/BohnBawerick/catallax/pull/697): reduce trace-rebuild and curve-fit work within unchanged caller controls.
- [Certified display triangulation](https://github.com/BohnBawerick/catallax/pull/696): fixed trim samples can still leave thin facets. Keep reporting aspect bounds rather than hiding or relaxing them.
- [Scoreboard cutter](https://github.com/BohnBawerick/catallax/pull/677): up to 98 OCCT common operations per part, with 44 for 00000632. Reading bore axes first could reduce this unmeasured setup cost. Moving the three phase-2 targets into `parts/` also waits until its CI job is ready for them.
- [Browser app](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/329): edited-body meshing took two thirds of worker time in the recorded edit; CPU picking over the shown mesh could avoid GPU readback. Later Cancel work has since landed, so those old readback numbers are not a current Cancel verdict. Retained unhealed faces still need a usable surface display representation.
- [Browser lazy reads](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/323): bracket opening, worker build time and first exact-body reads remain recorded costs. Its visibility-sensitive pick-readback proposal was not a valid fix merely because one edit became faster.

The P8 thread sweep and refused threaded union remain in the [phase-5 work cases](../phase5-work-cases.md). The 20M target and 100M application allowance have not changed.

## Release order

After the cut: Catallax first, then SolidDraft pinned to that latest Catallax, then the kernel site update. Release execution is outside this reconciliation. Do not publish draft numbers as final release measurements.

## Draft order of phase 3

The owner already chose drawings from PMI to follow this release. The rest is a proposal, not a new execution order.

1. Drawings from PMI. Put PMI on the part, carry AP242 in and out, recognize features, lay out dimensions, run the interview and produce STEP plus a drawing with one action. This is the chosen user outcome, rather than an isolated hidden-line algorithm.
2. Assembly and exchange completion, roadmap 3.2. Preserve source occurrences, names, colours and annotations through round trips before claiming whole-product exchange. Keep mixed prototype-unit refusal explicit until supported.
3. Clash, clearance and mass, roadmap 3.3. Every assembly pair needs clash, clear by a stated gap, or unresolved; bounded answers must remain distinguishable from proof.
4. Remaining SMLib-only tools, roadmap 3.4. Bend, twist, curvature-continuous blends and mesh operations follow the selected drawing and product-exchange work because they unlock fewer of those paths.

Unfinished phase-2 gates stay visible alongside this proposal. General optimisation remains phase 5; the accepted P8 and mixed-unit gaps do not become completed milestones at the checkpoint.
