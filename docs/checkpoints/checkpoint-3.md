# Checkpoint 3

This checkpoint measures [Catallax v0.5.0](https://github.com/BohnBawerick/catallax/releases/tag/v0.5.0), commit `a59f3aa9a07a217b0c963d68c97865a9804619ee`. The scoreboard and [Checkpoint 4 STEP census](checkpoint-4.md#step-census) share one cut. [SolidDraft v0.8.0](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/releases/tag/v0.8.0) pins that kernel in its dependency and lockfile. The phase delivery tables use the release's landed evidence.

## Rerun the cut

Run each command once from the clean v0.5.0 checkout, on the work PC or an equivalent Linux build host. Give the run its own build output directory and keep at least 20 GB free. Pass runtime paths for the comparison installation, frozen STEP sources and selection. The source files and raw outputs stay outside the checkout.

```sh
python3 -B tools/checkpoint/checkpoint.py scoreboard --tools "$TOOLS" --sources "$SOURCES" --selection "$SELECTION" --out "$RUN/scoreboard"
python3 -B tools/checkpoint/checkpoint.py census --tools "$TOOLS" --sources "$SOURCES" --selection "$SELECTION" --out "$RUN/census"
```

The wrapper and geometry came from the tag. The scoreboard ran on the shared Linux work PC with one pinned core, one warmup and one measured sample per engine, 100M work per Catallax phase and 120 seconds per engine. It used CadQuery 2.6.1 and smlibpy 0.4.0 as reference builders and OCP 8.0.1.1 for STEP validation. The Mac ran the checkpoint reporting regression and reproduced all 26 scoreboard rows from the raw observations.

## Edit-chain scoreboard

The [Checkpoint 2 definition](checkpoint-2.md#edit-chain-scoreboard) applies, with the P3 release checker changes below. Both the original and edit must build, and both STEP exports must open as valid closed solids with the producer's face count and no wrong checked value. Each census-part drill and radius edit must reduce the independently measured volume by more than one millionth of the source volume. Catallax name survival is reported separately. The curved housing's healing-off read is diagnostic. The command compares P0 through P5, P8, P9, the curved housing and ten census-part edit chains. P6 and P7 have separate field printing checks.

The scoreboard completed in 907.42 seconds (15.1 minutes), including its build. The host was an Intel Core i7-14700. The scoreboard binary SHA-256 is `d9ff8cd0294e578095c4ce4423f96c4368a4bc7d6b897f4bb23d66ba405e37ba`.

| Target | Catallax | OpenCASCADE | SMLib | Catallax work, build / edit | Names |
| --- | --- | --- | --- | --- | --- |
| P0 block-sphere | done | done | done | 633,934 / 397,086 | passed |
| P0 crossing-pipes | done | done | edit build command-failed | 3,169,337 / 1,997,714 | passed |
| P0 equal-cylinder-sphere | done | done | done | 2,813 / 76,393 | passed |
| P0 general-cone | construction build refused | done | done | 2,341 / 2,341 | not-run |
| P0 general-torus | construction build refused | done | done | 4 / 4 | not-run |
| P0 line-tangent | done | done | construction build command-failed | 35,038 / 39,208 | passed |
| P0 near-miss | done | done | construction build command-failed | 3,625 / 3,590 | passed |
| P0 seam | done | done | done | 296,359 / 280,564 | passed |
| P1-bracket | done | done | done | 13,618,847 / 13,479,812 | passed |
| P2-collar | done | done | construction 2 wrong checked values | 1,345 / 1,382 | passed |
| P3-filleted-boss | done | construction STEP raw read invalid | construction STEP raw read invalid | 20,773 / 21,591 | passed |
| P4-shelled-tray | done | done | done | 917,132 / 936,262 | passed |
| P5-flanged-housing | done | construction build command-failed | construction build command-failed | 1,063,025 / 1,062,466 | failed |
| P8-spring-threaded-boss | construction build refused | done | done | 55,626,285 / 55,626,285 | not-run |
| P9-lofted-duct | construction build refused | done | done | 4,793,088 / 4,793,088 | not-run |
| curved-housing | construction build refused | done | done | 1,511,266 / 1,511,274 | not-run |
| ABC 00000131 | done | done | drill raised the volume | 469,022 / 479,854 | passed |
| ABC 00000173 | construction build refused | done | done | 1,989,712 / 1,989,922 | not-run |
| ABC 00000190 | done | done | done | 31,083 / 30,963 | passed |
| ABC 00000405 | done | done | done | 415,413 / 410,612 | passed |
| ABC 00000511 | done | done | done | 30,403 / 30,283 | passed |
| ABC 00000605 | done | done | done | 27,077 / 27,250 | passed |
| ABC 00000632 | done | done | done | 697,371 / 1,128,070 | passed |
| ABC 00000633 | done | done | done | 112,187 / 113,026 | passed |
| ABC 00000695 | done | done | done | 86,583 / 83,250 | passed |
| ABC 00000821 | done | done | drill raised the volume | 218,037 / 218,237 | passed |

| Kernel | P0 to P5, of 13 | Phase-2 targets, of 3 | Census chains, of 10 |
| --- | --- | --- | --- |
| Catallax | 11 | 0 | 9 |
| OpenCASCADE, CadQuery 2.6.1 | 11 | 3 | 10 |
| SMLib, smlibpy 0.4.0 | 7 | 3 | 8 |

Healing-off STEP at 1e-7 mm for both variants.

- P3-filleted-boss. Catallax valid. OpenCASCADE, CadQuery 2.6.1 invalid. SMLib, smlibpy 0.4.0 invalid.
- curved-housing. OpenCASCADE, CadQuery 2.6.1 invalid. SMLib, smlibpy 0.4.0 invalid.

### Catallax gaps

- P0 general-cone. Construction build refused at `intersect.traced-pair-required`.
- P0 general-torus. Construction build refused at `intersect.traced-pair-required`.
- P8-spring-threaded-boss. Construction build refused at `intersect.parametric-surface`.
- P9-lofted-duct. Construction build refused at `loft.admissibility-unproven`.
- curved-housing. Construction build refused at `intersect.unsupported-pair`.
- ABC 00000173. Construction build refused at `boolean.no-rational-point`.

### What changed since v0.4.0

- P3 moved from a Catallax membership failure to a pass after the checker correction. The release checker addresses false inside classifications and adds four junction probes and a required healing-off STEP read for all three engines. Catallax has zero wrong values in both variants, with 192 classified samples and 28 boundary samples per variant. Both reference engines pass the normal reader. Their original and edited STEP files fail the new healing-off gate. Their P0 to P5 totals each drop by one under the tagged checker. [P3 checker change](https://github.com/BohnBawerick/catallax/pull/755)
- Catallax P5 passes the build and export gate and records 16 unnamed entities in the separate name check. P0's general cone and torus still refuse. On the census drill chains, ABC 00000173 still refuses before its cut, and SMLib's drills on ABC 00000131 and 00000821 raise the measured volume.
- The three phase-2 chains still refuse in Catallax. P8 stops at the threaded union before the spring, P9 at the outer loft before the inner cut, and the curved housing at the spline cut before its fillet.

## Coverage recount

The [line-by-line recount](checkpoint-3-coverage.md) awards 31 of 61 points (50.8%) under the adopted acceptance bars. All 61 lines have working evidence. Spline editing meets its full bar and 60 lines have partial credit. The older 55% and 61% figures use their historical scoring. The census and edit-chain rates have separate denominators.

## Phase 3 delivery

| Roadmap item | Delivered | Open at the cut |
|---|---|---|
| 3.0A tolerance-first core | The edge-gap checker admits a measured local tolerance or refuses. The export records its tolerance and the check-once path reopens it. [Gap checker](../milestones/3.0A-edge-gap-checker.md) [Check-once evidence](../milestones/3.0A-check-once.md) | General edge and import coverage remains open. The STEP census measures the net result under the existing limits. |
| 3.0B one model and predictor | Exact and tolerant entities coexist. A supported through hole after a fitted fillet builds and writes checked STEP. The 117-row predictor corpus exercised multiple worker counts on two machines. [One model](../milestones/3.0B-one-model.md) [Predictor](../milestones/3.0B-predictor.md) | Other later-feature families can still refuse `features.tolerant-input`. The 2M predictor tier remains unmeasured. The owner accepted the measured wall-time band at its limit. |
| Drawings and PMI | P1 bracket PMI travels in AP242 and drives a manufacturing package. The native app bore edit keeps its tolerance. Bare-import feature recognition and the interview have app evidence. [P1 specification](checkpoint-2-drawings-pmi.md) [SolidDraft acceptance](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/346) [Hole interview](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/355) | OpenCASCADE reads 13 of the 14 bracket items. The owner retained the NIST-accepted representation. The browser bore-edit proof remains open. |
| 3.1 drawing views | SolidDraft drawings use kernel section, silhouette and hidden-line curves. Fitted pieces are labelled approximate and unsupported exports report a refusal. [Drawing curves](../interface-changes/2026-11-21-sections-and-drawing-curves.md) [SolidDraft drawing acceptance](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/357) | Some fitted and unsupported views still lack exact curves. Mesh preview remains available. |
| 3.2 assemblies and exchange | The rooted STEP assembly round-trips shared occurrences, rigid placements, names, per-part colours and supported PMI. Whole spheres and ring tori export to IGES. [Exchange evidence](../milestones/3.2-assemblies-and-exchange.md) | Wider assembly encodings, per-occurrence and per-face styles, general IGES and complete annotation acceptance remain open. |
| 3.3 clash, clearance and mass | The 78-pair test assembly obtains 38 clear, two clash and 38 named unresolved answers. SolidDraft uses the three-result contract. [Assembly evidence](../milestones/3.3-clash-clearance.md#results-3) [SolidDraft checks](https://github.com/ZIGRAWBOHNBAWERICK/soliddraft/pull/357) | The later coverage bar requires resolving pairs outside tolerance. Field-only envelopes, fitted P5 pairs, exact contact and Boolean work exhaustion account for the unresolved results in the recorded exercise. |
| 3.4 SMLib-only tools | Bend, twist, G2 fillets, mesh Booleans and decimation have first builds. [Bend and twist](../RECIPE.md#420-bend-and-twist) [G2 fillets](../RECIPE.md#453-curvature-continuous-fillets) [Mesh Booleans and decimation](../interface-changes/2026-11-21-mesh-boolean-and-decimation.md) | General input remains open for all three coverage lines. |
| 3.5 one display mesher | Shared B-rep charts display analytic and spline faces, including the M11 loft fillets and the three-turn coil. [Display evidence](#one-display-mesher) | Corpus-wide correctness and the downstream SolidDraft M11 scenario remain open. The work and size gates move to [phase 5.1](../../ROADMAP.md#display-mesher-work). |

### One display mesher

| Acceptance line | Release evidence |
|---|---|
| Frozen corpus at 0.01 and 0.005 mm | The three M11 loft controls display at both chords. The complete frozen corpus has no passing receipt. [M11 display test](../../crates/catallax-features/tests/m11_spline_fillet.rs) |
| Shared edge polylines and closed, chord, angle and no-fold checks | The shared-edge path and focused closed and bounded controls landed. The complete corpus check remains open. [Shared chart delivery](https://github.com/BohnBawerick/catallax/pull/748) [Chart checks](https://github.com/BohnBawerick/catallax/pull/754) |
| Dense independent checks against each face bound | Landed spot checks use a 3 by 3 grid on the affected faces. The dense per-face checks remain open. [Chart checks](https://github.com/BohnBawerick/catallax/pull/754) |
| Point prediction equals production | The release has no passing prediction check for the corpus. [Release exits](../milestones/M5.md#050-phase-4-alpha-release) |
| Mesh bytes and work agree at 1, 2, 4 and 24 workers on two machines | Focused worker-count and differential checks landed. A complete corpus pass on two machines remains open. [Worker-count checks](https://github.com/BohnBawerick/catallax/pull/752) [Differential checks](https://github.com/BohnBawerick/catallax/pull/754) |
| SolidDraft displays the M11 loft fillet with one later feature | The release keeps this scenario open. The kernel pin and a successful display of the fillet by itself leave the downstream feature unverified. [Release exits](../milestones/M5.md#050-phase-4-alpha-release) |
| Whole-body forward compare on M11 and P5 inside 100M | Selected rows have focused receipts. Complete coverage at the cut remains open. [Shared mesher](https://github.com/BohnBawerick/catallax/pull/748) [Independent spline compare](../M42A-FIELD-FILLETS.md#forward-display-comparison) |
| Three-turn coil displays with its periodic seam | The retained kernel control builds checked display and matching previews, including the closed mesh check. Keep the SolidDraft scenario and its periodic seam in the acceptance. [Coil control](../../crates/catallax-app/tests/coil_display.rs) |
| Field charts display all built 4.2a fillets and patch shells | The named M11 loft patch shell displays. Other shell/fillet routes and complete two-sided dense checks remain open. [Field chart scope](../../crates/catallax-field/src/patch_shell/display.rs) |
| M11 display under 1M at 0.005 mm | The three loft tests assert this limit. It moves to phase 5.1 with the same limit. [M11 display test](../../crates/catallax-features/tests/m11_spline_fillet.rs) |
| Each analytic body at or below its baseline work | The differential run pins exceptions. A complete passing analytic-work table remains open in phase 5.1. [Differential checks](https://github.com/BohnBawerick/catallax/pull/754) |
| Three-turn coil inside 100M and a few turns under 20M | The focused test asserts the 100M allowance and prints work against 20M. Both figures move to phase 5.1. [Coil control](../../crates/catallax-app/tests/coil_display.rs) |
| At least 2,000 fewer `catallax-mesh` lines | The crate gained 2,990 lines (4,873 insertions and 1,883 deletions) between the parent of the shared-mesher delivery and the tag. The size goal moves to phase 5.1. [Shared mesher](https://github.com/BohnBawerick/catallax/pull/748) |

The release retains the broader spline fillet and sweep gaps. A successful numerical loft fillet does not establish radius-edit survival or the complete curved housing. P8's refused threaded union and work above 20M remain in the [phase-5 work cases](../phase5-work-cases.md).

Continue with [Checkpoint 4](checkpoint-4.md) for the census, field delivery and speed cases.
