# Patch-first fillet ownership, round 2

## Preregistration (before code or measurements)

**Question.** Can contact rails, rational circular sections and a small fitted B-spline cover replace round 1's finite-sphere unions while meeting the merged error policy and the 1,000 ms cached-input feature/display budget? Baseline: pinned `main` **321e416e7ed740beb349b6e87a18c9f1a7141705** (verified against live remote). Read PLAN §§8.2/8.4, ARCHITECTURE's merged policy, RECIPE §4.5 and [round 1](fillet-ownership.md). No AGENTS.md exists in this checkout or its ancestors; `.agents/skills/verify-catallax` is the applicable verification skill. No matching patch-first branch/PR exists; round 1 remains historical.

**Candidates.** A: four cubic longitudinal spans per octant, quadratic rational sections (32 tensor patches for G1). B: two degree-seven longitudinal spans per octant, the same sections (16 patches). Interpolate homogeneous section controls at Chebyshev-Lobatto stations; use each fitted surface's boundary columns as its boundary rails, never independently fit the rails. These are construction candidates, not certified B-rep answers. Keep elementary P3 carriers exact. No sphere unions, production changes, measuring-tool edits or added kernel carrier family.

**Cases, in order.** P3 unchanged `parts/P3-filleted-boss.json`, then junction radius 2→5/2 and block-top radius 1→3/2 with IDs/selections retained. Next G1: union of radius-8, height-32 z/x cylinders centred at the origin, selected wall/wall edge, fillet radius 1 then 3/2. Its offset-centre branches satisfy x=±z and x²+y²=(8+r)². For each signed octant let t∈[0,1], x=(8+r)(1−t²)/(1+t²), y=(8+r)2t/(1+t²), |z|=|x|. Contacts are C−r(x,y,0)/(8+r) and C−r(0,y,z)/(8+r). Sections use the shorter exact circle between these contacts, including its algebraic middle weight. Patches must cover all eight octants; collapsed sections at |y|=8+r require explicit end/topology treatment. If G1 lacks certified, valid trimmed solid coverage, stop expansion: G2 oblique cone/plane and G3 offset sphere/cylinder from round 1 are recorded unmeasured.

**Frozen measurement.** One measured round only; no repeat is authorized. Unchanged release verification CLI and driver, 100,000,000 request work, face options width 1/10,000,000 mm, 4,096 cells, 64 path segments; default display/export settings. Use nice 10, one doctor-selected core, at most two build jobs/test threads; retain actual priority/affinity/load, commands, hashes, exit codes and artifacts outside the repository. P3 baseline kernel-loop recipe/field/compare requests precede candidate work; use unchanged face checks for the five analytic blend faces. The prototype harness constructs each fixed G1 layout once per radius, with no fit tuning, measuring-tool changes or tolerance/budget increases. Fixed diagnostics: 33×17 grid per patch, 17 edge points; these are samples, never certificates. Record degree/spans/control counts, sampled definition-correspondence error in nm, boundary reuse, minimum sample Jacobian, actual OCCT face/mesh/STEP outcomes and topology/coverage limitations. Attempt unchanged measurement admission; unsupported raw spline input is an explicit prerequisite gap, never a zero deviation.

**Separate phase clocks, preregistered.** Contact tracing (centres/rails), section construction, homogeneous fitting, face trimming/topology, unchanged whole-face certification/admission, and display meshing; record milliseconds independently, including skipped phases and failures. Record input construction, whole request, export and independent import separately. Partial/open patch display is not first display of a valid filleted solid. No speedup claim without the skill's warmup and three alternating matched pairs; these single construction observations cannot decide click performance. Application observations are unmeasured unless the actual application loop is available.

**Decision.** Promote only if complete trimmed solid, two-sided whole-face certificates ≤100 nm (hard cap 1 µm and one tenth of the nearest feature), shared incidence/orientation, edit identity, unchanged earlier answers and independent valid STEP all pass, with canonical matched-pair click evidence. Any larger sample error can falsify the target but cannot certify it. Missing raw-spline measurement alone is a prerequisite gap, not a geometric failure. Park a layout that demonstrably misses a required gate or remains an incomplete construction at this bound; distinguish that from impossibility. `repeat once` is allowed as a verdict only for a wrong measurement and authorizes no further work here. Prototype implementation stays on `prototype/fillet-patch-first`, never merged; deliver exactly one docs-only draft PR from a separate branch based on the pinned main. Do not change PLAN or production files.

**Pre-measurement source addendum.** The requested updated §8.2 is in the still-open roadmap [PR #58](https://github.com/BohnBawerick/catallax/pull/58), read at head `5e5adade34ce141e580a4b68513bf0a075aa22e8` before measurement. It explicitly adds strip/trimmed-neighbour shared-curve incidence and collapsed-boundary treatment. These remain gates above; the baseline remains pinned main and PLAN is untouched. The report distinguishes reuse of a raw surface boundary from sharing that edge with a trimmed input wall. No budget, tolerance, layout or round count changes.

## Results

### Verdict: park

**One round completed; neither layout is promoted. No repeat was run or requested.** The compact degree-seven cover gives small sampled correspondence errors, but the constructed artifacts are open compounds with **zero solids**, no trimmed input walls, no shared neighbour incidence and no whole-face certificate. That incomplete construction independently prevents promotion. The missing general-spline measuring route is a prerequisite gap, not a measured geometric failure. These findings do not prove that patch-first fillets are impossible or too slow.

### P3 exact control, before the G1 probes

The unchanged kernel-loop CLI built both P3 recipes as exact 13-face solids (five cylinders, seven planes, one torus). The existing face-check wrapper measured all five blend faces at each radius set. Each returned the whole-face interval **[0, 0] nm**, `converged: true`, by the measuring stick's analytic identity route; this is a certificate, not a grid maximum. The separate two-radius rebuild retained all five measured blend origins, with five resolving once, zero missing and zero ambiguous. The wrapper's additional scalar edit also resolved all five references. This establishes only the elementary control.

| P3 radius set, junction / top | Plain recipe request, ms | First B-rep mesh from face-check request, ms | Whole face-check request, ms | Five-face measurement spans, ms | B-rep triangles / vertices |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 / 1 | 36.624 | 213.055 | 231.855 | 0.091 | 10,584 / 5,294 |
| 5/2 / 3/2 | 19.143 | 270.481 | 291.922 | 0.109 | 18,912 / 9,458 |

These are single observations, not paired latency claims. First-mesh clocks include prefix construction and precede certification/edit/export. Internal contact/section/fitting subphases are not separately exposed for the unchanged elementary builder. Default B-rep chord tolerance is 0.01 mm, distinct from the zero construction-deviation certificates.

Independent OpenCASCADE 7.8.1.1 imports of both plain-recipe and both face-check STEP artifacts passed without repair, each valid with 13 faces. Plain and face-check STEP bytes match for each radius set. Kernel-volume agreement was not requested or measured. The independent field-display and comparison requests **all exhausted the frozen 100,000,000-unit allowance**: original/edited field requests took 9,156.868 / 9,956.268 ms; comparisons took 16,168.273 / 15,256.625 ms, all until refusal. No OBJ or STEP remained for those refusals. P3 therefore does not establish full target-part acceptance or two-kernel agreement.

### Contact sections and fitted G1 patches

The G1 prefix built an eight-face input solid with four elliptic and four circular edges (207.471 ms single request). The selected wall/wall run reaches the elementary builder's `fillet.unsupported-topology` refusal, “elementary fillets need a line or circle”; the input Boolean itself succeeded. The scratch constructor uses the preregistered closed-form offset-contact branches, not ball unions or an independently fitted boundary rail.

For centre C and unit vectors nA, nB from C toward the feet, set w=√((1+nA·nB)/2). The exact rational quadratic section has homogeneous controls `(C+r nA,1)`, `(w C+r(nA+nB)/(2w),w)`, `(C+r nB,1)`. The implemented station coordinates/weights and fits use binary64 arithmetic; naming the analytic formula does not make these numerical patches exact circular sections or exact contact traces. Each longitudinal fit interpolates these three homogeneous controls, with boundary columns reused by `VIso(0)` and `VIso(1)`. This supplies one surface/boundary construction, but the input-wall trims and their shared topological edges were not constructed.

| Candidate | Radius, mm | Patches | Homogeneous controls | Sampled correspondence error, nm | Sampled contact-rail correspondence error, nm | Whole-face certificate |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| A, degree (3,2) | 1 | 32 | 384 | 1,383.617 | 1,251.533 | unavailable |
| B, degree (7,2) | 1 | 16 | 384 | 30.051 | 3.033 | unavailable |
| A, degree (3,2) | 3/2 | 32 | 384 | 1,628.743 | 1,283.779 | unavailable |
| B, degree (7,2) | 3/2 | 16 | 384 | 36.628 | 3.080 | unavailable |

These diagnostics compare corresponding parameters on a fixed 33×17 grid per span in the positive octant; sign reflections preserve Euclidean error. They neither bound the unsampled face nor give a certified lower bound on closest-point deviation. In particular, A's correspondence errors above 1 µm are not reported as a proven violation of the policy's distance cap, and B's sub-100 nm samples are not a target pass. No local tolerance was issued. All sampled surface/own-boundary discrepancies were zero at the 17 boundary stations; this does not prove incidence with the absent trimmed neighbours. Minimum sampled finite-difference Jacobian was zero in every layout at the collapsed end. Regularity, branch orientation, collapsed-end topology and reverse coverage remain unproved.

### Preregistered phase observations

| Candidate / radius | Contact, ms | Sections, ms | Fitting, ms | Rectangle-face construction, ms | Whole-face certification | Open-patch display, ms |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| A / 1 | 0.492 | 0.152 | 0.882 | 13.215 | not run | 61.317 |
| B / 1 | 0.097 | 0.064 | 0.211 | 6.387 | not run | 85.000 |
| A / 3/2 | 0.119 | 0.068 | 0.240 | 9.861 | not run | 61.878 |
| B / 3/2 | 0.085 | 0.044 | 0.231 | 6.655 | not run | 88.266 |

The preregistered trimming clock reached only bounded parameter-rectangle faces. It did not reach input-wall trimming, sewing or a publishable solid. Contact is station evaluation of the analytic branches; reflection reuses that construction. Certification is null/unmeasured, not a near-zero successful phase. Diagnostics, imports and exports are excluded from these phase clocks. No end-to-end click result or speedup over round 1 is claimed; paired warmups/repeats were not run because the construction gate already failed.

All 32/16 individual faces were OCCT-valid. Each open compound exported and reimported without repair as a valid 32/16-face compound with **zero solids**. A compound validity check is not a valid-solid or fillet-definition check.

| Candidate / radius | Open-patch triangles / vertices | STEP export, ms | Independent read-back, ms |
| --- | ---: | ---: | ---: |
| A / 1 | 6,480 / 3,656 | 19.937 | 34.645 |
| B / 1 | 7,248 / 3,928 | 7.084 | 17.192 |
| A / 3/2 | 7,296 / 4,096 | 10.966 | 29.027 |
| B / 3/2 | 7,728 / 4,184 | 6.216 | 17.071 |

OCCT was only the prototype patch constructor, mesher and artifact reader. Its chord setting and validity checks are not the Catallax whole-face measuring stick. No candidate was substituted into the production fillet builder or reported through the CLI as an answered solid.

### Unchanged kernel-loop gates and proof gap

Both G1 radii refused ordinary B-rep requests with `fillet.unsupported-topology` and field display with `field.fillet-contact`. The actual existing certified-recipe entry point also refused both G1 recipes during elementary fillet construction, before any spline certification: 178.670 / 211.071 ms and 3,432,721 work units each. These were recipe admission probes, not measurements of the raw OCCT patches. Comparison requests returned **unresolved**, with zero witnesses; their reported sampled bound of zero is not agreement. Refused requests left no STEP/OBJ artifacts.

At pinned main, [`SplineSupport`](../../crates/catallax-shapes/src/spline_surface.rs) admits only checked `EllipticOffset` surfaces. Raw B-spline controls cannot be installed as a general fillet surface. The unchanged [`measure_face`](../../crates/catallax-field/src/measurement.rs) supplies analytic chart covers and speed bounds; it has no general fillet-spline chart route. The prototype did not extend either interface, change the measuring tool or relabel an elliptic-offset proof.

A later, separately authorized implementation would need bounds on the *homogeneous* section interpolation over every knot span, a strictly positive denominator bound, positional correspondence in both directions and a field measurement against the retained selected centre-path branch. It must also prove the exposed envelope, shared trimmed-neighbour incidence, end-collapse treatment, material orientation, nearest-feature caps and complete body coverage. Those are proof obligations, not results of this round. An interval bound on station fitting alone would not satisfy them.

G1 blend/corner edit identity, full-model STEP and complete first display are unmeasured because no candidate solid exists. G2 and G3 were not expanded after the G1 gate. The SolidDraft application loop, earlier-answer corpus/failure-pack replay, and one/two/four-worker promotion proofs were not run. No production implementation changed; nothing here establishes general-fillet ownership. Stop at this bound and park the scratch construction. A future trial needs new authorization and the missing construction/proof support; “repeat once” is not the verdict.

### Evidence and provenance

- Pinned baseline main: `321e416e7ed740beb349b6e87a18c9f1a7141705`. Built in a separate detached checkout, release, nice 10, two jobs; preserved CLI SHA-256 `9db4617f9858808406d67fac509ac9d11d54d0bf7aae7e064a16c21887c924bd`, with doctor `about.request = catallax-request/1` and embedded commit equal to the pin.
- Original preregistration commit: `ebab36d96f225fe3961c520099e4048c5b9731e3`, committed at **2026-10-24 05:20:17 UTC**; original file SHA-256 `61e306e4433652aa285c3644f88d47863c1ce11727b830b046a877d060a3b4de`. The pre-measurement roadmap-source addendum is `04ba013eb4220a342fa24766eba63e3388a7f1ab`; frozen protocol SHA-256 `fc371096f40860feb61e9990ddf5e823c6dd70cddfb3977c0e539b3ebac68154`.
- Measured local prototype head: `34635575d0e2cc524112cbc40210aadc8b293cd0`, branch `prototype/fillet-patch-first`. Harness `tools/prototypes/fillet_patch_first.py`, SHA-256 `79418b45566c013f4456b26e6442541a8597f2eb6d05b7790b10a4a7994f5fdd`. Code remains local, never proposed for merge.
- Unchanged driver SHA-256 `a05d0d4251ad2a326bc05bc503066286d2f295f0200b48b2ca7642454da06d91`; field measuring-stick source SHA-256 `e33e1a0e052194da09b2433360fc7bdf1789594dc86af408d4bf438d1372d2dc`. Geometry, CLI, driver, tolerances and budgets match pinned main throughout.
- Doctor-selected core 0. Harness and an actual verification worker were observed at nice 10, affinity `{0}`, without nested nice increments. One/five/fifteen-minute host load: `(1.465, 0.812, 0.341)` before, `(1.232, 0.861, 0.387)` after. No idle-host claim. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, OCP 7.8.1.1.
- Frozen public runtime input SHA-256: P3 `6691c703b7c31bb723528d1fd607e625bcfda84c5be76b682f21ef5fd3d6add4`; P3 edit `13b96162c2a2be93d04c8a73338f4993d5fe54265504eab64dd6f7a38752d031`; G1 `39054c537c577a58026fb0927134e02f1aa18abc85abf79cd2546807bc12adde`; G1 edit `b5e5caf6672e938d95056bdd82a81247b5f4201c6c5bdee7b345997e56c1edce`. P3 recipe content is unchanged apart from the explicit edit; runtime JSON is canonically serialized.
- Strict P3 STEP SHA-256: original `c623e24e8e81192284d1b6ea573e1151ab369349fbda0ea445b145896e618c37`; edited `cfb4a29b392d1ee78c07cbde5759d2b203221f7a64e0613def8b82680e9205b4`. Both equal their face-check export counterpart.

The round-start sentinel filesystem birth time is **2026-10-24 05:24:23.981862057 UTC**; the first P3 run configuration was created at **05:24:24.127222310 UTC**, after preregistration and its source addendum. Original and frozen preregistrations are preserved verbatim in the evidence. The evidence SHA-256 index is `a748067da04ecd33b738db881746eb85320df984c64f2ec11387462aa6b90417`. These timestamps record the observed local chronology.

The saved cloud workspace retains the external evidence directory `fillet-patch-first-evidence/round1`: doctor/provenance, exact commands and exit codes, frozen inputs, complete kernel requests/responses/manifests/operation logs/portable questions/STEP, four control nets, four open-patch STEP artifacts, phase results, import checks and refusal-artifact checks. The parent retains the full path; private inputs were never read or copied. Only public recipes and generated experiment results were used. Build logs and the preserved executable are beside that directory. The separate docs branch starts at the pinned main and copies only this final document; no prototype source, private corpus or unrelated work belongs in its draft PR.
