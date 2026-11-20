# Catallax and SolidDraft: Checkpoint 2 decisions

Written 20 Nov 2026 from the owner's grill session of 19 to 20 Nov, which worked through the Checkpoint 2 grill pack, a private working source. This file records every decision, the specification each new milestone needs, the roadmap edits, who builds what, and the findings behind them.

Research notes behind every claim are in [docs/research/cp2-grill/](../research/cp2-grill/). Section 11 lists them. Where a fact came from a search snippet or memory rather than a read source, the note says so.

## 0. What firstmate does with this file

1. Copy this file and [docs/research/cp2-grill/](../research/cp2-grill/) into the Catallax repository, for example `docs/checkpoints/checkpoint-2-decisions.md` and `docs/research/cp2-grill/`.
2. Write the ROADMAP.md update described in section 6 as one PR, plus the small edits in section 7. The owner approves that PR before any worker starts (decision Q29).
3. Then dispatch the work in section 8, in the order of section 5.
4. Nothing in this file reorders the existing roadmap except where the owner said so here (Q18, Q20, Q30 to Q32). Everything else is in addition.

## 1. The direction in one paragraph

Catallax becomes a tolerance-first kernel, in the way Parasolid and OpenCASCADE are. The goal is the most efficient way to build geometry that works downstream, not the strongest proof. The kernel builds fitted geometry in ordinary floating point, measures the gap between each edge's copies, stores that gap times 1.1 as the edge's tolerance within fixed caps, and checks once. Exact rational arithmetic stays where it is exact and cheap (planes, quadrics, closed-form blends) and keeps design intent such as built-in coincidence. Exact and tolerant pieces live in one model. A quick look at each operation's inputs sets a work cap for the exact attempt, and the fitted path takes over above it. The signed distance field stays: it owns lattices and field-driven design, checks the B-rep in tests, and will build every fillet on its own in time. To the user every result is approximate until they press a proof button, which comes late in the roadmap.

## 2. Decision log

Numbers follow the session. Q1 to Q4 came first; later numbers were added as questions arose.

| # | Question | Decision |
|---|---|---|
| Q1 | Direction: tolerance-first (T), measure and prove on demand (M), today's checker looser (L), S1 proof by shape kind (R) | **T**, plus an optional button that gives a user an exact proof for one part |
| Q2 | Add the 12 changes for 100% coverage | **Yes**, after the current milestones, nothing reordered; write `docs/COVERAGE.md` |
| Q3 | Drawings and PMI target | **Yes**, the P1 bracket test, checks C1 to C5 |
| Q4 | Working rules R1 to R5 | **Yes** to all five; R2 tuned per Q24 |
| Q5 | Skip checks on edge copies made from another copy | **Yes**, by tagging each copy; no storage change; CI audit of tags |
| Q6 | Checker design | **Setup A+** (section 4.2) |
| Q7 | Tolerance values | **Rules T0 to T8** (section 4.3) |
| Q8 | Watertight by construction | **Dropped**, except the exact spline fillet path, which is already watertight |
| Q9 | Which gap to measure | Own builds: same-parameter gap. Imports: geometric gap first; if the parametric gap is much larger, reparametrise instead of raising the tolerance |
| Q10 | The proof button | **B**, a real whole-edge proof, built late (item 5.2), low priority |
| Q11 | 14-edge known-answer test set | **Yes**, it is 3.0A's checker done-when. The owner's condition: the checker is general research-based code, never tuned to pass the test |
| Q12 | Exact and tolerant together | **A**, one model, exact means tolerance zero (Parasolid structure), CGAL-style filtered decisions for tolerant pieces |
| Q13 | Exact or fitted per operation | **A**, a predictor: a quick look sets the work cap tier, the cap protects; the two work-meter fixes come first; part of 3.0B |
| Q14 | Approximate labels | To the user, everything is approximate until the proof button. Inside the kernel, per-entity exact and tolerant flags as Q12 needs |
| Q15 | Field (SDF) plan | Keep the 22 Oct plan. The field may take over a B-rep step's result only when that is the only way to build the part at all and it causes no downstream problem. Long-term goal: the field builds every fillet on its own |
| Q16 | Field as a yardstick | **Tests only** (and behind the proof button). Never at user build time |
| Q17 | Field fillet milestone | **Yes**, sub-milestone 4.2a, steps 1 to 6 for analytic pairs |
| Q18 | When T work happens | **Now, before anything else**. 3.0B follows. Both before phase 3. T is design work, not speed work, so the phase 5 speed rule does not block it |
| Q19 | Numbering | **3.0A** tolerance-first core, **3.0B** one model, **4.2a** field fillets. No old number changes |
| Q20 | Drawings and PMI work | Runs mostly in parallel, but sits just behind 3.0A and 3.0B in the roadmap order |
| Q21 | Spline faces in the field | **A**, step 7 of 4.2a, including a spline-face building block for whole patches; trimmed faces in 4.4. Owner's rule: add as many building blocks as field parity needs, in the order it needs them |
| Q22 | Guards against tuning the checker to its test | **Yes**, three guards (section 4.2) |
| Q23 | What a field distance query returns | **A**, inside or outside, plus a guaranteed range for the distance near the surface |
| Q24 | App test (R2) tuning | **Yes**, kernel bumps run only the affected scenarios plus the bracket smoke test; full set nightly and at release |
| Q25 | Owner of 3.0A and 3.0B | **Fable**, one owner for both, in order |
| Q26 | If Fable quota runs short | Drawings move to **Claude Opus**; the core keeps Fable |
| Q27 | App lane | **Claude Opus**, or GPT 6.1 or 6.2 if available |
| Q28 | Reviews | **The dot (GPT Astra) reviews every PR** of 3.0A and 3.0B, from the code and research notes, not the author's summary |
| Q29 | First step | **ROADMAP.md update PR first**, owner approves, then work starts |
| Q30 | Obsolete proof work | **Remove** the "cheaper proofs" work and the 2,000 ms proof target from 5.1. Rewrite the "proof behind each result" claim and the "numerical first" rule to match T |
| Q31 | Proof lights in the app | **A**, remove the proved and bounded lights; keep the refused message and a working indicator |
| Q32 | Field preview | **Yes**, narrow 4.1: field preview only for field-only features (lattices, gyroids); the B-rep previews everything else |

Superseded earlier answers: board 1 "A, keep sampling, adaptive" lives on inside the checker; the Q1=C "predictor picks the path" answer of round 2 is replaced by Q13's specific predictor; board 2's provisional S1 is not adopted (Q1).

Withdrawn during the session, with reasons: a derivative-bound sampler ("Idea 2", too loose, about 10^5 samples per edge); watertight-by-construction everywhere (no precedent, section 9.4); "narrow the field to three roles" (it conflicted with the owner's 22 Oct rulings).

## 3. Rules the owner set or confirmed during the session

- Choose direction by logic and research, not by trying things. Measurements of existing code are fine.
- Prove nothing for proving's sake. Accept a gap and a tolerance.
- Own parts stay tight (1 nm); 0.01 mm is for poor imports only.
- The roadmap is the plan. New work is in addition. The owner reorders only by explicit decision (here Q18, Q20, Q30 to Q32).
- Everything is a trade-off; never assume a Catallax advantage.
- Checker and predictor code must be general and research-based, never tuned to a test set.
- Who builds what is decided last.

## 4. Specifications

### 4.1 Item 3.0A: tolerance-first core

Owner Fable; every PR reviewed by the dot. Starts now.

Contents, in order:

1. The 14-edge known-answer test set and its three guards (below). Built first, so the new checker is measured before it replaces the old one.
2. The A+ checker (section 4.2).
3. Copy tagging (section 4.4).
4. Tolerance rules T0 to T8 (section 4.3), including the export receiver profile.
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

### 4.2 The A+ checker

The job: find the largest gap on each freehand edge copy, reliably and cheaply, and store it times 1.1 with its location.

1. Simple pairs (line or circle on a plane, conic pairs) use a closed form. OpenCASCADE does the same (`ElementaryMaximum`).
2. Imports whose copies run at different speeds: measure the geometric gap both ways, reparametrise, then continue with the same-parameter gap (Q9).
3. Split the edge at every break: knots of the 3D curve, knots of the face-curve, every place the face-curve crosses a surface knot line, seams, poles, and joints of traced curves.
4. Per piece, take max(17, 3 x composed degree) samples, Chebyshev-spaced, offset from the fitter's own check points (the fitter uses eighths today, `numerical_embedding.rs:30-35`). Add the builder's hint location if it reports one.
5. Turn the samples into Chebyshev coefficients (a cheap transform). If the high-order coefficients are still large, double the points on that piece and repeat. This signal may only add points. It may never stop the search before step 4's count.
6. Polish every local maximum, not only the best one, with bounded Brent on the true gap, in f64.
7. Record a noise floor from the f64 error estimate instead of escalating precision. Escalate only when the gap is within a rounding error of a cap.
8. Store the maximum times 1.1 and its location. A hard per-edge work cap refuses by name.

Why the Chebyshev part helps, compared with plain sampling plus polish (the OpenCASCADE master pattern). Plain sampling picks its point count in advance from the formula's degree. Blends, sweeps, offsets and traced curves have no degree, and plain sampling cannot tell whether it took enough points. The coefficient tail answers "did I see enough": quiet tail, done at the same cost; loud tail, look harder on that piece only. Chebyshev point sets nest (17, 33, 65), so doubling reuses every point. On polynomial pieces, degree + 1 Chebyshev points reproduce the gap exactly. OpenCASCADE already trusts this kind of tail estimate in `AdvApprox` when building curves, but not in its checker. In the experiment (section 9.6), Chebyshev alone missed a narrow defect in 3 of 4 runs because it used a quiet tail to stop early. A+ forbids that. The gain of A+ over plain A is reasoned from those results, not yet measured.

### 4.3 Tolerance rules T0 to T8

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

Where this does better than OpenCASCADE, Parasolid and SMLib: own prismatic and sheet-metal parts, long operation chains (no growth, where SMLib's 2x margin can compound), mixed-size assemblies, MEMS parts, near-tangent fillets (refused clearly, not broken silently), export to tighter receivers. The same: brackets, gears, good imports, car bodies and wings against Parasolid. Worse: more refusals and opt-in clicks on poor and thin-walled imports, bigger files and slower fitting on very large freeform models.

Weakest evidence: T1's ratio (borrowed from Parasolid) and T7's numbers (no source gives one).

### 4.4 Copy tagging

Each edge copy carries `Independent` or `Derived { from, rule, bound }`, on `NumericalTrim` and `NumericalNativeTrim` (`native.rs:13-26`). Constructors must set it.

- Set by: fillet blend-face trims (the rail is a row of the blend surface), analytic line and circle embeddings, pcurves built by exact composition. Plane-projected trims stay checked, because that check is the only proof the 3D curve lies on the plane.
- Read by: the checker, which records the bound for derived copies and runs A+ on independent ones.
- On edit: any change to one copy resets it to `Independent`, unless the same exact operation applies to both.
- Export: STEP carries no tag; both copies are rounded and the export check covers rounding.
- Reopen: documents are recipes, so tags regenerate on rebuild.
- Import: everything is `Independent`.
- Wrong-tag risk: an unchecked gap ships silently. A CI audit runs the full check on every derived copy and fails if the gap exceeds its bound plus rounding, with one property test per producer.

Honest note: the checks skipped are cheap. The larger gain is that the code always knows which copy is authoritative.

### 4.5 Item 3.0B: one model, and the predictor

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

### 4.6 Item 4.2a: the field builds every fillet on its own

A numbered sub-milestone of 4.2, "the field catches up". Constant-radius rolling-ball fillets only; variable radius, setback vertex blends and conic or G2 sections keep their typed refusals.

Rule: the field derives its own ball paths. It does not reuse the B-rep's `rolling_pipe`, `rolling_saddle`, `rolling_sphere` or `rolling_mitre`, because a shared bug would fool both kernels. The milestone document names any shared dependency (the M10 pair search is one).

Field queries return inside or outside plus a guaranteed distance range near the surface (Q23). Note that `rolling_region` today answers inside or outside only; its own doc says it implies no exact distance. Steps 1 and 7 must extend it.

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

Measured on this machine (release build, i7-9750H; logs in [docs/research/cp2-grill/viability-logs/](../research/cp2-grill/viability-logs/)):

- M10 spine search 2,332,939 / 2,666,401 / 6,793,048 work on the three M11 edges, identical at 1, 2 and 4 workers, 0.14 s to 0.40 s per edge;
- full B-rep fillet 5.25M / 6.97M / 16.43M work, within 78 / 5 / 80 nm of the exact fillet.

The blocker comes first: the field cannot hold a spline face. `import_spline.rs` admits bilinear patches only. An implicit route does not scale (a bicubic patch has an implicit of degree up to 18). The model to copy is the closed-mesh leaf (`mesh.rs`): distance to the nearest face, sign by counting ray crossings. Estimate 1,500 to 3,000 lines for whole patches, about double with trimming.

Step 7 sub-steps:

- 7.0 The query contract is already set by Q23; no new decision is needed.
- 7.1 Spike: certified point-to-patch distance on the M27 patch pair, in catallax-field, timed over 10^4 points. This retires the biggest risk.
- 7.2 The spline-face building block for whole, untrimmed patches, with ray-crossing sign. Trimmed faces go to 4.4 with a typed refusal until then.
- 7.3 The field's own spine proof: prove the two end points on the face borders, then a parametric Krawczyk proof span by span. Tube width eps <= min(100 nm, part tolerance / 2). A build that cannot reach it refuses and reports the eps reached.
- 7.4 The swept-ball distance query on that tube.
- 7.5 Near-tangent sweep on a synthetic edge: record the smallest face angle proved inside 20M, then decide the tail policy (deflation per Leykin et al. 2006, or a labelled bound).
- 7.6 Spline-with-analytic pairs: M10 needs both sides as rational patches. Write the analytic offset (same family) as a rational patch with offset 0. Untested.
- 7.7 Typed refusals, including a singular-contact control.

Step 7 done when, on the three M11 edges plus one spline-plane edge:

- the field proves its own eps <= min(100 nm, tolerance / 2) inside 20M work, identical at 1, 2 and 4 workers;
- every B-rep display vertex of each fillet face lies inside the field range widened by the B-rep's reported bound;
- a dense-sampling third check in tests agrees (catches a bug shared through M10);
- median query time is recorded;
- a near-tangent control refuses or returns a labelled bound.

Why a proven tube and not an estimate: the field is the yardstick for B-rep bugs (Q16). An estimated path cannot tell whether the B-rep path is wrong. It runs in tests only, and it proves only down to the tolerance, never past it.

Caveat: the three M11 test edges are single-span bilinear or biquadratic patches, so they say little about real imported parts.

### 4.7 Item 5.2: the proof button

A real whole-edge proof on demand for one part (Q10 B). Late and low priority. Built from parts Catallax already has: exact range arithmetic (`evaluate_box`, `jet_box`, Ball and Interval) and a one-dimensional branch and bound that splits wherever an upper bound beats the best sample. Nothing in Catallax computes a certified maximum of |C(t) - S(P(t))| today. Expect minutes per part and "could not prove" on blends and sweeps (solidcraft's experience: one edge took 153 s, some bounds came out 60x loose, blends were not provable). VISION.md's "no proof jobs or proof services" still holds: the button is one computation on request, with nothing stored and no service.

### 4.8 Drawings and PMI (the P1 bracket test)

Spec unchanged and stored in [checkpoint-2-drawings-pmi.md](checkpoint-2-drawings-pmi.md), copied from the review report and grill pack. It covers manufacturing information written into the STEP (AP242 semantic PMI), automatic drawings from it, feature recognition on bare imports, and the interview. Owner Fable; if Fable quota runs short, it moves to Claude Opus (Q26). Runs mostly in parallel with 3.0A and 3.0B, but sits just behind them in the roadmap order (Q20). The bracket is all planes and cylinders, so T does not change it.

Checks C1 to C5, unchanged: C1 round trip of all 14 items (NIST analyzer, OCCT reads all 14, re-import, save, reopen); C2 bore edit 8 to 10 mm keeps the tolerance attached; C3 one button gives STEP plus PDF and DXF drawing with every item once, no overlaps, no hand moves, title block; C4 nothing on the drawing absent from the file; C5 bare import, recognition finds bore, boss, 4 holes, slot, the interview asks only for what is missing.

### 4.9 Working rules R1 to R5, with the R2 tuning

- R1: one feature, one owner across kernel and app; parallel lanes only on different features.
- R2: a kernel change to an app-used feature passes that app acceptance, with one later feature on top, before merge. Tuning (Q24): each SolidDraft scenario lists the kernel parts it uses; a Catallax pin bump runs only the scenarios whose kernel parts changed, plus the bracket scenario as a smoke test; the full set runs nightly and before every release; shared core changes or scenarios with no list run everything. Today any pin bump runs all eight Linux units and all 85+ scenarios (`soliddraft/docs/ci.md`, line 9).
- R3: a target that fails two different approaches is parked and brought to the owner.
- R4: once a release is called, only wrong-result, crash and broken-feature fixes go in.
- R5: scoreboard and census once per checkpoint, no dry runs; coverage counted from `docs/COVERAGE.md`.

### 4.10 Coverage: the 12 changes for 100%

Write `docs/COVERAGE.md` with the full bar for each of the 61 lines, and make every done-when line name the coverage lines it closes and at which level. Then:

1. Widen 2.2: least-squares fitting, G1 N-sided and Gordon fills, curve networks with tolerance, lofts with guide curves (C7, C8, C9, C10, F8).
2. New line: offsets of any curve, surface or solid, and projection of curves onto surfaces (C6, C11, E6).
3. Widen 3.4 to general input: bend and twist, curvature-continuous blends, mesh Booleans and decimation (C12, E8, Q8).
4. New line: intersections complete, every analytic pair, spline against spline, curve against surface, tangent contact under the bounded-answer rule (I1, I3, I4, I5).
5. Widen 1.3 and 2.1: Booleans and fillets on spline-faced solids pass a generated set and the census parts (S3, E3).
6. Widen 2.3 and add split and imprint (S6, S7, X7, S5, F10).
7. Widen 2.4: taper, guide, twist, helix and up-to-face on any profile along any path (F4, F5, F6, F7, F9, F11).
8. New line: edge and wall features in general, variable-radius fillets, chamfer roll-over, shell and draft on any face (E1, E2, E4, E5, E7).
9. Widen 2.5: every direct edit on any face type (E9).
10. New line: STEP import parity against OpenCASCADE's 90.03% (X2, X4).
11. Widen 3.2: colours, names and annotations round-trip (X8).
12. New line: geometry queries, closest point and curvature analysis, and widen 3.3 so no pair answers "unresolved" (Q4, Q5, Q6).

Full text in the grill pack, Part 3, "Smallest set of additions for 100% of both kernels". The count today: 61% recounted on 19 Nov; about 67% if every milestone meets its done-when, 85% at full planned scope.

## 5. Order of work

1. ROADMAP.md update PR, owner approval (Q29).
2. 3.0A tolerance-first core (Fable).
3. 3.0B one model and predictor (Fable).
4. Drawings and PMI: runs alongside from the start where it does not touch 3.0A/3.0B code; sits just behind 3.0B in the roadmap order.
5. App lane alongside (Opus or GPT 6.1/6.2): R2 CI tuning, removal of proved and bounded lights, 4.1 preview narrowing, then the app side of 3.0B (mixed exact and tolerant bodies, no refusal after a fillet).
6. Phase 3 items 3.1 to 3.4, after 3.0B.
7. Phase 4, with 4.2a inside 4.2.
8. Phase 5, then 5.2 the proof button.

## 6. ROADMAP.md edits

Add:

- 3.0A and 3.0B before 3.1, with the contents and done-when lines of sections 4.1 and 4.5.
- 4.2a inside 4.2, with steps 1 to 7 and done-when lines of section 4.6. Any extra building block field parity needs is added in the order it is needed (Q21).
- 5.2 the proof button (section 4.7).
- The drawings and PMI item, just behind 3.0B (section 4.8).
- The 12 coverage changes (section 4.10). Checkpoint 3's bar moves from "every line has something built" to the new done-when lines.
- Rules R1 to R5 with the R2 tuning (section 4.9).
- Tolerance rules T0 to T8 (section 4.3) under "Rules that stay fixed between checkpoints".
- The new lane assignment (section 8).

Remove (Q30):

- From 5.1, the "cheaper proofs" bullet (proofs from construction facts, prepared bounds, reuse across edits, coarse proofs refined where they fail) and the "2,000 ms for proof" click target.

Rewrite (Q30, Q32):

- The "numerical first" rule under "How we run it", to T: build fitted geometry in f64, measure the gap, store 1.1 x gap within the caps of T0 to T8, check once, write STEP marked approximate; the exact path is chosen by the predictor; the proof button comes at 5.2.
- The line "Twelve things neither reference kernel does ... a proof behind each result": drop "a proof behind each result". Keep "the same answer on every machine" (work units keep this true) and "exact numbers" for exact pieces.
- 4.1: the app draws field previews only for field-only features (lattices, gyroids); the B-rep previews everything else. Drop the proved and pending light wording (Q31).
- 4.2: add "the field may take over a B-rep step's result only when that is the only way to build the part, and only with no downstream problem" (Q15).
- "Delivery until phase 5 ... no speed work": add that 3.0A and 3.0B are design work and not blocked by it (Q18).
- The rule "Proven, bounded and refused answers stay separate" stays inside the kernel; add that the user sees results as approximate until the proof button (Q14).

Move to history:

- The phase 2 lane setup ("three kernel lanes", "Kernel lanes run on Codex models ... Codex cloud reviews and advises").

Keep unchanged: phase 2 open items 2.1, 2.4 (P8 parked), 2.5; phase 3; 4.3 to 4.5; the 20M target and 100M allowance; "kernel changes are general"; "no tolerance loosened to make one case pass"; "clock time never fails a change before phase 5".

## 7. Other file edits

- `data/captain.md`: "tolerant import default, certified the higher label" becomes "tolerant import default; the proof button gives the proof". Add the rules of section 3.
- VISION.md: no change needed. "Numerical where needed, honest everywhere", the field checking "on demand", and "no proof jobs or proof services" all fit the decisions.
- SolidDraft: remove the proved and bounded lights (Q31); keep the refused message and a working indicator; R2 scenario tagging and selection (Q24); narrow field preview (Q32).
- Catallax PRD: the rule "every recipe operation has a field version" stays as the long-term goal (Q15), with 4.2a as the route; it never blocks a B-rep feature.

## 8. Who builds what

| Work | Owner | Review |
|---|---|---|
| 3.0A and 3.0B (one owner, in order) | Fable | the dot (GPT Astra), every PR |
| Drawings and PMI | Fable; moves to Claude Opus if Fable quota runs short | the dot |
| App lane | Claude Opus, or GPT 6.1 or 6.2 if available | the dot |
| ROADMAP.md update | firstmate | owner approval |

Outside CLIs were unavailable at the time of the session. Reviews therefore use the three in-house models.

## 9. Findings worth keeping

### 9.1 What the tolerances really were

- Exact parts have no tolerance (not 1 nm). The 1 nm figure only limits an unproved coincidence decision on the certified-offset path.
- Own-built approximate parts are held to 1e-7 mm (0.1 nm), ten times tighter than the owner remembered (`evaluate_numerical_fillet.rs:115`); the fillet API refuses above 1 um.
- 0.01 mm applies to imports, plus an unused library default (`NumericalOptions::default()`, PR #603, no reason recorded).
- No receiving CAD system works as tight as 1e-7 mm. Parasolid's docs warn that precision of 1e-6 to 1e-8 m makes SP-curves with many segments and slows things down.

### 9.2 Why checking dominated cost

"Numerical" only replaced construction. The fitted result was stored as exact fractions and checked in 128-bit software arithmetic at 9 points per edge at build, again at export and again at reopen. On the loft fillet, checking was 96.7M of 97.4M work. Rim export cost about 22M against about 6M to build. OpenCASCADE and SMLib check in floating point, on demand, once, and export only writes numbers. Today's 9-point check also proves nothing between points, and the fitter aims at the same points the checker samples.

### 9.3 Why no kernel proves the gap along a whole edge

- OpenCASCADE: the plain sampling check alone took 17.9% of a Boolean (issue #1372); its exact checker became opt-in in 2021 (0032448).
- SMLib: a TODO in `SmGap.cpp` says solving for the true maximum "would probably be much more expensive".
- Parasolid's tolerances exist to absorb the gap (Jackson 1995).
- solidcraft (open-source Rust CAD, SMLib port) tried proven bounds: one edge took 153 s, some bounds came out 60x loose, blends were unprovable.
- A derivative bound from control nets bounds each copy separately and loses: about 10^5 samples on a 100 mm edge at 1e-3 mm margin. The bound that works composes the gap into one polynomial (degree 18 for bicubic with cubic) and loses about 9 of 16 digits to cancellation.

### 9.4 Watertight by construction

The exact intersection of two bicubic patches has degree 324; a trim curve on a bicubic patch has degree at most 18. So no low-degree edge lies on both surfaces, and every watertight method reshapes surfaces near every edge. In the one industrial test found (318 faces), all 312 analytic faces became splines, control points grew 23.8%, and fillets lost tangency. No shipping solid kernel is watertight by design: Fusion Form (T-splines), Rhino SubD and nTop all convert to tolerant B-rep for STEP. Autodesk bought T-Splines in 2011 and later dropped the SolidWorks and Rhino plugins. STEP AP242 edition 2 has a locally refined spline entity, but no translator was found that reads or writes it. The research was for isogeometric analysis. Catallax's exact spline fillet is already watertight by construction, which is the useful case.

### 9.5 One authoritative copy per edge

Parasolid and ACIS keep one authoritative copy per edge. An exact edge stores only its 3D curve (Parasolid evaluates intersection curves to any precision). A tolerant edge stores only SP-curves on each face, and the edge tolerance is the gap between them. Catallax's exact bodies already store only the 3D curve; its numerical bodies keep independently fitted copies that consumers read, which is why tagging (Q5) and not deletion was chosen.

### 9.6 How sampled checkers fail, and the experiment

Recorded misses: OCCT 0025597 (2015, Boolean face-curves never measured), 0025613 (2015, one global search found a local maximum), 0028030 (2016, above 101 knots fell back to 101 even points, which later let a peak through), 0029611 (2018, the Hessian was the gradient, wrong for about 3 years with no known-answer test), PR #1640 (Oct 2026, peaks set by surface knots, high-degree peaks, rational end weight, origin at 1e12), PR #1683 (open, a peak exactly between two equal samples); SMLib solidcore #3258 (tolerance 3x a 20-point sample still below the true gap) and #570 (an end-point edit stretched the parameter by 1.0000057); solidcraft PR #18166 (face-curve outside its domain lifted 5.6e13 off). Lessons: more points fixed none of these; splitting at all breaks plus polishing every local maximum did.

Synthetic experiment ([docs/research/cp2-grill/exp/](../research/cp2-grill/exp/), Python, evaluation counts not time): 23 uniform points found 16% to 84% of the true maximum; split at breaks plus 54 points per piece plus polish of every local maximum caught a defect 1e-4 of the edge wide; Chebyshev alone was 1.3 to 2.4x cheaper but missed that defect in 3 of 4 runs. Every method misses below some width (3e-5 here).

### 9.7 Tolerance numbers in the wild

PDQ standards (VDA 4955, SASIG, ISO 10303-59) name the checks and leave every value to the company. Receivers: NX, Creo, Rhino and Solid Edge work at 0.01 mm; CATIA merges at 0.001 mm and heals up to 0.1 mm; CATIA STEP files declare 0.005 mm. Kernels: OpenCASCADE Confusion 1e-7 mm, Approximation 1e-6 mm; Parasolid session precision 1e-5 mm, exact edges count as half; ACIS resabs 1e-6 mm, resfit 1e-3 mm. SMLib sets edge tolerance to 2x the gap while the gap is under the bounding-box diagonal / 800, with no stated reason; 2x compounds over operations, and /800 allows 1.25 mm on a 1 m part. No source shows a receiver rejecting a file for being too tight; the documented failures are real gaps larger than the declared tolerance.

### 9.8 The field (SDF), honestly

- It compiles the recipe, never the B-rep. It supports stock, extrudes, revolves, holes, patterns, Booleans, simple fillets (plane-plane, plane-cylinder, three-plane corners), chamfers, analytic shells. It refuses curved-curved fillets and the P5 radius-2 fillet. After a Boolean it gives a residual, not a true distance. It builds P1 in 710 work where the B-rep runs out of budget.
- It is 76.6k lines, about 13% of the code, and fills phase 4.
- Why fillets are hard in the field: a rolling-ball fillet is everything at distance r from the ball-centre path, and that path is where the two faces' r-offsets meet. Plane-plane gives a line, plane-cylinder a circle or ellipse, two cylinders a degree-4 curve, spline faces no formula. The hard part of a fillet is the same in both kernels. Smooth-min blends are a different shape; they match a rolling ball only between perpendicular planes.
- Copying the B-rep fillet into the field is possible but loses that step's independent check and makes each query a nearest-point search.
- Under T, the field is not needed at build time, as Parasolid and OpenCASCADE show. Its yardstick value is differential testing during development: the kernel can build a wrong shape that passes every gap check, and big vendors catch that with decades of test suites and millions of users that Catallax does not have.
- Meshless FEA needs inside or outside plus where the surface is; for a machined part that field can come from the B-rep, so FEA does not wait for field fillet parity. Intact.Simulation and SimScale run on nTop implicit geometry (search snippets).
- Parasolid v38 added implicit lattices and v39 describes B-reps, facets and implicit geometry in one model (Siemens blog, search snippet only).
- Fornjot is a pure B-rep kernel; its author did SDF experiments before it and hit a cliff on tangent, coincident and tolerance cases. Catallax's two-kernel plan came from the roadmap written after that retrospective (22 Oct).

### 9.9 The CGAL way, explained

Geometry code builds (new curves, points) and decides (yes or no: inside, touching, order). Small building errors give slightly wrong shapes; wrong decisions break topology. CGAL builds in floats and makes every decision exact on the stored numbers with a filter: float with an error bound, then interval, then exact only on close calls. Typical cost 8% to 35% over plain floats, up to 11x on near-degenerate input (Shewchuk 1997); about 4x for CGAL's lazy exact kernel. Catallax today does both building and deciding in exact rationals, and rationals grow with each step. Q12 keeps exact rationals where they hold design intent and uses CGAL's filter for every decision with a tolerant piece.

## 10. Open items, not verified, and measurements to run

- 3.0A: measure the gain of A+ over plain sampling plus polish on the 14 edges (reasoned, not measured).
- 3.0A: cost per f64 surface evaluation in Rust (100 to 300 ns is a guess).
- 3.0B: exact entities with no tolerance or a 5e-6 mm floor; measure the M8 cuts under both.
- 3.0B: the 2M and 20M predictor tiers; the fitted cost of the P1 fillet-end vertex; refunded exact cost on the rim; STEP work once `fit.rs` charges it.
- 4.2a step 7: per-query cost of spline-face distance (no code exists yet); how often near-tangent fillet ends occur on real parts; the analytic-to-patch wrapper.
- The dense re-check the S1 condenser proposed (rim, loft fillet, P5) was never run. It would show whether today's 9-point check ever passed a bad edge. The 14-edge set in 3.0A partly replaces it.
- Several sources were search snippets or memory (marked in the notes): Parasolid and ACIS internals, Sederberg 2008, Urick 2019 full text, the Parasolid v38/v39 implicit claims, Intact and SimScale.
- solidcraft PR #18166 merge status unknown.

## 11. Research notes (in [docs/research/cp2-grill/](../research/cp2-grill/))

| File | What it holds |
|---|---|
| [t-smarter-research.md](../research/cp2-grill/t-smarter-research.md) | Smarter-than-plain-T survey: one authoritative copy, derivative bounds, caching, repair before growth, watertight, theory |
| [edge-copies-code.md](../research/cp2-grill/edge-copies-code.md) | Which edge copies Catallax stores, derived or independent, who reads them (code, file:line) |
| [edge-bound-research.md](../research/cp2-grill/edge-bound-research.md) | Derivative and coefficient bounds, applicability by geometry type, why no kernel proves the gap |
| [tolerance-numbers-research.md](../research/cp2-grill/tolerance-numbers-research.md) | Kernel defaults, receiver tolerances, PDQ standards, SMLib rules |
| [tolerance-scenarios.md](../research/cp2-grill/tolerance-scenarios.md) | T0 to T8 against every part scenario, better, same or worse than others |
| [watertight-research.md](../research/cp2-grill/watertight-research.md) | Watertight-by-design precedents, cost, STEP, verdict |
| [fillet-watertight-code.md](../research/cp2-grill/fillet-watertight-code.md) | How Catallax builds fillets; watertight fillets per support kind; tagging design |
| [checker-menu.md](../research/cp2-grill/checker-menu.md) | Every checker method, the experiment, setups A, B, C |
| [checker-failures.md](../research/cp2-grill/checker-failures.md) | Recorded checker misses in OCCT, SMLib, solidcraft; Catallax risky geometry; reusable tools; the 14-edge set |
| [exp/gapcheck.py](../research/cp2-grill/exp/gapcheck.py), [exp/variants.py](../research/cp2-grill/exp/variants.py) | The synthetic checker experiment |
| [occt1640.diff](../research/cp2-grill/occt1640.diff) | The OpenCASCADE Oct 2026 checker fix, for reference |
| [exact-tolerant-coexist.md](../research/cp2-grill/exact-tolerant-coexist.md) | How exact and numerical bodies interact today; options A to E for one model |
| [sdf-role.md](../research/cp2-grill/sdf-role.md) | What the field does, what works, roles, cost |
| [sdf-fillet-research.md](../research/cp2-grill/sdf-fillet-research.md) | Field rolling-ball fillets for analytic pairs; milestone steps |
| [sdf-spline-fillet-research.md](../research/cp2-grill/sdf-spline-fillet-research.md) | Routes for spline-face field fillets |
| [sdf-spline-viability.md](../research/cp2-grill/sdf-spline-viability.md) | Viability study with measured runs; risks; plan |
| [viability-logs/](../research/cp2-grill/viability-logs/) | Raw logs of the measured M10 and spline fillet runs |
| [predictor-research.md](../research/cp2-grill/predictor-research.md) | Predictor designs, how each chooses on the known cases, the meter holes |
