# How sampled edge-gap checkers fail, and what Catallax has to build on

Study date 2026-11-20. Catallax worktree Catallax repository at bb2f8375.
Labels: [read] = I read the source, commit or PR text; [snippet] = search-result summary only;
[unread/memory] = not checked this session.

## Part 1: failure evidence

### OCCT, from commit messages and diffs (gh api, [read])

| Case | Geometry that caused it | Fix | Lesson |
|---|---|---|---|
| 0025597 (Jan 2015), commit 1b7ae95146 "Invalid curve on surface in the result of General Fuse" | Booleans built 2D curves and intersection curves (incl. plane/B-spline sections) whose lifted copy left the 3D curve; nothing measured it | Moved the max-distance tool into `BRepLib_CheckCurveOnSurface` and ran it on every pcurve the Boolean built; set edge tolerance from it | The producer must measure what it built; a pass/fail check run later is too late |
| 0025613 (Jul 2015), commit 769fb6a4ab "Wrong distance found by xdistef" | B-spline edges: one global search over the whole range found a local, not global, maximum | Split the range at knots; run particle swarm (math_PSO) then local refine per knot interval; parallelised | Search per knot span; one global optimiser misses peaks |
| 0025819 (Feb 2015), e002f1ce5b "Bad result of BOP cut on valid shapes" | Parameters adjusted at surface boundaries (periodic); check ran outside First/Last | Clamp First/Last in `Compute()` | Range handling on periodic and bounded carriers is part of the checker |
| 0026918 (2016), 309bad284d | Periodic 3D curve: checker could not compute a distance at all | Test case added (fix in earlier work) | Closed/periodic curves need explicit tests |
| 0028030 (Nov 2016), e841c38c71 | B-spline with very many knots: per-span search too slow | Above 101 knots, drop the knot grid and use 101 uniform sample points | A cost cap silently changed the method. 10 years later (below) this cap was the hole a narrow peak fell through |
| 0029611 (Mar 2018), 8e45500e21 Hessian fix | Any curve: Newton refinement used `Hessian = gradient` (diff shows `theHessian(1,1) = theGrad(1)`), so local refinement after PSO converged badly | Real second derivative `-2(|dC2-dC1|^2 + (P2-P1).(ddC2-ddC1))` | A refinement step can be wrong for years if no test has a known maximum to compare against |
| 0032448 (Oct 2021), 0ffecc2fc7 | Default `checkshape` samples a finite number of points; exact method was optional | Added `-exact` option (BRepLib_ValidateEdge). New tests on real files (bug27814.brep: 88 faulty; a sewn file: 8234 faulty) | OCCT itself labels the default "faster, possibly not correct". Default stayed sampled |
| PR #1640 (merged 18 Oct 2026), "preserve narrow extrema" | GTest `DenseNonuniformKnotsPreserveNarrowMaximum`: degree-1 spline, 100/101/102/128 knots, a 1e-6 peak spanning 2^-30 in parameter, periodic and shifted by a period. Also `SurfaceKnotsPreserveNarrowMaximum` (peak 2^-12..2^-42 wide set by **surface** knots, not curve knots), `SurfaceDegreePreservesInteriorMaximum` and `SeparateSurfaceMaximaAreRefined` (max-degree Bezier/B-spline surface where a narrow Bernstein peak beats the central one by 1e-6), `CurveResidualPreservesEndSpanContinuation` (rational weight 2 on an end pole, origin at 1e12) | Rewrote the search: extract polynomial pieces per span, cut at **surface** knot crossings too, sample each interval, bracket and refine every distinct local maximum, report an evaluated distance above a limit | Peak can be set by curve knots, surface knots, high surface degree, rational weights or a far origin. Sample count alone does not fix any of these |
| 8.0.1 #1373 (Jul 2026), d8b85a6626 | Line/line, circle and ellipse compositions on elementary surfaces | Closed-form maxima (endpoints for lines, trig roots for conics) instead of PSO | Exact formulas for the common elementary cases are cheaper and correct |
| PR #1683 (open, Oct 2026) surface boxes | Best sample exactly halfway between two equal samples; rounding decided whether to refine. Seam: the two sides evaluate 1e-14 apart, so a tie was refined on the wrong side | `< d + eps` instead of `<` | Tie-breaking and seams decide refinement; "refine only where the sample is best" loses symmetric peaks [read] |

Not reached: tracker.dev.opencascade.org did not resolve (DNS). Ticket pages unread; commit texts read instead.

### BRepCheck default behaviour [snippet]
- `ShapeAnalysis_Edge::CheckSameParameter` and BRepCheck use 23 equally spaced points by default. Deviation between points can be missed.
- unmesh PR #148: `BRepLib::SameParameter` absorbed large deviations by growing edge tolerance to 1.05e-2. Shape then passes. Same "tolerance = measured gap" design as Catallax's plan; pass says nothing about size.
- OCCTSwift #2929 / #2743: flags stale or ignored, so checker skipped measurement; 1.07 on R5 circle, 92 on a 10 box edge.
- FreeCAD via solidcraft #8502: `isValid()` True, `check(True)` raises InvalidCurveOnSurface on circular edges written as rational B-splines with no pcurve on a cylinder (8.1e-5 vs 1e-7).
- build123d #1512: UnifySameDomain on a sphere seam concatenated pcurves (ConcatC1); edge tolerance inflated to 0.8.

### SMLib / solidcore / solidcraft
- solidcore #570 [read]: thread flank trim (helical sweep thread cut by Boolean, cylinders R8 and R5) lifted 2.54e-5 off its edge vs tolerance 2e-5. First blamed on 5 points per span. Real cause: `sm_FixUVTrimCurves` trimmed then `EditParameterization()` stretched the trim by knot ratio 1.0000057, sliding every interior parameter (5.05e-5 at same t). Fix: move only the end point. Lesson: an end-point edit that reparametrises shows up as a same-parameter gap; the checker must sample after every edit, and a later cache-building step can create the gap.
- solidcore #3258 [read]: projected pcurves on the NURBS form of cylinder and sphere faces wander 8.376e-5 from the edge vs tolerance 8.266e-5. SMLib samples 20 points; `cover_gaps` then sets tolerance to 3x that sampled gap. Tightening the requested tolerance 8x returned the same pcurve. Found with even points + 4 per knot span of both curve and pcurve + golden-section at peaks. Lesson: tolerance set from a sampled max is under by the between-sample excess, and the x3 factor did not cover it here.
- solidcraft #13630 [read]: traced intersection edges store measured deviation (5 interior samples per span) as tolerance; no proved bound; proposal is grid error + interpolation remainder from derivative enclosures + Krawczyk.
- solidcraft PR #18166 [read, first half]: old rule raised edge tolerance to 2x reported error per refit, so tolerance ratcheted (1.37e-3 vs proved gap 2.7e-4). Out-of-domain trim (u = 71.94 on a B-spline face) lifted through surface **extrapolation**, 5.6e13 off. Blend-face and plane-line trims rest on samples alone. A proved bound 0.229 vs sample 3.6e-3: bound too loose to use. Bisection stalled at 0.098 for a 2.2e-4 gap.
- solidcraft PR #19050 [snippet]: an equation curve 25x its tolerance off passed silently when the fit could not be proven.
- No paper found that reports a miss rate for sampled checks [snippet]. Evidence is case-based.

## Part 2: Catallax inventory (bb2f8375, [read] unless marked)

### 2.1 Geometry that reaches the numerical checker
Checker accepts curves Line, Circle, BSpline (`catallax-topology/src/numerical/source.rs:62-76`) and surfaces Plane, Cylinder, Cone, Sphere, Torus, BSpline (`numerical/analytic.rs:39-48`); source path accepts Plane and RationalPatch / one-chart RationalAtlas (`source.rs:92-100`). Ellipse, Intersection, PatchSection, PatchPair, TorusSheet, EllipticCylinder are refused `numerical.unsupported-carrier` or converted before (embedding fits PatchSection/PatchPair to cubics, `numerical_embedding.rs:1-13`).

| Producer | 3D copy | Face copy | Risk properties |
|---|---|---|---|
| Embedding (sweeps, lofts, holes, numerical STEP) `catallax-features/src/numerical_embedding.rs:688-696` | fitted cubic, up to `MAX_SPANS`=256 spans (`:38`) | fitted native UV B-spline per face | many knots; fitter tuned at `CHECK_SHARES` eighths (`:39`) = same points admission checks |
| Numerical fillets `catallax-shapes/src/spline_fillet/numerical_carrier.rs:179-237` | rational rail rows of blend net (weights `:71-73`, `:195`) | non-rational cubic support pcurves, endpoint-edited (`fillet_numerical.rs:836-848`) | rational vs non-rational copies; end-point edits (solidcore #570 pattern); R3.99 rim: 324 stations, span 4.73e-4 |
| Numerical Booleans `catallax-boolean/src/numerical/partition.rs:118-163` | traced section fit | independent UV fits | traced curves near tangent faces (no proof of delta) |
| Strict STEP import `catallax-import/src/numerical.rs:399-460` | file Line/Circle/BSpline | linear UV between inverted end points, or projection on planes | any degree/knot count from file; circles on cylinder/cone/torus use angle parameters; linear-UV assumption only established by the 9 points |
| Tolerant import (separate checker) `tolerant/admission/fit.rs:21` | file | file PCURVE or rebuilt | 17 binary64 points; tolerance grows to cover measured |

Risky properties present: periodic faces with paired seams (`analytic.rs:128-155`, `analytic/periodic.rs:39-105`); sphere poles refused (`analytic.rs:460-468`, so no pole case today); closed circles carry a rational approximation of one period (`analytic.rs:23-27`); rational weights (fillet rails, RationalPatch faces); many knots (embedding up to 256 spans; fitter densities to 256); planar projected trims are an on-plane check only.

### 2.2 Reusable tools
| Tool | What it computes | Cost | Where |
|---|---|---|---|
| `BSplineCurve::evaluate_box` / `jet_box` | Certified enclosure of position (and 1st/2nd derivatives) over an exact rational parameter interval | One span pass in Ball arithmetic at chosen precision | `catallax-shapes/src/bspline/curve.rs:475-499` |
| B-spline surface `evaluate_box` / `jet_box` | Same over a rational u-v box, covering every crossed knot span | Per covered span | `catallax-shapes/src/bspline.rs:628-659` |
| Ball, Interval | Midpoint-radius balls, intervals, rationals | 128-bit used by admission | `catallax-numbers/src/ball.rs`, `interval.rs` |
| `catallax-shapes` distance branch-and-bound | Certified **minimum** distance, cell cap 4096, returns attained enclosure on exhaustion | Cells x surface eval | `catallax-shapes/src/distance.rs:1-5, 57-67`; **refuses B-spline and TorusSheet faces** (`:163`) |
| `catallax-measure` PreparedFace | point-to-face, face-to-face certified distances | as above | `catallax-measure/src/distance.rs:20-240` |
| Mesher Hessian / chord bound | Upper bounds on S_uu, S_uv, S_vv per chart; chord error <= Hessian-weighted variance / 2 | Built once per chart | `catallax-mesh/src/chart.rs:291-312, 405-431` |
| STEP curve writer | Halve while proven Taylor remainder > tol; stop when rounding dominates | Reserves work per span | `catallax-step/src/fit.rs:1-2, 164` |
| Import UV trace | Error-driven halving, cap 65,536 points, named refusal | binary64, cheap | `catallax-import/.../fit/trace.rs:24-65` |
| Fitter density doubling | 4..256 per piece until binary64 fit meets tol/4 | build side | `numerical_embedding.rs:36-39, 520` |
| SDF / field distance | Distances from the field's own recipe predicates, no B-rep input | n/a | `catallax-field/src/distance.rs:1`; not reusable for edge copies |

Key gap: no tool computes a certified **maximum** over t of |C(t) - S(P(t))|. All parts exist: interval of t -> `evaluate_box` on C and P -> surface `evaluate_box` on P's box -> upper bound of distance; split where the upper bound exceeds the current best sampled max. That is a 1D branch-and-bound, simpler than the 2D min search.

### 2.3 Known-answer tests
None found. The closest is `linear_native(..., spike)` in `catallax-topology/tests/numerical_trim/native.rs:84-108`: a 1e-4 spike placed **on a knot** (t = 1/16), which knot-span sampling hits by construction. No test places a gap between sample points or states a true maximum.

Proposed set (each with a closed-form or exact-rational true max; run the checker and assert reported >= true max and <= 1.1 x true max):
1. Degree-1 spline, narrow peak between eighths (e.g. at 5/16 + 2^-30), knot count 100-300 (OCCT #1640 case).
2. Peak set by surface knots, not curve knots (curve has one span).
3. High-degree (e.g. 25) Bezier surface with a narrow Bernstein peak near an end beating the central bump by 1e-6.
4. Rational: end pole weight 2, and a fillet-style rational rail vs non-rational cubic pcurve.
5. Far origin (1e6..1e12 mm) with small gap.
6. Peak exactly midway between two equal samples (tie case, PR #1683).
7. Periodic: edge across a cylinder or torus seam; same edge shifted by one period; full closed circle with the rational period.
8. Parameter-speed mismatch: pcurve reparametrised by a 1.0000057 knot stretch (solidcore #570), geometric gap small but same-t gap 5e-5.
9. Pcurve leaving the face domain (extrapolation, PR #18166).
10. Near-tangent faces: traced section where faces meet at under 5 degrees.
11. Fitted cubic that passes at every eighth but fails at sixteenths (fitter/checker coupling).
12. Many knots: 256-span fitted edge (R3.99 rim) for work-cost budget.
13. Zero-gap controls: line on plane, circle on cylinder, iso-line rail (must report ~0, catch false alarms).
14. Analytic exact answers: circle vs its rational-B-spline copy with known chord sag.
