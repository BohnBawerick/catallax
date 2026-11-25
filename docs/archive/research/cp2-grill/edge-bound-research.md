Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Certified edge gap bounds ("Idea 2"): research notes

Date: 2026-11-19. Labels: [read] = I read the abstract, page, source file or PR text in this session.
[snippet] = only a search-engine summary, not the page itself. [unread/memory] = from memory, not checked.

Local copies of source files read: `public upstream source files, not bundled here`
(BRepLib_ValidateEdge.cxx, GeomLib_CheckCurveOnSurface.cxx from OCCT master; SmEdgeuse.cpp, SmGap.cpp from NVIDIA-Omniverse/usd-brep main).

## Correction to the brief first

- OCCT BRepLib_ValidateEdge: `myControlPointsNumber(22)` and a loop `i = 0..22`, so 23 points. Confirmed [read].
- SMLib `SmEdgeuse::CalcMaxEdgeFaceGap`: `lSampleCnt = 20`, builds an `SmCrvSrfGapFunction` from the 3D curve to the *surface* (not through the pcurve), then takes the max sample [read]. `SmGapFunction::BuildSampleSet` notes that B-spline curves also get mid-knot samples "since many BSplineCurves (e.g. XSectCurves) are built to maximize their gap sizes at the mid-knot points" [read].
- OCCT master changed `GeomLib_CheckCurveOnSurface` in Oct 2026 (PR #1595, #1640): PSO + Newton replaced by "bounded Brent", plus a `CompensatedReal` double-double type and a `BSplineResidual` class that evaluates C(t) - S(p(t)) "without rounding intermediate points" because "this is needed when its final difference is small compared with either point" [read]. Still a search, still no bound between samples.
- "No kernel bounds the gap between samples" is not true any more. tangency-labs/solidcraft (Apache-2.0 Rust CAD app, ~17k stars, created Oct 2026, contains an SMLib port `crates/smlib`) already computes proved upper bounds on the edge/trim gap with rounded interval arithmetic and uses them to set edge tolerances (PR #18166, open as of 2026-11-16, and related merged-looking PRs) [read]. Details in Q4.

## Q1. Is it based on published research?

Yes. Every component is published; the specific combination (certify a stored edge's curve-on-surface gap, used as a kernel validity check) is not in any paper I found. It shows up only in an open-source kernel's PRs.

### Derivative-bound sampling and Lipschitz optimisation
- Filip, Magedson, Markot, "Surface algorithms using bounds on derivatives", CAGD 3(4):295-311, 1986 [snippet: abstract via ResearchGate summary]. Abstract: generalises piecewise linear approximation, min-max boxes and surface/surface intersection to all parametric C2 surfaces using "theorems from approximation theory which give information about the maximum deviation an approximation to a surface can have if bounds on partial derivatives are known." Note: the Filip bound uses *second* derivatives (linear interpolation error <= M2 h^2/8), not first. This matters for sample counts (see Q3).
- Tookey, Cripps, "Improved surface bounds based on derivatives", CAGD 14(8):787-791, 1997 [snippet].
- Piyavskii 1972 (USSR Comp. Math. Math. Phys. 12(4):57-67) and Shubert 1972 (SIAM J. Numer. Anal. 9(3):379-388), "A sequential method seeking the global maximum of a function" [snippet, citations confirmed]. The algorithm samples where the Lipschitz upper bound f(x_i) + L|x - x_i| is largest. This is the adaptive version of "max|g| <= max sample + L h/2".
- Hansen, Jaumard, "Global optimization of univariate Lipschitz functions I: survey", Math. Prog. 55, 1992 [snippet]. Extensions to Lipschitz derivatives exist (Gokcesu, arXiv 2108.10859) [snippet].
- Sergeyev, Kvasov (diagonal approach, Lipschitz derivative methods), Lyness, Novak (information-based complexity of optimisation) [unread/memory].

### Hodograph and derivative bounds
- Floater, "Derivatives of rational Bezier curves", CAGD 9(3):161-174, 1992 [snippet: abstract]. Bound scales with (max weight / min weight) times control point differences.
- Saito, Wang, Sederberg, "Hodographs and normals of rational curves and surfaces", CAGD 12(4):417-430, 1995 [snippet: abstract]. Scaled hodograph of a degree m x n rational surface is a polynomial of degree (2m-2) x 2n; gives direction, not magnitude.
- Wang, Sederberg, Saito, "Partial derivatives of rational Bezier surfaces", CAGD 14(4):377-381, 1997 [snippet].
- Mao Shi, "Counterexamples to a conjecture on first derivative bounds of rational Bezier curves", arXiv 2510.17300, v2 Mar 2026 [read: abstract]. A degree-7 counterexample shows the Li et al. 2013 conjectured bound is false; proposes a bound by degree elevation that is valid at every elevation and converges to the true sup as elevation grows. Lesson: published closed-form rational bounds can be wrong; elevation-based bounds are safe but loose at low elevation.
- A 2020 MDPI Mathematics paper improving derivative bounds for NURBS surfaces cites Filip's tessellation algorithm (8(8):1382) [snippet; page returned 403].

### Composition
- DeRose, Goldman, Hagen, Mann, "Functional composition algorithms via blossoming", ACM TOG 12(2):113-135, 1993 [snippet: partial abstract]. Computes H = F o G in Bezier or B-spline form; product algorithm more efficient, blossom algorithm more intuitive. DeRose 1988 "Composing Bezier simplexes" TOG 7(3) is the precursor.

### Curve-on-surface, Hausdorff distance
- Chen, Ma, Xu, Paul, "Computing the Hausdorff distance between two B-spline curves", CAD 42(12):1197-1206, 2010 [snippet]. Geometric pruning with upper/lower bounds. (The brief's "Chen et al. 2010" is this paper; Elber is not an author.)
- Elber, Grandine, "Hausdorff and minimal distances between parametric freeforms in R2 and R3", GMP 2008, LNCS pp. 191-204 [snippet]. Exact distance by solving polynomial systems; a citing paper reports cost "at least O(n^2) for two curves composed of O(n) coefficients".
- Kim, Oh, Yoon, Kim, Elber, "Efficient Hausdorff distance computation for freeform geometric models in close proximity", CAD 45(2):270-276, 2013; Barton, Hanniel, Elber, Kim, CAGD 27(8), 2010 [snippet]. "Close proximity" is exactly the edge-gap case and is called out as the hard case.

### Interval and Taylor models in CAD
- Hu, Patrikalakis, Ye, "Robust interval solid modelling" Part I (CAD 28(10):807-817) and Part II (819-830), 1996 [snippet: abstract]. Interval B-rep in rounded interval arithmetic.
- Neumaier, "Taylor forms - use and limits", Reliable Computing 9:43-79, 2003 [snippet: abstract]. Emphasis on overestimation and the wrapping effect.
- Berz, Makino Taylor models [unread/memory].
- Sherman, Michel, Carbin, "Sound and robust solid modeling via exact real arithmetic and continuity", ICFP 2019 [snippet: abstract]. Sound Hausdorff distance via exact reals; no performance numbers found.
- Shewchuk 1997 adaptive predicates [unread/memory; well known].

### Certified tessellation with the same kind of bound
- Lutterkort, Peters slefes (2001); Peters, Wu "On the optimality of piecewise linear max-norm enclosures based on slefes" (cubics within 7% of optimal) [snippet].
- Yeo, Bhandare, Peters, "Efficient pixel-accurate rendering of (animated) curved surfaces", 2012 [snippet]. Tessellation factor guaranteeing half-pixel error; reports the AABB / min-max bound is "looser by an order of magnitude" than slefe width, and loose bounds cause over-tessellation. Directly relevant: the tightness of the bound sets the sample count.
- Patent US5377320A (trimmed NURBS rendering) cites Filip et al. [snippet].

## Q2. Applicability by geometry type

"Cheap bound" means from stored data with no search.

| Geometry | Cheap derivative / difference bound? | Tightness | What breaks it |
|---|---|---|---|
| Polynomial B-spline curve on polynomial B-spline surface | Yes. Hodograph nets; surface partials over the pcurve's uv box; or compose to Bezier and difference nets | Bounding C' and (S o p)' separately: very loose (see Q3). Difference net of g: tight, O(h^2) under subdivision | Degree blow-up (18 for bicubic with cubic pcurve), many knots, near-tangent pcurve |
| Rational (NURBS, weighted) | Yes, via Floater or Saito-Wang-Sederberg, or by homogeneous composition | Loose by factor ~ max w / min w, worse with w spread | Weight ratios, published-bound errors (arXiv 2510.17300); need to divide by weight after composing |
| Plane, cylinder, cone, sphere, torus in angle parameters | Yes. Closed-form derivatives; sin/cos bounded by 1 and by radius; interval evaluation of trig | Tight for plane/cylinder; sphere and cone fine away from apex/pole | Sphere poles and cone apex (degenerate partials, pcurve speed unbounded in u); seams. Alternatively these are exact rational quadrics, so the NURBS route applies |
| Surface of revolution, extrusion | Yes, if the profile curve is a spline (product of profile bound and rotation/translation bound) | Moderate | Axis crossing the profile (pole), rational profile |
| Offset surface | Only via interval enclosure of the normal: |S_u x S_v| in the denominator | Loose; blows up | Near-zero |S_u x S_v| and curvature radius near the offset distance |
| Swept / blend / procedural surfaces | No control net. Need interval evaluation of the defining procedure | Loose or none | solidcraft PR #18166 leaves these OPEN ("a trim on a procedural blend has no exact net") [read] |
| Procedural intersection curves | No closed form; need enclosure of the marching approximant plus proof it is near the true intersection | None today | solidcraft issue #13630: "neither proves the B-spline between its points"; planned fix needs Krawczyk root existence [read] |
| Imported STEP entities | Depends on entity: B-spline and analytic yes; trimmed, offset, swept, surface_curve with no pcurve, not same-parameter | Mixed | Missing pcurves, different parameterisations (parametric gap overestimates the geometric gap; solidcraft builds a reparametrisation psi), out-of-domain pcurves (solidcraft gen0642: trim runs at u=71.94 on a [-0.05,1.05] domain, extrapolated surface 5.6e13 away) [read] |

## Q3. How it would work, in plain English

1. Split the edge at every knot of C, of p, and every surface knot line the pcurve crosses. Inside each piece every function is one polynomial (or one rational).
2. On each piece, get a bound on how fast the gap can change. Two ways:
   - Naive (the brief's version): bound |C'| from C's hodograph net, bound |(S o p)'| <= |S_u||u'| + |S_v||v'| using the surface's differenced nets over the uv box. Then L = bound(|C'|) + bound(|(S o p)'|).
   - Difference form: build g itself as one polynomial (compose S with p, elevate C to the same degree, subtract control points). Then max|g| <= max |control point of g| directly, with no samples. Its hodograph gives a small L too.
3. Sample g at points h apart; take the largest sample; add the slack L h/2 (first derivative) or M2 h^2/8 (second derivative, Filip style).
4. If the bound fits under the tolerance margin, the piece passes. If the largest sample is already over, it fails. If unsure, halve the piece and repeat (Piyavskii-Shubert style adaptive refinement), up to a budget.
5. Arithmetic: compute in doubles with a forward error term; only near the threshold redo in double-double or intervals. OCCT now evaluates the residual in double-double already; solidcraft uses outward-rounded intervals throughout.

Sample counts. With first-derivative slack L h/2 <= m: n = L H / (2m) samples on a span of parameter length H. With second-derivative slack: n = H sqrt(M2/(8m)).

The trap: the naive L does not bound g', it bounds |C'| + |(S o p)'| ~ 2 x edge speed, because the triangle inequality throws away the cancellation. For a 100 mm edge over t in [0,1], L ~ 200 mm. With m = 1e-6 mm (Parasolid-scale resolution) that is 1e8 samples; with m = 1e-3 mm, 1e5 samples. That is unusable. The second-derivative version with M2 ~ 100 mm and m = 1e-3 needs ~ 110 samples, better but still bounding the curve's own curvature, not the gap's. Only the difference form gets a slack proportional to the gap's own wiggle, which for a good edge is near the gap size itself, so tens of samples or a few subdivisions suffice.

Where it gets expensive:
- Composition degree: bicubic S with cubic p gives degree 18 per knot span (3x3 + 3x3), times the number of surface patches the pcurve crosses. Rational surfaces double the degree (numerator and weight).
- Cancellation: g is ~1e-7 relative to coordinates ~1e2, so 9 of 16 digits cancel. Bounds need compensated or interval arithmetic everywhere.
- Slack that never shrinks: solidcraft hit a bound stalled at 0.098 for a true gap of 2.2e-4 because the uv box overshoots the net end; the overshoot term does not shrink with bisection [read].
- Wall time: solidcraft reports one case at 153 s for the full bound and 43 s with an early-exit cutoff [read].
- Loose cases remain: solidcraft lists edges where the sample says 3.6e-3 and the bound says 0.229 [read].

## Q4. Why has no kernel done it?

### Evidence (each read in this session)
1. Cost, OCCT. Issue #1372 (2026): `GeomLib_CheckCurveOnSurface::Perform` was 17.9% of a whole boolean operation by callgrind; ~620 evaluations per sub-interval to measure a gap of ~1e-16. Fix #1373 skips search when samples are below `Precision::SquareConfusion()`. Even the heuristic search is too expensive; a certified bound is more work.
2. Cost, OCCT 2016. Commit 0028030 "GeomLib_CheckCurveOnSurface takes too much time for Bspline curves with big number of knots": above a knot count it stops working per knot interval and uses fixed sample points.
3. Cost, OCCT 2021. Commit 0032448 made the exact method an opt-in `-exact` flag for `checkshape`; "Default mode is calculating in finite number of points."
4. Origin of the OCCT checker. Commit 0025597 (Jan 2015) "Invalid curve on surface in the result of General Fuse operation" moved the max-distance tool into BRepLib_CheckCurveOnSurface and applied it to pcurves built during booleans. Then 0025819, 0026323 "Tolerance computing unification", 0026918 (periodic curves failed). This is the post-2015 history the brief mentions.
5. Cost, SMLib. `SmGap.cpp` BuildSampleSet: "TODO: Gap Functions should be represented by functions and not by sampling. Max and Min gaps should be found by calling the Global solver ... That code is not written yet and would probably be much more expensive to run since we expect the Gap Function size to commonly have one maximum within each and every BSpline knot interval." `SmEdgeuse.cpp` GetMaxEdgeFaceGap: "TODO: extend SmSurface::GlobalCurveSolve to find the max gap between a curve and a surface"; "after a 1st expensive call - it's cheap to get" (cached).
6. Looseness and coverage. solidcraft PR #18166 and issue #13630: bounds are loose on some edges (60x), stall on non-shrinking slack, cannot cover procedural blends or traced intersection curves, and cost tens of seconds on one case.
7. Tolerances are already a model. Jackson 1995 "Boundary representation modelling with local tolerances" (Parasolid): edges carry their own tolerance that "reflect the extent to which surfaces don't meet" [snippet]. Parasolid PK_EDGE_check has a fault "SP-curves not within edge's tolerance of nominal geometry" [snippet]; how it measures is not documented in what I found.
8. Hausdorff papers call close proximity the hard case and report O(n^2) cost (Elber-Grandine, Kim et al. 2013) [snippet].

### Inference (not stated by any source)
- No customer demand: a missed gap between 23 samples rarely changes a boolean result, and the edge tolerance is often inflated afterwards anyway (solidcraft #18166 found SMLib-derived code doubling tolerance on each refit, so stored tolerances were already far above the true gap).
- Procedural geometry (offsets, blends, intersection curves) is the common case in real models, and the method gives nothing there without a separate enclosure for each surface type. A partial certificate is hard to sell.
- Implementation complexity: composition, rational handling, rounding control, pole handling and per-type enclosures are a large code surface for a check users do not see.
- Kernel tolerances (1e-7 to 1e-5) are loose compared with what samples miss on a smooth gap function; a smooth gap with one bump per knot span is mostly caught by mid-knot samples (SMLib's own comment).

## What I could not verify
- Full text of any paper; ACM, MDPI and ScienceDirect returned 403. Filip et al. abstract only via search summary.
- Parasolid and ACIS checking internals; the q-solid.com Parasolid docs were unreachable.
- OCCT forum and dev mailing list threads from 2015; I used the commit log instead.
- Any peer-reviewed paper that certified B-rep edge gaps and reported cost. None found.
- Sergeyev/Kvasov, Lyness, Novak, Berz/Makino specifics.
- Whether solidcraft PR #18166 merged.
