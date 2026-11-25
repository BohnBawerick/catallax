Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Edge-gap checker menu under contract T (measure, store gap x 1.1)

Date 2026-11-20. Read-only research. Builds on grill pack Parts C-E, edge-bound-research.md and
t-smarter-research.md; it does not repeat their findings.

Labels: [read] = I read the source, page or code this session. [snippet] = search summary only.
[unread/memory] = from memory, not checked. [exp] = my synthetic experiment (section 3).
[inference] = my reasoning, not stated by any source.

## 0. What changes when proof is dropped

Under T the checker has no pass/fail threshold except the tolerance cap. Its job is to return
the largest gap it can find, and a location. Two things follow [inference]:

- A miss does not cause a wrong pass/fail on its own. It causes a stored tolerance that is too
  small, which later shows up as a boolean or meshing fault at that edge. So the error that
  matters is *under-measurement*, and the only way to under-measure is to never look near the
  true peak.
- f64 is enough for the measurement almost everywhere. At 100 mm coordinates, f64 rounding in
  g = C(t) - S(p(t)) is about 1e-13 to 1e-12 mm. Kernel tolerances are 1e-7 to 1e-3 mm. Higher
  precision is needed only when the gap itself is near 1e-12 mm, and then the right answer is a
  noise floor, not a better number (section 2.10).

## 1. Menu summary

| # | Method | Handles well | Misses when | Cost (evals per edge, f64) | Used by |
|---|---|---|---|---|---|
| 1 | Fixed uniform samples | Lines, circles, single-span curves | Peak between samples; any fit whose spans do not line up with the grid. [exp] 23 points found 16-84% of the true gap | 20-23 | OCCT default (23) [read], SMLib (20) [read], Catallax today (9) [read] |
| 2 | Split at knots, seams, poles; k per piece | All piecewise geometry; fitted edges (gap peaks once per span) | Feature narrower than piece/k not tied to a knot. [exp] k=8 found 94-100% on clean fits, 29-36% on narrow defects | (pieces) x k, ~100-300 | OCCT exact-mode front end [read], Catallax native trims (8/span) [read], SMLib mid-knot samples [read], truck 6 x degree [read in Part D] |
| 3 | Samples plus local maximisation | Everything smooth per piece; pins the reported value | Peak with no local max visible in the samples. [exp] split + 54 samples/piece + Brent on every local max found 100% down to bump width 1e-4 of edge, missed 3e-5 | ~1,000-3,300 at 3 x degree samples/piece | OCCT master GeomLib_CheckCurveOnSurface (PSO + bounded Brent per bracket) [read] |
| 4 | Bend-driven spacing; error-driven halving | Long gentle edges; near-limit regions | Bend is not gap: a flat region with a fit error. Halving needs a point to see something first | Data dependent | OCCT GCPnts, Approx_SameParameter, BRepMesh [read in Part C]; Catallax import trace [read in Part C]; GAIL funmin (guaranteed under "non-spiky" cone) [read: abstract] |
| 5 | Chebyshev / spectral approximation of g | Polynomial pieces (exact at degree+1 points); analytic surfaces; procedural surfaces you can only evaluate | Narrow feature that aliases below the chop level. [exp] missed bump width 1e-4 on one case that #3 caught | ~400-1,500 in [exp] (1.3-2.4x fewer than #3) | Chebfun [read]. No CAD checker found. OCCT AdvApprox uses Jacobi-series tail as max-error estimate for fitting, not checking [read] |
| 6 | Projection (geometric) gap vs same-parameter gap; Hausdorff | Imports with mismatched speed; Parasolid-style tubes | Projection jumps to the wrong branch on near-closed or folded geometry; it hides parameter drift that same-parameter code downstream will see | 3-10 x #1 (each sample is a local projection) | OCCT ValidateEdge falls back to two-way local projection when not same-parameter [read]; Approx_SameParameter reparametrises [read]; Chen et al. 2010, Kim et al. 2013 [snippet] |
| 7 | Builder hints | Own builds (fits, traces) | The hint is wrong or correlated with the checker's own points | ~0 extra | OCCT Approx reports MaxError per span [read]; SMLib TestSpanAccuracy [read in Part C] |
| 8 | Reuse mesh points / SDF | Free extra seeds where the mesh is dense | Mesh is coarse where flat; display deflection (1e-2 mm) is 1e4 x coarser than edge tolerance; SDF cells far coarser still | 0 if mesh exists | OCCT BRepMesh splitSegment tests curve-on-surface agreement while meshing [read in Part C] |
| 9 | Random / low-discrepancy / extreme-value statistics | Giving points the builder did not aim at | Smooth 1D functions: worse than deterministic grids. [exp] 64 golden-ratio points found 26-98% | as chosen | de Haan 1981 JASA 76:467 [snippet]; Chebfun explicitly avoids randomness [read] |
| 10 | Floating-point filter | Gaps near the noise floor | Never needed for the T measurement except below ~1e-11 mm | x1 normally, x10-20 at noise floor | OCCT master double-double residual (Oct 2026) [read]; Shewchuk/CGAL pattern [unread/memory] |
| 11 | Closed form for elementary pairs (added) | Line on plane, circle on plane/cylinder, conic pairs | Nothing (exact) | O(1) | OCCT `ElementaryMaximum` [read] |
| 12 | Lipschitz search with estimated constant (added) | Gives a "how much could I be missing" estimate | Estimate of L is wrong on spikes | 1.5-3 x #3 | Strongin 1978, Sergeyev 1995 local tuning [snippet] |
| 13 | Coefficient/interval bound (the "prove" button) | Polynomial and rational pairs (coefficients); everything with a formula (intervals) | Loose, slow; procedural blends | 10-1000 x | solidcraft PR #18166 [read in edge-bound]; Catallax balls [read in Part C] |

## 2. Each method in detail

### 2.1 Fixed uniform samples
How: evaluate g at n equally spaced parameters, take the max.
Evidence: OCCT `BRepLib_ValidateEdge` `myControlPointsNumber(22)`, i = 0..22 [read, local src].
SMLib `lSampleCnt = 20` [read]. OCCT header calls it "faster, but possible not correct" [read in Part D].
Handles: lines, circles, single-span curves on analytic faces.
Misses: fitted geometry. A fit's error peaks once per knot span (SMLib comment: "built to
maximize their gap sizes at the mid-knot points" [read]). With more spans than samples, most
peaks are never seen. Equal parameter spacing also bunches where the curve is slow.
[exp]: ../../../research/cp2-grill/on 4 clean fitted edges uniform 23 returned 0.49 to 0.84 of the true max. With a narrow
defect added, 0.16 to 0.28.
Cost: n evaluations. Complexity: trivial.

### 2.2 Split at all breaks, k samples per piece
How: cut the parameter range at every knot of C, every knot of p, every parameter where p
crosses a surface knot line, every seam crossing, every pole approach, every C0 join of a
procedural curve (trace chart points). Inside each piece every input is one smooth function.
Sample each piece.
Evidence: OCCT `FillSubIntervals` merges 3D and 2D knots and surface intervals; sample count per
piece is `3 * max(nbParticles, composed degree)` where composed degree = pcurve degree x (surface
u degree + v degree) [read, local src lines 1739-1835]. That rule is principled: a polynomial of
degree d has at most d-1 interior extrema, so 3d samples give a few samples per extremum.
Handles: everything piecewise. Fitted edges become the easy case because their peaks are tied
to spans.
Misses: a feature narrower than piece length / k that is not tied to a break. Sources of such
features [inference]: a pcurve passing near a surface pole; a near-tangent intersection trace;
an imported pcurve with a local kink that has no knot; an offset surface near its curvature
radius.
[exp]: ../../../research/cp2-grill/k = 8 per piece found 0.94-1.00 of the true max on clean fits, 0.29-0.36 with a narrow
defect.
Cost: pieces x k. A 324-station fitted strip with k = 8 is ~2,600 evals [read in Part C]; in
f64 that is a millisecond or two [inference].
Complexity: low. Finding surface knot-line crossings needs a 1D root find on u(t) - knot.

### 2.3 Samples plus local maximisation
How: after the samples, bracket every local maximum of the sampled |g| (each sample that is at
least as large as both neighbours) and run a bounded 1D maximiser (Brent, golden section) on the
true g inside the bracket.
Evidence: OCCT master `MinComputing`: samples the piece, collects every local-extremum bracket,
runs PSO (3+ particles), then `MathOpt::Brent` with `XTolerance = PConfusion / bracket length`
on each bracket [read, local src lines 1446-1611]. Oct 2026 history: fix to "preserve narrow
extrema" [read in Part D]; issue #1372: 17.9% of a boolean, ~620 evals per sub-interval [read
in edge-bound].
Key finding [exp]: polishing *every* local max, not just the best one, is what catches narrow
features. A narrow bump on top of a sloped smooth gap makes a small local max in the sample
sequence even when the bump's own peak is far between samples. Dense 54 samples per piece
without polishing found 0.82-0.88 of a width-5e-4 bump; 17 samples per piece plus polishing
every local max found 1.00. Polishing every local max at 54/piece caught width 1e-4; nothing
tested caught width 3e-5.
Handles: all smooth geometry. The answer is a real local max, accurate to ~1e-10 in t.
Misses: a peak that leaves no local max in the samples (sits on a monotone stretch and is
narrower than the spacing); a piece with many shallow maxima only the best of which is polished.
PSO adds little in 1D once every bracket is polished [inference]; skip it for determinism.
Cost: ~10-20 evals per polished bracket. [exp] 1,000-3,300 evals per edge for the OCCT-like rule.
Complexity: low to medium. Brent is ~60 lines.

### 2.4 Bend-driven spacing and error-driven refinement
How: (a) set spacing h from curvature or second-derivative size so h^2 M / 8 is a share of the
tolerance; (b) halve intervals where the measured gap uses a large share of the current max.
Evidence: Part C (OCCT GCPnts, Approx_SameParameter IncreaseNbPoles, BRepMesh splitSegment,
Catallax import trace). New: Choi, Ding, Hickernell, Tong, "Local adaption for approximation and
minimization of univariate functions", J. Complexity 40:17-33, 2017 [read: abstract]. Their
`funmin` estimates the second derivative from divided differences of the samples and is
guaranteed for a cone of "non-spiky" functions; cost O(sqrt(||f''||_1/2 / eps)). This is the
closest published thing to a "guarantee under a stated smoothness assumption" without a proof
of the bound.
Handles: long gentle edges, regions close to the cap.
Misses: bend in C or S is not bend in g (two curves that bend together have a flat gap); refines
only where a sample already saw something.
Cost: data dependent; worst case set by a cap. Complexity: medium (thresholds, caps, refusal).

### 2.5 Chebyshev / spectral approximation of g
How: on each piece, evaluate g at n+1 Chebyshev points (cos(j pi / n)), get Chebyshev
coefficients by FFT, double n (17, 33, 65, 129; values are reused because the grids nest) until
the coefficients fall to a plateau, chop there. The tail size estimates the interpolation error.
The max is then found from the interpolant (roots of its derivative via a colleague matrix, or
dense evaluation of the cheap polynomial), and polished on the true g.
Evidence:
- Chebfun construction: grids of 17, 33, 65, ... points, stop when coefficients reach a plateau
  [read, Aurentz and Trefethen, "Chopping a Chebyshev series", arXiv 1512.01803, ACM TOMS 2017].
  standardChop: upper envelope of |coeffs|, search for a "sufficiently long, sufficiently flat
  plateau" no higher than tol^(2/3), then chop at the elbow [read].
- Aliasing guard: `sampleTest` evaluates at two arbitrary off-grid points. T_128 equals 1 at every
  point of the 17, 33 and 65 grids, so a coefficient test alone would call it a constant. The
  authors: "in principle a function might arise that outwits sampleTest, but we have never seen
  this happen" [read]. Footnote: randomised tests would defeat any fixed adversary, but "Chebfun
  avoids randomization, because it is too disturbing to users" [read].
- Narrow spikes: Chebfun Guide section 8 minSamples. Default 17 (9 before v5); a spike near
  x = 0.47 between grid points 0.383 and 0.556 is missed; fix by raising minSamples. The
  "Spike integral" example: default settings miss the narrowest spike in splitting mode [snippet].
- Root finding: Battles and Trefethen 2004 (SISC) and Boyd 2002 (SINUM) colleague-matrix roots
  with recursive subdivision for high degree [unread/memory].
- CAD use: I found no CAD checker that uses it [snippet: search turned up none]. The nearest
  production use is OCCT `AdvApprox_ApproxAFunction`, which fits functions in a Jacobi
  (Legendre-family) basis and estimates max error as sum over dropped coefficients of
  |coeff| x max|basis| (`PLib_JacobiPolynomial::MaxError`) [read, local src]. It is used for
  Approx_CurveOnSurface, Approx_SameParameter, sweeps and curve conversion [read: grep of
  callers]. So spectral tail error is trusted in a shipping kernel for *building* curves.
Why it fits edge gaps [inference]:
- After splitting at all breaks, g on a polynomial piece *is* a polynomial of known degree
  (bicubic surface + cubic pcurve + cubic C: degree 18). Chebyshev interpolation at degree+1
  points reproduces it to rounding. No aliasing is possible. The max then comes from a degree-36
  polynomial |g|^2. This is a near-exact max for polynomial pieces at ~19 evals.
- Rational pieces and angle-described analytic faces give analytic g; Chebyshev converges
  geometrically, so 17-33 points per piece usually suffice.
- Procedural surfaces (offsets, sweeps, blends, traced intersections) only need an evaluator.
  Convergence speed then depends on how close complex singularities are (offset near curvature
  radius, near-tangent traces), and the coefficient decay shows it: slow decay is a built-in
  "this piece is hard or has a hidden kink, split it" signal.
- Chop level must be absolute, tied to the tolerance (say 1e-3 of the measured gap or of the
  cap), not relative to machine epsilon. That keeps degrees low [inference].
Robustness [exp]: on clean fits it matched the truth at 1.3-2.4x fewer evals than the OCCT-like
rule. It caught narrow bumps of width 2e-3 and 5e-4. It missed width 1e-4 on 3 of 4 runs: once
the smooth part was resolved, the bump's coefficients sat under the chop level. Tightening the
chop to 1e-5 caught it on the cylinder but not the B-spline case.
Misses when: a narrow feature aliases below the chop level; non-smooth points not split out
(algebraic decay, n grows to the cap); noise floor (gap ~1e-12, coefficients are noise).
Cost: n per piece, nested doubling; FFT and colleague-matrix work are negligible next to evals.
Complexity: medium. DCT, chop rule, off-grid test, root finding or dense interpolant eval.
Verdict [inference]: excellent *estimator* and diagnostic; on its own slightly less robust than
#3 with local-max polishing, because chopping is a decision to stop looking.

### 2.6 Geometric (projection) gap vs same-parameter gap
Definitions: same-parameter gap max_t |C(t) - S(p(t))|. Geometric gap: one-sided Hausdorff,
max over points of one copy of the distance to the other copy (point to surface, or curve to
curve), taken both ways.
Evidence: OCCT ValidateEdge switches to projection when the edge is not flagged same-parameter
or the parameter ranges differ by more than PConfusion; it projects each sample both ways with
`Extrema_LocateExtPC` (a local search seeded at the proportional parameter) [read, local src
lines 118-222]. SMLib CalcMaxEdgeFaceGap measures C to the *surface*, not through the pcurve
[read in edge-bound]. ACIS TEDGE: max distance between equiparametric positions [read in
t-smarter]. Parasolid tolerant edge: a geometric tube [read in t-smarter]. Approx_SameParameter
fits a reparametrisation from projected correspondences, then measures at 2 x NbSamples, x1.05,
and accepts within aMult = 250 of the best achieved [read, local src lines 193, 367-466].
Hausdorff literature: Chen, Ma, Xu, Paul 2010 (B-spline curves, pruning); Elber and Grandine
2008 (exact, polynomial systems, ~O(n^2)); Kim et al. 2013 (close proximity is the hard case)
[snippet, from edge-bound].
When each is right [inference]:
- Same-parameter gap is right when downstream code evaluates the copies at a shared parameter
  (boolean edge splitting, meshing along the edge, OCCT's SameParameter contract, Catallax's
  constructed shared parameter). Then the stored tolerance must cover the parametric gap.
- Geometric gap is right for a tube-style tolerance (Parasolid) and for judging whether an
  import's copies describe the same curve.
- For imports with mismatched speeds: measure the geometric gap; if parametric >> geometric, the
  problem is parameterisation, not geometry. Reparametrise (Approx_SameParameter style) and
  re-measure same-parameter. Store the post-repair parametric gap. Do not inflate the tolerance
  to cover speed mismatch (SMLib rebuilds the trim instead of raising tolerance [read in
  t-smarter]).
Failure modes: local projection lands on the wrong branch on near-closed curves, at seams, or
when the gap exceeds the local curvature radius; projection near a pole is ill-conditioned.
Geometric gap <= parametric gap always, so using it alone under-measures for same-parameter
consumers.
Cost: each sample is a Newton projection, 3-10 evals [inference]. Complexity: medium (seeding,
seams, both directions).

### 2.7 Builder-supplied hints
How: the fitter, tracer or importer reports its own measured residual and where it peaked. The
checker (1) seeds a polish there, (2) treats the hint as a lower bound on the result, (3) flags
a bug if it measures less than the builder claimed at the same t.
Evidence: OCCT approximators expose MaxError per span (AdvApprox, Approx_CurveOnSurface) [read];
SMLib TestSpanAccuracy measures 5 interior points per traced span [read in Part C]; Catallax's
fitter tests at the interior eighths and logs the largest sampled residual (2.96e-8 mm on R3.99)
[read in Part C].
Misses: when the builder's points and the checker's points coincide, both see the same thing
(Part C 1.3). Hints never replace independent points.
Cost: ~0. Complexity: low, but it adds a data contract between every builder and the checker.

### 2.8 Reuse of tessellation or SDF
How: wherever the mesher already evaluated both copies at mesh vertices, use those values as
free extra samples or seeds.
Evidence: OCCT BRepMesh_CurveTessellator `splitSegment` compares the surface point at the UV
midpoint with the chord, recursively to depth 10 [read in Part C]. No kernel found that feeds
mesh points into its tolerance check.
Limits [inference]: mesh points are placed by display deflection (typically 1e-2 to 1e-1 mm), so
they are sparse where the shape is flat, which is where fit drift hides. Meshing may not exist at
check time (import, export). The SDF/field print grid resolution is far coarser than 1e-7 mm, so
it can only screen for gross errors (gap > cell size).
Cost: 0. Complexity: low, but couples checker results to mesh settings, which hurts determinism.
Verdict: optional seeds only.

### 2.9 Statistical methods
How: random or low-discrepancy samples; extreme-value statistics on the top order statistics to
estimate how far the true max lies above the best sample.
Evidence: de Haan, "Estimation of the minimum of a function using order statistics", JASA
76(374):467-469, 1981 [snippet, citation confirmed]; a confidence interval for the optimum
assuming the function behaves like |t - t*|^alpha near it, exponent d/2 for smooth functions
[snippet via arXiv 2402.11168]. Zhigljavsky's global random search literature [snippet].
Chebfun rejects randomness for user trust [read].
Verdict [inference + exp]: for a 1D smooth function, random points are strictly worse than a
deterministic grid plus polishing (64 golden-ratio points found 0.26-0.98 of the true max). The
extreme-value interval assumes the max looks like a smooth peak, which is exactly the assumption
that fails on narrow features. Two useful pieces survive: (1) a deterministic low-discrepancy
*offset* gives points the builder did not aim at; (2) "distance from best sample to estimated
max" is a cheap confidence indicator, but Chebyshev tails and Brent already give it.

### 2.10 Floating-point filter
How: run in f64, attach a forward error bound, escalate to double-double or 128-bit only where
the answer is not resolved.
Evidence: OCCT master `CompensatedReal` and `BSplineResidual` evaluate C(t) - S(p(t)) "without
rounding intermediate points" because the difference is small compared with either point [read
in edge-bound; class present in local src lines 92-320]. OCCT #1372: Perform was 17.9% of a
boolean while measuring gaps of ~1e-16; fix #1373 skips the search below
Precision::SquareConfusion [read in edge-bound].
Under T [inference]: there is no threshold to resolve. The only place f64 fails is when the gap
is below its own rounding error (~1e-13 to 1e-12 mm at 100 mm scale). Do not escalate there;
record max(measured, noise floor) where the floor comes from the f64 error estimate. Escalate
only when the measured gap is within a rounding error of the cap, or when the user presses the
128-bit "exact at points" button. This also removes the OCCT #1372 cost trap: searching noise.

### 2.11 to 2.13 (added)
- Closed forms for elementary pairs: OCCT `ElementaryMaximum` solves line/conic cases with no
  search [read]. Catallax already has exact closed forms for analytic pairs [read in t-smarter].
- Lipschitz search with estimated constant (Strongin 1978; Sergeyev, SIAM J. Optim. 5(4):858-870,
  1995, local tuning) [snippet]: samples where an estimated upper envelope is highest. Not a
  proof (L is estimated), but a principled refinement order. Similar role to 2.4b.
- Coefficient or interval bounds: the optional "prove" button; see edge-bound-research.md.

## 3. Experiment [exp]

Script: `exp/gapcheck.py`
and `variants.py` in the same folder. Python, numpy 2.5 / scipy 1.18.

Setup: bicubic B-spline surface (7x7 net, 100 mm square, random heights) or cylinder R = 10 mm
in angle parameters; cubic pcurve with 5 spans; C = least-squares cubic fit of S(p(t)) with 5 to
30 interior knots, ends pinned. So the gap is a real fit residual. "Narrow bump" adds a Gaussian
defect of 3 x the residual at t = 0.4137 (not on a knot), width w in edge-parameter units.
Truth: 400,001 uniform samples plus polish. Numbers are found max / true max (evals).

| Case | uniform23 | golden64 | split+8 | split+8+polish best | OCCT-like (split, 54/piece, polish every local max) | split+Chebyshev+polish |
|---|---|---|---|---|---|---|
| B-spline, 5 knots | 0.66 (23) | 0.98 | 1.00 (161) | 1.00 (566) | 1.00 (1621) | 1.00 (747) |
| B-spline, 20 knots | 0.49 | 0.78 | 1.00 (281) | 1.00 (921) | 1.00 (2979) | 1.00 (1257) |
| Cylinder, 6 knots | 0.67 | 0.98 | 0.94 (89) | 1.00 (298) | 1.00 (965) | 1.00 (413) |
| Cylinder, 30 knots | 0.84 | 0.85 | 1.00 (281) | 1.00 (889) | 1.00 (3190) | 1.00 (1256) |
| B-spline 20 + bump w = 5e-4 | 0.16 | 0.26 | 0.33 | 1.00 | 1.00 | 1.00 |
| Cylinder 30 + bump w = 5e-4 | 0.28 | 0.28 | 0.33 | 0.34 | 1.00 | 1.00 |
| B-spline 20 + bump w = 1e-4 | 0.16 | 0.26 | 0.33 | 0.33 | 1.00 | 0.33 |
| Cylinder 30 + bump w = 1e-4 | 0.28 | 0.28 | 0.33 | 0.34 | 1.00 | 0.34 (1.00 with chop 1e-5) |
| Either + bump w = 3e-5 | 0.16-0.28 | | 0.33 | | 0.33 | 0.33 |

From variants.py: split + 17 uniform per piece + polish every local max found width 5e-4 at
~1,760-2,020 evals; split + 54 per piece with no polishing missed it (0.82-0.88).

Lessons [exp + inference]:
1. Splitting at all breaks is what makes clean fitted edges easy (split+8 already 0.94-1.00).
2. Polishing every local max, not only the best sample, extends reach to features ~5x narrower
   than the sample spacing.
3. Chebyshev is the cheapest at full accuracy on smooth pieces but its chop can stop before a
   narrow feature shows.
4. Every method has a width below which it misses (D1 in Part D). The detectable width scales
   with sample spacing.
Caveats: synthetic, one random surface, Python; eval counts, not time; the bump is an artificial
defect; no rational surfaces, poles, seams or procedural surfaces were tested.

## 4. Geometry kinds

| Kind | Front-end split points | Best engine | Notes [inference unless marked] |
|---|---|---|---|
| Line/circle on plane, cylinder, cone, sphere | none | closed form (#11) | OCCT ElementaryMaximum [read] |
| Spline curve on plane | C knots | polynomial pieces: Chebyshev at degree+1 is exact; or samples+polish | |
| Spline/rational on B-spline surface | C, p, surface knot-line crossings | samples (3 x composed degree) + polish every local max; Chebyshev tail as check | composed degree 18 for bicubic+cubic [read: OCCT rule] |
| Rational (weights) | as above | same; g is rational, so no exact degree; adaptive n | weight spread slows convergence |
| Torus, sphere near poles, cone near apex | seam crossings; pole approach (split where pcurve speed blows up) | samples+polish; geometric gap near poles | pcurve parameter speed unbounded near poles |
| Offsets | inherit base breaks; split where offset distance nears curvature radius | adaptive (Chebyshev decay flags trouble) | no control net |
| Sweeps, blends | profile and rail knots; spine breaks | evaluator-only methods (#3, #5) | solidcraft cannot bound these [read in edge-bound] |
| Traced intersection curves | every chart/trace point (C0 joins) | #3 with all chart points as breaks; SMLib tightens near tangency [read in Part C] | near-tangent traces are where narrow features live |
| Imported STEP, mismatched speeds | knots of all copies; parameter range ends | geometric gap first, reparametrise, then same-parameter (#6) | ISO 10303-42 promises no parameter tie [read in Part E] |

## 5. Complete setups

### Setup A: split, sample, polish everything (OCCT master pattern, cleaned up)
1. Elementary pairs: closed form.
2. Imports with mismatched parameters: geometric gap both ways; reparametrise; continue with
   same-parameter.
3. Split at every break (section 4).
4. Per piece: n = max(17, 3 x composed degree) samples, offset from the builder's check fractions
   (Catallax fitter uses eighths) so they are independent; add builder hint location.
5. Polish every local max of the samples with bounded Brent on the true g, f64. No PSO.
6. Noise floor from f64 error estimate; double-double only near the cap or on request.
7. Store max x 1.1 and the location. Hard per-edge eval cap with a named refusal (Part C).
Cost: ~1-3k f64 evals for a 20-30 knot edge [exp], ~0.3-1 ms [inference: 100-300 ns per
B-spline surface point, unverified]. Complexity: low.

### Setup B: split plus Chebyshev
Same front end as A. Per piece: polynomial pieces get degree+1 Chebyshev points (exact);
others adaptive 17/33/65/129 with an absolute chop tied to the tolerance and two off-grid
sample tests; max from all local maxima of the interpolant, each polished on the true g; slow
coefficient decay forces a split. Cost ~1.3-2.4x fewer evals than A [exp]; returns an error
estimate per piece. Complexity: medium. Weakness: the chop can stop before a narrow feature
shows [exp].

### Setup C: A plus estimated-envelope refinement
A, then a GAIL/Sergeyev-style pass: from divided differences estimate a local second derivative
per interval, inflate it (say x4), and add points wherever sample + inflated M h^2/8 exceeds the
current max by more than 10%. Stops with an estimated headroom that is reported, not proven.
Cost 1.5-3x A. Complexity: medium-high. Guarantee holds only for "non-spiky" gaps (Choi et al.
2017).

### Pick
A, with two pieces of B grafted on: use Chebyshev-spaced points (so the sample set nests and
doubles without waste) and use the coefficient tail as a trigger to *double* n, never to stop
below the floor. Reasons: f64 evals are cheap, so spending 2x evals over B buys the robustness
that the experiment showed (A caught a width-1e-4 defect B missed). Every step is plain code a
reviewer can follow. It matches what the only kernel that has iterated on this for a decade
(OCCT 2015-2026) converged to. Confidence: medium. The experiment is synthetic, and the cost per
eval in Rust is unmeasured.

## 6. Could not verify
- Full text of Battles-Trefethen 2004, Boyd 2002, de Haan 1981, Choi et al. 2017, Sergeyev 1995
  (abstracts or summaries only).
- Any CAD kernel that uses Chebyshev or spectral tails for *checking* (none found).
- Parasolid and ACIS check internals.
- Per-eval cost in Rust for Catallax evaluators; numbers above are estimates.
- Behaviour on rational surfaces, poles, seams, offsets and traced curves: not in the experiment.
