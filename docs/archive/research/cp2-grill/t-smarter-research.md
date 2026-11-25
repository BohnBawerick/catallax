Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Smarter ways to reach T's goals: research notes

Question: is there an approach, in a shipping kernel or in the literature, that meets T's goals
better than baseline T? Baseline T means float build, adaptive gap sampling, stored per-entity
tolerance, one check, and exact arithmetic only where it is cheaper.

Source labels: [read] means I read the source or the code this session. [unread/memory] means
it comes from prior knowledge and I did not verify it here.

Catallax fit notes come from `crates/catallax-topology/src/numerical.rs` (doc comment on
`NumericalBody`, struct defs). Each `NumericalEdge` keeps a 3D `curve`. Each `NumericalTrim`
keeps uv endpoints (linear) or a retained native uv curve. The check evaluates 9 equally spaced
locations per coedge, plus 8 subdivisions per knot span for native uv curves. The doc comment
says outright that these are "pointwise distance enclosures, not bounds between samples". So the
96.7M work units buy proof at 9 points and nothing in between.

---

## 1. One authoritative copy per edge (Parasolid exact/tolerant edges, ACIS EDGE/TEDGE)

**What it is.** An edge never stores two copies that can disagree. Parasolid has two cases.
- An *exact* edge has only a 3-space curve and no stored pcurve. When a surface/surface
  intersection produces the curve, it is a procedural "intersection curve": the two surfaces,
  the two end limits, and a "chart" of points. Parasolid can evaluate it "precisely at any
  point".
- A *tolerant* edge has **no** 3D curve. Each fin (coedge) carries an SP-curve: a 2D B-curve
  in the face's parameter space, so it lies on its surface by construction. The edge tolerance
  is a tube around one chosen SP-curve, and every SP-curve must stay inside the tube.

ACIS mirrors this. In an exact EDGE the edge geometry is primary and the coedge pcurve is
secondary. In a TEDGE the TCOEDGE pcurve is primary, the 3D curve is "lazy", and the edge
tolerance is "the maximum distance between any two equiparametric positions on any of its
tolerant coedges".

**Who.** Parasolid since about 1994 (Jackson 1995) and ACIS (R10 docs and later).

**Evidence.** Claimed, not measured: Jackson says tolerant modelling avoids "performing
numerical computations to inappropriately high accuracy, which would be slow and unreliable".
Parasolid docs: "enquiries on exact geometry are usually faster than enquiries on approximate
SP-curves". They also say that precision in the range 1e-6 to 1e-8 makes SP-curves with many
segments and a "drop in performance".

**Cost/risk.** For an exact edge you must derive a pcurve whenever an algorithm needs uv. You
can cache it, but the cached copy is never the authority. For a tolerant edge, every 3D query
lifts a pcurve through its surface. For a manifold edge the two lifted pcurves differ by up to
the tolerance, so you pick one as canonical.

**Fit to T and Catallax.** This is the best structural fit. Today a fitted edge has three
geometric copies (3D curve plus two uv trims), so there are two curve-vs-trim gaps to check
and two to export. If a fitted or imported edge drops the stored 3D curve, the tolerance check
becomes one comparison per edge: lifted pcurve A against lifted pcurve B, at shared parameters.
Neither copy can be off its own surface. Where a pcurve is an isoline, such as the rail of a
rolling-ball fillet on the fillet surface, that side is exact for free. For exact edges
(analytic pairs, line/circle on plane), keep only the 3D curve and treat pcurves as derived
caches that never need checking. That matches T's "exact where cheaper" rule.

STEP export still wants an EDGE_CURVE. Write a lifted SP-curve (or a fit of it) as the 3D curve
and attach the tolerance. Parasolid's STEP path has to do the same.

**Sources.**
- Jackson 1995, https://ftp.cs.wisc.edu/pub/users/prem/jackson-SM-95.pdf [read]
- Parasolid V35 docs, Session and local precision, ch. 16,
  http://www.q-solid.com/Parasolid_Docs_V35/chapters/fd_chap.017.html [read]
- Parasolid XT format reference: intersection curve chart/limits, null edge curve when
  tolerant (via search-result excerpts) [partly read]
- ACIS R10 Tolerant Modeling ch. 6 pages 0001/0002 (www-isl.ece.arizona.edu mirror) [read]

## 2. Continuous gap bound from samples plus a derivative bound, in float, with exact fallback (filter pattern)

**What it is.** Let g(t) = C(t) - S(p(t)) be the gap vector. If |g'(t)| <= L on a span of
width h, then the maximum of |g| on the span is at most max|g(samples)| + L*h/2. Bound L
cheaply from control nets. A B-spline's derivative lies in the hull of its differenced control
points. For S(p(t)), the chain rule gives
|dS/dt| <= |S_u|max*|p_u'|max + |S_v|max*|p_v'|max, using differenced surface nets restricted
to the uv box of the pcurve span. Choose the sample count per span so that L*h/2 fits in the
margin. That gives adaptive sampling with a real between-sample bound.

A tighter variant composes S with p into Bezier form by blossoming (DeRose et al. 1993) and
subtracts C after degree elevation. By the convex hull property, the largest control-point
difference bounds the gap. The degree is high, though: a bicubic composed with a cubic gives
degree 18.

Run all of this in f64 with a standard forward error bound. Call Catallax's 128-bit balls only
when the float result lands within its error band of the tolerance or cap. That is the
Shewchuk / CGAL filtered-predicate pattern applied to a measurement.

**Who.** Filip, Magedson, Markot 1986 ("Surface algorithms using bounds on derivatives",
CAGD) give second-derivative bounds for flatness and sample-count control [unread/memory].
DeRose, Goldman, Hagen, Mann 1993, "Functional composition algorithms via blossoming", ACM TOG
[unread/memory]. Shewchuk 1997 adaptive predicates, CGAL `Lazy_exact_nt` / filtered kernels,
and Cherchi et al. 2020/2022 indirect predicates for mesh booleans are the filter pattern
[unread/memory].

What shipping kernels do instead, for comparison:
- OCCT `BRepLib_ValidateEdge` defaults to 22 control points (23 evaluations), evenly spaced,
  with a 1.00001 multiplier on the found distance and early exit once the tolerance is
  exceeded. It has an optional "exact" mode, `GeomLib_CheckCurveOnSurface`. That mode splits at
  the knots of the 3D curve, the pcurve and the surface intervals, then runs `math_PSO` (3
  particles) plus Brent per sub-interval, in parallel. It is a better optimiser but gives no
  guarantee. [read]
- SMLib `SmEdgeuse::CalcMaxEdgeFaceGap` uses `SmCrvSrfGapFunction` with 20 evenly spaced
  samples. A TODO says to replace it with a global solve. [read]

No kernel I read gives a between-sample bound.

**Evidence.** The mathematics is standard. I found no measured kernel benchmark for this exact
use.

**Cost/risk.** Derivative bounds from control nets can be loose on high-curvature or rational
spans and force many samples there. Rational weights need homogeneous-space care (see arXiv
2510.17300 on rational Bezier derivative-bound pitfalls). The cost is bounded and runs in f64.

**Fit.** Strong. It does the adaptive sampling T already wants. It also improves on today's
guarantee: today's 9 points are proven, but nothing between them is. The float measurement plus
a derivative bound gives a bound everywhere, far cheaper than balls at 9 points. Balls become
the near-threshold fallback and the "prove this part" button.

**Sources.** OCCT `src/ModelingAlgorithms/TKTopAlgo/BRepLib/BRepLib_ValidateEdge.cxx` and
`src/ModelingData/TKGeomBase/GeomLib/GeomLib_CheckCurveOnSurface.cxx` [read]. SMLib
`source/SMLib/src/SmEdgeuse.cpp` (`CalcMaxEdgeFaceGap`, `lSampleCnt = 20`) [read]. Others as
labelled above.

## 3. Lazy, cached gap per entity with invalidation; persist the tolerance

**What it is.** Compute the max gap once, on first request, cache it on the entity, and clear
the cache when geometry changes. Export writes the stored number. Reopen trusts the stored
tolerance, with an optional re-check.

**Who.**
- SMLib: `GetMaxEdgeFaceGap(bForceCalc)` notes "MaxFaceEdgeGap is Cached - after a 1st
  expensive call - it's cheap to get". Mates share the cache. [read]
- ACIS: `TEDGE::set_update(logical)` marks that the tolerance needs recalculating, plus
  `api_update_tolerance`. [read via search excerpt of the ACIS reference]
- Parasolid XT stores the edge and vertex tolerance in the file. [read]

**Evidence.** Design claims only.

**Cost/risk.** Cache invalidation has to be correct. With Catallax's immutable `rebuild`
model, that is easy: a new body means a new cache.

**Fit.** This directly removes two of today's three checks (export and reopen). T already
plans it. The kernels confirm it is standard, so this is a confirmation rather than a new idea.

## 4. Repair before you grow; shrink when you can; cap growth

**What it is.** When a gap is too big, first try to remove it geometrically, and raise the
tolerance only if that fails. Lower tolerances again when the geometry allows.

- **SMLib** `SmBrep::ValidateAndUpdateTolerances` [read]:
  - If the edge-vs-surface gap is good but the edge-vs-uv-trim gap is bad, rebuild the uv trim
    instead of raising the tolerance.
  - Otherwise set edge tolerance = 2 x max gap, but only while the gap < bbox diagonal / 800.
    Larger gaps are left as errors ("symptoms of more serious Brep problems").
  - Vertex tolerance = 2 x max(vertex-edge gap, vertex-face gap).
- **OCCT** `Approx_SameParameter` [read]: when the 3D curve and pcurve are not parameterised
  alike, it projects samples to find the parameter correspondence. It then fits a cubic
  B-spline reparameterisation, recomposes the pcurve so it shares the 3D curve's parameter,
  and only then takes tolerance = 1.05 x deviation at 2 x NbSamples points. `ShapeFix` /
  `ShapeFix_ShapeTolerance::LimitTolerance` caps tolerances. [memory for the ShapeFix part]
- **Parasolid** [read]:
  - `PK_EDGE_optimise` lowers a tolerant edge's precision toward the true SP-curve deviation.
  - `PK_EDGE_reset_precision_2` deletes the SP-curves and re-intersects the surfaces, making the
    edge exact again (`int_only` option).
  - The docs warn that setting precision smaller without replacing the SP-curves fails checking.
  - Short edges, "all positions within vertex bubbles", are reported (`PK_ERROR_edge_too_short`).
  - Loop gaps can be closed by extending/trimming SP-curves or by reparameterising the surface.
- **ACIS** [read]: tweak re-intersects non-tangent tolerant edges to precise geometry. Tangent
  and near-tangent ones stay tolerant.

**Fit.** This adds two mechanisms to baseline T:
- The parameter-mismatch part of a gap is repairable by refitting the pcurve, so the tolerance
  only has to cover geometric mismatch.
- A "re-exactify" pass for edges whose two surfaces are analytic, using Catallax's cheap exact
  closed forms.

SMLib's bbox/800 cap and Parasolid's short-edge rule are concrete numbers for T's cap and
"no short edge" requirements.

## 5. Watertight by construction (change the surfaces so they share the edge)

**What it is.**
- Urick et al. 2019 post-process a conventional Boolean result in three stages: parameter-space
  analysis, reparameterisation, and model-space update. Surfaces are updated so that they
  "geometrically join explicitly at their approximated intersection", using the SSI tolerance
  the CAD system already uses. The result is gap-free and editable, built from standard
  splines. Follow-up: Urick et al. 2020 JCISE, "Reconstruction of trimmed NURBS surfaces for
  gap-free intersections".
- Sederberg et al. 2008 convert each trimmed NURBS to an untrimmed T-spline so that
  neighbours match exactly along the former trim. The geometry changes within tolerance.
- Hybrid: "Making trimmed B-spline B-reps watertight with a hybrid representation" (about 2019)
  seals narrow bands around trims with Bezier triangles.

**Evidence.** The authors claim watertight output, and IGA simulation results support analysis
suitability. I found no speed numbers against tolerant modelling.

**Cost/risk.** These methods rebuild or reparameterise surfaces near every edge. Control-point
counts grow, especially for T-splines, and STEP has no T-spline entity. They also change
surfaces the user may want left alone.

**Fit to T.** Low as a general mechanism. The ideas are useful in one place: Catallax's own
fillet and loft output, where the kernel builds the surface and can choose it so the boundary
is an isoline matching a pcurve of the neighbour. In that case one side of the gap goes to zero
by construction.

**Sources.** IGA2019 abstract [read]. CAD 115 (2019) 147-160, abstract and highlights via
search [partly read]. Sederberg 2008 ACM TOG 27(3) [unread/memory]. Hybrid paper via search
snippet [unread].

## 6. Gap-tolerant downstream queries (reduce what the tolerance must protect)

**What it is.** Generalised winding numbers for collections of trimmed NURBS (arXiv 2504.11435,
2025) give inside/outside classification that degrades smoothly with gaps. A 2025 ACM TOG paper
("Boolean operation for CAD models using a hybrid representation", doi 10.1145/3730908) maps
the B-rep to a mesh with a bijection and conservative intersection, comparing itself against
OCCT, ACIS and Rhino.

**Fit.** Optional. It makes point classification and tessellation robust to small leftover
gaps, so fewer downstream operations fail. It does not replace tolerances. [unread, search
snippets only]

## 7. Tolerance theory (Segal 1990, Fang/Bruderlin/Zhu 1993, Stewart, Qi/Shapiro)

These are mainly polyhedral. The shared idea is that each entity carries a tolerance region,
decisions must stay consistent with the regions, and on conflict you grow the regions and merge
topology. Jackson cites this "adaptive tolerance-based" line as the precursor of Parasolid's
scheme. Qi and Shapiro (epsilon-topology) formalise when tolerant solids are valid. These papers
justify T's vertex >= edge >= face ordering and its grow-then-merge rule. They add no cheaper
check mechanism. [unread/memory, except Jackson's summary of refs 6-7, which I read]

## 8. Hybrid B-rep plus SDF for gap handling

I found no paper that uses an SDF to repair or certify B-rep edge gaps. DualBrep (arXiv
2606.31579, 2026) uses SDF+UDF fields for generation and reverse engineering, not for gaps.
Catallax's SDF evaluator does not appear to have a role here beyond containment, where
winding numbers (section 6) are the better-studied tool. [search snippets only]

---

## What I could not verify
- Sederberg 2008, Filip 1986, DeRose 1993, Shewchuk, Cherchi, Segal, Fang, Qi/Shapiro and
  Stewart: from memory only.
- OCCT ShapeFix cap specifics: from memory.
- Urick 2019 full text: paywalled. I read only the abstract, highlights and the IGA abstract.
- Parasolid's procedural intersection-curve evaluation details: from XT reference snippets,
  not the full PDF.
- I did not run any experiment measuring derivative-bound sample counts on Catallax's loft
  fillet. That is the obvious next check before adopting approach 2.
