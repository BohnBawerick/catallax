# Spline-face rolling-ball fillets in the Catallax distance field

Worktree: Catallax repository at bb2f8375 (CX below).
Follows sdf-fillet-research.md in this scratchpad.
Labels: [read] = read this session; [snippet] = search-result summary only;
[unread/memory] = my own knowledge, not checked. (inference) marks my reasoning.

## 0. The previous study is out of date on one point

The previous note says M27's six-variable search refused at depth 24 and that nobody had a
certified spline spine. Both statements are now stale inside CX itself [read]:

- docs/M11-OFFSET-PROOF.md, "The spine by the M10 engine": since 12 November the M10 pair
  search takes exact offsets (`RationalPatch::offset(d)`, `OffsetPatch::offset_meeting`). It
  carries L = |N| as an unknown with L^2 = N.N and certifies cells only where the offset is
  proved regular. All three test edges (M27 loft edge, next loft edge, biquadratic edge) prove
  as one open branch of three certified arcs, 2.3M / 2.7M / 6.8M of 20M work units, depth 3,
  ends exactly isolated on support borders, identical at 1/2/4 workers.
- `SplineFillet::fit` (crates/catallax-shapes/src/spline_fillet.rs, 6.5k lines with submodules)
  then proves, per span of tau, a parametric Krawczyk certificate of a unique smooth spine point
  for every tau in the span, with interval Taylor jets of the implicit solution to fourth order.
  It fits a rational (D,2) carrier with Hermite rows and reports a two-sided bound: 78 nm, 5 nm,
  80 nm on the three edges; total work 5.3M / 7.0M / 16.4M including the spine. Target 100 nm,
  hard cap 1 um, refusal `fillet.spline-deviation` with the bound reached.
- A numerical path exists too (spline_fillet/numerical.rs): Newton plus quarter-point checks,
  "not a complete branch or Hausdorff proof" [read].
- The M27 and M11 depth-24 refusals were the six- and four-variable regular-graph experiments.
  They remain baselines, not the production route [read].

So "fit the ball-centre path once with a proven eps" is not hypothetical in CX. The B-rep does it
today with a sound a priori bound (Krawczyk existence/uniqueness plus Taylor remainder).

## 1. Prerequisite the previous study did not name

The field does not hold general spline faces. crates/catallax-field/src/import_spline.rs [read]
admits only bilinear equal-weight patches: flat ones as planes and curved ones as saddle roofs by
their exact quadric. Everything else is `field.import-spline`. Field lofts are listed under 4.2
("field versions of fillets, sweeps and lofts", ROADMAP.md:196 [read]); imported spline parts
in the field are 4.4 (ROADMAP.md:198 [read]). A spline fillet in the field needs distance to its
spline supports first, for the body around the fillet and for the trimmed contact rails.
Step 7 therefore depends on the field's loft step (or whichever step lands certified
point-to-rational-patch distance).

## 2. Literature by route

### Exact offsets
- Offsets of polynomial/rational surfaces are not rational in general because of |N|; only PN
  (Pythagorean-normal) surfaces have rational offsets. Pottmann, "Rational curves and surfaces
  with rational offsets", CAGD 12(2) 1995; Peternell and Pottmann, Laguerre approach, CAGD 15,
  1998. Finding all polynomial PN surfaces is still open (2012 statement). [snippet]
  https://www.geometrie.tuwien.ac.at/peternell/laggeo2.pdf
- Maekawa, "An overview of offset curves and surfaces", CAD 31(3) 165-173, 1999: survey of exact
  representation, approximation, self-intersection. [snippet, content unread]
- Consequence (inference): no closed-form spine for general B-spline pairs. CX supports are
  arbitrary rational patches, so the PN route is closed.

### Approximate offsets with error bounds
- Elber and Cohen, "Error bounded variable distance offset operator for free form curves and
  surfaces", IJCGA 1(1) 1991; Elber, Lee, Kim, "Comparing offset curve approximation methods",
  IEEE CG&A 17(3) 1997; Lee, Kim, Elber, circle approximation, CAD 28(8) 1996; Piegl and Tiller,
  "Computing offsets of NURBS curves and surfaces", CAD 31(2) 1999. [snippet]
- Their bounds are on each offset separately [unread/memory]. Intersecting two approximate
  offsets adds the angle amplification 1/sin(theta) between offset normals (inference), and still
  needs a certified intersection. Strictly worse than solving the exact offset system, which CX
  already does.

### Subdivision / Bernstein solvers for the offset intersection
- Hu, Maekawa, Patrikalakis, Ye, "Robust interval algorithm for surface intersections", CAD 29(9)
  1997: interval projected polyhedron (IPP), Bernstein subdivision with rounded interval
  arithmetic, handles ill-conditioned tangential cases in theory. [snippet]
  https://www.sciencedirect.com/science/article/abs/pii/S0010448596000991
- Elber and Kim 2001: zero-dimensional solutions certified; one-dimensional sets fitted to
  samples (previous note). [snippet]
- Mourrain and Pavone, Bernstein subdivision solver (JSC 2009) [unread/memory].
- CX equivalent: M10 pair search with offset inputs. Proven on the three edges [read].

### Certified tracing of the spine (one-manifold)
- Martin, Goldsztejn, Granvilliers, Jermann, "Certified parallelotope continuation for
  one-manifolds", SINUM 51(6) 2013: each parallelotope is crossed by exactly one curve piece,
  proved with interval Newton; code IbexCont. [snippet]
  https://pagesperso.ls2n.fr/~jermann-c/dw/_telech.php?pdf=mggj-sinum13.pdf
- Mukundan, Ko, Maekawa, Sakkalis, Patrikalakis, "Tracing surface intersections with validated
  ODE system solver", SM 2004: validated tube around the intersection curve. [snippet]
  https://web.mit.edu/deslab/pubs/files/SM04.pdf
- Duff and Lee, certified Krawczyk homotopy tracking, ISSAC 2024 (arXiv 2402.07053); a priori
  bounds, arXiv 2512.01355; Burr et al., certified surface approximation with the interval
  Krawczyk test, ISSAC 2026 (arXiv 2602.07718), Hausdorff bound from the Krawczyk boxes. [snippet]
- "Topology guaranteed and error controlled curve tracing for parametric surface-surface
  intersection", CAGD 2025 (S0167839625000214). [snippet]
- CX equivalent: the parametric Krawczyk plus Taylor jets inside spline_fillet [read]. This is
  the same family as ParCont and the validated ODE tube.

### Distance to a spline curve with a certified bound
- Standard (unread/memory): point-to-polynomial-curve distance is a 1D minimisation; Bernstein
  subdivision of |x - C(t)|^2 with hull exclusion, or an interval panel cover, gives an
  enclosure. CX rolling_region already does the panel cover form for a closed-form centre curve
  (96 bits, <= 4096 panels, typed refusals) [read header].
- Distance to the fillet face = |dist(x, spine) - r| where the nearest spine point is interior
  and the direction lies in the contact sector, while r is below the local reach (previous note,
  inference). With the spine inside a tube of radius eps the enclosure widens by eps each side.

### Closing/opening limited to a window
- Rossignac and Requicha 1984/1986 (previous note). Needs exact distance to trimmed spline
  supports (which the field lacks, section 1) inside a 3-unknown branch-and-bound per query.
  Each inner evaluation is itself a 2-parameter spline closest-point solve (inference).

## 3. Routes compared

| Route | Guarantee | Build per fillet | Query | Reuse in CX |
|---|---|---|---|---|
| A. Exact offset closed form (PN) | exact | n/a | n/a | none; general supports are not PN |
| B. Approximate each offset, intersect approximations | bounded, eps grows ~ 1/sin(angle) | moderate | cheap | none useful |
| C. Direct Bernstein search per query (no stored spine) | certified | none | very high: a 4-5 unknown search per point (inference) | M10 engine; too slow per point |
| D. Spine tube built once: M10 branch topology + field's own parametric Krawczyk with Taylor jets per tau span; query by panel cover | certified enclosure; width floor 2 eps, eps proven | ~ B-rep spine cost: 2-7M units topology plus span proofs, under the B-rep's 5-16M total (inference) | 1-param panel cover over spans: tens to hundreds of panels, est. 10-100 us (inference, not measured) | M10 `offset_meeting`, `BSplineSurface::jet_box`/`evaluate_box`, Ball/Decider, rolling_region panel cover and sectors |
| D'. Same as D, refine spans at query time | certified, width shrinks on demand | as D | higher near rails/ends (re-runs Krawczyk) | as D |
| E. Numerical spine (Newton samples, fitted) | estimate only, measured deviation | low | cheap | numerical.rs pattern (in B-rep, not to be shared) |
| F. Windowed closing/opening with exact trimmed spline distances | certified if inner distances certified | none | very high: 3D B&B with nested 2D spline minimisations, est. 1e3-1e5 inner solves | none yet; needs field spline face distance first |
| G. Voxel offsets (PicoGK style) | voxel size, rounds all edges | grid build | cheap | none; wrong semantics |

Spline-with-analytic pairs use D unchanged; the analytic side's offset is exact and closed form,
and the M10 slices for plane, cylinder and sphere against a patch already exist [read M27 doc].
One unknown pair fewer (inference).

## 4. Recommended contract

Under the owner's rulings (proof-cost.md [read]: 13 Nov "accept a measured threshold or change
strategy"; 14 Nov minimum options: "report what the checks establish", numerical B-reps under
stated tolerances), and FIELD.md's value kinds [read]:

Contract: **the field fits each spline fillet's ball-centre path once, with its own proven
Hausdorff bound eps, and returns `Value::TrueDistance` enclosures that contain the exact rolling-
ball distance with width at least 2 eps. eps must be no larger than min(100 nm, part tolerance /
2). A build that cannot prove that refuses with the eps it reached.** Optionally a second,
labelled tier: when the certificate refuses, fall back to a numerical spine with a measured
(not proved) deviation, returned as an approximate value, never as TrueDistance.

Why not a tighter "certified enclosure that shrinks on demand" (D'): a proven 2 eps width below
tolerance already is a certified enclosure of the true distance. Refining past the part tolerance
is proving for proving's sake. Why not numerical only (E): the field exists partly to check the
B-rep, and a numerical field spine cannot catch a B-rep spine error larger than its own unknown
error. The proven eps costs the same order as the B-rep's own fit, which CX already pays.

Independence (inference): sharing the M10 engine is sharing generic intersection code, like
offset_spline.rs sharing `certify_elliptic_offset_sweep`. Do not share spline_fillet's Krawczyk,
jets or carrier. The field writes its own span certificate (different variables are fine, e.g.
contact parameters on support A only plus the analytic offset equation) so a bug in the B-rep
fit cannot pass both. M10 branch counting stays a common-mode dependency; say so in the doc.

## 5. Step 7 of the field-fillet sub-milestone (4.2)

Prerequisite: the field step that gives certified distance to a trimmed rational spline face
(the 4.2 loft step). Step 7 does not start before it.

7a. Spine tube type: per tau span, enclosures of C(tau) as a Taylor model (order 3 coefficients
    plus a remainder ball), contact parameters on both supports, and a curvature bound from the
    jets. Built from M10 `offset_meeting` topology and the field's own parametric Krawczyk over
    `BSplineSurface::jet_box`. Halve spans until eps <= min(100 nm, tol/2); refuse with reached eps.
7b. Query: plug the tube into the step-1 generalised rolling_region panel cover as a new piece
    kind. Sector test from the contact enclosures. Reach check r * kappa_max < 1 from the jets.
7c. Spline-with-analytic: analytic offset in closed form, same tube and query.
7d. Ends: reuse step 4's edge-tube and vertex-sphere pieces where the spline spine meets a border.
7e. Typed refusals: `field.spline-fillet-spine` (topology), `field.spline-fillet-eps` (with
    reached eps), `field.spline-fillet-reach`, plus the singular-contact control.
7f. Optional numerical tier, labelled approximate, only if the owner wants spline parts to build
    in the field when the certificate refuses.

Done when: on the three M11 edges (M27 loft edge, next loft edge, biquadratic edge) and one
spline-plane edge, the field builds each fillet with no fillet refusal, its own eps is at most
min(100 nm, part tolerance / 2), and catallax-compare finds every B-rep display vertex of each
fillet face inside the field's TrueDistance enclosure widened by the B-rep's reported bound;
build work per edge stays inside the 20M allowance; records match at 1, 2 and 4 workers; a
singular-contact control refuses with a typed code; the median query time on each fillet is
measured and recorded in the milestone doc.

## Not verified
- Field build cost and query time for route D: estimates from the B-rep's measured work, not runs.
- Whether the field's loft step will provide point-to-spline-face distance, or how.
- Content of Maekawa 1999, Elber-Lee-Kim 1997 bounds, Mourrain-Pavone details (citations only).
- Whether M10's offset search handles an analytic-vs-patch offset pair directly (the doc shows
  patch-vs-patch offsets; plane/cylinder/sphere slices exist for non-offset inputs).
- Whether spline_fillet's span certificate covers tangent (G1) support pairs; M10 doc says an
  offset pair with a tangent contact refuses at the split limit.
