Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Roadmap](../../../../ROADMAP.md).

# Watertight fillets and copy tagging in Catallax

Worktree: Catallax repository, commit bb2f8375. Read-only study.
Paths below are relative to `crates/`. "Inferred" marks anything not read directly.
Builds on edge-copies-code.md in the same scratchpad.

## 0. Headline finding

The exact M11 spline fillet already is a watertight fillet in the sense of proposal 2.
`catallax-shapes/src/spline_fillet.rs:11-13`: "its two boundary rows are the supports composed
with cubic Hermite curves of their contact parameters, so each rail lies exactly on its support
and is the trim curve there". The code does it in `spline_fillet/carrier.rs:403-518`: cubic
Hermite rail parameters per span (`cubic`, carrier.rs:300-309), then `compose(poly, params)`
(carrier.rs:279-298) substitutes them into each support's exact homogeneous polynomial
(`support.exact`, carrier.rs:434-439), all in exact `Rational`. Every row is degree-elevated to
the largest degree (carrier.rs:454-466) and the carrier is a rational B-spline of degree
`[degree, 2]` (carrier.rs:498-517). Its STEP export writes the rail pcurve as "its exact cubic
parameter curves" (`catallax-step/src/export/spline_fillet.rs:1-9`).

So the open question is the numerical fillet, not the exact one.

## 1. How rolling-ball fillets are built today

There are three families.

### A. Exact elementary strips (fillet.rs)
- `catallax-features/src/fillet.rs:154-290`: match on (support, support, edge curve). Plane-plane
  with a line edge gives a `Cylinder` carrier and two `Line` rails (fillet.rs:171-213).
  Plane-cylinder circle gives a torus (fillet.rs:215-219). Oblique plane-cylinder/cone ellipse
  gives `RollingPipe` with `Ellipse` rails (fillet.rs:710-757). Cylinder-cylinder line gives
  `RollingSaddle` (fillet.rs:261, 369). Sphere-cylinder ring (fillet.rs:290-330).
- The rail is a closed-form 3D curve on both surfaces. The exact body stores only that 3D curve
  (entity.rs:44-62), so there is one copy and no gap. Watertight already.

### B. Exact M11 spline fillet (catallax-shapes/src/spline_fillet.rs)
1. Spine: centre where offsets `P + side r N/|N|` of two `RationalPatch` supports meet, a graph
   over one support parameter `tau` (spline_fillet.rs:3-5). Krawczyk certificate per span
   (spline_fillet.rs:17-20).
2. Contact parameters at nodes, Hermite cubic per span (carrier.rs:300-309, 415-422).
3. Rails = support polynomials composed with those cubics (carrier.rs:434-439). Exact.
4. Cross section: rational quadratic, middle row = weight cubic times middle-point cubic
   (carrier.rs:440-449), so degree 6 before elevation.
5. Carrier: B-spline `[degree, 2]` with knots of multiplicity `degree` at nodes
   (carrier.rs:494-517).
6. Rails proved inside support domains (`rails_inside`, carrier.rs:370-399). Rail border nodes
   pin a support parameter to a rational value at a node (`RailBorder`, spline_fillet.rs:144-157).
7. Error bound two-sided by Hermite remainders (spline_fillet.rs:20-24); halve spans until met.

### C. Numerical rolling fillet (fillet_rolling.rs, fillet_numerical.rs, spline_fillet/numerical*.rs)
1. Supports wrapped as `NumericalFilletSupport` over `Surface` (plane, cylinder, cone, sphere,
   torus, elliptic cylinder, rational atlas): numerical_support.rs:7-60, analytic proposals in
   numerical_analytic.rs. Quadrics use angle parameters (`u.sin_cos()`,
   numerical_analytic.rs:115-160).
2. Spine: Newton continuation across atlas charts (spline_fillet/numerical.rs:1-4).
3. Per station, `packed` (numerical_carrier.rs:31-81): binary64 contact points and normals from
   `contact_f64`, rounded to rationals, middle point and weight `sqrt((1+cos)/2)`, plus the
   contact parameters `node.parameters` in slots 12-15.
4. Slopes by a local 5-point Lagrange derivative of each packed coordinate separately
   (`station_slopes`, numerical_carrier.rs:84-121); cubic Hermite controls (`hermite`, :123-148).
5. `candidate` (numerical_carrier.rs:157-243): carrier `[3,2]`; rails = rows 0 and 2 of the
   carrier (rational cubic, weight 1 at contacts so polynomial); `parameter_rails` = a separate
   non-rational cubic of slots 12-15.
6. Fit loop checks centre, section and `rail_error` at quarter points and halves spans
   (numerical_carrier.rs:~263-385).
7. Topology: `constrained_rail` edits pcurve end controls to hit cut uv values
   (fillet_numerical.rs:423, 836-848); border edges get uv from source native curves via `uv_at`
   or plane inversion (fillet_numerical.rs:806-808); 3D edge = `fillet.rail(side)`
   (fillet_numerical.rs:849-851). Body rechecked by `rebuild_with_sampled_carriers`.

Why the numerical gap exists (inferred from 3-5): rail point and pcurve value at a node are
rounded independently (contact_f64 vs node.parameters), and slopes are differenced per
coordinate, so the rail slope is not `S_u u' + S_v v'`. Gap is about rounding at nodes and
O(h^4) between, then the endpoint edit moves the pcurve only.

## 2. Exact boundary per support kind (maths)

Write the support homogeneously, pcurve degree d (cubic: 3). Composition degree in tau =
d times the largest total degree of a monomial in the homogeneous numerator/weight.

| Support | Exact boundary possible | Rail degree in tau (d=3) | Notes |
|---|---|---|---|
| Plane | Yes, trivially | 3 (no increase) | S affine. Today's numerical carrier stays [3,2]; rail controls = S(pcurve controls). |
| Cylinder, angle param | Only for iso-line contacts (generator line, coaxial circle) | 1 in angle | cos(poly) is transcendental. Elementary native strips already cover these. |
| Cylinder, rational param (Weierstrass s=tan(u/2) or quarter-arc NURBS) | Yes | 9 | Monomial h s^2 has degree 3. Nodes needed at arc joints. |
| Cone, rational | Yes | 9 | Same structure as cylinder. |
| Sphere, rational biquadratic | Yes | 12 | Stereographic (total degree 2) gives 6, but a nonstandard chart. |
| Torus, rational biquadratic | Yes | 12 | s^2 t^2 term, degree 4. |
| Rational bicubic B-spline | Yes (M11 does it) | 18 | Biquintic gives 30. Nodes where the pcurve crosses support knot lines (RailBorder pattern). |
| Earlier fillet blend [n,2] as support | Yes | 3(n+2): 15 for n=3, 60 for n=18 | Fillet over fillet compounds degree. |

Blend degree = max(rail degree, middle row degree). The numerical carrier fits the homogeneous
middle row directly as a cubic, so the middle row adds nothing; M11 uses weight x point = 6.
Control points per span per row = degree, so a bicubic support costs 6x today's numerical count.

Watertight boundary removes only the copy gap. Ball tangency, radius error and G1 to supports
are still fit errors and still need checks.

## 3. Interference

- Later fillets: they read support native trims (fillet_numerical.rs:408-420, 806). An exact
  preimage pcurve is a better input. But the endpoint edit at fillet_numerical.rs:836-848 breaks
  exactness. A watertight build must cut by exact subdivision at `cuts[k].tau_ends()` and define
  the vertex as `S(p(tau_cut))`, reversing today's authority (vertex first, pcurve snapped).
  Fillet over fillet blows degree up (table).
- Booleans: numerical Booleans are scoped to loft vs planar stock and cone rings
  (catallax-boolean/src/numerical.rs header, partition.rs:1). They do not take fillet bodies
  today (inferred from scope). Splitting a watertight edge stays exact if both copies are split
  by knot insertion at the same tau (composition commutes with restriction). A refit breaks it.
- Quadric pcurves: the stored native trim is a polynomial B-spline in the angle chart
  (native.rs:13-26). A rational-chart pcurve does not fit that slot. Options: store the support
  as a rational B-spline surface (loses analytic type), or keep angle charts and accept
  approximation there.
- STEP: numerical writer emits both copies (curved.rs:146-188) and rechecks rounded pcurves
  (analytic.rs:620-656). Decimal rounding reopens a gap at ulp scale times the weight
  condition; the self-check must stay. STEP CYLINDRICAL_SURFACE etc. use angle parameters, so an
  exact rational-chart pcurve cannot be written on an analytic face; the exact path solves this
  with pcurves derived at write time (catallax-step/src/pcurve.rs:1-2). Degree 18 is already
  exported for M11. OCCT's B-spline degree cap of 25 (from memory, not read) rules out
  biquintic supports and cascades without degree reduction, which is approximation again.
- Exact pipeline: unaffected; it already works this way (sections 0, 1A, 1B).
- Tolerance: the numerical body has one body-wide `tolerance_mm` (default 0.01,
  topology/src/numerical.rs:107-113). Per-entity tolerances exist in tolerant import only. Under
  a per-entity model with 1.1x margin, 1.1 x measured 0 = 0, which fails any receiver after
  rounding. A watertight edge needs tolerance from a rounding bound, not a measurement:
  `max(1.1 * measured, rounding_bound)`.

## 4. Do the two ideas combine

Yes. A watertight support pcurve is a Derived copy with bound 0 (exact rationals), the strongest
case. Nothing breaks if:
- the tag is cleared by any edit to one copy (constrained_rail, uv endpoint rewrite,
  vertex repair at topology/src/numerical.rs:283-305);
- vertex, loop and cap checks still run on tagged edges;
- the STEP self-check still runs (rounding is not covered by the tag);
- tolerance uses the floor in section 3.
The payoff of skipping is small (the check is cheap next to the fit). The real value of the tag
is recorded authority, so consumers know which copy to trust and edit.

## 5. Idea 1 in Catallax

- Tag: a field on `NumericalTrim` and `NumericalNativeTrim` (native.rs:13-26), for example
  `Independent` or `Derived { from: Edge3d | Carrier(face), rule, bound: Rational }`. Constructors
  require it; no default to Derived.
- Setters: numerical fillet (blend-face iso-line trims: Derived from carrier; support pcurve:
  Derived only if built by exact composition, else Independent), embedding of analytic lines
  and coaxial circles, Boolean pieces split by knot insertion. Plane projected trims stay
  checked: that check is the only proof the 3D curve lies on the plane.
- Readers: numerical checker (numerical.rs:828-904) records `bound` instead of sampling for
  Derived; samples Independent copies at points the fitter never used (seeded random or
  irrational offsets, not eighths, see numerical_embedding.rs:30-35). STEP writer reads the
  bound for its tolerance.
- On edit: any mutation of one copy resets the tag to Independent, unless the same exact map is
  applied to both copies (knot insertion, rigid transform in rationals).
- Export: STEP carries no tag. Both copies are rounded and self-checked as today.
- Reopen: documents are recipes and handles are not persisted (catallax-app/src/lib.rs:3-8), so
  tags regenerate on rebuild. Any future body cache must include the tag.
- Import: all file copies are Independent. Tolerant admission unchanged.
- Wrong tag: an unchecked gap ships to STEP and to later fillets, silently. Catch it with a
  debug/CI audit that runs the full check on Derived copies and asserts measured <= bound plus
  rounding, plus a property test per producer.
