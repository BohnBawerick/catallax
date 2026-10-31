# Piece 4: exact intersections (`catallax-intersect`)

Status: complete for blind test writing and blind building.

This is the input the test writer and the creator get, together with `docs/DESIGN.md`, `docs/ARCHITECTURE.md`, piece 1's description `docs/pieces/01-numbers-and-decisions.md` and piece 2's description `docs/pieces/02-shapes.md`. This piece stores everything it makes in piece 2's `GeomStore` and decides through piece 1's `Decider`, exactly as those descriptions define them. Nothing else is needed.

## Approved conic and algebraic migration

The 2026-10-22 architecture decision supersedes the rational-only types, line-only open-curve operations and closed-form restrictions below. The existing pair tables remain the rational baseline, not a claim of coverage for newly representable surfaces. Implementation lands in two steps.

### Conic preparation

First widen `CurveForm::Circle::radius_sq` and `CurveForm::Ellipse::{major_sq, minor_sq}` to `Real`. Add these forms, with rational directions in this preparation step:

```rust
Parabola {
    vertex: Point3, normal: RVector3, axis_dir: RVector3, focal: Real,
}
Hyperbola {
    center: Point3, normal: RVector3, axis_dir: RVector3,
    transverse_sq: Real, conjugate_sq: Real, branch: StrictSign,
}
```

For either form, `e1` is the normalised `axis_dir`, `e3` the normalised `normal`, and `e2 = e3 x e1`; those directions are nonzero and perpendicular. A parabola requires `focal > 0` and has `C(t) = vertex + focal (t^2 e1 + 2t e2)`, for every real `t`. Its equations are its plane and `y^2 - 4 focal x = 0`, with `x,y` measured from `vertex` in this frame. A hyperbola requires both squared dimensions positive. With `a = sqrt(transverse_sq)`, `b = sqrt(conjugate_sq)` and `s` the sign of `branch`, use `C(t) = center + s a sqrt(1+t^2) e1 + b t e2`, for every real `t`. Its equations are its plane and `x^2/a^2 - y^2/b^2 - 1 = 0`, with the additional selector `sign(x) = branch`. Two hyperbola branches are different components even though their equations agree. Increasing `t` follows the natural direction; choose the frame accordingly from the oriented intersection. These are forms of stored intersection components, not new piece 2 curve variants.

Add `curve_is_closed(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<bool, Refusal>`. Lines, parabolas and hyperbola branches are open; complete circles and ellipses are closed. For `Branch`, certify the actual component's closure, including singular endpoints. Unknown closure returns a specific existing numeric or unsupported-curve failure, never a guessed boolean. `Curve::is_closed()` remains the cheap `Option<bool>` query on the stored definition. Consumers needing a definite answer use the piece 4 query.

Extend `order_on_curve` to open conics using their increasing parameter. `cyclic_order` requires proven closure. Open-curve `arc_point` and `trace` require two finite endpoints in natural order; `None, None` is valid only for a proven closed curve. Branch membership includes the component/nappe selector, not just the implicit equations. At a conic endpoint or singularity preserve the stored component and one-sided tangent contract. An unresolved selector cannot supply a marked membership or topology choice.

Keep every existing successful rational answer byte-identical through this preparation. Adapt matches in boolean, measurement, export and CLI together. The unchanged version 2 exchange has only `line`, `circle`, `ellipse` and `intersection` edge tags, with `intersection` restricted to branches without a closed form. It cannot truthfully report a parabola or hyperbola form. Keep the new forms available through library/recipe evaluation; publishing those results through the CLI requires a separate versioned exchange decision. Do not invent a tag or silently relabel a known conic as a branch. This is an exchange coverage gap, not impossible geometry.

### Shared algebraic foundation

After conic preparation, widen every `CurveForm` direction to `Vector3`, including the new forms. Positions remain `Point3`; all dimensions are `Real`. Stored analytic curves follow piece 2's matching field changes. Configuration quantities become `Real`; use rational fast paths when available, otherwise the decider. `Branch` now means a certified intersection component without a recognised analytic form, not only a quartic curve. Its coefficients and equations use `RealPoly` and piece 1's certified algebraic solver path. Never clear an algebraic coefficient by replacing it with a decimal or by accepting all conjugate roots.

The required new coverage includes cone/plane open conics, torus/plane circles with irrational squared radii, and elliptic-cylinder/plane sections. `supports` identifies kind pairs with an implemented path, not a guarantee for every pose. General cone, torus and elliptic-cylinder coverage still requires component isolation, singular handling, membership, ordering and certified traces. A surface enum variant alone is not coverage. Missing paths retain specific temporary failures; every coverage change updates the affected baseline table and tests.

Preserve nappe selectors, shared root correlations and construction incidences through all operations. Surface-pair and hit ordering must use exact comparisons in the documented order. Pair reversal changes orientation as specified, never the selected point sets. Traces enclose coefficient uncertainty and chart derivatives at the caller's precision and retain finite trims on open curves.

Use independent section equations to check torus cuts with irrational squared radii, parabola orientation and both hyperbola branches. Test that a point on the opposite branch fails membership and that an open component cannot become a ring edge. Add irrational axis-length ratios to the CC6 cases when their coverage lands. A high-degree algebraic expression must take the shared proof path, not a rational-only special-case refusal.

## Finite spline curves and explicit coverage

The `BSpline` curve form retains real finite spline geometry, not the original
ellipse's implicit equations. Membership recovers a checked unique finite
preimage. Ordering and finite traces use its actual parameter domain, including
closed seam identification; endpoint queries preserve both arriving and
leaving tangents at a knot. A zero polynomial cannot stand for missing curve
equations.

The first spline support's producer implements certified axial planar trimming
by retaining the exact target and rebuilding the finite product solid. This
does not implement general spline/spline or spline/analytic intersection
isolation: unsupported pairs retain an explicit coverage refusal.

`intersect_curve_surface` answers an exact, unmarked `Disjoint` for a line that
`SplineSupport::misses_line` proves clear of the complete finite carrier, after
curve validation and recorded facts. Any other line against a spline support
keeps `intersect.parametric-surface`, including a line that really meets the
patch. There is no finite spline solver yet, so an empty hit list always means
a proved miss.

`trace_elliptic_offset_quarter` traces the retained exact normal-offset law on
the principal quarter `0 <= theta <= pi/2` for nominal property integration.
It is not a trace of the fitted boundary. The fitter's centred quarter charts
and their finite spline parameters have a different phase; callers must not
substitute one for the other.


## What this piece is for

Piece 4 answers "where do these meet?" exactly, for the surfaces of the first milestone: planes, cylinders and spheres, in every pair and every pose.

- Where two surfaces meet: nothing, the same point set, or a list of components (curves and isolated points). Each curve is added to the store as a definition: a line, circle or ellipse with exact rational data where one exists, and an intersection curve (piece 2's `Curve::Intersection`) otherwise.
- Where a stored curve meets a surface: nothing, the whole curve, or a list of points along the curve, each added as a meeting point (`PointDef::Meet`).
- Questions about stored curves: the closed form a curve has, whether a point lies on it, whether two curves are the same, the order of points along it, its tangents, a point strictly inside an arc.
- Certified traces of arcs: parametrised pieces with proven `Ball` enclosures of the point and its Taylor coefficients, for measurement (piece 6) and export (piece 7).

It works on geometry only. It never sees faces or edges; trimming curves to faces is piece 5's job.

## Number types

| What | Type | Meaning |
| --- | --- | --- |
| Configuration quantities: dot and cross products, squared distances, squared radii | `Rational` | Exact. Every configuration of a plane, cylinder and sphere pair is decided by the sign of a rational, with no `Decider` call |
| Fields of the lines, circles and ellipses this piece adds | `Rational`, `RPoint3`, `RVector3` | Exact rationals, as piece 2's analytic curves require |
| Positions of the points this piece adds (witnesses, singular points, tangent points, hits) | `Point3`, whose coordinates are `Real` | Exact real numbers built from rationals, square roots and `Algebraic` leaves, or certified system roots |
| Closed forms returned by `curve_form` | `Point3` for points, `RVector3` for directions, `Real` for squared lengths and focal length | Exact; points and dimensions may be irrational, directions remain rational |
| Traces | `Ball` | Proven enclosures at the caller's precision |
| Parameter values inside traces | `Rational` ranges, `Ball` arguments | Caches for measurement and export; no decision reads them |

A conic's squared dimensions are exact `Real` expressions. Primitive analytic curves in piece 2 still have rational squared dimensions. Extract a `Rational` with `Real::as_rational` only when storing such a primitive; an irrational dimension stays in its intersection definition.
Directions are never normalised; unit vectors appear only inside `Real` and `Ball` expressions.
For a vector `v`, `|v|` means `sqrt(v . v)` as a `Real`, and `v^` means `v / |v|`.

## Coverage

`supports(a, b)` is true when both kinds are among `Plane`, `Cylinder`, `Cone`, `Sphere` and `Torus`, and for plane/elliptic-cylinder pairs. Plane/cylinder/sphere paths retain their existing coverage. Cone/torus paths cover the closed-form configurations below, not arbitrary intersecting poses.

| Pair | Current cone/torus coverage |
| --- | --- |
| Plane–cone | Circles, ellipses, parabolas, the selected positive-nappe hyperbola branch, apex points and generator rays, subject to the native compatibility boundary below |
| Cone–cylinder, cone–sphere, cone–cone | Coaxial circles, tangent circles, isolated apex points, coincidence and separation |
| Plane–torus | Perpendicular circular sections, meridian circles through the axis, supporting-plane separation and isolated supporting tangency |
| Cylinder–torus, sphere–torus, cone–torus, torus–torus | Coaxial circles, tangencies and separation |
| Sphere–torus off axis | Certified containment/separation bounds, including a sphere inside the tube, and the meridian tangency circle of a sphere centred on the core |

Unimplemented intersecting poses return `intersect.traced-pair-required`. The rational fast path retains the established `intersect.algebraic-direction-required` cases so existing native answers stay unchanged. Algebraic definitions use `Real` directions and certified coefficient decisions rather than that rational-only direction restriction. These temporary coverage codes do not mean impossible geometry.

The Boolean caller can exclude finite parallel cone bands, or a cone band and
a parallel cylinder, without tracing their unbounded supports. It first
certifies each cone band between its complete coaxial rings, then bounds its
radius by its two ring radii. A band that reaches its apex has a smaller ring
of radius zero. Strict external tube separation, a tube strictly inside the
band's smaller ring, or strict containment inside a cylinder, excludes a
meeting. Axis distances and exact cone probe frames retain these proofs under
rational rigid placements. Touching or overlapping tubes do not authorize an
exclusion; truly meeting cone pairs still require tracing.

Curve queries include closed conics from these pairs and finite open-conic arcs. Membership retains the cone nappe, generator-ray and component selectors. A polynomial carrier alone does not establish membership. An unresolved cone nappe or new conic divider makes hit construction refuse with the originating numeric code; it cannot retain the other component with a mark. Legacy bounded divider marks are unchanged.

## Where two surfaces meet

`intersect_surfaces(geom, a, b, dec)` returns `SurfacePair`:

- `Disjoint`: the surfaces have no common point;
- `Coincident { same_sense }`: they are the same point set; `same_sense` is true when their normals agree (for two planes, `normal_a . normal_b > 0`; two cylinders and two spheres always agree);
- `Meets(components)`: the meeting set is the union of the components, in the order of the tables below.

A component is either `Component::Curve { curve, contact, singular }` or `Component::Point { point }`.

- `curve` is a new curve in the store.
- `contact` is `Crossing` when the surfaces cross along the curve and `Tangent` when they touch along it without crossing.
- `singular` lists the singular points of the meeting set that lie on the curve, in the order of the tables. A singular point is a point of both surfaces where their gradients are parallel (piece 2). A point shared by two components is one `PointId`, listed by both.
- `point` is an isolated touching point, a new point in the store.

### Rules shared by every configuration

**Exact configuration.** The plane/cylinder/sphere tables decide configuration boundaries with exact rational arithmetic, without a `Decider` call, statistics or marks, even after the deadline. Their existing exception remains the CC6 equal-radius crossing-cylinder pair with an irrational axis-length ratio: it refuses with `intersect.unsupported-pair` before adding geometry. Question operands do not meet that exception because their axes are rational unit vectors. New cone/torus paths use the caller's decider where an algebraic sign is needed; a refused pair adds no partial components and returns no marks.

**When there is a rational closed form.** A curve component that is a line, circle or ellipse whose fields, as given by the tables, are all rational is added as that analytic curve. Otherwise it is added as an intersection curve `Curve::Intersection { surfaces: [a, b], witness }` with the witness the tables give, and `curve_form` returns its closed form. A table says when a field can be irrational: a line origin `F + sqrt(h2) v` is rational exactly when the rational `h2` is the square of a rational, and similarly for the other square roots. Tangent curves always have rational fields.

**Natural direction.** An intersection curve runs along `grad f_a x grad f_b` (piece 2), with `f_a` the implicit of the first argument. Every analytic curve this piece adds for a crossing component runs the same way, so the tables give each line's `dir` and each circle's and ellipse's `normal` with the sign that makes this true. Two exceptions, where `grad f_a x grad f_b` vanishes or changes sign: tangent curves, whose direction the tables fix by an axis; and the two ellipses of configuration CC6, which pass through singular points.

**Points.** Every point this piece adds in `intersect_surfaces` is `PointDef::Meet(Meet { of: MeetOf::Surfaces(vec![a, b]), position: Position::Closed(p) })` with the exact position `p` of the tables.

**What a call adds.** One curve per curve component, one point per point component, one point per distinct singular point, and one witness point per intersection curve. Nothing else. A second call on the same pair adds a second, equal set of definitions; nothing is shared between calls.

**Facts.** Piece 2 records the incidences of every intersection curve on its surfaces and of every meeting point on its surfaces. This piece also records, with origin `FactOrigin::Construction("intersect.surfaces")` and `at: None` unless stated:

- `Incident(curve, a)` and `Incident(curve, b)` for every analytic curve it adds;
- `Incident(witness, curve)` for every intersection curve, after the curve is added;
- `Incident(point, curve)` for every singular point and every curve that lists it;
- `Tangent(Surface(a), Surface(b))` with `at: Some(Curve(curve))` or `at: Some(Point(point))` for every tangent component;
- `Coincident(Surface(a), Surface(b))` when the pair construction proves coincidence.

Circular components also acquire any structurally guaranteed coaxial or concentric relations from the geometry store. Pair construction consumes recorded coincidence, concentric sphere centres and coaxial cylinders before recomputing placement. The remaining radius or orientation comparison is still exact. Tangent queries skip the gradient-cross zero proof only at the recorded point, on the recorded component, or at a point recorded incident on that component. A tangent fact about one component never settles another component's contact.

**Order of the arguments.** Swapping `a` and `b` gives the same configuration and the same meeting set, and each intersection curve's `surfaces` swapped. The fields, the order and the directions of the components are those the tables give for the new order, so they can differ: a line's origin (PP1, CP4, CP5); the `dir` of the tangent line of CP4 and CP5 when the axes point opposite ways; the order of the two lines of CP3 when `w_a . w_b > 0`; the `ref_dir` of the SS3 circle, which becomes `perp(-u) = -perp(u)` and so moves the parameter origin to the opposite point of the circle; a witness of CC3 when the radii are equal; and the second ellipse of CC6, whose normal `n+` is symmetric in the two cylinders while its `major_dir` changes sign. Tests of a swapped call compare against the tables evaluated with the swapped arguments, never against the unswapped answer with directions reversed. Where a table needs a role ("the plane", "the cylinder", "the larger cylinder"), it says how the role is chosen.

### Notation for the tables

Planes have `origin o`, `normal n`, `ref_dir`. Cylinders have `origin q`, `axis w`, `radius r`, `ref_dir`. Spheres have `center c`, `radius R`. Subscripts `a` and `b` name the arguments when both have the same kind.
`perp(v)` is `v x e_k`, where `e_k` is the coordinate unit vector of the first axis (`x`, then `y`, then `z`) at which `|v_k|` is smallest. It is rational and non-zero for every non-zero `v`.

### Plane and plane

`m = n_a x n_b`.

| Configuration | Condition | Answer |
| --- | --- | --- |
| PP1 crossing | `m != 0` | One line: `origin = o_a - (((o_b - o_a) . n_b) / \|m\|^2) (n_a x m)`, `dir = m`. `Crossing` |
| PP2 apart | `m = 0`, `(o_b - o_a) . n_a != 0` | `Disjoint` |
| PP3 coincident | `m = 0`, `(o_b - o_a) . n_a = 0` | `Coincident { same_sense: n_a . n_b > 0 }` |

The PP1 origin is the point of the line nearest `o_a`.

### Plane and cylinder

The plane and the cylinder are whichever arguments have those kinds. `s = n . w`. When `s = 0`, `F = q - (((q - o) . n) / |n|^2) n` is the foot of the axis on the plane, `e2 = ((q - o) . n)^2 / |n|^2` the squared distance from the axis to the plane, `v = n x w`, and `h2 = (r^2 - e2) / |v|^2`.
`sigma` is `+1` when the plane is argument `a` and `-1` when it is argument `b`.

| Configuration | Condition | Answer |
| --- | --- | --- |
| PC1 circle | `w x n = 0` | One circle: `center = q + (((o - q) . n) / s) w`, `normal = sigma n`, `radius_sq = r^2`, `ref_dir` = the cylinder's `ref_dir`. `Crossing` |
| PC2 ellipse | `s != 0`, `w x n != 0` | One ellipse: `center` as PC1, `normal = sigma n`, `major_dir = w - ((w . n) / \|n\|^2) n`, `major_sq = r^2 \|n\|^2 \|w\|^2 / s^2`, `minor_sq = r^2`. `Crossing` |
| PC3 two lines | `s = 0`, `e2 < r^2` | First the line through `F + sqrt(h2) v` with `dir = -sigma w`, then the line through `F - sqrt(h2) v` with `dir = sigma w`. Each `Crossing`. Rational exactly when `h2` is a rational square; otherwise each is an intersection curve whose witness is that point |
| PC4 tangent line | `s = 0`, `e2 = r^2` | One line: `origin = F`, `dir = w`. `Tangent` |
| PC5 apart | `s = 0`, `e2 > r^2` | `Disjoint` |

### Plane and sphere

`k = (c - o) . n`, `F = c - (k / |n|^2) n` (the foot of the centre), `e2 = k^2 / |n|^2`, and `sigma` as above.

| Configuration | Condition | Answer |
| --- | --- | --- |
| PS1 circle | `e2 < R^2` | One circle: `center = F`, `normal = sigma n`, `radius_sq = R^2 - e2`, `ref_dir` = the plane's `ref_dir`. `Crossing` |
| PS2 tangent point | `e2 = R^2` | One point at `F` |
| PS3 apart | `e2 > R^2` | `Disjoint` |

### Sphere and sphere

`u = c_b - c_a`, `D2 = |u|^2`.

| Configuration | Condition | Answer |
| --- | --- | --- |
| SS1 coincident | `D2 = 0`, `R_a = R_b` | `Coincident { same_sense: true }` |
| SS2 concentric | `D2 = 0`, `R_a != R_b` | `Disjoint` |
| SS3 circle | `(R_a - R_b)^2 < D2 < (R_a + R_b)^2` | One circle: `center = c_a + t u` with `t = (D2 + R_a^2 - R_b^2) / (2 D2)`, `normal = u`, `radius_sq = R_a^2 - t^2 D2`, `ref_dir = perp(u)`. `Crossing` |
| SS4 outside tangent | `D2 = (R_a + R_b)^2` | One point at `c_a + (R_a / (R_a + R_b)) u` |
| SS5 inside tangent | `D2 > 0`, `D2 = (R_a - R_b)^2` | One point at `c_a + (R_a / (R_a - R_b)) u` |
| SS6 apart or nested | `D2 > (R_a + R_b)^2`, or `0 < D2 < (R_a - R_b)^2` | `Disjoint` |

### Cylinder and cylinder, parallel axes

`w_a x w_b = 0`. `u = (q_b - q_a) - (((q_b - q_a) . w_a) / |w_a|^2) w_a` is the offset from axis `a` to axis `b`, `D2 = |u|^2`. For two lines, `t = (D2 + r_a^2 - r_b^2) / (2 D2)`, `v = w_a x u` and `h2 = (r_a^2 - t^2 D2) / |v|^2`.

| Configuration | Condition | Answer |
| --- | --- | --- |
| CP1 coincident | `D2 = 0`, `r_a = r_b` | `Coincident { same_sense: true }` |
| CP2 coaxial | `D2 = 0`, `r_a != r_b` | `Disjoint` |
| CP3 two lines | `(r_a - r_b)^2 < D2 < (r_a + r_b)^2` | First the line through `q_a + t u + sqrt(h2) v` with `dir = w_a`, then the line through `q_a + t u - sqrt(h2) v` with `dir = -w_a`. Each `Crossing`. Rational exactly when `h2` is a rational square; otherwise each is an intersection curve whose witness is that point |
| CP4 outside tangent line | `D2 = (r_a + r_b)^2` | One line: `origin = q_a + (r_a / (r_a + r_b)) u`, `dir = w_a`. `Tangent` |
| CP5 inside tangent line | `D2 > 0`, `D2 = (r_a - r_b)^2` | One line: `origin = q_a + (r_a / (r_a - r_b)) u`, `dir = w_a`. `Tangent` |
| CP6 apart or nested | otherwise | `Disjoint` |

### Cylinder and cylinder, crossing or skew axes

`w_a x w_b != 0`. Call `B` the cylinder with the larger radius and `S` the other; when the radii are equal, `B` is argument `a`. So `r_B >= r_S`.

- `m = w_B x w_S`, `d = q_S - q_B`.
- The feet of the common perpendicular: `f_B = q_B + ((d x w_S) . m / |m|^2) w_B` on axis `B` and `f_S = q_S + ((d x w_B) . m / |m|^2) w_S` on axis `S`. Both are rational.
- `delta2 = (d . m)^2 / |m|^2` is the squared distance between the axes, and `delta = sqrt(delta2)`.
- The frame: `e3 = w_B^`; `e2 = (f_S - f_B) / delta` when `delta2 > 0`, and `m^` when `delta2 = 0`; `e1 = e2 x e3`. In it, `w_S^ = alpha e1 + gamma e3` with `alpha = w_S^ . e1` (never zero) and `gamma = w_S^ . e3`.
- A point `f_B + x e1 + y e2 + z e3` lies on both cylinders exactly when `x^2 + y^2 = r_B^2` and `(alpha z - gamma x)^2 = r_S^2 - (y - delta)^2`. So the configurations depend on `r_B`, `r_S` and `delta` only, never on the angle between the axes.

| Configuration | Condition | Answer |
| --- | --- | --- |
| CC1 apart | `delta2 > (r_B + r_S)^2` | `Disjoint` |
| CC2 outside tangent point | `delta2 = (r_B + r_S)^2` | One point at `f_B + (r_B / (r_B + r_S)) (f_S - f_B)` |
| CC3 one branch | `(r_B - r_S)^2 < delta2 < (r_B + r_S)^2` | One intersection curve, `Crossing`, no singular point. Witness: `f_B + r_B e2 + (sqrt(r_S^2 - (r_B - delta)^2) / alpha) e3` |
| CC4 two branches through a node | `r_B > r_S`, `delta2 = (r_B - r_S)^2` | Two intersection curves, each `Crossing`, each with `singular = [node]`, `node = f_B + (r_B / (r_B - r_S)) (f_S - f_B)`. First the branch with witness `W+`, then the one with witness `W-` |
| CC5 two branches | `r_B > r_S`, `delta2 < (r_B - r_S)^2` | Two intersection curves, each `Crossing`, no singular point. First the branch with witness `W+`, then the one with witness `W-` |
| CC6 two ellipses | `delta2 = 0`, `r_B = r_S` | Two ellipses crossing at two singular points; see below |

In CC4 and CC5, with `X = sqrt(r_B^2 - delta2)`, the witnesses are `W+ = f_B + X e1 + (f_S - f_B) + ((gamma X + r_S) / alpha) e3` and `W- = f_B - X e1 + (f_S - f_B) + ((r_S - gamma X) / alpha) e3`. Each lies where `y = delta`, so `(alpha z - gamma x)^2 = r_S^2`. The branch with witness `W+` is the one on the side `(p - f_B) . e1 > 0`.

In CC6 the axes meet at `O = f_B = f_S`. Let `k = |w_B| |w_S|`. When `k^2` is not the square of a rational, the pair is refused with `intersect.unsupported-pair`. Otherwise `k` is rational, `n- = |w_S|^2 w_B - k w_S` and `n+ = |w_S|^2 w_B + k w_S` are rational and perpendicular, and `cos = (w_B . w_S) / k`. With `r = r_B = r_S`:

1. the ellipse with `center = O`, `normal = n-`, `major_dir = n+`, `major_sq = 2 r^2 / (1 - cos)`, `minor_sq = r^2`;
2. the ellipse with `center = O`, `normal = n+`, `major_dir = n-`, `major_sq = 2 r^2 / (1 + cos)`, `minor_sq = r^2`.

Both are `Crossing`, and both list `singular = [O + sqrt(r^2 / |m|^2) m, O - sqrt(r^2 / |m|^2) m]`. These two ellipses lie in the planes `(p - O) . n- = 0` and `(p - O) . n+ = 0`, which hold the meeting set because `(p . w_B^)^2 = (p . w_S^)^2` there. Their normals are fixed as listed, not by `grad f_a x grad f_b`, which changes sign at the singular points.

### Cylinder and sphere

The cylinder and the sphere are whichever arguments have those kinds. `F = q + (((c - q) . w) / |w|^2) w` is the foot of the centre on the axis, `d2 = |c - F|^2`, `d = sqrt(d2)`. `sigma` is `+1` when the cylinder is argument `a` and `-1` when it is argument `b`.

Centre on the axis, `d2 = 0`, with `k2 = (R^2 - r^2) / |w|^2`:

| Configuration | Condition | Answer |
| --- | --- | --- |
| CS1 two circles | `R > r` | First the circle with `center = c + sqrt(k2) w`, `normal = -sigma w`; then the one with `center = c - sqrt(k2) w`, `normal = sigma w`. Both `radius_sq = r^2`, `ref_dir` = the cylinder's `ref_dir`, `Crossing`. Rational exactly when `k2` is a rational square; otherwise each is an intersection curve whose witness is its circle's parameter origin, `center + sqrt(r^2 / \|ref_dir\|^2) ref_dir` |
| CS2 tangent circle | `R = r` | One circle: `center = c`, `normal = w`, `radius_sq = r^2`, `ref_dir` = the cylinder's. `Tangent` |
| CS3 inside | `R < r` | `Disjoint` |

Centre off the axis, `d2 > 0`. The frame is `e3 = w^`, `e1 = (c - F) / d`, `e2 = e3 x e1`. A point `F + r cos(theta) e1 + r sin(theta) e2 + z e3` of the cylinder lies on the sphere exactly when `z^2 = K + 2 r d cos(theta)`, with `K = R^2 - r^2 - d2`. The comparisons below are rational, because `R > r + d` is `R - r > 0` and `(R - r)^2 > d2`, and so on.

| Configuration | Condition | Answer |
| --- | --- | --- |
| CS4 apart | `R < \|r - d\|` | `Disjoint` |
| CS5 tangent point | `R = \|r - d\|` | One point at `F + (r / d) (c - F)`. Here `d` is `r + R` or `r - R`, so the point is rational |
| CS6 one branch | `\|r - d\| < R < r + d` | One intersection curve, `Crossing`, no singular point. Witness: `F + (r / d) (c - F) + sqrt((R^2 - (r - d)^2) / \|w\|^2) w` |
| CS7 two branches through a node | `R = r + d` | Two intersection curves, each `Crossing`, each with `singular = [node]`, `node = F - (r / d) (c - F)`. First the upper branch, witness `F + (r / d) (c - F) + sqrt(4 r d / \|w\|^2) w`; then the lower, with `- sqrt(...)` |
| CS8 two branches | `R > r + d` | Two intersection curves, each `Crossing`. First the upper, witness `F + (r / d) (c - F) + sqrt((R^2 - (r - d)^2) / \|w\|^2) w`; then the lower, with `- sqrt(...)` |

The upper branch is the one on the side `(p - F) . w > 0`. In CS7, `d = R - r` is rational, so the node and the witnesses have one square root at most.

## Closed forms

`curve_form(geom, c, dec)` returns the exact shape of the carrier of a curve, as `CurveForm`:

- `Line { origin, dir }`, `Circle { center, normal, radius_sq, ref_dir }` or `Ellipse { center, normal, major_dir, major_sq, minor_sq }`, with the same meaning as piece 2's analytic curves and the same natural direction, except that `origin` and `center` are `Point3` and the squared dimensions are `Real`;
- `Parabola { vertex, normal, axis_dir, focal }`: the open conic `y² = 4 focal x`, with positive exact `Real` focal length. The unit `x` direction is `axis_dir`, the unit `y` direction is `normal x axis_dir`, and the origin is `vertex`;
- `Hyperbola { center, normal, axis_dir, transverse_sq, conjugate_sq, branch }`: one branch of `x²/transverse_sq - y²/conjugate_sq = 1`, with positive exact `Real` squared dimensions and `sign(x) = branch`. The unit `x` direction is `axis_dir`, the unit `y` direction is `normal x axis_dir`, and the origin is `center`. Opposite `StrictSign` selectors name different components;
- `Branch`: an intersection branch with no conic closed form.

For both open conics, the two defining directions are nonzero and perpendicular. Increasing `y` fixes the natural direction. Plane–cone sections now return these forms for their selected positive-nappe component; finite trims remain part of every open-arc query.

For a line, circle or ellipse it returns that curve's own fields and makes no decision. For an intersection curve on surfaces `[s0, s1]` with witness `p` it returns the form of the component of the tables that contains `p`. The configuration is exact rational arithmetic; which component holds `p` (and so the sign of a line's `dir` or a circle's `normal`) is one `Decider::sign` decision on the dividing quantity named below. When that decision ends `NotExact`, `curve_form` refuses with the `NotExact`'s code; a `Refused` is returned unchanged:

| Surfaces | Form |
| --- | --- |
| two planes | `Line { origin: p, dir: n0 x n1 }` |
| plane and cylinder with `s != 0` | the circle of PC1 or the ellipse of PC2, with `sigma` from the order of `s0` and `s1` |
| plane and cylinder with `s = 0` | `Line { origin: p, dir }`, with `dir = +w` or `-w`, whichever is the natural direction at `p`: the decision is the sign of `(p - F) . v`, and PC3 gives the direction for each sign |
| plane and sphere, two spheres | the circle of PS1 or SS3 |
| two parallel cylinders | `Line { origin: p, dir: +w0 or -w0 }`, by the sign of `(p - q_a - t u) . v` and CP3 |
| two cylinders with crossing or skew axes | `Branch`; on a CC6 pair the curve is refused (below) |
| cylinder and sphere, centre on the axis | `Circle { center: c + z w, normal, radius_sq: r^2, ref_dir }`, with `z = ((p - c) . w) / \|w\|^2`, the cylinder's `ref_dir`, and `normal` as CS1 gives it for the sign of `(p - c) . w` |
| cylinder and sphere, centre off the axis | `Branch` |

An intersection curve with a `Circle` form is the whole circle. A `Line` form is a carrier: cone generators remain rays from the apex, and `on_curve` checks the ray domain. Such a ray is not globally the same curve as the whole line. Parabolas and hyperbola branches likewise retain their selected domain. An intersection curve built by hand on a CC6 pair would be an open arc of an ellipse between singular points, with a direction the CC6 normals do not fix; it still refuses with `intersect.unsupported-curve`. `intersect_surfaces` stores those CC6 ellipses analytically instead.

Rational curve–surface elimination may skip finite-coordinate projections when
the carrier basis already proves both transverse coordinates affine in the
existing RUR separator. It retains that root polynomial and those coordinate
expressions: equal positions from another line parameter need not give identical
measurement or STEP enclosures. Public `curve_form` origins and natural
directions, component-domain checks, contacts and hit ordering are unchanged.

The front end counts edges by `curve_form`, so an intersection curve with a line form counts as a line in the manifest's `edges_by_kind`, and piece 7 writes it as a `LINE` (`JUDGE-INTERFACE.md`).

## Where a curve meets a surface

`intersect_curve_surface(geom, c, s, dec)` returns `CurveSurface`:

- `Contained`: every point of the curve lies on the surface;
- `Disjoint`: no point of the curve lies on the surface;
- `Meets(hits)`: the curve and the surface share exactly these points, each a `CurveHit { point, contact }`.

The points of the curve are those of its closure: for an intersection curve that ends at a singular point, the singular point counts.

**Facts first.** When `holds(Incident, Curve(c), Surface(s))` is true, the answer is `Contained`, settled by one `Decider::note_fact_decision` call and no arithmetic. This covers every intersection curve against either of its own surfaces.

**Circles far from a cone.** When a circle and a cone are not both rational, the circle's axis-aligned box, taken relative to the apex, is enclosed and the cone's equation is evaluated over it. A strictly signed enclosure answers `Disjoint` before any elimination: no point of the circle lies on either nappe. An enclosure that contains zero, including every touching or crossing circle, or one that cannot be settled, keeps the solver below. Rational pairs keep the solver and its accounting. This is what lets two blind drill points with irrational tip slopes share one hole feature.

**Contact.** A hit is `Crossing` when the sign of `f_s` (the surface's implicit) changes along the curve through the hit, and `Tangent` when it does not. So a hit of odd multiplicity along the curve is `Crossing` and one of even multiplicity is `Tangent`.

An exact tangent contact with proven branch membership records `Tangent(Curve(c), Surface(s))` at its new point, with origin `Construction("intersect.curve-surface")`. Bounded contact guesses and uncertain merges never become facts. On a later query, the same curve, surface and exact point definition settle contact with one `note_fact_decision` and no arithmetic. Other hits of that pair still get their own contact decisions.

**Order.** Hits are listed along the curve's natural direction. On a curve whose form is a line, in increasing `(p - origin) . dir` of that form. On a closed curve, starting from its parameter origin (piece 2's `parameter_origin`: the `t = 0` point of a circle or ellipse, the witness of an intersection curve); a hit at the parameter origin comes first. On an intersection curve that ends at a singular point, the singular point is passed once, between the two ends of the arc.

**Points.** Each hit is a new point `PointDef::Meet(Meet { of: MeetOf::CurveSurface(c, s), position })`. Piece 2 records its incidences on `c` and `s`. The position is `Position::Closed` with exact `Real` coordinates, or `Position::Isolated` with a root of the three-variable system made of the curve's equations (piece 2's `Curve::equations`) and `f_s`, certified by `Decider::certify_root`. A `Tangent` hit is always `Closed`, because `certify_root` refuses a tangent root (piece 1), except a hit merged on a `NotExact` (below). A hit on an analytic line is always `Closed`: it is a root of a polynomial of degree 2 at most in the line's parameter, with rational coefficients. A hit kept on the branch by a `NotExact` (below) is added as `MeetOf::Surfaces(vec![s0, s1, s])`, the curve's two surfaces and `s`, instead of `MeetOf::CurveSurface`, so that no fact puts it on a branch it may not be on (DESIGN.md, R5); it is a root of the system, so it lies on all three surfaces.

**Exactness.** For a line, circle or ellipse with rational data and a supported surface, the hits are the real roots of polynomials with rational coefficients, found by `Decider::real_roots` and other exact methods, so the answer is exact under the default budget. For an intersection curve, it also needs exact decisions about which branch a solution of the system lies on; they are signs of the frame quantities of the tables (for example the sign of `(p - F) . w`), decided by `Decider::sign`.

**When a decision ends `NotExact`.** The answer takes the degenerate branch and adds a mark built with `Mark::from_not_exact`:

| Undecided question | Degenerate branch | Mark kind and bound |
| --- | --- | --- |
| whether two candidate hits are the same point | one hit, merged as below | `Distance`: a proven upper bound on the distance between the two candidates |
| whether a candidate is on the curve's branch | on the branch, added as a point of the three surfaces (above) | `Distance`: a proven upper bound on the candidate's distance to the dividing plane of the branch (listed under `on_curve`) |
| whether a hit is crossing or tangent | `Tangent` | `Value`: the undecided quantity's bound |
| the order of two hits along the curve | one hit, merged as below | `Distance`: a proven upper bound on the distance between the two hits |

A merged hit is one new point. Its position is that of one of the two candidates, chosen deterministically by the creator, and may be `Isolated`; its contact is `Tangent`. Tests check it with `same_point` against both expected positions and accept either.

No tie is ever broken by cached coordinates.

A refusal of piece 1 (`numbers.deadline`, `numbers.root-not-certified`, `numbers.budget-exhausted` from `real_roots`) is returned unchanged, and nothing is added to the store on a refusal.

This function records no fact of its own.

## Questions about stored curves

These take `&GeomStore` and never add to it or record facts. A caller that wants to keep a proven zero records it with `GeomStore::record_fact` and `FactOrigin::Proof(proof)`. Tests check the outcome and, for a `NotExact`, its bound kind; they do not check which proof kind a `Proven` or `ProvenCoincident` answer carries.

**`on_curve(geom, c, p, dec) -> Outcome<Apart>`**: whether the point lies on the closure of the curve itself. For a line, circle or ellipse this is piece 2's `point_on_curve`. For an intersection curve it is `point_on_curve` (both surfaces) and then the branch: the point must be on the component the witness is on. That holds when the point and the witness give the same sign of the dividing quantity of the configuration, or the point's quantity is zero, decided with `Decider::sign`:

| Configuration | Dividing quantity |
| --- | --- |
| PC3 | `(p - F) . v` |
| CP3 | `(p - q_a - t u) . v` |
| CS1 | `(p - c) . w` |
| CC4, CC5 | `(p - f_B) . e1` |
| CS7, CS8 | `(p - F) . w` |

The quantity is zero only at a singular point, which lies on both branches. The configurations with one curve component (CC3, CS6) need no branch decision. It is settled by fact when `holds(Incident, Point(p), Curve(c))` is true, with one `note_fact_decision`. A `NotExact` carries a `Distance` bound: piece 2's `point_on_curve` bound when a surface was undecided, and otherwise a proven upper bound on the point's distance to the dividing plane (piece 1's bound on the quantity divided by a rational lower bound on the length of its vector). The surfaces are decided first.

**`same_curve(geom, a, b, dec) -> Outcome<Apart>`**: whether the closures of the two curves are the same point set. Settled by fact when `holds(Coincident, Curve(a), Curve(b))` is true, which includes `a == b`. Two analytic curves are compared with exact rational arithmetic, like piece 2's `same_surface`: two lines when their directions are parallel and one origin lies on the other line; two circles when the centres, the squared radii and the normal lines agree; two ellipses when the centres, the normal lines, `major_sq`, `minor_sq` and the major axis lines agree. A circle is never the same as an ellipse, and kinds of `curve_form` that differ are never the same. Otherwise it compares the closed forms with decisions, or for two branches, checks that both lie on the same two surfaces (up to `same_surface`) and that one witness is `on_curve` of the other. A `NotExact` carries the `Distance` bound of the first comparison that ended `NotExact`: a proven upper bound on the distance between the two points it compared, or the bound of `on_curve` of one witness on the other curve.

For circular carriers, recorded `Coaxial` settles parallel normals and recorded `Concentric` settles centre equality. Each used relation counts as one fact decision without arithmetic. Neither relation settles radius equality, and shared axes alone do not settle centre equality.

**`curve_is_closed(geom, c, dec) -> Result<bool, Refusal>`**: whether a curve closes on itself. Analytic lines return false; analytic circles and ellipses return true, without a decision. For an intersection curve it validates the definition and queries its form. Line, parabola and hyperbola forms are open; the currently supported circle, ellipse and quartic branch forms are closed. Branch closure follows from the certified cylinder/sphere or cylinder/cylinder frame. Unsupported definitions or unknown closure refuse with `intersect.unsupported-curve`, never a guessed boolean. Callers use this query rather than assume that every non-line is closed.

**`order_on_curve(geom, c, p, q, dec) -> Outcome<StrictOrder>`**: on a curve whose form is a line, the sign of `(q - p) . dir`: `Proven(Less)` when `p` comes before `q` in the natural direction, `ProvenCoincident` when the two points are at the same place along the line. It is one `Decider::sign` decision, and a `NotExact` keeps piece 1's bound on `(q - p) . dir` (a `Value` bound).

**`cyclic_order(geom, c, [p, q, r], dec) -> Outcome<Cyclic>`**: on a closed curve, `Proven(Forward)` when travelling from `p` in the natural direction meets `q` before `r`, `Proven(Backward)` when it meets `r` first, and `ProvenCoincident` when two of the three points are the same point. On a circle or an ellipse it is the sign of `((q - p) x (r - p)) . normal`, one decision, and a `NotExact` keeps piece 1's bound on that quantity (a `Value` bound). On a branch it is decided from the frame of the tables without angles (DESIGN.md, R8), and a `NotExact` carries the `Value` bound of the first undecided quantity. The three points must lie on the curve.

**`tangent_at(geom, c, p, dec) -> Result<[Vector3; 2], Refusal>`**: exact tangent vectors, unnormalised, in the natural direction: `[arriving, leaving]`. For an analytic curve no decision is made and both entries are the tangent: `dir` for a line; `normal x (p - center)` for a circle; for an ellipse, `normal x g`, with `d = p - center`, `V = normal x major_dir` and `g = ((d . major_dir) / (|major_dir|^2 major_sq)) major_dir + ((d . V) / (|V|^2 minor_sq)) V`. For an intersection curve it first decides whether `p` is a singular point, by deciding with `Decider::sign` whether each component of `grad f0 x grad f1` at `p` is zero; nothing in the store marks a point as singular. At a regular point both entries are `grad f0 x grad f1`. At a singular point, which is always an end of the curve's arc, they are the one-sided tangents with which the curve arrives at and leaves the point, each a positive multiple of the true one-sided tangent. A decision that ends `NotExact` makes it refuse with that code. Tests compare the tangents of lines, circles and ellipses exactly and those of intersection curves as positive multiples.

**`arc_point(geom, c, from, to, dec) -> Answered<Point3>`**: an exact point strictly inside an arc. On a curve whose form is a line, `from` and `to` are both `Some`, with `from` before `to`, and the point lies strictly between them. On a closed curve, `from` and `to` both `None` mean the whole curve and the point is any point other than the parameter origin; `Some(p), Some(p)` means the whole curve starting at `p` and the point is not `p`; `Some(p), Some(q)` with `p != q` mean the arc from `p` to `q` in the natural direction and the point lies strictly inside it. Which point is chosen is the creator's choice; it is deterministic. `arc_point` returns no marks: a decision that ends `NotExact` makes it refuse with that code.

## Traces

`trace(geom, c, from, to, dec) -> Answered<CurveTrace>` gives the arc with the same meaning of `from` and `to` as `arc_point`, as a list of `TracePiece` values in order along the natural direction.

Each piece is a parametrised curve `C(t)` on a rational range `t0 < t1`:

1. `C` is real-analytic on an open interval containing `[t0, t1]`, one-to-one there, and its derivative is never zero there;
2. `C` increases along the natural direction of the curve;
3. the first piece starts at the start of the arc (`from`, or the parameter origin for a whole curve given by `None`), the last ends at its end, and each piece ends at the point where the next starts;
4. the union of the pieces' images is exactly the arc, and two pieces share only their joining points.

How the arc is cut into pieces, and which parametrisation each piece uses, is the creator's choice. Parameter values, angles included, are caches here; they never enter a decision. An intersection curve that ends at a singular point is still regular along each of its branches, so its pieces satisfy rule 1 up to the singular end.

A `TracePiece` encloses:

- `range()`: `(t0, t1)`;
- `evaluate(t)`: three balls containing `C(s)` for every `s` in the ball `t`;
- `derivative(t)`: three balls containing `C'(s)` for every `s` in `t`;
- `taylor(t, order)`: `order + 1` entries; entry `k` holds three balls containing `C^(k)(s) / k!` for every `s` in `t`.

The result's precision is that of `t`, as for piece 2's `evaluate`. These enclosures are what measurement integrates and what export fits; both need proven bounds, not samples, and high Taylor orders make validated quadrature cheap. A ball `t` outside `[t0, t1]` is allowed as long as it lies in the open interval of rule 1; outside it the balls may be of any width, but never a false enclosure. Tests evaluate only inside `[t0, t1]`.

A trace piece encloses its non-rational constants when a call needs them, at that call's precision (through `Decider::enclose_at`), and keeps each enclosure for later calls at that precision or lower. `trace` itself encloses no constant to a precision no call has asked for, so evaluating at `p` bits costs work for about `p` bits, whatever `max_precision_bits` is.

`trace` makes decisions to locate `from` and `to` on the curve. When one ends `NotExact`, or a refusal comes from piece 1, `trace` refuses with that code; a trace of the wrong arc would make measurement wrong. It returns no marks.

## Decisions, outcomes and marks

- Every numeric decision goes through the `Decider` or is exact rational arithmetic. No decision reads a cached `f64` coordinate.
- A `Proven` or `ProvenCoincident` outcome is never wrong (DESIGN.md, R2). An `intersect_surfaces` answer is always exact.
- Repeating any call on equal stores with a `Decider` in the same state gives equal results, the same new handles in the same order, and the same statistics, on every machine and thread. The order of new handles and of recorded facts within a call is not otherwise part of the contract; tests find new handles through the answer, never by index.

## How this piece uses pieces 1 and 2

- **Outcomes.** Operations that add to the store return `Answered`; questions return `Outcome`, mapped from piece 1's and piece 2's as stated above. A `NotExact` from piece 1 or 2 inside an operation becomes the degenerate branch and a mark; inside a question it is returned as `NotExact` with the bound stated.
- **Budgets.** Every decision runs on the caller's `Decider` and its `BudgetPolicy`. This piece never creates a `Decider` or a policy. The deadline reaches it only through `Decider` calls that do work (ARCHITECTURE.md, "The time budget"): `intersect_surfaces` makes none, and a fact-settled decision makes none.
- **Facts.** This piece records the construction facts listed above and uses the store's structural identities. Fact-settled incidence, coincidence, shared-axis, shared-centre and local tangency decisions call `note_fact_decision` without numeric work.
- **Marks.** `intersect_curve_surface` builds marks with `Mark::from_not_exact`, with the kinds of its table. `intersect_surfaces` and `trace` never build one.
- **Refusals.** This piece adds the codes below and passes every `numbers.*` refusal through unchanged.

## Refusal codes

| Code | When |
| --- | --- |
| `intersect.unsupported-pair` | the CC6 pair whose axis lengths have an irrational ratio, or an unknown surface kind |
| `intersect.unsupported-curve` | an unsupported manually defined intersection, including the CC6 definition described above |
| `intersect.traced-pair-required` | a valid cone/torus pose needs the general certified branch tracer |
| `intersect.algebraic-direction-required` | a valid closed form needs an irrational direction that the current rational frame cannot store |

`numbers.*` refusals pass through. This piece adds no `NotExact` code.

## Programming errors

These may panic. None can come from a question file:

- a handle from another store, or out of range;
- `order_on_curve` on a curve whose form is not a line, or `cyclic_order` on a line;
- `arc_point` or `trace` on a curve whose form is a line, with `from` or `to` missing, or with `from` not before `to`;
- points that do not lie on the curve, passed to `order_on_curve`, `cyclic_order`, `tangent_at`, `arc_point` or `trace`;
- an intersection curve whose witness is not a regular point of its surfaces (piece 2's rule).

## Exact public interface

The crate re-exports these types at its root. Everything listed is `pub`; fields not listed are private. The compiled crate in `crates/catallax-intersect` is authoritative if a signature listing becomes stale; update this contract with any API change.

```rust
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Contact { Crossing, Tangent }
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Component {
    Curve { curve: CurveId, contact: Contact, singular: Vec<PointId> },
    Point { point: PointId },
}
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum SurfacePair { Disjoint, Coincident { same_sense: bool }, Meets(Vec<Component>) }
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CurveHit { pub point: PointId, pub contact: Contact }
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum CurveSurface { Disjoint, Contained, Meets(Vec<CurveHit>) }
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Cyclic { Forward, Backward }
#[non_exhaustive]
#[derive(Clone, Debug)]
pub enum CurveForm {
    Line { origin: Point3, dir: RVector3 },
    Circle { center: Point3, normal: RVector3, radius_sq: Real, ref_dir: RVector3 },
    Ellipse { center: Point3, normal: RVector3, major_dir: RVector3, major_sq: Real, minor_sq: Real },
    Parabola { vertex: Point3, normal: RVector3, axis_dir: RVector3, focal: Real },
    Hyperbola { center: Point3, normal: RVector3, axis_dir: RVector3, transverse_sq: Real, conjugate_sq: Real, branch: StrictSign },
    Branch,
}

#[derive(Clone, Debug)]
pub struct TracePiece;
impl TracePiece {
    pub fn range(&self) -> (Rational, Rational);
    pub fn evaluate(&self, t: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal>;
    pub fn derivative(&self, t: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal>;
    pub fn taylor(&self, t: &Ball, order: usize, dec: &Decider) -> Result<Vec<[Ball; 3]>, Refusal>;
}
// Enclosures check the caller's deadline during series evaluation and subdivision.
#[derive(Clone, Debug)]
pub struct CurveTrace;
impl CurveTrace { pub fn pieces(&self) -> &[TracePiece]; }

pub fn supports(a: SurfaceKind, b: SurfaceKind) -> bool;
pub fn intersect_surfaces(geom: &mut GeomStore, a: SurfaceId, b: SurfaceId, dec: &mut Decider) -> Answered<SurfacePair>;
pub fn intersect_curve_surface(geom: &mut GeomStore, c: CurveId, s: SurfaceId, dec: &mut Decider) -> Answered<CurveSurface>;
pub fn curve_form(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<CurveForm, Refusal>;
pub fn curve_is_closed(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<bool, Refusal>;
pub fn on_curve(geom: &GeomStore, c: CurveId, p: PointId, dec: &mut Decider) -> Outcome<Apart>;
pub fn same_curve(geom: &GeomStore, a: CurveId, b: CurveId, dec: &mut Decider) -> Outcome<Apart>;
pub fn order_on_curve(geom: &GeomStore, c: CurveId, p: PointId, q: PointId, dec: &mut Decider) -> Outcome<StrictOrder>;
pub fn cyclic_order(geom: &GeomStore, c: CurveId, points: [PointId; 3], dec: &mut Decider) -> Outcome<Cyclic>;
pub fn tangent_at(geom: &GeomStore, c: CurveId, p: PointId, dec: &mut Decider) -> Result<[Vector3; 2], Refusal>;
pub fn arc_point(geom: &GeomStore, c: CurveId, from: Option<PointId>, to: Option<PointId>, dec: &mut Decider) -> Answered<Point3>;
pub fn trace(geom: &GeomStore, c: CurveId, from: Option<PointId>, to: Option<PointId>, dec: &mut Decider) -> Answered<CurveTrace>;
```

## What defines correct behaviour

The test writer builds tests from the cases below. Every case has an answer that maths gives, so tests compare against that answer, never against a run of the code.
Independent oracles are allowed: a third-party big-rational crate used only as a dev-dependency, or values generated by a separate program, as long as nothing of this piece's own code computes the expected answer.

Tests compare an added analytic curve field by field with `==` against the table's formula evaluated in exact rationals. They compare points with piece 2's `same_point` against the expected exact position, built as a `Position::Closed` point, and check membership with `side`. Where a case gives a curve "through" a point, the test checks `on_curve` of that point.

Unless stated, surfaces use these definitions. `ref_dir` is `(1, 0, 0)` for planes with normal along `z` and for cylinders with axis along `z`; `(0, 1, 0)` for planes with normal along `x`; `(0, 0, 1)` for planes with normal along `y` and for cylinders with axis along `x`; `(1, 0, 0)` for the plane of case 4 and the cylinders with axis `(0, 1, 1)` of case 30; `(0, 1, 0)` for the cylinder with axis `(3, 0, 4)` of case 22; `(0, 0, 1)` for the cylinder with axis `(1, 1, 0)` of case 22. Every sphere has axis `(0, 0, 1)` and `ref_dir (1, 0, 0)`. Where these rules and a case disagree, the case's `ref_dir` must be perpendicular to its axis or normal, and the expected answers do not depend on it unless the case names it.

- `Z5`: the cylinder with origin `0`, axis `(0, 0, 1)`, radius 5;
- `S5`: the sphere with centre `0`, radius 5, axis `(0, 0, 1)`;
- `plane x = k`: origin `(k, 0, 0)`, normal `(1, 0, 0)`; similarly for `y` and `z`.

### Plane and plane

1. The plane with origin `(0, 0, 1)` and normal `(0, 0, 2)` and the plane `x = 3`: one line, `origin (3, 0, 1)`, `dir (0, 2, 0)`, `Crossing`. Swapped: `origin (3, 0, 1)`, `dir (0, -2, 0)`. Facts `Incident(line, a)` and `Incident(line, b)` hold.
2. `z = 0` and the plane with origin `(0, 0, 10^-30)` and normal `(0, 0, -3)`: `Disjoint`. With origin `(5, 7, 0)` and normal `(0, 0, -3)`: `Coincident { same_sense: false }`; with normal `(0, 0, 4)`: `Coincident { same_sense: true }`. A plane with itself (the same handle): `Coincident { same_sense: true }`.

### Plane and cylinder

3. `z = 2` (normal `(0, 0, 1)`) and `Z5`: one circle, `center (0, 0, 2)`, `normal (0, 0, 1)`, `radius_sq 25`, `ref_dir (1, 0, 0)`. Swapped: `normal (0, 0, -1)`.
4. The plane through `0` with normal `(0, 3, 4)` and `Z5`: one ellipse, `center 0`, `normal (0, 3, 4)`, `major_dir (0, -12/25, 9/25)`, `major_sq 625/16`, `minor_sq 25`. Its major vertex `(0, -5, 15/4)` lies on both surfaces.
5. `x = 3` and `Z5`: two lines (`h2 = 16`, rational). First `origin (3, -4, 0)`, `dir (0, 0, -1)`; then `origin (3, 4, 0)`, `dir (0, 0, 1)`. Swapped: the same order with `dir (0, 0, 1)` then `(0, 0, -1)`.
6. `x = 1` and `Z5`: `h2 = 24` is not a rational square, so two intersection curves on `[plane, Z5]` with witnesses at `(1, -2 sqrt(6), 0)` and `(1, 2 sqrt(6), 0)`, in that order. `curve_form` gives `Line { origin: witness, dir: (0, 0, -1) }` and `Line { origin: witness, dir: (0, 0, 1) }`. The call adds 2 curves and 2 points.
7. `x = 5` and `Z5`: one line `origin (5, 0, 0)`, `dir (0, 0, 1)`, `Tangent`, and the fact `Tangent(plane, Z5)` at that line. `x = 5 + 10^-30`: `Disjoint`. `x = 5 - 10^-30`: two intersection curves (`h2 = 10^-29 - 10^-60` is not a square).
8. Along the seam: `y = 0` and `Z5`: two lines (`h2 = 25`), first `origin (5, 0, 0)` (the seam line of `Z5`), `dir (0, 0, -1)`; then `origin (-5, 0, 0)`, `dir (0, 0, 1)`.

### Plane and sphere

9. `z = 3` and `S5`: one circle, `center (0, 0, 3)`, `normal (0, 0, 1)`, `radius_sq 16`, `ref_dir (1, 0, 0)`.
10. `z = 5` and `S5`: one point at `(0, 0, 5)`, the north pole, and the fact `Tangent(plane, S5)` at that point. `z = 5 + 10^-30`: `Disjoint`. `z = 5 - 10^-30`: one circle with `radius_sq = 10^-29 - 10^-60`.
11. `y = 0` and `S5`: the great circle through the seam, `center 0`, `normal (0, 1, 0)`, `radius_sq 25`, `ref_dir (0, 0, 1)`.

### Sphere and sphere

12. `S5` and the sphere with centre `(6, 0, 0)` and radius 5: one circle, `center (3, 0, 0)`, `normal (6, 0, 0)`, `radius_sq 16`, `ref_dir (0, 0, 6)` (`perp((6, 0, 0)) = (6, 0, 0) x (0, 1, 0)`). Its parameter origin is `(3, 0, 4)`.
13. `S5` and centre `(10, 0, 0)`, radius 5: one point `(5, 0, 0)`. `S5` and centre `(2, 0, 0)`, radius 3: one point `(5, 0, 0)` (inside tangent).
14. `S5` with centre `(10 + 10^-30, 0, 0)`, radius 5: `Disjoint`. With centre `(1, 0, 0)`, radius 3: `Disjoint` (nested). With centre `0`, radius 4: `Disjoint`. With itself as another handle with equal fields: `Coincident { same_sense: true }`.

### Cylinder and cylinder

15. `Z5` and the cylinder with origin `(8, 0, 0)`, axis `(0, 0, -2)`, radius 5: two lines (`h2 = 9/64`). First `origin (4, 3, 0)`, `dir (0, 0, 1)`; then `origin (4, -3, 0)`, `dir (0, 0, -1)`.
16. `Z5` and the cylinder with origin `(10, 0, 0)`, axis `(0, 0, 1)`, radius 5: one line `origin (5, 0, 0)`, `dir (0, 0, 1)`, `Tangent`. `Z5` and the cylinder with origin `(2, 0, 0)`, axis `(0, 0, 1)`, radius 3: one line `origin (5, 0, 0)`, `dir (0, 0, 1)`, `Tangent` (inside). `Z5` and origin `0`, axis `(0, 0, 3)`, radius 2: `Disjoint`; radius 5: `Coincident { same_sense: true }`.

In cases 17 to 21, `A2` is the cylinder with origin `0`, axis `(0, 0, 1)`, radius 2, and `X1(y0)` is the cylinder with origin `(0, y0, 0)`, axis `(1, 0, 0)`, radius 1. In cases 17 to 20, `B = A2`, `m = (0, 1, 0)`, `e1 = (1, 0, 0)`, `e2 = (0, 1, 0)`, `e3 = (0, 0, 1)`, `alpha = 1` and `gamma = 0`.

17. `A2` and `X1(0)` (CC5, `delta = 0`): two intersection curves with witnesses `(2, 0, 1)` then `(-2, 0, 1)`. At `(2, 0, 1)`, `tangent_at` gives a positive multiple of `(0, -1, 0)` (`grad f_A x grad f_X = (4, 0, 0) x (0, 0, 2) = (0, -8, 0)`).
18. `A2` and `X1(2)` (CC3): one intersection curve with witness `(0, 2, 1)`.
19. `A2` and `X1(1)` (CC4): two intersection curves with witnesses `(sqrt(3), 1, 1)` then `(-sqrt(3), 1, 1)`, each with `singular = [(0, 2, 0)]`, the same `PointId` in both. The call adds 2 curves and 3 points.
20. `A2` and `X1(3)` (CC2): one point `(0, 2, 0)`. `X1(3 + 10^-30)`: `Disjoint`. `X1(3 - 10^-30)`: one intersection curve (CC3).
21. Steinmetz: `A1` (origin `0`, axis `(0, 0, 1)`, radius 1) and `X1(0)`: two ellipses. First `center 0`, `normal (-1, 0, 1)`, `major_dir (1, 0, 1)`, `major_sq 2`, `minor_sq 1`; then `normal (1, 0, 1)`, `major_dir (-1, 0, 1)`, `major_sq 2`, `minor_sq 1`. Both list `singular = [(0, 1, 0), (0, -1, 0)]`. The facts `Incident(p, ellipse)` hold for both points and both ellipses.
22. Steinmetz at another angle: `A1` and the cylinder with origin `0`, axis `(3, 0, 4)`, radius 1 (`k = 5`, `cos = 4/5`): first `normal (-15, 0, 5)`, `major_dir (15, 0, 45)`, `major_sq 10`, `minor_sq 1`; then `normal (15, 0, 45)`, `major_dir (-15, 0, 5)`, `major_sq 10/9`, `minor_sq 1`. With axis `(1, 1, 0)` instead (`k^2 = 2`): refused with `intersect.unsupported-pair`.

### Cylinder and sphere

23. `Z5` and the sphere with centre `(0, 0, 1)`, radius 13 (`k2 = 144`): two circles, first `center (0, 0, 13)`, `normal (0, 0, -1)`; then `center (0, 0, -11)`, `normal (0, 0, 1)`; both `radius_sq 25`, `ref_dir (1, 0, 0)`. With radius 6 (`k2 = 11`): two intersection curves with witnesses `(5, 0, 1 + sqrt(11))` then `(5, 0, 1 - sqrt(11))`, whose `curve_form` is `Circle { center: (0, 0, 1 + sqrt(11)), normal: (0, 0, -1), radius_sq: 25, ref_dir: (1, 0, 0) }` and the same with `1 - sqrt(11)` and normal `(0, 0, 1)`. With radius 5: one circle `center (0, 0, 1)`, `normal (0, 0, 1)`, `Tangent`. With radius 4: `Disjoint`.

In cases 24 to 27, `Z1` is the cylinder with origin `0`, axis `(0, 0, 1)`, radius 1, and every sphere centre is on the positive `x` axis, so `F = 0` and `e1 = (1, 0, 0)`.

24. Viviani's curve: `Z1` and the sphere with centre `(1, 0, 0)` and radius 2 (CS7, `d = 1`): two intersection curves with witnesses `(1, 0, 2)` then `(1, 0, -2)`, each with `singular = [(-1, 0, 0)]`. On the first, `tangent_at` at the witness is a positive multiple of `(0, -1, 0)`, and at the node it is `[arriving, leaving]` = positive multiples of `[(0, 1, -1), (0, 1, 1)]`.
25. `Z1` and the sphere with centre `(2, 0, 0)` and radius 2 (CS6): one intersection curve with witness `(1, 0, sqrt(3))`.
26. `Z1` and the sphere with centre `(1, 0, 0)` and radius 3 (CS8): witnesses `(1, 0, 3)` then `(1, 0, -3)`.
27. `Z1` and the sphere with centre `(3, 0, 0)` and radius 2 (CS5, outside): one point `(1, 0, 0)`. The cylinder with origin `0`, axis `(0, 0, 1)`, radius 3, and the sphere with centre `(1, 0, 0)`, radius 2 (CS5, inside): one point `(3, 0, 0)`. Near misses: in the first pair, moving the centre to `(3 + 10^-30, 0, 0)` gives `Disjoint` (CS4) and to `(3 - 10^-30, 0, 0)` one branch (CS6); in the second, moving it to `(1 + 10^-30, 0, 0)` gives one branch (CS6) and to `(1 - 10^-30, 0, 0)` `Disjoint` (CS4).

### Curves and surfaces

28. The line with origin `(-10, 0, 0)` and `dir (1, 0, 0)`: with `S5`, hits `(-5, 0, 0)` then `(5, 0, 0)`, both `Crossing`; with `Z5`, the same; with `x = 3`, one `Crossing` hit `(3, 0, 0)`; with `y = 0`, `Contained`. The line with origin `(-10, 5, 0)` and `dir (1, 0, 0)`: with `Z5` and with `S5`, one `Tangent` hit `(0, 5, 0)`. Moved to origin `(-10, 5 + 10^-30, 0)`: `Disjoint` for both.
29. The circle with centre `0`, normal `(0, 0, 1)`, `radius_sq 25`, `ref_dir (1, 0, 0)` (call it `K5`): with `x = 3`, hits `(3, 4, 0)` then `(3, -4, 0)`; with `x = 5`, one `Tangent` hit `(5, 0, 0)`, the parameter origin; with the sphere of centre `(5, 0, 0)` and radius 5, two `Crossing` hits `(5/2, 5 sqrt(3)/2, 0)` then `(5/2, -5 sqrt(3)/2, 0)`; with the cylinder of origin `0`, axis `(1, 0, 0)`, radius 3, four `Crossing` hits `(4, 3, 0)`, `(-4, 3, 0)`, `(-4, -3, 0)`, `(4, -3, 0)` in that order; with `z = 0`, `Contained`; with `Z5`, `Contained` (found by numbers, since no fact links them).
30. A tilted cylinder: `K5` with the cylinder of origin `(4, 0, 0)`, axis `(0, 1, 1)`, radius 2. In the plane `z = 0` that cylinder's implicit is `y^2 + 2 (x - 4)^2 - 8`, which meets `x^2 + y^2 = 25` where `x^2 - 16x + 49 = 0`. The root with `|x| <= 5` is `x = 8 - sqrt(15)`, so there are two `Crossing` hits `(8 - sqrt(15), ±sqrt(25 - (8 - sqrt(15))^2), 0)`, the one with `y > 0` first. With radius 2 and origin `(3, 0, 0)` instead, the equation is `x^2 - 12x + 35 = 0`, its root `x = 5` gives one `Tangent` hit at `(5, 0, 0)`, the parameter origin.
31. The ellipse of case 4 with `x = 0`: two `Crossing` hits at its major vertices, `(0, -5, 15/4)` (the parameter origin, so first) then `(0, 5, -15/4)`. With `y = 0`: on the ellipse `y = 0` forces `z = 0` and `x^2 = 25`, so two `Crossing` hits, `(5, 0, 0)` then `(-5, 0, 0)` (counter-clockwise about `(0, 3, 4)` from the parameter origin).
32. Branch and plane: the first curve of case 24 (the upper branch of Viviani's curve) with `z = 1`: two `Crossing` hits, `(-1/2, -sqrt(3)/2, 1)` then `(-1/2, sqrt(3)/2, 1)`. With `z = 0`: one `Tangent` hit at the node `(-1, 0, 0)`. With `z = 3`: `Disjoint`. The second curve of case 24 with `z = 1`: `Disjoint`.
33. Branch and its own surfaces: every intersection curve of cases 17 to 26 with either of its surfaces is `Contained`, settled by fact with no work.

### Questions about stored curves

34. `on_curve`: `(3, 4, 0)` on `K5` is `ProvenCoincident`; `(3, 4, 10^-30)` is `Apart`. On the first curve of case 26 (upper branch, witness `(1, 0, 3)`), the point `(1, 0, -3)` lies on both surfaces but on the other branch: `Apart`. The node of case 24 is `ProvenCoincident` on both of its curves.
35. `same_curve`: the circle of case 3 and the circle with the same fields but `normal (0, 0, -5)` and `ref_dir (0, 1, 0)` are `ProvenCoincident`; with `radius_sq 25 + 10^-30`, `Apart`. The two curves of case 6 are `Apart`. For the first intersection curve of case 23 with radius 6, the test adds by hand a second intersection curve on the same two surfaces whose witness is another point of the same circle, `(-5, 0, 1 + sqrt(11))` as a closed meeting point: the two are `ProvenCoincident`. With the witness `(-5, 0, 1 - sqrt(11))` they are `Apart`.
36. `order_on_curve` on the line of case 28: `(-5, 0, 0)` before `(5, 0, 0)` is `Proven(Less)`; the reverse is `Proven(Greater)`; a point with itself is `ProvenCoincident`.
37. `cyclic_order` on `K5`: `[(5, 0, 0), (0, 5, 0), (-5, 0, 0)]` is `Forward`, `[(5, 0, 0), (-5, 0, 0), (0, 5, 0)]` is `Backward`, and two equal points give `ProvenCoincident`. On the upper branch of case 24, from the witness `(1, 0, 2)`, the hits of case 32 come in the order listed: `[(1, 0, 2), (-1/2, -sqrt(3)/2, 1), (-1/2, sqrt(3)/2, 1)]` is `Forward`.
38. `tangent_at` on `K5` at `(3, 4, 0)`: both entries are `(-4, 3, 0)`. On the line of case 28: `(1, 0, 0)`.
39. `arc_point`. The test never adds the returned `Point3` to the store, since a meeting point would record facts that make the check pass by default. It checks `p` with piece 2's `side_of(s, &p, dec)`, which must be `ProvenCoincident` on each surface the curve lies on, and decides the order with `Decider::sign` on quantities it builds with `Point3` and `Vector3` arithmetic. On `K5` from `(5, 0, 0)` to `(0, 5, 0)`: `side_of` coincident on `z = 0` and `Z5`, and `((p - (5, 0, 0)) x ((0, 5, 0) - (5, 0, 0))) . (0, 0, 1) > 0`. On the line of case 28 between its two `S5` hits: coincident on `y = 0` and `z = 0`, and `p_x + 5 > 0`, `5 - p_x > 0`. On the whole upper branch of case 24 (`from` and `to` both `None`): coincident on `Z1` and the sphere, `p_z >= 0` (the dividing quantity), and `p` differs from the witness `(1, 0, 2)` in some coordinate.

### Traces

The test checks a trace with its own validated arithmetic: each piece's `evaluate` at its range ends encloses the expected end points; `evaluate` of random sub-balls gives balls in which each surface's implicit polynomial, evaluated in `Ball`, contains 0; `derivative` has a positive dot product with `tangent_at` at the ends; and `taylor` agrees with `evaluate` and `derivative` on its first two entries.

40. The whole of `K5`: the arc length `10 pi`, computed by the test from the Taylor enclosures with any validated quadrature, encloses `10 pi`. The arc of `K5` from `(5, 0, 0)` to `(0, 5, 0)`: length `5 pi / 2`.
41. The segment of the line of case 28 between its `S5` hits: length 10.
42. The whole of the first curve of case 25: its length, an elliptic integral, contains the value `L25` given under "Values checked for this description".
43. The upper branch of case 24, which ends at the node, traced with `from` and `to` both the node `(-1, 0, 0)` (the whole curve starting there): the first piece starts and the last ends at the node, the pieces pass through `(1, 0, 2)`, and the length contains `L24`. Traced with both `None`, it starts and ends at the witness `(1, 0, 2)` and has the same length.

### Scale, placement and determinism

44. Scale (DESIGN.md, R9): cases 1 to 33 with every point, origin, radius and length multiplied by `10^k` for `k` in `-30..=30` (squared radii by `10^(2k)`, directions unchanged) keep their configurations, their component and hit counts, contacts and orders, and every expected position scales by `10^k`. Under the default budget each scaled answer carries no mark, as the unscaled answer does; a mark there is a failure, not a permitted outcome.
45. Placement: the same cases moved by 20 random rational poses (integer quaternion components up to 100, translations up to `10^6`) keep their configurations, counts, contacts and orders, with no mark under the default budget; every expected field and position is the moved one (directions moved by `apply_vector`), except the `ref_dir` of an SS3 circle, which is `perp` of the moved `u`, because `perp` picks a coordinate axis; the parameter origins and hit orders on those circles follow it. A generated family repeats each configuration of the tables at random rational parameters chosen on and near every boundary, posed by random poses, with the expected configuration computed by the test's own rational oracle from the table's conditions.
46. Determinism: repeating any call on a clone of the store with a fresh `Decider` gives equal results and statistics. Results compare with `==` where the type has it and otherwise by their `Debug` text, and so do the definitions the call added. Eight threads, each with its own clone and `Decider`, give identical results.
47. Budgets: under `BudgetPolicy { max_precision_bits: 64, max_work: 10_000 }` (this budget admits piece 1's `f64` filter, which costs one work unit; with `max_work = 0` piece 1 decides nothing, not even by the filter), `intersect_surfaces` gives its default-budget answer (it makes no decision), and cases 28 to 39 give their default answer, a marked answer whose marks have the table's kinds (for a question, a `NotExact` with the stated bound kind), or a refusal; never a wrong proven answer. A `Decider` whose deadline has passed: `intersect_surfaces` and every fact-settled answer (the `Contained` of case 33, `on_curve` of the node in case 34, `same_curve` of a curve with itself) still answer, since they make no `Decider` call that does work; every other call either refuses with `numbers.deadline` or gives its default-budget answer (ARCHITECTURE.md, "The time budget"), and never any other answer.

### Volume and running time of the tests

The generated families (cases 44, 45 and 46) must together run at least 30,000 maths-answered questions per test run, from fixed seeds.
Keep the default test run under about 10 minutes in release mode on one core. `intersect_surfaces` is exact rational arithmetic and cheap. Decisions on square roots took about 0.3 to 0.9 ms each in piece 1's build, and `real_roots` of degree 4 is a few milliseconds, so run a few thousand curve-surface questions at most and the large sweeps on `intersect_surfaces`.

## Values checked for this description

The architect checked these with exact rational arithmetic and 60-digit decimals in an independent program (`sympy`):

- every witness, node, tangent point and analytic curve of cases 1 to 27 lies on both surfaces, every witness is a regular point, every node is singular, and every stated natural direction agrees with `grad f_a x grad f_b`;
- the general formulas of the tables on 40 random rational poses per configuration, including non-unit axes;
- the hit positions and orders of cases 28 to 32 and the one-sided tangents of case 24;
- the lengths `L24` and `L25` below, by direct numerical integration at 50 digits.

- `L24 = 7.64039557805542403580952416434288658` (the upper branch of Viviani's curve, the closed arc from the node through `(1, 0, 2)` back to the node: the integral of `sqrt(1 + sin(t/2)^2)` over `t` in `[-pi, pi]`);
- `L25 = 9.53883361768658472990916036073734440` (the loop of case 25: twice the integral of `sqrt(1 + 4 sin(t)^2 / (4 cos(t) - 1))` over `t` in `[-acos(1/4), acos(1/4)]`).

The test writer's own oracle restates these to the digits it needs.

## Speed: measured, never gated

No test fails for speed. A benchmark target measures and prints:

1. time of `intersect_surfaces` per pair and configuration, at identity and at random poses;
2. time of `intersect_curve_surface` per curve kind and surface kind, and the share of its decisions settled by facts and by piece 1's filter;
3. the largest precision reached in the near-miss families (`10^-30` from every tangent and coincident boundary);
4. time of `on_curve`, `same_curve`, `order_on_curve` and `cyclic_order` per curve kind;
5. time of `trace` per arc kind, and of `taylor` at orders 4 and 16 at 64 and 1,000 bits.

The numbers are written to the benchmark output so later runs can be compared.

## Out of scope for this piece

- Cones and tori (next), free-form surfaces (later).
- Faces, edges and trimming to faces (piece 5).
- Meetings of two curves other than through a surface; piece 5 finds where two curves on one face meet as the meeting of one curve with the other's second surface.
- Any tolerance. This piece has none.
