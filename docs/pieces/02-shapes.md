# Piece 2: shapes (`catallax-shapes`)

Status: complete for blind test writing and blind building.

This is the input the test writer and the creator get, together with `docs/DESIGN.md`, `docs/ARCHITECTURE.md` and piece 1's description `docs/pieces/01-numbers-and-decisions.md`. This piece is built on piece 1's numbers and outcomes and uses them exactly as that description defines them. Nothing else is needed.

## Approved algebraic-definition contract

The shared foundation implements the 2026-10-22 architecture decision below. This section supersedes the rational-only number table, pose/surface/curve signatures, equality guarantees and relevant tests below. The remaining formulas and rational examples specify the rational fast path; the old signature listing is its historical baseline, not the current API.

### Definition types and validation

| Definition | Stored type |
| --- | --- |
| Surface and analytic-curve point fields | `Point3` with `Real` coordinates |
| Directions, normals, axes, references and translations | `Vector3` with `Real` coordinates |
| Radii, squared radii, squared semi-axes and slopes | `Real` |
| Rotation matrix entries | `Real` |
| Implicit polynomials and curve equations | Piece 1 `RealPoly`, in three variables |
| Rational inputs and exact rational fast paths | Existing `Rational`, `RPoint3`, `RVector3`, `MPoly` |

`Real` is algebraic, not an arbitrary computable-real interface. Keep expressions and shared certified roots as definitions; cached coordinates or enclosures do not replace them. Directions remain unnormalised. Squared dimensions remain squared, even when they are irrational. Conversion from the rational types is exact and must not require a proof.

Analytic constructors take the target field types and a final `&mut Decider`, returning `Result<Self, Refusal>`. First validate every field's expression domain in declaration order, including position coordinates. Then preserve the invariant-check order below: directions, positive dimensions and kind-specific inequalities, then perpendicularity. Within each group, visit fields in declaration order. Success requires proofs of all invariants. `NotExact` cannot validate a nonzero axis, positivity, perpendicularity or a divisor; return its original numeric reason as a temporary failure. Proven invalid geometry uses the existing `shapes.*` code. No approximate constructor is added. Public fields, if retained during migration, carry the same caller obligation as the existing raw structs: invalid raw definitions are programming errors.

Remove `PartialEq`, `Eq` and `Hash` from definitions containing `Real`, including `Rotation`, `Pose`, analytic structs, `Surface` and `Curve`. Handles and kind tags retain their traits; rational input structs retain theirs. Geometry equality belongs to `same_surface`, `same_point` and piece 4's `same_curve`, using facts then exact decisions. An expression or construction key can prove identity for cache reuse, but unequal keys cannot prove different geometry. Do not compare independent algebraic roots by their isolating intervals or polynomial text.

`Surface::implicit()`, `GeomStore::implicit()` and `Curve::equations()` return `RealPoly` or `Vec<RealPoly>`. Preserve the exact scaling of every formula below. For rational defining data, construct the existing polynomial without avoidable square roots, so it converts back to the identical `MPoly`. An algebraic expression that happens to be rational need not be detected by the cheap conversion. Gradients and rigid-transform identities hold for all validated algebraic points. `same_surface` retains the listed geometric criteria, evaluated through the decider when quantities are non-rational. Its rational fast path retains the existing no-decision result. For elliptic cylinders, compare axis lines and the transverse quadratic forms, including their eigenvalues and eigendirections; reference directions are geometric for an ellipse. Equality of proportional implicit polynomials alone does not prove equality of selected nappes or trimmed domains.

Chart evaluation must enclose algebraic constants through piece 1 under the caller's policy. Add a final `&Decider` argument and an outer `Result<_, Refusal>` to `Surface::{evaluate, derivatives, parameters_of}` and `Curve::evaluate`; preserve each method's existing inner return type and input-precision/containment rules. A `parameters_of` seam/domain ambiguity remains `Ok(None)`, distinct from an enclosure failure. `implicit`, `gradient`, `transformed` and `parameter_origin` still only build exact expressions and need no decision context. This prevents unbudgeted coefficient refinement or hidden rationalisation in measurement, meshing and export.

### Poses and explicit points

Keep `Rotation::from_quaternion(Rational, Rational, Rational, Rational)` as the checked rational-input convenience constructor. Add `Rotation::from_real_quaternion(Real, Real, Real, Real, &mut Decider) -> Result<Rotation, Refusal>` with the same matrix formula and a proven nonzero norm. `matrix()` returns `&[[Real; 3]; 3]`. The algebraic constructor validates expression domains before division. Composition and inversion preserve certified orthogonality and determinant one by construction; do not re-prove them after every application.

`Pose::translation` uses `Vector3`. The checked constructor is `Pose::new(rotation: Rotation, translation: Vector3, dec: &mut Decider) -> Result<Pose, Refusal>`; it validates the translation expressions. `identity`, `then` and `inverse` remain infallible on validated poses. `Rotation::apply` and `Pose::apply_vector` take `&Vector3` and return `Vector3`; `Pose::apply_point` takes `&Point3` and returns `Point3`. The existing `apply_real` and `apply_real_point` may remain aliases. Rational callers convert once at this boundary. Posing an entity preserves its exact dimensions and construction provenance. Tests compare posed coordinates through exact decisions, not structural `==`.

`Point3::rational_components`, `Vector3::rational_components` and `Rotation::rational_matrix` expose borrowed rational entries when all components are cached exact rationals. They do not eliminate, refine or round algebraic expressions. `same_definition` only proves positive construction identity for reuse; a false result is not a geometric inequality proof.

A recipe angle in rational degrees uses piece 1's certified sine/cosine constructor. A rotation about an algebraic axis can use the half-angle quaternion with the exactly normalised axis. A 45-degree revolve therefore has exact end-cap and edge positions. Normalising `[1, 0, 1]` gives `[1/sqrt(2), 0, 1/sqrt(2)]`; an extrude along it retains that expression rather than refusing or rounding it.

Widen `PointDef::Explicit` to `Explicit(Point3)`. Its coordinates are exact input or construction expressions, not rounded measurements. Rational callers use `to_real()`. `Meet` continues to store the defining objects and the selected solution; do not replace a meeting definition with an explicit coordinate cache. The caller must validate explicit expressions before insertion, just as it guarantees meeting-point incidence. All explicit-point approximations become caches of the exact expressions with the same enclosure rule as other points.

### Circular cylinders

`Cylinder` stores `radius_sq: Real`, not a metric radius. `from_radius_sq` validates the positive squared radius directly. The radius-input convenience constructor `new` validates before squaring. Implicit equations, surface equality and parallel circular-carrier intersections use the stored square. `radius()` extracts the positive root only for metric evaluation or exchange. `radial_reference()` scales the seam to a point on the circle, so an extruded rational arc retains a rational parametrisation even when its radius is irrational.

Rational cylinder witnesses select their stereographic pole from the exact sign of the target's seam projection. This avoids the singular chart for rounded targets near the opposite seam point. The resulting witness still needs a distance bound and a trimmed-face membership check.

### Elliptic cylinders

Add `EllipticCylinder` to `Surface` and `SurfaceKind` with fields `origin: Point3`, `axis: Vector3`, `major_dir: Vector3`, `major_sq: Real`, `minor_sq: Real`. Its constructor follows the validation rules above, requiring nonzero perpendicular directions and `0 < minor_sq < major_sq`. Invalid squared semi-axis values use `shapes.ellipse-axes`, as for an ellipse. Zero direction vectors still use `shapes.zero-vector`. Equal axes reduce to `Cylinder` with `radius_sq = major_sq` before this constructor; negative or zero axes do not silently reduce.

Let `e3 = axis / |axis|`, `e1 = major_dir / |major_dir|`, `e2 = e3 x e1`, `d = x - origin`, `a = sqrt(major_sq)` and `b = sqrt(minor_sq)`.

- The implicit is `f(x) = (d . e1)^2 / major_sq + (d . e2)^2 / minor_sq - 1`, with outward gradient and no rescaling. Expand normalised squared projections using unnormalised dot products divided by squared direction lengths, avoiding unnecessary square roots in rational coefficients.
- `S(u,v) = origin + a cos(u) e1 + b sin(u) e2 + v e3`; domain `[Periodic, Line]`.
- `dS/du = -a sin(u) e1 + b cos(u) e2`; `dS/dv = e3`.
- `parameters_of` uses `u = angle((d . e1)/a, (d . e2)/b)` and `v = d . e3`, with the existing seam/axis enclosure rules. The liveness scale is `b`.
- Changing `major_dir` changes geometry except for its sign. Axis and reference signs choose parameter orientation and seam, not a different point set. Transformation preserves both squared semi-axes.

For a circle of radius `r` in a plane with unit normal `n`, extruded along unit direction `e`, nonzero `n . e`, the transverse ellipse has squared semi-axes `r^2` and `r^2 (n . e)^2`. A major direction is `n x e` when nonzero. If they are parallel, use `Cylinder`. The original start/end planes trim this carrier obliquely; replacing them by perpendicular caps would change the solid. Build the carrier and trims from the same profile/sweep definition and retain its recipe origins.

No unrestricted `Quadric`, `Extruded` or `Revolved` surface variant is added in this foundation. The analytic reductions cover line extrusion, circular-arc extrusion, meridional line revolution and ring circular-arc revolution. Other valid patches remain named implementation gaps until their domain-aware carrier and full consumer support land. In particular, `shapes.torus-self-intersecting` describes an invalid whole ring torus; it is not evidence that every trimmed circular revolution with that carrier is impossible.

### Acceptance of the migration

Preserve the successful rational corpus and add independent cases for algebraic poses/dimensions, an arc with radius `sqrt(2)`, 45-degree and non-quadratic rational-degree rotations, extrusion along `[1,0,1]`, and oblique circular extrusion. Check exact incidence, normal orientation, component selection and pose/inverse identities. An elliptic-cylinder check must use its defining swept circle as an independent oracle. Compare cold/warm caches and thread counts. Failure to prove a constructor invariant must never produce a successful geometry definition.

## What this piece is for

Piece 2 holds the exact definitions of geometry and answers exact questions about one or two objects at a time.

- Rational points and vectors, and poses: rigid motions with exact rational rotations.
- The analytic surfaces (plane, cylinder, cone, sphere, torus) and the analytic curves (line, circle, ellipse), each with exact rational parameters.
- Geometry stored by definition: an intersection curve ("this branch of where these two surfaces meet") and a point where objects meet ("where these meet, this one").
- The write-once geometry store, `GeomStore`, with typed handles and the fact book of construction facts.
- Decisions about single objects: which side of a surface a point is on, whether a point is on a curve, whether two points are the same point, whether two surfaces are the same point set, and whether two surface normals point the same way.
- Exact parametrisations, evaluated as `Ball` enclosures, for measurement and export.

It does not intersect objects with each other (piece 4) and has no topology (piece 3). It holds no tolerance.

## Number types

Every value in this piece has one of the exact types below. The table is part of the contract.

| What | Type | Meaning |
| --- | --- | --- |
| Defining parameters of surfaces and curves: points, directions, radii, squared radii, slopes | `Rational`, `RPoint3`, `RVector3` | Exact rationals, given by the caller |
| Rotation matrices and translations | `Rational`, `RVector3` | Exact rationals |
| Positions of points, gradients, normals | `Point3`, `Vector3`, whose coordinates are `Real` | Exact real numbers. Building them decides nothing and never fails; their signs are decided through the `Decider` |
| Implicit polynomials and curve equations | `MPoly` in 3 variables with `Rational` coefficients | Exact. Variable 0 is `x`, 1 is `y`, 2 is `z` |
| Parameter values and evaluated points | `Ball` | Proven enclosures at the caller's precision |
| Cached coordinates | `f64` | For display and search only. No decision reads them |

A length that is the square root of a rational is stored squared, as a `Rational`: the radius of a circle cut from a sphere by a plane, and the semi-axes of an ellipse.
Directions are never normalised. A direction is any non-zero rational vector; only its direction carries meaning (and its sign, where the table says so).
Unit vectors such as `ref_dir / |ref_dir|` appear only inside `Real` and `Ball` expressions, where the square root is exact or enclosed.

## Poses

A rotation is given by a non-zero rational quaternion `(w, x, y, z)`, not necessarily of unit length. With `n = w^2 + x^2 + y^2 + z^2`, the rotation matrix is

```
        | w^2 + x^2 - y^2 - z^2   2(xy - wz)              2(xz + wy)            |
(1/n) * | 2(xy + wz)              w^2 - x^2 + y^2 - z^2   2(yz - wx)            |
        | 2(xz - wy)              2(yz + wx)              w^2 - x^2 - y^2 + z^2 |
```

which is exactly orthogonal with determinant 1. It acts on column vectors: `v' = R v`.

- `Rotation::from_quaternion(w, x, y, z)` refuses the zero quaternion with `shapes.zero-quaternion`. Otherwise it holds the matrix above with entries in lowest terms.
- `Rotation` equality and hashing compare the matrix only. So `(w, x, y, z)`, `(-w, -x, -y, -z)` and `(c w, c x, c y, c z)` for any non-zero rational `c` give equal rotations.
- `identity()` is the identity matrix. `matrix()` returns the matrix row by row. `apply(v)` is `R v`; `apply_real(v)` is the same on a `Vector3`.
- `a.then(&b)` applies `a` first, then `b`. Its matrix is `B A`. In quaternions: `from(q1).then(&from(q2)) == from(q2 q1)`, with the Hamilton product `(w1, v1)(w2, v2) = (w1 w2 - v1 . v2, w1 v2 + w2 v1 + v1 x v2)`.
- `inverse()` is the transpose, and equals `from_quaternion(w, -x, -y, -z)`.

A `Pose { rotation, translation }` maps a point `p` to `R p + t`.

- `identity()`: the identity rotation and a zero translation.
- `a.then(&b)` maps `p` to `b(a(p))`: rotation `B A`, translation `B t_a + t_b`.
- `inverse()`: rotation `R^T`, translation `-(R^T t)`.
- `apply_point(p)` is `R p + t`. `apply_vector(v)` is `R v`, without the translation. `apply_real_point(p)` is `R p + t` on a `Point3`.

All pose arithmetic is exact rational arithmetic. It decides nothing and cannot fail, except `from_quaternion` on the zero quaternion.

## Surfaces

### Fields, frame and validation

Every surface carries a frame, computed when needed and never stored:

- `e3 = axis / |axis|` (for a plane, `normal / |normal|`);
- `e1 = ref_dir / |ref_dir|`;
- `e2 = e3 x e1`.

Because `ref_dir` is exactly perpendicular to `axis` (or `normal`), `(e1, e2, e3)` is a right-handed orthonormal frame.

| Surface | Fields | Implicit `f` (positive on the normal side) | Parametrisation `S(u, v)` | Domain `[u, v]` |
| --- | --- | --- | --- | --- |
| Plane | `origin`, `normal`, `ref_dir` | `(x - origin) . normal` | `origin + u e1 + v e2` | `[Line, Line]` |
| Cylinder | `origin`, `axis`, `radius_sq`, `ref_dir` | `\|(x - origin) x axis\|^2 - radius_sq \|axis\|^2` | `origin + sqrt(radius_sq) (cos u e1 + sin u e2) + v e3` | `[Periodic, Line]` |
| Cone | `apex`, `axis`, `slope`, `ref_dir` | `\|(x - apex) x axis\|^2 - slope^2 ((x - apex) . axis)^2` | `apex + v (e3 + slope (cos u e1 + sin u e2))` | `[Periodic, Ray]` |
| Sphere | `center`, `radius`, `axis`, `ref_dir` | `\|x - center\|^2 - radius^2` | `center + radius (cos v (cos u e1 + sin u e2) + sin v e3)` | `[Periodic, Latitude]` |
| Torus | `center`, `axis`, `major`, `minor`, `ref_dir` | `(\|x - c\|^2 + major^2 - minor^2)^2 - 4 major^2 (\|x - c\|^2 - ((x - c) . axis)^2 / \|axis\|^2)` | `c + (major + minor cos v)(cos u e1 + sin u e2) + minor sin v e3` | `[Periodic, Periodic]` |

In the torus row, `c` is `center`.

- A cone's `slope` is the radius gained per unit of distance along the axis, and its `v` is the distance from the apex along `e3`. Its surface is the nappe on the side `axis` points to (`v >= 0`).
- The sign of `axis` matters for a cone (it picks the nappe). For a cylinder, sphere and torus only the line of the axis matters for the point set; the sign of `axis` and the choice of `ref_dir` only move the parametrisation. For a plane, the sign of `normal` picks the orientation, not the point set.
- `ParamKind` describes each parameter: `Line` (all reals), `Periodic` (period `2 pi`, seam at `u = 0`, on the `e1` side), `Latitude` (`[-pi/2, pi/2]`), `Ray` (`[0, infinity)`).
- The seams follow the judge contract (`JUDGE-INTERFACE.md`, seams and poles): for the frames the builders of piece 3 use (`ref_dir = +x` and an axis along the local `z` axis, pointing either way), every `u` seam is on the local half-plane `y = 0, x > 0`, a sphere's poles are on the local `z` axis, and a torus's `v` seam (`v = 0`) is its outer equator.

Each surface struct has public fields and a constructor `new(...) -> Result<Self, Refusal>`. The constructor runs these checks in order and returns the refusal of the first that fails:

1. `shapes.zero-vector`: a direction field (`normal`, `axis` or `ref_dir`) is the zero vector.
2. `shapes.non-positive-radius`: `radius <= 0`, or for a torus `major <= 0` or `minor <= 0`. `shapes.non-positive-slope`: `slope <= 0`. Then `shapes.torus-self-intersecting`: `minor >= major`. So a torus with `major = 0` and `minor = 1` gives `shapes.non-positive-radius`.
3. `shapes.not-perpendicular`: `ref_dir . axis != 0` (for a plane, `ref_dir . normal != 0`).

All checks are exact rational arithmetic. A constructor never panics on any input.

### Orientation and the implicit polynomial

`Surface::implicit()` returns exactly the polynomial `f` of the table, expanded, with no scaling. Tests may therefore compare `implicit().eval(p)` with the table's formula evaluated by hand in exact rationals.
Note the factors this fixes: a cylinder's `f` scales with `|axis|^2`, a cone's with `|axis|^2`, a plane's with `|normal|`.

The normal side of a surface is where `f > 0`; the surface normal at a point is the gradient `grad f`, not normalised.
A plane's normal is its `normal` field. A cylinder's and a cone's normals point away from the axis, a sphere's away from the centre, and a torus's away from its core circle (the circle of radius `major` about the axis).

`Surface::gradient(p)` returns the three partial derivatives of `implicit()` evaluated at `p`, as exact `Real` values. It decides nothing.

The zero set of `f` is the surface's carrier, and every decision of this piece about a surface is about the sign of `f`.
For every kind but the cone, the carrier is the surface, and `grad f` is non-zero at every point of it.
For a cone, the carrier is the double cone: both nappes and the apex. A point on the opposite nappe (where `(p - apex) . axis < 0`) has `f = 0`, so `side` calls it coincident. The apex is the one singular point: `grad f` is zero there.
A piece that needs to know the nappe tests the sign of `(p - apex) . axis` itself.

### Parametrisation, derivatives and parameters

`Surface::evaluate(u, v)` returns three balls that contain `S(u, v)` for every `u` and `v` in the input balls.
Any real input is accepted; the domain is advisory. A cone at `v < 0` evaluates to the opposite nappe, and a sphere at `|v| > pi/2` continues the formula.
The result's precision is that of the inputs (the larger of the two).

`Surface::derivatives(u, v)` returns `[dS/du, dS/dv]`, each three balls with the same containment rule. They are the derivatives of the table's formulas:

| Surface | `dS/du` | `dS/dv` |
| --- | --- | --- |
| Plane | `e1` | `e2` |
| Cylinder | `radius (-sin u e1 + cos u e2)` | `e3` |
| Cone | `v slope (-sin u e1 + cos u e2)` | `e3 + slope (cos u e1 + sin u e2)` |
| Sphere | `radius cos v (-sin u e1 + cos u e2)` | `radius (-sin v (cos u e1 + sin u e2) + cos v e3)` |
| Torus | `(major + minor cos v)(-sin u e1 + cos u e2)` | `minor (-sin v (cos u e1 + sin u e2) + cos v e3)` |

`Surface::parameters_of(p)` takes three balls enclosing a point and returns balls enclosing its parameters, computed by these formulas, with `d = p - origin` (`apex`, `center`):

| Surface | `u` | `v` |
| --- | --- | --- |
| Plane | `d . e1` | `d . e2` |
| Cylinder | `angle(d . e1, d . e2)` | `d . e3` |
| Cone | `angle(d . e1, d . e2)` | `d . e3` |
| Sphere | `angle(d . e1, d . e2)` | `atan2(d . e3, rho)` |
| Torus | `angle(d . e1, d . e2)` | `angle(rho - major, d . e3)` |

Here `rho = sqrt((d . e1)^2 + (d . e2)^2)`, and `angle(a, b)` is the unique `theta` in `[0, 2 pi)` with `(a, b) = |(a, b)| (cos theta, sin theta)`.
For a point on the surface, `S(parameters_of(p)) = p`. For a point off it, the same formulas apply.

- Correctness: when the result is `Some`, its balls contain the formula's values for every point in the input box.
- `None` is allowed when the box meets a place where a formula is undefined or jumps: the closed seam half-plane `{d . e2 = 0, d . e1 >= 0}` and the axis `{d . e1 = d . e2 = 0}` for every `u`; for a torus `v`, also the closed half-plane `{d . e3 = 0, rho >= major}` and the core circle; for a sphere `v`, the centre. A plane never gives `None`.
- Liveness, for tests: at precision 64 or more, a box of radius zero whose point is at a distance of at least `L / 100` from every set listed above must give `Some`. `L` is the `radius`; for a torus, `minor`; for a cone, `slope` times the point's distance from the apex along the axis.

`Surface::domain()` returns the domain column of the table.

### Moving a surface

`Surface::transformed(pose)` moves every point field with `apply_point` and every direction field with `apply_vector`, and keeps radii and slopes. The result is exact and field-for-field comparable with `==`.
Because rotations preserve dot products, cross products and lengths, the implicit commutes with posing: `transformed(P).implicit().eval(P(p)) == implicit().eval(p)` exactly, for every rational point `p`.

`Surface` equality (`==`) is structural: two planes with different origins on the same plane are not `==`. Whether two surfaces are the same point set is the decision `same_surface`.

## Curves

| Curve | Fields | Natural direction | Parametrisation `C(t)` | Closed |
| --- | --- | --- | --- | --- |
| Line | `origin`, `dir` | `dir` | `origin + t dir` (`t` is not arc length) | no |
| Circle | `center`, `normal`, `radius_sq`, `ref_dir` | counter-clockwise seen from the tip of `normal` | `center + sqrt(radius_sq) (cos t e1 + sin t e2)` | yes |
| Ellipse | `center`, `normal`, `major_dir`, `major_sq`, `minor_sq` | counter-clockwise seen from the tip of `normal` | `center + sqrt(major_sq) cos t e1 + sqrt(minor_sq) sin t e2` | yes |
| Intersection | `surfaces: [SurfaceId; 2]`, `witness: PointId` | `grad f0 x grad f1` at regular points | none here; piece 4 traces it | piece 4 decides |

For circles and ellipses, `e1` is `ref_dir` (or `major_dir`) normalised, `e3` is `normal` normalised, and `e2 = e3 x e1`. The parameter `t` is periodic with period `2 pi`, and `t = 0` is the point in the `e1` direction.

Each analytic curve struct has public fields and a constructor `new(...) -> Result<Self, Refusal>` that runs these checks in order:

1. `shapes.zero-vector`: `dir`, `normal`, `ref_dir` or `major_dir` is the zero vector.
2. `shapes.non-positive-radius`: a circle's `radius_sq <= 0`. `shapes.ellipse-axes`: an ellipse without `0 < minor_sq < major_sq`.
3. `shapes.not-perpendicular`: `ref_dir . normal != 0` or `major_dir . normal != 0`.

An ellipse with equal semi-axes is a circle and must be built as a `Circle`, so `minor_sq == major_sq` is refused with `shapes.ellipse-axes`.

`Curve::equations(store)` returns exactly these polynomials, in this order:

| Curve | Equations |
| --- | --- |
| Line | the components of `(x - origin) x dir`, in `x`, `y`, `z` order, leaving out any that is the zero polynomial (two or three equations) |
| Circle | `(x - center) . normal`, then `\|x - center\|^2 - radius_sq` |
| Ellipse | `(x - center) . normal`, then `((x - center) . d1)^2 / (\|d1\|^2 major_sq) + ((x - center) . d2)^2 / (\|d2\|^2 minor_sq) - 1`, with `d1 = major_dir` and `d2 = normal x major_dir` |
| Intersection | `store.implicit(surfaces[0])`, then `store.implicit(surfaces[1])` |

For a line, a circle and an ellipse, the common real zero set of the equations is exactly the curve. For an intersection curve it is the whole meeting set of the two surfaces, which contains the curve.

- `Curve::is_closed()`: `Some(false)` for a line, `Some(true)` for a circle and an ellipse, `None` for an intersection curve.
- `Curve::transformed(pose)`: moves point fields with `apply_point` and directions with `apply_vector`, keeps `radius_sq`, `major_sq` and `minor_sq`. `None` for an intersection curve, whose surfaces are handles: to move one, move its surfaces and witness and add a new curve.
- `Curve::evaluate(t)`: three balls containing `C(t)` for every `t` in the ball. `None` for an intersection curve.
- `Curve::parameter_origin(store)`: the exact point at `t = 0` as a `Point3`. A line's `origin`; a circle's `center + sqrt(radius_sq / |ref_dir|^2) ref_dir`; an ellipse's `center + sqrt(major_sq / |major_dir|^2) major_dir`; an intersection curve's witness position. Piece 4 starts closed-curve orderings there, and piece 7 puts a ring edge's vertex there.

### Intersection curves

An intersection curve is the connected component, of the meeting set of its two surfaces with its singular points removed, that contains `witness`. The witness must be a regular point of that set: a point of both surfaces where their gradients are not parallel.
Piece 4 creates these curves and guarantees the witness. This piece checks neither the witness nor which branch a point lies on; it treats the curve through its two surfaces.
The two surfaces must be different handles.

## Points

A point is stored by definition:

- `PointDef::Explicit(RPoint3)`: an exact rational point.
- `PointDef::Meet(Meet)`: a point where objects meet. `Meet::of` names them: `MeetOf::Surfaces` (two or more distinct surfaces), `MeetOf::CurveSurface(curve, surface)`, or `MeetOf::Curves(a, b)`. `Meet::position` says which point of the meeting it is: `Position::Closed(Point3)`, exact closed-form coordinates in `Real`; or `Position::Isolated(IsolatedRoot)`, a certified root of a polynomial system in the three variables `x, y, z`.

The creator of a `Meet` point (pieces 3 and 4) guarantees that its position lies on every object it names. This piece does not check it.

- `GeomStore::position(p)` returns the exact position as a `Point3`: the rational coordinates of an explicit point, a clone of a closed position, or `coordinate(0)`, `coordinate(1)`, `coordinate(2)` of an isolated root.
- `GeomStore::approx_point(p)` returns cached `f64` coordinates for display and search. For an explicit point each value is `to_f64_nearest` of the coordinate. For a meeting point each value lies inside `Real::approx()` of the exact coordinate. No decision reads these values.
- `GeomStore::enclose_point(p, max_width, dec)` returns a proven ball for each coordinate. It makes one `Decider::enclose` call per coordinate, `x`, then `y`, then `z`, with that `max_width`, so its statistics are those of those calls. When a call refuses, it returns that refusal at once and makes no further call. When the budget runs out first, the balls are the tightest reached; the caller checks the width.

## The store

`Handle<T>` is a typed index; `Arena<T>` is an append-only list.

- `Arena::push` returns handles in order from index 0. `get` returns the entity; `len`, `is_empty` and `handles` (all handles in order) describe the arena. Nothing in an arena ever changes once added. There is no `get_mut`.
- `GeomStore` holds three arenas (surfaces, curves, points) and one `FactBook<GeomKey>`. `GeomKey` names a surface, a curve or a point. Its order is by variant (surface, then curve, then point), then by index.
- `add_surface`, `add_curve` and `add_point` never deduplicate. Adding the same definition twice gives two handles; whether they name the same object is a decision.
- `GeomStore::implicit(s)` equals `surface(s).implicit()`. It may be cached.
- `GeomStore::normal(s, p)` is `surface(s).gradient(&position(p))`.

Two constructions record facts, because the definition makes them true:

- `add_curve` of an intersection curve records `Incident(curve, surfaces[0])`, then `Incident(curve, surfaces[1])`.
- `add_point` of a `Meet` records `Incident(point, x)` for every object `x` it names, in the order named (`CurveSurface(c, s)` records the curve, then the surface).

Their origin is `FactOrigin::Construction("shapes.definition")` and their `at` is `None`. Adding an analytic curve, a surface or an explicit point records nothing.
`record_fact` passes a fact to the fact book and returns what `FactBook::record` returns. Recording a fact that is not true is a programming error of the caller.

## Decisions

Six decisions answer questions about one or two objects. Each takes a `&mut Decider` and returns `Outcome`.

| Decision | Question | `Proven` | `ProvenCoincident` |
| --- | --- | --- | --- |
| `side(s, p)`, `side_of(s, &point)` | sign of `f_s` at the point | `Side::Front` (`f > 0`), `Side::Back` (`f < 0`) | the point is on the carrier (`f = 0`) |
| `point_on_curve(c, p)` | is the point on the common zero set of `equations()` | `Apart` | on it |
| `same_point(a, b)` | are the two positions equal | `Apart` | equal |
| `same_surface(a, b)` | are the two surfaces the same point set | `Apart` | the same point set |
| `normals_agree(a, b, at)` | sign of `grad f_a . grad f_b` at the point | `StrictSign` of the dot product | the dot product is zero (perpendicular normals, or a zero gradient at a singular point) |

### Facts first

Each decision looks in the fact book before it does any arithmetic.

| Decision | Settled by a fact when |
| --- | --- |
| `side(s, p)` | `holds(Incident, Point(p), Surface(s))` |
| `point_on_curve(c, p)` | `holds(Incident, Point(p), Curve(c))`; or, for an intersection curve, `holds(Incident, Point(p), Surface(s))` for both of its surfaces |
| `same_point(a, b)` | `holds(Coincident, Point(a), Point(b))`, which is true when `a == b` |
| `same_surface(a, b)` | `holds(Coincident, Surface(a), Surface(b))`, which is true when `a == b` |
| `side_of`, `normals_agree` | never: `side_of` has no point key, and no fact kind says which way two normals point |

The fact book closes `Incident` chains and crosses `Coincident` links (piece 1, facts), so a point incident on a curve that is incident on a surface is on that surface by fact.

A decision settled by a fact returns `ProvenCoincident` with `ZeroProof::ConstructionFact`, calls `Decider::note_fact_decision` exactly once, and calls no other `Decider` method. So it adds one to `decisions` and one to `coincident_fact`, and adds no work.
It answers even when the `Decider`'s deadline has passed, because it does no work.

### Then numbers

`same_surface` compares the defining data of two analytic surfaces with exact rational arithmetic, without the `Decider`. The sign of a rational is known exactly (piece 1, sign), so this needs no decision. Off the fact path it returns only `Proven(Apart)` or `ProvenCoincident` with `ZeroProof::ExactArithmetic`, and it adds nothing to the statistics. It never returns `NotExact` or `Refused`.

- Different kinds are never the same point set: `Proven(Apart)`.
- Two planes: same point set exactly when `normal_a x normal_b = 0` and `(origin_b - origin_a) . normal_a = 0`.
- Two cylinders: `radius_sq_a = radius_sq_b`, `axis_a x axis_b = 0` and `(origin_b - origin_a) x axis_a = 0`.
- Two cones: `apex_a = apex_b`, `slope_a = slope_b`, `axis_a x axis_b = 0` and `axis_a . axis_b > 0`.
- Two spheres: `center_a = center_b` and `radius_a = radius_b`.
- Two tori: `center_a = center_b`, `major_a = major_b`, `minor_a = minor_b` and `axis_a x axis_b = 0`.

`ref_dir`, the sign of a plane's normal and the sign of a cylinder's or torus's axis never matter here.

The other decisions build exact `Real` values from the definitions and decide them with `Decider::sign` (or `compare`):

- `side` and `side_of` decide the sign of `f_s(position)`.
- `point_on_curve` decides the value of each equation at the position, one `Decider::sign` per equation in the order of `equations()`. `same_point` decides the coordinate differences `a - b`, one `Decider::sign` each, in the order `x`, `y`, `z`.
  The first of these that ends `Proven` non-zero gives `Proven(Apart)`, and the first that ends `Refused` gives that refusal, whichever comes first; no further decision is made. When neither happens, the outcome is `NotExact` if any ended `NotExact`, and `ProvenCoincident` if all were proven zero.
- `normals_agree` decides the sign of the dot product of the two gradients at `position(at)`. The point need not lie on either surface.

So `side`, `side_of` and `normals_agree` make exactly one piece 1 decision off the fact path, and `point_on_curve` and `same_point` make one per equation or coordinate until they stop. Which proof path each decision takes, and its work, are piece 1's.

### When a decision cannot prove

| Piece 1 returned | This piece returns |
| --- | --- |
| `Proven` or `ProvenCoincident` | the matching outcome, with the `Coincidence` of piece 1 |
| `NotExact` | `NotExact`, keeping piece 1's `code` and `detail`, with the bound converted as below |
| `Refused` | the same `Refused`, unchanged (for example `numbers.deadline`, `numbers.division-by-zero` for an undefined position) |

A `NotExact` from `side`, `side_of`, `point_on_curve` or `same_point` carries a **distance bound**: a rational `b` with a proof about distances, in model units, as listed below. Callers that mark it use `BoundKind::Distance`.

- `side`, `side_of`: the distance from the point to the carrier of `s`.
- `point_on_curve` on a line, circle or ellipse: the distance from the point to the curve.
- `point_on_curve` on an intersection curve: the distance from the point to each of the two surfaces' carriers. This is not a bound on the distance to the curve: where the surfaces meet at a small angle, a point close to both can be far from their meeting. No sharper bound follows from the two implicit values alone.
- `same_point`: the distance between the two points.

A `NotExact` from `normals_agree` keeps piece 1's bound on the dot product itself. Callers mark it with `BoundKind::Value`.

Converting a bound on an implicit value into a distance is the creator's job. Each conversion must be proven, and any rational upper bound of the right side of these inequalities is valid. With `b` the bound on `|f|`:

| Object | Distance from the point to the object is at most |
| --- | --- |
| Plane | `b / \|normal\|` (an equality) |
| Cylinder | `b / (\|axis\|^2 radius)`, since `f = \|axis\|^2 (rho^2 - radius^2)` |
| Cone | `sqrt(b / \|axis\|^2)`, since `f = \|axis\|^2 (rho^2 - slope^2 h^2)` and `\|rho - slope \|h\|\| <= sqrt(\|rho^2 - slope^2 h^2\|)` |
| Sphere | `b / radius`, since `f = (d - radius)(d + radius)` |
| Torus | `b / (minor (major^2 - minor^2))`, since `f = (D^2 - minor^2)(D'^2 - minor^2)` with `D'^2 >= major^2` |
| Line, bounds `b_i` on its equations | `sqrt(sum b_i^2) / \|dir\|` |
| Circle, bounds `b1`, `b2` | `sqrt(b2 + 2 h^2)` with `h = b1 / \|normal\|` |
| Ellipse, bounds `b1`, `b2` | `sqrt(major_sq) b2 + b1 / \|normal\|` |
| Two points, bounds `bx`, `by`, `bz` on the coordinate differences | `sqrt(bx^2 + by^2 + bz^2)` |

Here `rho` is the distance from the axis, `h` the signed distance along it from the apex, `d` the distance from the centre, `D` the distance from the torus's core circle and `D'` the distance from its point opposite. An equation proven zero has bound 0.
The `code` and `detail` of the converted `NotExact` are those of the first piece 1 decision, in the order above, that ended `NotExact`.

### Rules shared by every decision

- The decisions take `&self` and never record facts. A caller that wants to keep a proven zero records it with `record_fact` and `FactOrigin::Proof(proof)`.
- When the deadline has passed, every decision that reaches the `Decider` returns `Refused("numbers.deadline")`. Fact-settled decisions and `same_surface` still answer.
- A `Proven` or `ProvenCoincident` outcome is never wrong (DESIGN.md, R2). A point off the carrier never comes back `ProvenCoincident`; a point on it never comes back `Proven`.
- Repeating a decision on the same store with a `Decider` in the same state gives the same outcome and adds the same statistics, on every machine and thread.

## How this piece uses piece 1

- **Outcomes.** Every decision returns `Outcome`, mapped from piece 1's as in the table above. Constructors return `Result<_, Refusal>`.
- **Budgets.** Every numeric decision runs on the caller's `Decider` and its `BudgetPolicy`. This piece never creates a `Decider` or a policy, and holds no budget of its own.
- **Facts.** The store holds one `FactBook<GeomKey>`. This piece records only the incidences that a definition makes true, consults facts before numbers, and reports fact-settled decisions through `note_fact_decision`.
- **Marks.** This piece never builds a `Mark`. Its callers build one from each `NotExact` they accept, with `BoundKind::Distance` for `side`, `side_of`, `point_on_curve` and `same_point`, and `BoundKind::Value` for `normals_agree`.
- **Refusals.** This piece adds only the constructor codes below. Every `numbers.*` refusal passes through unchanged.

## Refusal codes

| Code | When |
| --- | --- |
| `shapes.zero-quaternion` | `Rotation::from_quaternion(0, 0, 0, 0)` |
| `shapes.zero-vector` | a direction field is zero |
| `shapes.not-perpendicular` | `ref_dir` or `major_dir` is not exactly perpendicular to the axis or normal |
| `shapes.non-positive-radius` | a radius, `major`, `minor` or `radius_sq` is zero or negative |
| `shapes.non-positive-slope` | a cone's `slope` is zero or negative |
| `shapes.torus-self-intersecting` | `minor >= major` |
| `shapes.ellipse-axes` | an ellipse without `0 < minor_sq < major_sq` |

This piece adds no `NotExact` code; it passes piece 1's on.

## Programming errors

These may panic. None can come from a question file, because pieces 3 and 4 build every definition from validated data:

- a handle from another store, or out of range;
- a surface or curve struct built directly with fields that its `new` would refuse;
- `MeetOf::Surfaces` with fewer than two surfaces or a repeated surface, or an intersection curve on one surface twice;
- `Position::Isolated` whose system is not in exactly three variables;
- a `Real` position or a recorded fact that is not true.

## Exact public interface

The crate re-exports these types at its root. Everything listed is `pub`; fields not listed are private. The compiled crate in `crates/catallax-shapes` is authoritative if a signature listing becomes stale; update this contract with any API change.

```rust
// Handles and arenas
pub struct Handle<T>;                      // Copy, Eq, Ord, Hash, Debug ("#index")
impl<T> Handle<T> { pub fn index(&self) -> u32; }
#[derive(Clone, Debug)]
pub struct Arena<T>;                       // Default
impl<T> Arena<T> {
    pub fn new() -> Self;
    pub fn push(&mut self, value: T) -> Handle<T>;
    pub fn get(&self, handle: Handle<T>) -> &T;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    pub fn handles(&self) -> Vec<Handle<T>>;
}
pub type SurfaceId = Handle<Surface>;
pub type CurveId = Handle<Curve>;
pub type PointId = Handle<PointDef>;
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum GeomKey { Surface(SurfaceId), Curve(CurveId), Point(PointId) }

// Points and vectors
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct RPoint3 { pub x: Rational, pub y: Rational, pub z: Rational }
impl RPoint3 {
    pub fn new(x: Rational, y: Rational, z: Rational) -> Self;
    pub fn origin() -> Self;
    pub fn minus(&self, other: &RPoint3) -> RVector3;
    pub fn plus(&self, v: &RVector3) -> RPoint3;
    pub fn to_real(&self) -> Point3;
}
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct RVector3 { pub x: Rational, pub y: Rational, pub z: Rational }
impl RVector3 {
    pub fn new(x: Rational, y: Rational, z: Rational) -> Self;
    pub fn zero() -> Self;
    pub fn is_zero(&self) -> bool;
    pub fn plus(&self, other: &RVector3) -> RVector3;
    pub fn minus(&self, other: &RVector3) -> RVector3;
    pub fn scaled(&self, c: &Rational) -> RVector3;
    pub fn dot(&self, other: &RVector3) -> Rational;
    pub fn cross(&self, other: &RVector3) -> RVector3;
    pub fn norm_squared(&self) -> Rational;
    pub fn to_real(&self) -> Vector3;
}
#[derive(Clone, Debug)]
pub struct Point3 { pub x: Real, pub y: Real, pub z: Real }
impl Point3 {
    pub fn new(x: Real, y: Real, z: Real) -> Self;
    pub fn minus(&self, other: &Point3) -> Vector3;
    pub fn plus(&self, v: &Vector3) -> Point3;
}
#[derive(Clone, Debug)]
pub struct Vector3 { pub x: Real, pub y: Real, pub z: Real }
impl Vector3 {
    pub fn new(x: Real, y: Real, z: Real) -> Self;
    pub fn plus(&self, other: &Vector3) -> Vector3;
    pub fn minus(&self, other: &Vector3) -> Vector3;
    pub fn scaled(&self, c: &Real) -> Vector3;
    pub fn dot(&self, other: &Vector3) -> Real;
    pub fn cross(&self, other: &Vector3) -> Vector3;
    pub fn norm_squared(&self) -> Real;
}
// From<RPoint3> for Point3, From<RVector3> for Vector3

// Poses
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Rotation;
impl Rotation {
    pub fn from_quaternion(w: Rational, x: Rational, y: Rational, z: Rational) -> Result<Self, Refusal>;
    pub fn identity() -> Self;
    pub fn matrix(&self) -> &[[Rational; 3]; 3];
    pub fn then(&self, next: &Rotation) -> Rotation;
    pub fn inverse(&self) -> Rotation;
    pub fn apply(&self, v: &RVector3) -> RVector3;
    pub fn apply_real(&self, v: &Vector3) -> Vector3;
}
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Pose { pub rotation: Rotation, pub translation: RVector3 }
impl Pose {
    pub fn identity() -> Self;
    pub fn new(rotation: Rotation, translation: RVector3) -> Self;
    pub fn then(&self, next: &Pose) -> Pose;
    pub fn inverse(&self) -> Pose;
    pub fn apply_point(&self, p: &RPoint3) -> RPoint3;
    pub fn apply_vector(&self, v: &RVector3) -> RVector3;
    pub fn apply_real_point(&self, p: &Point3) -> Point3;
}

// Surfaces
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum SurfaceKind { Plane, Cylinder, Cone, Sphere, Torus }
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum ParamKind { Line, Periodic, Latitude, Ray }
// Each struct: #[derive(Clone, PartialEq, Eq, Hash, Debug)], public fields, validating new(...)
pub struct Plane { pub origin: RPoint3, pub normal: RVector3, pub ref_dir: RVector3 }
pub struct Cylinder { pub origin: Point3, pub axis: Vector3, pub radius_sq: Real, pub ref_dir: Vector3 }
pub struct Cone { pub apex: RPoint3, pub axis: RVector3, pub slope: Rational, pub ref_dir: RVector3 }
pub struct Sphere { pub center: RPoint3, pub radius: Rational, pub axis: RVector3, pub ref_dir: RVector3 }
pub struct Torus { pub center: RPoint3, pub axis: RVector3, pub major: Rational, pub minor: Rational, pub ref_dir: RVector3 }
impl Plane { pub fn new(origin: RPoint3, normal: RVector3, ref_dir: RVector3) -> Result<Self, Refusal>; }
impl Cylinder {
    pub fn new(origin: Point3, axis: Vector3, radius: Real, ref_dir: Vector3, dec: &mut Decider) -> Result<Self, Refusal>;
    pub fn from_radius_sq(origin: Point3, axis: Vector3, radius_sq: Real, ref_dir: Vector3, dec: &mut Decider) -> Result<Self, Refusal>;
    pub fn radius(&self) -> Real;
    pub fn radial_reference(&self) -> Vector3;
}
impl Cone { pub fn new(apex: RPoint3, axis: RVector3, slope: Rational, ref_dir: RVector3) -> Result<Self, Refusal>; }
impl Sphere { pub fn new(center: RPoint3, radius: Rational, axis: RVector3, ref_dir: RVector3) -> Result<Self, Refusal>; }
impl Torus { pub fn new(center: RPoint3, axis: RVector3, major: Rational, minor: Rational, ref_dir: RVector3) -> Result<Self, Refusal>; }
#[non_exhaustive]
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Surface { Plane(Plane), Cylinder(Cylinder), Cone(Cone), Sphere(Sphere), Torus(Torus) }
impl Surface {
    pub fn kind(&self) -> SurfaceKind;
    pub fn implicit(&self) -> MPoly;
    pub fn gradient(&self, p: &Point3) -> Vector3;
    pub fn transformed(&self, pose: &Pose) -> Surface;
    pub fn domain(&self) -> [ParamKind; 2];
    pub fn evaluate(&self, u: &Ball, v: &Ball) -> [Ball; 3];
    pub fn derivatives(&self, u: &Ball, v: &Ball) -> [[Ball; 3]; 2];
    pub fn parameters_of(&self, p: &[Ball; 3]) -> Option<[Ball; 2]>;
}

// Curves
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum CurveKind { Line, Circle, Ellipse, Intersection }
// Each struct: #[derive(Clone, PartialEq, Eq, Hash, Debug)], public fields
pub struct Line { pub origin: RPoint3, pub dir: RVector3 }
pub struct Circle { pub center: RPoint3, pub normal: RVector3, pub radius_sq: Rational, pub ref_dir: RVector3 }
pub struct Ellipse { pub center: RPoint3, pub normal: RVector3, pub major_dir: RVector3, pub major_sq: Rational, pub minor_sq: Rational }
pub struct IntersectionCurve { pub surfaces: [SurfaceId; 2], pub witness: PointId }
impl Line { pub fn new(origin: RPoint3, dir: RVector3) -> Result<Self, Refusal>; }
impl Circle { pub fn new(center: RPoint3, normal: RVector3, radius_sq: Rational, ref_dir: RVector3) -> Result<Self, Refusal>; }
impl Ellipse { pub fn new(center: RPoint3, normal: RVector3, major_dir: RVector3, major_sq: Rational, minor_sq: Rational) -> Result<Self, Refusal>; }
#[non_exhaustive]
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Curve { Line(Line), Circle(Circle), Ellipse(Ellipse), Intersection(IntersectionCurve) }
impl Curve {
    pub fn kind(&self) -> CurveKind;
    pub fn equations(&self, store: &GeomStore) -> Vec<MPoly>;
    pub fn is_closed(&self) -> Option<bool>;
    pub fn transformed(&self, pose: &Pose) -> Option<Curve>;
    pub fn evaluate(&self, t: &Ball) -> Option<[Ball; 3]>;
    pub fn parameter_origin(&self, store: &GeomStore) -> Point3;
}

// Points
#[derive(Clone, Debug)]
pub enum PointDef { Explicit(RPoint3), Meet(Meet) }
#[derive(Clone, Debug)]
pub struct Meet { pub of: MeetOf, pub position: Position }
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum MeetOf { Surfaces(Vec<SurfaceId>), CurveSurface(CurveId, SurfaceId), Curves(CurveId, CurveId) }
#[derive(Clone, Debug)]
pub enum Position { Closed(Point3), Isolated(IsolatedRoot) }

// Decision answers
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Side { Front, Back }
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Apart;

// The store
#[derive(Clone, Debug, Default)]
pub struct GeomStore;
impl GeomStore {
    pub fn new() -> Self;
    pub fn add_surface(&mut self, s: Surface) -> SurfaceId;
    pub fn add_curve(&mut self, c: Curve) -> CurveId;
    pub fn add_point(&mut self, p: PointDef) -> PointId;
    pub fn surface(&self, id: SurfaceId) -> &Surface;
    pub fn curve(&self, id: CurveId) -> &Curve;
    pub fn point(&self, id: PointId) -> &PointDef;
    pub fn surface_count(&self) -> usize;
    pub fn curve_count(&self) -> usize;
    pub fn point_count(&self) -> usize;
    pub fn facts(&self) -> &FactBook<GeomKey>;
    pub fn record_fact(&mut self, fact: Fact<GeomKey>) -> bool;
    pub fn implicit(&self, s: SurfaceId) -> MPoly;
    pub fn position(&self, p: PointId) -> Point3;
    pub fn approx_point(&self, p: PointId) -> [f64; 3];
    pub fn enclose_point(&self, p: PointId, max_width: &Rational, dec: &mut Decider) -> Result<[Ball; 3], Refusal>;
    pub fn normal(&self, s: SurfaceId, p: PointId) -> Vector3;
    pub fn side(&self, s: SurfaceId, p: PointId, dec: &mut Decider) -> Outcome<Side>;
    pub fn side_of(&self, s: SurfaceId, p: &Point3, dec: &mut Decider) -> Outcome<Side>;
    pub fn point_on_curve(&self, c: CurveId, p: PointId, dec: &mut Decider) -> Outcome<Apart>;
    pub fn same_point(&self, a: PointId, b: PointId, dec: &mut Decider) -> Outcome<Apart>;
    pub fn same_surface(&self, a: SurfaceId, b: SurfaceId, dec: &mut Decider) -> Outcome<Apart>;
    pub fn normals_agree(&self, a: SurfaceId, b: SurfaceId, at: PointId, dec: &mut Decider) -> Outcome<StrictSign>;
}
```

## What defines correct behaviour

The test writer builds tests from the cases below. Every case has an answer that maths gives, so tests compare against that answer, never against a run of the code.
Independent oracles are allowed: a third-party big-rational crate used only as a dev-dependency to evaluate the table's formulas exactly, or values generated by a separate program, as long as nothing of this piece's own code computes the expected answer.

Exact points on curved surfaces come from rational parametrisations. With `c = (1 - k^2) / (1 + k^2)` and `s = 2k / (1 + k^2)` for a rational `k`, `c^2 + s^2 = 1` exactly, so these points are rational and lie exactly on the surfaces with axis `z` and `ref_dir` `x`: `(r c, r s, v)` on a cylinder; `(slope h c, slope h s, h)` on a cone with apex at the origin; `(r c1 c2, r c1 s2, r s1)` on a sphere centred at the origin; `((major + minor c2) c1, (major + minor c2) s1, minor s2)` on a torus. Posed by a rational pose they stay rational and stay on the posed surface.

Where a case says "the outcome kinds are unchanged", a `Proven` answer keeps its value.

### Poses

1. `(1, 0, 0, 0)`, `(2, 0, 0, 0)`, `(-3, 0, 0, 0)` and `(1/7, 0, 0, 0)` all give `Rotation::identity()`.
2. `(1, 0, 0, 1)` gives `[[0, -1, 0], [1, 0, 0], [0, 0, 1]]`, a quarter turn about `z` taking `x` to `y`.
3. `(1, 1, 1, 1)` gives `[[0, 0, 1], [1, 0, 0], [0, 1, 0]]`, taking `x` to `y`, `y` to `z` and `z` to `x`.
4. `(2, 1, 0, 0)` gives `[[1, 0, 0], [0, 3/5, -4/5], [0, 4/5, 3/5]]`.
5. `(1, 2, 3, 4)` gives `[[-2/3, 2/15, 11/15], [2/3, -1/3, 2/3], [1/3, 14/15, 2/15]]`. `(0, 0, 0, 1)` gives `diag(-1, -1, 1)`.
6. `q`, `-q`, `5 q` and `q / 7` give rotations that are `==` and hash the same.
7. Composition: with `q1 = (1, 2, 3, 4)` and `q2 = (2, -1, 0, 5)`, `from(q1).then(&from(q2)) == from(-16, -12, 20, 10)`, the Hamilton product `q2 q1`. For random rational quaternions, `then` equals the rotation of the Hamilton product.
8. For random rational quaternions with components of up to 64 bits: `R R^T = I` and `det R = 1` exactly; `inverse()` equals the transpose and `from_quaternion(w, -x, -y, -z)`; `r.then(&r.inverse())` is the identity.
9. The zero quaternion is refused with `shapes.zero-quaternion` and never panics.
10. The pose with rotation `(1, 0, 0, 1)` and translation `(1, 2, 3)` maps `(1, 0, 0)` to `(1, 3, 3)`; its inverse maps `(1, 3, 3)` back. For random poses `a`, `b` and points `p`: `a.then(&b).apply_point(p) == b.apply_point(&a.apply_point(p))`, `a.then(&a.inverse())` is the identity, and `apply_vector` ignores the translation.

### Construction and implicit polynomials

11. Each refusal code with a minimal example, and none panics: a plane with normal `(0, 0, 0)` gives `shapes.zero-vector`; a plane with normal `(0, 0, 1)` and `ref_dir` `(1, 0, 1)` gives `shapes.not-perpendicular`; a cylinder of radius 0 or -1 gives `shapes.non-positive-radius`; a cone of slope 0 gives `shapes.non-positive-slope`; a torus with `minor = major` or `minor > major` gives `shapes.torus-self-intersecting` and with `minor = 0` gives `shapes.non-positive-radius`; a circle with `radius_sq = 0` gives `shapes.non-positive-radius`; ellipses with `minor_sq = major_sq`, `minor_sq > major_sq` or `minor_sq = 0` give `shapes.ellipse-axes`. Order: a cylinder with a zero axis and radius -1 gives `shapes.zero-vector`; one with radius -1 and a `ref_dir` not perpendicular gives `shapes.non-positive-radius`.
12. Exact implicit polynomials:
    - plane, origin `(0, 0, 1)`, normal `(0, 0, 2)`, `ref_dir` `(1, 0, 0)`: `f = 2z - 2`;
    - cylinder, origin `0`, axis `(0, 0, 1)`, radius 3: `f = x^2 + y^2 - 9`; with axis `(0, 0, 2)`: `f = 4x^2 + 4y^2 - 36`, so `f(1, 1, 1) = -28`;
    - cone, apex `0`, axis `(0, 0, 1)`, slope `1/2`: `f = x^2 + y^2 - z^2/4`;
    - sphere, centre `(1, 2, 3)`, radius 2: `f = x^2 + y^2 + z^2 - 2x - 4y - 6z + 10`;
    - torus, centre `0`, axis `(0, 0, 1)`, major 3, minor 1: `f = (x^2 + y^2 + z^2 + 8)^2 - 36(x^2 + y^2)`, with `f(4, 0, 0) = f(2, 0, 0) = f(3, 0, 1) = 0`, `f(0, 0, 0) = 64` and `f(3, 0, 0) = -35`.
13. For random rational poses `P`, random surfaces of every kind and random rational points `p`: `S.transformed(P).implicit().eval(P(p)) == S.implicit().eval(p)` exactly.
14. The cylinder with origin `0`, axis `(0, 0, 1)`, radius 5, `ref_dir` `(1, 0, 0)`, moved by the pose with rotation `(2, 1, 0, 0)` and translation `(10^6, -7, 1/3)`, is `==` to the cylinder with origin `(10^6, -7, 1/3)`, axis `(0, -4/5, 3/5)`, radius 5, `ref_dir` `(1, 0, 0)`. The local point `(3, 4, 7)` moves to `(1000003, -51/5, 116/15)`, where the moved cylinder's `f` is exactly 0.
15. Gradients: the sphere with centre `0` and radius 3 at `(1, 2, 2)` has gradient `(2, 4, 4)`; a plane's gradient is its `normal` everywhere; the cone of case 12 has gradient `(0, 0, 0)` at its apex.

### Side

16. Points on the carrier give `ProvenCoincident`: the plane of case 12 at `(5, 7, 1)`; the cylinder of radius 5 about the `z` axis at `(3, 4, 7)`; the cone of case 12 at `(3, 4, 10)` and at its apex; the sphere of radius 3 at the origin at `(1, 2, 2)`; the torus of case 12 at `(54/25, 72/25, 4/5)`, `(4, 0, 0)`, `(2, 0, 0)` and `(3, 0, 1)`; the moved cylinder and point of case 14.
17. The opposite nappe: the cone of case 12 at `(3, 4, -10)` gives `ProvenCoincident`, because `side` answers about the carrier.
18. Near points, with `e = 10^-30`, then again with `e = 10^-7` and `e = 10^-3`: the cylinder at `(3(1 + e), 4(1 + e), 7)` is `Front` and at `(3(1 - e), 4(1 - e), 7)` is `Back`; the sphere at `(1, 2, 2)(1 + e)` is `Front` and at `(1 - e)` times is `Back`; the plane at `(5, 7, 1 + e)` is `Front` and at `(5, 7, 1 - e)` is `Back`; the cone at `(3, 4, 10(1 - e))` is `Front` and at `(3, 4, 10(1 + e))` is `Back`; the torus at `(54/25, 72/25, 4/5 (1 + e))` is `Front` and at `(54/25, 72/25, 4/5 (1 - e))` is `Back`. The same with every surface and point translated by `(10^6, 10^6, 10^6)`. Each must be `Proven` with that side under the default budget. A generated family repeats this on random exact points of every kind, posed by random poses, with the expected side the sign of the table's formula evaluated by an exact rational oracle.
19. Points with square-root coordinates, as `Position::Closed`: on the cylinder of radius 1 about the `z` axis, `(sqrt(1/2), sqrt(1/2), 0)` is `ProvenCoincident`; `(sqrt(1/2) + 10^-30, sqrt(1/2), 0)` is `Front` (its distance from the cylinder is about `7.071e-31`); `(sqrt(1/2) - 10^-30, sqrt(1/2), 0)` is `Back`. On the sphere of radius 1 at the origin, `(1/2, 1/2, sqrt(1/2))` is `ProvenCoincident`.
20. Points from certified roots, as `Position::Isolated`: the system of three spheres of radius 5 centred at `(0, 0, 0)`, `(6, 0, 0)` and `(0, 6, 0)` has the roots `(3, 3, sqrt(7))` and `(3, 3, -sqrt(7))` (piece 1, case 51). Certified in the boxes `[2, 4] x [2, 4] x [2, 3]` and `[2, 4] x [2, 4] x [-3, -2]`: each root is `ProvenCoincident` on all three spheres through `side_of(position)`; on the plane through the origin with normal `(0, 0, 1)` the first is `Front` and the second `Back`; on the sphere centred at `(3, 3, 0)` with radius 3 both are `Back` (`|p - c|^2 = 7 < 9`).
21. Budget: under `BudgetPolicy { max_precision_bits: 64, max_work: 10_000 }`, every point of cases 18 to 20 gives either its default-budget outcome or `NotExact`. A `NotExact` bound is at least the true distance from the point to the carrier (for case 18 at `e = 10^-30`: `5e` for the cylinder, `3e` for the sphere, `e` for the plane with normal `(0, 0, 2)`). A point off the carrier never comes back `ProvenCoincident`, a point on it never comes back `Proven`, and no point comes back on the wrong side.
22. A `Decider` whose deadline has passed: `side` on a point with no fact gives `Refused("numbers.deadline")`; `side` on a point settled by a fact still gives `ProvenCoincident` with `ConstructionFact`.
23. `side_of` on the sphere with centre `0` and radius 1, at a closed position whose `x` is `1 / (sqrt(2) sqrt(2) - 2)`, gives `Refused("numbers.division-by-zero")`, passed up from piece 1 unchanged. Whether an undefined coordinate that the implicit polynomial does not use is refused is not part of the contract; tests do not ask it.

### Points on curves and same points

24. The line with origin `(1, 2, 3)` and `dir` `(2, -1, 0)` has equations `[z - 3, 2z - 6, -x - 2y + 5]`. `(5, 0, 3)` (at `t = 2`) is `ProvenCoincident`; `(5, 0, 3 + 10^-30)` and `(5 + 10^-30, 0, 3)` are `Proven(Apart)`. The line with origin `0` and `dir` `(1, 0, 0)` has equations `[z, -y]`.
25. The circle with centre `0`, normal `(0, 0, 1)`, `radius_sq` 25 and `ref_dir` `(1, 0, 0)` has equations `[z, x^2 + y^2 + z^2 - 25]`. `(3, 4, 0)` is on it; `(3, 4, 10^-30)` and `(3, 4 + 10^-30, 0)` are `Apart`. A circle with `radius_sq` 2 contains `(1, 1, 0)`; one with `radius_sq` `2 + 10^-30` does not.
26. The ellipse with centre `0`, normal `(0, 0, 1)`, `major_dir` `(2, 0, 0)`, `major_sq` 25 and `minor_sq` 9 has equations `[z, x^2/25 + y^2/9 - 1]`. `(3, 12/5, 0)`, `(5, 0, 0)` and `(0, 3, 0)` are on it; `(3, 12/5 + 10^-30, 0)` is `Apart`.
27. Intersection curves: the curve of the sphere with centre `0` and radius 5 and the plane with origin `(0, 0, 3)` and normal `(0, 0, 1)`, with explicit witness `(4, 0, 3)`: `(0, 4, 3)` and `(-4, 0, 3)` are `ProvenCoincident` by numbers; `(0, 4, 3 + 10^-30)` is `Apart`. A point added as a `Meet` of both surfaces is `ProvenCoincident` by fact. Branch membership is not this piece's question: the curve of the cylinder of radius 1 about the `z` axis and the plane `x = 1/2` (normal `(1, 0, 0)`), with closed witness `(1/2, sqrt(3)/2, 0)`, contains `(1/2, -sqrt(3)/2, 5)` for this piece, although that point is on the other line of the meeting.
28. `same_point`: explicit `(1/3, 0, 0)` and closed `(Real(1) / Real(3), 0, 0)` are `ProvenCoincident`; explicit `(1, 1, 0)` and closed `(sqrt(2) sqrt(1/2), 1, 0)` are `ProvenCoincident`; explicit `(1 + 10^-30, 1, 0)` and closed `(sqrt(2) sqrt(1/2), 1, 0)` are `Apart`; the root of `[x^2 + y^2 + z^2 - 1, x - y, z]` certified in `[1/2, 1] x [1/2, 1] x [-1/2, 1/2]` and closed `(sqrt(1/2), sqrt(1/2), 0)` are `ProvenCoincident`. Adding explicit `(1, 2, 3)` twice gives two handles that are `ProvenCoincident` by numbers; a point with itself is `ProvenCoincident` by fact.
29. Under the tiny budget of case 21, cases 24 to 28 give their default outcome or `NotExact`, and a `NotExact` bound is at least the true distance: `10^-30` for `(3, 4, 10^-30)` from the circle of case 25; about `8e-31` for `(3, 4 + 10^-30, 0)` from the same circle. For the intersection curve of case 27, the bound is at least the larger of the true distances from the point to the sphere and to the plane (`10^-30` for `(0, 4, 3 + 10^-30)`).

### Same surfaces and normals

30. The plane with origin `(0, 0, 1)`, normal `(0, 0, 1)`, `ref_dir` `(1, 0, 0)` and the plane with origin `(5, 7, 1)`, normal `(0, 0, -4)`, `ref_dir` `(0, 1, 0)` are `ProvenCoincident` with `ExactArithmetic`; moving the second origin to `(5, 7, 1 + 10^-30)` makes them `Apart`.
31. The cylinders with origin `0`, axis `(0, 0, 1)`, radius 2 and with origin `(0, 0, 5)`, axis `(0, 0, -3)`, radius 2 are `ProvenCoincident`. Origin `(10^-30, 0, 0)`, radius `2 + 10^-30`, or axis `(10^-30, 0, 1)` with `ref_dir` `(1, 0, -10^-30)` each makes them `Apart`.
32. The cones with apex `0`, slope `1/2` and axes `(0, 0, 1)` and `(0, 0, 2)` are `ProvenCoincident`; with axes `(0, 0, 1)` and `(0, 0, -1)` they are `Apart` (the other nappe); with slopes `1/2` and `1/2 + 10^-30` they are `Apart`.
33. Spheres with the same centre and radius but different `axis` and `ref_dir` are `ProvenCoincident`.
34. Tori with axes `(0, 0, 1)` and `(0, 0, -7)` and equal centres and radii are `ProvenCoincident`; a different `minor` makes them `Apart`.
35. Different kinds are `Apart`, even when they touch along a whole curve: the cylinder of radius 1 about the `z` axis and the sphere of radius 1 at the origin. `same_surface` changes no statistics except on the fact path, and answers after the deadline.
36. `normals_agree`: the sphere of radius 2 at the origin and the plane with origin `(0, 0, 2)` and normal `(0, 0, 1)`, at `(0, 0, 2)`: `Proven(Positive)`; with normal `(0, 0, -1)`: `Proven(Negative)`. The planes `x = 0` (normal `(1, 0, 0)`) and `y = 0` (normal `(0, 1, 0)`) at `(0, 0, 5)`: `ProvenCoincident`. The cone of case 12 and any plane, at the apex: `ProvenCoincident`.

### Evaluation and parameters

The `Ball` inputs below are exact rationals or multiples of `Ball::pi`. Tests check containment of the exact answer, never closeness.

37. The cylinder with origin `0`, axis `(0, 0, 1)`, radius 2 and `ref_dir` `(3, 0, 0)` (not unit): `evaluate(0, 3)` contains `(2, 0, 3)`; `evaluate(pi/2, 3)` contains `(0, 2, 3)`; `evaluate(pi, -1)` contains `(-2, 0, -1)`. `derivatives(0, v)` contains `[(0, 2, 0), (0, 0, 1)]`.
38. The sphere with centre `(1, 2, 3)`, radius 2, axis `(0, 0, 1)`, `ref_dir` `(1, 0, 0)`: `evaluate(u, pi/2)` contains `(1, 2, 5)` for every `u`; `evaluate(0, 0)` contains `(3, 2, 3)`; `derivatives(0, 0)` contains `[(0, 2, 0), (0, 0, 2)]`.
39. The torus of case 12 with `ref_dir` `(1, 0, 0)`: `evaluate(0, 0)` contains `(4, 0, 0)`, the outer equator on the seam; `evaluate(0, pi)` contains `(2, 0, 0)`; `evaluate(pi/2, pi/2)` contains `(0, 3, 1)`; `derivatives(0, 0)` contains `[(0, 4, 0), (0, 0, 1)]`.
40. The cone of case 12 with `ref_dir` `(1, 0, 0)`: `evaluate(0, 10)` contains `(5, 0, 10)`, and `evaluate(u, 0)` contains the apex for every `u`. The cone with apex `0`, axis `(0, 0, -1)`, slope `1/2` and `ref_dir` `(1, 0, 0)` has `e2 = (0, -1, 0)`, so `evaluate(pi/2, 2)` contains `(0, -1, -2)`.
41. The plane with origin `(1, 1, 1)`, normal `(0, 0, 1)` and `ref_dir` `(0, 2, 0)` has `e1 = (0, 1, 0)` and `e2 = (-1, 0, 0)`, so `evaluate(1, 1)` contains `(0, 2, 1)`.
42. Curves: the circle with centre `0`, normal `(0, 0, 1)`, `radius_sq` 4 and `ref_dir` `(0, 2, 0)`: `evaluate(0)` contains `(0, 2, 0)` and `evaluate(pi/2)` contains `(-2, 0, 0)`, which is counter-clockwise seen from `+z`. The circle with `radius_sq` 2 and `ref_dir` `(1, 1, 0)`: `evaluate(0)` contains `(1, 1, 0)`. The ellipse of case 26: `evaluate(0)` contains `(5, 0, 0)` and `evaluate(pi/2)` contains `(0, 3, 0)`. The line of case 24: `evaluate(2)` contains `(5, 0, 3)`. An intersection curve gives `None`.
43. `parameters_of`: on the cylinder of radius 5 about the `z` axis with `ref_dir` `(1, 0, 0)`, `(3, 4, 7)` gives `u` containing `atan2(4, 3)` and `v` containing 7, and `(-5, 0, 0)` gives `u` containing `pi`; `(5, 0, 0)`, on the seam, gives `None` or a `u` containing 0. On the torus of case 12, `(0, 2, 0)` gives `u` containing `pi/2` and `v` containing `pi`. A sphere's pole may give `None`; a plane never does. For random exact points from the rational parametrisations above, posed by random poses and at least `L / 100` from the sets of the liveness rule, the result is `Some`, and each ball overlaps the enclosure of the true parameter computed with `Ball::atan2` (and `Ball` arithmetic) from the exact coordinates in the local frame. For random `u` in `[0, 2 pi)` and `v` in the domain, `parameters_of(evaluate(u, v))` is `None` or contains `(u, v)`.
44. `parameter_origin` is exact: the circle of case 42 with `ref_dir` `(0, 2, 0)` gives a point `same_point` with `(0, 2, 0)`; the circle with `radius_sq` 2 and `ref_dir` `(1, 1, 0)` gives `(1, 1, 0)`; the ellipse of case 26 gives `(5, 0, 0)`; a line gives its origin; an intersection curve gives its witness position.

### Store and facts

45. Handles are issued from 0 in order; `len`, `handles` and `get` agree with the order of `push`.
46. `add_curve` of an intersection curve on `[s0, s1]` records `Incident(curve, s0)` then `Incident(curve, s1)`; `add_point` of `Meet` of `Surfaces([s0, s1, s2])` records three incidences in that order; `CurveSurface(c, s)` records `Incident(p, c)` then `Incident(p, s)`; `Curves(a, b)` records `Incident(p, a)` then `Incident(p, b)`. Each origin is `Construction("shapes.definition")`. Analytic curves, surfaces and explicit points record nothing.
47. Fact paths: for a point added as a `Meet` of spheres `a` and `b`, `side(a, p)` is `ProvenCoincident` with `ConstructionFact`, and it adds exactly one `decision`, one `coincident_fact` and no work. After `Incident(p, c)` and `Incident(c, s)` are recorded, `side(s, p)` is settled by fact; after `Coincident(s, s2)` is recorded, `side(s2, p)` is too. `same_surface(a, a)` is settled by fact.
48. `position` returns the defining coordinates; `approx_point` of explicit `(1/3, 0, 0)` is `[to_f64_nearest(1/3), 0.0, 0.0]`; for meeting points each value lies in the coordinate's `Real::approx()`.
49. `enclose_point` with `max_width = 2^-100` on the root `(3, 3, sqrt(7))` of case 20, under the default budget: three balls of width at most `2^-100`, containing 3, 3 and `sqrt(7)` (checked against `sqrt(7) = 2.6457513110645905905016157536392604257102...`).

### Scale, placement and determinism

50. Scale (DESIGN.md, R9): cases 16 to 36 repeated with every point, origin, radius and length multiplied by `10^k` for `k` in `-30..=30` (squared radii by `10^(2k)`, slopes and directions unchanged) keep their outcome kinds.
51. Placement: the same cases moved by 20 random rational poses (integer quaternion components up to 100, translations up to `10^6`) keep their outcome kinds.
52. Determinism: repeating any decision with a fresh `Decider` gives the same outcome and the same statistics. Running the same list of decisions on eight threads, each with its own clone of the store and its own `Decider`, gives identical outcomes and statistics.

### Volume and running time of the tests

The generated families (cases 7, 8, 13, 18, 43, 50, 51 and 52) must together run at least 30,000 maths-answered questions per test run, from fixed seeds.
Keep the default test run of this piece under about 10 minutes in release mode on one core. Questions whose inputs are rational are cheap. Hidden zeros with square roots took about 0.3 to 0.9 ms each in piece 1's build (comparison 01), so use at most a few hundred of those per run; run the larger scale and pose sweeps on rational inputs.

## Values checked for this description

The architect checked these with exact rational arithmetic and 120-digit decimals in an independent program:

- the matrices of cases 2 to 5, the Hamilton product `q2 q1 = (-16, -12, 20, 10)` of case 7, and that `from(q1).then(&from(q2))` equals its matrix;
- the implicit values of case 12, the torus point `(54/25, 72/25, 4/5)`, the ellipse point `(3, 12/5, 0)`, and the moved cylinder and point of case 14;
- case 19: the distance of `(sqrt(1/2) + 10^-30, sqrt(1/2), 0)` from the unit cylinder is about `7.071068e-31`, with `f` about `1.414214e-30`;
- the distance conversions for the torus and the cone, on 2,000 random rational points each, against the exact distance; the ellipse conversion against a dense parameter scan.

## Speed: measured, never gated

No test fails for speed. A benchmark target measures and prints:

1. time per `side` on random rational points in general position, by surface kind, and the share settled in piece 1's filter;
2. time and precision reached for `side` on the near points of case 18 at `e = 10^-30`, and on the square-root and certified-root points of cases 19 and 20;
3. time of a fact-settled decision, and the share of decisions settled by facts in a mixed workload;
4. time of `point_on_curve`, `same_point` and `normals_agree` per curve and surface kind;
5. time of `evaluate`, `derivatives` and `parameters_of` per kind at 64 and 1,000 bits;
6. time of `implicit`, `transformed` and `Rotation::from_quaternion` per kind, and of `enclose_point` at widths `2^-64` and `2^-1000`.

The numbers are written to the benchmark output so later runs can be compared.

## Out of scope for this piece

- Intersecting surfaces and curves with each other, finding meeting points, choosing branches, and ordering points along curves (piece 4).
- Whether a point on a cone's carrier lies on the nappe; callers test the sign of `(p - apex) . axis`.
- Faces, edges and every other topology (piece 3).
- Free-form surfaces, fillet and offset surfaces. They come later as new variants of the non-exhaustive enums.
- Any tolerance. This piece has none.
