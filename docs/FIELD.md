# Analytic field kernel

`catallax-field` compiles `catallax-recipe/3` stock solids, line/circle/arc sketch profiles, signed extrusions, two-sided through-all cuts, full and rational-frame partial revolutions, and regularised union, subtraction and intersection. It calls no B-rep builder. Sketch validation calls `catallax_features::Profile::new`, so both kernels use the same intersection, tangency and nesting rules and refusal codes. Every parameter and rigid placement retains its exact definition from [RECIPE.md](RECIPE.md).

This is original MIT code. Fidget supplied design ideas, notably fixed-order tapes, independent point evaluation and reusable evaluator storage. No Fidget code or dependency is included. Regional pruning, meshing, fillets, field-only features, FEA and application integration are not implemented.

## Shared algebraic definitions

The evaluator uses the 2026-10-22 shared algebraic foundation, including the existing profile-sweep consumers. Rational definitions retain the existing arithmetic path. This migration adds no recipe operation or elliptic-cylinder metric coverage; the public distance/result distinctions below continue to apply.

Store the same exact `Surface`, `RealPoly`, frame and dimension definitions as B-rep. Widen polynomial constraints, gradients and cached ray coefficients together. Keep the rational sign path for rational coefficients; use the shared decider for algebraic signs, zero coefficients and support comparisons. Normalising a polynomial by its first syntactic coefficient is invalid if that expression might be zero. Certify the first nonzero coefficient before division. Equal support equations do not remove nappe or trim constraints. Definition keys may permit safe work sharing, but unequal expression keys cannot prove different sets.

Retain exact coefficients, reuse compiled rational enclosures, and enclose algebraic coefficients under the caller's policy for proven queries. Keep separate rounded caches for fast queries. Preserve input order and complete result-bit equality across thread counts. The fixed 192-bit baseline is an initial enclosure precision, not an equality proof or a universal accuracy promise. Algebraic predicates need deterministic refinement/work accounting. An unfinished proof returns the existing unresolved query outcome or a specific propagated failure, never `Boundary` from an interval containing zero. Compilation must validate defining invariants before publishing a plan. Checked compilation and query entrypoints must accept the shared deterministic policy; convenience entrypoints may use its documented default. Do not introduce a separate field-only proof budget.

Algebraic rigid transforms preserve metric distance. An affine stretch used to turn an elliptic cylinder into a circle does not. New elliptic-cylinder/extrusion queries need a certified closest-point calculation on their actual trimmed boundary, or conservative distance bounds with unresolved upper distance. A raw implicit polynomial has no automatic 1-Lipschitz property: either prove a Lipschitz bound on the query region and scale the residual bound accordingly, or return a zero distance lower bound. Do not label that residual `TrueDistance`, and do not reuse the round-solid metric formula for oblique caps.

Add cross-kernel fixtures for 45-degree revolve, extrusion along `[1,0,1]`, oblique circle extrusion and algebraic-radius arcs as both evaluators gain those operations. Compare exact supports, trims, orientations and origin sets, then certified metric enclosures. A mesh sample agreement remains a display check and cannot certify 1 nm.

## Query contract

| API or result | Meaning |
|---|---|
| `Plan::compile`, `compile_checked` | Check the recipe, construct exact supporting surfaces, then compile a separate evaluation tape. The checked variant uses the caller's decider; unsupported operations and impossible dimensions return named refusals. |
| `Plan::fast` | Tentative scalar display/filter value. `true_distance_formula` identifies the underlying analytic formula, not a floating-point accuracy certificate. |
| `Plan::proven`, `proven_checked` | Exact rational/algebraic predicates, outward-rounded binary-ball evaluation and explicit result semantics. Input coordinates are exact `RPoint3` values in millimetres; the checked variant propagates numeric refusals under the caller's policy. |
| `Plan::transformed` | Compose a validated exact rigid pose into supports, constraints, metric placements and face records using the caller's decider. |
| `Value::TrueDistance` | Enclosure of true signed Euclidean distance to a primitive boundary. Negative inside. |
| `Value::Residual` | Enclosure of a sweep support residual or min/max Boolean expression. It is not Euclidean distance. |
| `DistanceBound` | Bounds on unsigned geometric distance to the result boundary, separate from the field value. Infinity means unbounded, not a measured tolerance. |
| `Location` | Proven interior, exterior or boundary; `Indeterminate` when an exact trim or the local Boolean proof cannot settle the query. |
| `Direction` | Componentwise enclosure of a unique outward unit normal at a certified boundary point, `Ambiguous` at a sharp or unresolved configuration, `Undefined` elsewhere. |

Convenience compilation and queries use `Decider::default()`. If a convenience proven query cannot finish a numeric proof, it returns `Indeterminate`, an entire residual interval, unbounded distance and an ambiguous direction. `proven_checked` and `face_query_checked` instead preserve the specific numeric refusal. No unresolved sign becomes a boundary classification.

The empty solid has no boundary. Its unsigned distance is positive infinity and its location is `Outside`. A Boolean zero is not automatically a boundary. In particular, two blocks meeting only on a plane have an empty regularised intersection, even though their max residual is zero on the contact plane. Compilation proves that intersection empty when opposite coincident primitive halfspace constraints exclude any common interior. Exact expression identities also fold empty results. The contact plane is never a face of that empty result.

Compilation supports every listed Boolean combination without a geometric construction failure. That does not make this a complete classifier for all degenerate regularised zero sets. `Indeterminate` is an unresolved query, not a geometric refusal, a tolerance-based coincidence or a proven answer. Callers must not promote it to one.

## Why the distance bounds hold

For a box, let `q_i = max(-p_i, p_i - size_i)`. Its signed distance is the norm of the positive components of `q`, plus `min(max(q_i), 0)`. Thus a box corner uses Euclidean distance, not the maximum plane residual.

A sphere uses `norm(p) - radius`. A ring torus uses `sqrt((sqrt(x²+y²) - major)² + z²) - minor`.

For a capped round solid, reduce to `(rho, z)`, where `rho = sqrt(x²+y²)`. The meridian side is the segment joining `(base_radius, 0)` and `(top_radius, height)`. The nearest side point is its clamped orthogonal projection. The cap candidates are the nearest points on the two radial disk segments. The minimum of these three distances is the true unsigned boundary distance, including points beyond cap rims and cone apices. Exact cap and side predicates supply its sign. Equal radii give the same definition as a capped cylinder. Rotational symmetry justifies the meridian reduction, and exact rigid motions preserve distance.

Each primitive signed distance is 1-Lipschitz. Minimum, maximum and negation preserve that property, so the compiled Boolean residual `f` is also 1-Lipschitz. Every genuine boundary point `q` of the regularised result satisfies `f(q) = 0`, even though regularisation may remove other zeros. Consequently `abs(f(p)) <= norm(p-q)`. Taking the infimum gives a proven lower bound on distance to the result boundary.

For an enclosed residual `[lo, hi]`, the lower distance bound is zero if the enclosure contains zero, otherwise `min(abs(lo), abs(hi))`. The Boolean upper distance bound is positive infinity. We do not compose primitive upper bounds through booleans: hidden or removed source faces make that unsound. At a certified boundary point the distance is exactly zero. Primitive queries return both finite distance bounds when representable as doubles. An overflowed endpoint remains an honest outward infinity.

These lower bounds support certified exclusion and filtering. They do not certify proximity, clearance below a threshold, offsets, or a 1 nm hand-off. Such uses need finite upper bounds from certified closest points on surviving trimmed faces. There is no claim that a small residual proves closeness.

### Profile sweep certificates

Sweeps use conservative residuals, even when a particular profile could reduce to a stock solid. Each supporting plane, cylinder, sphere or torus contributes its unsigned distance to the untrimmed surface. A cone contributes `abs(rho - slope * abs(z)) / sqrt(1 + slope²)` in its unit-axis frame. This is a 1-Lipschitz lower bound to the double-cone support; it does not select a nappe or replace the profile trim. Take the minimum over all supports and give it the exact solid-membership sign. Every genuine sweep boundary belongs to at least one support, so this signed residual is continuous and 1-Lipschitz, including across a change of membership. It can vanish on inactive supports inside or outside the solid. Its absolute value is a certified distance lower bound, never evidence of proximity. The upper distance bound is infinity except at a certified boundary point.

Extrusion projects a query along its defining vector back into the sketch plane, then tests exact odd-loop membership and the signed blind interval. Through-all omits both interval caps and subtracts that infinite tool from the target. It has no construction caps or arbitrary extent. Revolution projects to the meridian using the exact square root of squared radial distance, tests both radial sides of the profile, and applies the cardinal angular sector for partial turns. Profile boundary predicates retain circular arcs rather than replacing them with chords. Shared `Decider` signs settle algebraic meridian and ray-crossing expressions; an exhausted exact decision produces `Indeterminate`, not a guessed sign.

Full turns omit cap names and axis segments that sweep no area. At 180 degrees, caps connected along an axis segment expose both `start` and `end` names over the merged domain. Consecutive collinear or cocircular segments expose every merged segment origin over the maximal side domain. Face queries check that domain as well as the supporting equation, including the correct branch of a revolved surface.

The supported rational subset and `features.unsupported-*` refusals match [RECIPE.md section 7.14](RECIPE.md#714-b-rep-feature-geometry-bugs). Legal algebraic placements, non-cardinal cap frames, elliptic cylinders and generated irrational surface radii remain explicit unfinished geometry, not impossible inputs.

## Exact signs and surface directions

Each primitive is an intersection of the exact polynomial inequalities of its supporting surfaces and caps. Tori use the ring-torus quartic, not a squared equation interpreted without its ring-radius restriction. A cone's cap inequalities select the intended nappe. Rational evaluation settles rational primitive signs, including exact boundary points, without a numerical tolerance; algebraic coefficients use certified shared decisions.

For a Boolean point on a source boundary, strict interior/exterior logic settles unambiguous cases. Remaining cases use exact polynomial ray germs. Substituting `p + epsilon*d` into each support polynomial yields degree at most four. The first nonzero coefficient gives its sign for every sufficiently small positive epsilon. Stock-solid predicates use that coefficient directly. Sweep trims evaluate at an exact rational displacement proved to stay in every support polynomial's initial sign cell: if `c_k` is the first nonzero coefficient, choose `0 < epsilon <= min(1, abs(c_k)/(abs(c_k) + sum_{j>k} abs(c_j) + 1))` across supports. The leading term then strictly dominates the tail throughout that interval. This displacement is a certified germ, not a sampled geometric tolerance. A ray confined to a true boundary is not counted as interior.

When every vanishing constraint of every boundary leaf is planar, the classifier covers the complete local plane arrangement. Nonzero constraints keep their signs near the point. Coincident planes are one plane, regardless of normal orientation. One plane has two open sides; multiple planes divide the direction sphere into chambers. Each chamber meets a vertex formed by two nonparallel planes. Certified perturbations to both sides of every incident line visit all adjacent sectors at both ends of that vertex axis. Each perturbation preserves the first nonzero coefficient's sign under the shared decider. Interior and exterior chambers prove boundary; uniformly interior or exterior chambers prove `Inside` or `Outside`. This removes buried meeting faces and cancelled common supports, including their edges and corners, without interpreting a zero residual as boundary.

For other configurations, a fixed order of active face normals and 26 rational coordinate directions searches for both an interior and exterior germ. Finding both proves the point is a genuine boundary point. Failure to find both proves nothing and returns `Indeterminate`. In particular, outside straight rays alone cannot remove a curved cusp.

Normals come from exact supporting-surface gradients, not finite differences. Interior/exterior germs orient them in the final body, including subtraction's reversal. Distinct active directions, zero gradients, cap rims and corners return `Ambiguous`, never a tie-selected branch normal. Normalisation uses certified binary square-root and division enclosures.

## Face names and trimming

`Plan::faces` is an immutable source-face registry outside the optimised tape. It holds each recipe `Origin` and exact world-space `Surface`. `select` resolves an origin, not an instruction index. Geometrically identical leaves may share evaluation work without sharing or losing their names. A union of identical solids keeps both sets of selectable origins. Subtraction can fold to empty without deleting its source names from the registry.

Registry entries are candidates, not a claim that all their geometry survives. `face_query` checks the supporting equation, source primitive domain, result classification and orientation. It returns `NotOnFace` for removed or out-of-trim points. Distinct active supporting surfaces at a Boolean trim intersection can return `Indeterminate`; mere contact with a support is not proof of a surviving face patch. Primitive edge/corner incidence remains selectable with an ambiguous normal. Coactive names on a common surviving support provide the merged face's origin set. This crate does not construct maximal B-rep faces, edge graphs or split components.

A `FaceId` is local to the plan that returned it. Persist the `Origin`, not this handle, and resolve it after rebuilding.

## Whole-face definition measurements

`catallax_field::measurement::measure_face` takes a B-rep `Model`, a topology
`FaceId`, an optional `FaceDefinition`, `MeasurementOptions` and the caller's
`Decider`. This checker is separate from `Plan` and its model-wide field.
`FaceDefinition::Fillet` holds a `CentrePath` and rolling radius.
`FaceDefinition::Offset` holds the original analytic supporting surface and
unsigned offset distance. A missing record declares an exact face and returns
`[0, 0]`. Callers retain these records alongside their named faces; the checker
does not infer a feature definition from the finished B-rep.

For a fillet, the checked function is `abs(distance(p, centre_path) - radius)`.
The two outward-offset rolling supports define the centre path. Points on the
fillet face are not themselves at the rolling radius from both supports.
For an offset the function is `abs(distance(p, original_support) - distance)`.
These are Euclidean distances to the declared definitions, not polynomial
residuals, Boolean min/max values, display chord errors or distances to the
finished part's boundary.

`FaceMeasurement.lower <= sup(face deviation) <= FaceMeasurement.upper`.
The upper endpoint is the face's local tolerance in millimetres. The record also
reports `converged`, live covering `cells` and exact trimmed-face `witnesses`.
`converged` means the enclosure width met the requested width, not that the
geometry met a manufacturing target. A 100 nm target is `0.0001` mm.

### Whole-face proof

Certified traces of every boundary arc bound the linear chart coordinates of
compact planar, cylindrical, conical and elliptic-cylindrical faces. A linear
chart coordinate cannot have a strict interior extremum on these supports.
Periodic coordinates cover a full turn. Spheres and tori use complete supporting
charts. This rectangle covers the whole face, including boundaries and holes.
It may include trimmed-away regions, so its upper bound can be conservative
even after refinement. The report does not claim trim-optimality.

Each parameter cell has a certified centre point and a uniform spatial radius
`E = L_u * width_u / 2 + L_v * width_v / 2`. The constants bound the Euclidean
norms of the chart derivatives everywhere on the covering domain:

| Supporting surface | `L_u` | `L_v` |
| --- | --- | --- |
| Plane | 1 | 1 |
| Circular cylinder | radius | 1 |
| Cone | slope times maximum absolute axial coordinate | `sqrt(1 + slope²)` |
| Sphere | radius | radius |
| Ring torus | major plus minor radius | minor radius |
| Elliptic cylinder | major semi-axis | 1 |

Angles use radians in this formula. The implementation stores turns and encloses
the conversion with a proven pi ball. An orthonormal surface frame preserves
these derivative bounds. Integrating the two derivatives proves that every
point of the cell lies within `E` of its centre.

Distance to a set and `abs(distance - radius)` are 1-Lipschitz. The centre's
outward-rounded deviation enclosure plus `E` therefore bounds every point in
the cell. The maximum cell upper endpoint bounds the whole face. Samples alone
never supply this upper certificate. Exact support identities instead prove a
constant deviation over the whole support, including every possible trim.
They cover rolling cylinders, rolling tori, parallel plane offsets and
concentric circular, spherical and toroidal offsets.

The lower endpoint comes only from exact points on the actual trimmed face.
Refinement projects candidate witnesses onto the support and checks exact trim
membership before using their distance bounds. A point in a hole or on a
trimmed-away support cannot raise the lower endpoint. Refinement visits the
largest upper bound first, keeps a fixed tie order and intersects each child's
upper certificate with its parent's. More cells cannot worsen the bound.

Lines and circles have analytic path distances. `CentrePath::Arc` retains the
model's exact carrier, selected branch and endpoints. Certified trace-piece
boxes cover the entire arc, including non-conic degree-four intersections.
Minimum box-distance lower bounds and distances to certified path witnesses
bound closest distance from both sides. Elliptic-cylinder support distance uses
the same argument after orthogonal projection onto its transverse ellipse.
`path_segments` controls this independent path refinement.

The checker uses 128-bit outward binary balls and exact rational endpoints.
The default requested width is `1e-7` mm, with at most 4096 live face cells and
64 subdivisions per path trace piece. A cell limit returns the attained bound
with `converged=false`, never an invented precision. Non-compact or degenerate
chart covers and invalid controls return `field.face-measurement`; numeric,
intersection and exact trim-query refusals propagate. No B-rep fillet or offset
is constructed, and no tolerance cap or merge policy is applied.

## Determinism and storage

Compilation follows recipe order. Equivalent leaves and elementary Boolean identities are folded deterministically. Unreachable instructions are removed while the face registry remains intact. Evaluation follows the resulting topological tape order. Hash iteration, point reduction, scheduling and random seeds do not change a proven answer. Checked operations honor the caller's wall-clock deadline and may return a numeric refusal; default convenience queries have no wall-clock deadline.

The proven evaluator uses exact rational/algebraic predicates and initial 192-bit `catallax-numbers::Ball` metric enclosures. Algebraic coefficient enclosures and signs use the caller's numeric policy. Endpoints pass through `Interval::from_rational` for outward rounding. It calls no platform square root or transcendental function. Enclosures can widen at extreme scales; there is no hidden 1 nm assertion. Exact predicates still distinguish features smaller than a double's spacing.

The fast evaluator uses fixed scalar arithmetic and a bit-seeded, fixed-iteration software square root. It has no JIT, SIMD-dependent path, fused operation or platform maths call. NaN and infinite inputs are rejected. Arithmetic outside its finite range reports `field.fast-range`; callers can use the proven path. Returned zeros and interval endpoint zeros are canonical positive zero.

`Workspace` reuses tape registers. Fast scalar evaluation allocates nothing after the registers have warmed up. Proven arithmetic still allocates exact integers and binary balls as needed. Polynomial ray coefficients use fixed-size arrays over cached terms.

`Batch` reuses registers and output capacity through `fast_batch` and `proven_batch`. Output index always matches input index. Workers share an immutable `Plan`, each with its own storage. There is no internal thread pool. Empty batches are valid, and a fast batch error exposes no partial result.

## Running and checking

```sh
cargo run -p catallax-field --example query -- part.json 3 2 2
cargo test -p catallax-field -p catallax-recipe
```

The example reads a real recipe, prints tentative and proven answers, resolves named faces and checks scalar/batch agreement. Exact coordinate arguments accept the recipe number grammar.

Acceptance tests cover every stock solid against independent exact `catallax-shapes` surface definitions, known rational and irrational distances, capped and expanding frusta, rigid placements, all booleans, surviving and removed names, reversed tool normals, sharp edges and apex ambiguity, degenerate regularisation and large-coordinate cancellation. Sweep tests cover rectangles, circles, slots, holes and islands, oblique line extrusions, revolved hollow rings and cones, 90/180/270-degree turns, cap merging, infinite through-all cuts, B-rep origin availability and specific shared refusals. Both proven queries and tentative sweep values have complete result-bit equality across one, two and four threads. Cross-architecture execution and representative end-to-end performance have not been measured here. There are no throughput or interactive-budget claims.
