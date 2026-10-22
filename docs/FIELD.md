# Analytic field kernel

`catallax-field` compiles `catallax-recipe/3` blocks, capped cylinders, spheres, capped cone frusta, ring tori and regularised union, subtraction and intersection. It calls no B-rep builder. Every parameter and rigid placement retains its exact rational definition from [RECIPE.md](RECIPE.md). The new through-all operation and the sketch-based features still report `field.unsupported-operation`; a schema version does not imply every evaluator implements every operation.

This is original MIT code. Fidget supplied design ideas, notably fixed-order tapes, independent point evaluation and reusable evaluator storage. No Fidget code or dependency is included. Regional pruning, meshing, fillets, field-only features, FEA and application integration are not implemented.

## Approved shared-definition migration

The 2026-10-22 architecture decision requires field to migrate with the algebraic foundation, before either kernel claims widened recipe coverage. The rest of this document describes the current rational stock-solid evaluator. Its public distance/result distinctions continue to apply.

Store the same exact `Surface`, `RealPoly`, frame and dimension definitions as B-rep. Widen polynomial constraints, gradients and cached ray coefficients together. Keep the rational sign path for rational coefficients; use the shared decider for algebraic signs, zero coefficients and support comparisons. Normalising a polynomial by its first syntactic coefficient is invalid if that expression might be zero. Certify the first nonzero coefficient before division. Equal support equations do not remove nappe or trim constraints. Definition keys may permit safe work sharing, but unequal expression keys cannot prove different sets.

Compile certified coefficient enclosures for the proven evaluator and separate rounded caches for fast queries. Preserve input order and complete result-bit equality across thread counts. The fixed 192-bit baseline is an initial enclosure precision, not an equality proof or a universal accuracy promise. Algebraic predicates need deterministic refinement/work accounting. An unfinished proof returns the existing unresolved query outcome or a specific propagated failure, never `Boundary` from an interval containing zero. Compilation must validate defining invariants before publishing a plan. Checked compilation and query entrypoints must accept the shared deterministic policy; convenience entrypoints may use its documented default. Do not introduce a separate field-only proof budget.

Algebraic rigid transforms preserve metric distance. An affine stretch used to turn an elliptic cylinder into a circle does not. New elliptic-cylinder/extrusion queries need a certified closest-point calculation on their actual trimmed boundary, or conservative distance bounds with unresolved upper distance. A raw implicit polynomial has no automatic 1-Lipschitz property: either prove a Lipschitz bound on the query region and scale the residual bound accordingly, or return a zero distance lower bound. Do not label that residual `TrueDistance`, and do not reuse the round-solid metric formula for oblique caps.

Add cross-kernel fixtures for 45-degree revolve, extrusion along `[1,0,1]`, oblique circle extrusion and algebraic-radius arcs as both evaluators gain those operations. Compare exact supports, trims, orientations and origin sets, then certified metric enclosures. A mesh sample agreement remains a display check and cannot certify 1 nm.

## Query contract

| API or result | Meaning |
|---|---|
| `Plan::compile` | Check the recipe, construct exact supporting surfaces, then compile a separate evaluation tape. Unsupported operations and impossible dimensions return named refusals. |
| `Plan::fast` | Tentative scalar display/filter value. `true_distance_formula` identifies the underlying analytic formula, not a floating-point accuracy certificate. |
| `Plan::proven` | Exact rational predicates, outward-rounded binary-ball evaluation and explicit result semantics. Input coordinates are exact `RPoint3` values in millimetres. |
| `Value::TrueDistance` | Enclosure of true signed Euclidean distance to a primitive boundary. Negative inside. |
| `Value::Residual` | Enclosure of the min/max Boolean expression. It is not Euclidean distance. |
| `DistanceBound` | Bounds on unsigned geometric distance to the result boundary, separate from the field value. Infinity means unbounded, not a measured tolerance. |
| `Location` | Proven interior, exterior or boundary; `Indeterminate` when a degenerate Boolean zero cannot be settled by the fixed local proof. |
| `Direction` | Componentwise enclosure of a unique outward unit normal at a certified boundary point, `Ambiguous` at a sharp or unresolved configuration, `Undefined` elsewhere. |

The empty solid has no boundary. Its unsigned distance is positive infinity and its location is `Outside`. A Boolean zero is not automatically a boundary. In particular, two blocks meeting only on a plane have an empty regularised intersection, even though their max residual is zero on the contact plane. This implementation reports `Indeterminate` at that degenerate zero unless an exact plan identity already establishes emptiness. It never certifies the contact plane as a face of that empty result.

Compilation supports every listed Boolean combination without a geometric construction failure. That does not make this a complete classifier for all degenerate regularised zero sets. `Indeterminate` is an unresolved query, not a geometric refusal, a tolerance-based coincidence or a proven answer. Callers must not promote it to one.

## Why the distance bounds hold

For a box, let `q_i = max(-p_i, p_i - size_i)`. Its signed distance is the norm of the positive components of `q`, plus `min(max(q_i), 0)`. Thus a box corner uses Euclidean distance, not the maximum plane residual.

A sphere uses `norm(p) - radius`. A ring torus uses `sqrt((sqrt(x²+y²) - major)² + z²) - minor`.

For a capped round solid, reduce to `(rho, z)`, where `rho = sqrt(x²+y²)`. The meridian side is the segment joining `(base_radius, 0)` and `(top_radius, height)`. The nearest side point is its clamped orthogonal projection. The cap candidates are the nearest points on the two radial disk segments. The minimum of these three distances is the true unsigned boundary distance, including points beyond cap rims and cone apices. Exact cap and side predicates supply its sign. Equal radii give the same definition as a capped cylinder. Rotational symmetry justifies the meridian reduction, and exact rigid motions preserve distance.

Each primitive signed distance is 1-Lipschitz. Minimum, maximum and negation preserve that property, so the compiled Boolean residual `f` is also 1-Lipschitz. Every genuine boundary point `q` of the regularised result satisfies `f(q) = 0`, even though regularisation may remove other zeros. Consequently `abs(f(p)) <= norm(p-q)`. Taking the infimum gives a proven lower bound on distance to the result boundary.

For an enclosed residual `[lo, hi]`, the lower distance bound is zero if the enclosure contains zero, otherwise `min(abs(lo), abs(hi))`. The Boolean upper distance bound is positive infinity. We do not compose primitive upper bounds through booleans: hidden or removed source faces make that unsound. At a certified boundary point the distance is exactly zero. Primitive queries return both finite distance bounds when representable as doubles. An overflowed endpoint remains an honest outward infinity.

These lower bounds support certified exclusion and filtering. They do not certify proximity, clearance below a threshold, offsets, or a 1 nm hand-off. Such uses need finite upper bounds from certified closest points on surviving trimmed faces. There is no claim that a small residual proves closeness.

## Exact signs and surface directions

Each primitive is an intersection of the exact polynomial inequalities of its supporting surfaces and caps. Tori use the ring-torus quartic, not a squared equation interpreted without its ring-radius restriction. A cone's cap inequalities select the intended nappe. Rational evaluation settles primitive signs, including exact boundary points, without a numerical tolerance.

For a Boolean point on a source boundary, strict interior/exterior logic settles unambiguous cases. Remaining cases use exact polynomial ray germs. Substituting `p + epsilon*d` into each polynomial yields degree at most four. The first nonzero coefficient gives its sign for every sufficiently small positive epsilon. No finite displacement or sampled tolerance is substituted for that proof. A ray confined to a supporting surface remains unresolved rather than being counted as interior.

A fixed order of active face normals and 26 rational coordinate directions searches for both an interior and exterior germ. Finding both proves the point is a genuine boundary point. Failure to find both proves nothing and returns `Indeterminate`.

Normals come from exact supporting-surface gradients, not finite differences. Interior/exterior germs orient them in the final body, including subtraction's reversal. Distinct active directions, zero gradients, cap rims and corners return `Ambiguous`, never a tie-selected branch normal. Normalisation uses certified binary square-root and division enclosures.

## Face names and trimming

`Plan::faces` is an immutable source-face registry outside the optimised tape. It holds each recipe `Origin` and exact world-space `Surface`. `select` resolves an origin, not an instruction index. Geometrically identical leaves may share evaluation work without sharing or losing their names. A union of identical solids keeps both sets of selectable origins. Subtraction can fold to empty without deleting its source names from the registry.

Registry entries are candidates, not a claim that all their geometry survives. `face_query` checks the supporting equation, source primitive domain, result classification and orientation. It returns `NotOnFace` for removed or out-of-trim points. Distinct active supporting surfaces at a Boolean trim intersection can return `Indeterminate`; mere contact with a support is not proof of a surviving face patch. Primitive edge/corner incidence remains selectable with an ambiguous normal. Coactive names on a common surviving support provide the merged face's origin set. This crate does not construct maximal B-rep faces, edge graphs or split components.

A `FaceId` is local to the plan that returned it. Persist the `Origin`, not this handle, and resolve it after rebuilding.

## Determinism and storage

Compilation follows recipe order. Equivalent leaves and elementary Boolean identities are folded deterministically. Unreachable instructions are removed while the face registry remains intact. Evaluation follows the resulting topological tape order. No hash iteration, point reduction, scheduling decision, random seed or deadline affects an answer.

The proven evaluator uses exact rational arithmetic and `catallax-numbers::Ball` at a fixed 192-bit precision. Endpoints pass through `Interval::from_rational` for outward rounding. It calls no platform square root or transcendental function. Enclosures can widen at extreme scales; there is no hidden 1 nm assertion. Exact predicates still distinguish features smaller than a double's spacing.

The fast evaluator uses fixed scalar arithmetic and a bit-seeded, fixed-iteration software square root. It has no JIT, SIMD-dependent path, fused operation or platform maths call. NaN and infinite inputs are rejected. Arithmetic outside its finite range reports `field.fast-range`; callers can use the proven path. Returned zeros and interval endpoint zeros are canonical positive zero.

`Workspace` reuses tape registers. Fast scalar evaluation allocates nothing after the registers have warmed up. Proven arithmetic still allocates exact integers and binary balls as needed. Polynomial ray coefficients use fixed-size arrays over cached terms.

`Batch` reuses registers and output capacity through `fast_batch` and `proven_batch`. Output index always matches input index. Workers share an immutable `Plan`, each with its own storage. There is no internal thread pool. Empty batches are valid, and a fast batch error exposes no partial result.

## Running and checking

```sh
cargo run -p catallax-field --example query -- part.json 3 2 2
cargo test -p catallax-field -p catallax-recipe
```

The example reads a real recipe, prints tentative and proven answers, resolves named faces and checks scalar/batch agreement. Exact coordinate arguments accept the recipe number grammar.

Acceptance tests cover every primitive against independent exact `catallax-shapes` surface definitions, known rational and irrational distances, capped and expanding frusta, rigid placements, all booleans, surviving and removed names, reversed tool normals, sharp edges and apex ambiguity, degenerate regularisation, large-coordinate cancellation, and complete query-bit equality across one, two and four threads. Cross-architecture execution and representative end-to-end performance have not been measured here. There are no throughput or interactive-budget claims.
