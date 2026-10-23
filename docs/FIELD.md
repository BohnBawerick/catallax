# Analytic field kernel

`catallax-field` compiles `catallax-recipe/3` stock solids, line/circle/arc sketch profiles, signed extrusions, two-sided through-all cuts, full and rational-frame partial revolutions, linear and circular body patterns, plane mirrors, and regularised union, subtraction and intersection. It calls no B-rep builder. Sketch validation calls `catallax_features::Profile::new`, so both kernels use the same intersection, tangency and nesting rules and refusal codes. Every parameter and rigid placement retains its exact definition from [RECIPE.md](RECIPE.md).

This is original MIT code. Fidget supplied design ideas, notably fixed-order tapes, independent point evaluation and reusable evaluator storage. No Fidget code or dependency is included. Definition-native display meshing is described below. General fillets, field-only features, FEA and application integration are not implemented.

## Shared algebraic definitions

The evaluator uses the 2026-10-22 shared algebraic foundation, including the existing profile-sweep consumers. Rational definitions retain the existing arithmetic path. That foundation migration added no recipe operation or elliptic-cylinder metric coverage; the trimmed-distance API below now supplies certified metric bounds on valid analytic definitions.

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
| `Plan::prepare_distances` | Prepare reusable compact analytic chart covers from recipe bounds and trim predicates. No B-rep is built. |
| `distance::Prepared::{point_to_face, face_to_face, point_to_boundary}` | Outward-rounded minimum-distance bounds with exact admitted witnesses, an explicit convergence flag and deterministic refinement limits. |
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

These lower bounds alone support certified exclusion and filtering, not proximity or clearance below a threshold. Prepared trimmed-face queries supply finite upper bounds from certified points on surviving faces. A small residual still does not prove closeness, and neither API constructs offsets or certifies a 1 nm hand-off.

### Profile sweep certificates

Sweeps use conservative residuals, even when a particular profile could reduce to a stock solid. Each supporting plane, cylinder, sphere or torus contributes its unsigned distance to the untrimmed surface. A cone contributes `abs(rho - slope * abs(z)) / sqrt(1 + slope²)` in its unit-axis frame. This is a 1-Lipschitz lower bound to the double-cone support; it does not select a nappe or replace the profile trim. Take the minimum over all supports and give it the exact solid-membership sign. Every genuine sweep boundary belongs to at least one support, so this signed residual is continuous and 1-Lipschitz, including across a change of membership. It can vanish on inactive supports inside or outside the solid. Its absolute value is a certified distance lower bound, never evidence of proximity. The upper distance bound is infinity except at a certified boundary point.

Extrusion projects a query along its defining vector back into the sketch plane, then tests exact odd-loop membership and the signed blind interval. Through-all omits both interval caps and subtracts that infinite tool from the target. It has no construction caps or arbitrary extent. Revolution projects to the meridian using the exact square root of squared radial distance, tests both radial sides of the profile, and applies the cardinal angular sector for partial turns. Profile boundary predicates retain circular arcs rather than replacing them with chords. Shared `Decider` signs settle algebraic meridian and ray-crossing expressions; an exhausted exact decision produces `Indeterminate`, not a guessed sign.

Full turns omit cap names and axis segments that sweep no area. At 180 degrees, caps connected along an axis segment expose both `start` and `end` names over the merged domain. Consecutive collinear or cocircular segments expose every merged segment origin over the maximal side domain. Face queries check that domain as well as the supporting equation, including the correct branch of a revolved surface.

The supported rational subset and `features.unsupported-*` refusals match [RECIPE.md section 7.14](RECIPE.md#714-b-rep-feature-geometry-bugs). Legal algebraic placements, non-cardinal cap frames, elliptic cylinders and generated irrational surface radii remain explicit unfinished geometry, not impossible inputs.

### Patterns and mirrors

Each copy composes an exact inverse isometry into the seed leaves' metric and trim evaluators. Supporting surfaces, constraints, gradients, display caches and names move together. A circular pattern uses `k * angle` about the authored axis, not a total angle or a rounded quaternion. Reflection uses the rational plane formula for normals of any nonzero length. Isometries preserve primitive distance formulas and sweep residual bounds; a union of copies still has Boolean residual semantics.

Copy zero keeps its origins. Later pattern copies and mirror images wrap each source origin with their creating step and `copy.instance`. Every copied subgraph reads a snapshot of its seed's face candidates, independent of evaluation-leaf sharing. The registry retains removed candidates, but selectable ownership follows the actual result body. Thus coincident copies keep all authored names, while a symmetric mirror without its original exposes only the image's wrapped names.

Feature-pattern editing saves ordinary operations through `catallax_app::pattern::expand_feature`; the field compiler needs no feature-pattern opcode. The rational candidate-sampling APIs still report their algebraic-placement limitation when a transformed candidate has no native rational coordinates. That limitation does not round the field's geometry or change its proven query semantics.

### Elementary fillet regions

Plane-plane line contacts reduce to exact cylinder strips. A plane meeting a perpendicular cylinder at a ring reduces to a torus strip. The field resolves the source supports from recipe origins and classifies their material sector with exact positive-infinitesimal rays. It constructs no B-rep. Each rolling-ball section retains its circular arc, then extrudes or revolves the enclosed edge region. Convex regions subtract material; concave regions add it.

Finite straight runs use independently resolved source-edge intervals and perpendicular end planes. Every requested reference must resolve a sharp input component. The blend origin uses the input edge name, not the radius or a generated handle. The retained definition includes both supports, their orientations, the selected carrier, radius, material side and contact rails; rigid transforms move that definition with the field.

The sweep residual and conservative distance bounds above still apply. This is not a new true-distance formula. Filtered references, algebraic section coordinates, oblique line-strip ends, horn or spindle tori, neighbour roll-over and general corner patches remain explicit unfinished coverage. An unproved contact continuation returns `field.fillet-contact` or `fillet.unsupported-overflow`, not an impossible-geometry claim. Independent P3 point and face queries cover convex removal, concave addition and blend origins. Application display and sampled whole-part comparison remain gated by periodic trimmed-face meshing.

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

## Trimmed minimum distances

`catallax_shapes::distance` bounds minimum distance on planes, cylinders, cones,
spheres, tori and elliptic cylinders. Rational quarter-circle charts cover seams
and poles without an angle approximation. Frames and dimensions may be algebraic.
Elliptic distances use the actual orthonormal frame and semi-axes, not a circular
distance after an affine stretch. Cones retain only their declared nappe.

Every live parameter cell has an outward-rounded spatial box. The minimum box
separation is a lower certificate. Exact carrier projections also provide support
lower bounds, with a scaled Lipschitz residual for the ellipse. Upper certificates
come only from exact points admitted by the trim adapter. The search refines low
distance cells and admitted witness cells so outside cells touching a trim edge
cannot indefinitely prevent better witnesses on its surviving side.
For ellipse semi-axes `a >= b`, the function
`b * (sqrt((x/a)^2 + (y/b)^2) - 1)` is 1-Lipschitz in its orthonormal transverse
frame and vanishes on the support. Its absolute value is therefore a distance
lower bound, not an ellipse distance formula.

`catallax_measure::distance::PreparedFace` uses analytic face bounds, certified
boundary-loop traces and exact `face_contains`, including holes and closed edges.
A connected cell with an excluded centre is discarded only if its whole spatial
cover avoids every boundary loop. The field adapter instead uses its own recipe
profiles, support constraints and Boolean ray germs. If only one source leaf is
on a boundary, all other source memberships are locally constant. Evaluating the
Boolean tape with that leaf inside and outside proves whether its face survives,
without rebuilding radial roots for normal germs. Joint boundaries retain exact
ray proofs. Interval constraint signs prove constant trim membership over cells.
Buried tangent sources cannot supply witnesses at isolated contacts; nearby
open-face witnesses approach their surviving closed trim boundaries.
Neither adapter uses display triangles as its distance proof.
Only face entries owned by the compiled result receive distance domains.
Reusing a coincident evaluation leaf does not reactivate discarded copy names;
such registry-only handles report `field.face-id` from the distance API.
Canonical public names retain the plan's generated-domain role for trim
predicates; renaming a revolved cone cap does not turn it into an angular cap.

`Options::width` requests an enclosure width, not geometric acceptance.
`max_cells` limits refinement. Exhaustion returns the attained bounds with
`converged=false`; `upper=None` means no admitted witness was established.
`face_to_face` measures the minimum distance, not Hausdorff distance. These queries
do not replace `measurement::measure_face` or its whole-face definition metric.

Analytic query coverage does not add missing recipe emitters. Oblique circular
extrusion still reports `features.unsupported-elliptic-cylinder`. The elliptic
distance regressions construct exact B-rep circle loops and an independent
native profile/sweep definition directly; they check the Euclidean metric and
the excluded extension beyond the end trim without enabling that recipe feature.

The sampled comparison uses prepared field bounds when the residual alone has
no finite upper, and prepared B-rep bounds when a support projection misses its
trim. `Report::sampled_distance_upper` includes forward mesh error. It is finite
only when every checked distance has a finite upper. Agreement remains sampled;
the comparison accuracy, mesh controls and numeric budgets are unchanged.
After a required trimmed witness settles proximity, later mesh vertices try
incident witnesses before repeating the signed field query. This changes proof
order, not acceptance: every early upper still comes from an exact surviving
boundary point, and a miss retains the original distance and gap checks.


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

## Definition-native display

`catallax_field::display::mesh_plan` consumes a compiled `Plan`, `Options` and
the caller's `Decider`. It never constructs or meshes a B-rep. The returned
`Mesh` contains binary64 positions, triangles tagged with field `FaceId`s,
visited-cell work and a proven two-sided `max_deviation` in millimetres.
`Mesh::write_obj` writes the display artifact. It is not an exportable solid,
and adaptive trim cells and differing refinement levels need not sew watertight.

`mesh_plan_reported` performs the same work with non-charging `Diagnostics`.
The report survives a refusal and records full-request work per shared support,
per-face query costs, cell counts and maximum depth. It is an executed prefix,
not partial geometry. Unreported meshing allocates no diagnostic arrays.
`FaceId::index` maps a display handle to this plan's registry; it is not a
persistent name. See [VERIFY-INTERFACE.md](VERIFY-INTERFACE.md#field-display)
for support sharing and cost aggregation. Work counters require an active
`with_operation_budget` scope and are zero outside it.


The default tolerance is 0.01 mm, with at most one million distinct vertices
and two million visited parameter cells. Plane, cylinder, cone, sphere and
ring-torus charts use rational unit frames and dimensions. Unsupported frames
and surfaces return `field.display-unsupported`. A result without a certified
compact domain returns `field.display-unbounded`; invalid controls return
`field.display-options`. Limits and unresolved refinement return
`field.display-budget`, with visited cells and maximum depth per support.
Coordinate range or rounding failures return `field.display-precision`.
No refusal switches to the B-rep.

### Coverage and deviation proof

Primitive dimensions, profile bounds, sweep extents and exact rigid transforms
give a result bounding sphere. Union takes the larger finite bound,
intersection takes either available bound, and subtraction keeps its target's.
The display covers every candidate support within that sphere. Coincident
supports share a chart cover, so coactive recipe names do not duplicate area.

Each source leaf also supplies a conservative spatial box from its dimensions,
profile and sweep. Cells outside every coactive source box are rejected before
trim queries. These are defining bounds, not fitted boxes or sampled extents.

Project these boxes onto a planar chart's two axes, or a cylinder chart's
height axis, before subdivision. The hull of the coactive source projections
still covers the entire named support domain. Every leaf also uses its own box
to prove regional absence, including leaves unrelated to the candidate face.

Chart coordinates are enclosed directly with outward binary64 interval
arithmetic and cached rational-frame enclosures. Quarter-circle cosine and
sine use their exact monotonicity on each parameter interval, including the
sphere's `[-1,1]` latitude chart. An exact centre point is constructed only
when a membership or boundary witness requires it.

Circular coordinates use rational quarter-circle charts
`((1-t²)/(1+t²), 2t/(1+t²))`. Their speed is at most 2 and their acceleration
norm at most 8. Products give uniform first- and second-derivative bounds for
each supported surface. For a parameter rectangle with spans `du,dv`, linear
triangle interpolation differs from its exact chart by at most
`(Huu*du² + 2*Huv*du*dv + Hvv*dv²)/8`. Exact conversion residuals add the
binary64 vertex-rounding error. These bounds cover whole triangles, not samples.

Planar caps bounded by coaxial circular supports and parallel planes use
four radial quarter-circle charts. A conservative source-box L1 radius covers
their entire domain. A squared-radius cylinder may have a rational nonunit
radial frame. For its squared norm `n` and L1 norm `s`, `n/s` is a positive
lower metric bound and `s` is an upper bound. Exact rational norms replace
both bounds when available. The parameter extent divides the source radius
by the lower bound; physical chord bounds multiply by the upper bound.
A constraint proved equal to a coaxial cylinder restricts
to `A*r²+B`; exact radial endpoint signs classify its open cells, and rational
roots become subdivision cuts. This equality is established from the complete
support definitions, not fitted from samples.

For a radial cell, the emitted triangles form the chord trapezoid of an
annular sector, or a fan when the inner radius is zero. Every point of that
trapezoid has the form `r*((1-s)*c(u0)+s*c(u1))` for an intermediate radius
and angle fraction; the converse construction maps each sector point into the
trapezoid. The quarter-circle acceleration bound therefore gives two-sided
set distance at most `r_max*du²`, where `r_max` bounds the physical radius,
including radial-frame scale. No radial subdivision is required merely
to reduce this bound. This is a set-distance certificate, not an interpolation
claim at matching chart parameters. Exact coordinate-rounding residuals are
added as for other charts. Degenerate fan triangles are not published.

Other sign-stable cells split according to curvature contributions rather
than physical edge length, so straight cylinder generators need not subdivide
with their curved direction. Unsettled trims retain witness-radius refinement.

Outward polynomial ranges determine which leaf constraints can change sign
over a chart cell. Identical supporting constraints remain exactly zero.
Exact leaf membership, field face membership and result-oriented normal queries
settle sign-stable cells. A sample is used only after this regional proof.

For a rational support reference point `a`, cache the exact polynomial
`Q(d) = P(a + d)` once and evaluate it on outward enclosures of `box - a`.
This preserves every sign and zero while removing the dependency cancellation
of expanded terms such as `x² - 2*a*x + a²`. A nonrational reference uses
the original coordinates. Coincidence checks and exact point predicates keep
the original definitions; no surface or trim is moved by this evaluation frame.


Boundary queries reuse regional strict interior/exterior classifications,
which remain constant under a sufficiently small ray from the cell's interior
witness. Only boundary or unresolved registers need the field's exact
infinitesimal ray predicates. Opposite proven interior/exterior rays establish
a genuine boundary point and orient its normal. No metric calculation on
unrelated leaves is needed to prove that boundary point.
For an unresolved leaf, exact strict membership at the witness can likewise
settle both infinitesimal directions without a polynomial ray expansion.
This witness-local state is discarded at the next cell and never promoted to
a regional trim or coverage proof.

A leaf's proved location is cached with its open parameter rectangle.
Descendants reuse it only when their open interiors lie inside that rectangle.
This reuses interval proofs and sign-stable membership, never unproved point
samples. Fixed traversal order makes this regional simplification independent
of worker scheduling.

On a planar chart, exact coincident support identities also define two open
normal half-regions. When every other support of a leaf has a fixed nonzero
sign on the open cell, its two infinitesimal material sides are constant there.
Certify that pair with the exact ray predicate and cache it with the rectangle.
The evaluation tape represents both side predicates as reduced ordered Boolean
diagrams. Variable order is operand identity, then germ kind. Strict material
states are constants; unresolved predicates are variables.

An unresolved operand shares one variable across the two sides only when every
support polynomial has an exact, nonzero restriction to the candidate plane.
The complete affine pullback is computed in rational arithmetic. Its zero set
has no open planar patch, so any remaining boundary there belongs to that
operand's own carrier charts. A zero or unproved restriction retains separate
germ variables; coincident membership is never equated without a proof.

Union, intersection and subtraction combine canonical predicates, including
negation for a cut. Equal canonical side expressions prove that this plane has
no exposed open patch, even when the common material predicate varies over the
cell. This prunes the patch without claiming that its entire region is inside
or outside. Nodes and apply results are shared across neighbouring cells in a
request-local arena. Construction and apply operations check the unchanged
shared numeric work policy.

Opposite constant side expressions prove a wholly exposed planar patch.
Unequal unresolved expressions still require trim refinement or a certified
witness; inequality alone is not an exposure proof.
The chart already lies exactly on the common plane, so the split is by its
two infinitesimal normal directions, not by a rounded spatial offset.
Closed-inside membership alone is insufficient: an exterior cap also has
`f = 0` on its exposed patch. No `sup(f) <= 0` result alone removes
a patch. Shared cell edges retain the closure coverage described below.

Display witnesses for rational extrusions of polygon loops and full circles
use native rational projection, exact odd-loop membership, cap and segment
tests. Their ray displacement uses the same all-support coefficient-dominance
proof as the shared exact-real feature predicate. These caches borrow the
defining loops and do not change field query APIs. Arcs, revolutions and
nonrational placements retain the shared exact-real predicates.

A plane constraint on a planar chart is affine. Its four exact corner values
prove a strict sign on the open rectangle when they are weakly of one sign
and not all zero. Zero endpoint pairs classify aligned edges exactly.
On cylindrical and radial cap charts, exact Bernstein signs of the angular
derivative prove when the projected circular coordinate is monotone on a cell.
Height or radial position is affine, so the same corner argument then applies.
Cells spanning an unproved angular extremum retain interval refinement.
No arbitrary near-zero value is treated as an aligned edge.
Axis-aligned affine trim roots become exact chart cuts, rather than being
approached indefinitely by midpoint subdivision. The same rule applies to
axial plane trims of cylinders.

Open cell interiors are disjoint. Their shared edges have no area, and the
closure of an exposed open patch covers its limiting edge points without
adding another area patch. A regularised solid's boundary is the closure of
its exposed two-dimensional support patches. Thus rejecting an absent open
cell does not lose an isolated trim edge, and emitting the closed triangles
does not double count area at an aligned cell edge.

For an unsettled trim cell, enclose its distance from the exact chart centre
by `R = sqrt(sum_i max(abs(lo_i-centre_i),abs(hi_i-centre_i))²)`, rounded
outward. Refine until `R` is at most half the requested tolerance. When that
centre is a certified surviving boundary point, publish it as a mesh vertex
and triangulate the cell as a fan.

Every true boundary point in the cell is within `R` of this mesh vertex.
Conversely, every fan triangle lies in the convex hull of corners and centre,
each within `R` of the true boundary witness. Convexity of Euclidean balls
therefore bounds the whole displayed fan, not just its vertices. Vertex
rounding adds its exact residual as before. Using this radius rather than
the full cell diameter sharpens the proof without changing the requested
tolerance or its half-tolerance rounding reserve.

A cell without such a witness may be omitted only when
an existing mesh vertex is proved within the requested tolerance of its
entire spatial cover. Otherwise it is refined or refused. This also prevents
an unseen small component from silently disappearing. The union of accepted
cells and certified vertex covers gives the two-sided boundary distance bound.
This display certificate is distinct from whole-face feature-definition
measurements and does not establish a 100 nm fillet construction tolerance.

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
