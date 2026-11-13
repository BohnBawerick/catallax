# Known SMLib misses in the M9 gate

This is the R1 exception list for the intersection port at SMLib
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`. Matching the pinned C++ output
on a listed case does not count as a correct geometric answer. The differential
gates check candidate parity. The separate regression tests below check the
known exceptions, and Catallax's exact proof layer still decides acceptance.
Raw evidence retained privately.

## Missing components and contacts

| Case | Pinned SMLib result | Required Catallax result and evidence |
| --- | --- | --- |
| Perpendicular segments of lengths `1e-10` and `1e-8`, meeting at their exact midpoint | Success with no candidates because the closest-line helper rejects the small cross product | The exact control points share the point at `t1 = t2 = 1/2`. Eight `known-miss` rows cover both curve orders and flags. [Line-pair witnesses](../crates/catallax-trace/tests/curve_curve.rs). |
| One requested line interval collapsed to `[1/2,1/2]` at a common midpoint | Success with no candidates because one evaluated segment has zero length | The closed interval contains the meeting. Eight `known-miss` rows check both curve orders and flags. [Line-pair exceptions](../crates/catallax-trace/tests/fixtures/curve-curve/README.md#r1-known-misses-on-line-pairs). |
| Hill cut 1e-6 mm below the peak at `225/256` | No curve, missing the closed loop | Complete plane-section search finds one loop, 12 arc and 76 empty cells, depth 21, at 174,443 work units. [Plane proof record](M10-PLANE-SECTIONS.md#verification). |
| Hill cut exactly at the peak | Success with no curve, omitting the isolated contact | The plane proof names the touch point `(5/2, 5/2, 225/256)` in its tangent refusal at 60,998 work units. An empty candidate list cannot prove emptiness. [Scope case table](M25-INTERSECTION-PORT-SCOPE.md#historical-proof-at-scale-check). |
| Nearly parallel finite line crossing a plane between its ends | The 0.0001-degree parallel test can report no meeting | Keep the crossing when the ends have opposite signs. [Plane negative gate](../crates/catallax-trace/tests/plane.rs), `near_parallel::a_near_parallel_interior_crossing_is_found`. |
| Short line direction against a surface-normal polar box | `HasPerpendicularToVector` can discard a true perpendicular direction because its sign test uses the unnormalized vector | Normalize before the sign test; a direction too short to normalize cannot prune. [Surface-cache gate policy](../crates/catallax-trace/tests/fixtures/surface-cache/README.md#polar-boxes). The saved differential vectors do not expose this difference. |
| Near-parallel line through a large finite box | A nonzero direction component below `SM_EFF_ZERO` can disappear from the slab test | Keep the crossing. The regression uses `2^44` and a `2^-44` direction component, with crossing interval `[2^42, 3*2^41]`. [Line-clip negative gates](../crates/catallax-trace/tests/line_clip.rs). |
| Left-handed pseudo-box basis | Ordering crossings by basis handedness can discard a true meeting | The reflected unit-box case crosses at `[1,2]`. [Line-clip negative gates](../crates/catallax-trace/tests/line_clip.rs), `a_left_handed_pseudo_box_is_crossed_where_the_line_meets_it`. |

## Invalid answers and classifications

| Case | Pinned SMLib result | Required Catallax result and evidence |
| --- | --- | --- |
| Reparameterization to a single value | Succeeds with an empty active knot interval in 24 saved cases | Return `Fault::InvalidInput`. [Construction exceptions](../crates/catallax-trace/tests/fixtures/curve-build/README.md), `reparameterisation_to_a_single_value_is_refused`. |
| Refinement with descending knot inputs | Succeeds with decreasing knots in 16 saved cases | Return `Fault::InvalidInput`. [Construction exceptions](../crates/catallax-trace/tests/fixtures/curve-build/README.md), `refinement_by_unordered_knots_is_refused`. |
| Missing second solution object or ray-fire ordering on a leaf | Reads a null object or missing child | Return a typed refusal. [Solver behavior tests](../crates/catallax-trace/tests/solver_behaviour.rs). These paths are not credited as successful oracle calls. |
| Curve/surface endpoint fallback outside a requested sub-domain | Keeps a point found in the natural domain, such as `(0,0.75,0.5)` outside the requested half-plane | Drop the out-of-domain answer. [Curve/surface gate exceptions](../crates/catallax-trace/tests/fixtures/curve-surface/README.md). |
| Local surface/surface side search outside requested domains | Keeps meetings from the whole iso-curve outside the requested parameter rectangle | Keep both parameter pairs inside the requested domains. Four saved cases translate the invalid C++ answer to no answer and retain the original privately. [Surface/surface gate exceptions](../crates/catallax-trace/tests/fixtures/surface-surface/README.md). |
| Failed evaluator in a degenerate curve/surface or local surface/surface solve | Can report a found answer despite the evaluator failure | Return the typed failure. [Curve/surface](../crates/catallax-trace/tests/curve_surface.rs) and [surface/surface](../crates/catallax-trace/tests/surface_surface.rs) regression tests. |
| Quartic whose knot samples lie on a parameter line but whose interior bends away | `DropIsoCurve` accepts the straight line from knot samples alone | Require the whole-span control-point bound before taking the iso-line shortcut. Otherwise use the general drop. [Drop-curve negative gate](../crates/catallax-trace/tests/drop_curve.rs), `review::a_curved_curve_is_not_dropped_as_a_parameter_line`. |
| Saddle `z = u*v` with parallel normals | Reads `N.S_uu` and `N.S_vv` but omits `N.S_uv`, so calls the point planar and the meeting tangent | Classify the mixed-curvature saddle as singular. [Intersector negative gate](../crates/catallax-trace/tests/intersector.rs), `review::a_mixed_curvature_saddle_is_not_planar`, in both surface orders. |
| Singularity at a start point with a previous direction | Can keep the first branch found even when it leaves the domain | Select the inside branch. [Intersector negative gate](../crates/catallax-trace/tests/intersector.rs), `review::a_start_singularity_takes_an_inside_branch`. |

The three formerly recorded flat-patch crashes were probe faults. The old
probe dereferenced the absent UV curves returned for natural-boundary grazing
answers. On identical inputs, the old probe exits with signal 11 and the
null-safe probe returns success with one curve in all three cases. They are
restored to the candidate gate and do not count as SMLib misses.

This list records the misses and invalid results already identified by the
port and M25 study. It is not a completeness claim about SMLib. A crash,
unsupported path or `SolverRefusal::NotPorted` is not a parity success.
Numeric ULP allowances compare candidates only; they never establish a proof.
The existing tolerances, deterministic work allowances and 477 saved replay
answers stay unchanged.

## Port status and proof boundary

Main already contains steps 1 to 9, the range-answer follow-ups, analytic
polynomial two-point lines and analytic circular arcs. In particular,
[analytic lines](https://github.com/BohnBawerick/catallax/pull/498),
[circular arcs](https://github.com/BohnBawerick/catallax/pull/508) and
[Hermite and line polar boxes](https://github.com/BohnBawerick/catallax/pull/563)
are merged. This follow-up adds line/circle analytic pairs, rational
two-point lines passing the sampled line test, natural-boundary grazing
curves and recovery through singular side points.

Remaining typed refusals include a curve intersected with itself, a circle
SMLib rebuilds through `MakeNurb`, a zero tangent where SMLib reads past its
array, the projected side search for curves not known to lie on the surface,
closed-surface seam starts and the general tracer fallback. Two of the
132 end-to-end drop cases retain `SolverRefusal::NotPorted` for the projected
side search and tracer fallback. They are not hidden or counted as completed
geometric answers. A ported candidate unit never replaces Catallax's exact
proof layer, and this list makes no claim that every SMLib operation is ported.
