## 2026-10-23: certified trimmed analytic minimum distances

Pieces: shapes, measurement, field, comparison, verification CLI
Change: `catallax_shapes::distance` exports compact `Domain` covers, `Trim`,
`Membership`, `Options`, `Bounds`, `project`, `point_distance` and `face_distance`.
`catallax_measure::distance::PreparedFace` adapts exact B-rep loops.
`Plan::prepare_distances` returns field-native `distance::Prepared`, with
`point_to_face`, `face_to_face` and `point_to_boundary`. Comparison reports add
`sampled_distance_upper`; verification exposes it as the additive
`comparison.sampled_distance_upper_mm` rational string or null.
Why: curved support residuals alone cannot prove proximity to a surviving trimmed
face. Exact admitted witnesses give finite upper bounds; interval chart covers
prove lower bounds without replacing trims or the whole-face definition metric.
Migration: distance callers must inspect `converged` and optional `upper`, retain
the prepared model or plan's definitions, and supply a certified compact box when
using the shapes-level API. Default query width is 1e-6 mm with 4096 parameter
cells. Invalid covers or controls use `shapes.distance-domain`; unmarked B-rep
bounds and traces are required by `measure.distance-marked`. A field result
without a certified compact recipe domain returns `field.distance-unbounded`;
absent field handles retain `field.face-id`. Numeric and exact trim refusals
propagate. Comparison accuracies, numeric budgets, recipe construction, STEP
output and whole-face definition measurements are unchanged. Previous curved
extrusion, revolution, round-Boolean and small planar-Boolean perturbation
limitation tests now require agreement with finite sampled bounds.
