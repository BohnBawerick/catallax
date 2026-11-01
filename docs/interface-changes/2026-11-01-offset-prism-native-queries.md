## 2026-11-01: native queries on the checked inward offset prism

Pieces: numbers, shapes, intersect, topology, boolean

Change: `classify_point`, `face_contains` and `intersect_curve_surface` with a line now
answer on the complete checked elliptic-offset prism, including its offset edit and a
placed or rescaled definition. `SplineSupport::line_hits` returns `Some` for
`EllipticOffset`, through the new `SplineSurface::line_hits`. `SplinePrism` gains
`of_body`, `of_face`, `locate` and `locate_on_face`, `Model` gains `bodies_with_face`, and
`RealPoly` gains `simple_root_in`, which certifies a simple common root whose coefficients
lie in one quadratic tower.

Why: the body already built and exported, while 24 of 28 native point queries and every
transverse line refused, so nothing could rely on its material or finite faces.

Migration: callers that matched `intersect.parametric-coincidence`,
`intersect.unsupported-pair` or `intersect.parametric-surface` for this family now receive
an answer. A tangent line and a line overlapping a generator keep
`intersect.parametric-surface`. A face or body that does not complete the checked prism
stays on its earlier path and keeps its refusals. General spline Booleans, spline-curve
equivalence, fits, tolerances, budgets, schemas, STEP bytes and certified output are
unchanged. Certified bodies stay `exact:false`.
