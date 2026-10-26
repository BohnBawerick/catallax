## 2026-10-25: exact finite ellipse application edge queries

Pieces: application (`catallax_app::curves`)
Change: `CurveSpan` adds `FullEllipse` and `EllipticArc { interior, crosses_seam }`. `edge_curve` now answers ring and coincident-endpoint ellipse edges as `FullEllipse`, and finite ellipse edges as `EllipticArc` with an exact interior point from `arc_point` and a strict seam-crossing flag from `cyclic_order` against the major-axis point. The original definition, ordered endpoints, traversal sense, origins and construction pairs are retained as in [0170](0170-2026-10-25-exact-trimmed-application-edge-queries.md). Line and circle answers are unchanged, unresolved intersection definitions still refuse with their original code, and parabolas, hyperbolas, splines and general branches stay `Unsupported`. This is a query capability, not silhouette or visibility support.
Why: ellipse edges of represented bodies were reported as `Unsupported` although their exact carrier and trim are known.
Migration: callers with an exhaustive `match` on `CurveSpan` add the two variants; callers that treated ellipse edges as `Unsupported` read the new spans.
