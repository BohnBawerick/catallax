## 2026-10-24: reuse exact display guards and boundary enclosures

Pieces: 4, display mesh, application benchmark

Change: face charts retain their exact Hessian coefficients. Cones still scale the axial term by the current parameter extent. Edge UV chord guards first compare outward-rounded intervals, then use the unchanged rational predicate when inconclusive. Exact fallback results are local to the chart and keyed by binary64 parameter bounds and the exact threshold. Strict UV orientation uses the same filter/fallback rule. Incident planar faces share each boundary vertex's working-precision enclosure.

Exact traced boundary points evaluate their scalar DAG directly instead of constructing order-zero Taylor coefficient vectors. Shared DAG nodes are evaluated once per parameter, and the original expression order, zero additions and reciprocal definitions are retained.

Why: repeated coefficient construction, rational UV predicates and enclosure evaluation consume the edit display budget. The standalone application display benchmark records build and display time separately and captures complete meshes without instrumenting application or kernel builds.

Migration: none. Default 0.01 mm accuracy, refinement order, vertex and arithmetic allowances, exact geometry, picking handles, manifests and STEP output are unchanged. Caches live only for one mesh request; no answer or enclosure crosses model definitions or precision settings.

