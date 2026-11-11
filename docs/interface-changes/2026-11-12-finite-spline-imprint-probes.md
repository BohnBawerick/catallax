## 2026-11-12: finite spline queries for section clipping

Pieces: 2, 4, 5.

Change: Finite rational patches can invert a nonsingular affine pair of world coordinates when every knot span has constant nonzero rational weight. The candidate must lie in its closed span and satisfy the remaining homogeneous coordinate equation exactly. Unsupported frames retain the existing solver. Projected-parameter and side queries retain their definitions.

A nonzero coefficient of a plane's homogeneous pullback proves that a finite spline support differs from the plane. This proves nonidentity, not disjointness. An identically zero pullback retains the existing support-equality refusal.

Native finite B-spline curves can meet rational planes through complete knot-span polynomial roots. Only transverse interior roots produce crossings. The existing Bezier handler retains precedence. The new `intersect.spline-plane` refusal names missing contact proofs at tangencies, knots, endpoints or repeated spatial points, partial-span overlaps, and unsupported coefficients.

Boolean section clipping can obtain an exact interior probe from a retained section branch. Cell witnesses are checked against the selected arc. Otherwise, a rational coordinate plane strictly between the arc endpoints supplies a complete meeting set, filtered by branch membership and arc order. `boolean.patch-section-probe` names an unresolved probe proof.

`RationalPatch::embedding_projection` exposes the patch's existing checked embedding axes. `RationalAtlas::affine_projection` checks whether every chart has the same nonsingular affine world-coordinate pair in original source parameters, after dividing by constant rational weight. Equal projected positions then imply equal source parameters and world points. The returned axes preserve parameter orientation. Unsupported maps return `None`.

Face membership uses certified winding when a spline has such a proved injective projection. Other supports retain their existing query route. `boolean.spline-chart` names a missing global atlas projection when winding needs one. This projection does not change the spline surface's definition.

Accounting: New coefficient arithmetic and rational probe construction reserve the caller's existing exact-arithmetic cost before execution. Exhaustion retains `numbers.operation-budget`.

Migration: Existing signatures, formats, tolerances and work allowances are unchanged. Callers must accept the new refusal codes. These queries do not supply section trace enclosures, complete trimmed-face classification, a sewn loft Boolean or section-branch STEP output.
