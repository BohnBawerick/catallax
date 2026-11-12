## 2026-11-12: finite spline queries for section clipping

Pieces: 2, 4, 5.

Change: Finite rational patches can invert a nonsingular affine pair of world coordinates when every knot span has constant nonzero rational weight. The candidate must lie in its closed span and satisfy the remaining homogeneous coordinate equation exactly. Unsupported frames retain the existing solver. Projected-parameter and side queries retain their definitions.

A nonzero coefficient of a plane's homogeneous pullback proves that a finite spline support differs from the plane. This proves nonidentity, not disjointness. An identically zero pullback retains the existing support-equality refusal.

Native finite B-spline curves can meet rational planes through complete knot-span polynomial roots. Only transverse interior roots produce crossings. The existing Bezier handler retains precedence. The new `intersect.spline-plane` refusal names missing contact proofs at tangencies, knots, endpoints or repeated spatial points, partial-span overlaps, and unsupported coefficients.

For a native spline and a rational line, a nonzero homogeneous pullback of a line constraint proves that an off-line spline point exists. This proves nonidentity even when the two curves share endpoints. Containment in the line retains the missing finite-equivalence refusal.

`native_spline_line_points` supplies every isolated meeting in either input order. It solves a complete univariate cut on every knot span, clips roots to the finite domain, verifies both independent line constraints exactly, and deduplicates shared-knot points. Tangencies remain points without a crossing claim. Unsupported coefficients return `intersect.spline-line`; an overlap interval returns `intersect.spline-line-overlap` with its span and the missing ownership proof. The caller receives no partial point set on failure. Existing Bezier queries and known coincidence facts keep precedence.

`BSplineCurve::proves_control_hull_separation` reuses the patch-pair engine's exact closed-box exclusion certificate for native curves. Positive weights keep the complete curve in its Cartesian control hull. Rational control boxes must be strictly disjoint; touching boxes and unsupported controls provide no certificate. Curve comparison consumes a successful certificate as nonidentity. It makes no equality claim when the hulls overlap.

`native_spline_normal_acceleration` computes a native spline's outgoing normal acceleration at a unique finite preimage. It includes both weight derivatives, selects the finite knot side from the supplied nonzero velocity, and checks both ends of a closed seam. A zero knot-side speed or unproved selector returns `intersect.spline-germ`.

At an exposed clamped endpoint, a strict control-coordinate certificate excludes every interior preimage. Other points retain the complete polynomial preimage query. Positive algebraic weights are supported; weak or mixed control bounds retain the existing solver.

`PatchSectionBranch::normal_acceleration` computes plane-section acceleration from the retained homogeneous patch derivatives. It checks every retained chart preimage before selecting outgoing cell sides, including closed seams. The plane equation and constant-speed condition determine the parameter acceleration and are verified exactly. Unequal selected accelerations, a singular frame or an unsupported cutter return `shapes.patch-section-germ`. Nonplane cutters name the missing exact Hessian query. Atlas regularity and joins do not establish global embedding.

The full Boolean can reuse the existing complete apart-body proof when a fitted operand first needs a native spline germ. Its private certificate binds the unchanged input bodies and is paid once from the caller's allowance. Each speculative worker owns its context; only the main call carries a proof to the Boolean fallback. Standalone imprint queries retain their public signature and run their own geometry queries.

Boolean section clipping can obtain an exact interior probe from a retained section branch. Cell witnesses are checked against the selected arc. Otherwise, a rational coordinate plane strictly between the arc endpoints supplies a complete meeting set, filtered by branch membership and arc order. `boolean.patch-section-probe` names an unresolved probe proof.

`RationalPatch::embedding_projection` exposes the patch's existing checked embedding axes. `RationalAtlas::affine_projection` checks whether every chart has the same nonsingular affine world-coordinate pair in original source parameters, after dividing by constant rational weight. Equal projected positions then imply equal source parameters and world points. The returned axes preserve parameter orientation. Unsupported maps return `None`.

Face membership uses certified winding when a spline has such a proved injective projection. Other supports retain their existing query route. `boolean.spline-chart` names a missing global atlas projection when winding needs one. This projection does not change the spline surface's definition.

When the face-piece search cannot produce a finite spline probe, `boolean.spline-face-probe` names the missing exact parameter-space path and boundary-clearance certificate. Analytic faces retain their existing refusal. This slice does not claim a sewn box-with-loft solid.

Accounting: New coefficient arithmetic and rational probe construction reserve the caller's existing exact-arithmetic cost before execution. Exhaustion retains `numbers.operation-budget`.

Migration: Existing signatures, formats, tolerances and work allowances are unchanged. Callers must accept the new refusal codes. These queries do not supply section trace enclosures, complete trimmed-face classification, a sewn loft Boolean or section-branch STEP output.
