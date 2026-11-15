# Numerical finite cone/cylinder bands

`catallax_intersect::numerical_native_finite_intersection` accepts the original
cone/cylinder carriers in either operand order, a strictly positive increasing
cone-height interval and an increasing cylinder-axis interval. Distances are
millimetres along unit axes from the original analytic origins. This finite
query does not require the global support window to exclude an apex meeting or
a parallel cone generator. It does require the requested cone band to exclude
the apex and opposite nappe.

The result is `NumericalNativeIntersection`. Its original analytic carriers,
rounded positive-weight source definitions, all reported 3D/UV fits, preset
1e-7 mm tolerance, actual checks and stage work remain available. Sixteen pairs
of native quarter rectangles are searched through the existing metered SMLib
intersector. Fits keep their native intervals and knots. Shared chart contacts
are retained in every reported piece. `NumericalIntersectionPiece::source_domains`
identifies each piece's requested source rectangles in operand order. Independent
samples check those rectangles and physical source/carrier residuals. Empty
finite results remain numerical, with no exact emptiness claim.

`NumericalNativeWindows::cylinder_origin` identifies the source frame origin.
The existing global-support query keeps its cone-apex projection origin. The new
finite query uses the original cylinder origin, so the same cylinder and axis
interval produce the same cylinder source for different cone faces. The existing
`cylinder_parameters` getter now documents its coordinates relative to this
explicit origin. No existing global source or fit is reparameterized.

Invalid or collapsed bands refuse before calling the tracer, whose empty-domain
fallback must not turn a finite request into a natural-domain search. The query
shares the caller's operation allowance and deadline. Before tracing, the
converted binary64 control net must preserve distinct axial controls at every
angular pole. A far-origin band whose distinct retained rows collapse during
proposal conversion refuses with `intersect.numerical-native-proposal-collapse`.
The common source origin and native coordinates are not recentered. Successful
queries record this actual proposal check separately from source/carrier checks. The standalone allowance
stays 20M. No numerical fit becomes an exact pair, incidence fact or exact result.
The result does not contain Boolean topology, a complete branch-count proof or
a whole-curve/body deviation bound.

The `numerical_hole_carriers` example accepts `--finite-bands MIN MAX` to run this
query on actual circular cone-band boundaries in E07 recipe prefixes. The caller
supplies one cylinder interval for both faces. Its reports retain original
analytic double inputs for OCCT comparisons and each piece's requested rectangles.
`--check-supports` independently admits every retained quarter as a rational patch.

This follows the existing SMLib requested-domain intersector and OCCT's finite
face intersection ordering. OCCT constructs finite curves after domain boundary
classification in [IntTools_FaceFace::MakeCurve](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/IntTools/IntTools_FaceFace.cxx#L690)
and [GeomInt_LineConstructor::Perform](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/GeomInt/GeomInt_LineConstructor.cxx#L112).
The P04 port map and known-miss list still apply. Returning all reported pieces
and matching independent finite controls do not establish general completeness.
