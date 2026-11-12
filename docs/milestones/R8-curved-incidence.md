# Circular spline repair on capped cylinders

A cubic approximation of a circular STEP edge now imports when its complete neighbourhood consists of one unchanged cylinder and two perpendicular planar caps. The synthetic quarter-circle case previously stopped at measured `spline_curve_surface_incidence`. It now builds one certified body and records one local tolerance. The attained pointwise displacement is enclosed by 0.000000000097 to 0.000850114026 mm. The upper bound is below the unchanged 0.01 mm import cap.

This is a bounded first step. Physical generator edges, larger shells, spline supports, moved vertices, holes and other analytic supports still need their own complete neighbourhood certificate. The refusal names the missing support, endpoint, ring, displacement or feature step. No source model geometry is committed. Raw evidence retained privately.

## Finite correspondence and feature bounds

The source must have exact authored endpoints on the unchanged cylinder and exactly in its cap plane. Charged rational blossoming extracts every active homogeneous Bezier span, including nonclamped knots. Positive weights and exact C0 joins retain the entire source interval.

In one exact rational orthonormal circle frame, write the homogeneous source as `X,Y,0,W`, radius `R`, and `Q = R W + X`. A positive Bernstein hull for `Q` gives one finite stereographic chart. A strictly signed Bernstein numerator for the derivative of `Y/Q` proves a one-to-one correspondence onto the native finite arc. The source polar derivative must have the same strict interior sign. Thus the source boundary also follows the same angular sector without folding or reaching the axis.

The target homogeneous circle is `R(Q²-Y²), 2R QY, 0, Q²+Y²`. Exact Bernstein products bound the Euclidean displacement at every shared parameter. Its upper bound bounds both finite-curve Hausdorff directions. The lower bound measures the maximum chosen pointwise correspondence displacement, not a Hausdorff lower bound. Midpoint subdivision stops at 4096 visited cells, depth 24 or the caller's existing work allowance.

Complete angular ring coverage is required on both caps and the cylindrical band. Matching source edge identities, unchanged vertices and exact circle carriers prevent skipped sectors, overlapping arcs or a different seam being inferred from proximity. Radius, cap separation and distinct ring-vertex chords give a conservative lower bound `L` for neighbouring face widths and edge lengths. Source polar order makes each cap a radial graph. Boundary displacement at most `delta` reduces neighbouring separations by at most `2 delta`. The test `12 delta <= L` therefore proves `10 delta <= L - 2 delta` for the source feature limit. Bounds are refined against the smaller of `0.01 mm` and `L/12` before deciding.

The derivative then passes the existing strict solid builder, intersection checks and topology checks. Each changed physical source edge retains its source geometry hash, measured bounds, approved tolerance and proved source feature bound. Supports and shared vertices remain exact. No work allowance changes.
