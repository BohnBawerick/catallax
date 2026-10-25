## 2026-10-18: ND5 export deviation, B-spline branches, no pcurves

Pieces: 7, judge interface
Change: `max_deviation` is at most `curve_tolerance`, default 1e-7 mm, unless the answer carries `step.tolerance-not-reached`, which piece 7 also adds when rounding at `significant_digits` alone exceeds the tolerance. Branches with no closed form are one non-rational `B_SPLINE_CURVE_WITH_KNOTS` each. Edges carry no pcurves.
Why: a pcurve is a second approximation per face that would need its own proven distance, and readers build it from the 3D curve and the analytic surface.
Migration: piece 7 writes no `PCURVE`, `SURFACE_CURVE` or `SEAM_CURVE`, and checks rounding against the tolerance.

