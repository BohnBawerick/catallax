## 2026-10-23: checked finite spline supports retain their rigid definitions

Pieces: shapes and geometric consumers
Change: `SplineSurface` owns a checked regular convex offset sweep, its exact target and finite axial domain. `SplinePointParameters` provides exact finite incidence, chart tangents and outward normals. Cross-section side queries use a bounded Bernstein interval filter before the existing exact path. Rigid transports retain canonical source data and the composed exact pose; inverse queries reuse the source frame instead of re-solving world-coordinate equations.
Why: separate interval boxes lose parameter correlation in tangent half-planes; a rigid coordinate permutation can also make the first two world-coordinate equations dependent. Neither warrants a larger work allowance or a geometry refusal disguised as impossible input.
Migration: the local support remains distinct from a valid recipe-bound B-rep and its complete region certificate. Cross-section side is not finite-body membership. Uncertain filtering never proves contact; poses preserve source domains and error bounds by exact isometry. Existing analytic kinds and legacy exchange formats remain unchanged.

