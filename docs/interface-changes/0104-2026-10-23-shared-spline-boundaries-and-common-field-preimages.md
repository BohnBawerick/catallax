## 2026-10-23: shared spline boundaries and common-field preimages

Pieces: numbers, shapes and field
Change: `RealPoly::gcd` computes a certified monic univariate gcd without eliminating the existing algebraic coefficient field. `BSplineCurve` adds finite positive-weight curve data, shared `SplineBoundary` views, interval jets and checked finite point preimages. `SplineCurveParameter` retains its verified span. The whole-domain sweep checker now lives in shapes and returns an opaque `OffsetSweepProof`; the field measurement API delegates to it without changing its public report.
Why: real fitted edge controls are algebraic. Recovering their parameter through a global rational elimination exhausted the fixed request allowance; common-field coordinate factors recover actual curve incidence without manufacturing a global implicit.
Migration: exact preimages are not uniqueness or body proofs. Surface-boundary identity is a positive shared-data witness, not a negative geometric comparison. Publication still requires valid shared trims, orientation, complete body correspondence and source scope. No existing `Surface`/`Curve` kind, legacy request or exchange output changes in this addition.

