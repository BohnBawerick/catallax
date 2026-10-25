## 2026-10-22: restrict straight contact quadrics using their source proof

Piece: 5
Change: a verified straight planar probe with proven source contact uses `f(t) = t*c1 + t²*(f(q)-c1)`, where `c1 = grad(f)(q)·p + 2*f(q) - grad(f)(q)·q`. The coefficients are affine in the algebraic source point.
Why: expanding quadratic source-coordinate products reconstructed a zero already supplied by the contact proof.
Migration: only exact carrier checks and existing incidence proofs enable the reduction. Contact polynomials and endpoints use `RealPoly` and `Point3`; the shared decider certifies degree and endpoint contact without rational-only coefficient assumptions. Other paths retain their original expansion, and whole-path separation and component checks remain mandatory. Regressions independently derive simple, double and contained contacts, reject extraneous algebraic carrier candidates, and retain non-rational quadric contact/interior roots. A standalone smoke checks nine coefficients and fifteen exact path values.

