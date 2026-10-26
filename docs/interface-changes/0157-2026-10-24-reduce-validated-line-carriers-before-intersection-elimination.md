## 2026-10-24: reduce validated line carriers before intersection elimination

Pieces: 5

Change: curve–surface intersections use the existing univariate line solver for
validated line carriers against rational as well as algebraic implicit targets.
Selected-component, cone-nappe, contact and hit-ordering checks are retained.
Other carriers keep the existing multivariate solver.

Why: a constructed line already supplies a separating parameter; eliminating
its world coordinates again discards that degree-reduction certificate.

Migration: none. Public signatures, branch acceptance rules and work allowances
are unchanged.

