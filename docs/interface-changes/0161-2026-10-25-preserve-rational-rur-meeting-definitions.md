## 2026-10-25: preserve rational RUR meeting definitions

Pieces: 4, 6, 7

Change: this corrects “reduce validated line carriers before intersection
elimination” and “keep construction-frame elimination private”. Rational targets
retain the existing RUR separator, primitive-root polynomial and coordinate
expressions. When the carrier basis proves both transverse coordinates affine,
the solver omits redundant finite-coordinate projections. Algebraic targets
retain their existing witness-frame line construction; the private solver form
is removed.

Why: forward substitution through a different exact line parameter could loosen
centroid enclosures and STEP deviation bounds even when the meeting positions
were proven coincident and the STEP bytes were unchanged. Preserving the
published certificate requires preserving the construction, not just its value.

Migration: none. Public signatures, component guards, work allowances and
tolerance rules are unchanged. The earlier line-carrier entry concerns piece 4,
not piece 5.

