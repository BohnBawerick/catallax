## 2026-10-24: keep construction-frame elimination private

Pieces: 4

Change: the univariate line solver keeps the selected carrier's construction
frame instead of importing a possibly larger witness field. Public
`curve_form` line origins remain the regular witness required by piece 4.
Component-domain, cone-nappe, contact and ordering checks are unchanged.

Why: an offset along a known line can inflate polynomial coefficients without
changing the line or its intersections.

Migration: none. Public forms, signatures, work allowances and tolerance rules
are unchanged. The earlier entry, “reduce validated line carriers before
intersection elimination”, also describes piece 4, not piece 5.

