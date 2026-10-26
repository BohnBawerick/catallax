## 2026-10-26: nearest-face draft selection on a hollow annulus

Pieces: draft, field draft
Change: a saved `near` point on a draft face reference now selects among the four trimmed faces of a coaxial annulus (two cylindrical walls, two perpendicular annular caps) in both kernels, instead of refusing with `draft.unsupported-near` or `draft.unsupported-reference`. Cap distance clamps the radial projection to the annulus, so the hole does not count as part of the cap; wall distance clamps height to the finite interval. Ties remain `recipe.ambiguous-reference`. A selected cap reaches the existing annular construction, which refuses it as before.
Why: close the remaining gap in the annular draft ledger without adding drafted geometry.
Migration: none. Non-coaxial walls, equal-radius walls and other mixed shells keep their nearest-face refusals.
