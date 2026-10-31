## 2026-11-01: planar draft between a preserved boss and bore

Pieces: features, field

Change: `draft` now answers for a convex planar body with one complete outward axial
cylinder standing on a stationary cap and one separate complete axial through bore,
made before the draft in either order, in both evaluators. Both walls, the boss's end,
all four rings and every name stay unchanged, and the named cap keeps both exclusions.
`bossed::BossedCell::new` takes a `Vec<Bore>` and `BossedCell` gains `bored`. The detail
text of the `bored::certify_disjoint` refusals now says "walls".

Why: SolidDraft had to suppress a boss or a bore before drafting a part's walls.

Migration: callers of `BossedCell::new` pass an empty list for a boss alone. Callers that
matched `draft.unsupported-arrangement` or `draft.unsupported-source-boundary` for this
family now receive an answer. A boss touching or crossing a bore, a bore through a boss,
several bores beside a boss, several bosses, and a selected wall, end or cap keep those
refusals. This corrects the boss-and-bore line of
[the preserved boss entry](2026-11-01-draft-preserved-boss.md). Recipes, budgets,
tolerances and the answers for a boss or bores alone are unchanged.
