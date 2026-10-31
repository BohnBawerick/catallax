## 2026-11-01: planar draft under a preserved boss

Pieces: features, field

Change: `draft` now answers for a convex planar body with one complete outward axial
cylinder standing on a stationary cap, made before the draft, in both evaluators. The
boss's wall, end disk, both rings and names stay unchanged, and the named cap keeps its
exclusion. `catallax_features::draft::bossed` adds `Boss`, `BossedCell` and
`certify_footprint`, and `catallax_features::draft::bored` adds `ring_clear`. The detail
text of the `bored::check_axial` refusal now says "wall".

Why: SolidDraft had to suppress a boss before drafting a part's walls, or make the boss
after the draft.

Migration: callers that matched `draft.unsupported-arrangement` or
`draft.unsupported-source-boundary` for this family now receive an answer. A boss
touching or crossing the edge of its cap, a selected wall, end or cap, a boss on an
oblique neutral or pull, several bosses, a boss combined with a bore, and a curved or
concave base keep those refusals. The field also keeps its refusal for a tool sunk into
the target. Recipes, budgets, tolerances and undrafted or boss-after-draft answers are
unchanged.
