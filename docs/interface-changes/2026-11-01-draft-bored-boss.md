## 2026-11-01: planar draft of a body with a bored boss

Pieces: features, field

Change: `draft` now answers for a convex planar body with one complete outward axial
cylinder standing on a stationary cap and one complete axial bore, newer than that
cylinder, strictly inside its footprint and passing through both it and the base, in
both evaluators. Both walls, the boss's annular end, all four rings and every name stay
unchanged. `bossed::BossedCell::through` certifies this arrangement, `BossedCell` gains
`through`, and `bored` gains `certify_nested`, `axis_inside` and
`BoredCell::clear_of_bores`.

Why: SolidDraft had to suppress a bored boss before drafting a part's walls.

Migration: callers that matched `draft.unsupported-arrangement` or
`draft.unsupported-source-boundary` for this family now receive an answer.
`BossedCell::new` is unchanged and still proves the bore apart from the boss. A bore
touching or crossing the boss's wall, a blind or stepped bore, a bore older than the
boss, several bored bosses, and a selected wall, end or cap keep those refusals. This
corrects the bore-through-a-boss line of
[the boss and bore entry](2026-11-01-draft-boss-and-bore.md). Recipes, budgets,
tolerances and the answers for a boss or bores alone are unchanged.
