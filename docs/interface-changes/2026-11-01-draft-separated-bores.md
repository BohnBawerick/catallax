## 2026-11-01: planar draft around separated through bores

Pieces: features, field

Change: `draft` now answers for a convex planar body with several mutually separated
axial through bores between the same two caps, in both evaluators. The bores may
differ in radius, come from repeated hole or subtraction steps, or come from one step
whose tool cuts several of them. Every bore wall, support, cap ring and name stays
unchanged. `catallax_features::draft::bored` gains `Bore`, `certify_disjoint` and
`certify_bores`. `BoredCell` now holds `bores` in place of `cylinder` and `origins`,
its `caps` are ascending, and `BoredCell::new` takes the bore list.

Why: SolidDraft had to suppress all but one hole before drafting a part.

Migration: callers that matched `draft.unsupported-arrangement` or
`draft.unsupported-source-boundary` for this family now receive an answer. A selected
bore or cap, a bore touching or crossing an outer wall or another bore, bores that
do not share both caps, a blind or non-axial tool, an oblique neutral plane, a
concave envelope and any curved outer wall keep those refusals. Recipes, budgets,
tolerances, zero-bore and one-bore answers are unchanged.
