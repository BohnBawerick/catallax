## 2026-11-01: planar draft of a body with a countersunk through hole

Pieces: features, field

Change: `draft` now answers for a convex planar body with one complete axial through
hole opened on one cap by a coaxial cone, in both evaluators. The wall, the cone, the
ring they share, both cap openings and every name stay unchanged. `bored` gains `Head`,
`certify_head`, `seam` and `BoredCell::countersunk`, and `BoredCell` gains `head`.

Why: SolidDraft had to move a countersink after the draft or suppress it to draft a
part's walls.

Migration: callers that matched `draft.unsupported-arrangement` or
`draft.unsupported-source-boundary` for this family now receive an answer. `Bore` and
`BoredCell::new` are unchanged. A cone touching or crossing an outer wall, a blind hole,
a counterbore, a countersunk hole beside another bore or a boss, and a selected wall,
cone or cap keep those refusals. Recipes, budgets, tolerances and the answers for a
plain bore are unchanged.
