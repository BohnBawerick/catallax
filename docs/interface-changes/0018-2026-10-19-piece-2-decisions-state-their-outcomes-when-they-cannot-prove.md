## 2026-10-19: piece 2 decisions state their outcomes when they cannot prove

Pieces: 2, 4, 5
Change: a `NotExact` from `side`, `side_of`, `point_on_curve` or `same_point` carries a proven distance bound in model units (callers mark it with `BoundKind::Distance`); on an intersection curve it bounds the distance to each of the two surfaces, not to the curve. `point_on_curve` and `same_point` decide equations or coordinates one at a time in a fixed order and stop at the first proven non-zero or refused one. One from `normals_agree` carries the dot-product bound (`BoundKind::Value`). `same_surface` uses exact rational arithmetic without the `Decider` and is never `NotExact` or `Refused`. A fact-settled decision calls only `note_fact_decision` and answers after the deadline. `side` answers about the implicit's zero set, so a cone's opposite nappe is coincident. The decisions take `&self` and never record facts.
Why: piece 1's comparison showed that blind tests need the exact outcome of every unproven case, and marks from these decisions must have a geometric meaning.
Migration: callers of these decisions build marks with the stated `BoundKind`. A caller that needs the cone's nappe tests `(p - apex) . axis` itself. A caller that wants a proven zero kept records it with origin `Proof`.

