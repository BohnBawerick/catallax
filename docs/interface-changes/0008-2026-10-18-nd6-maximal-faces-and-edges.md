## 2026-10-18: ND6 maximal faces and edges

Pieces: 5, judge interface
Change: the kernel promises maximal faces and edges: no edge between two faces on the same surface with the same sense, no vertex joining exactly two edges on the same curve that bound the same two faces, and a ring edge where an edge covers a whole closed curve alone. This entry's original marked-merge rule is superseded by the later `boolean.merge-ambiguous` rule.
Why: the manifest counts then match the counts maths gives, and the pieces merged are proven or marked like every other decision.
Migration: piece 5's sew removes redundant vertices as well as redundant edges.

