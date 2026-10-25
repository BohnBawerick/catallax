## 2026-10-18: ND4 touching material is separate lumps sharing no entity

Pieces: 3, 5, 7, judge interface
Change: where a result's material touches itself, the two sides share no face, edge or vertex. Pieces that touch only along curves or at points are separate lumps, each its own `MANIFOLD_SOLID_BREP`; a connected piece touching itself at isolated points is one lump; a connected piece touching itself along a curve is refused with `boolean.self-touching-shell`. This replaces the rule that touching lumps may share edges and vertices.
Why: lumps with their own entities export as self-contained manifold solids, and `result` counts are then the sum of the lumps' counts; a non-manifold solid has no `MANIFOLD_SOLID_BREP` form.
Migration: piece 3 and piece 5 creators never share topology between lumps; piece 7 writes each lump's entities on their own.

