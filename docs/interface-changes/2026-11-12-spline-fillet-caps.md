## 2026-11-12: spline fillet faces cut at their neighbours

Pieces: 2
Change: `RationalPatch::extended` gives a one-span patch's exact map over a larger rectangle. `SplineFilletDefinition::with_borders` takes `RailBorder` lines where the rails must meet their supports' borders, and `SplineFillet::borders` and `SplineFillet::nodes` return the fitted nodes. `SplineFillet::plane_cut` returns a `FilletCut`, or refuses with the new code `fillet.spline-cut`.
Why: a fillet on a loft edge runs past the loft's caps and is cut by them ([prerequisite](2026-11-12-spline-fillet-faces.md)).
Migration: none; the change is additive.
