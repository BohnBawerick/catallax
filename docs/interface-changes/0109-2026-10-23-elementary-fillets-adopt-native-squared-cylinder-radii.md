## 2026-10-23: elementary fillets adopt native squared cylinder radii

Pieces: features, shapes
Change: elementary fillet definitions construct and compare `Cylinder::radius_sq` directly. Metric bounds and torus major-radius calculations use `Cylinder::radius()`.
Why: integrating exact fillets with native squared-radius carriers must preserve both definitions without extracting roots for incidence predicates.
Migration: cylinder consumers use the squared field for algebraic support equations and the positive radius accessor for metric expressions. No recipe schema, budget, tolerance or STEP-option change.

