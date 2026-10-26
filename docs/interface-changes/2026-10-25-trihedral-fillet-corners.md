## 2026-10-25: exact trihedral fillet corners

Pieces: features, field

Change: `catallax_features::fillet_corner::CornerDefinition` reduces three
convex line strips of equal radius between mutually perpendicular planes to an
exact sphere sector with three quarter-circle boundaries. Both evaluators use it:
the B-rep sews a `corner.vertex` sphere face to the three cylinder strips, and the
field subtracts the corner cube less the ball. Other three-strip junctions now
refuse with the new code `fillet.unsupported-corner` instead of the general
`fillet.unsupported-topology` or `field.fillet-contact` corner refusal.

Why: block corners with all three edges filleted are a defined RECIPE.md section
4.5.1 case with an exact analytic answer.

Migration: callers matching the old three-strip refusal text should match
`fillet.unsupported-corner`. `Built::fillets` is unchanged and does not retain
corner definitions. Recipes, budgets, tolerances and STEP options are unchanged.
