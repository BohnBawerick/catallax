## 2026-10-31: exact skew fillet corners

Pieces: features, field

Change: `catallax_features::fillet_corner::CornerDefinition` no longer requires
mutually perpendicular supports. Three convex line strips of equal radius between
independent planes reduce to the sphere about the three inward offsets, with one
sector in the positive cone of the contact normals. Its boundaries are the shorter
great-circle arcs between those normals, so skew corners have arcs that are not
quarter circles. The definition gains `vertex` and `dual_basis`, and `in_sector`
now tests the cone rather than the sign of each normal's dot product. Both
evaluators use it: the B-rep sews the sector to the three cylinder strips, and the
field subtracts a six-plane cell less the ball in place of the corner cube. The
field also accepts an oblique planar end on a finite line strip, widening the
strip's span to the furthest contact rail station.

Why: an extruded triangle or any other convex polygon with every edge filleted is
a defined RECIPE.md section 4.5.1 case with an exact analytic answer. Its volume
has a closed form, and the app's existing multi-edge fillet action can now reach it.

Migration: corners with dependent normals, mixed concave strips, unequal radii or
nonplanar supports, and corners whose ball centre lies at or past the vertex along
an edge, still refuse with `fillet.unsupported-corner`. Callers that matched that
code for a skew convex corner now receive an answer. Recipes, budgets, tolerances
and STEP options are unchanged, and perpendicular corners keep their frame and
export encoding.
