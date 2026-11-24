## 2026-11-24: face-bounded extrusion recipes

Pieces: recipe, features, app, judge interface
Change: `catallax-recipe/6` adds `extrude_from`. It enables selected planar and Bézier trims, cylinder cap crossings and circular-profile enclosure proofs for the B-rep `extrude_to` and `extrude_next` paths.
Why: a sketch can start on an existing trimmed face or use it as a limit. Versions 3 through 5 keep their original extrusion paths and work charges.
Migration: retain the saved format to keep prior evaluation behavior. Upgrade to version 6 before authoring `extrude_from` or using the new limit paths. See [the recipe specification](../RECIPE.md#412a-blind-extrusion-from-a-face).
