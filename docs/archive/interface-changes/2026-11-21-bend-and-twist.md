Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-21: bend and twist recipe operations

Pieces: recipe, features, application, field
Change: `catallax_recipe::Operation` gains `Twist { body, axis, angle, length }` and `Bend { body, axis, toward, radius }`, written `"op": "twist"` and `"op": "bend"` and accepted only in `catallax-recipe/5` ([RECIPE.md section 4.20](../../RECIPE.md#420-bend-and-twist)). Both build only through `evaluate_incremental_outcome` and the application's no-display numerical outcome, as fitted geometry admitted by `build_with_sampled_carriers`. The exact evaluators refuse them with `deform.requires-numerical-evaluation`. New refusal codes: `deform.definition`, `deform.bend-axis`, `deform.bend-overlap`, `deform.unsupported-face`, `deform.unsupported-edge`, `deform.fit`, `deform.source` and `deform.name`. Neither operation makes faces; references keep naming the input's steps.
Why: roadmap 3.4, tools only SMLib has: a first build of bend and twist, after SMLib's `SmBendVolume` and `SmTwistVolume`.
Migration: `Operation` is `#[non_exhaustive]`, so callers that match it already have a fallback arm. A caller that lists operations which make no faces adds the two. The field kernel refuses both with `field.unsupported-operation`, and the two-kernel comparison reports them as not comparable.
