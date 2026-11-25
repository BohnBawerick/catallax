Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-15: one planar side opening in a bored prismatic shell

Pieces: catallax-features, ordinary B-rep shell evaluation.

Change: `shell` accepts one named planar side opening in a convex prismatic envelope pierced by separated parallel cylindrical through bores. The two end caps remain perpendicular to the bore axes, and every side plane must be parallel to those axes. Existing closed and cap-open bored shells retain their established route. No recipe schema, public type, numerical policy or work allowance changes.

Proof: the shared draft bore inventory remains read-only. The shell wall code separately proves the prismatic side-plane condition before constructing a cavity. The convex cavity check still proves that minimum kept-plane projections stay inside their finite planar trims. Cap projections preserve radial coordinates, side supports strictly miss the source cylinders, and radial bore witnesses remain between the perpendicular caps. Expanded cylinders must remain strictly clear of every cavity side plane, including the unshifted opening plane, and of every other expanded bore. Equality and crossing refuse before any model writes; nonpositive thickness also refuses before any writes.

Coverage: the public `bored-side-shell.json` fixture opens `stock.side.x_min` in a 40 × 30 × 12 block with a diameter-6 through bore. Thickness edits from 2 to 1.5 preserve the hole, finite rectangular opening rim, inner-cylinder origin and radius, and all result-face origins. Incremental edits match fresh STEP output. Regressions also cover multiple bores, rational rigid placement, touching/crossing expanded envelopes, and zero or negative thickness. Side openings in tapered envelopes, multiple openings including a side, bore-wall openings, adjacency filters and ambiguous nearest references remain typed refusals. The independent field compiler is unchanged and still refuses side-open bored shells.

Migration: use the existing `Operation::Shell` / recipe `shell` operation with a resolved planar `FaceRef` and positive thickness. Application edit/save/reopen integration remains a separate SolidDraft gate.
