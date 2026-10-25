## 2026-10-24: physical draft references and axial separation

Pieces: features, field

Change: `Skin::resolve_faces` applies face patterns and shared-edge adjacency before comparing exact distances to closed planar trims. Projection inside a face and nearest points on its outer or hole edges use the same exact squared-distance order. Ties return `recipe.ambiguous-reference`; missing faces and swallowed-cap picks return `recipe.missing-reference`. Physical selection precedes alias expansion. Both kernels reuse their extracted skin for selection and construction.

`Skin::is_convex` retains the direct support-intersection construction for convex Boolean inputs and filtered references. This avoids imposing moved-fan restrictions on their non-simple vertices.

`RoundCell::certified_disjoint` proves axial or radial separation from exact height intervals and affine radii. The B-rep adapter no longer assembles intersecting drafted round cells as independent untrimmed lumps. Such cells need a reconstructed curved arrangement and retain `draft.unsupported-arrangement`; this is not an impossible-input classification. Separated axial cells retain their existing construction. Mixed curved field inputs and filtered curved references remain unfinished.

Why: coactive names can label disconnected pieces, and adding draft can make originally disjoint cylindrical components intersect. Neither shared ancestry nor source separation proves the replacement topology.

Migration: no new origins, recipe version, approximation tolerance or work-allowance changes. Face-reference errors use the recipe namespace. This supersedes the planar filtered-reference limitation recorded above, not the remaining curved or non-simple nonconvex limits.

