## 2026-10-28: another sphere frame when a pole lies on a trim

Pieces: 7, STEP consumers

Change: A trimmed spherical face whose parameter pole lies on one of its trims
now joins the partial-chart path of
[the regular periodic trims](2026-10-26-partial-periodic-step-trims.md) and
[the interior pole caps](2026-10-28-step-sphere-pole-caps.md). Neither a chart
nor a cap of the stored frame admits such a pole. The frame only parameterizes
the sphere, so every trimmed face on that sphere is written on one new
`SPHERICAL_SURFACE` with the same centre and radius and another exact frame.

Candidates come in a fixed order: the stored axis `a` and reference `r` swapped,
then `a × r` as the axis, then the diagonals `a ± r`, `a ± a × r` and
`r ± a × r`. Each pair is perpendicular because the stored pair is. The first
candidate whose two poles are certified outside every face on the sphere, or
strictly inside one where a cap then applies, is used. Membership is the same
exact query as before; it depends only on the trims, not on the stored frame.
If no candidate is admitted, the whole scene keeps its legacy encoding. The
retry on perpendicular frames applies the same rule to its own poles.

Scenes whose stored frames are already admitted keep their bytes. Kernel
topology and geometry are unchanged. Report fields and manifest formats do not
change.

Why: A ball cut by planes through its axis puts a pole on a trim corner. The
legacy encoding of such a scene fails the pinned raw translator.

Migration: None for callers. Consumers see a different `SPHERICAL_SURFACE` axis
and reference direction for a sphere whose stored pole lies on a trim.
