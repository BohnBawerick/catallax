Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-20: edge-gap checker in piece 3

Pieces: 3
Change: `catallax_topology::numerical::gap` is new. `gap::measure` takes the two copies of one edge use (`EdgeCopies`: a 3D curve, edge parameters, a `GapSurface` and a `FaceCopyRef`), `Hints` and `Limits`, and returns a `GapMeasure`: the largest same-parameter distance found, its measurement error (rounding plus what the samples left unresolved), a domain exit if the face copy leaves a spline surface, the first sample plan and every later doubling of a piece's points with its reason. `gap::plan` takes the same copies, hints and limits and returns the first plan alone. New refusal codes are `numerical.gap-unsupported` and `numerical.gap-work-cap`.
Why: roadmap item 3.0A replaces the nine-point and knot-span sampling of numerical admission with a checker that finds the largest gap. See [the checker note](../milestones/3.0A-edge-gap-checker.md).
Migration: none. Nothing calls the module yet, and no existing type, check name or refusal changes.
