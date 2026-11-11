## 2026-11-11: Atlas sections, branch identity and branch meetings

Pieces: shapes, intersect; Boolean callers of curve identity and curve-surface hits.

Change: `RationalAtlas::sections(&SectionSurface)` returns an `AtlasSection` whose `PatchSectionBranch` components run over every chart of a joined face, joined across chart seams and periodic rows. `PatchSectionBranch::patch()` becomes `charts()`, a slice holding one patch for a single-patch section. `PatchSectionBranch::same_component` decides whether two branches have the same point set when the proof is cheap and exact. `PatchSectionBranch::surface_hits(&SectionSurface)` returns `SectionHits::Contained` or the exact transverse meetings in traversal order.

`catallax_intersect::supports` admits a `RationalAtlas` with a plane, cylinder or sphere, and the pair planner returns its branches. `same_curve` now proves patch branches equal or apart when `same_component` decides, and proves a branch apart from a line or B-spline when one of its cell-end witnesses is off that curve. `intersect_curve_surface` answers a patch branch against a plane, cylinder or sphere.

New refusal codes: `intersect.atlas-section-row` (a cut contains a whole chart boundary row) and `intersect.patch-section-hit` (a tangent or partly contained meeting, or a meeting that cannot be isolated within the unchanged allowance).

Why: M10's second slice needs loft faces, which are atlases, and the Boolean of a box with a lofted body needs branch identity and third-surface meetings (items 3, 4 and 5 of the M25 scope study's section 4.6).

Migration: replace `branch.patch()` with `branch.charts()[0]` for single-patch branches. Code that matched `intersect.unsupported-pair` for an atlas with a plane, cylinder or sphere now receives branches or the section's own typed refusals. Plane, cylinder and sphere sections of single patches keep their answers and recorded work units. See [M10 atlas sections and meetings](../M10-ATLAS-SECTIONS.md).
