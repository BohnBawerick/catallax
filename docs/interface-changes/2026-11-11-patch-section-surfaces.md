## 2026-11-11: Rational patch sections take a cylinder or sphere

Pieces: shapes, intersect; consumers that match on `Curve` or `CurveForm`.

Change: `RationalPatch::plane_sections(&Plane)` becomes `RationalPatch::sections(&SectionSurface)`. `SectionSurface` is `Plane`, `Cylinder` or `Sphere`; `SectionSurface::of` admits a `Surface`. The plane-only names become surface-neutral: `PatchPlaneBranch` is `PatchSectionBranch`, `PlaneSectionCandidate` is `SectionCandidate`, `PlaneBranchApproximation` is `SectionBranchApproximation`, `PlaneSectionSearch` is `SectionSearch`, `Curve::PatchPlane` and `CurveKind::PatchPlane` are `PatchSection`, and `CurveForm::PatchPlane` is `CurveForm::PatchSection`. `PatchSectionBranch::plane` and `RationalPatchSection::plane` become `surface`, and `follows_plane_first` becomes `follows_surface_first`. `catallax_intersect::supports` now admits a rational patch with a cylinder or sphere, and the pair planner returns their exact branches.

New refusal codes: `intersect.patch-cylinder-tangent`, `intersect.patch-sphere-tangent`, `intersect.patch-cylinder-selector`, `intersect.patch-sphere-selector`, and `intersect.patch-{cylinder,sphere}-candidate-{domain,endpoint,cap}`. Plane codes are unchanged.

Why: M10's second slice cuts a patch by a cylinder or sphere with the plane slice's search and tube. A cylinder branch named `PatchPlane` would misstate its definition.

Migration: replace the old names. Pass `&SectionSurface::Plane(plane)` where a plane was passed. The plane slice's answers, refusal codes and recorded work units are unchanged. See [M10 cylinder and sphere sections](../M10-QUADRIC-SECTIONS.md).
