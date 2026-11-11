## 2026-11-11: Rational patch sections take a cone or torus

Pieces: shapes and intersect; consumers matching `SectionSurface`.

Change: `SectionSurface` adds `Cone(Cone)` and `Torus(Torus)`.
`SectionSurface::of`, `RationalPatch::sections` and
`catallax_intersect::supports` admit these families. The pair planner returns
exact retained branches with the existing membership, order, tangent and
orientation queries. The torus must be a ring torus. A cone span must have a
proved strict nappe before its polynomial is used.

New refusal codes: `intersect.patch-cone-nappe`,
`intersect.patch-{cone,torus}-{tangent,nonregular,selector}` and
`intersect.patch-{cone,torus}-candidate-{domain,endpoint,cap}`.
These families refuse multiple roots in rational reductions as well as
nonregular bivariate contacts. Existing families keep their decisions and
work accounting.

Migration: handle the two added enum variants in exhaustive matches.
See [M10 cone and torus sections](../M10-CONE-TORUS-SECTIONS.md) for the
complete-search proof, measurements and remaining limits.
