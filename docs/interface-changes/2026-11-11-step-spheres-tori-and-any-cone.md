## 2026-11-11: STEP import builds spheres, tori and cones of any semi-angle

Pieces: 2, 3, 8, field, compare, mesh

Change: The shared STEP importer builds spherical and toroidal faces and cones at any semi-angle. It no longer refuses them as `unsupported`. `catallax_import::Surface` gains `Cone(Placement, Rational, SemiAngle)` for cones other than exactly 45 degrees in a degree context, and `FaceLoop` gains `vertex: Option<Vertex>` for a `VERTEX_LOOP` bound. The field builds the parts these admit, with the file's face names. It extends [the earlier subset](2026-11-11-field-step-import-subset.md), whose closing paragraph no longer describes the importer. `field::mass::volume` also integrates sphere and torus primitives.

Why: M21 asks for imported parts made of planes, cylinders, cones, spheres and tori in both kernels, with the two kernels agreeing on request.

Details:

- Spheres and tori. A face whose bounds reduce to seams, or to one `VERTEX_LOOP` on the surface, is the complete closed surface, as `catallax_topology::sphere` and `torus` build it. Circle edges on spheres, and torus parallels and meridians, are checked exactly. Other circles on a torus are `unsupported`. A `VERTEX_LOOP` beside other bounds or on an open surface is `unsupported`; one off its surface is `not_solid`.
- Cones. tan(semi-angle) is irrational at every writable angle except 45 degrees, so a decimal rim away from the placement plane never lies exactly on the declared cone. The importer builds the cone through the face's two coaxial rims, or its rim and apex. The apex is a `VERTEX_LOOP` or a seam end on the axis. That cone has a rational slope. The whole-face deviation from the declared cone is certified, and it is accepted only within the existing 0.01 mm automatic cap and a tenth of the face's nearest feature. The apply route records it as a `FaceTolerance` with `approved_tolerance_mm` `0.01`. Beyond either limit the import is `not_solid`, with both bounds in the detail. 45-degree degree-unit cones keep their exact declared build and record nothing.
- Field. A complete outward sphere or torus becomes the field primitive of the same size and placement. A torus axis of irrational length returns `field.import-domain`. A pointed cone with one disk becomes a round cell with an apex. Trimmed spheres and tori return `field.import-domain`; they no longer return `field.import-unsupported`.
- Mesh preview draws a rim-free declared cone, unproved, from a 256-bit enclosure of its slope.

No existing answer changes. Every file these rules now admit was refused before. Canonical derivative text for existing surfaces is unchanged. No tolerance, saved replay answer or work allowance changes. Raw evidence retained privately.

Migration: Exhaustive matches on `catallax_import::Surface` add the `Cone` arm. Struct literals of `FaceLoop` add `vertex: None`. Callers that expected `field.import-unsupported` for spheres or tori now receive a built field, or `field.import-domain` for trimmed ones. Readers of `RepairRecord::face_tolerances` will see entries for rim-built cones in otherwise clean imports.
