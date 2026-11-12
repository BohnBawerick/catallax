## 2026-11-12: the field builds imported parts whose spline faces are flat bilinear patches

Pieces: field, compare

Change: `Plan::compile_with_files` builds an imported STEP part with spline faces when the shared importer admits the solid and every spline face is a bilinear equal-weight patch whose four controls lie in one plane. The field holds that face as the plane it lies in, oriented by the patch normal and the face sense, and the spline edges that are exactly straight (every control on the line of its ends) as straight sides. The file's face names are kept. A new `catallax_field::flat_support` returns that plane for a B-rep face, so `catallax_compare` can compare the flat spline faces of a mesh with the field's planes.

Why: M21 asks for imported parts with spline faces to build in the field, with the two kernels agreeing on request. The shared importer admits a spline face today only as a complete embedded rational patch with straight boundary rows, so the field takes those it can hold exactly.

Details:

- The check is exact: a bilinear patch is `a + u b + v c + u v d`, and its controls are coplanar if and only if `b x c . d` is zero. A flat patch lies in the plane through its controls, whatever its parameter ranges. The orientation is the sign of `du x dv` at the centre of the unit square; the importer's embedding proof fixes that sign over the whole patch.
- `field.import-spline` remains the refusal for spline parts the field does not hold, now with the reason: a spline solid the shared importer does not certify (the original `M10` reason, with the importer's own code), a degree above one, unequal weights, controls that are not rational, and a bilinear patch that is not flat ("a curved spline face"). The check for these runs after the shared build instead of on the source faces; files the importer refused keep their refusal code and `M10` wording.
- Parts with no flat-spline face are unchanged. The planar route is the existing one; only the plane of a face and the straightness test of an edge widen.
- Curved bilinear patches, patches of higher degree and rational patches are not built; they need an exact membership and trimming route for a curved support, which waits for the field side of M10.

No existing answer changes: every part this admits was refused before. No tolerance, saved replay answer or work allowance changes. Raw evidence retained privately.

Migration: Callers that treated `field.import-spline` as final for every spline face now receive a built field for flat patches. Exhaustive uses of the refusal detail should read it as free text.
