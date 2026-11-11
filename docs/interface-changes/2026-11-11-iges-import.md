## 2026-11-11: IGES manifold solids read into the shared importer

Pieces: 8, STEP and IGES consumers

Change: `catallax_import` gains `read_iges`, `read_iges_with_decider` and
`import_iges_bodies`. An IGES 5.3 manifold solid B-rep object (type 186) is read
into the same `StepFile` a STEP file gives, then built by the unchanged
`build_file`. Face keys are directory-entry sequence numbers such as `D17`;
edges and vertices are named by list and index, such as `D3:5`.
`import_iges_bodies` checks the saved hash and self-intersection and builds
without repairs, as `import_step_bodies` does for STEP. The supported entities,
conventions and refusals are listed in
[piece 8](../pieces/08-step-import-repair.md#iges-source-recognition).

Why: PLAN M24 asks for IGES import. Sharing the recognized source model keeps
one builder, one set of exact incidence checks and one set of refusal codes for
both formats.

Coverage, raw evidence retained privately: files from Catallax's IGES writer
read back as the same body for a block, a cylinder, a 45-degree frustum, a
cylinder turned a quarter about an axis and a bored plate, with overlapping
certified volumes. Files written by OpenCASCADE 8.0 in B-rep mode read exactly
for a box (volume 30) and a bored plate. Its cylinder, cone, sphere and torus do
not: their seam vertices and arc matrices carry binary sine and cosine noise of
order 1e-16, so they refuse `not_solid` or `Unsupported` (a matrix that is not
exactly orthonormal) rather than being moved. Hand-made files cover an inch
cube (exactly 16.387064 cm³), a sphere closed through its seam and degenerate
pole uses, voids, an unknown unit, a surface model, a scaling transform and an
arc of irrational radius.

Migration: none for existing callers. No STEP behaviour, recipe format, saved
answer, tolerance or work allowance changes. Exhaustive matches on
`ImportRefusal` are unchanged.

Known gaps: there is no `import_iges` recipe operation yet, so the application,
field and comparison routes do not see IGES files. Near-exact files from other
writers need the curved repairs with recorded local tolerances that M24 lists
separately. Trimmed-surface models (143, 144), surfaces of revolution, ruled and
offset surfaces, composite and conic curves remain typed refusals.
