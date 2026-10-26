## 2026-10-22: B-rep display meshes and OBJ export

Pieces: display mesh, consuming 1 to 4; the judge interface is unchanged
Change: new crate `catallax-mesh` exports `mesh_body`, `MeshOptions`, `Mesh`, `Triangle`, `FaceMesh` and `EdgeMesh`. A mesh has shared binary64 positions, face-tagged triangles, per-face proven chord bounds and per-edge certified polylines. `Mesh::write_obj` writes shared vertices and `face_<id>` groups. Options specify an exact positive millimetre chord tolerance and a deterministic vertex limit. Invalid options return `mesh.invalid-options`; resource, chart-resolution or coordinate-precision failures return `mesh.precision-or-budget`. Numeric trace/deadline errors propagate and input marks remain attached.
Why: PRD section 6.5 requires meshes from the B-rep for display and later background comparison with the field kernel. Display approximations must not replace the exact body or claim the 1 nm agreement guarantee.
Migration: new library consumers may mesh an existing `Model` and `BodyId`, then pick by `Triangle::face`. Arena face IDs are not recipe origin names. No existing caller, STEP byte format, question format or manifest format changes.

