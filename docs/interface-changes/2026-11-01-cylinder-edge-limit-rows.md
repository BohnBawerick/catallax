## 2026-11-01: an edge limit on a cylinder face seeds rows

Pieces: mesh

Change: when `max_edge_length` is set and a cylinder or elliptic cylinder face is too long
for one grid cell to meet it, `mesh_body` seeds interior rows at that face's rail columns
before certification. A rail sampled finer than the face needs keeps only the columns the
face's chord, angle and length ask for. Chord, rounding, normal, aspect, length and closure
proofs, every cap and work allowance, the geometry and the public types are unchanged.
`ruled.rs` still declines a strip that exceeds the limit, and `finish` still measures every
returned edge. Preview takes no edge limit and is unchanged.

Why: `edges.rs` splits a long ruling at the limit, so `open_strip` declines the strip and the
general triangulator takes the face. Its grid sized cells from chord and angle alone, which
gives a cylinder fine columns and no rows. `finish` then split long edges one midpoint at a
time, and Delaunay refinement in a chart of turns by millimetres spread toward square chart
cells. A block with one rounded edge reached `max_vertices` on that face at the drawing chord
and 2 degrees with an edge limit of 8 or less.

Migration: callers that matched `mesh.precision-or-budget` for such a face now receive a mesh.
A closed or partial cylinder that failed with "periodic or trimmed triangulation did not
close consistently" under a short limit can now answer too. Output that already answered with
an edge limit shorter than the face changes, usually to fewer vertices and sometimes to more:
P3 at a 1.5 mm limit grew from 76,809 to 95,130 vertices, because four cylinder faces with
finely sampled rails now take rows. Calls with no edge limit, and bodies whose faces already
meet the limit in one cell or have no cylinder face, keep their exact output. In a changed
body the cone, sphere, torus and plane faces keep their triangles, though compacted indices
can differ. Exact properties and STEP bytes are unchanged.
