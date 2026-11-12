# M16 distance, loads and material regions

The opt-in solver adapter adds `boundary_distance` and `boundary_distance_with`,
with distance enclosures, admitted witnesses and explicit direction/owner ties.
The P1 exterior probe now completes within its unchanged work allowance.

`face_patch` and `face_patch_with` integrate supports, traction, total force and
pressure on named planar and supported curved faces. Unsupported curved trims
refuse `solver.patch-coverage`, including a partly surviving spherical band.
Both `_with` methods reject a service from another definition instance with
`solver.definition-mismatch` before assigning a revision to its result.

Optional recipe `materials` entries select body-valued steps as material regions.
`Definition::from_recipe` compiles them through `Materials::Recipe`. Box regions
must partition the complete certified bounds; gaps and overlaps refuse
`solver.material-partition`. Curved regions retain conservative per-region
moment enclosures and explicit unresolved partition cells.

Comparison sampling retains the original typed refusal on failure and clears
its incomplete sampled distance bound. The 750,000-unit P3 regression remains
a resource refusal. Its allowance and assertions are unchanged.

No solver, tolerance, work allowance or saved replay answer changes. See
[the query contract](../SOLVER-QUERIES.md) for coverage and numeric evidence.
