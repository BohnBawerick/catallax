# Isolated plane contacts in complete spline sections

Pieces: shapes, intersection and Boolean consumption.

`RationalPatch::sections` and `RationalAtlas::sections` can now retain isolated
plane contacts alongside regular branches. `AtlasSection::points()` exposes the
contacts. The intersection adapter publishes them as `Component::Point`, so the
ordinary Boolean graph can consume them. `SectionSearch::contact_cells` counts
cells with corner contacts separately from empty cells and regular arcs.

The new check uses exact Bernstein coefficient signs. A cell with coefficients
of one weak sign and at least one nonzero coefficient has no interior zero.
Every edge must also have a nonzero coefficient. Its only possible zeros are
then corners, which are evaluated on the retained patch and kept explicitly.
An identically zero edge fails this check. Opposed coefficient signs also fail
it. The search still accounts for every cell and every regular branch.

At an exact stationary midpoint, a positive Hessian determinant permits further
subdivision instead of the previous immediate tangent refusal. This local
classification does not discard the cell or certify its whole zero set. The
complete Bernstein checks must finish before the section can return. Duplicate
contacts at chart joins are merged by exact point equality. Points covered by a
regular branch are excluded from the isolated-contact list.

This is a cheap exact path. It does not create a numerical result or change a
tolerance. The existing work, cell and depth allowances still apply. Contacts
that these checks cannot resolve keep their refusal. Curve contacts and higher
degeneracies need their own finite classification. Cylinder and sphere contact
handling is unchanged.

The regression uses the recorded bowl, border and hill fixtures in their five
exact placements. It also checks the tiny hill loop, reversed intersection
operands, worker counts 1/2/4, and both subtraction orders at the M10 roof peak.
The Boolean checks retain topology and strict interior/exterior material
witnesses. STEP export is a separate remaining gate for the rebuilt atlas body.

The approach follows OCCT's separation of isolated points from intersection
curves in [`IntPatch_ImpPrmIntersection`](https://github.com/Open-Cascade-SAS/OCCT/blob/master/src/ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_ImpPrmIntersection.cxx).
Catallax keeps its complete cell search and uses exact signs for this cheap
case. The [P04 port map](../phase1/P04/PORT-CONSUMERS.md) and the
[known SMLib misses](../M9-KNOWN-SMLIB-MISSES.md) explain why an empty candidate
trace cannot establish an empty section, especially beside a tiny loop.
