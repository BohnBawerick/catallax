## 2026-10-25: sharp selected-face thickening

Pieces: features, field, recipe names and comparison

Change: version 4 thickening accepts proved convex planar face sets and complete circular planar/cylindrical/conical face sets. Front images meet through exact cylindrical, toroidal and spherical projection sectors. Straight open rims and planar sector-end fans close the selected sheet without rounding its boundary. Back images regularise by union. Front/back unions remove their shared source sheets.

Generated far patches use `offset.of`; interior edge and vertex sectors use `round.edge` and `round_corner.vertex`. Boundary end fans retain both incident `rim.edge` origins. The version 4 name checker accepts the new round references. The field builds its own finite images and clipped sectors, without reading B-rep topology.

Why: separate smooth-face sweeps leave a gap at a convex sharp edge, while extending their supports produces a mitre rather than the specified thickness. Tangent round/plane contacts additionally need closed-patch incidence proofs. Exact second-order normal graphs and polynomial parabolic germs now certify surviving tangent patches without admitting isolated buried donors.

Change: generated rotational sectors retain the boundary predicate of their clipped node before union with neighboring images. Exact strict point locations of historical boundary nodes can reject a buried name before contact analysis. Alternative ownership paths require all branches to be excluded; this negative test never admits a face. The symmetric-cylinder regression preserves the exposed round and rejects its buried meridians without changing comparison sampling.

Migration: version 3 keeps its prior sharp-selection refusal. Single-face paths, outward cylinder offsets, tolerances and deterministic work allowances are unchanged. New images reduce proved rational expressions before constructing their arrangements, avoiding needless algebraic elimination in two-sided unions. Nonconvex source skins, branched boundary fans, conical folds, tips and non-ring carriers retain typed coverage refusals.

Prior art: [OpenCASCADE arc joins](https://dev.opencascade.org/doc/refman/html/class_b_rep_offset_a_p_i___make_offset_shape.html) distinguish round joins from intersection mitres. [Normal cones](https://macaulay2.com/doc/Macaulay2-1.24.11/share/doc/Macaulay2/Polyhedra/html/_normal__Cone_lp__Polyhedron_cm__Polyhedron_rp.html) determine the selected projection sectors. [Libfive's offset discussion](https://github.com/libfive/libfive/discussions/545) explains why a generic Boolean residual cannot replace a true distance construction.

