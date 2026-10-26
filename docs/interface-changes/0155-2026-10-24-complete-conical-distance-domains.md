## 2026-10-24: complete conical distance domains

Pieces: features, field and comparison tests

Change: `cone_skin::ConeSkin` retains convex conical cells with exact full circular or elliptic caps and optional coaxial cylindrical neighbours. Inward offsets rebuild shared native conics and retain the positive nappe. `cone_sweep::ConeSweep` constructs finite conical band normal images with analytic sloping rims. Complete-section recognition checks the source topology before accepting these definitions.

Why: translating a cone support alone neither proves a finite offset body nor supplies the rims of a normal sweep. The field compiles its own constraint inventory and excludes the other algebraic nappe without registering an invented physical face. Additional rational cone candidates and finite bounds are query-only.

Migration: no recipe, exchange format, tolerance or work allowance changes. Unproved outward joins, concave skins, intersecting trims and conical folds retain typed refusals. The property-integration coverage boundary is unchanged.

Prior art: [OpenCASCADE offset joins](https://dev.opencascade.org/doc/refman/html/class_b_rep_offset_a_p_i___make_offset_shape.html) distinguish extended supports from arc joins. [Truck's normal offsets](https://github.com/ricosjp/truck/blob/master/truck-geometry/src/decorators/offset/surface.rs) retain a parameter domain. These distinctions motivate separate finite-domain and join proofs here; no source code is copied.

Change: conical distance comparisons add a separate planar polar candidate chart. Each ray splits at the field's own quadratic constraints, including circular and elliptic trim crossings. Every point still requires independent finite-trim, origin and direction certification. This resolves the narrow open conical shell rim without increasing work allowances or adding part-construction work.

Change: outward perpendicular-capped frustum offsets directly assemble conical, planar and toroidal sectors. The field evaluates a separate nearest-face and circular-edge decomposition through `ConeSkin::outward_pieces`. Circular cap thickening retains only the selected disk's normal image.

Why: an outward support mitre is not a distance offset. Direct assembly avoids repeated B-rep Boolean unions while preserving shared exact meeting circles.

Migration: outward cone tips, cylindrical neighbours, oblique noncircular edge carriers and horn or spindle tubes retain typed refusals. Cone-tip display refinement remains unproved at the existing 1/200 mm chord tolerance; no mesh, geometric tolerance, work allowance or property-integration path is changed to bypass it.

Change: conical torus candidate charts partition meridians at the field's independent plane and coaxial quadric crossings. Native and distance-feature visits share this query-only chart and sample open parameter cells rather than their G1 contact wires. Every emitted point still needs the unchanged finite-trim, origin and direction certificates. No finding is filtered or comparison predicate weakened.

