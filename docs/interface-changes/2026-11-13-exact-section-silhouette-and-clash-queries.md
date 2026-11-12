## 2026-11-13: sections, silhouettes and exact body relations

Piece: application (`catallax_app::queries`). Additive API; no caller migration.

`section_curves` returns an owned model, exact trimmed `SectionEdge` branches and
isolated contact points for a rational plane. It uses the existing intersection
graph against a cutting block whose other five sides lie strictly beyond a
certified enclosing box. It reads no display positions. Endpoints and curve
components keep their original exact definitions. The supplied `Accuracy` controls
the enclosing box, not a section fit or geometry tolerance. A section containing a
planar face refuses with `app.section-coplanar` until planar region extraction is
covered. Other intersection and measurement refusals propagate unchanged.

`silhouette` returns orthographic outline generators in model coordinates, for a
rational direction toward the viewer. A convex planar solid returns exact finite
edges between a front-facing face and a back-facing or edge-on face. The query
proves convexity against every outward face halfspace and checks all boundaries
are finite lines with two face uses. Whole spheres return an exact great circle.
Multiple shells, nonconvex bodies and other curved surfaces refuse with
`app.silhouette-occlusion`, `app.silhouette-nonconvex` or
`app.silhouette-surface`. Consumers project the returned definitions; these are
not display polylines or a list of potentially hidden edges.

`overlap_volume` intersects two retained bodies and measures the result with the
existing volume integrator. It returns an enclosure in cubic millimetres, with
measurement and construction marks intact. Empty overlap is successful zero.
`contact_faces` intersects, imprints and classifies the boundaries, returning an
owned model with clipped patches from both operands. `OnSame` and `OnOpposite`
distinguish equally oriented shared boundaries from opposing contact boundaries.
Edge and point touches have no contact area. Patch origins remain available in
the returned model. Lists are not assumed to pair one-to-one. Volume and contact
are separate requests so a refusal in one does not discard an answer to the other.

Free functions accept any valid retained model and body handles, including an
incremental checkpoint. `ApplicationBody` has matching methods; its two-body
methods take the other handle from the same build's model. Intermediate bodies
are available through `Built::bodies`. Foreign-model handles are programming
errors, as in the existing topology queries.

Every new query uses the caller's active arithmetic allowance or the unchanged
20,000,000-unit application default. `QueryRefusal` retains the original `Refusal`,
including any undecided numeric bound, plus each operand's face, curved-face,
edge and vertex counts. No partial result replaces a refusal. The source model,
recipe, display and picking names remain untouched.

Already present on main: regularised empty Boolean bodies and their application
display. `empty_boolean_is_a_valid_empty_display` in `tests/application.rs`
exercises self-subtraction through `build_soliddraft`; `catallax_boolean::boolean`
also returns `BodyKind::Empty` for apart intersections. This change reuses those
answers and does not alter Boolean, tolerance or replay rules.

Remaining work: smooth contour trimming for cylinders, cones, tori and spline
faces; visibility splitting for nonconvex bodies and several shells; and section
boundaries on coplanar faces. These cases retain typed refusals. Exact sections,
contact and overlap on freeform bodies inherit the Boolean engine's current
coverage and work limits.

Focused release-part coverage: P2's mid-height section has two exact circles with
squared radii 36 and 144 mm². Its overlap with the retained unchamfered ring is
between 3998 and 3999 mm³, and each operand has four clipped, equally oriented
shared patches. Its silhouette refuses with six faces, four curved. P4's section
at z = 6 mm has eight line branches and four circular branches, preserving both
the outer boundary and the cavity. Unit coverage also checks an oblique cube
section, a tangent sphere's isolated point, axis and oblique convex silhouettes,
partial-face contact, edge-only contact, empty results, retained construction
marks and exhausted decision allowances.
