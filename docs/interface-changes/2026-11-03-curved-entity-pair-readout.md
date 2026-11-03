## 2026-11-03: curved face and edge pair readout

Pieces: 2 (shapes), 6 (measurement), SolidDraft Catallax adapter

Change: Added `catallax_measure::entity_pair` with `PairEntity` (a face or an edge) and `PairRequest` (`MinimumDistance`, `CentreDistance`, `AngleDegrees`). It answers faces on planes, cylinders, cones, spheres and ring tori, and finite straight or circular edges, against each other. Minimum distance is measured on the closed trimmed entities from exact stationary pairs of their faces, edges and vertices; centre distance joins axes, centres, straight-edge carriers and supporting planes; the angle is the unoriented acute angle between axes and normals. New refusal codes are `measure.distance-uncovered` (a pair of strata with no closed form that cannot be ruled out) and `measure.centre-domain` (a line or plane that meets a plane). Spheres have no direction and refuse `measure.angle-domain`. Added `catallax_shapes::distance::Finite::contains_carrier_point`, the exact trim test of a recognized finite face for a point already on its carrier. `edge_pair`, `PreparedFace` and their refusals are unchanged.

Why: A real-use pass of SolidDraft picked two hole walls in the Measure panel and got `measure.nonplanar-face`. No certified readout existed between curved faces or between faces and edges.

Migration: Existing callers need no changes. A front end that wants curved pairs calls `entity_pair` with `Accuracy::default()` and treats marks as an incomplete proof, as it does for `edge_pair`. Parallel planar faces return the same bounds as `PreparedFace::parallel_face_to_face`, converted to an `Enclosure` with an accuracy mark when the search did not reach the requested width.
