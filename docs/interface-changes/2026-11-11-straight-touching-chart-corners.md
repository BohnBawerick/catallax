## 2026-11-11: straight corners and one-edge rings at touching spherical trims

Pieces: 7, STEP consumers

Change: A touching spherical chart now accepts a corner whose two branches leave the shared vertex on opposite sides, including a one-edge ring whose corner pairs with itself. Such a corner gets no pcurve repair. Before this change the writer refused it and fell back to the unseamed legacy encoding. Hole nesting also accepts a hole whose first point is a vertex on its outer boundary.

Why: A ring tangent to another trim, such as a bore's circle touching a latitude cut on a ball, has one vertex. At that vertex both wires pass straight through, so the pinched-corner joint of [certified joint pcurves at touching spherical trims](2026-10-30-certified-touching-sphere-pcurves.md) does not apply. The legacy encoding of these faces fails the direct no-healing read.

The certificate is exact. Every control point of each corner span after the vertex lies strictly on that branch's side of a line through the vertex. Positive Bernstein weights keep the span on that side, so the two corner spans meet only at the vertex. Every other span of either curve clears each corner span: it stays at or beyond that span's far end along the line, strictly behind the vertex, or strictly outside that span's control box. In a ring the first and last spans are the corner spans and only the spans between them are tested. Nothing moves, so the UV error and the STEP deviation bound are unchanged by the corner. Pinched corners still try the joint repair first; a corner that passes neither certificate keeps the legacy fallback.

A hole's containment query used the hole's first point. For a one-edge ring that point is the touching vertex, which lies on the outer boundary, and the query refused with a budget error. The writer now tries the hole's span starts in walk order and uses the first that the exact winding test does not find on an outer boundary. A hole whose first point lies on no outer boundary is queried exactly as before.

Migration: No recipe, question, native topology or property interface changes. Native answers keep their exactness, marks, properties and topology. Answers of this class get new STEP bytes: charted faces with seams and pcurves in place of the raw-invalid legacy faces. Only the STEP-report fields `added_edges`, `added_vertices`, `approximated_curves` and `max_deviation` can differ. The public regression `test_one_edge_ring_touching_a_latitude_trim_is_raw_valid` checks raw and standard validity and an independent volume.
