## 2026-11-11: nearest named face (F8)

Pieces: field

Change: `catallax_field::distance::Prepared` gains `nearest_face(point, options, dec)` and `prepare_nearest(dec)`. The module adds `NearestOptions`, `Ownership`, `Owner`, `Ambiguity` and `Candidate`. `nearest_face` returns `Ownership::Proved(Box<Owner>)` when one named face's certified upper distance bound is strictly below every other name's lower bound. Otherwise it returns `Ownership::Refused(Ambiguity)`, whose `refusal()` has the new stable code `field.ownership-ambiguous`. Invalid options refuse with the new code `field.ownership-options`; a result without a surviving named face refuses with `field.ownership-empty`.

Every answer is proved under the [M14 decision](../ARCHITECTURE.md#unproved-hardware-preview-evaluator-decision-of-2026-11-10). The query reads exact definitions and certified bounds only; no hardware preview value, clock or worker count selects an answer or a proof branch. An `Owner` carries the distance enclosure, an exact witness, the enclosed unit direction toward that witness, the runner-up name and a margin. The margin is a proven lower bound on the distance to the nearest ownership border. An `Ambiguity` carries the contenders' bounds and their spread. [FIELD.md](../FIELD.md#nearest-named-face) gives the method and its limits.

Why: PLAN M14 F8 asks for the nearest named face from any point, for loads and supports on named faces, picking and the proof light. `face_query` answers only whether a point is on one named face.

Migration: none. Existing `point_to_face`, `point_to_boundary`, `face_to_face`, `face_query` and display answers, refusals and work counts are unchanged. Callers that share one `Prepared` between workers call `prepare_nearest` once beforehand, so per-query logical work does not depend on the partition. It fills the same planar-region cache that `point_to_face` already fills on demand.
