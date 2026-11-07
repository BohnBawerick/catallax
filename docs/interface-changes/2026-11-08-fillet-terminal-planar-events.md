## 2026-11-08: straight fillets that end at a convex boundary event

Pieces: features.

Change: a single straight plane/plane roll-over run whose foot reaches a whole straight boundary edge of its face before the requested radius now builds. The event radius is the smallest positive fraction of the request at which a contact rail meets a straight edge of its face parallel to the selected edge. At that radius the run must pass the face-limit proof from [the face-limit regions](2026-11-07-fillet-face-limit.md), and every rail that ends on a neighbour must lie along one whole straight edge shared with it. That edge must be convex and not tangent. The B-rep then builds the exact region at the event radius. The retained definition keeps the authored request in `requested_radius`, the field added by [the rim consumption change](2026-11-07-fillet-rim-consumption.md). The consumed face is reported with the blend origin.

Why: section 4.5.1 stops a contact sheet when no continuation keeps the centre/foot limit, and keeps the union of the regions already swept. At a convex, non-tangent boundary edge the old foot normal is not normal to the next face, and a ball inside the material cannot rest on a convex edge. So the region stops at the event radius. This is the planar counterpart of the closed-rim terminal event, not a smaller-radius retry: the declared radius stays in the definition.

Refusals: an event at a concave edge keeps its refusal, because the ball could roll onto the next face. So does an event where the boundary edge covers only part of the rail, a neighbour that fails the hull proof at the event radius, or any run outside the single straight plane/plane stratum.

Evidence: `past-face-limit` (radius 10.001) and `swallowed-face` (radius 12) build the same six-face solid as `at-face-limit`, the region `W_10`, with volume `2000 + 500*pi`. OCCT refuses all three. Material witnesses, the recorded request and the consumed side face are checked natively. A concave-edge control keeps its refusal. Full-corpus results belong to the clean published head's report.
