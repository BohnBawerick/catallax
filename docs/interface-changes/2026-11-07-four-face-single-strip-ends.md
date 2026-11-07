## 2026-11-07: a single planar strip ends at a four-face vertex

Pieces: features, recipe CLI

One selected straight edge of a convex polyhedron can now end at a four-face
vertex. It previously refused with `fillet.unsupported-end`. The selected
strip remains its exact cylinder. Two conic arcs on the unselected end walls
close it, meeting at the unique hit of those walls' shared input edge in the
strip's shorter normal sector. The other input edges shorten to these exact
contacts; connection alone never starts another blend run.

The input skin must have planar faces and straight edges, with every vertex
in every outward halfspace. This proves the convex material sector. Both end
walls must be transverse to the cylinder axis and face the remaining run.
Their conics project one-to-one onto the cylinder's normal circle. Excluding
the other rail selects the correct sector interval for each arc. Every
plane crossing of each complete conic is enumerated. No event may lie in the
retained arc's interior, and an interior point fixes its material side.
Thus the proof covers each complete finite arc, not sampled endpoints alone.
The retained contact rails must remain positively ordered between their ends.

The cylinder and its cuts are exact. No fitted model face or approximation
allowance is introduced. Each shared boundary is one edge, strip naming
stays `blend.edge`, and the surviving input walls keep their origins. A
single strip end creates no finite-area corner face. There is no public API,
recipe schema or overflow-policy change.

`vertex-4-one` exports a raw-valid six-face solid at an explicit 20M allowance.
It uses 4,167,270 work units, including 2,658,585 for the final fillet. At the
ordinary 1M scoreboard allowance it still exhausts construction work. The
snapshot therefore moves from `fillet.unsupported-end` to
`numbers.operation-budget`, not to a checked export at 1M. Work optimization
is deliberately deferred.

Seven material witnesses and an independent clipped-section volume integral
check the actual CLI export. Its reported file bound is
4.564144149380011770722525218291367494769e-9 mm. Original, translated,
radius-edited and shallower-pyramid controls check closed shared topology and
material membership.

Two selected strips at a four-face vertex, curved or nonconvex end fans and
unequal earlier mitres remain separate construction work. The independent
field evaluator is unchanged.
