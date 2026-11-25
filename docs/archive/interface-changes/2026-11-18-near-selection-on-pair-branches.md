Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

## 2026-11-18: a fillet's near point picks between straight or circular branches of one surface pair

Pieces: 7, features; SolidDraft adapter
Change: when a fillet edge reference with a `near` point resolves to several edges, an edge whose curve is stored as a branch of a surface pair is now measured exactly where `catallax_intersect::curve_form` shows the branch to be a line or a circle and both edge vertices lie on it: as the segment between its vertices, or as the arc the branch's own tangent at the start vertex runs along. The edge at strictly the smallest distance is chosen. Equal smallest distances refuse `recipe.ambiguous-reference`, as for stored lines and circles.
Why: a cylinder meets a plane along its axis in two lines, and the kernel keeps each as a branch of the pair. Two such edges between the same faces refused `fillet.unsupported-topology` ("nearest-edge selection among curved edges needs an exact trimmed-edge distance query") whatever the near point, so SolidDraft could not round either.
Migration: none. Callers that send `near` get the edge it is nearest to where they got that refusal. The refusal remains for branches with no line or circle form: ellipse, parabola and hyperbola arcs, splines and patch sections.
