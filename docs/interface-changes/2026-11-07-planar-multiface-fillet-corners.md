## 2026-11-07: exact spherical polygons at fully selected planar vertices

Pieces: features, recipe CLI

The B-rep now closes a fully selected planar vertex of more than three faces
when all incoming strips have the same radius and material side, and their
axes meet at one ball centre. It previously refused these junctions with
`fillet.unsupported-corner`.

The constructor proves the ball contacts every input plane at the requested
radius, the contact normals lie in one open hemisphere, and each selected
pair bounds an exposed edge of their positive cone. Every support has two
incident strips. These predicates identify the complete spherical convex hull
required by RECIPE section 4.5.1. Its shorter great-circle arcs are shared
edges with the strips, not separately fitted trim curves. Non-incident input
supports must be clear of the sphere. The centre must lie within each input
edge's axial interval, and the retained rails must remain positively ordered
between their two corner ends. Consumed strips still require continuation.

The sphere is exact over the whole polygon. There is no fitted carrier or
approximation allowance. The patch keeps `corner.vertex`, including all of
that vertex's input face origins. Strip names remain `blend.edge`. Existing
three-contact construction and the public `CornerDefinition` and `SphereSector`
types are unchanged. There is no recipe, schema or public API change.

The four-valence pyramid constructs at the ordinary 1M allowance, then
exhausts that allowance during STEP export. An explicit 20M replay exports
one raw-valid solid at 14,012,254 work units, including 257,925 for the final
fillet. Work-limit optimization is not part of this change. The native body
has 18 faces, including one four-sided spherical apex patch; STEP has 20
faces after chart splitting. Its composed file bound is
5.569808867159556034910916092483513545794e-9 mm.

Coverage includes the original pyramid, a translated radius edit and a
shallower pyramid whose contact normals are not perpendicular. Exact centre,
radius, shared great-circle edges, material witnesses and the independent
inner-parallel volume formula are checked. Raw STEP import uses no healing.

Partial selections at four-valence vertices, non-co-spherical contacts,
mixed convex/concave fans and unequal earlier mitres remain separate
construction gaps. The independent field evaluator is unchanged and does not
acquire this multiface construction. General valence-five coverage and
SolidDraft interaction are not claimed.
