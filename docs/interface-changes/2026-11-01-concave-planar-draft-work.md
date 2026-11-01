## 2026-11-01: concave planar draft at ordinary angles

Pieces: features, field

Change: `draft` now answers for a concave line-bounded extrusion at ordinary angles within
the 20M application allowance, in both evaluators. A hole-free L bracket drafted at 3 and 5
degrees, and a half-scale copy in an oblique frame, answer with eight planar faces and
exact raw-valid STEP. `catallax_features::draft::arrangement` adds `Foot` and
`Polyhedron::foot`, which read a point on a face's carrier against its trim from enclosures.
Each polyhedral edge's line now takes its direction from its two carriers, and edges on the
same pair of carriers record a `Coincident` construction fact. The field's distance to a
polyhedral cell is now the minimum of each candidate's bounds. It stays a rigorous enclosure
of true distance, and the query's exact location is unchanged.

Why: each failure was one exact zero sent to a general proof, whose reservation exhausted the
request. Arrangement vertices are now direct plane triples instead of chains of interpolations,
and the zeros that hold by construction are read from it. A query point's own ties and trim
contacts need no proof, because a distance enclosure does not depend on which candidate is
nearest. No allowance, work charge or guard changed.

Migration: callers that matched `numbers.operation-budget` for this family now receive an
answer. A preserved bore or boss on a concave base keeps its refusal, as do touching or
crossing tools and curved concave sources. Recipes, budgets, tolerances and every earlier
answer are unchanged.
