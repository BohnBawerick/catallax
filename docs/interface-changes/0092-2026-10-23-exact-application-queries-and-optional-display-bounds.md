## 2026-10-23: exact application queries and optional display bounds

Pieces: measurement, Boolean point queries, meshing and application
Change: `catallax-measure` adds `face_centroid`, `edge_length`, `edge_midpoint`,
`bounding_box`, `face_bounding_box` and `BoundingBox`. `ApplicationBody` borrows its
existing build for each query. Bounds select analytic extrema with exact trim
membership and use budgeted second-derivative enclosures for general branches.
`face_contains` accepts exact `Point3` through a source-compatible generic argument;
its rational caller and classifier are unchanged. Measurement now depends on the
Boolean crate without introducing a dependency cycle.

Display options add `angular_deviation`, `max_edge_length` and an optional
`max_aspect_ratio` request. Facet normal certificates include analytic normal
variation across each parameter triangle and coordinate rounding. Length checks
use the exact rational values of returned binary64 positions. Face reports give
an aspect-ratio bound and the requested certified angular bound. Unmet aspect
requests appear in `Mesh::unavailable_controls`, with the face and reason.
Exact planar trims retain their incidence through shared dyadic edge refinement.
`TracePiece::evaluate_exact` retains the algebraic positions of boundary samples
before display rounding, including samples introduced by length controls.
The aspect bound is `None` when rounding collapses a facet to zero area; requested
aspect and angular controls then appear as unavailable rather than finite bounds.
Planar triangulation proposes indexed CDT connectivity, then certifies the exact
convex hull, every triangle orientation and the oriented trim-boundary incidence.
Collapsed or uncertified proposals retain the exact trim clipper. Outward-rounded
binary64 interval filters accept only proven signs; inconclusive predicates use
`Real` arithmetic.
Sine/cosine range reduction uses a stored 4096-bit Machin enclosure of pi with
outward truncation and a series fallback above that precision. The constant has
no numeric work charge or per-request cache state. Direct pi queries are unchanged.
Edge subdivision reuses its certified endpoint enclosures down the recursion.

Why: applications need exact geometric summaries and explicit display guarantees,
not quantities inferred from triangles or controls silently discarded.
Migration: mesh option literals should use `..Default::default()` for optional
controls. An edge midpoint is halfway along arc length, including closed edges
with their start vertex or carrier seam. Bounds contain no display chord padding;
an empty body has no box. Existing mass-property requests, question defaults,
answer formats and STEP definitions are unchanged.

