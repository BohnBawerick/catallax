## 2026-10-23: filtered trace Taylor enclosures and curved display meshing

Pieces: 4, display mesh, application
Change: `TracePiece::taylor_interval` returns optional outward-rounded binary64 Taylor enclosures through the requested order. `None` requests the existing ball path; deadline errors still propagate. Immutable trace constants retain a binary64 enclosure; nonrational constants derive it at 128-bit working precision, and ball enclosures remain precision-specific. Display charts now include elliptic cylinders and filter evaluation, chord and aspect bounds before exact-rational or ball fallback. Fillet cylinder carriers retain exact squared radii, and display coefficients obtain radii through the definition accessor.
Why: repeated high-precision evaluation and incremental UV insertion blocked display latency. Shared-rail ears also created coincident sliver triangles between exact fillet strips.
Migration: existing trace and meshing callers need no API changes. Filter callers must handle `None` with the existing certified path. Display triangle counts, indices and conservative reported bounds may change within the requested certificate. The 0.01 mm default, vertex and work limits, exact B-rep, question formats and STEP definitions are unchanged.

