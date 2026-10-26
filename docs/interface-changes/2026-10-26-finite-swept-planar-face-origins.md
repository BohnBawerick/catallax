## 2026-10-26: finite ownership of swept planar face origins

Pieces: catallax-field, named-face query callers.
Change: post-Boolean planar sweep origins are associated only with material
terms whose finite boundary facet exactly matches an owning physical face tile.
Disconnected coplanar rims no longer answer each other's named-face queries.
Why: support and outward sense did not distinguish two holes on the same plane.
Migration: no signatures, formats, tolerances or allowances change. Explicit
origin selectors retain correct-hole membership and reject the other hole.
Every physical tile must have exactly one proved owner before publication;
unproved correspondence fails closed. This corrects the ownership behavior in
[the original planar sweep fragment](2026-10-25-post-boolean-polygon-face-thicken.md).
