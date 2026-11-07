# Exact rolling-ball bands on oblique cone sections

## Scope

A fillet on the closed ellipse where a plane cuts a cone obliquely now builds
one named band instead of refusing with `fillet.unsupported-contact`. The band
is the exact rolling-ball surface, a rational carrier of degree (6, 2): the
plane and oblique cylinder band of degree (4, 2) with one more factor for the
cone's centre path, `1 + e cos t` for the section's polar eccentricity `e`.
Its plane rail stays an exact ellipse. Its cone rail is the carrier's last row,
exact but no conic, so the fillet's cone contact and the cone face's new edge
are that `bspline` curve.

## Consumers

- Point classification proves points whose rays cross the band. The band's
  line crossings now solve over algebraic coefficients, which also covers
  cylinder bands with irrational laws that previously had no line equations.
- `RealPoly::real_roots` lowers a polynomial whose coefficients lie in one
  quadratic tower into the existing tower arithmetic. It takes the rational
  norm one generator at a time, as `(A + sqrt(r) B)(A - sqrt(r) B)`, isolates
  the norm's roots, and keeps a candidate exactly when the polynomial's sign
  there is zero. That sign comes from the signs of `A`, `B` and
  `A^2 - r B^2`, down to rational polynomials. Other algebraic coefficients
  keep the general eliminant.
- Volume, area, centroid and mass properties answer. Spline curve traces now
  expand each knot span about its own start, which keeps interval enclosures
  of spans far from parameter zero near the curve.
- Display meshes close.
- STEP writes the band as its written carrier and the cone rail as the same
  row, with exact rational power-form pieces for the rail's chart curve on the
  `CONICAL_SURFACE`. Only a cone face trimmed by such a rail takes the chart
  path; every other cone keeps its legacy topology and bytes.
- The CLI returns `catallax-manifest/5` for a result whose spline faces are all
  exact rolling-pipe bands with at least one on a cone.

## Observed coverage

The corpus case `oblique-cone` (frustum of radii 10 and 4 over height 24, cut
at `2 atan(7/24)` through `(0, 0, 8)`, radius 1/2) builds at 23.3M work in the
fillet and writes STEP at 9.9M more. At the corpus's 1M allowance it refuses
with `numbers.operation-budget` inside the fillet. Its raw STEP imports into
OCCT as one valid solid of volume 1957.4896872686613 mm^3, against the exact
Steiner value 1957.48968726864861 mm^3; OCCT's own fillet gives
1957.4901903688324 mm^3.
