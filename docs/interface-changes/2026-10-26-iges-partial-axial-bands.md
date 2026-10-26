## 2026-10-26: IGES finite axial cylinder and cone sectors

Piece: catallax-iges

Change: `write_iges` accepts one finite four-edge axial cylinder/cone band bounded
by two coaxial circular arcs and two generators, when the angular endpoints are
cardinal in the retained support frame and the existing exact numeric encodings
apply. Signed 90/180/270-degree spans, reversed uses and seam crossings produce
one continuous parameter-space rectangle. Shared input edges and face counts
are retained; these sectors add no seam or vertex. Generators use exact bounded
degree-one type 126 curves. Finite type 100 arcs start at zero in an exactly
reparameterized orthonormal circle frame.

The report retains its whole-domain positive-similarity proof and adds the
sector arc-domain endpoint error to the conical slant error before common
scaling: `abs(scale - 1) * coordinate_bound + scale * trim_deviation`.
`trim_deviation` includes angular compatibility error as well as slant rounding,
so it can be nonzero for cylinders. Inclusive tolerance checks and all geometric
refusals still precede any sink write. Direct unhealed OpenCASCADE transfer and
sampled pcurve comparisons remain independent numerical evidence, not that proof.

Migration: no signature, option, modelling tolerance, work allowance or STEP
change. Callers that previously received `UnrepresentableDefinition` for the
covered bands can receive a certified IGES file. No application export command,
importer, sphere/torus support, apex cap or arbitrary angle encoding is added.
Other periodic trim shapes and unencodable definitions retain typed refusals.
