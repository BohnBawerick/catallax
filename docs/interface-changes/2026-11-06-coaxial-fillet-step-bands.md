# Complete coaxial STEP bands for fillet scenes

The complete cylindrical, conical and toroidal band encoder now also accepts
eligible scenes containing a toroidal fillet. Previously, a trimmed torus
needed a cone neighbour to enter this route. Torus-only fillet scenes were
sent to generic chart preparation even when every boundary was a complete
coaxial ring.

A torus face's recorded `blend` origin selects the new route. It does not
certify the geometry. The existing whole-scene checks still validate every
periodic face, shared ring, coaxial relation, incidence and band orientation.
An unhandled trimmed sphere or arbitrary periodic patch declines the entire
plan. It must not become legacy output just because another face is a band.
Unblended torus scenes retain their previous chart route and STEP bytes.

The corpus's radius-0.2 rounds on a torus cap now export within 1M work as one
raw-valid, four-face solid. Both material sides, the central hole and the
retained side of the cap have independent witnesses. The volume also matches
a separate washer integral of the rounded meridian section.

The native fillet definition, radii, senses, origins, overflow semantics,
17 significant digits and 1/10000000 mm file tolerance do not change. The
encoder uses analytic carriers and certified isoparameter curves, not a new
fitted carrier. Its file-deviation report remains subject to the caller's
requested tolerance.

No public API signature, recipe format, chamfer/draft logic, display mesh,
SolidDraft UI or dependency pin changes. Sphere-cap export, arbitrary partial
bands and general corner construction remain separate work.
