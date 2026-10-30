# Countersink placement invariance

Separated countersinks and neighbouring cylindrical faces no longer require
branch tracing merely because a part is translated or rotated. The Boolean
fallback retains its existing support-box proof, then tries parallel-axis tube
bounds for certified finite cone bands. It also excludes two separated finite
cone bands. The bounds require strict separation, so tangencies and actual
meetings cannot pass this exclusion.

Cone face probes and their adjacency paths now share an exact normalized frame
even when a rational seam direction has an irrational length. Other surface
probe frames keep their existing compatibility boundary.

The public P5 regression builds all four countersinks at the identity, under a
translation, under an axial rotation, and under a general rational rotation
plus translation. Each placement checks the independently derived volume
`11932*pi/3`, the surviving housing origins, and one countersink and wall face
per hole. Neighbour tests include the application's corner placement and its
countersinks-then-boss operation order. The moved P5 STEP regression checks raw
and standard OCP 8 validity, a single solid, and the independently derived volume.
Boundary tests preserve touching and crossing cone/cylinder and cone/cone pairs.

This does not implement general cone pair tracing, the full P5 junction fillet,
or spline exchange. The SolidDraft pin, scenario extension and raw-reader gate
remain a separate companion change.
