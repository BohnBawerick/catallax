# Finite oblique annular separation bound

Hollow coaxial oblique draft now checks both finite radial laws and retains four
explicit cap intersections. `ObliqueAnnularBound` exposes the two checked ruled
supports, each cap plane with its outer and inner exact closed spline row, an
outward rational enclosure of the signed radial gap in millimetres, and an exact
separation relation. It contains no body or assembled cap faces.

In a common axial and angular frame, each radius is
`A(z) + b*cos(u) + c*sin(u)`. Their difference is affine in height. Its exact
minimum and maximum on the entire finite domain are the minimum/maximum endpoint
amplitudes minus/plus `sqrt((b_outer-b_inner)^2+(c_outer-c_inner)^2)`. The decision
uses exact signs; enclosure width never turns touching into separation. This is
a radial-order certificate, not a shortest Euclidean wall-distance measurement.
The constructor retains the existing positive, convex-section requirements on
each individual support. Unsupported folds and numeric exhaustion keep their
own refusals.

For the annular fixture with source radii 8 and 4 mm, caps at z=0 and z=2,
neutral plane `x+8(z-1)=0`, and both walls drafted by 45 degrees, the laws are
`rho_outer=9-z-cos(u)` and `rho_inner=3+z+cos(u)/2`. The full signed radial gap
is enclosed by `[1/2,15/2]` mm, with both extrema attained on cap intersections.
Changing the neutral-axis rate from 8 to 6 gives touching with gap `[0,8]` mm;
rate 4 gives crossing with gap `[-1,9]` mm. Zero rate has no closed transverse
neutral curve and is diagnosed before any division by that rate.

B-rep evaluation and field compilation still publish nothing for these hollow
oblique cases. A positive bound returns `draft.unsupported-oblique-annulus`,
with `separation=Separated`, `certified_radial_gap_mm=[lower,upper]`, four exact
cap intersection rows, `shared_trimmed_cap_rings=unbuilt`, and
`body_published=false`. Cap intersection curves do not certify assembled trimmed
cap faces or a closed material body. Touching, crossing, reversed radial order,
and zero neutral-axis rate respectively return
`draft.oblique-annulus-touching`, `draft.oblique-annulus-crossing`,
`draft.oblique-annulus-reversed`, and
`draft.oblique-annulus-zero-neutral-rate`. Existing source dimensions, wall
selection, cap span, tangent and neutral tilt diagnostics remain available.

Zero-angle identity and supported axial annular construction retain their
existing paths. No recipe version, modelling tolerance, work allowance or
publication contract changes.
