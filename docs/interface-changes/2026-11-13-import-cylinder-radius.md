## 2026-11-13: a full cylindrical import band repaired to its unchanged rims

Pieces: 8, STEP and IGES repair records

After both the exact import and the established binary-noise cleanup refuse,
`apply_source` can refit the radius of one complete capped cylinder. The file
must contain one solid with one cylindrical face and two perpendicular planar
caps. Each cap has one complete circular rim. The two rims must be coaxial,
have the same positive radius and use the same source edges as the cylinder's
two physical boundaries. A repeated straight seam is allowed. The cap planes,
curves, vertices, face names and source bytes stay unchanged.

An inconsistent source has no exactly incident trimmed cylindrical face.
`catallax-import-repair/4` states its interpretation explicitly as
`cylinder_between_cap_planes/1`: the complete source cylinder at its declared
radius between the two source cap planes. Its equal-angle, equal-height
correspondence with the repaired band has distance exactly the radius
change everywhere. This proves both Hausdorff directions over the entire
band, including its rims. The one changed face gets that nonzero local bound.
Caps have zero deviation and get no added exception.

The radius change must be at most the unchanged 0.01 mm automatic cap and at
most a tenth of the minimum of the source radius, result radius and cap
separation. Those are conservative feature bounds for these complete disks
and band. A breach returns `ValidationFailed` with the measured radius change,
cap and squared feature bound. Arbitrary curved trims, unequal rims, several
solids, off-axis rings, cones, spheres, tori and spline faces retain their
existing refusals through this route. A missing bound never admits a repair.

The candidate must pass the ordinary exact incidence, self-intersection and
solid checks, and have positive certified volume. It builds in a scratch model
and commits only on success. `/4` records add `cylinder_radius`, containing the
source interpretation, source and repaired radii and unchanged rim edge keys.
The derivative hash and whole-face tolerance bind the repair to the immutable
source. `rebuild_source` remeasures and compares the interpretation as well as
the bounds and derivative. Forged records leave the caller's model unchanged.
This automatic route creates no acceptance proposal or looser-tolerance override.

Existing exact and binary-cleaned successes keep their records, work and
precedence. Their `/2` and `/3` JSON omits the empty new field. Callers persisting
or reopening `/4` records need this reader version; older readers reject its
format rather than silently applying another source interpretation.

Public synthetic checks cover inward and outward changes at exactly 10 µm,
an over-cap 10.001 µm change, a 6 µm change on a cylinder only 50 µm high,
a posed STEP copy, IGES replay, forged records and identical records,
derivatives and checked STEP bytes at 1, 2 and 4 workers. Raw evidence retained
privately. This qualifies a complete-cylinder family, not general curved gap
closure or spline-face import. No modelling tolerance or work allowance changes.
