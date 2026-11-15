# Native cylinder spans and numerical Hole cap curves

`catallax_features::hole::HoleTool::cylinder_span()` returns one retained cylinder
and signed millimetre interval for a contiguous, constant-radius drill. Its origin,
compiled unit axis and seam come from the original Hole definition. Entry and wall
bands retain their individual intervals and `HolePart` provenance. Changing-radius
heads and conical tips return `None`.

`hole::cap_curve` accepts that span, one quadratic/linear cylinder quarter, a cap
plane and preset `NumericalOptions`. It returns a `CapCurve` containing a 3D edge,
plane UV curve, cylinder-quarter UV curve and numerical evidence. The caller's key
and the source quarter's finite U interval remain unchanged. Parallel projection
returns `None`; invalid sources, cancellation and failed checks return a refusal.

The projection moves each source pole along the compiled cylinder axis to the cap
plane. The 3D and plane UV curves keep the quadratic degree, native knots and
weights. The cylinder UV curve uses homogeneous degree elevation to represent
native U and the projected axial height in one cubic rational curve. Its V maps
to the supplied source interval. Clamped endpoint controls avoid interval widening
of a homogeneous quotient beyond an integer native U boundary. UV coordinates are
not clamped or reparameterized.

## Prior art

The pole projection follows OCCT's
[`ProjLib_ProjectOnPlane`, B-spline case](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_3/src/ProjLib/ProjLib_ProjectOnPlane.cxx#L951):
project the poles and retain the curve parameterization. The SMLib port map in
[THIRD-PARTY.md](../THIRD-PARTY.md) identifies `SmCrvOnSurf` evaluation in
`catallax-trace/src/curve_on_surface.rs`. Its UV-then-surface composition is the
model for the retained-parameter checks. No new SMLib code is ported here.

## Evidence and limits

Every cap records six source axial-control checks and 17 samples for each of:
distance to the original plane, distance to the original cylinder, plane UV
composition and cylinder-quarter UV composition. All use the caller's preset
physical tolerance. These samples establish no complete curve or body deviation
bound and grant no exact B-rep incidence.

The tests exercise 32 templates across top/bottom caps, axial/skew Hole frames,
four source quarters and a rigid placement. Independent surface-section evaluations
check the 3D and UV curves. Wrong radius, wrong or shifted source interval,
parallel-plane and cancelled-request controls cover refusal behavior.

A closed skew-cylinder fixture assembles two cap loops and four wall quarters
with 24 retained native trim bindings. A diameter edit keeps the stable face keys;
numerical handover preserves the bindings, tolerance and edit history. Live and
restored numerical STEP bytes agree. This fixture qualifies the cap/wall contract;
the collar Hole still needs the finite cone-face partitioner and app consumer.

The ordinary exact Hole evaluator is unchanged. Numerical assembly and STEP remain
separate typed results. The cap helper neither inserts fitted facts into the exact
model nor converts numerical evidence into an exact answer.
