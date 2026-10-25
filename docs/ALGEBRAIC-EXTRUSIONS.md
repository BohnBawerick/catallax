# Signed extrusion with algebraic offsets

A finite extrusion accepts a rational model-space direction with a nonzero
component normal to the sketch. Its signed distance is a physical distance:
for direction `d` and distance `h`, the end offset is exactly `h*d/sqrt(d.d)`.
The unit vector need not be rational. Line, circle and circular-arc profiles
retain their analytic carriers, shared boundary edges and construction origins.

The B-rep builder keeps the existing rational path when `sqrt(d.d)` is rational.
The separate algebraic builder retains exact-real offsets for vertices, boundary
curves and cap planes. Side planes use the unnormalised rational direction;
circular sides reduce to cylinders or elliptic cylinders, with squared minor
radius `r²*(n.d)²/(d.d)` for the profile's unit normal `n`.

The independent field implementation uses a rational projective sweep direction
and an exact algebraic finite end parameter. Profile projection cancels the
normalisation. Its analytic surfaces and finite cap predicates do not consult
B-rep topology. Bounds enclose exact offsets; display remains a separate
approximation at the existing requested accuracy.

Native display of an isolated algebraic extrusion meshes the same upright
profile using independent field charts, then applies the exact invertible affine
sweep map and its final isometry. The map's Frobenius norm bounds amplification
of the source mesh's two-sided distance certificate. A tighter internal source
accuracy reserves half of the original requested tolerance for coordinate
rounding; complete vertex enclosures and facet-orientation checks certify the
transport. Caller tolerances, work and mesh limits remain unchanged. This path
does not mesh a B-rep or widen the general rational chart implementation.

Signed distance edits preserve `start`, `end` and segment `side` origins. Zero
direction retains `extrude.zero-direction`; a nonzero direction in the sketch
plane retains `extrude.zero-thickness`. Zero distance retains
`extrude.zero-distance`. Bounded first-hit extrusion is a separate implementation.

The rational candidate sampler cannot emit points on a general algebraic cap.
A sampled kernel comparison can therefore remain unresolved with the existing
algebraic sampling-coverage refusal. Exact membership and independently controlled
volume checks establish their stated properties; sampling is not a whole-body
certificate. General downstream draft/Boolean coverage is not implied.

Public regression inputs are in
`crates/catallax-features/tests/fixtures/algebraic-extrusions`. Rust tests check
analytic volume and containment, application properties and distance-edit names.
`tools/step-regressions/test_algebraic_extrusions.py` independently imports signed
line, circle and arc outputs with OpenCASCADE. Import validity and double-precision
mass checks are separate from the exact property-enclosure tests.
