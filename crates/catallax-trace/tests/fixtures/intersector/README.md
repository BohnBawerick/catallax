# Surface intersector gate fixtures

`cases.txt` holds outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the stages of
`SmAdvSurfaceIntersector::DoIntersection` ported so far, read through a probe
subclass of the intersector: the constructor's tolerance, closed flags and
clipped domains; the root-box test; the `SM_SO_INTERSECTION_TEST` search with
the advanced intersector's overrides (the union of the meeting boxes and
domains of the leaf pairs, and their count); and `DoBoundaryPlaneTest`. The
library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

The cases are the eleven surface cases of the M25 scope study's coverage
driver (a bicubic bump and its rational copy against level and tilted planes,
with the default and a tight tolerance; a hill against planes through its
plateau, at its top and 1e-6 below it; a rational cylinder patch against a
wall along its ruling and across it; two bumps; a loft against a sloping plane
with two tolerances), and variations: planes missing the surfaces and touching
them at a corner, an edge and the top, sub-domains, a requested domain beyond
the bump's natural one, a smaller plane domain, swapped surfaces, walls 1e-6
off the cylinder's ruling, the patch with a pole against planes through and
below its pole, more levels, and flat patches side by side with a gap of 5e-4,
2e-2 and 1e-9, sharing an edge, overlapping, at an angle and turned 0.005 and
0.02 radians about z. Each case runs the search alone and with the
boundary-plane test.

When two side planes graze, SMLib intersects the side curves with the curve
intersector, which is not ported; the port refuses with
`SolverRefusal::NotPorted` there and the driver records `-1`, as for disjoint
root boxes.

Binary64 values match within 16 units in the last place; every one matched bit
for bit when the fixture was made. Swapping the subdivision order of
`sm_FindLargestNonLeaf`, dropping the tolerance expansion of the pseudo-box
test, changing the default tolerance, removing the disjoint outcome of the
plane test or the clipping of a requested domain to the natural one changes
cases. Widening the opposed-normal test of the side planes from 1e-4 to 1e-2
leaves every output unchanged: the turned patches' planes do not separate the
surfaces either way.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/intersector/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/intersector/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
