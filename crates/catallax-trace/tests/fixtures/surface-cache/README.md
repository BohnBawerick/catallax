# Surface cache gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the trimmed-surface cache that
`SmCacheMgr::GetOrCreateObjectCache(SM_OC_TRIMSRF, surface)` builds for a
B-spline surface without a face (its tolerances from
`SmSurface::CacheMakeOrValidate`), for `sm_ComputeNetConstants` on the control
net as the cache asks for it, for `SmSurface::CalculateContinuities` in u and
v, for `SmBSplineSurface::CalculateBoundingBox` with its pseudo box, and for
the side planes of `SmPatchBoundaryBoundingPlanes` over the natural domain
(each side's plane, degeneracy, planarity and samples). The
library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

A cache case lists the knot continuities in u and v and then the whole tree in
preorder: each node's auxiliary data type and box, then for a Bezier patch node
its parameter domain, pseudo box, side continuities, split direction and
tolerance flags, and for an organising node its domain and split direction.

The surfaces are those of the surface point solver gate (bumps, a rational
loft-like patch, a pole, a bilinear patch, a closed cylinder) and a bicubic
with a triple interior knot whose net folds there, planes whose
parametrisation kinks along a knot with the normal unchanged and with the
normal turned by half a degree (SMLib raises the first kink to G1 and not the
second), biquadratic waves with nine knots across (SMLib samples five), a strip
a thousand times longer than wide, where subdivision stops at the smallest side
ratio, a doubly folded patch, a net whose diagonal control points coincide, a
net with control points 1e-7 apart (below the zone tolerance the turning angle
skips), and flat zigzag nets of degrees 5 to 8 (in u, in v and in both, rational,
with a small height, and across two patches).

The zigzags turn more than 540 degrees, where SMLib's angle test falls back to
the surface's own turning, `sm_ComputeSurfConstants`: it samples the whole
surface at five by five points over the patch's domain and sums the angles
between successive chords along each row and column. A second family of nets goes
forward and back along a nearly straight line: the net turns 180 degrees at every
vertex and the surface hardly turns, so the sampled angles are small where the
net's are large and the tree differs from the one the net's angles would give.
Without the sampled angles eight of the cases fail. Only the angle outputs of
`sm_ComputeSurfConstants` are carried, which is all the cache asks for.

Binary64 values match within 16 units in the last place. Every value matched
bit for bit when the fixture was made except one control-net turning angle on
the pole, 1 ulp apart: the sum of arc cosines is the same in both, and the C
library's `acos` that SMLib calls and the `libm` one the port calls differ in
the last place for one of the terms. Mutating the split continuity, the
normal-agreement threshold, the box pairing direction, the chord-height
fraction, the zone tolerance in the turning angle, the rough turning-angle
scale or the 540-degree limit changes cases. The sampling fractions across
knots and the exhausted diagonal searches only act on surfaces these cases do
not include; for a net whose whole control diagonal coincides the exhausted
search gives the same answer either way.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/surface-cache/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/surface-cache/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
