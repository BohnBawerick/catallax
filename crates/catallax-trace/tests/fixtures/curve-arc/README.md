# Circle-arc and end-point gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmBSplineCurve::CreateCircleSegment`
(through `N_CreateCircArc`: the quadratic arc and the quintic circle) and
`SmBSplineCurve::EditEndPoint`, and for the NLib routines the quintic arc cuts its
circle with: `N_CrvDerivs`, `N_CrvClosestPt` and `N_CrvExtractCrvSeg`. The library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles
or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. A curve is printed as stored: dimension, degree, counts,
every knot and each control point's four stored coordinates. The frame of a
circle is made by `SmAxis2Placement::SetCanonical`, as the library's callers do.

The circle cases use eight frames (axis-aligned, rotated, non-unit axes, axes
0.001 and 1e-13 off perpendicular), radii from 1e-12 to 1e12, two and three
dimensions, and sweeps on both sides of each span-count threshold (91, 181 and
271 degrees) and of the closure tolerance (359.99), negative and wrapped angles,
equal start and end, sweeps of 1e-7 to 2e-6 degrees (nearly parallel end
tangents, which take SMLib's flat-span branch), the angle limits of plus and
minus 360 and values beyond them, and refused radii and dimensions. The end-point cases
edit the first and last control point of polynomial and rational curves in two
and three dimensions with and without a weight.

Every binary64 output is bit-identical except 23 values in 5 cases, at most 3
units in the last place apart. Their cause is the sine and cosine: the port
takes them from `libm` for repeatable last bits on every machine, SMLib from
the C library. Swapping the C library's functions into the port makes all 7,356
outputs bit-identical, which is how the cause was established. The test pins
the count, so no other difference can join them unseen, and the gate allows
4 ulp.

The NLib cases run `N_CrvDerivs` on curves of degrees 1 to 4, polynomial and
rational, planar and spatial, at the ends, inside, from both sides, to the
derivative orders SMLib allows (four) and beyond; `N_CrvClosestPt` from random
points and starts, to tolerances 1e-7, 1e-12 and 1e9, and on closed curves;
and `N_CrvExtractCrvSeg` over the whole curve, over parts, at knots, 1e-13
from a knot, and over empty, reversed and out-of-range intervals. The quintic
circle cases include partial arcs (0 to 90, 0 to 180, -90 to 90, 30 to 300,
0 to 359, 350 to 20, 5 to 5.5, 100 to 99 and others), which are cut from the
circle through these routines.

One group of cases is not compared value for value. It is pinned by its own test
in `../../curve_arc.rs`:

- `edit-end-invalid` (6 cases): SMLib writes a zero weight or a non-finite
  coordinate into an end control point and succeeds; a checked `BSplineCurve`
  refuses both with `InvalidInput`.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/curve-arc/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/curve-arc/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
