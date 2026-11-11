# Circle-arc and end-point gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmBSplineCurve::CreateCircleSegment`
(through `N_CreateCircArc`: the quadratic arc and the quintic circle) and
`SmBSplineCurve::EditEndPoint`. The library build is the one described in
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

Two groups of cases are not compared value for value. Each is pinned by its
own test in `../../curve_arc.rs`:

- `circle-quintic-partial` (7 cases): a quintic arc short of a whole circle is
  cut from the circle by `N_CrvClosestPt` and `N_CrvExtractCrvSeg`, which are
  not ported; the port refuses with a degenerate fault where SMLib returns the
  arc. The whole quintic circle is compared.
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
