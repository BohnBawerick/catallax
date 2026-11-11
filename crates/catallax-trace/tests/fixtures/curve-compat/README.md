# Curve-compatibility gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the routines that make curves
share one knot vector before they are skinned or ruled: `N_CrvElevateDegree`,
`N_CrvRemoveKnot`, `N_GetCompatibleKnotVectorToTol` and
`SmBSplineCurve::SyncronizeKnotsOfCurves`, which runs
`N_CrvsMakeCompatibleAdjKnots` (rational and dimension compatibility, scaling to
the unit interval, degree elevation, removal of a last span no longer than the
tolerance, knot merging and refinement). The library build is the one described
in [`../foundation/README.md`](../foundation/README.md). Catallax never
compiles or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. A curve is read and printed as stored: dimension,
degree, counts, every knot and each control point's four stored coordinates.
This gate accepts no unit-in-the-last-place budget: every binary64 output is
bit-identical.

The cases cover degrees 1 to 4, polynomial and rational curves in two and three
dimensions, Bezier curves, interior knots of every multiplicity up to the degree
(`tangent`) and elevation up to degree 17, the limit of SMLib's tables. Knot
removal runs on curves built by inserting a knot (removable), on genuine knots
at tolerances from 1e-12 to 1e300, on parameters that are not knots, on the end
knots and a parameter 1e-13 from a knot, on a curve whose weights would leave
`[1e-3, 100]`, and with a request for more copies than the knot has. Knot
merging runs on two and three vectors with shared, repeated, missing and
near-equal knots (within 1e-9 and 0.02), on gaps below the tolerance, and on
end knots that differ. Synchronisation runs on curves of different degree,
rationality, dimension and parameter range (scaled to the unit interval), on
ranges that differ by 1e-9 and 1e-7, on identical curves, on a single curve, and
on a last span of 1e-9, where SMLib removes the knot or, when
the removal reports its basis search as failed, refuses the whole call.

One group is not compared value for value. SMLib keeps a planar curve's
dimension label of 2 when synchronisation makes its data three-dimensional;
a checked `BSplineCurve` takes its dimension from whether its points carry z.
The three `sync-mixed-dimension` cases are pinned by their own test in
`../../curve_compat.rs`, which checks that every value is bit-identical and the
port reports 3 where SMLib keeps 2.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/curve-compat/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/curve-compat/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
