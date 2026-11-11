# Curve-fit gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmBSplineCurve::CreateInterpolatingCurve`
through points with a vector at each (`N_FitCrvFirstDeriv`, with the banded
LU solve, the parameter values and the knots under it), the two-point line
shortcut, and `SmBSplineCurve::FairBlendCurveDerivatives`. The library build is
the one described in [`../foundation/README.md`](../foundation/README.md).
Catallax never compiles or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. A curve is printed as stored: dimension, degree, counts,
every knot and each control point's four stored coordinates, then the input
vectors as SMLib leaves them (with `bIgnoreVectorMagnitude` it normalises the
caller's array and scales it by the data polygon's length). The fit is plain
arithmetic and a banded solve, so the gate allows no ulp budget: every binary64
output is bit-identical.

The interpolation cases cover the uniform, chord-length and centripetal
parameterisations, degrees 2 and 3, two and three dimensions, three to thirty
points, vectors as derivatives and as directions, straight data, vectors along
the path, data scaled from 1e-6 to 1e8, repeated and nearly repeated points,
data with no length (SMLib leaves the interior parameters unwritten and its
allocator returns zeroed memory in the pinned build, which the port copies),
zero and tiny vectors, refused degrees (0, 1, 4, 5, 17) and dimensions, and the
line shortcut. The fairing cases cover random end conditions with and without
inflection prevention, parallel, opposite and chord-aligned tangents, S-shaped
blends, coincident points, zero and tiny tangents, and two, four and odd counts
of higher derivatives.

Two groups of cases are not compared value for value. Each is pinned by its own
test in `../../curve_fit.rs`:

- `interpolate-planar` (48 cases): SMLib builds a two-dimensional curve from the
  three-dimensional fit and labels it two-dimensional, keeping the points' z. A
  `BSplineCurve` takes its dimension from whether its points carry z, so the port
  stores the planar marker. Every other output, including the vectors, is
  bit-identical.
- `interpolate-other` (19 cases): the other branches of `CreateInterpolatingCurve`
  (end vectors, a start vector, no vectors, two points) use `N_FitCrvHighDerivs`,
  which is not ported. The port refuses where SMLib returns a curve.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/curve-fit/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/curve-fit/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
