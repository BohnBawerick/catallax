# Curve-construction gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the curve construction and
knot-editing routines the blend and fillet surface generators build on:
`SmBSplineCurve::CreateLineSegment`, `CreateDegenerateCurve` (and
`CreatePointCurve`), `InsertOneKnot` (through `N_CrvInsertKnot`), `InsertKnots`
(through `N_CrvRefine`) and `SmCurve::EditParameterization` (through
`N_BasisReparam`). The library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles
or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. A curve is read and printed as stored: dimension,
degree, control-point and knot counts, every knot and each control point's four
stored coordinates, so the weight and two-dimensional markers are compared too.
This gate accepts no unit-in-the-last-place budget: every binary64 output is
bit-identical.

The cases cover degrees 1 to 4, polynomial and rational curves in two and three
dimensions, Bezier curves and curves whose control points all coincide.
Insertions run at fresh parameters and at existing knots up to the multiplicity
limit (`tangent`), one past it (`singular`), within 1e-13 and 3e-12 of a knot
(`near-singular`, where the constructor's `sm_FixupKnotVector` merges or keeps
the knot), and on or outside the end knots. Refinements run with one to four
new knots, repeated and existing knots, knots a hair from existing ones, and
lists that touch or pass the ends. Reparameterisations run to new intervals,
reversed intervals, intervals within SMLib's scaled zero of the old one (kept
as they are) and intervals of width 1e-9. Lines run at the 1e-12 length
threshold from both sides.

Three groups of cases are not compared value for value, because SMLib returns
a result the checked `BSplineCurve` cannot hold. Each is pinned by its own test
in `../../curve_build.rs`, which checks the port's refusal or value and the
SMLib result it differs from:

- `reparam-nonincreasing` (24 cases): a reparameterisation to a single value
  succeeds in SMLib and returns a curve whose active knot interval is empty;
  the port refuses with `InvalidInput`.
- `refine-unordered` (16 cases): `N_CrvRefine` documents an ascending list. A
  descending pair succeeds in SMLib and returns decreasing knots; the port
  refuses with `InvalidInput`.
- `point-planar` (2 cases): `CreateDegenerateCurve` with dimension 2 goes
  through `CreateCanonical`, which keeps the point's stored z; the port stores
  the two-dimensional marker. Every other value is bit-identical.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/curve-build/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/curve-build/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
