# Offset curve-surface gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for SMLib's local curve-surface
solvers on offset surfaces: `SmSurface::LocalCurveIntersect` and
`SmSurface::LocalCurveSolve` called on an `SmOffsetSurface` with an `SmIsoCurve`
of another `SmOffsetSurface` (the solve `SmFilletSolver::StepSolve` makes for a
clipped step), and `SmSurface::IsSingularity` and `SmSurface::GetSingularities`
on an offset surface. The library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. A surface is read as stored (`du dv nu nv rational
knots_u... knots_v... (x y z w)...`), each after its offset distance. The
operations are:

- `oc-intersect`: offset surface A's iso-curve (axis, parameter, whether the
  parameter is clamped to the natural domain; from the left, as `StepSolve`
  builds it) against offset surface B over a domain, with a tolerance, a curve
  parameter guess and a surface parameter guess. The output is the status, found,
  and the curve parameter, surface parameters and deviation when found.
- `oc-solve`: the same curve and surface with `LocalCurveSolve` minimising or
  normalising from the guess.
- `oc-sing`: `IsSingularity` at a parameter pair with a tolerance, both forms, then
  `GetSingularities` with its pole points. `SmOffsetSurface` hides that virtual with
  a cached flag of the same name, so the driver calls `SmSurface`'s, as
  `IsSingularity` does.

The surfaces are a floor and a wall (flat, bumped by 0.06 and 0.08, and rational)
whose offsets meet as a fillet's do, an arch whose offset crest the floor offset's
iso-line grazes (tangent), touches within 1e-7 (touching) or crosses twice, and
fans with a collapsed side at `u = 0`, `v = 0` and `u = 1` (singular, near-singular
at 1e-9). The offset distances are positive, negative and zero (an offset of a pole
is a curve, not a point, unless the distance is zero). The boundary cases put the
iso-curve on and past the edges, clamped and not (the offset surface clamps the
evaluation, which is SMLib's behaviour with surface extension off), with guesses
outside the domain and sub-domains that hold the meeting and that miss it.

Every binary64 output is bit-identical; the gate accepts no unit-in-the-last-place
budget.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/offset-curve-surface/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/offset-curve-surface/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
