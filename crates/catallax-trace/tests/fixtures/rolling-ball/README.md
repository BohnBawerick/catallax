# Rolling-ball gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the constant-radius rolling-ball
equations: `SmConstantRadiusFS` with `SmFilletSolver::ComputeSurfaceValues`,
`ComputePointValues` (which runs `ComputeUVVectors`) and `StepSolve`. The library
build is the one described in [`../foundation/README.md`](../foundation/README.md).
Catallax never compiles or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. The head of the inputs is `r1 r2 outward1 outward2
tolerance` and the two surfaces, each read as stored (`du dv nu nv rational
knots_u... knots_v... (x y z w)...`). The operations are:

- `rb-surface`: the surface values at a parameter pair on each surface. The
  output is the point kind, then per surface the position, `S_u`, `S_v`, unit
  normal, parameters and parameter tangent, then the curve point, curve tangent
  and tangent-plane angle. SMLib's undefined marker is kept where it does not
  write a field.
- `rb-point`: the same point after `ComputePointValues`, with a trace direction,
  both parameter domains and an optional previous point whose tangent fixes the
  direction at a tangent point.
- `rb-step`: `StepSolve` from a point with parameter tangents, a step and both
  domains. The output is the found, clipped and boundary-hit flags, the step
  SMLib leaves in its reference argument, and the next parameters when found.

The surfaces are a flat floor and a wall, polynomial and rational, flat and with
control-point perturbations of 0.06 and 0.08, with radii equal and unequal and
all three orientation pairs. The tangent cases roll between a floor and a ceiling
at radii that sum to the gap, where the offset surfaces are parallel and SMLib
steps to find the tangent. The singular cases use zero and 1e-13 radii and two
coincident floors; the touching case is the coincident pair at `rb-point`. The
boundary cases step from the ends of the domain.

Every binary64 output is bit-identical; the gate accepts no unit-in-the-last-place
budget. One difference of representation remains: SMLib's default point holds
undefined tangents where the port's `TsectPoint::default()` holds zeros, so the
gate's harness builds points with undefined tangents, as SMLib does.

**Pinned refusal.** A step that a domain clips and that does not end on the
boundary goes on, in SMLib, to `LocalCurveIntersect` on an offset surface. That
needs the M9 local solvers generalised over a surface trait, which this port does
not carry. The four such cases are `rb-step-clipped`: the gate checks that the
port refuses each with `StepRefusal::ClippedStep` and that SMLib's answer for it
is found with the clipped flag set. They are not excluded; the refusal is the
pinned difference. Clipped steps that end as boundary hits do not reach that
branch and match bit for bit.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/rolling-ball/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/rolling-ball/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
