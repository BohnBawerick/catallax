# Fillet-intersector gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the stepping half of
`SmFilletIntersector`: its constructor's tolerances, `ComputePointValues`,
`ComputePointValuesDelta`, `TestSpanAccuracy`, `ComputeNextPoint` and
`EvaluateLawPoint`, over a `SmConstantRadiusFS`. The library build is the one
described in [`../foundation/README.md`](../foundation/README.md). Catallax never
compiles or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. The head of the inputs is `r1 r2 outward1 outward2
solver-tolerance`, the two surfaces (each read as stored: `du dv nu nv rational
knots_u... knots_v... (x y z w)...`) and the two domains. Every operation but
`fi-ctor` then sets the intersector's members, which are protected in SMLib:
`m_dThisApproxTol3d`, `m_dThisAngTolRad`, `m_dCurveTraceDirection` and
`m_bClippedByThroughPoint`. The operations are:

- `fi-ctor`: the tolerance, angle tolerance, trace direction and through-point flag
  after construction.
- `fi-values`, `fi-delta`: `ComputePointValues` with an optional previous point,
  and `ComputePointValuesDelta` with the two curve parameters that fix the chord's
  direction.
- `fi-law`: `EvaluateLawPoint` for a plane and guesses.
- `fi-span`: `TestSpanAccuracy` between two points that `ComputePointValues` makes
  from parameter pairs, with an optional middle point.
- `fi-next`: `ComputeNextPoint` from such a point and a step; the output includes
  the current point's tangents, which the through-point branch can change.

The surfaces are a flat floor and wall, polynomial and rational, flat and with
control-point perturbations of 0.06 and 0.08; a floor and wall stretched by a
square root, so the parameter speed changes along the curve; and a floor and wall
cubic in the curve's direction, where it changes more than threefold. Radii are
equal and unequal. The tangent cases use a floor and ceiling whose offset surfaces
are parallel; the singular cases use zero and 1e-13 radii and a zero chord; the
touching cases use two coincident floors; the boundary cases step from the ends of
the domains. Tolerances run from 1e-3 to 1e-12 and angle tolerances from twenty
degrees to 1e-4 radian, to reject spans by angle and by distance. Six `fi-next`
cases take the through-point branch.

Every binary64 output is bit-identical; the gate accepts no unit-in-the-last-place
budget. SMLib's default point holds undefined vectors where the port's
`TsectPoint::default()` holds zeros, so the gate's harness builds points with
undefined vectors, as SMLib does. `ComputePointValues` of the intersector
reports any failure as `SM_ERR`, and the port does the same.

**Clipped steps.** A step that a domain clips and that does not end on the boundary
ends on the clipping side, where SMLib re-solves it with `LocalCurveIntersect` on an
offset surface and the rest of the span uses the rescaled step. There are 31 such
`fi-next` cases; 30 are compared bit for bit (a test keeps them in the gate) and one
is pinned below.

**Surface extension.** An offset surface in SMLib extends its base surface when an
evaluation lies outside the base surface's domain (`sm_ExtendSurface`). In five
`fi-next` cases a Hermite parameter curve leaves the domain (three by about 2e-16,
two by a little more on the cubic surfaces) and SMLib's rail points differ from
clamped ones: one reports a deviation of 6.7e-16 where a clamped value is 8.9e-16.
These five were pinned as refusals while the extension was not ported; with it
ported ([`../offset-extend/README.md`](../offset-extend/README.md)) they are
compared like the rest and match bit for bit. So are six clipped steps whose rail
curves leave the domain by rounding after the step ends on the clipping side.

**Pinned difference, an evaluation the extension cannot reach.** In one clipped
`fi-next` case (fixture line 605) the next point carries parameter tangents of about
3.3e5, and the rail curve's Hermite parameter curve reaches `v = -49.9`, fifty
domain lengths outside. SMLib's offset surface fails to evaluate there after
extending (a fresh `SmOffsetSurface` with extension on answers `SM_ERR` at that
parameter), and `SmCrvOnSurf::Evaluate` then falls back to no derivatives and leaves
the tangent unset, which `ComputeNextPoint` reads uninitialised; it reports no good
point. The port refuses with `StepRefusal::Failed` instead, as its `lift` does for
that fallback. The test checks the refusal and SMLib's answer.

The tracing loop that calls these methods is gated in
[`../fillet-trace/README.md`](../fillet-trace/README.md). `FlushCurve` is not
carried.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository. The driver is
compiled with `-fno-access-control` to reach the intersector's protected members.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 -fno-access-control \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/fillet-intersector/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/fillet-intersector/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
