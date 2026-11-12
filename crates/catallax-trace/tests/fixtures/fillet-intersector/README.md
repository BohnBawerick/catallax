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
degrees to 1e-4 radian, to reject spans by angle and by distance. The last pair
is curved: its floor's y runs as the cube of its parameter, its wall bends as
x = 1 + v^2, and their offset surfaces meet on a parabola. Its inputs are the
offset intersection at a wall parameter, so the two chords a `1e-5` step and a
full step give are not collinear.

The through-point repair changes the current tangent in six `fi-next` cases:
two on the cubic pair and four on the curved one. SMLib fills the short step's
point values at the full step's UVs, not at the short step's own UVs
(`ComputeNextPoint` passes `rNextTSP.UVPos` to `ComputePointValues`), so the
repaired chord is the full step's direction, divided by `1e-5` and unitized;
the small step supplies no direction of its own. The port keeps that behavior
for parity, and these six cases pin it: computing the values at the short
step's UVs instead makes them differ from SMLib.

Every binary64 output is bit-identical; the gate accepts no unit-in-the-last-place
budget. SMLib's default point holds undefined vectors where the port's
`TsectPoint::default()` holds zeros, so the gate's harness builds points with
undefined vectors, as SMLib does. `ComputePointValues` of the intersector
reports any failure as `SM_ERR`, and the port does the same.

**Pinned refusal, clipped steps.** A step that a domain clips and that does not end
on the boundary goes on, in SMLib, to `LocalCurveIntersect` on an offset surface,
which the port does not carry. The 31 such cases are `fi-next-clipped`: the gate
checks that the port refuses each with `StepRefusal::ClippedStep` and that SMLib's
answer for it carries the clipped flag.

**Pinned refusal, surface extension.** An offset surface in SMLib extends its base
surface when an evaluation lies outside the base surface's domain
(`sm_ExtendSurface`), which the port does not carry; it would clamp instead. In
five `fi-next` cases a Hermite parameter curve leaves the domain (three by about
2e-16, two by a little more on the cubic surfaces) and SMLib's rail points differ
from the clamped ones: one reports a deviation of 6.7e-16 where the clamped value is
8.9e-16. The port refuses these with `StepRefusal::SurfaceExtension`. The test lists
the five fixture lines, checks SMLib's answer for each, and fails if any other case
refuses for this reason.

The tracing loop that calls these methods, `FlushCurve`, and the extension of the
offset surfaces' domains are not carried.

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
