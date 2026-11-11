# Fillet-trace gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmSurfaceIntersector::DoPointIntersection`
as `SmFilletIntersector` runs it: the through points' values, the trace through
them with the virtual `ComputePointValues`, `ComputeNextPoint` and
`TestSpanAccuracy` of the fillet intersector, and the optional trace before the
start. The driver subclasses the fillet intersector and replaces `FlushCurve`
(which builds the fillet surface and its rails) with a dump of the traced points
that it would flush. The library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. The inputs are the head of
[`../fillet-intersector/README.md`](../fillet-intersector/README.md) (radii,
orientations, solver tolerance, two surfaces and two domains), then the request:
the two extension flags, an optional start direction and end direction, an
optional distance tolerance and angle tolerance, the extension distance and the
through points as surface parameter pairs. The outputs are the status, a flag for
a clipped step (below), whether `FlushCurve` ran, the traced points (kind,
parameters and parameter tangents on both surfaces, curve point and tangent, curve
parameter, deviation and tangent-plane angle), whether the curve closed, and each
through point's kind and curve parameter, which SMLib writes into the caller's
points.

The surfaces are a flat floor and wall, a floor and wall stretched by a square
root, and a floor and wall cubic in the curve's direction, each with the exact
parameters of the offset surfaces' meeting for radii equal and unequal. Requests
trace two points with a start direction (forward and reversed, so the trace
direction flips), several points with and without directions and with each pair of
extension flags, one point, an end direction, looser and tighter tolerances, a
duplicated through point and an extension distance. The tangent cases use a floor
and ceiling whose offset surfaces are parallel; the touching cases two coincident
floors; the singular cases no points, a first point off the surfaces, zero and
1e-13 radii and a zero start direction. Spans inside the domains end on their end
point.

Every binary64 output of the compared cases is bit-identical; the gate accepts no
unit-in-the-last-place budget. SMLib's default point holds undefined vectors where
the port's holds zeros, so the harness builds points with undefined vectors.

**Clipped steps.** A trace that runs on to a side of a domain steps, near it, a
step that the domain clips. SMLib re-solves that step on the clipping side with
`LocalCurveIntersect` on an offset surface, as the port does. The driver records
whether any step was clipped without ending on the boundary, and the port reports
the same (`FilletTrace::clipped_step`). All 131 such traces are compared bit for
bit (a test keeps them in the gate), including every boundary case. In most of them
a rail curve after the clipped step leaves the domain by rounding (a parameter of
`1.0000000000000002` or `-2.8e-17`), where the offset surfaces extend their bases as
SMLib's do.

The tracing loop is the one of `surface_trace.rs`, which now takes any stepper
that implements `TraceSteps`; the fillet intersector is the second. `FlushCurve`
is not carried.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository. The driver is
compiled with `-fno-access-control` to reach the intersector's protected members.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 -fno-access-control \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/fillet-trace/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/fillet-trace/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 20 s CPU cap and a 40 s timeout. Raw evidence retained privately.
