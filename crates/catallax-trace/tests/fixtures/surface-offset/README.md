# Offset-surface gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmOffsetSurface`:
`EvaluatePoint`, `EvaluateNormal`, `Evaluate` and `EvaluateSimple` (through
`EvaluatePlainOrSimple`) on offsets of B-spline surfaces. The library build is
the one described in [`../foundation/README.md`](../foundation/README.md).
Catallax never compiles or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. A surface is read as stored (`du dv nu nv rational
knots_u... knots_v... (x y z w)...`). An evaluation prints every cell of the
`(hu + 1) x (hv + 1)` matrix, SMLib's undefined marker included where it does not
write a cell. This gate accepts no unit-in-the-last-place budget: every binary64
output is bit-identical.

The surfaces are fifteen: flat, gently curved and strongly curved patches of
degrees 1 to 3 in each direction, polynomial and rational, patches with a
collapsed control row at the start of U, at both ends of U and at the start of V
(poles, where an offset has no normal and SMLib extrapolates), and patches scaled
by 1e-6 and 1e6. The offsets are 0, 0.5, -0.5, 2 and 1e-9. The parameters are
interior, on every edge and corner, 1e-9 and 1e-12 off an edge, and outside the
domain (SMLib clamps; for the normal it then uses the offset surface's own
derivatives). Evaluations run to every order SMLib allows (a total of three),
the full and the upper-half matrix, simple and repaired tangents, and both sides
of each parameter, and include the orders above three that SMLib refuses.
The pole cases run at each pole, with positive and negative offsets, in both
evaluation kinds.

SMLib's evaluation has two quirks that the port reproduces and the gate pins:
the mixed derivative `S_uvv` and the third derivative `S_vvv` store the base
derivative `B_uuv` where `B_uvv` and `B_vvv` belong; and in upper-half mode it
writes cells beyond a non-square matrix, which the port keeps in scratch and
drops. SMLib's surface extension and nested offsets are not carried: the port
clamps to the domain, which is SMLib's behaviour with extension off.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/surface-offset/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/surface-offset/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
