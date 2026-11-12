# Offset-surface extension gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmOffsetSurface` with
extension allowed (`SetAllowExtension(TRUE)`, as the fillet solver sets it):
`EvaluatePoint`, `EvaluateNormal`, `Evaluate` and `EvaluateSimple` at
parameters outside the base surface's domain, which extend the base surface
(`sm_ExtendSurface`, over the NLib routines of
[`../surface-extend/README.md`](../surface-extend/README.md)), wrap the
parameter where the base is periodic and stop it at a singular side; the
extension kept between evaluations; and what the constructor works out of the
base surface (`IsPeriodic` in each direction and `GetSingularities`). The
library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code. The driver is compiled with `-fno-access-control` to read
the offset surface's members.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. The only operation is `ox-seq`: an offset distance, the
extension flag, a base surface (read as stored: `du dv nu nv rational knots_u...
knots_v... (x y z w)...`) and a sequence of steps on that one offset surface.
A step is a point, a normal or an evaluation of given orders and flags at a
parameter pair. The output is the periodic flags and singular sides the
constructor found, then for each step its status and values, the surface's
natural domain afterwards (which is the extension's once one is held) and
whether an extension is held, and at the end the extension itself.

The base surfaces are a plain bilinear patch, curved Bezier patches of degrees 1
to 3 (polynomial and rational), a patch with interior knots, a rational quadratic
cylinder that is smoothly periodic in U and its transpose periodic in V, and
patches with a collapsed edge (poles, in U and in V). Distances are 0, 0.5 and
-0.5. Parameters are outside the domain by a third of its range and by more, on
one side and on two, and inside it. The sequences evaluate outside, come back
inside, go further outside (so SMLib extends its extension), and use each
combination of evaluation orders and flags. The same steps run with extension
off, where the parameter is clamped, and as seeded random sequences.

Every binary64 output is bit-identical; the gate accepts no unit-in-the-last-place
budget. 269 of the 1,080 cases include a failing step, which the port matches as
a status, and 712 hold an extension at some step. The zero-offset views that the
fillet code evaluates share the offset surface's extension, as SMLib's one object
does after `SetSurfacesToZeroOffset`.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 -fno-access-control \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/offset-extend/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/offset-extend/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
