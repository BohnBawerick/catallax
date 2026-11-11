# Ruled-surface gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmBSplineSurface::CreateRuledSurface`
(through `N_CreateRuledSrf` and `N_CrvsMakeCompatibleKnotTol`) and
`SmBSplineSurface::Reparameterize` (through `N_SrfReparamToInterval`). The
library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles
or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. A surface is printed as stored: degrees, control-net and
knot counts, every knot of both directions and each control point's four stored
coordinates in the order SMLib stores the net. This gate accepts no
unit-in-the-last-place budget: every binary64 output is bit-identical.

The ruled cases join curves of degrees 1 to 3 in each combination, polynomial
and rational, planar and spatial (including a planar curve against a spatial
one), in both directions. The curves have different interior knots (so the
merge refines one or both), identical knots, knots 1e-10 and 1e-7 apart,
different parameter ranges (scaled to the unit interval), a range differing by
1e-13 and by 1e-9, and control points that all coincide. The reparameterisation
cases take surfaces of degrees 1 to 3 in each direction, polynomial and
rational, to new ranges, to the range they already have (left alone), to ranges
1e-13 and 1e-9 off, to reversed ranges and to a range of width 1e-9.

One group is not compared value for value. SMLib reparameterises a surface to a
single-valued range and returns a surface whose active knot interval is empty;
a checked `BSplineSurface` has a positive-length domain, so the port refuses
with `InvalidInput`. The 30 `reparam-single` cases are pinned by their own test
in `../../surface_build.rs`, which checks the port's refusal against SMLib's
success.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/surface-build/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/surface-build/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
