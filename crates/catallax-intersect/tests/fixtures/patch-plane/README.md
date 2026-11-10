# Patch-plane oracle fixtures

These are untrusted parameter-curve controls from NVIDIA-Omniverse/usd-brep,
commit `c9c5979f04cea8c38213004d4f04df0dcfb40d05`, using the rust-brep oracle
preprocessor definitions with floating-point contraction disabled. Catallax never compiles or links this C++
code. Cargo tests read the saved text only.

`oracle.cpp` cuts the rational loft-like patch with `x+y+z=5/2`, the wavy patch
with `z=0.35`, and the hill with `z=0.3`, `225/256-1e-6`, and `225/256`.
`lofttight` asks the tracer for `1e-5` accuracy. `loft` keeps its default.
The driver checks that each parameter curve is two-dimensional and
nonrational before emitting unit weights. The unused SMLib weight sentinel is
never treated as a curve weight.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against that library. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` names an
output executable outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-intersect/tests/fixtures/patch-plane/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-intersect/tests/fixtures/patch-plane/generate.py "$ORACLE"
```

Each `POINTS` row stores exact round-trippable binary64 `(u,v,weight)` triples.
The test reader checks the cubic Bezier knot decomposition and retains every
binary64 coordinate as an exact rational. It accepts no SMLib accuracy claim.
`wavy-heights.txt` freezes the surface heights from the same driver.

`tiny.txt` contains no curve, although the exact section has a closed loop.
`tangent.txt` contains no curve or contact point. These omissions are part of
the regression inputs. The complete Catallax search uses neither file to
choose what geometry exists. Raw evidence retained privately.
