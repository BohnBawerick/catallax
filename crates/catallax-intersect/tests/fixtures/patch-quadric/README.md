# Patch-cylinder and patch-sphere oracle fixtures

These are untrusted parameter-curve controls from NVIDIA-Omniverse/usd-brep,
commit `c9c5979f04cea8c38213004d4f04df0dcfb40d05`, built with the same flags
as the [plane fixtures](../patch-plane/README.md). That build reproduces the
plane fixtures byte for byte. Catallax never compiles or links this C++ code.
Cargo tests read the saved text only.

`oracle.cpp` cuts the rational loft-like patch with the cylinder
`x² + z² = 25/4` and the sphere `x² + y² + z² = 25/4`. It cuts the bicubic
hill with a radius-1 cylinder along `x` through `y = 5/2`, and a radius-1
sphere centred over `(5/2, 5/2)`. Their lowest points are at `0.3`
(`hill`), `225/256-1e-6` (`tiny`) and `225/256` (`tangent`). `lofttight` asks
the tracer for `1e-5` accuracy; `loft` keeps its default.

To regenerate, build the pinned SMLib outside the repository as the plane
notes describe, then compile this driver against it:

```sh
c++ -std=c++17 -O1 -ffp-contract=off -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-intersect/tests/fixtures/patch-quadric/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-intersect/tests/fixtures/patch-quadric/generate.py "$ORACLE"
```

The SMLib library was built from all of `source/SMLib/src` with
`-O2 -fno-fast-math -ffp-contract=off` and the same preprocessor definitions.

`cylinder-tiny.txt` and `sphere-tiny.txt` contain no curve, although each
exact section has a closed loop. The two `tangent` files contain no curve or
contact point. The tight sphere trace starts `4.3e-19` off the `u = 0`
border. These defects are regression inputs. The complete Catallax search uses
none of these files to choose what geometry exists. Raw evidence retained
privately.
