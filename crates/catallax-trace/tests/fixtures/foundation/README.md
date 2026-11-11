# Foundation gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the vector, extent and
`SmGeomUtility` routines that `catallax-trace` adapts. The library was built
on x86-64 Linux from every `source/SMLib/src/*.cpp` with the rust-brep oracle
preprocessor definitions, `-O2 -fno-fast-math -ffp-contract=off` and system TBB.
Catallax never compiles or links this C++ code. Cargo reads the saved text only.

Each line is `operation kind | inputs | outputs`. The kinds are `regular`,
`tangent`, `touching`, `singular`, `near-singular` and `boundary`; the test
refuses to pass without tangent, touching and singular cases. The comparison
rules are in `../../common/mod.rs`: integers, statuses and classes match
exactly, and binary64 values match within 4 units in the last place, or 16 for
routines that call the inverse sine or cosine, which the port takes from `libm`
and the oracle from the C library.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/foundation/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/foundation/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
