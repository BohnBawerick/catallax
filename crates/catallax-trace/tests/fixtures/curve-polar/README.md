# Line and Hermite polar-box gate

The oracle uses NVIDIA-Omniverse/usd-brep SMLib at
`c9c5979f04cea8c38213004d4f04df0dcfb40d05` on x86-64 Linux. It calls the
polar part of `SmLine::CalculateBoundingBox` and
`SmHermiteCurve::CalculateBoundingBox` directly, with the binary64 protocol
and constructors of the [output-curve gate](../output-curves/README.md).
Catallax does not compile or link C++.

`generate.py` selects the 898 line and Hermite inputs from that gate's
seeded generator. The line cases cover bounded constructors, segments,
inside-out parameters, planar projection, edited parameter ranges, short
and collapsed intervals, and directions around `SM_EFF_ZERO`. A line box
uses the evaluated interval chord, including after reparameterization.
Hermite cases cover general and planar cubics, coincident ends, zero and
small end derivatives, and stationary points. Each Hermite runs over
`[0,1]`, `[0.2,0.7]` and `[0.5,0.5]`; SMLib uses all three control-polygon
edges regardless of that interval.

The gate compares the state exactly and the basis and angular domain within
16 ULP, matching the existing polar-box candidate policy for inverse
trigonometry. An uninitialized box has no defined C++ basis, so only its
state is read. Every run must include tangent, touching and singular cases.
The Rust test also splits the same fixture across 2 and 4 workers and
compares the candidate bits and order to the 1-worker result.

The measured gate has 23,415 of 23,439 binary64 outputs bit-identical.
The other 24 values are angular bounds, one ULP apart, in eight Hermite
cases repeated over three intervals. Candidate bits and order agree at
1, 2 and 4 workers.

These boxes remain unproved candidates. No kernel proof calls the new
methods, and the port does not promote their tolerance tests to exact
bounds. Any future proof consumer must establish its own enclosure before
using a polar box to reject an intersection. The exact acceptance routes,
tolerances and work allowances stay unchanged. Known C++ misses have a
separate [R1 exception list](../../../../../docs/M9-KNOWN-SMLIB-MISSES.md).

Build the pinned library as in the
[foundation instructions](../foundation/README.md), then compile the driver
outside the repository:

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/curve-polar/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/curve-polar/generate.py "$ORACLE" "$RAW_DIR"
cargo test -p catallax-trace --test curve_polar -- --nocapture
```

Raw evidence retained privately.
