# Conic gate fixtures

`cases.txt` holds outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmPolynomial::SolveCubicEqn`
and `SolveQuarticEqn` and for the planar conics of the Dupin indicatrix
(`SmConic`): points along a conic and their inversion, the implicit equation
of a moved conic, and `IntersectConic`. The library build is the one
described in [`../foundation/README.md`](../foundation/README.md). Catallax
never compiles or links this C++ code.

The polynomial cases have three, two, one and no real roots, double, triple
and quadruple roots, the quartic's `r = 0` branch, negligible and zero leading
coefficients, and seeded random coefficients. The conic cases pair ellipses
turned against each other, aligned, tangent and identical; an ellipse with
hyperbola branches; hyperbolas; lines (parallel, crossing, tangent); a
parabola; moved origins; and seeded random pairs.

The residual distance of each conic meeting is rounding noise that no caller
reads; the driver records only that it lies within the tolerance.

Binary64 values match within 16 units in the last place, except in the two
`conicx` cases below. When the fixture was made 773 of 793 values matched bit
for bit. The other differences, at most 8 units, come from the port taking
its cube root, sine, cosine and arc cosine from the `libm` crate (repeatable
last bits on every platform), which rounds some arguments one unit
differently from the C library that SMLib links; a double root of the cubic
magnifies one such unit to 8.

**Exceptions (`conicx`, 32 units).** Two pairs, an ellipse against a moved
circle and the second seeded random pair, differ by up to 25 and 29 units in
the meeting parameters (and 16 and 15 in the second parameter). The cause was
traced to the cube roots of the quartic's resolvent cubic: of the eight
`cbrt` calls the two pairs make, six return values one unit apart in `libm`
and in the C library (all four in the first pair, two in the second). With
the C library's cube root substituted, both pairs matched SMLib bit for bit.
The build has no fused multiply-add (SMLib compiled with
`-ffp-contract=off`; Rust does not fuse), so operation order and fusion are
not the cause. The quartic's near-coincident roots magnify the one-unit
differences. The test `conic_exception_cube_roots` in `tests/conic.rs` pins
the six arguments with both libraries' results. Matching the C library would
mean copying its cube-root code, which is under the LGPL, so the port keeps
`libm`. Whether to accept this exception or to adopt a correctly rounded,
permissively licensed cube root is an open decision.

Mutation checks: dropping the half-ellipse test of the implicit intersection,
the quartic's duplicate-root filter or its zero test of `u`, the hyperbola's
clamp or its second-root retry, or reversing a transform's angle changes
cases.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/conic/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/conic/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
