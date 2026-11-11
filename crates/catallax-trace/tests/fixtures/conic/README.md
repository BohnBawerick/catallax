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

Binary64 values match within 16 units in the last place. When the fixture was
made 723 of 737 values matched bit for bit and the rest were at most 8 units
apart. The port takes its cube root, sine, cosine and arc cosine from the
`libm` crate, kept for repeatable last bits on every platform; the C library
that SMLib links rounds some of them one unit differently, and a double root
of the cubic magnifies one such unit to 8.

Two pairs are left out of `cases.txt`, their outputs kept with the raw
evidence: an ellipse against a moved circle and the second seeded random
pair. There the resolvent cubic's cube root differs by one unit between the
two libraries and the quartic's roots magnify it to 25 and 29 units. Writing
out the C library's cube root would close that gap, but its code is under the
LGPL, so the port does not copy it.

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
