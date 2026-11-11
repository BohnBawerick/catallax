# Local solver gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmLocalSolve1d::SolveIt` and
`SmLocalSolveNd::SolveIt` on the test functions of `oracle.cpp`, which
`../../local_solve.rs` repeats operation for operation. The library build is
the one described in [`../foundation/README.md`](../foundation/README.md).
Catallax never compiles or links this C++ code.

The one-variable functions are a cubic with simple roots, a parabola with a
double root or a near miss, a signed square root that makes Newton steps
oscillate, a line that cannot be evaluated past a point, and a double well
solved with each of SMLib's algorithms. The several-variable functions are a
circle with a crossing, tangent or missing line (with an exact or a
finite-difference Jacobian), two circles that cross, touch or separate, and a
rank-deficient linear system. Cases also cover parameter boxes, periodic
wrapping, the three boundary handlers and step cutting. The comparison rules are
in `../../common/mod.rs`; no inverse trigonometry is involved, so the budget is
4 units in the last place.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/local-solve/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/local-solve/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
