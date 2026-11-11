# Surface point drop gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for `SmSurface::DropPoint` (its
success flag and, on success, the parameters, gap and multivalued flag) on
surfaces without a face. The library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code.

The surfaces are those of the surface point solver gate (the bicubic bumps,
the loft, the patch with a pole, the bilinear patch and the cylinder) and two
planes, one unbounded and one bounded with scaled parameters. The cases drop
points on, above and beyond the surfaces, over whole domains and sub-domains,
without a guess (`DropPointFast` on planes, then the global search) and from
good, poor and far guesses (the local solve, kept or refused by
`sm_GoodDropResult`), minimising, normalising and meeting within SMLib's
default intersection tolerance of 2e-5. They include the pole and points
turned towards it, the cylinder's axis and seam (where the search returns two
seam answers), points on a plane, 1e-6 and 1e-4 off it and beyond a plane's
sub-domain, and a seeded sweep of 75 points with guesses on and near the sides
of sub-domains.

Binary64 values match within 16 units in the last place; every one matched bit
for bit when the fixture was made. Removing the zero-gap acceptance, halving
the step limit, removing the small-gradient acceptance or the comparison of
the guess's and the result's distances, or the multivalued flag, changes
cases. The side tests of `sm_GoodDropResult` (a result on a side whose
gradient points inward, an interior result with a gradient above 0.001) only
act when the local solver stops short, which no case reaches. Picking the
answer nearest the guess among several global answers only acts when the local
solve from that guess fails, which no case reaches either. A plane's
`DropPointFast` clamps into the domain, so its containment test always passes.

To regenerate, build the driver as in the foundation README with
`crates/catallax-trace/tests/fixtures/surface-drop/oracle.cpp`, then run

```sh
python3 crates/catallax-trace/tests/fixtures/surface-drop/generate.py "$ORACLE" "$RAW_DIR"
```

Raw evidence retained privately.
