# Surface-extension gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the NLib routines that
`SmOffsetSurface` extends a base surface with: `N_SrfRemoveKnotConditional`,
`N_SrfExtendByParamDist` (through `N_SrfEvalPtDerivsUnbounded`, for `NL_G1R`,
the continuity `sm_ExtendSurface` asks for) and `N_SrfIsClosedSmooth`. The
library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles
or links this C++ code.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. A surface is read as stored (`du dv nu nv rational
knots_u... knots_v... (x y z w)...`) and printed as stored, with its degrees,
net and knot counts. The operations are:

- `srf-remove`: a surface, a direction (0 for U, 1 for V), a parameter, a count
  and a tolerance. The output is NLib's error flag, the number of copies of the
  knot removed and the surface, edited in place as NLib's two surfaces are the
  same.
- `srf-extend`: a surface and a sequence of parameter pairs, evaluated with one
  `esur` carried through the calls. After each call the output is the error
  flag, the point when there was none, whether `esur` holds a surface, and then
  its closure flag and the extended surface. A parameter outside the domain
  extends the surface (or wraps, where it closes smoothly), and a later call
  extends the surface extended so far.

The surfaces are Bezier patches of degrees 1 to 3 in each direction, polynomial
and rational, with interior knots inserted by Boehm's algorithm so that they
can be removed exactly, perturbed by 1e-9 to 1e-1 so that they can only be
removed to a tolerance, with multiple knots, weights that a removal would take
out of range, a collapsed edge, and a rational quadratic cylinder that closes
smoothly in U. Targets are inside, on the edges and outside the domain, by less
and by more than the range, in both directions and in sequences that extend
twice.

Every binary64 output of the compared cases is bit-identical; the gate accepts
no unit-in-the-last-place budget. Three of NLib's habits are reproduced: the
Euclidean point of a control point is formed with the weight's reciprocal; a
removal that fails on a later row or column leaves the earlier ones edited; and
a basis maximum that cannot be found is reported as an error although the
removal went on with a maximum of 1.

**Pinned refusal, scratch overrun.** `N_SrfRemoveKnotConditional` keeps its
working rows in arrays that run to the index `2 * max(p, q)`, and removing a knot
of high multiplicity from a surface of degree 1 (20 cases) or 2 (5 cases) reaches
further. SMLib then writes past them and dies, so there is no answer to compare;
the 25 such cases are `srf-remove-overrun`, and the port refuses each. The test
checks that the port refuses these, and every other case is compared.

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/surface-extend/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/surface-extend/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
