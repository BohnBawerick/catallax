# Curve-on-surface gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the part of `SmCrvOnSurf` that the
global solver reads (`EvaluatePoint`, `EvaluateUVPoint` and
`EvaluateSurfaceNormal`) and for the branch of
`SmGlobalSolver::BranchMayContainAnswers` that samples a curve on a surface
(`SM_SO_SURFACE_PROJECTED_INTERSECT`). The library build is the one described in
[`../foundation/README.md`](../foundation/README.md). Catallax never compiles or
links this C++ code. The driver is compiled with `-fno-access-control` to reach the
solver's protected branch test.

Each line is `operation kind | inputs | outputs`; the comparison rules are in
`../../common/mod.rs`. The inputs begin with a surface (`du dv nu nv rational
knots_u... knots_v... (x y z w)...`) and a parameter-space curve: `0` and eight
numbers for a Hermite curve (`p1 d1 p2 d2`), or `1` and a planar B-spline curve as
stored (`dim degree count knots... (x y z w)...`). The operations are:

- `cos-eval`: a parameter and a side; the output is the surface point, the surface
  parameters, and the parameter-curve position and tangent with the unit surface
  normal, each with SMLib's status.
- `cos-branch`: the node of the curve's cache tree that SMLib's tree walk
  numbers `index`, the box of a node of the other tree, and the tolerance; the
  output is whether the branch may contain answers. The driver builds the curve's
  cache, so the second node's box and span interval come from SMLib; the generator
  reads them with a `cos-node` call and appends them to the inputs, where the
  port reads them as the tree's node data (the port has no cache for curves on
  surfaces).

The surfaces are a bumpy and a rational bicubic patch, a rational loft, a flat
bilinear patch, a cylinder, a patch with a pole and a wave patch of many knots. The
parameter-space curves are Hermite and B-spline curves of degrees 1 to 3 that run
across each surface, one along the diagonal of a flat patch's domain and one along
the edge where a patch has its pole. Evaluations run at both ends and three interior parameters from each side.
Branches use the root, a middle and the last node of each cache tree, a box that
overlaps the node's, and boxes that touch it (tangent), are the tolerance away
(touching), 0.01, 0.2 and 3 away along x, z and the diagonal, on both sides, for
tolerances 1e-3 and 0.05.

Every binary64 output is bit-identical and every branch answer equal; the gate
accepts no unit-in-the-last-place budget. Mutating the normal's flip changes 330
cases and the sample count of ten changes two.

Not carried: the rest of `SmCrvOnSurf` (derivative lifting, bounding boxes, the
B-spline approximation and cache, dropping, trimming, transforms, database I/O).

To regenerate, build the pinned SMLib outside the repository, then compile the
driver against it. `SMLIB_ROOT` names its `source/SMLib` directory;
`SMLIB_LIB_DIR` names the built shared library directory. `ORACLE` and
`RAW_DIR` name an executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 -fno-access-control \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/curve-on-surface/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/curve-on-surface/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 2 GiB address-space cap, a
20 s CPU cap and a 40 s timeout. Raw evidence retained privately.
