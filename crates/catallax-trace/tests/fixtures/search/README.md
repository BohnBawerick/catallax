# Search gate fixtures

`cases.txt` holds the outputs of NVIDIA-Omniverse/usd-brep SMLib, commit
`c9c5979f04cea8c38213004d4f04df0dcfb40d05`, for the routines of port step 2
that `catallax-trace` adapts: the general search solver (`SmGlobalSolver`)
and the solution types and array (`SmSolutionEnd`, `SmSolution`,
`SmSolutionArray`). The library was built on x86-64 Linux from every
`source/SMLib/src/*.cpp` with the rust-brep oracle preprocessor definitions,
`-O2 -fno-fast-math -ffp-contract=off` and system TBB. Catallax never compiles
or links this C++ code. Cargo reads the saved text only.

Each line is `operation kind | inputs | outputs`, in the format of the
foundation gate (`../foundation/README.md`). The comparison rules are in
`../../common/mod.rs`: integers, statuses and classes match exactly, and
binary64 values match within 4 units in the last place. Every binary64 output
of the saved run is bit-identical.

The kinds are `regular`, `tangent`, `touching`, `singular`, `near-singular`
and `boundary`; the test refuses to pass without tangent, touching and
singular cases. Here the kinds mean: `tangent` is boxes that share a face or
ranges that share an end; `touching` is a gap or a parameter difference at
the tolerance an SMLib test uses; `boundary` is a value one unit in the last
place either side of a comparison, or exactly on it; `singular` is a
zero-sized box, a zero vector, a solution without objects, or a tree with no
node to split.

| Operation | What it runs |
| --- | --- |
| `solcmp` | `SmSolution::operator==`, `ClassifySolutionEnd` for the start and the end, `IsEqual` and `SmGlobalSolver::IdenticalSolutions` on pairs of points and ranges |
| `solidx` | `GetIndex`, `GetObjectIndexForParameterIndex`, `GetParameterInterval` and `GetObjectTolerance` for object mixes of curves, surfaces, a vertex, an edge, a face and an unclassified object |
| `addsorted` | `SmGlobalSolver::AddSortedSolution` for every sort key, request and operation category: ordering, duplicates, merging and splitting of ranges, best-answer tracking and trimming |
| `array` | `SmSolutionArray` `Add`, `InsertAt`, `RemoveAt`, `SetSize`, `SetAt`, `SetAtGrow`, `Swap`, `RemoveDuplicateSolutions`, `AddSortedSolution`, `Append`, `Copy` and `ReSet` |
| `branch` | `BranchMayContainAnswers` for every operation SMLib handles, with boxes apart, touching, overlapping and degenerate |
| `tree` | `AreReadyForLocalSolve`, `FindBestSubdivisionNode` and `Subdivide` on random trees, equal-size ties, flat boxes and Bezier-surface markers |
| `solve` | `SolveIt` with SMLib's own `LocalSolve`, which records the leaf node pairs, over one to three trees |

## What the gate does not run

- The `SM_SO_SURFACE_PROJECTED_INTERSECT` branch of `BranchMayContainAnswers`
  that samples a curve on a surface needs `SmCrvOnSurf`, which a later step
  ports. The port returns `SolverRefusal::NotPorted` there. Only the
  overlapping-box answer, which SMLib gives before that code, is gated.
- Edge, face and volume objects. The oracle builds an edge and a face; it has
  no volume. `tests/solver_behaviour.rs` covers the volume rules, and SMLib's
  null-pointer reads (a missing second object, a ray-fire order on a leaf), which
  return typed refusals here.
- `SmSolution::GetPoint`, which evaluates geometry that later steps port.
- Derived solvers. SMLib's derived solvers override `LocalSolve` and a few
  other virtual methods; the port exposes them as the `SolverHooks` trait.
  Step 2 gates only the base behaviour that every derived solver inherits.

`BranchMayContainAnswers` for `SM_SO_PERSPECTIVE_INTERSECT` projects the second
sphere's offsets from the first sphere's centre and sizes both spheres with the
first sphere's radius. The port keeps both, so the gate's perspective cases
exclude a branch in 4 of 41 runs.

## Regenerate

Build the pinned SMLib outside the repository, then compile the driver
against it. `SMLIB_ROOT` names its `source/SMLib` directory; `SMLIB_LIB_DIR`
names the built shared library directory. `ORACLE` and `RAW_DIR` name an
executable and a folder outside the repository.

```sh
c++ -std=c++17 -O1 -ffp-contract=off -w -DNDEBUG -DNOMINMAX \
  -DSM_GRAPHICS_CALLBACKS=1 -D_GLIBCXX_USE_CXX11_ABI=1 \
  -I"$SMLIB_ROOT/inc" crates/catallax-trace/tests/fixtures/search/oracle.cpp \
  -L"$SMLIB_LIB_DIR" -Wl,-rpath,"$SMLIB_LIB_DIR" -lSMLib -ltbb -ltbbmalloc -o "$ORACLE"
python3 crates/catallax-trace/tests/fixtures/search/generate.py "$ORACLE" "$RAW_DIR"
```

The generator runs every case in a fresh process under a 1 GiB address-space
cap, a 10 s CPU cap and a 20 s timeout. Raw evidence retained privately.
