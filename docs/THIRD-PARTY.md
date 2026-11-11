# Third-party source code

Catallax's original code remains MIT licensed. The following adapted files are Apache-2.0 licensed and retain their NVIDIA copyright and modification notices. The affected crate manifests use `MIT AND Apache-2.0`.

| Catallax file | Donor unit | Adaptation |
| --- | --- | --- |
| `crates/catallax-shapes/src/spline_tools.rs` | rust-brep `crates/kernel/src/curve/construction.rs`, `N_CrvExtractCrvSeg` | Exact `Rational` weighted-coordinate endpoint insertion and tensor restriction. Separate Euclidean controls and positive weights, finite-domain checks, typed refusals and caller work checkpoints |
| `crates/catallax-mesh/src/preview/imported.rs` | rust-brep `crates/kernel/src/tess_subdivide_mesh.rs`, `ST_TessSubdivideSrf` rectangle stack and midpoint bisection | Immutable source knots, Catallax interval evaluation, sampled flatness, bounded work and memory, shared STEP edges and unproved `PreviewMesh` output |
| `crates/catallax-trace/src/scalar.rs` | SMLib `SmMath`, `SmTol::GetScaledZero`; rust-brep `vector2d.rs`, `sm_geom_utility.rs` | Constants and clamped scalar functions only, without SMLib's error state; inverse sine and cosine from `libm` |
| `crates/catallax-trace/src/vector.rs` | rust-brep `crates/kernel/src/vector2d.rs` and `vector3d.rs` (`SmVector2d`, `SmVector3d`) | The operations the intersection chain enters, `Fault` results and named cross products |
| `crates/catallax-trace/src/extent.rs` | rust-brep `crates/kernel/src/extent1d.rs`, `extent2d.rs`, `extent3d.rs`; SMLib `SmExtent2d::ClipLine2d` | The interval and box operations the intersection chain enters, `Fault` results and returned optional outputs |
| `crates/catallax-trace/src/geometry.rs` | rust-brep `crates/kernel/src/sm_geom_utility.rs`, `sm_geom_intersect.rs`; SMLib `smgu_EvaluateCurveContinuity` | Returned values and `Fault`s instead of status codes and output references |
| `crates/catallax-trace/src/solution.rs` | SMLib `SmSolutionArray.h`, `SmSolutionArray.cpp` (`SmSolutionEnd`, `SmSolution`, `SmSolutionArray`); rust-brep `crates/kernel/src/solution_array.rs` | Object and tree-node pointers become identities, object kinds and parameter domains come from a caller-supplied table in place of a fixed object registry, and reads through a null object return `Fault`s |
| `crates/catallax-trace/src/global_solver.rs` | SMLib `SmGlobalSolver.h`, `SmGlobalSolver.cpp`; `SmCoreTypes.h` (operation, category, request and sort enumerations); `SmTree.h` (the node fields the solver reads) | Translated directly, rust-brep has no port. Virtual methods become the `SolverHooks` trait with SMLib's behaviour as the default, trees become an index arena, the traversal keeps its own stack in SMLib's order, and the curve-on-surface branch returns a typed refusal until its step ports |
| `crates/catallax-trace/src/search_math.rs` | SMLib `SmExtent3d` (box distances) and `SmVector3d` (rotated and perspective projection of a point) | Translated directly; assertion macros dropped |
| `crates/catallax-trace/src/spline.rs` | SMLib `SmBSplineCurve.cpp`, `SmBSplineSurface.cpp`, `SmNurbsSrf.cpp`, `SmNurbsCrv.cpp`, `NL_Math.cpp`, `NL_SrfTool.cpp`; rust-brep `curve/sm_nurbs.rs`, `surface/sm_bspline_net.rs`, `surface/partial_net_decompose.rs`, `basis/basic.rs` | Owned candidate curves and control nets, checked storage, object services, evaluations and splits; no contexts, inheritance or allocation emulation. [Differential gate and supported operations](../crates/catallax-trace/tests/fixtures/spline/README.md) |

The donor units are pinned to rust-brep commit `4c83df59d777ccfcb14f441b1c8d984d1b5d93ab`, translated from NVIDIA USD BRep SMLib commit `c9c5979f04cea8c38213004d4f04df0dcfb40d05`. Units rust-brep has not ported are translated directly from that SMLib commit. See the repository [NOTICE](../NOTICE) and the preserved [Apache-2.0 licence](../licenses/rust-brep-APACHE-2.0.txt). NVIDIA does not affiliate with or endorse Catallax or rust-brep. No trademark rights are granted.

`dist-workspace.toml` explicitly includes `NOTICE` and `licenses/` in platform archives. The universal macOS assembler retains those assets, and the release smoke checks their bytes against the source checkout before exercising the packaged programs.

No donor topology, error-state machinery, status codes, oracle implementation, C++ runtime emulation, TBB or fitting routines are included. The adapters use existing Catallax dependencies. The donor's floating-point differential records do not qualify these changed exact-arithmetic algorithms or certify preview triangles.

Surface fitting was deliberately deferred. Current Catallax spline families already have checked carriers and writers. Donor fitting becomes useful when a new definition-held family, such as a loft or general sweep, needs a carrier and its boundary and whole-domain proof contract exists. Exact surface definitions and the existing certified STEP intersection-curve fitter are unchanged.

The [M10 plane fixtures](../crates/catallax-intersect/tests/fixtures/patch-plane/README.md)
and [cylinder and sphere fixtures](../crates/catallax-intersect/tests/fixtures/patch-quadric/README.md)
are parameter-curve control data generated by the pinned C++ SMLib above.
A separate test-only driver and regeneration script record their provenance.
Cargo reads the saved data and never builds or links SMLib. The complete cell
search and tube certificate are Catallax code; the tracer supplies no proof.

`catallax-trace` holds the spline-intersection port that M9 extracts straight
into Catallax under the owner's 11 November ruling, in the M25 unit order. It
produces unproved floating-point candidates only; a candidate becomes part of
an answer only through Catallax's own proofs. Its [foundation gate
fixtures](../crates/catallax-trace/tests/fixtures/foundation/README.md) are
outputs of the pinned C++ SMLib above, recorded by a test-only driver, as are
the [search gate fixtures](../crates/catallax-trace/tests/fixtures/search/README.md)
for the general search solver and the solution array. Cargo reads the saved
data and never builds or links SMLib.
