# Numerical native intersections

`catallax_intersect::numerical_native_windows` builds finite cone/cylinder windows in either operand order. It retains the original analytic carriers, rounded native B-spline sources, positive cone-height interval, cylinder interval and actual source/control checks. Algebraic axes use 128-bit enclosures. The preset sample tolerance is 1e-7 mm.

`numerical_native_intersection` searches those windows and retains every returned 3D fit and both UV fits. All three definitions keep the same native parameter interval. Each nonzero span receives endpoint and quarter-span checks. UV sample locations must remain inside the finite source domains, without clamping. Independent checks record distances to the native patches and original analytic carriers.

The result is `NumericalNativeIntersection`, with immutable pieces, `NumericalEvidence`, tracer visit counts and measured work stages. It publishes no exact incidence, complete branch-count proof, whole-curve bound, Boolean topology or STEP result. Exact pair plans retain their existing proof requirements. Windows reaching a parallel cone generator or apex still need finite-face or singular-point handling.

The search follows the P04 SMLib `SmAdvSurfaceIntersector` consumer through finite native patches, subdivision, local refinement and piecewise cubic output. Checked output uses each span's interior-side UV derivatives at native circle knots. Ordinary unproved candidate consumers retain their existing output. OCCT also keeps separate walking and curve-approximation controls in [`GeomInt_IntSS::MakeCurve`](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/GeomInt/GeomInt_IntSS_1.cxx#L241). This producer uses that distinction to fit at one quarter of its preset sample tolerance and then check the retained definitions independently.

`Intersector::do_intersection_with_work` enables thread-local checkpoints only during its call. Searches, local iterations, surface/curve evaluations, source items and trace/output visits call the budget callback. Cancellation discards partial output even when a ported helper catches a failure. `WorkReport` counts checkpoint visits. The numerical query's work stages include the callback's additional deadline checks and exact validation work. Standalone queries use 20M; queries inside an operation share its allowance.

The portable example reads actual recipes ending in one hole and checks both cone faces in both operand orders:

```sh
cargo run --release -p catallax-cli --example numerical_hole_carriers -- \
  /tmp/native-pairs-1 1 E07-axial.json E07-skew.json
```

Repeat with fresh output directories and worker counts 2 and 4. Reports retain the native definitions, every sampled location, evidence and work stages. They are numerical query diagnostics and cannot be reused as exact manifests. The regression tests cover the native-knot derivative defect, both operand orders, off-carrier controls, cancellation and unchanged ordinary candidate output.
