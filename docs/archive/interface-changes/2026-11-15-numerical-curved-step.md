Archived on 2026-11-26. Historical record only. Not current rules or plans. Current: [Interface changes](../../INTERFACE-CHANGES.md).

# Numerical STEP for native spline boundaries

`catallax_step::numerical::build` accepts checked numerical bodies with finite line or native B-spline edges and planar or single rational-patch faces. It serializes every retained physical loop and shared edge/vertex identity. Every native 3D spline carries an explicit `TRIMMED_CURVE` over its admitted parameter interval, with STEP parameter master representation. Repeated endpoint preimages outside that interval cannot select another traversal. Patch edges carry finite linear UV pcurves in the same native parameter direction. The writer indexes face/coedge uses in one charged walk, including planar uses it skips for pcurve output. The original planar line path keeps its existing behavior.

The curved path reconstructs the decimal controls, knots, weights, finite parameters, line magnitudes and STEP plane frames, then admits that model under the body's numerical tolerance. The report keeps this emitted evidence separate from the generating definition, source evidence and edit history. Planar spline walks retain their sampled check names. Patch trim polygons retain their declared linear UV interpretation. Shared topology stays unchanged. Vertex coordinates change only through recorded decimal rounding; the writer performs no healing.

`NumericalStepReport::spline_rounding` is `Some(NumericalSplineRounding)` for this path and `None` for the existing planar line path. It records the maximum parameter error in native parameter units and the maximum dimensionless weight error. These values are separate from the length-scalar error and normalized direction-ratio error. Unchecked native curve closure and self-intersection flags are written as unknown.

The artifact and STEP description remain numerical and approximate. The tolerance governs the recorded checks; sample maxima and scalar rounding errors do not establish a whole-body deviation bound. Internal file checks cover byte identity and the seam representation contract. Raw OpenCASCADE import and edit acceptance require the separate reference check on the actual file.

The `numerical_loft_cut` example builds the actual near-contact loft and recorded M10-008 through this API at an explicit 20M allowance and 1e-7 mm tolerance. It saves the generating recipe, complete source/emitted numerical checks and stage work beside the STEP. Its reversed- and missing-hole controls remove a prior output before admission, and publish a refusal response with no STEP.

```sh
cargo run --release -p catallax-step --example numerical_loft_cut -- /tmp/numerical-curved/tiny tiny
cargo run --release -p catallax-step --example numerical_loft_cut -- /tmp/numerical-curved/m10-008 m10-008
cargo run --release -p catallax-step --example numerical_loft_cut -- /tmp/numerical-curved/reverse-hole m10-008 reverse-hole
cargo run --release -p catallax-step --example numerical_loft_cut -- /tmp/numerical-curved/missing-hole m10-008 missing-hole
python3 tools/step-regressions/check_numerical_curved_step.py /tmp/numerical-curved
```

The two invalid controls exit 1. The raw checker requires OCP 8 and checks finite topology, physical void/material witnesses, volume, healing-off import and a later OCCT cut on both actual files. It also verifies both refused outputs contain no stale STEP. The main exact `catallax-step/1` request adapter retains its existing contract; this example exposes the separate numerical Rust path.

The repeated-endpoint regression uses a degree-one carrier A, B, C, A, B with C off both incident planes and an admitted interval [0,1]. Set `CATALLAX_N03_CAPTURE` while running the `numerical_curved` tests to save `repeated-endpoint.step`. Pass that file as `--repeated-endpoint` to the raw checker to verify the physical volume and finite edge samples. The pcurve index unit test checks exact visit charges and exhaustion one unit below them.
