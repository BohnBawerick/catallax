# Numerical STEP for native spline boundaries

`catallax_step::numerical::build` accepts checked numerical bodies with finite line or native B-spline edges and planar or single rational-patch faces. It serializes every retained physical loop and shared edge/vertex identity. Patch edges carry finite linear UV pcurves in the edge's parameter direction. The original planar line path keeps its existing behavior.

The curved path reconstructs the decimal controls, knots, weights, finite parameters, line magnitudes and STEP plane frames, then admits that model under the body's numerical tolerance. The report keeps this emitted evidence separate from the generating definition, source evidence and edit history. Planar spline walks retain their sampled check names. Patch trim polygons retain their declared linear UV interpretation. Vertex positions and topology are retained without healing.

`NumericalStepReport::spline_rounding` is `Some(NumericalSplineRounding)` for this path and `None` for the existing planar line path. It records the maximum parameter error in native parameter units and the maximum dimensionless weight error. These values are separate from the length-scalar error and normalized direction-ratio error. Unchecked native curve closure and self-intersection flags are written as unknown.

The artifact and STEP description remain numerical and approximate. The tolerance governs the recorded checks; sample maxima and scalar rounding errors do not establish a whole-body deviation bound. Internal file checks cover byte identity and the seam representation contract. Raw OpenCASCADE import and edit acceptance require the separate reference check on the actual file.
