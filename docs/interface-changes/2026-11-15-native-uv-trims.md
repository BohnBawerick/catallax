# Retained native UV trims

`numerical::build_with_native_trims` admits a `NumericalDefinition` with optional `NumericalNativeTrim` sidecars.
A binding names a stable face and edge plus its loop/use address.
Its B-spline x/y coordinates are surface parameters at the 3D edge's native parameter; z is unused.
The complete finite edge interval must lie in the UV carrier's domain.
Duplicate, unknown and stale bindings refuse.
Existing definition and trim literal layouts stay unchanged, and empty sidecars preserve ordinary numerical builds.

The producer follows SMLib's `SmCrvOnSurf` composition in `catallax-trace`: evaluate the UV curve at t, then evaluate the surface there.
It merges the 3D/UV knot boundaries clipped to the finite edge interval and checks endpoints and eight subdivisions per resulting span.
Each actual UV sample must lie in the retained patch domain without clamping.
The checks record 3D-edge/surface distance and physical endpoint displacement from the declared closed-loop endpoints, in millimetres.
They do not treat parameter-coordinate differences as millimetres.
Every join on a face with native bindings compares the two actual pcurve endpoints in physical space, including the closing join.
These checks run on the source and again on the actual decimal definitions, before declared loop endpoints can hide a gap.
The sampled loop uses declared endpoints and native interior locations; simplicity and hole placement remain qualified as sampled.
Existing analytic periodic-seam and continuous circular-hole clearance checks still run before native loop sampling.
Continuous containment of general native curves, separation and whole-curve deviation remain unresolved.
All results retain their numerical category and stated tolerance, including zero sampled residuals.

`NumericalBody` retains original and current native trims.
Rebuild, surface replacement and vertex repair recheck them and preserve earlier bindings in `NumericalHistory`.
`rebuild_with_native_trims` explicitly replaces the current bindings when an edit changes their addresses or definitions.
The generating source and previous checked body remain unchanged on refusal.
The body handover codec retains original, current and historical bindings and their native curve graphs.
Numerical body and history records start with separate layout tags, `CXNUMB02` and `CXNUMH02`.
Earlier positional payloads refuse with `handover.version-mismatch`; they are never interpreted as the new layout.

The numerical STEP writer emits retained rational 2D B-splines through PCURVE on rational patch and plane faces.
Analytic exchange keeps its existing empty-sidecar path; a native binding on a currently unsupported analytic carrier refuses instead of dropping the binding.
The 2D and 3D curves use the same finite native parameter interval and traversal sense.
Before decimal rounding, the writer restricts each UV basis exactly to that interval with positive homogeneous blossom extraction, retaining every clipped span and its original parameter values.
This follows OCCT's [native 2D B-spline representation](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/GeomToStep/GeomToStep_MakeBSplineCurveWithKnots.cxx#L47).
The finite basis is necessary because [OCCT's 2D trimmed-curve reader returns a B-spline basis before reading the trim bounds](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/StepToGeom/StepToGeom.cxx#L2229).
A TRIMMED_CURVE wrapper alone does not preserve that finite basis through this reader.

UV controls and knots use native parameter rounding; weights have their own rounding record.
The actual decimal sidecars are re-admitted with the decimal 3D geometry.
`NumericalStepReport` retains the original, current source and emitted native trims and their separate evidence.
Exact restriction of retained definitions adds no exact relation to an original analytic source and no STEP approximation bound.

The `numerical_loft_cut` example has a `native-uv` control on M10-008.
It reparameterizes one existing generatrix with `g(t)=t/2+3*t²/2-t³`, keeping the physical segment unchanged.
The 3D carrier and both UV carriers cover `[-1,2]`, while the edge uses `[0,1]`.
Missing, reversed, partial, wrong, stale and duplicate controls must refuse without stale STEP output.
An optional final argument selects 1, 2 or 4 workers.
The raw STEP checker accepts `--native-uv-root` to check the physical hole, material, topology, volume, native parameter correspondence and a later OpenCASCADE cut with healing off.
