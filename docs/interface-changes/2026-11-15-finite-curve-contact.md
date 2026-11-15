# Finite curve contact and numerical spline-line queries

The recorded quadratic loft now builds at a finite tangent interval and at
placements 1 micrometre above and below it. The tangent and above cases retain
material. The below case removes the thin intersecting region. The finite
queries retain closed knot domains instead of extending a straight spline
into an unbounded line.

The cheap exact path reuses a univariate patch section when its active atlas
chart has finite physical ends. Complete point sets clip these sections
without changing the existing two-sided contact labels. Straight spline
comparison proves the complete finite interval, including opposed senses.
Linear and quadratic algebraic cuts use the existing exact root solver.
Spline-line comparison needs only the bounded spline domain and a valid
unbounded line to prove different complete point sets. Its point query checks
both transverse constraints and filters every root to its finite span.
Partial overlaps still require interval ownership and refuse as isolated points.

`numerical_spline_line_points` handles higher-degree algebraic coefficients
through a separate numerical result. Its preset tolerance is 1e-7 mm. It
examines closed knot-span intervals using 128-bit position and derivative
enclosures, rejects separated intervals, subdivides the remaining intervals
and refines candidates with clamped Newton steps. Each returned rational point
is checked against the retained spline and line. The result records both
parameters, residual samples and the number of examined cells. Overlap or
unresolved cells refuse. Adjacent candidate intervals merge at the numerical
resolution; this is not an exact root count or a topology certificate.

Numerical points do not enter the exact Boolean meeting adapter. The return
type carries `NumericalEvidence`, has no exact-answer conversion and publishes
no incidence facts. Numerical body construction and approximate STEP still
need their own consumers and checks. Cheap exact cases keep their exact status.

The numerical approach follows the finite domains, subdivision, Newton
refinement and distance checks in OCCT's
[`IntCurve_IntPolyPolyGen` at V7_9_0](https://github.com/Open-Cascade-SAS/OCCT/blob/V7_9_0/src/IntCurve/IntCurve_IntPolyPolyGen.gxx).
Catallax uses its existing interval evaluator rather than copying OCCT code.
The [P04 port map](../phase1/P04/PORT-CONSUMERS.md) places candidate search and
curve services in separate units and identifies the limits of candidate parity.
The new query does not infer completeness from a tracer's empty output.

Focused regressions cover three isolated components, algebraic line
coefficients, reversed operands, finite endpoints, tangency, two-sided off-line
controls and partial overlap. `loft_curve_contact` checks the actual recipe
Boolean and material witnesses under the unchanged 20M operation allowance.
These checks do not complete milestones 1.1 or 1.3.
