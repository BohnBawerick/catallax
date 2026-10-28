## 2026-10-26: regular partial-cylinder and pole-free sphere STEP charts

Pieces: 7, STEP consumers, CLI verification

Change: Regular partial circular-cylinder and pole-free spherical faces now use
bounded parameter charts. Exact meridian cuts split shared physical edges once
and propagate those splits to every incident face. Source loop-corner incidence
is preserved at touching vertices. Disconnected chart components become separate
STEP faces with distinct surface identities; a cut shared by two components is
not mislabeled as a same-face seam.

Selection first tries exterior meridians so contractible partial walls and caps
do not acquire unnecessary face fragments. Spherical membership aims rational
query planes toward enclosed trim chords instead of relying only on blind
carrier rays that can all miss a small cap. The enclosures choose candidates,
not geometric approximations: arc membership and crossing orientation still
require exact proofs, and the existing request budgets are unchanged.

Each incident chart has a bounded degree-seven pcurve using the written space
curve's native parameter. Interval Taylor remainders, endpoint replacement and
decimal rounding contribute to the existing certified distance bound. The
surface frame is the frame actually written to STEP, not an idealized source
frame. Branch fitting reserves part of the requested distance budget for carrier
projection and pcurve fitting; this is not a tighter mass-property requirement.
The source carrier and kernel topology are unchanged. No geometry regularization
is introduced, and existing regularization records and displacement bounds are
preserved through the export-only face/edge remapping.

A vertex on a chart cut writes the cut's exact parameter, so every pcurve along
one side of the cut shares one `u` value. The pinned raw translator pairs a
seam's two pcurves by the exact direction of each pcurve's end-to-end chord.
Inverting rounded vertex points could leave `u` values one binary64 step apart,
and the translator then paired a cylinder band's seam the wrong way round. The
endpoint replacement term in the certified bound covers the change.

`StepReport.faces` counts emitted chart fragments. `added_vertices` includes
chart-cut vertices; `added_edges` includes the net increase from splitting
existing edges as well as new seams. `approximated_curves` counts emitted 3D
branch fits, so splitting a branch can increase it. Report fields and manifest
formats do not change. Complete-band and loopless exports retain their previous bytes.

An exact analytic 3D `CIRCLE` remains exact beneath its `SURFACE_CURVE`
association. Its auxiliary pcurve can still contribute a nonzero certified
`max_deviation`. Consumers compare the actual writer-returned rational bound
with the requested linear tolerance rather than assuming an analytic master
curve implies a zero bound; the 3D curve is not replaced by an approximation.

Eligibility is whole-export, not per-face. Only regular cylinders and proven
pole-free spheres are rebuilt as charts. Uncharted planar caps must retain
native analytic space boundaries: explicit lines/circles/ellipses or intersections
of a plane with planes, cylinders or spheres. Opaque boundaries on other faces
reject the entire cutover rather than inheriting a different branch fit.
Cone faces, arbitrary torus trims, elliptic cylinders and spherical domains
containing a parameter pole keep the whole scene on its legacy encoding.
The earlier complete-band path remains unchanged; conical bands no longer enter
the partial-chart path alongside partial faces. Structural rejection happens
before chart membership, cuts or fitting and uses no runtime import or measurement.
Before cutover, distinct endpoints of
every nonclosed written edge, charted or not, must remain separated after the
existing decimal rounding by more than `max(abs(curve_tolerance), 1e-7)` plus
their certified rounding errors. Each written vertex's rounding error must also
meet the requested distance tolerance; otherwise an existing over-tolerance
low-digit answer stays on its legacy path. The `1e-7` floor matches the pinned
raw translator's model-unit precision. A spherical chart face whose wires pass
the same vertex on the same chart side twice, as trims touching at a tangency
do, also keeps the legacy encoding; touching cylinder trims keep their charts.
Unresolved scenes keep their original bytes; this guard neither heals them nor
turns an existing answer into a new refusal.

Pole membership is a bounded search over chord-aimed and 32 blind query planes.
A pole without a certificate keeps the legacy encoding instead of refusing.
Charts are an optional encoding: any refusal raised while preparing or writing
them, other than an exhausted request allowance or deadline, falls back to the
legacy encoder.

Known follow-up: Public fixture `boolean-13` intersects a radius-two cylinder
with a box at a height offset of `1e-30`. Its certified volume is approximately
`2.666666666666667e-44`, while the existing seventeen-digit placement rounds that
offset away. It remains byte-identical and outside this repair. Its parent raw
import is invalid; the default standard reader reports valid topology but a
converged volume near `125.66370614359145`, not the source volume. Correcting
that rounding/regularization problem requires separate work.

Verification: The pinned OpenCASCADE regression retains the non-repairing raw
translator and also reports the standard reader separately. Partial-trim checks
use adaptive Gauss–Kronrod volume convergence rather than the default fixed
quadrature, which can disagree on fitted pcurves. Diagnostic volume agreement
allows the exported solid's total face area times its requested linear curve
tolerance, plus the spread of the converged integration sequence, around the
certified source interval. This first-order tube comparison is not a proof near
corners or contacts, not a certification and not a kernel tolerance. Exact
residuals, intervals, areas, tolerances, spreads and margins remain qualification evidence. The
public cylinder/sphere regression exercises the unchanged 100,000,000-work
request allowance. Qualification must retain byte identity for already
raw-valid fixtures and record every supported invalid-to-valid change.
