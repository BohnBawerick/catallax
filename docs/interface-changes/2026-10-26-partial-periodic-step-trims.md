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

Eligibility: A selected scene must use supported regular charts throughout.
Arbitrary cone/torus trims, elliptic cylinders and spherical domains containing
a parameter pole keep their legacy encoding. Complete non-apex conical bands
can accompany supported partial faces. Before cutover, distinct endpoints of
nonclosed chart edges must remain separated after the existing decimal rounding
by more than `max(abs(curve_tolerance), 1e-7)` plus their certified rounding
errors. Each chart vertex's rounding error must also meet the requested distance
tolerance; otherwise an existing over-tolerance low-digit answer stays on its
legacy path. The `1e-7` floor matches the pinned raw translator's model-unit precision.
Unresolved scenes keep their original bytes; this guard neither heals them nor
turns an existing answer into a new refusal.

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
