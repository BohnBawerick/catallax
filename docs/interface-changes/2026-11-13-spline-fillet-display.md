# Spline fillet display and bounded rendering

Pieces: shapes, mesh, app

`mesh_body` prepares a checked cut chart for one spline fillet between two
polynomial bilinear patch faces with planar ends. Its trimmed side faces use
exact ribbon parameterizations. The shared-station mesher covers the remaining
complete patch faces, including single-chart atlases. This entry still requires
chord, embedding and closure proofs and returns `Answered<Mesh>`.

The fillet grid follows the two cap-cut graphs. Derivative queries reuse exact
Bezier span views translated by a constant origin. Complete derivative boxes
bound barycentric chord error across carrier knots. A fixed projection checks
positive symmetric secant matrices and output rounding. Cut-graph parameter
error and projection onto the exact cap planes contribute to these bounds.

`display_mesh_body` is a separate rendering entry returning
`OperationResult<DisplayMesh>`. It tries certification first. The fixed policy
assigns one quarter of the remaining caller allowance to that attempt, charges
its actual work, then permits the positional method to use what remains. No
attempt starts a fresh allowance. Without a caller allowance the whole request
has the unchanged 20M guard. Other tolerances and optional controls use the
certified entry. A failed direct attempt returns `mesh.spline-display-budget`
with the rule level and total request work, or a typed measured-bound refusal.

A successful certified attempt returns `Proven`. Otherwise supported spline
boundaries can return `Bounded` under `spline-display-evaluator-v1`, with kind
`spline display mesh` and level `1/100` mm. The evidence identifies the immutable
recipe/model revision and every covered original face. The result reports the
largest measured face/edge deviation and the measurement method. Unsupported
trims, a bound above the level, insufficient vertices or exhausted work return
a typed refusal. Bounded rendering data has no conversion to `Mesh` or export
writer. It supplies no solid embedding or exact-property proof.

The positional method evaluates a 3x3 stencil on each cell against the exact
retained carrier. It encloses the residual derivatives over every intersected
source knot span, then adds their axis contributions between stencil nodes to
the largest sampled residual. The L1 distance bound also bounds Euclidean
distance. The derivative cover uses outward homogeneous control enclosures and
positive weight bounds. Affine subtraction precedes division, so common weight
factors do not inflate the bound for a nearly straight direction. The larger
residual derivative contribution chooses the next refinement axis.

The bounded path omits the new ribbon and fillet embedding preparations. It
retains the finite-domain and shared-row checks needed to measure the right
boundary. Cut-domain displacement, shared-edge replacement and binary64 rounding
contribute to the final bound. Samples alone never establish acceptance. Body
and trace fit marks remain available in `DisplayMesh.source_marks`; the distance
measurement concerns the represented surface, rather than its ideal feature law.

`EvaluatedApplicationBody::display_result` exposes the rendering query without
consuming the retained checkpoint. `build_incremental_with_display_result`
builds the checkpoint, names and rendering answer within one allowance. Existing
certified application entries retain their types and contracts. SolidDraft must
handle all three rendering variants when adopting this entry.

`CutCurve.parameter_upper` reports the uniform error in the linear graph's tau
coordinate before spatial-control rounding. The STEP curve, spatial bound and
output controls retain their definitions. Callers constructing `CutCurve`
literals must supply this bound.

`RationalPatch::display_grid_adaptive` is an opt-in grid entry that retains
prepared span/cell bounds and bisects a failed cell's limiting interval. It
preserves exact knots and shared stations. `display_grid` retains its existing
refinement and output contract.

`BSplineSurface::derivatives_box` encloses first derivatives over an exact
rational box. It evaluates prepared derivative nets on each intersected owned
knot span and joins the enclosures. At a knot it includes both one-sided
derivatives. Existing binary-box entries keep their behavior.

`BSplineSurface::boundary_orientation` exposes the complete boundary-carrier
identity check without an interior embedding claim. The existing
`RationalPatch::boundary_orientation` delegates to it with the same behavior.
