# Target parts

Compare the same edit chains across kernels. The [roadmap](../ROADMAP.md) sets the order of work. The saved recipes under [`parts/`](../parts/) keep their dimensions so results can be compared across revisions.

| Part | What it exercises |
| --- | --- |
| P1 Bracket. A plate with a flange along one edge, a boss with a through bore, four bolt holes on a rectangular pattern and a slot through the flange. | Extrusions and booleans, including a union whose operands share a face, and a chain of about a dozen subtractions. |
| P2 Turned collar. A rectangle revolved a full turn into a ring, chamfered on both outer edges. | Revolution, torus and cone surfaces, chamfer on circular edges. |
| P3 Filleted boss. A block with a cylindrical boss on top, the junction edge filleted, the block's top edges filleted. | Plane-to-plane and plane-to-cylinder fillets, concave and convex, with transversal ends. |
| P4 Shelled tray. A block shelled to a thin wall with its top face open, then intersected with a cylinder to round its corners. | Shell with an open face, then a boolean on the shelled body. |
| P5 Flanged housing. A revolved housing with a pattern of countersunk holes and fillets where the flange meets the body. | Pattern, countersink cones in booleans, fillet against a revolved face. |

The next five are named so that no near-term decision closes them off: P6 lattice-filled bracket, P7 field-driven bracket with gyroid infill, P8 helical spring and threaded boss, P9 lofted duct, P10 topology-optimised bracket. Read the corresponding recipe or scoreboard row before claiming a part is qualified.

P0 is the hard-case set: equal-diameter pipes crossing at right angles, cylinders tangent along a line, a cylinder and a sphere with equal radii, a block and sphere union, poses along a seam, near-miss poses, and the general cone and torus poses that refuse today. Measure P0 at checkpoints. Fix refusals in the roadmap step that owns their cause. Every P0 case also runs as translated, scaled and seam-crossing copies of itself; a copy that answers differently is a gauge failure. Wrong successful results are caught by point-set membership as well as by volume, because a closed body can still be wrong.

A target part **builds** when:

- the owning kernel builds it from the recipe with no refusal, and background proof establishes every local tolerance under its target;
- its background proof completes within the applicable verification-request allowance; a postponed or failed proof does not count;
- after export proof, its STEP file opens in OpenCASCADE with the same face count and a volume within the exported deviation;
- a parameter edit such as plate thickness, boss radius or fillet radius builds the same way, with every recipe name resolving to the same face;
- the application loop builds and edits it in SolidDraft and passes, including unproved-preview, proof-completion and typed-failure behaviour.

Measure the complete selectable display and requested outputs, not just command return. Use the [roadmap measurement rules](../ROADMAP.md); the old stock, Boolean, feature and background proof times are not acceptance targets.

Field comparison is an on-demand check, not a universal "counts as built" gate. When a user or verification run requests it, record its outcome, budget and independent/taken-over coverage. No request means not requested, never implied agreement. Field-only parts still build in the field, with their applicable export or stand-in checks.

Use the same part and edit in each reference kernel. Keep reference comparisons in verification, not on every application edit.
