# Target parts

Moved on 14 November 2026 from section 9 of the plan now archived as [PLAN-2026-11.md](archive/PLAN-2026-11.md). [ROADMAP.md](../ROADMAP.md) orders the work. The "first stage" column refers to the staged path in [section 6 of the archived plan](archive/PLAN-2026-11.md#6-the-staged-path).

Target parts are the scoreboard and a tie-breaker. Work proceeds in the build order whether or not a target part needs it yet; when two items compete for the same hands, the one an earlier unbuilt part needs goes first. P1 to P10 below are confirmed ([section 15 of the archived plan](archive/PLAN-2026-11.md#15-decided-on-22-october-2026)). Each part is a `catallax-recipe` file under `parts/` with fixed dimensions, so timings compare across revisions.

| Part | What it exercises | First stage that can build it |
|---|---|---|
| P1 Bracket. A plate with a flange along one edge, a boss with a through bore, four bolt holes on a rectangular pattern and a slot through the flange. | Extrusions and booleans, including a union whose operands share a face, and a chain of about a dozen subtractions. | 1 |
| P2 Turned collar. A rectangle revolved a full turn into a ring, chamfered on both outer edges. | Revolution, torus and cone surfaces, chamfer on circular edges. | 3 |
| P3 Filleted boss. A block with a cylindrical boss on top, the junction edge filleted, the block's top edges filleted. | Plane-to-plane and plane-to-cylinder fillets, concave and convex, with transversal ends. | 3 |
| P4 Shelled tray. A block shelled to a thin wall with its top face open, then intersected with a cylinder to round its corners. | Shell with an open face, then a boolean on the shelled body. | 3 |
| P5 Flanged housing. A revolved housing with a pattern of countersunk holes and fillets where the flange meets the body. | Pattern, countersink cones in booleans, fillet against a revolved face. | 4 |

The next five are named so that no near-term decision closes them off: P6 lattice-filled bracket, P7 field-driven bracket with gyroid infill, P8 helical spring and threaded boss, P9 lofted duct, P10 topology-optimised bracket. Their recipes are written when their turn comes.

P0 is the hard-case set: equal-diameter pipes crossing at right angles, cylinders tangent along a line, a cylinder and a sphere with equal radii, a block and sphere union, poses along a seam, near-miss poses, and the general cone and torus poses that refuse today. P0 is both a gauge, measured at every milestone, and a work queue for the refusal-repair lane. Every P0 case also runs as translated, scaled and seam-crossing copies of itself; a copy that answers differently is a gauge failure. Wrong successful results are caught by point-set membership as well as by volume, because a closed body can still be wrong.

Under the 11 November ruling ([section 1.2 of the archived plan](archive/PLAN-2026-11.md#12-capability-first-speed-last)) there are two levels. A target part **builds** when:

- the owning kernel builds it from the recipe with no refusal, and background proof establishes every local tolerance under its target;
- its background proof completes within the applicable unchanged deterministic allowance; a postponed or failed proof does not count;
- after export proof, its STEP file opens in OpenCASCADE with the same face count and a volume within the exported deviation;
- a parameter edit such as plate thickness, boss radius or fillet radius builds the same way, with every recipe name resolving to the same face;
- the application loop builds and edits it in SolidDraft and passes, including unproved-preview, proof-completion and typed-failure behaviour.

A target part **counts as built at speed** when it builds and also:

- construction preview and complete selectable display, including mesh, edge polylines and names, meet the 100 ms stock, 500 ms Boolean or 1,000 ms feature-click target, for the build and for the parameter edit;
- its background proof completes within the separate 2,000 ms proof budget.

The capability milestones ask that a part builds. Milestone 5.1 of [the roadmap](../ROADMAP.md), speed and simplification, asks that every target part counts as built at speed. Timings that count come from a machine to itself ([section 10.5 of the archived plan](archive/PLAN-2026-11.md#105-jobs-that-need-a-machine-to-themselves)).

Field comparison is an on-demand check, not a universal "counts as built" gate. When a user or verification run requests it, record its outcome, budget and independent/taken-over coverage. No request means not requested, never implied agreement. Field-only parts still build in the field, with their applicable export or stand-in checks.

The scoreboard also carries the head-to-head columns of item 5 of the 22 October list ([section 11.4 of the archived plan](archive/PLAN-2026-11.md#114-the-22-october-next-items-on-11-november)): the same parts built and edited in a reference open-source CAD tool and in Catallax, in the scoreboard run, never in the application.
