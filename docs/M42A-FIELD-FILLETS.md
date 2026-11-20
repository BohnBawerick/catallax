# 4.2a: the field builds every fillet on its own

This is the milestone note for ROADMAP 4.2a and [Checkpoint 2 decisions, section 4.6](checkpoints/checkpoint-2-decisions.md#46-item-42a-the-field-builds-every-fillet-on-its-own) (Q15 to Q17, Q21, Q23). The decisions and the roadmap win over this note.

Scope: constant-radius rolling-ball fillets only. Variable radius, setback vertex blends and conic or G2 sections keep their typed refusals.

## The independence rule and the shared dependencies

The decisions say: "the field derives its own ball paths. It does not reuse the B-rep's `rolling_pipe`, `rolling_saddle`, `rolling_sphere` or `rolling_mitre`, because a shared bug would fool both kernels." No 4.2a code reads those modules.

Code the field shares with the B-rep, and what a bug in it would do:

| Shared code | Used by the field for | Risk |
| --- | --- | --- |
| `catallax-numbers` (`Real`, `Ball`, `Decider`, allowances) | all exact signs and interval enclosures | A bug here breaks many unrelated tests at once, so it is unlikely to hide. |
| `catallax-shapes` points, vectors and surface records | inputs and outputs | Plain data, no fillet logic. |
| `catallax_shapes::RollingEnvelopeDefinition` | the P5 rim piece's validated dimensions, frame and normal-cone regime | The field computes the rim centre curve itself from these numbers. A wrong regime check would admit or refuse the same rims in both kernels. Step 4 replaces it with a field-side check or records why not. |
| `catallax_features::fillet::FilletDefinition` and `CornerDefinition` | today's compiled plane-plane and plane-cylinder field fillets (edge, contacts, section) | These elementary definitions predate 4.2a. The step 1 path pieces do not use them; the step 1 tests compare the two, so a bug in either shows. Production keeps them until a later step moves the compile route onto the path pieces. |
| `catallax-intersect` support intersection | finding which support pairs meet at a sharp edge | Edge discovery only. The ball path is not computed there. |
| M10 pair search (`offset_meeting`) | step 7 only, spline spine topology | Named by the decisions. Step 7 adds the dense binary64 third check that catches a bug shared through M10. |

## Query contract (Q23)

A field query returns inside or outside plus a guaranteed range for the distance near the surface. The rolling region of M18 answered inside or outside only. Step 1 extends it.

## Step 1: one general routine for a path piece

Code: `crates/catallax-field/src/rolling_region/path.rs`, `rolling_region/elementary.rs` and `rolling_region.rs`. Tests: `crates/catallax-field/tests/rolling_paths.rs`, with the M18 tests `rolling_region`, `rolling_sectors` and `m18_p5_verdict` unchanged.

A path piece is a finite ball-centre curve `c(u)` for `u` in `[-1, 1]`, a ball radius `r` and the two contact feet. Its face is the union of the short arcs of the characteristic circles between the feet. Its region is its fan, given by exact clip walls, less every ball. Each path supplies its local frame, its walls and interval enclosures of its centre, tangent and feet over a parameter panel. A path with closed forms also supplies its separation from the balls and its nearest station.

- `path::classify` decides the walls exactly. With a closed-form separation it decides the balls exactly, `Boundary` included. Otherwise it runs the M18 whole-interval panel cover: a panel is cleared only when every ball of the panel misses the point, and a midpoint inside a ball proves exclusion. Sampling never proves inclusion. The limits are M18's: 96 bits, 4096 panels, depth 64.
- `path::distance` returns `[lower, upper]` for the distance to the face. For one characteristic circle the distance is `sqrt((|v_n| - r)^2 + |v_t|^2)`, with `v` the offset from the centre split into its normal-plane and tangent parts. It is the distance to the arc when `v_n` points into the arc's sector, and otherwise the nearer foot is nearest. The lower bound of a panel uses the panel's enclosures. The upper bound comes from exact face points at panel midpoints. Panels are refined lowest bound first, ties by parameter, until the range is at most the requested width or the limits are met; `converged` says which. The range is guaranteed either way.
- `query_batch` runs either answer on the shared field pool with the M18 index-based allowance split, so records and work are the same at any worker count.

Path kinds in step 1:

| Kind | Derived from | Closed forms |
| --- | --- | --- |
| `Region` (M18) | P5 rim: body cylinder offset and mouth-circle tube | none; panel cover, lower-envelope semantics kept |
| `elementary::Strip` | two planes offset by `r` toward the ball meet in a line | projection onto the line; squared distance to the line less `r^2` |
| `elementary::Ring` | a plane and a perpendicular cylinder offset by `r` meet in a circle of radius `C = R + r` or `R - r` | the meridian through the point; the ring-torus polynomial, valid for `C > r` |

`Region::classify` now runs through `path::classify`. Its decider calls are in the same order, so every answer it gave before, and its work, are unchanged; the M18 tests pass unchanged. One refusal can now become an answer: a panel at the depth limit (a point tangent to one ball of the family) used to end the query with a refusal, and now the search for an exclusion witness continues over the rest of the cover. Inclusion is still refused if any panel stays unresolved. Production field fillets do not use the new pieces yet, so every compiled field answer is unchanged by construction.

The step 1 tests compile the box-edge, oblique (126.87 degree), concave L, convex ring and concave boss-rim fillets and check, at every grid point the fillet alone decides:

- the piece's membership equals the compiled field location (bit-identical, exact on both sides);
- the closed-form range contains the exact distance from an independent section-plane formula, and is converged at 1e-12 mm;
- the same piece run through the general panel cover gives the same membership off the boundary, and its range contains the same exact distance.

The P5 rim range is checked against a dense binary64 sampling of the face written from the M18 chart formulas.

How the reference kernels do this: OpenCASCADE projects a point on a line or circle in closed form (`Extrema_ExtPElC`; an ellipse through `math_TrigonometricFunctionRoots`) and on a general curve by sampling the curve and refining where `(C(u)-P).C'(u)` changes sign (`Extrema_GExtPC`, 17 samples per interval, degree + 1 for a B-spline). Its rolling-ball blend (`BlendFunc_ConstRad`) puts each section in the plane normal to the guide at the guide point. Step 1 follows the closed forms for lines and circles. For general paths it replaces the samples with a whole-interval cover, because a sampled minimum is not a guaranteed range. SMLib's source is not available to this lane, so no SMLib approach was checked.

New refusal: `field.rolling-path-reach` names a path the field cannot derive: parallel planes, a plane not perpendicular to the cylinder axis, a non-positive radius, or a centre circle no larger than the ball (horn or spindle tube).

Not covered by step 1: the straight generator between a cylinder and a parallel plane (a line piece whose section follows the cylinder arc) goes with the other cylinder poses in step 2. Three-plane corners and planar corner fans are vertex pieces and stay on their current route until step 4.

## Frozen case list for the done-when lines

Frozen at `945f2504a92b82267e2c615c3f122f99e8714d69` (20 Nov 2026 main).

From [the fillet inventory](../parts/fillets/cases.jsonl), the constant-radius rolling-ball rows. `variable-linear`, `variable-zero-end` (variable radius) and `explicit-setbacks` (setback) are out of scope.

B-rep builds, so catallax-compare checks the field against it:

| Step | Rows |
| --- | --- |
| 1 (plane-plane strips, plane-perpendicular-cylinder rings) | `box-single`, `near-face-limit`, `at-face-limit`, `past-face-limit`, `swallowed-face`, `acute-5deg`, `acute-halfdeg`, `obtuse-174deg`, `obtuse-179deg`, `l-concave`, `cylinder-rim`, `cylinder-collapse`, `cylinder-past-collapse`, `boss-rim`, `boss-seam-rotated`, `hole-mouth`, `r12-cube-at-contact`, `r12-cube-past-contact` |
| 2 (plane with cylinder in any pose, cone or sphere) | `cone-top`, `cone-base`, `sphere-dome`, `oblique-cylinder`, `narrow-cone` |
| 4 (ends at an edge or vertex) | `box-all`, `l-mixed-all`, `vertex-4-one` |
| 5 (fillet over fillet) | `box-mitre`, `tangent-boss`, `tangent-pocket`, `unequal-ordered`, `roll-over-mitres`, `mixed-ordered` |
| 6 (torus pairs) | `torus-slice`, `torus-intersection` |

B-rep refuses today, so the field result is recorded with its reason and there are no display vertices to compare: `cross-equal`, `cross-unequal`, `oblique-cone`, `sphere-cylinder`, `vertex-4-two`, `vertex-4-all`, `partial-oblique`, `unequal-mitres`, `narrow-mitres`, `vertex-4-all-unequal-slopes`, `r12-legacy-plane-bound`, `r12-rim-rational-bound`, `r12-generator-bound`.

From the fillet test files the decisions name, every constant-radius recipe in these files at the frozen commit, under `crates/catallax-features/tests/`: `oblique_fillets.rs`, `oblique_cone_fillets.rs`, `generator_fillets.rs` (step 2); `saddle_fillets.rs`, `coaxial_fillets.rs`, `fillet_sphere_ring.rs`, `equal_crossing_fillets.rs` (step 3); `rim_contact_fillets.rs`, `hole_ring_fillets.rs`, `slot_fillets.rs` and P5 at radius 2 (`parts/P5-flanged-housing.json`) (step 4); `fillet_over_fillet.rs` (step 5).

The list does not change during the milestone. A case added later joins the next milestone's list.

## Steps

| Step | Content | State |
| --- | --- | --- |
| 1 | General path routine; line and circle pieces match today's answers | this change |
| 2 | Plane with cylinder in any pose, cone or sphere | next |
| 3 | Quadric pairs (closed forms per Dupont, Lazard, Lazard, Petitjean 2008) | |
| 4 | Ends at an edge or vertex, P5 at radius 2 | |
| 5 | Fillet over fillet, nested path | |
| 6 | Torus pairs | |
| 7 | Spline faces (7.0 to 7.7 of the decisions) | |
