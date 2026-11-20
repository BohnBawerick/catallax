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
- Corrected in step 2, after a post-merge review: a path's nearest station gives the distance range directly only when the path declares itself uniform (`Path::uniform`): a straight path whose feet keep their offsets, or a circle whose feet turn with its meridian. Only then does the face point nearest a point lie on that station's arc. Along a straight centre line with turning feet, the projection's arc can be far while another station's arc nearly touches the point; `rolling_paths.rs` keeps that case as a test. Other paths take the panel cover. Strip and Ring are uniform, so their answers did not change.
- `query_batch` runs either answer on the shared field pool with the M18 index-based allowance split, so records and work are the same at any worker count. `distance_batch` runs the range alone, for points on the face whose membership a cover cannot settle.

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

## Step 2: a plane with a cylinder in any pose, a cone or a sphere

Code: `crates/catallax-field/src/rolling_region/revolved.rs`, `rolling_region/generator.rs`, `rolling_region/piece.rs` and `crates/catallax-field/src/fillet_rolling.rs`. Tests: `crates/catallax-field/tests/rolling_conics.rs`, and the step 1 tests.

Two new path kinds, derived by the field from the two supports, the ball side of each and the radius. Neither reads a B-rep fillet module or a fillet definition's section, contacts or carrier.

| Kind | Supports | Path | Closed forms |
| --- | --- | --- | --- |
| `revolved::Revolved` | a plane with a cylinder, cone or sphere, any pose | offset plane meets the coaxial offset surface: a circle when the plane is normal to the axis (always for a sphere), else an ellipse | circles: meridian station, separation `(rho - C)^2 + h^2 - r^2` (ring-torus polynomial for `C > r`, a square root for a spindle) |
| `generator::Generator` | two planes or cylinders parallel to one line | offsets meet in a line, solved in the section through the edge | projection clamped to the run; squared distance to the line less `r^2` |

The conic is written about the pivot where the axis meets the offset plane: `c = p + (G / D) (m - a cos (n.u1) / n.a)` with `m` the meridian, `D = n.a + k cos (n.u1)` and `G` constant. `(n.a)^2 (1 + k^2) > k^2` keeps the conic closed and `G D > 0` keeps its radius positive. An ellipse whose smallest radius is not above `r` is refused, so every ball stays off the axis. Panel enclosures of the meridian use the mean-value form at the panel's midpoint; a Taylor series over the whole panel was a millimetre loose near `pi`.

The fan is the edge's Voronoi zone in the solid of ball centres, as exact walls: the ball side of both supports; points whose projection along the plane normal onto the offset plane falls beyond the offset conic; and, in the point's own meridian, the edge side of the offset surface's normal line through the station. A straight piece has rail planes through each foot and end planes, which must be normal to its line. Bands keep far zones of these walls out and never cut the region: within the normal section of a station the region lies in the quadrilateral of the edge point, the feet and the centre, so it lies within `2r` of the plane, and within `2r` of a support the ball lies outside of. Between parallel supports the rail planes form an unbounded wedge, which crosses a cylinder the ball lies outside of and leaves it on the far side; that cylinder adds a chord wall through the edge and its rail, on the centre's side. This case was found by the test of the union of two parallel cylinders, whose other cusp the wedge reached.

The compile route keeps every pair and edge the elementary route sections today, so its answers and work are unchanged. The rolling route takes what it skipped or refused: ellipses, cones, spheres, parallel cylinder pairs, straight runs that share one name, and a straight strip whose section the elementary route cannot write exactly (an irrational centre) that meets no other run. The ball side of each support comes from the elementary route's exact material quadrant. The piece's region is the leaf `Leaf::Rolling`: membership from `path::classify`, `Indeterminate` where a cover cannot settle it, and a residual that is a lower bound of the distance to the region's boundary (the supports' 1-Lipschitz residuals, `|m - r|` with `m` the distance to the centre path, and outside the region the distance to its box). Every other support must be clear of the piece's box, proven by interval evaluation of its polynomial; one that reaches it is refused with `fillet.unsupported-overflow`, as the elementary route refuses it.

The pieces register no named face yet. Names, rays from boundary points and face witnesses for catallax-compare are open; see below.

How the reference kernels do this. OpenCASCADE builds plane-plane, plane-cylinder (line or circle spine) and plane-cone (circle spine) fillets in closed form (`ChFiKPart_ComputeData_FilPlnCyl`, `_FilPlnCon`, `_FilPlnPln`) and walks every other pair, the oblique ellipse included, with `BRepBlend_Walking` on `BlendFunc_ConstRad`. SMLib intersects the two offset surfaces (`SmFilletSolver`, `SmOffsetSurface` with an orientation per side) by stepping, then sweeps a cross-section. The field follows the offset-surface definition both use, solved in closed form for these pairs, and replaces walking and stepping with exact walls and a whole-interval cover, because a traced spine gives no guaranteed range.

Results (`rolling_conics.rs`): the field builds the step 2 corpus rows `cone-top`, `cone-base`, `sphere-dome`, `oblique-cylinder`, `narrow-cone` and `oblique-cone`, and agrees with every corpus witness (checked against OpenCASCADE). It builds the oblique file's cut rim (radius 1 and 2), boss and bore, the oblique-cone file's cut cone and its mirror, and the generator file's parallel-cylinder intersection, union and crescent, the flat's two generators and the single generator with a rational and an irrational centre, and agrees with those files' B-rep witnesses. Every display vertex of every B-rep blend face lies within the field's guaranteed range of the field's own face at 1 nm: about 2e-15 mm on circles and lines, and within 1e-7 mm on ellipses when a 1e-7 mm range is requested. The test asks only for a range as wide as the tolerance, which a face witness settles quickly; a narrow range on an ellipse costs thousands of panels, because a panel's lower bound is first order in its width. Records and work are the same at 1, 2 and 4 workers.

New and kept refusals, each named: `field.rolling-path-reach` (an open conic, a plane parallel to the axis for a revolved piece, an ellipse within `r` of its axis, a ball no smaller than the cylinder or sphere it rolls inside, offsets that do not meet), `fillet.unsupported-overflow` (another support reaches the piece), `field.fillet-contact` (a straight piece's end plane not normal to its line, a closed edge without a rational point for its material test, such as the coaxial file's dome rim, a circle of radius `sqrt 21` at height 2, which has no rational point because 21 is not a sum of two rational squares; the material probe needs algebraic points to clear it). Filtered edge references (`near`, `touching`) keep their elementary refusal.

## Step 3: quadric pairs

Code: `crates/catallax-field/src/rolling_region/coaxial.rs`, `rolling_region/saddle.rs`, the compile route in `crates/catallax-field/src/fillet_rolling.rs` and the algebraic material probe in `fillet.rs`. Tests: `crates/catallax-field/tests/rolling_quadrics.rs`, and the updated field assertions in `saddle_fillets.rs` and `coaxial_fillets.rs`.

Two new path kinds, derived by the field from the two supports, the ball side of each and the radius. Neither reads a B-rep fillet module or a fillet definition's section, contacts or carrier.

| Kind | Supports | Path | Closed forms |
| --- | --- | --- | --- |
| `coaxial::Coaxial` | two cylinders, cone nappes or spheres about one axis | the offset profiles meet in the meridian half-plane: Cramer's rule for two lines, a quadratic along a line for a line and a circle, the radical line for two circles | the circle's: meridian station, ring-torus separation |
| `saddle::Saddle` | a cylinder with a second cylinder or a sphere off its axis | along each ruling of the offset cylinder, at angle `theta`, the other offset fixes the axial station by `A s^2 + 2 B s + C = 0`, so `s = (-B + sigma sqrt(D)) / A`: one ring (`D > 0` everywhere, one sign), one loop (both signs where `D >= 0`), or, for two equal cylinders whose axes meet, two ellipses in the bisecting planes (`D` a square, `s` linear) | none; whole-interval panel cover |

Both follow the closed form the decisions name. Dupont, Lazard, Lazard and Petitjean (2008) parametrize the meeting of two quadrics through a ruled quadric of their pencil, with at most one square root for a smooth quartic; when one support is a cylinder it is that ruled quadric, and its rulings are the angle parameter. OpenCASCADE's `IntAna_Curve` (`SetCylinderQuadValues`, read at V7_8_0) writes a cylinder meeting a quadric the same way: the cylinder's angle, two roots of a quadratic in the axial coordinate (`TwoZForATheta`, `ZIsPositive`) over the domain where they are real. It then fillets such pairs by walking (`BRepBlend_Walking` on `BlendFunc_ConstRad`); SMLib intersects the offset surfaces by stepping (step 2). The field keeps the closed form and replaces walking with exact walls and a whole-interval cover, because a traced spine gives no guaranteed range.

Coaxial pieces work in the meridian half-plane, where the fan is the sector at the centre between the rays to the two feet. A sphere the ball lies outside of adds its chord wall, as the cylinders of step 2 do. Of several meetings of the offsets the piece keeps the one whose feet lie ahead of the edge along both faces; exactly one must. Spindles (`C <= r`) are refused.

Saddle pieces: a ring keeps the root sign its edge has at the edge's witness point (`A s + B` there); a loop needs `D` to change sign on one run of the circle, proven on 256 panels, with zeros only at the run's two ends, for the edge and for the offsets alike. Otherwise the ball changes the topology (offsets that touch, a loop that pinches to a node) and the piece refuses. Near a branch point `ds/dtheta` grows without bound, so the frame's tangent is `sqrt(D)` times the derivative. Panels where `D < 0` hold no station: the shared cover now asks `Cover::exists` before a midpoint is used as a witness and skips panels `Cover::empty` proves empty (both default to the old answer, so earlier paths are unchanged). The fan's walls are the ball side of both supports and the two projection walls (a point's radial projection onto either offset falls outside the other offset), all exact. A box keeps far zones out, such as the other ring of the offsets: the hull of the centres, feet and edge points over the cover, with the chord cube of each support the ball lies inside, padded by `r`. It is a band (`Path::bands`, new, default none): positive in the region, never a face. The crossing's single piece covers both ellipses, so every ball of the crossing clips all four bands, which meet at the poles.

No closed form gives a saddle's nearest station, so each distance query first takes a binary64 guess (`Path::guess`, new, default none): the centre nearest the point among 1024 samples, polished by a golden-section search. The query evaluates the exact face point at that station before the cover; the guess enters a bound only through that exact witness. A display vertex then settles at once: the saddle file's large loop has 58,630 blend vertices, which took 872 s on the Mac without the guess and 13 s with it.

Clearance: a ring's whole box holds its axis, so a support near the ring's inside reached it. When the whole box fails, every panel's box is tested; a support clear of all of them is clear. This also clears the coaxial file's chamfer, whose top round's box reached the cylinder wall, and keeps every earlier clearance, which the whole box already proved, unchanged.

Algebraic edge points. Edges with no rational point (the dome rim, a circle of radius `sqrt 21`; every saddle ring) are probed for their material at an algebraic point: the edge's witness, or the end of a closed edge's axis. The probe is the rational one through the field's algebraic membership and ray queries (`distance::location`, `distance::ray`, `Leaf::on_face_algebraic`); rational points take the old probe unchanged.

Scope of the admission. Only named edges admit the step 3 pairs. An all-edges selection skips them as it did before, so its answers do not change; whether it should round them is left to a later step.

Results (`rolling_quadrics.rs`). The field builds every case of the four files that the B-rep builds, except the countersink's mouth and the ring beside a clear neighbour (below), and agrees with every witness of those files: the coaxial file's frustum, chamfer, dome, seat, throat and rod cap; the saddle file's boss, bore, larger ball, either order, lean, through bore, mirror, and its loops (boss, notch, lens, bite, larger ball, either order, the far side, a leaning boss and notch, equal radii, almost past, a mirror), with the loop witnesses from that file's own ball model; the sphere ring and its mirror; the equal crossing as union and intersection, oblique and spatial, reflected, and the offset and unequal crossings, which are a loop and two rings. Every display vertex of every B-rep blend face that displays lies within the field's guaranteed range of the field's own face, at the part's tolerance: 1 nm on exact faces, and on fitted saddle bands their certified bound, as rules T0 to T8 store it (for the saddle boss 2.26e-5 mm; the field's ranges put the B-rep's vertices at up to 2.24e-5 mm from the exact face). Records and work are the same at 1, 2 and 4 workers. The whole test file takes 98 s on the Mac.

A fix to step 2. A plane with a cone whose ball path is a circle did not check the ball's foot on the cone. In the coaxial file's spike (a pointed cone on a block, radius 12) the foot lies past the apex, on the other nappe, so no ball touches the face; the B-rep refuses and the field built a piece. `Revolved::new` now refuses it (`field.rolling-path-reach`). Every step 2 case keeps its answer.

Refusals, each named. Where the B-rep builds and the field refuses, with `fillet.unsupported-overflow`: the countersink's mouth and both rounds of it, because the hole's entry cylinder passes through the mouth circle into the region and only its face's extent keeps it clear (the field has no face-extent proof for a support that meets a piece along its edge); and the saddle file's ring beside a clear neighbour, whose band stops about 0.05 mm short of the neighbour while the field's panel boxes, which hold the ball centres, reach it. Where both refuse, `the_field_names_what_it_refuses` pins the field's code for each: the bore at radius 3 and 10, the short narrow and short wide cylinders, the neighbour in the way, the loop past a cylinder's end, offsets that touch, the clipped and open crossings, a crossing with the ball on different sides, the dome at radius 2 and 3/2, the narrow frustum at radius 2, the spike, and the torus groove and slice (step 6). `field.rolling-path-reach`: offsets that change the edge's topology, a crossing whose ball lies inside one cylinder and outside the other, a spindle, two balls or none in a coaxial corner. `field.fillet-contact`: a crossing clipped to fewer than four bands, a fillet face's own classifier at an algebraic point, and the torus groove (step 6).

Residual: a saddle piece has no nearest station, so its residual near the path is the supports' and, outside its box, the box's distance; inside the box it is zero.

## Open for the done-when line

The done-when line asks catallax-compare to find every B-rep display vertex within the field's bound. Compare's forward pass needs a field face witness to bound a vertex from above, and a rolling piece registers no named face. Until the pieces' faces join the field's face registry (names, directions, prepared witnesses from `path::distance`), compare leaves these vertices unresolved; `rolling_conics.rs` runs the same check directly against the pieces. This building block is needed before the milestone's done-when can be shown through compare.

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
| 1 | General path routine; line and circle pieces match today's answers | done ([#716](https://github.com/BohnBawerick/catallax/pull/716)) |
| 2 | Plane with cylinder in any pose, cone or sphere | [#724](https://github.com/BohnBawerick/catallax/pull/724) |
| 3 | Quadric pairs (closed forms per Dupont, Lazard, Lazard, Petitjean 2008) | this change |
| 4 | Ends at an edge or vertex, P5 at radius 2 | |
| 5 | Fillet over fillet, nested path | |
| 6 | Torus pairs | |
| 7 | Spline faces (7.0 to 7.7 of the decisions) | |
