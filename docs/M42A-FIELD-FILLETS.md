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
| `catallax_shapes::RollingEnvelopeDefinition` | the M18 `Region` building block's input, in the M18 tests only | Replaced in step 4: the field's rim piece (`rolling_region::rim`) finds each mouth from its own leaves and supports and checks the regime itself. The M18 `Region` still takes the B-rep's definition as its input format; no recipe compiles through it. |
| `catallax_features::fillet::FilletDefinition` and `CornerDefinition` | today's compiled plane-plane and plane-cylinder field fillets (edge, contacts, section) | These elementary definitions predate 4.2a. The step 1 path pieces do not use them; the step 1 tests compare the two, so a bug in either shows. Production keeps them until a later step moves the compile route onto the path pieces. |
| `catallax-intersect` support intersection | finding which support pairs meet at a sharp edge | Edge discovery only. The ball path is not computed there. |
| `RationalPatch::bezier_spans` | step 7, exact homogeneous Bézier controls for patch distance and the spine's polynomial charts | Both kernels use the same blossom extraction. The step 7.1 distance check uses binary64 points from the original support polynomials as its oracle. |

M10's pair search informed the design. The field computes and certifies its
own seed, trace and border roots. The dense binary64 check solves the contact
equations from the original support polynomials independently.

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

## Step 4, part 1: rim contact

Code: `crates/catallax-field/src/rolling_region/rim.rs` and the compile route `crates/catallax-field/src/fillet_rim.rs`. Tests: `crates/catallax-field/tests/rolling_ends.rs`, and the updated M18 tests `m18_p5_verdict.rs` and `m18_clearance.rs`.

One new path kind, derived by the field from its own supports and leaves. It reads no B-rep fillet module, no fillet definition's section, contacts or carrier, and not the B-rep's `RollingEnvelopeDefinition`.

| Kind | Supports | Path | Closed forms |
| --- | --- | --- | --- |
| `rim::Rim` | a plane and a perpendicular cylinder, the ball outside it, with circular mouths in the plane | the ring's centre circle of radius `C = R + r`; over a mouth of radius `m` the ball leaves the plane for the mouth's edge and its centre drops to `H = sqrt(r^2 - (m - D)^2)`, `D` the horizontal distance to the mouth's axis | the ring's, when no mouth is crossed |

The path is one piece for the whole ring, mouths included, so no seam between a mouth and the ring needs its own proof. Its region is M18's: in the ring's frame, `0 <= y <= r` and `R <= rho <= C`, below the lower envelope of every ball. The interval cover writes the dip as the gap `g = sum of max(0, m_k - D_k)` over the mouths, so a panel that straddles a mouth's sector holds both the ring's and the dipped stations; no station comes nearer a mouth's axis than `o - C`. The edge contact is the edge point nearest the centre, which lies on the ball's characteristic circle, so the face's arc runs from the cylinder foot to it.

The route runs only where the elementary ring refuses with `fillet.unsupported-overflow`, so every earlier answer is unchanged. It proves every other source face either misses the ring's fan or is the side of a mouth:

- a support misses the fan when its exact second-order form about each of 64 panel centres, in the ring's frame, excludes zero on the panel; a support of higher degree is evaluated on the panel's world box;
- a face misses it when its leaf does: a Round primitive parallel to the ring's axis is compared with the fan exactly (heights against `[0, r]`, its disc against the annulus `[R, C]`), any other Round primitive by its box;
- a mouth is a Round primitive parallel to the axis, wholly on the material side of the plane, with one cap in it. Its cap circle is the edge; its side face's depth and flare give the edge's normal cone. A ray probe just below the plane, inside the mouth and the fan, separates a hole (air) from a buried cap (material), and one above shows the ring's side open. Mouth discs must keep apart;
- the ball must cross a mouth on its near side (`R < o - m < C < o`), stay above it (`C - (o - m) < r`), stay on the edge rather than the side face (`depth H + flare g > 0` at the deepest station), stay above the side face's far end, and the mouths' sectors must keep apart (the test `Sectors::new` uses).

A mouth the centre circle only touches (P5 at radius 1) or clears leaves the ring unchanged: the piece is the ring with its closed forms.

How the reference kernels do this. OpenCASCADE rolls a ball on a surface and an edge with `BRepBlend_SurfRstConstRad`, its constant-radius function with the edge's curve on the other face as the second support, walked by `BRepBlend_SurfRstLineBuilder` (read at V7_8_0). SMLib (the Apache-2.0 tree in `NVIDIA-Omniverse/usd-brep`) names the case a cliff rollover (`SM_FG_CLIFF_ROLLOVER` in `SmFilletGeom.h`), found where a rail crosses an edge (`SM_FV_RAIL_X_EDGEUSE`), and steps it like its other blends. Both take the edge as the second support, as the field does; the field keeps the closed form of that contact, the tube of radius `r` about the edge circle, and replaces walking and stepping with exact walls and a whole-interval cover.

Results (`rolling_ends.rs`). The field builds every P5 row of the rim-contact file that the B-rep builds: radius 2 and 5/2, translated, scaled by 1/2 and 2, turned (RZ) and turned in space (RQ), with and without a move; the four rows whose rail touches every mouth (radius 1); and the four whose rail clears them (the elementary torus, as before). It agrees with every membership witness of that file except the points in the flange plane under the round and at the contacts, which lie on a face the piece shares with the housing: a rolling piece registers no face to settle them. Every display vertex of every B-rep blend face lies within the field's guaranteed range of the rim piece's face, at the part's tolerance: 1 nm on the exact tori, the certified bound on the fitted envelopes. On P5 at radius 2 the worst of 3,192 vertices is 9.2e-8 mm from the field's face; on the touching rows 8e-15 mm. Records and work are the same at 1, 2 and 4 workers.

The hole-ring file's 17 rings already built on the elementary route. The test checks each ring's two side points from that file and every B-rep blend display vertex against the field's own `elementary::Ring` path, within 3.3e-15 mm.

Refusals, each named, where both kernels refuse: P5 at radius 7/2 (the B-rep: the ball would leave the mouth's edge for the cone; the field: the drill's top, 1 mm above the flange, reaches the wider fan first), and the three rows with a wall lower than the ball (its top plane reaches the fan), all `fillet.unsupported-overflow`.

## Step 4, part 2: tangent junctions and runs past the rails

Code: `crates/catallax-field/src/fillet_chain.rs`, the junction clips in `rolling_region/path.rs`, junction ends in `rolling_region/generator.rs`, ring sectors in `rolling_region/revolved.rs`, the face-extent proof `Feature::face_above` in `features.rs`, and the clearance in `fillet_rolling.rs`. Tests: `crates/catallax-field/tests/rolling_ends.rs`.

The route runs only where the other routes refuse a named selection with `fillet.unsupported-overflow`, `fillet.unsupported-contact` or `field.fillet-contact`. It starts again from the plan as it stood before them, so every earlier answer is unchanged. It reads no B-rep fillet module and no fillet definition. It never takes a body that an earlier fillet or chamfer rounded; that is step 5.

The field finds each selected pair's sharp runs from its own supports. A line is cut where other supports cross it, as before. A circle is cut at every point where another support crosses or touches it, each arc is probed for its material at its midpoint, and neighbouring arcs of one material are joined. The recipe version decides what happens at a run's end, as the B-rep's builder reads it: version 4 on continues a run across a tangent neighbour, version 3 runs a straight strip past its rails.

A tangent junction. Two runs share one face with the same material side and are both convex or both concave. Their other faces are tangent at the shared end, with the same material side, and the runs leave that point in opposite directions. The ball there is then one ball, and the two runs' sections in the plane normal to them at that point are one section. That plane is a junction: a straight piece (`Generator`) ends there with an end plane, a circle piece (`Revolved::with_sector`) with a meridian half-plane. A junction clips the region but is no face (`Path::junctions`). A point on it inside the fan and outside every ball lies inside the union of the two pieces, so a zero there does not make it a boundary point. A ring sector spans at most half a turn. Its nearest station is the meridian through the point, or the nearer end outside the sector, because the distance to a station's centre, and to its arc, grows with the angle between the two meridians. Version 4 extends the selection as the B-rep does: from each free end, along the face the run keeps, onto every face through that end that is tangent to the face it leaves.

A run past the rails (version 3). A straight convex run whose two ends meet cylindrical end walls becomes a strip over the whole line, past every leaf. The end walls are exempt from its clearance, so the body cuts the strip. Each end wall's face must lie beyond its run's end. The ball at each run end must lie inside the end wall's solid, on the side that holds its centre (`rho + r <= R` inside, `rho - r >= R` outside). That second condition is the slot file's stated reason for its radius 6 refusal: the ball does not fit inside the radius 5 end wall.

Face extent. At a tangent junction the neighbour's support passes through the piece's box, so the support test refuses. The field clears the neighbour's face from its leaf instead. An extrude's side face lies in its profile segment swept along the extrusion. An affine function's least value over that face lies at one end of the sweep, and on the segment at an end or, on an arc, where the arc's radius points against the function's gradient. A face on the far side of a junction meets the region at most in the junction plane, which holds no face of the region. The slot's end wall, a half circle swept up, lies beyond the junction plane through its two tangent points, and the slot's sides lie before it.

How the reference kernels do this. OpenCASCADE's `ChFi3d_Builder::PerformElement` (read at V7_8_0) collects "all mutually tangent edges" into one spine, and requires the faces opposite two adjacent edges to be tangent; the blend is then walked along that spine. SMLib's `SmFilletExecutive::ProcessSmoothEdgeSequences` collects smoothly meeting filleted edges (`CollectTangentEdgeSequence`), joins their curves into one, and adjusts the tangent corners between them (`AdjustTangentCorner`). The field keeps the same chain rule. It cuts the chain at each tangent junction into pieces with closed forms, a straight strip and a ring sector, whose sections agree there, so nothing is walked. For the version 3 run past the rails the field follows the recipe's rule, the strip cut by its end walls; neither reference kernel has recipe versions.

Results (`rolling_ends.rs`). The field builds every case of the slot file that the B-rep builds: the strip along the top at radius 1, 2 and 4 and the strips along the bottom and top at radius 1 and 4 (version 3), and the tangent loop of two strips and two half rings at radius 1, 2 and 4, named whole in version 3 or 4, continued from one arc in version 4, and named in another order with a repeat. At radius 4 the half rings are spindle sectors (`C = 5 - 4 < r`), which `Revolved` holds with its square-root separation. The field agrees with every membership witness of that file, and with points derived for the strips past their rails, each of which the test first checks against the B-rep. Every display vertex of every B-rep blend face lies within the field's guaranteed range of a piece's face at 1 nm; the worst is 5.7e-15 mm. Records and work are the same at 1, 2 and 4 workers.

Refusals, each named, where both kernels refuse: the loop at radius 6 (`field.rolling-path-reach`: the ball is no smaller than the end wall it rolls inside) and the strip at radius 6 (`fillet.unsupported-overflow`: the ball at a run end crosses its end wall).

Shared points. Points that lie on a face a rolling piece shares with the body (P5's flange plane under the round and its four contacts, the slot's rail ends) are kept per case with the B-rep's answer. A rolling piece registers no named face yet, so the field leaves them `Indeterminate`; the test checks that it never puts one on the wrong side. They are the control cases for the face registry below.

## Step 4, part 3: a mixed vertex and a strip end at a four-face vertex

Code: `crates/catallax-field/src/fillet_mixed.rs`, the run merge and the buried-face clearance in `crates/catallax-field/src/fillet.rs`, leaning end planes in `fillet_rolling.rs` (`Plan::oblique_stations`), and `Generator::section_at` and `Generator::enclose_cover` in `rolling_region/generator.rs`. Tests: `crates/catallax-field/tests/rolling_ends.rs`.

`vertex-4-one` rounds one ridge of two intersected triangular prisms, which ends at the four-face apex. The two apex walls lean across the ridge. The B-rep keeps its exact cylinder and cuts it by those walls (interface note of 7 Nov, "a single planar strip ends at a four-face vertex"). The elementary route widens a strip's span for a leaning end plane, but it cannot write this strip's section, whose frame along the diagonal is irrational, so the strip became a straight piece (`Generator`), which needed end planes normal to its line. A straight piece between two planes now ends at a normal plane past a leaning end plane. Each section of its region lies in the triangle of its edge point and two feet, so past the station where one end plane has passed all three, the section lies beyond that plane; each end takes the plane that passes them first. The body must lie on the run's side of every leaning end plane: it is an intersection of certified convex planar leaves (`convex_planar_input`, as the four-contact fans use) whose planes include the end plane. Then no material lies past an end plane, and the region removes the strip as the walls cut it. An irrational station is rounded outward to a dyadic rational, which adds only air. Both supports must be planes; with a cylinder the refusal stays.

`l-mixed-all` rounds every edge of an L of two blocks. Two things stopped it. First, the end faces of the two blocks are coplanar, so each edge along them came back from several face pairs and five or more runs met at a corner. An all-edges selection now keeps one run per edge: runs with the same supports, convexity and material sides whose intervals meet are one run, named by the first pair. Second, a strip along one block crosses the other block's face, buried in the union, such as the wall's front face below the base's top. A union with a child that holds the strip's hull equals that child inside it, and an intersection with such a child equals the other child, so the left-out child's faces cannot bound the body there (`Plan::bounding_leaves`; a child holds the hull when it is a certified convex leaf whose closure holds every hull point). This runs only where the strip would otherwise refuse.

Then the two ends of the concave edge are mixed vertices: one concave and two convex runs. RECIPE section 4.5.1 settles them: the concave blend forms first and the convex rounds roll over it, and the patch over the concave blend is the vertex's `corner.vertex`. The field builds that from its own pieces. The concave strip keeps its region. Each convex strip stops where the concave strip's foot on its face crosses the end face `P`, at a junction. Between the two junctions the ball rolls on `P` and the concave blend's cylinder, whose axis is the concave strip's centre line, derived by the field's `Generator`. With `P` normal to that line this is a ring sector (`Revolved::with_sector`) between the meridians through the two feet. At each junction the ring's ball is the convex strip's ball. The ring's fan holds only material once the concave strip is added: every support except the concave strip's two clears the ring's box, eight rays from the vertex prove the body there is the concave corner cut by `P`, and the sector between the feet is the wedge between the concave strip's rail planes, so a fan point on the air side of both concave faces lies in the concave strip's region. The ring sector is applied after every strip.

How the reference kernels do this. OpenCASCADE's `ChFi3d_Builder::PerformIntersectionAtEnd` (V7_8_0) ends one fillet at a vertex with more than three sharp edges by intersecting it with the faces there, as the field cuts its strip by the apex walls. At three filleted edges with one inverted concavity, `ChFi3d_FilBuilder::PerformThreeCorner` takes the inverted fillet as the pivot and rolls a constant-radius ball (`BRepBlend_ConstRad`) on the face at the end and the pivot's fillet surface along a corner spine; the field's ring sector is that rolling contact, in closed form. SMLib closes a mixed corner with an N-sided patch that blends the adjacent fillets (`SmFilletConcaveNxNCorner` in `SmFilletCorner.h`), and its single-fillet end (`SmFilletNx1Corner`) is not implemented past three faces. The field follows the recipe's order and OpenCASCADE's contact, not SMLib's patch.

Query cost. A straight piece along the diagonal has an algebraic frame (square roots of 2 and 6), and each distance query enclosed its whole interval cover again, about 23,000 units against the vertex check's 20,000 cap per query. The cover depends only on the path, so a straight piece now encloses it once when it becomes a piece (`Generator::enclose_cover`, called from `Piece::new`). A rigid motion drops the cover, and the moved centre and feet, written through the motion, cost more to enclose on every query (a translated `vertex-4-one` exhausted the 20,000 cap on all 18 vertices even with its cover enclosed again). So a placement (`Plan::transformed`) or a pattern or mirror copy (`Plan::copy`) derives a straight piece again from its moved supports, which are the same points, and encloses its cover under the caller's decider (`Generator::rederive`). `a_placed_straight_piece_keeps_its_enclosed_cover` checks `vertex-4-one` moved by a rational translation at the same per-query cap. Ranges are bit-identical; a query on `vertex-4-one` costs about 16,000 units, and the four corpus rows that already had a straight piece (`acute-5deg`, `acute-halfdeg`, `obtuse-174deg`, `obtuse-179deg`) pay about 3,700 to 4,000 more units once at compile time, with unchanged answers.

Results (`rolling_ends.rs`). The field builds both rows and agrees with every corpus witness (confirmed by OpenCASCADE) and with points derived here, each checked against the B-rep first: in `l-mixed-all` each ring sector's ball centre, a point of the concave blend a sector removes and the same point beyond its reach; in `vertex-4-one` two points past the apex station, one between the edge and the blend (removed) and one under the blend (kept). `vertex-4-one`'s witness on the rounded edge itself lies on both of the piece's supports and is kept as a shared point; the field leaves it `Indeterminate`. Every display vertex of `l-mixed-all`'s two tori named `corner.vertex` lies within the field's guaranteed range of a ring sector's face (338 vertices, worst 3.0e-15 mm), and every display vertex of `vertex-4-one`'s blend lies within range of the straight piece's face (18 vertices, worst 9.1e-16 mm), with the same records at 1, 2 and 4 workers. The strips of `l-mixed-all` are elementary regions with registered faces, which compare already reaches.

Seam control (reviewer note on #738): on the slot loop's junction and the four junctions of `l-mixed-all`, a point on the junction plane and its neighbours 1/1000 either side take the same answer, inside the ball and in the material the pieces remove, where the plane is no face; no face of the fillet lies in a junction plane.

The frozen rows the B-rep refused when the list was frozen keep the field's refusal, now recorded with its reason: `vertex-4-two` (`field.fillet-contact`: two strips meet at the apex, the elementary section along a diagonal needs irrational profile coordinates, and runs that meet at a point with no tangent junction need a corner the field has not built) and `vertex-4-all` (`field.fillet-contact`: the two prisms' base faces lie in one plane but are separate supports with different records, so each base edge counts twice at the base corners and five runs meet there). Every other fillet corpus row and every part keeps its answer and its field work; a scratch survey compared both at the base and at this change.

## Step 5, part 1: rolling over an earlier straight round

Code: the chain route in `crates/catallax-field/src/fillet_chain.rs`, the clearance in `fillet_rolling.rs` (`Plan::outside_solid`), `Feature::face_above` in `features.rs`, and `Plan::refuse_over_rolling` in `fillet.rs`, called by the fillet and chamfer steps. Tests: `crates/catallax-field/tests/rolling_over.rs`, and the field assertion in `crates/catallax-features/tests/fillet_over_fillet.rs`.

Decisions 4.6 step 5: "Fillet over fillet, using the nested path." The research note gives the path: an earlier round is a pipe of radius `R` about its spine; its offset by the ball radius `r` is the coaxial pipe of radius `R - r` or `R + r`, and the new spine is where that pipe meets the other face's offset. This part takes the case where the earlier round is a straight elementary strip, the other face a plane normal to the strip's centre line, and the ball on the centre line's side of the round: a convex chain over a convex round or a concave chain over a concave one. The spine is then a circle of radius `C = R - r` about the centre line, a ring (`C > r`), horn (`C = r`) or spindle (`C < r`) sector, which `Revolved::with_sector` already holds. The field derives it from the strip's cylinder face and the plane, with the ball side from its own material probe; it reads no B-rep fillet module and no fillet definition.

The chain route of step 4 part 2 now takes a body an earlier straight strip rounded. Version 4 continues a named run across the strip's face, which is tangent to the faces it rounded, so the chain rolls over the round and on along the next face. Three things were added for it:

- Free ends. A chain over an earlier strip may end a straight run at a source plane normal to it through its end, as an elementary strip ends at its end plane: `unequal-ordered` runs from the y_max face over the front round and down to the bottom face. The plane must pass through the end and its face must hold the end point.
- Clearance at the round. The round's tangent neighbours, such as the block's x_max and y_min faces beside a corner round, pass through the sector's box, and their leaf is the block, so no face extent clears them. The sector's region lies inside the round's solid cylinder (its wall `rho^2 <= R^2` when the ball lies inside), so a plane parallel to the axis at least `R` from it meets the region at most on its boundary, where a tangent plane touches the cylinder. A plane that supports a convex set meets it only on its boundary. A straight piece next to the round clears the strip's face by its leaf, as step 4 cleared a slot's sides: a straight strip's face is its section's `blend` arc swept along the edge (`Feature::face_above`).
- Rational probes. A fillet face's own classifier takes rational points only, and the midpoint of an arc on a round's face is algebraic. An arc on an earlier strip's face is probed at a rational point of the arc: a line through a rational cut of the circle in a rational direction of its plane meets the circle again at a rational point, the rational parametrization of a conic from one of its points.

Each support pair's runs are now found once per chain, which lowers the work of five refused rows (`at-face-limit`, `past-face-limit`, `swallowed-face`, `r12-cube-at-contact`, `r12-cube-past-contact`) by 16 units each. The other changes run only on a chain over an earlier strip. A scratch survey of every fillet corpus row, every part and every `fillet_over_fillet.rs` recipe at the base and at this change shows no other change of answer, membership grid or work.

A rolling piece registers no face, so a later fillet or chamfer on its body could neither roll over it nor see it in its clearance proof. This part made that reachable: the third round of `fillet_over_fillet.rs`, which the B-rep continues over the second round's chain, came out as a plain strip. A fillet or chamfer on a body that holds a rolling piece now refuses with `field.fillet-contact` (`Plan::refuse_over_rolling`) until the face registry exists. No corpus row or part runs a fillet after a rolling piece.

How the reference kernels do this. Both treat the earlier round as an ordinary face. OpenCASCADE's `ChFi3d_Builder::PerformElement` (V7_8_0) collects the tangent chain, and a plane with a cylinder on a circular spine is built in closed form as a torus (`ChFiKPart_ComputeData_FilPlnCyl`, read for step 2). SMLib's `SmFilletExecutive::CollectTangentEdgeSequence` collects smoothly meeting filleted edges, and its solvers offset both faces (`SmOffsetSurface`) and intersect the offsets by stepping. The field takes the same offset definition, the nested path, in closed form.

Results (`rolling_over.rs`). The field builds the corpus rows `tangent-boss` (four horn sectors, `C = r = 1`, and four strips), `tangent-pocket` (the same, concave) and `unequal-ordered` (a horn sector between two strips that end at the y_max and bottom faces), and from `fillet_over_fillet.rs` the boss at radii 3/1 (ring sectors) and 3/2 (spindle), the 3/1 loop named edge by edge in version 3, the pocket at 3/1, the corner at 3/1 from the top edge or from the round's own arc and at 3/2, and the boss and corner at 3/1 in a rational pose turned a third of a turn about (1, 1, 1). It agrees with every witness: `unequal-ordered`'s from OpenCASCADE, the test file's, and points derived here near each strip and sector, each checked against the B-rep first. `box-mitre`, which the elementary route already built (its two strips meet in the mitre plane), agrees with points derived at its corner. Every display vertex of the B-rep's blend faces of the rolling step lies within the field's guaranteed range of a piece's face at 1 nm: 5,487 vertices over 12 cases, the worst 5.5e-15 mm, with the same records at 1, 2 and 4 workers. A query costs at most 5,016 units, and 7,276 in the rational pose, against the check's cap of 20,000; each case's pieces are compiled once and reused by every query. Compile work is 6,688 units for `tangent-boss`, 9,668 for `tangent-pocket`, 1,588 for `unequal-ordered` and 14,692 for the boss at 3/1. The pocket at 3/1 costs 543,250: each corner circle crosses the other three upright cylinders at algebraic points, so it is cut eight times and six arcs inside the material are probed too.

Seam control: at the eight junction points where a strip meets a sector over an earlier round, a point on the junction plane and its neighbours 1/1000 either side take the same answer, in the material the chain removes or adds and inside its ball.

Refusals. The B-rep builds every other recipe of `fillet_over_fillet.rs` and the field still refuses each, recorded with its code in `the_field_names_what_it_has_not_rolled_over_yet` for the later parts of step 5: a ball as wide as the round (`field.rolling-path-reach`, the ball no smaller than the cylinder it rolls inside; the B-rep's sphere corner), a ball wider than the round (the cliff of RECIPE 4.5.1, the same code), a ball outside the round's cylinder (`mixed-ordered`, the bracket and the boss base: the concave round's end cap has no registered face and the round's tangent faces cross the sector), an upright under an earlier mitre (`roll-over-mitres` and the mitre recipes) and every edge at once (`fillet.unsupported-overflow`), and a third round over a chain of rolling pieces (`field.fillet-contact`, above). Where both kernels refuse (a ball wider than a round named by its arc, and version 3 strips that do not roll over), the field's code is recorded too.

## Step 5, part 2: a ball as wide as the round, or wider

Code: `Revolved::new`, its walls and separation in `crates/catallax-field/src/rolling_region/revolved.rs`; the cliff in `fillet_chain.rs` (`Plan::narrower_round`, `Plan::cliff_end`, `Run::cut_by`); `Generator::plane_prism`, `Plan::beside_prism` and `Plan::beside_box` in `fillet_rolling.rs`. Tests: `crates/catallax-field/tests/rolling_over.rs`.

A sphere corner. RECIPE 4.5.1: "A round exactly as wide as the ball is not a cliff: the ball's centre stays on the round's axis." The round's cylinder offset by the ball toward its axis is the axis itself (`C = R - r = 0`), so the plane's offset meets it in one point and every ball of the sector is one ball. The ring sector of part 1 takes this as its limit, for a plane normal to the axis: the separation is the ball's own `|x - c|^2 - r^2`, and the fan is the corner cell between the plane, the plane through the centre and the round's cylinder, within the sector. The meridian wall degenerates on the axis, where the fan's points lie inside the one closed ball; only the plane contact there is a boundary point. OpenCASCADE's `ChFiKPart_ComputeData_FilPlnCyl` (V7_8_0) makes the same sphere when `cylrad - Radius` is zero (`c1sphere`). SMLib's `SmFilletStandardSolver.cpp` steps the offset intersection numerically and has no closed form for a centre that stays fixed; the field follows OpenCASCADE.

A cliff. RECIPE 4.5.1: where a foot reaches a round of the strip's convexity narrower than its ball, "The run takes the cliff rule instead. Both feet stay on their original supports, and the run continues along its spine's extension", and "A tangent chain stops at such a round instead of continuing over it." The chain no longer continues onto such a round (`narrower_round`: the round's support passes through the run end, the ball lies inside its cylinder, and `R < r`). A straight run between planes then runs on to the nearest plane normal to it, past the cliff, that holds the other rail's point there (`cliff_end`), as an elementary strip ends at its end plane. This follows the B-rep's cliff (`fillet_build.rs`, `cliff`), which walks its cut from the cliff foot across the faces past the round until it reaches the strip's other face at the other rail's end. The round's face must lie past the plane normal to the run at the cliff, proven from its leaf (`Feature::face_above`), and by the cliff rule it cuts the extension, so the clearance leaves it to that cut (`Run::cut_by`). Two cliffs at one end, as under an earlier mitre, keep the refusal. OpenCASCADE has no cliff: `ChFiKPart_ComputeData_FilPlnCyl` stops with "the fillet can't pass" when `cylrad - Radius` is negative.

Two clearance proofs were added, each tried only where the box test would refuse, so no earlier answer or work changes. A straight piece's one box is its axis-aligned hull; in a rotated pose a wide strip's hull reaches the part's far faces (the boss at radii 3/3 placed as in part 1 reached its bottom face). A plane face with every vertex of the piece's section prism (edge point, feet and centre at both end stations) on one closed side meets the region at most on its boundary (`beside_prism`). And a curved face evaluated over a long box can seem to reach it by interval dependency alone: the strip of the boss at 1/3 runs past its cliffs along the whole edge, and the expanded polynomials of the far corner rounds straddle zero over its box. A face whose leaf puts it beyond one side of the box meets the region at most in that plane (`beside_box`).

Results (`rolling_over.rs`). The field builds, from `fillet_over_fillet.rs`, the sphere corners of the boss at radii 3/3 and 2/2 (four spheres and four strips each), the pocket at 3/3 and 2/2 (concave), the corner at 3/3, 2/2 and 2/2 from the round's arc, the boss at 3/3 in the rational pose of part 1, and the cliffs of the boss at 1/3, the pocket at 1/3 and the corner at 1/2 (one straight piece each). Each agrees with points derived here, checked against the B-rep first, and with the test file's points: kept inside each ball, cut or added outside it in each strip and corner cell, on the round's axis just short of the plane, on a cliff's extension where the round's material (or, in the pocket, its air) meets the strip, and edges the chain does not reach left sharp. Every display vertex of the B-rep's blend faces of the rolling step lies within the field's guaranteed range of a piece's face at 1 nm, with the same records at 1, 2 and 4 workers: 6,448 vertices over the eight sphere cases, the worst 6.1e-15 mm, at most 4,674 units a query and 6,934 in the rational pose; and 152 over the three cliffs, the worst 3.7e-15 mm, at 2 units a query from the straight piece's closed form. The sphere faces hold rational points `c + (1, 2, 2)`: the B-rep classifies them as boundary, and in the field exactly one piece, the sector, holds each on its region's boundary with a distance range from 0. The plan itself leaves them `Indeterminate`, as every point on a rolling piece's face: a point on a subtracted or added region's boundary needs that region's face registered (below).

Compile work: 14,692 units for the boss at 3/3, 6,688 at 2/2, 33,364 for the boss at 3/3 in the rational pose, 3,961 for the corner at 3/3, and 543,215 for the pocket at 3/3, whose corner circles the other uprights cut at algebraic points as at 3/1; the cliffs take 842 (boss), 1,565 (pocket) and 298 (corner).

A scratch survey of every fillet corpus row, every part and every `fillet_over_fillet.rs` recipe at the base and at this change shows no other change of answer or membership grid. Four refused recipes pay less, because the chain now stops at a narrower round instead of following it: the corpus row `unequal-mitres`, which both kernels refuse (11,185 to 1,286 units), the mitred loop at 1/2, the upright at 3/4 under two rounds and the mitred pocket at 1/(3/2). Version 3's boss at 3/1, which both kernels refuse, pays 5 units more (7,262 to 7,267) for the new clearance proofs it tries before refusing.

Refusals. Still refused, each recorded with its code in `the_field_names_what_it_has_not_rolled_over_yet`: a ball outside the round's cylinder (`mixed-ordered`, the bracket and the boss base), an upright under an earlier mitre, with a cliff at both feet or rolling over the mitre (`roll-over-mitres` and the mitre recipes), every edge at once (`fillet.unsupported-overflow`), and a third round over a chain of rolling pieces (`field.fillet-contact`). A ball wider than a round named by its own arc stays refused in both kernels.

## Step 5, part 3: an upright under an earlier mitre

Code: `Plan::mitre_corner`, `Plan::apply_corner` and `Plan::cliffs_end` in `crates/catallax-field/src/fillet_chain.rs`; `Plan::trihedral` in `fillet.rs`, which the elementary corner now calls too. Tests: `crates/catallax-field/tests/rolling_over.rs`, and the field assertion in `crates/catallax-features/tests/fillet_over_fillet.rs`.

An earlier mitre is two straight rounds meeting in the plane that bisects their faces' corner, as a block's top loop rounded at once leaves it. An upright whose straight run reaches such a corner reaches both rounds at once, one at each foot.

A ball as wide as both rounds. RECIPE 4.5.1: "A round exactly as wide as the ball is not a cliff: the ball's centre stays on the round's axis." Both rounds' axes pass through the upright's ball centre where its feet reach them, so the ball stops there, touching the run's two planes, both rounds and the plane the rounds roll on. The material-set rule of RECIPE 4.5.1 takes "the part of the input material inside the translated positive cone of those active normals and outside the common ball". The active normals are the two planes' and the rounds' contact normals, which turn from each plane's normal to the third plane's, so the cone is the one spanned by the three planes' normals at the centre, and those three planes bound the material in it. That is the trihedral corner of the three planes, as if the three edges had been rounded together, and it swallows the mitre. The field detects the end where the other routes leave it open (`mitre_corner`): two rounds of the run's convexity, as wide as the ball, each with its axis parallel to its foot's plane, normal to the run and through the ball's centre, and a source plane normal to the run one radius past the centre; the three directions must be orthogonal unit vectors, as at a block's corner. The run's straight piece ends at the plane normal to it through the centre, a junction, and the corner is the cell of the elementary corner, the cone of the three directions short of the three planes, less the ball (`trihedral`). Its sphere face is registered under a corner name of the step, so the plan answers `Boundary` on it. RECIPE 4.5.1 names a junction sector `corner.vertex` "for their original input vertex", the vertex of the step's input body where the run ends: its two planes and the two rounds meet there, so those four names make the corner's name, as the B-rep's does. Every source face other than the three planes and the two rounds must stay clear of the corner's cube, or the fillet refuses. The B-rep builds the same corner with a `CornerDefinition` of the three strips (`fillet_build.rs`); its test file's volume, "the corner cube less its ball octant", says the same. OpenCASCADE refuses the corpus row (`occt.fillet-not-done`). SMLib's corner code (`SmFilletCorner.cpp`, called from `SmFilletExecutive`) types a filleted vertex by how many of its edges one operation rounds (its Nx1 and Nx2 corners) and solves a blend patch there; the material-set rule makes this corner the trihedral cell of three planes, so the field needs no patch solver.

A ball wider than both rounds. RECIPE 4.5.1's cliff rule at both feet: the strip keeps both planes, runs on past the cliff, and both rounds cut it. Neither rail goes on, so part 2's end (the plane that holds the other rail's point) does not exist. The run ends at the nearest source plane normal to it past the cliff to which both rounds are tangent, the plane they roll on (`cliffs_end`); the rounds' faces end there, so every cut lies short of it. The B-rep bounds its cliff region the same way, by a plane normal to the edge past every cut curve (`Cliff::cut` in `fillet_build.rs`), and its test file's upright "runs on to the top and the rounds cut it". Each round's face must lie past the plane normal to the run at the cliff, proven from its leaf, as in part 2.

Results (`rolling_over.rs`). The field builds the corpus rows `roll-over-mitres` (a top loop and its uprights at radius 1) and `unequal-mitres` (the uprights at radius 2, which the B-rep refuses inside the corpus allowance and builds without one), and from `fillet_over_fillet.rs` the mitred corner at 3/3 and 2/2, the mitred loop at 3/3, the mitred pocket at 1/1 (concave), the mitred corner at 3/3 in the rational pose of part 1, and the cliffs at both feet of the mitred corner at 3/4 and the mitred pocket at 1/(3/2). Each agrees with the corpus witnesses, the test file's points and points derived here, each checked against the B-rep first: kept inside the strip's and the corner's balls, cut or added outside them, where the mitre stood (inside both rounds, outside the corner ball), under a narrower round outside the upright's disk, and inside both disk and round where the mitre's outer part survives. Every display vertex of the B-rep's blend faces of the rolling step lies within the field's guaranteed range of a piece's face at 1 nm, with the same records at 1, 2 and 4 workers: 630 vertices over the six swallowed-mitre cases, the worst 4.3e-15 mm, and 368 over the three with cliffs at both feet, the worst 3.5e-15 mm, each at 2 units a query from the straight piece's closed form. Every display vertex of the B-rep's corner faces, 3,306 over the six cases, lies within 1 nm of the field's registered corner sphere, at most 1 nm outside its cone, with the centres and radii derived here; the rational points `c + (1, 2, 2)` of the 3/3 corner are boundary points in both kernels.

Refusals. Still refused, each recorded with its code in `the_field_names_what_it_has_not_rolled_over_yet`: a ball outside the round's cylinder (`mixed-ordered`, the bracket and the boss base), a ball narrower than a mitre's rounds, which rolls over the mitre along the ellipse where the rounds' offsets meet (the mitred loop at 1/(1/2), the corner at 3/2, the pocket at 1/(1/2)), every edge at once (`fillet.unsupported-overflow`), and a third round over a chain of rolling pieces (`field.fillet-contact`).

## Step 5, part 4: a ball outside the round's cylinder

Code: `Plan::buried_rails`, `Plan::arc_material`, `Plan::strip_cap` and `Plan::same_face` in `crates/catallax-field/src/fillet_chain.rs`, and the coplanar continuation in its `Plan::continue_tangent`; `Feature::on_cap` in `features.rs`; `Plan::on_own_plane` in `fillet_rolling.rs`; `Plan::root_location` and `Plan::material_unless` in `fillet.rs`. Tests: `crates/catallax-field/tests/rolling_over.rs` (`outside_rounds`), and the field assertion in `crates/catallax-features/tests/fillet_over_fillet.rs`.

`mixed-ordered` rounds an L's inside corner, then the base's right top edge, which rolls over the inner round's end on the L's x_max plane and up the wall. The chain is convex and the round concave, so the ball lies outside the round's cylinder: its centre stays `r` from the plane and `R + r` from the round's axis, a ring sector about that axis whose centre circle has radius `R + r`. The test file's bracket does the same at other radii, and its boss base rolls a concave chain round a boss's foot over its convex upright rounds. RECIPE 4.5.1: "The result equals one fillet operation on the concave edges followed by another on the remaining selected convex edges and those arcs." The piece is step 2's plane with a cylinder, with the ball outside (`Revolved` with `outside`), so only the chain's rules and the clearance proof are new. OpenCASCADE's `ChFiKPart_ComputeData_FilPlnCyl` (V7_8_0) builds the same torus with spine radius `cylrad + Radius` when the ball is not inside the cylinder (`dedans` false). SMLib offsets both faces (`SmOffsetSurface`) and intersects them by stepping; the field takes the offset definition in closed form.

Three things stopped it.

- The arc on the plane. The inner round ends at the L's x_max plane, so its cap lies in that plane and holds the plane's face along the arc where the chain rolls onto the round; the B-rep merges that cap into the x_max face. The field's base x_max face ends at the base's top, so the arc's material probe found no face there. A strip registers only its blend face. Where a plane normal to an earlier strip's line passes through the arc point and the point lies on that strip's cap (`Feature::on_cap`: at its first or last station, inside or on its profile), the cap's test replaces the plane face's own (`arc_material`). The four material rays are unchanged.
- Coplanar faces. Up the wall the chain keeps the x_max plane, now as the wall's x_max face, a face of another leaf with the same support. The B-rep merges such faces into one. Where the kept face itself has no run through a chain end, the chain tries the coplanar faces of other leaves (`continue_tangent`), and a junction may share a face by its support (`same_face`).
- The round's two faces. The faces the earlier round rounded pass through the sector's box, and their leaves bound the body there, so no box, leaf or solid argument clears them. The fan is shown to be all material for a convex sector, or all air for a concave one, so removing or adding the fan less its balls is the rolling ball's result whatever those faces do inside it (`buried_rails`). The fan `Q` is the set of points between the plane and its offset by `r`, `R <= rho <= R + r` from the round's axis, between the meridians through the round's two rails in the plane. Every other face is cleared from the box as for any run, and the round runs at least `r` past the plane on the ball's side (its face holds a point there), so inside `Q` the body can change only across the round's two faces: the plane and the cylinder bound `Q`, the round's walls lie on its two faces and its caps lie outside the open slab. Those faces are tangent to the cylinder at the rails, so they cut `Q` into at most four cells: the round's own section between the arc and the corner, the cell beyond each face alone, and the cell beyond both, which is empty when the round's corner lies at least `R + r` from the axis. In cylindrical coordinates about the axis each cell is a slab interval times a planar cell; beyond either face the inner disc does not reach, so that planar cell is convex, and the round's own section is an interval of `rho` above `R` at each angle, so every cell is connected and the body is constant on it. One rational probe strictly inside each nonempty cell, classified exactly in the step's input body, settles it: dyadic points rounded from targets at fixed fractions of the angle and of the radial gap, each verified inside its cell by exact signs. This is the mixed vertex's argument (step 4, part 3) with the concave round already in the body instead of formed in the same step.

A boss's bottom face lies in the plate's top plane, a face of another leaf on the sector's plane. The region lies on the ball's closed side of each support it rolls on, so such a plane meets it at most on its boundary (`on_own_plane`, tried only where the box test would refuse).

Names. RECIPE 4.5.1: "Every exposed edge-run patch keeps `blend.edge` for the originally selected input edge, including after roll-over." A continued span now takes its name from the selected span it continues, as the B-rep names it: the three patches of `mixed-ordered` and of the bracket, the eight of the boss base, and the continued patches of `unequal-ordered`, `tangent-boss`, `tangent-pocket` and the test file's recipes of parts 1 and 2, which were named by their own faces before. Membership and distances do not change.

Results (`rolling_over.rs`). The field builds `mixed-ordered` (three pieces: the base's strip, the sector over the inner round, the wall's strip), the test file's bracket at radii 3/1, 3/4 and 2/(3/2), and the boss base at 2/1, 2/3 and 3/1 (four sectors and four strips each). Each agrees with points derived here, checked against the B-rep first: cut or added in each strip and sector outside its ball, kept or air inside it, the inner round's material past the sector's slab, and the air inside the inner round. Every display vertex of the B-rep's blend faces of the rolling step lies within the field's guaranteed range of a piece's face at 1 nm, with the same records at 1, 2 and 4 workers: 4,122 vertices over the seven cases (203 for `mixed-ordered`), the worst 4.1e-15 mm, at most 5,016 units a query. At the two junctions of `mixed-ordered`, a point on the junction plane and its neighbours 1/1000 either side take the same answer, cut outside the ball and kept inside it. Compile work is 9,348 units for `mixed-ordered`, 30,219 for the bracket at 3/1 and at 3/4, 9,628 at 2/(3/2), and about 41,000 for the boss base at 2/1 and 2/3; at 3/1 it is 350,644, because each corner circle crosses the neighbouring upright rounds at algebraic points, as in the pocket at 3/1.

Compare agrees on `mixed-ordered`: 1,367 forward vertices, 194 reverse samples and no finding, in 0.6M work. `roll-over-mitres` and `box-mitre` still agree with the new names (5,015 and 212 forward vertices, no finding).

Refusals. Still refused, each recorded with its code in `the_field_names_what_it_has_not_rolled_over_yet`: a ball narrower than a mitre's rounds, which rolls over the mitre along the ellipse where the rounds' offsets meet (the mitred loop at 1/(1/2), the corner at 3/2, the pocket at 1/(1/2)), every edge at once (`fillet.unsupported-overflow`), and a third round over a chain of rolling pieces (`field.fillet-contact`).

Reviewer note on #744: the walk behind `Plan::refuse_over_rolling` followed both inputs of every Boolean, and Booleans share their inputs' history, so a shared history was scanned many times outside the work meter. It now marks each live instruction once in one backward sweep, as the copy and compaction walks do (`Plan::rolling_below`); `the_rolling_walk_visits_shared_history_once` counts 41 visits for 40 levels of doubly shared Unions, and the third round over a chain of rolling pieces still refuses with `field.fillet-contact`.

## Step 5, part 5: a ball narrower than a mitre's rounds

Code: `crates/catallax-field/src/rolling_region/seam.rs` (new), `Plan::mitre_seam` and `Plan::apply_seam` in `crates/catallax-field/src/fillet_chain.rs`, and a run's solid cylinders in `fillet_rolling.rs` (`Run::inside`, `Plan::outside_solid`). Tests: `crates/catallax-field/tests/rolling_over.rs` (`mitre_seams`, the corpus row `narrow-mitres`, `mitre_seam_vertices_lie_on_the_field_seam_face`), and the field assertion in `crates/catallax-features/tests/fillet_over_fillet.rs`.

An upright whose ball is narrower than both rounds of an earlier mitre reaches them together, as in part 3, but they are not cliffs: the ball can roll onto them. RECIPE 4.5.1: "The ball rolls between the two rounds to the mitre's apex, with its centre on the ellipse where their offsets meet, and sweeps one exact rational band." The rounds' axes run along the run's plane normals `n_b` and `n_a`, each `R` from its plane, and meet at `q = p + R (n_a + n_b)` for the run end `p`. The ball lies inside both cylinders, so its centre lies `k = R - r` from both axes. The two offsets are equal cylinders on meeting axes, step 3's crossing, and meet in the planes that bisect the axes; the branch through the upright's last centre lies in the mitre plane:

```text
c(t) = q - k cos t (n_a + n_b) + k sin t d,      t in [0, pi/2],
```

from that centre to the apex, whose ball touches the plane the rounds roll on where their rails meet. With `tan(t/2) = (u + 1) / 2` every station is rational. The contact rays `-cos t n_a + sin t d` and `-cos t n_b + sin t d` are radial from each axis and square to the tangent, so the face is the arc between them in the normal plane.

The region. In the corner's frame `x - q = a n_a + b n_b + h d`, a point `c(t) + s_1 m_1 + s_2 m_2` of the fan has `a = -(k + s_1) cos t`, `b = -(k + s_2) cos t` and `h = (k + s_1 + s_2) sin t`. In the open octant `a, b < 0 < h`, eliminating `s_1` and `s_2` leaves `-(a + b) / cos t - h / sin t = k`, whose left side rises from minus to plus infinity on `(0, pi/2)`. So one station holds each point's section, and by that rise `s_1 >= 0` and `s_2 >= 0` are the polynomial walls `a^2 (b^2 + h^2) >= k^2 b^2` and `b^2 (a^2 + h^2) >= k^2 a^2`. With `a <= 0`, `b <= 0` and the ball side of both rounds' cylinders they bound the fan exactly. Two bands, outside both offsets, keep out the segment of the line through `q` along `d` below `q + k d`, where both polynomial walls, `a` and `b` vanish. `h >= 0` is a junction with the upright's straight piece, whose section is the band's at `t = 0`. The region lies inside both rounds' solid cylinders, so the run's two planes and the plane the rounds roll on, each tangent to one of them, meet it at most on its boundary: part 1's clearance, now for each of a run's solid cylinders. No closed form gives the station (a quartic in `tan(t/2)`), so membership near the face and distance take step 1's panel cover with step 3's binary64 guess, as the saddle does. The band lies on no elementary surface and registers no face, so the plan leaves points on it `Indeterminate`.

The route (`mitre_seam`) tries only run ends that every earlier route leaves open, so every earlier answer and its work stand. It needs a straight run between two planes with orthogonal unit normals, a round of the run's convexity, wider than the ball, tangent to each foot's plane along the line through the end normal to the run, both of one radius, and a source plane normal to the run `R` past its end. The other open ends may then end at a plane normal to the run, as beside a mitre corner. RECIPE 4.5.1 names a junction sector `corner.vertex` after the operation's input vertex; the band is named after the vertex where the run's planes meet the two rounds, as the B-rep names it. Rounds of different radii stay refused, as in the B-rep, and so do planes that are not orthogonal, which the B-rep builds.

How the reference kernels do this. OpenCASCADE builds the corpus row `narrow-mitres`, whose witnesses are its own; between two cylinders it walks `BlendFunc_ConstRad` with `BRepBlend_Walking` (step 3). SMLib intersects the offset surfaces by stepping (`SmOffsetSurface`, step 2). The B-rep builds the band as one exact rational patch over the same ellipse. The field takes the offsets' meeting in closed form, as for step 3's crossing, and replaces walking with exact walls and the panel cover.

Results (`rolling_over.rs`). The field builds the corpus row `narrow-mitres` (four uprights and four bands), which OpenCASCADE builds and the B-rep refuses inside the corpus allowance, and from `fillet_over_fillet.rs` the mitred corner at 3/2 (one upright, one band) and the mitred pocket at 1/(1/2) (concave, four and four). Each agrees with the corpus witnesses and points derived here, each checked against the B-rep first: kept (or air) inside a band's ball in the mitre plane, cut (or added) between that ball and the mitre, inside both rounds and outside every ball of the band, the strip's cut and kept points below the junction, and material the step leaves alone. Every display vertex of the B-rep's strip faces lies within the field's guaranteed range of a piece's face at 1 nm (174 vertices over the three cases, the worst 2.8e-15 mm, 2 units a query), with the same records at 1, 2 and 4 workers. The B-rep's band is a corner face of the step under the vertex's name, and the field's seam pieces carry the same names, one per band. Every display vertex of the B-rep's bands lies within the field's range of a seam piece's face at 1 nm: 859 vertices over the three cases, the worst 2.6e-15 mm, with the same records at 1, 2 and 4 workers. A vertex near the band it lies on settles at 4 units, from the exact face point at the binary64 station guess, which bisection on the slope of the squared distance finishes, since a golden-section search alone stops near the square root of the roundoff. A vertex of another corner's band, queried against this one, costs up to 374 units before its range closes to 1 nm.

The B-rep refuses the mitred corner at 3/2 in the rational pose of part 1 (`intersect.unsupported-pair` while it cuts its band), so there the field is checked against itself: each derived point, placed, takes its unplaced answer, and the unplaced B-rep's band vertices, placed exactly, lie within the placed seam piece's range at 1 nm (235 vertices, the worst 5.3e-15 mm).

Refusals. Still refused, each recorded with its code in `the_field_names_what_it_has_not_rolled_over_yet`: every edge at once (`fillet.unsupported-overflow`, `fillet.unsupported-contact`), and a third round over a chain of rolling pieces (`field.fillet-contact`).

## Step 5, part 6: every edge at once, and a third round over rolling pieces

Code: `Plan::chain_fillet`'s `"all"` selection, `Plan::block_corner` and `Plan::swallows` in `crates/catallax-field/src/fillet_chain.rs`, and `Plan::beside_bounds` in `fillet_rolling.rs`; for the third round, `Plan::earlier_round` and `Plan::hidden_pieces` in `fillet_chain.rs`, the route in `Plan::fillet` (`fillet.rs`), `Piece::end_ray` and `Piece::above` in `rolling_region/piece.rs`, `Generator::end_plane_at`, and the rolling arm of `Plan::face_above` in `fillet_rolling.rs`. Tests: `crates/catallax-field/tests/rolling_over.rs` (`every_edges`, the case `corner bottom 3/2/1`, vertex groups 18 to 20, and the block corners in `mitre_corner_vertices_lie_on_the_field_corner_face`).

Every edge at once. RECIPE 4.5.1: "Every chosen edge uses the same `r`, and all edges are blended at once", and "an `"all"` selection leaves smooth edges out". After an earlier round the elementary route refuses an `"all"` selection, because a chain must roll over the round's ends. The chain route now takes it where at least one source face is an earlier straight strip, with every sharp run of every support pair; every other `"all"` selection keeps the elementary route's answer, and one the chain cannot close keeps that route's reason. Tangent faces give no run, so the round's rails drop out. Three things were new.

- Block corners. Three runs that end at one point, each between two of three planes with pairwise orthogonal unit normals and all of one convexity, close with the trihedral corner of the three planes less the ball, as the elementary route closes a box's corner: at the vertex, RECIPE 4.5.1's material-set rule takes the input material inside the cone of the three active normals at the ball's centre, one radius from each plane. Each run's ball touches the third plane where its centre lies one radius short of the vertex, so each run ends there at a junction with the corner. The corner is named after its three planes, as the B-rep names it. Any other meeting of three or more runs keeps its refusal. OpenCASCADE closes three fillets of one radius at a vertex with a sphere (`ChFi3d_FilBuilder::PerformThreeCorner` through `ChFiKPart_ComputeData::ComputeCorner`, V7_8_0); SMLib fills a three-sided corner with a sphere through its boundary arcs where one fits (`SmFilletCorner.cpp`, `CreateSphereFromArcs`).
- The seam under a mitre corner. With every edge at radius 2 after a mitred top loop of radius 2, each upright swallows the mitre above it (part 3), and the selection also holds the four seams where the loop's rounds meet, ellipses the chain does not section. A ball one radius from both rounds' cylinders, on their material side, has its centre on both axes, which is the corner's centre, so a seam has no round of its own; the B-rep's test says "the mitre above it goes without a round of its own". Each round rounds one of the run's planes and the plane both roll on, and its face lies between its rails, so a point on both faces lies in the corner's cube and at least one radius from its centre: the seam lies in the cell the corner removes. A pair the chain cannot section is accepted only as the two rounds of a corner this step builds, each rounding those two planes.
- A far round over a corner's cube. A round's expanded cylinder equation, evaluated over a box far from its axis origin, can seem to reach the corner by interval dependency alone (the loop's y_max round over the block's corner at (0, 20, 0)). As for a run's box, the face's leaf then shows the whole face beyond one of the cube's sides.

Results (`rolling_over.rs`). The field builds the four recipes of `every_edge_rounds_over_an_earlier_round` in `fillet_over_fillet.rs`: after a front round at 3/1 and at 2/2 (thirteen pieces each: eleven strips and two sectors over the round's ends, one ball each at 2/2; six block corners), after rounded uprights (sixteen pieces: two loops of four strips and four ring sectors) and after a mitred loop (eight strips, four block corners at the bottom and four mitre corners at the top). Each agrees with points derived here, checked against the B-rep first: cut outside each strip's, sector's and corner's ball and kept inside it, the round left alone between corners, and uncut material. Every display vertex of the B-rep's blend faces lies within the field's range of a piece's face at 1 nm (720, 781, 1,416 and 370 vertices, the worst 3.8e-15 mm, at most 4,978 units a query), with the same records at 1, 2 and 4 workers. The B-rep's corner faces of the step are its block corners and, under the loop, its mitre corners; their display vertices lie within 1 nm of the field's corner spheres, inside each corner's cell (738 at 3/1, 1,165 at 2/2 and 1,535 under the loop), but for one at 2/2 and one under the loop.

A B-rep display finding. With every edge at radius 2, one display vertex of the B-rep's corner sphere about (2, 18, 2) lies 0.62 µm past the circle where that sphere meets the bottom strip along x_min: outside the B-rep's own corner face, and 96 nm inside the strip's exact cylinder, which is the field's surface there. The same vertex comes out at 2/2 and under the mitred loop, and the test records it, one per case. It is for the B-rep's owners; no B-rep code changes here.

A third round over rolling pieces. The right chain of `corner 3/2` is three rolling pieces: a strip along the top edge, a spindle sector over the front round's end and an upright strip down the front edge. The third round rounds the bottom right edge, which RECIPE 4.5.1 continues over the foot of that upright round and along the bottom front edge; the B-rep builds it as two strips and a torus sheet. Every route refused a fillet on a body with rolling pieces (`Plan::refuse_over_rolling`, part 1). Now only the chain route takes such a body, with three additions, and every other route and a chamfer still refuse it with the old reason.

- A straight rolling piece as an earlier round. Its registered cylinder is the same quarter cylinder between two rails as an elementary strip's face, so the chain rolls over it as over a strip (`earlier_round`), probing its arcs at rational points.
- Rays at a strip's end. The arc where the third round meets the upright round lies on that round's end plane, the bottom face. Near a point of the face on that plane, the piece's region is the face's positive side within the end plane's side, so a ray's side follows from the first nonzero coefficients of both equations along it (`Piece::end_ray`), as at an elementary strip's cap. Only points the plan left `Indeterminate` change.
- Clearance. A straight piece's region lies in the prism of its two end sections, so its face lies beyond a junction or a box side when every vertex of that prism does (`Piece::above`), as an elementary strip's leaf shows it. A piece that registers no face, such as the spindle sector, stands for no source face, so no clearance test sees it; every run of the chain must keep its box apart from such a piece's box, or the chain refuses.

Results. The field builds the third round, `corner bottom 3/2/1` (three pieces: two strips and a horn sector, `C = r = 1`, over the upright round's foot), and agrees with points derived here, checked against the B-rep first: cut outside each ball and kept inside it on the sector's meridian, the upright round's air and material above the chain, and uncut material. Every display vertex of the B-rep's blend faces of the third round lies within the field's range of a piece's face at 1 nm (147 vertices, the worst 2.5e-15 mm, at most 4,902 units a query), with the same records at 1, 2 and 4 workers.

Refusals. Every step 5 recipe of `fillet_over_fillet.rs` that the B-rep builds now builds in the field. The recipes both kernels refuse keep their codes in `the_field_names_what_it_has_not_rolled_over_yet`.

## Step 6: torus pairs

Code: `Revolution::Torus` and `Shape::Torus` in `crates/catallax-field/src/rolling_region/revolved.rs`; the route in `fillet_rolling.rs` (`admits`, `revolution`). Tests: `crates/catallax-field/tests/rolling_tori.rs`.

Both frozen rows cut a torus (core radius 6, tube 2) with a plane normal to its axis and round the two circles where they meet, with the ball inside the tube. The field derives the ball path as for a sphere in step 2: offset the plane by `r` toward the ball and the torus by `r` toward its core circle (the torus of tube radius `a - r`), and intersect. In the meridian half-plane the torus is the circle of radius `a` about `(R, 0)`, where a sphere's is about `(0, 0)`, so the offset plane at height `h0` meets the offset torus where `(rho - R)^2 + h0^2 = (a - r)^2`: two circles about the axis, `rho = R +- sqrt((a - r)^2 - h0^2)`, each edge taking the one on its side of the core circle. Each edge is a ring piece with the circle's closed forms (separation `(rho - C)^2 + h^2 - r^2`, the meridian as the nearest station). The fan's meridian wall is the offset tube's normal line through the station, which passes through the core circle's point `(R, 0)`; its edge side is fixed by the sign at the edge, which a cut below the torus's centre plane reverses for one branch of the formula. The station's normal to the plane bounds the fan on the other side: the point's meridian distance from the axis lies on the edge's side of the centre circle's, `sign(C - R) (rho - C) >= 0` toward the edge. (Written squared, as `(rho - R)^2 + h0^2 - (a - r)^2`, the wall also held on the far side of the core circle; with the cut below the centre plane the inner circle's meridian wall held there too, so its fan reached the outer edge and took kept material inside the outer ball, the point (7.7, 0, -0.85). The compare witness test found it.) The torus's tube point in the station's meridian gives the foot on the wall. Only a plane normal to the axis is taken; another plane meets the offset torus in a curve that is no circle, and refuses with `field.rolling-path-reach`. An all-edges selection still skips the pair, as it skips step 3's pairs, so its answers do not change.

`torus-intersection` (cut at z = 1, r = 1/2) is the case where the offset plane touches the offset tube's top: both edges share the centre circle `rho = 6` at z = 3/2, and their fans meet in the cylinder `rho = 6` over the core circle. What is left is the tube of radius 1/2 about that circle, the whole torus of major radius 6 and minor 1/2 that the row's reference volume describes.

How the reference kernels do this. OpenCASCADE's closed-form fillets (`ChFiKPart`, V7_8_0) cover a plane with a plane, a cylinder or a cone and sphere corners, not a torus, so it walks `BlendFunc_ConstRad`; it builds `torus-slice` and refuses `torus-intersection` (`occt.fillet-not-done`). SMLib offsets both faces (`SmOffsetSurface` in `SmFilletStandardSolver.cpp`), intersects the offsets to find the centre line and lofts the fillet surface, "in some special cases a cylinder or a torus". The field takes the same offset definition in closed form.

Two more things were needed. A plane-torus circle of irrational radius comes from the support intersection as an intersection curve, not a circle, so the route takes that curve with its witness point as the edge sample (`revolved_runs`, shared with the circle branch, whose calls are unchanged). And where the feet are opposite, as in `torus-intersection`, no short arc runs between them: a path frame now may carry a `side`, a direction from the centre into the face arc's half of the normal plane, and the arc's distance bound uses that half circle where the feet are parallel. Only torus frames set it (toward the edge point), so every other path's answers are unchanged.

Results (`rolling_tori.rs`). The field builds `torus-slice`, `torus-intersection` and a cut at z = -1 below the torus's centre plane (r = 1/5), two ring pieces each, and agrees with every witness: the corpus rows' (built by OpenCASCADE for `torus-slice`) and points derived here, each checked against the B-rep first; `torus-intersection`'s two boundary witnesses lie on a piece's region boundary at distance 0. Every display vertex of the B-rep's blend faces lies within the field's guaranteed range of a piece's face at 1 nm, with the same records at 1, 2 and 4 workers: 2,192, 3,136 and 1,288 vertices, the worst 4.7e-15 mm, at most 8,559 units a query. Compile work is 378 units for `torus-slice` and about 139,000 for the two cuts with irrational edges. The step 3 refusal `torus slice r=3/5` (`rolling_quadrics.rs`), where no ball wider than 1/2 fits the section and the B-rep refuses with that radius limit, now reaches the field's own reason: the offset plane misses the offset torus (`field.rolling-path-reach`).

## The face registry for rolling pieces

Code: `Piece::face`, `Piece::on_face` and `Piece::ray` in `crates/catallax-field/src/rolling_region/piece.rs`, `face` on `Revolved`, `Generator` and `Coaxial`, the registration in `Plan::apply_run` (`fillet_rolling.rs`), the rolling leaf's `ray` and `on_face` in `geometry.rs`, and its face test in `Prepared::member` (`distance.rs`). Tests: `crates/catallax-field/tests/rolling_over.rs` (`a_sphere_corner_face_holds_its_rational_points`, `a_straight_piece_face_is_a_registered_face`).

A piece whose face lies on an exact surface registers that surface as a named face, under the piece's name: the ring torus about a centre circle wider than the ball (`Revolved` circles and `Coaxial` pieces), the ball itself where the centre stays on the axis (the sphere corner of step 5), and the cylinder of the ball's radius about a straight centre line (`Generator`). They are the surfaces OpenCASCADE's closed-form fillets build for the same pairs (`ChFiKPart_ComputeData_FilPlnPln` a cylinder, `_FilPlnCyl` a torus, or a sphere when `cylrad - Radius` is zero, read at V7_8_0 for steps 2 and 5). Each equation has the sign of the path's closed-form separation: the cylinder's and the sphere's are the separation times a positive constant, and the ring-torus polynomial factors as `((rho - C)^2 + h^2 - r^2) ((rho + C)^2 + h^2 - r^2)`, whose second factor is positive for `C > r`. So the equation is negative in every ball and positive in the region near the face.

Two leaf answers make the face usable by the plan, which needs a side for each ray from a boundary point and a face test there.

- Rays (`Piece::ray`). Where every wall and junction is strictly positive at a point of the face, the region near the point is where the separation is positive. The side of `p + t d` for small `t > 0` is then the sign of the first nonzero coefficient of the face's equation along the ray, as for a primitive's constraint. Where a wall or junction passes through the point, the ray stays unresolved, as before.
- Face membership (`Piece::on_face`). A boundary point of the region lies on a ball of the family, which is the face, or outside every ball on a wall. The separation is zero exactly on the face.

The plan now answers `Boundary`, with a unique outward normal, at the rational points of these faces, which it left `Indeterminate`: the sphere corners' points `c + (1, 2, 2)` of the boss, pocket and corner at 3/3, and the points `c + 3 (3/5, 4/5, 0)` of the upright's straight piece under the mitred corner at 3/3, where points 1/1000 off the face along the normal are cut and kept as in the B-rep. A horn or spindle torus has no global equation (`Surface::TorusSheet`), and an ellipse path's pipe, a saddle and a rim lie on no elementary surface, so those pieces still register no face. The shared points that `rolling_ends.rs` keeps per case lie on a support face where a wall of the piece passes; the side of such a wall came later ("Compare on rolling faces" below), and it settles P5's flange-plane point under the round at radius 1, in each pose, as the B-rep does. The contacts and the slot's rail ends, where a junction or a second wall also passes, stay `Indeterminate`.

With every live leaf named, the solver's whole-boundary distance (`Queries::distances`) no longer refuses `solver.distance-coverage` for such plans. The region's boundary is its face and parts of the two supports, which the body already registers, so the named faces cover the result's boundary. Plans with a horn, ellipse, saddle or rim piece keep the refusal.

A later fillet or chamfer on a body that holds a rolling piece still refused (`Plan::refuse_over_rolling`): to roll over a piece's face, the chain must take it as a support. Part 6 of step 5 lets the chain route do that for straight pieces; every other route and a chamfer keep the refusal.

## Compare on rolling faces

Code: `Piece::face_point`, `Piece::ray`, `Piece::wall_side`, `Piece::planes_at` and `Piece::smooth_side` in `crates/catallax-field/src/rolling_region/piece.rs`; `Prepared::rolling_witness` and `Prepared::survives` in `distance.rs`; `Plan::smooth_location` in `lib.rs`; `Leaf::smooth_side`, the rolling arm of `Leaf::planes_at` and the primitive sampler's fillet roles in `geometry.rs`; the candidate visit in `samples.rs`; and one shared change in `crates/catallax-compare/src/check.rs`. Tests: `crates/catallax-field/tests/rolling_tori.rs` (`torus_blend_vertices_have_field_face_witnesses`), and the shared points of `rolling_ends.rs`, which the field now settles where the B-rep does.

The done-when line asks catallax-compare to find every B-rep display vertex within the field's bound. Compare bounds a vertex from above with a point of a field face that the field proves on the result's boundary, or from the field's classification at the vertex with a witness. Four gaps stopped it on rolling faces.

- A witness on an irrational face. Compare's witness is a rational point of the face's surface near the vertex (`catallax_boolean::surface_witness`), and a ring torus of irrational radius has none. The witness is now the surface's nearest point in closed form (`Piece::face_point`), as OpenCASCADE's `Extrema_ExtPElS` finds it for a sphere, cylinder or torus: along the ray from the sphere's centre, the perpendicular from the cylinder's axis, or the ray from the torus's core point in the vertex's meridian. It lies on the surface by construction, so the separation, which is the surface's equation up to a positive factor, is zero there without an exact zero test of an algebraic number. Where every band, wall and junction of the piece encloses strictly above zero, the point is a boundary point of the piece on its face, and the plan's other leaves must leave it on the result's boundary (`survives`, the second half of the face-membership test). A rational surface keeps its rational witness.
- A wall on a support. A point where a piece meets a support face, off every ball, lies on one wall of the piece, and the plan left its rays `Indeterminate`. The path's walls are faces there (`Path::walls_are_faces`): off the balls the region's boundary is its two supports. With one wall zero and one support through the point, that wall lies on the support, and the region is the support's side that holds the balls, since each ball touches both supports from the region's side (`Piece::wall_side`, read at a station's centre). The support's equation then settles each ray, as a primitive's constraint does. On a plane support the piece also gives the plan its plane (`planes_at`), so the plan reads the cells of the planes through the point.
- One smooth surface. On a curved support, the plan has no planes to read and found no rays both ways, so it stayed `Indeterminate`. Where every leaf with the point on its boundary is, near the point, one side of the same smooth surface (equal monic polynomials with a nonzero gradient there), the result is constant on each side, and the two rays along the normal read it (`Plan::smooth_location`, after the planar and witness routes). Two tangent spheres have different polynomials and stay out.
- Compare and sampling. A vertex the field cannot classify, such as one where several pieces meet on a plane, now tries the same face witness before compare reports it, as compare's witness-first path already does for classified vertices (the one compare change). The field's candidate visit stopped at the first rolling face; it now skips rolling faces, which have no rational chart, visits the rest, and reports that it did not cover every face. A primitive registered under a fillet's blend or corner name (a corner's sphere) is sampled on its curved face instead of failing the whole visit with `field.sample-role`. Every candidate is still proven by a face query.

Compare over the rolling corpus rows, each inside 100M work, at the face registry's commit and at this change:

| Row | Face registry | This change |
| --- | --- | --- |
| `box-single` (control) | agree | agree, unchanged |
| `torus-slice` | 100M exhausted | 55.4M work, 19,880 forward vertices, no forward finding |
| `sphere-dome` | 100M exhausted | 21.9M, 11,519 vertices, no forward finding |
| `cone-top` | 100M exhausted | 17.1M, 5,724 vertices, no forward finding |
| `cone-base` | 100M exhausted | 40.1M, 9,436 vertices, no forward finding |
| `l-mixed-all` | 2.4M, 74 vertices `Indeterminate`, sampling refused | 4.6M, 14,305 vertices, no forward finding, 540 reverse samples |
| `vertex-4-one` | one reverse finding | the same |
| `roll-over-mitres` | 100M exhausted | agree: 1.1M work, 5,015 forward vertices, 228 reverse samples, no finding |
| `torus-intersection` | not run | agree: 46.2M work, 12,416 forward vertices, no finding |

`roll-over-mitres` needed one more thing, a name. RECIPE 4.5.1 names a junction sector `corner.vertex` "for their original input vertex": the upright's mitre corner is named after the vertex of the step's input body where the run ends, its two planes and the two rounds, as the B-rep names it. The field named it from the three planes, so compare found no field face under the B-rep's name and fell back to whole classifications.

What remains is reverse. Each row keeps one or two `selected B-rep face lacks complete sampled origin/direction coverage` findings: compare wants field samples on each blend face, and a rolling face has none. `l-mixed-all` keeps 24 reverse `Origins` findings on its concave strip near the mixed corners, where the B-rep names parts of the strip's cylinder `corner.vertex`; RECIPE 4.5.1 keeps `blend.edge` on every exposed part of an edge's run and gives `corner.vertex` to the junction sector, so these read as a B-rep naming difference for the B-rep's owners, not a field change.

## Step 7.1: certified distance to a patch

Code: `crates/catallax-field/src/patch_distance.rs` (`PatchDistance`). Tests: `crates/catallax-field/tests/patch_distance.rs`.

The field cannot hold a spline face yet (decisions 4.6, step 7). The first piece it needs is the distance from a point to a whole rational patch, with bounds it proves. `PatchDistance` holds each Bezier span of the patch as outward binary64 interval enclosures of its homogeneous controls `(X, Y, Z, W)`. The squared distance from `p` to the span's point `X / W` is `N / D`, with `N = |X - W p|^2` and `D = W^2`, both polynomials of twice the span's degrees. On a parameter cell their Bernstein coefficients come from the cell's controls by the product rule (Farouki and Rajan, CAGD 5, 1988). Every `D_k` is positive, so `N / D` is a convex combination of the ratios `N_k / D_k`: the smallest ratio bounds the squared distance from below on the cell, and the corner ratios are values of it. A best-first search halves the cell with the smallest lower bound, across its wider side, drops cells whose lower bound passes the smallest corner value, and stops when the square roots of the two bounds differ by at most `WIDTH`, 2^-24 mm (about 0.06 nm, below a thousandth of step 7.3's 100 nm tube). Both bounds approach the minimum as the square of the cell size (Garloff 1986). A query stops at 4,096 cells, the rolling region's panel cap, or at 52 halvings, and then returns its wider bounds marked as not converged. A query is sequential, so its answer and work are the same at any worker count.

How the reference kernels do this. OpenCASCADE's `Extrema_GenExtPS` and SMLib's `SmSurfPtGlobalSolver::FindRoughGuess` (`SmSurfacePointSolve.cpp`: the control polygon, or 7 by 7, 5 by 5 and 3 by 3 grids zooming on the best sample) take a sampled guess and polish it with Newton steps. Neither proves that no other point is nearer. The field keeps the subdivision bound instead, which is a proof, and needs no Newton step: the corner values reach the minimum as the cells shrink.

The spike the decisions ask for (7.1): 10^4 points on the M27 loft edge's two bilinear supports, and the same on the biquadratic edge's pair, half within 0.3 mm of the shared edge and half in a box round both supports. Each enclosure is checked against an independent binary64 search (a 64 by 64 grid, then twenty 9 by 9 grids, each a quarter as wide, round the best sample), which shares no code with the field: its value lies inside every enclosure, to 1e-12 mm^2. Release build on the work PC:

| Supports | Queries | Median | 99th percentile | Slowest | Cells, median and largest |
| --- | ---: | ---: | ---: | ---: | --- |
| M27 loft edge, bilinear | 20,000 | 55 us | 226 us | 1.8 ms | 51, 281 |
| Biquadratic edge | 20,000 | 219 us | 586 us | 1.5 ms | 51, 183 |

Every query converged within `WIDTH`, and the answers and work are the same at 1, 2 and 4 threads. The risk the decisions named first (a slow certified distance per query) is retired for these supports: well under a millisecond a query. A point on a support gets the lower bound 0 and an upper bound within `WIDTH`. Multi-span patches are covered (each span is a root cell); trimmed faces are not, as the decisions put them in 4.4.

## Step 7.2: whole patches in a field solid

Code: `crates/catallax-field/src/patch_shell.rs` (`Shell`, `Face`, `PatchFace`, `Polygon`), `PatchDistance::crossings`, `hull`, `boxes` and `fast_nearest` in `patch_distance.rs`, the `Shell` leaf in `geometry.rs`, its print definition in `printing.rs`, and the ruled-loft route `Plan::patch_loft` in `sweep.rs`. Tests: `crates/catallax-field/tests/patch_shell.rs`.

The decisions' model is the closed-mesh leaf: distance to the nearest face, sign by counting ray crossings. A `Shell` is a closed solid bounded by whole rational patches and planar polygons, placed by a rational rigid frame as a mesh is. Its producer proves the faces a closed, outward shell with no crossings. A point is classified in four steps.

1. Distance. Each patch's distance comes from step 7.1, each polygon's exactly in rationals. A point exactly on a polygon is on the boundary. A point whose distance bounds reach zero on a patch stays `Indeterminate`: no exact zero test of a patch point is made.
2. A point outside the faces' bounding box is outside.
3. Winding number. The ray from the point along the first of eight fixed directions crosses each face. A polygon's crossing is exact. A patch's crossings are the roots of `F_k = (X - W p) . e_k`, `k = 1, 2`, for two directions `e_k` across the ray, isolated by Bernstein sign tests and proved one by one with the Krawczyk operator (Krawczyk 1969; Moore, Kearfott and Cloud 2009, section 8.2) on halved parameter cells. At a root the Jacobian's determinant is `W^2 |e1|^2 (S_u x S_v) . d`, so its sign over the proving cell is the crossing's orientation, and the sign of `(X - W p) . d` says whether the root lies ahead of the point. A ray that meets a face border, touches a face or leaves a root unproved after 40 halvings tries the next direction; a point no direction settles stays `Indeterminate`.
4. The value is the true signed distance: the nearest face's bounds, negative inside.

How the reference kernels do this. OpenCASCADE's `BRepClass3d_SolidClassifier` and SMLib's `SmBrep::RayFire` also cast a ray and count its hits, finding them by sampling and Newton steps and judging tangent or border hits by tolerance before trying another ray. The field counts only proved transverse crossings and takes another direction otherwise, so no tolerance enters the count.

The print and display classifier treats a shell as it treats a mesh: a box that meets none of the faces' boxes holds one winding number, read at its lower corner. The faces' boxes are each polygon's vertex box and the control hulls of each patch span's 64 cells after three halvings per side. Without this, a shell, which has no carrier polynomials, would have been settled by one corner at any box inside its bounds.

The first producer is the field's own ruled loft (`Plan::patch_loft`), tried where the translated-loft route refuses: two parallel sections of line segments, one loop each, with as many segments. Each pair of corresponding segments gives a bilinear side, the straight lines between points at equal parameters, which is the B-rep's ruled side (the M11 loft's patches in `catallax-shapes/tests/spline_fillet.rs`). The slice at fraction `v` of the height is the polygon with vertices `(1 - v) a_i + v b_i`, so the turn at a vertex is a quadratic in `v`; its Bernstein coefficients `a' x a''`, `(a' x b'' + b' x a'') / 2` and `b' x b''` (edges `'` before and `''` after the vertex) all of one strict sign at every vertex prove every slice a strictly convex polygon. Then the sides meet only along the lines through corresponding vertices and the caps only along the sections' edges, which is the closed shell the leaf needs. Any other loft keeps `field.loft-correspondence-unsupported`; a correspondence turned half a turn, whose middle slice collapses, is the test's case.

Results (`patch_shell.rs`, release build on the work PC, 4.6 s for the file). On 1,521 grid points round the M11 loft, the shell answers 294 inside and 1,216 outside, each as the B-rep does, and 3 points exactly on a side `Indeterminate`, where the B-rep says boundary. The B-rep refuses 8 points on the vertical line through the loft's corner `(1, -1)` (`shapes.rational-patch-unproven`, "the line is tangent to the finite rational patch"); the field answers outside there. Every distance agrees with a binary64 search of the six faces, and its bounds are within 2^-24 mm. Cap points are on the boundary, side points at rational parameters stay `Indeterminate`, a rotated and mirrored copy answers as the original, and answers and work are the same at 1, 2 and 4 threads. Through `Plan::compile`, the field builds the M11 loft and agrees with the B-rep at 505 points. The print classifier on 1,521 boxes of side 1/4 mm round it settles 159 inside and 533 outside, each as the shell answers at its corners and centre, and leaves 829 near the faces unknown, among them a box across a side. Ordinary compilation builds the M11 loft and stops at the fillet with `fillet.lost-edge`. Sampled comparison and test fixtures use the independent spline proof from steps 7.3 and 7.4.

## Step 7.3: an independent spline spine

`catallax_field::spline_spine` proves a centre path between two whole,
single-span rational patches. The caller supplies the radius, the side of
each patch normal, a continuation parameter, a seed station and the part
tolerance. The field uses the two homogeneous Bézier nets and solves directly
for the centre, both contacts and the two positive normal lengths. It shares
rational patch and interval primitives with the B-rep.

OpenCASCADE marches a guide with `BlendFunc_ConstRad` and
`BRepBlend_Walking`, using a Newton correction on each section.
SMLib's `SmConstantRadiusFS` solves for equal points on the two offset
surfaces, and `SmFilletIntersector` checks interpolated spans. The field uses
SMLib's equal-offset condition and OpenCASCADE's Newton iteration limit
for the trace. [OpenCASCADE solver](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_9_0/src/BRepBlend/BRepBlend_Walking.cxx)
[SMLib solver](https://raw.githubusercontent.com/NVIDIA-Omniverse/usd-brep/main/source/SMLib/src/SmFilletStandardSolver.cpp)
[SMLib span checks](https://raw.githubusercontent.com/NVIDIA-Omniverse/usd-brep/main/source/SMLib/src/SmFilletIntersector.cpp)

Outward-rounded Krawczyk tests isolate the two border roots and certify a
unique root at every station in each span. Taylor models enclose the
residual about an interpolated centre and contact curve. At a shared
station, a root enclosure must fit both uniqueness boxes. The solver
isolates the station's root separately if neither span supplies an
enclosure tight enough. The border and seed roots must also fit their
spans. Contact bounds and inward derivatives at border roots keep both
contacts on the patches. The centre enclosure must fit
`min(100 nm, part tolerance / 2)` in Euclidean distance.
[Certified continuation](https://doi.org/10.1137/130906544)

All three M11 fixtures passed with eight spans each. The spine work
recorded at `1a7321fc` predates charging for basis-table arithmetic.

| Edge | Work | Tube radius (nm, rounded up) |
| --- | ---: | ---: |
| M27 loft edge | 428,558 | 0.009859 |
| Next loft edge | 431,184 | 0.008290 |
| Biquadratic edge | 656,302 | 0.006276 |

The certificates and work match at 1, 2 and 4 workers and the default
setting. A separate Newton solver eliminates the centre and normal lengths
and checks 1,001 stations per edge. Every centre and contact lies within
the enclosure to within the oracle's binary64 test slack. A tighter
part tolerance forces subdivision, and an unattainable tolerance returns
`field.spline-spine-eps` with the width reached.

This certificate covers the component selected by the seed. The next step
ties it to the named edge of a convex ruled loft.

## Step 7.4: the swept-ball distance

The distance query covers the short circular sections between the two
contact rails. A panel supplies a lower bound. Only a certified station
inside the finite domain, or a certified endpoint, supplies an upper
witness. A binary64 Newton proposal can place a witness sooner without
entering the proof. OpenCASCADE refines sampled point-to-curve extrema
in `Extrema_GExtPC`. SMLib uses short circular sections in
`SmCircularCrossSectionFSG`. The field uses the interval panel cover
from step 1 on the spline's certified centre and contact enclosures.
[OpenCASCADE](https://raw.githubusercontent.com/Open-Cascade-SAS/OCCT/V7_9_0/src/Extrema/Extrema_GExtPC.gxx)
[SMLib](https://docs.nvidia.com/smlib/manual/smlib/fillets/)

Differentiate the span's offset system twice to bound the centre's speed
and acceleration. If the centre moves at least `m` per unit of parameter
along a fixed direction and its acceleration is at most `M`, a chord
for a parameter increment `dt` has length at least `m |dt|`.
Its projection on a unit normal is at most `M dt² / 2`.
The normal sections at radius `r` therefore clear every other ball when
`r M < m²`. This is the normal/chord reach inequality.
[Federer, *Curvature Measures*, Theorem 4.18](https://doi.org/10.1090/S0002-9947-1959-0110078-1)

The two rails of a loft can meet one cap at different stations. Follow the
same offset component onto the exact extended charts, tying it to the
whole-patch proof at the seed. Every short arc direction crosses each
cap once, as its cap-height derivative keeps one strict sign. A distance
witness on a clipped section comes from its nearest circle point, a
retained rail point, or a circle-plane intersection. This keeps the part
between the caps even where its first rail has left the original chart.

The loft route supports one named convex edge between adjacent sides.
For a ruled loft, each side is an exact quadratic inequality obtained
from its moving horizontal edge. A quadratic remainder bound proves
the rolling ball stays inside both touched sides and clears every other
side. The section's radial normal proves that the new face projects
one-to-one between the two rails along a common outward direction. That
graph replaces the sharp corner while preserving the remaining sides
and caps. Whole patches from other sources and subsequent rounds on
the trimmed loft keep their typed refusals.

Sampled comparison opts into the tube proof. The tentative evaluator
locates each rail with 24 height bisections, then solves height and across
together with at most 12 two-variable Newton steps. It seeds the solve by
interpolating between the rail cuts. The finite differences use `2^-16`
of the normalized parameter interval. The solve stops once both updates
are at most `2^-24`, the rail bracket resolution.

The focused release tests sampled 505 independent arc points per edge.
Median distance query times on the work PC were 4.388 us for M27, 4.474 us
for the next edge and 7.215 us for the biquadratic edge. These were
shared-machine runs. The cap regression also tests a point above the
first rail's last station and below the cap, using an independent solve
for the second rail.

The capped sweep timing test uses nine independent binary64 arc samples
per group. The near-cut controls sit 0.2 mm inside the cap on sections
clear of both cuts. The on-cap samples come from sections whose two
rails straddle that cap. These median query times are in microseconds
from a focused release run on the shared Blade with field optimization
level 1.

| Loft edge | Cap height (mm) | Near cut | On cut |
| --- | ---: | ---: | ---: |
| M27 | 0 | 16.817 | 21.201 |
| M27 | 2 | 15.505 | 22.734 |
| Next M11 | 0 | 15.542 | 21.390 |
| Next M11 | 2 | 15.261 | 22.808 |

### Forward display comparison

The comparison uses the existing complete-cell positional grid at
0.005 mm. The two caps identify the trimmed fillet's boundary. Only
the selected fillet face supplies a grid or edge stations. The
biquadratic fixture displays the complete fitted carrier, which has
no body cuts. The forward comparison includes the positional display
bound and the B-rep fit bound. A vertex passes when the field proves
its distance is at most the combined bound and the widened upper bound
fits the display accuracy. A field lower bound above the combined bound
proves a gap. A range that still straddles the combined bound after
refinement reports `field.spline-bound-unresolved`.
The bare carrier test in `catallax-compare` charges the field proof and
sweep, B-rep fit, grid and vertex queries to one request. The comparison
totals and field build work recorded at `e30b732a` predate charging for
basis-table arithmetic and cap normals.

| Edge | Forward vertices | Total work | Largest comparison upper bound (mm) |
| --- | ---: | ---: | ---: |
| M27 loft | 245 | 36,448,251 | 0.002888316 |
| Next M11 loft | 335 | 63,995,055 | 0.004464352 |
| Biquadratic | 104 | 72,210,076 | 0.004916114 |

Both lofts passed every forward vertex check. With fresh inputs at 1, 2
and 4 workers, the field build records and work matched for all four rows.
The capped field builds used 901,321 work for M27, 897,045 for the next
loft edge, and 887,445 for the spline-plane loft. The bare biquadratic
field proof and sweep used 667,187 work. On that row, the B-rep fit
took 16,438,284, and the 4 by 26 carrier
grid took 54,998,702. The grid's display bound was 0.004836298487 mm.
Vertex queries took the remaining 105,903 work. All three comparisons
used the existing 100M request allowance.

The whole-body B-rep display remains open for both lofts. Reverse
normal and origin coverage for the fitted blends also remains open.

## Steps 7.5 to 7.7: near tangency and a plane chart

The synthetic edge uses one parabolic side and one plane. The plane's
slope at the shared edge runs from `1` to `2^-24` in dyadic steps.
The smallest tested angle with a full tube and reach proof was
0.003497057 degrees (`2^-14`), at 261,376 work and a tube radius of
`8.77e-11 mm`. The near-tangent and plane-control work recorded at
`69b8a363` predates charging for basis-table arithmetic. Lower angles
failed at the seed, a border proof or the support-domain check. The exact
tangent control refuses by name.
When the feet point into the same half plane, the chord between them
bounds the distance from either rail to every point of the short arc.
The sweep labels this bound when it fits the spine target. At the smallest
tested angle that passed, the rail bound was 6.142 nm.

An analytic plane stays a plane at its offset. Its rational chart is
translated exactly along the normal and enters the field solver at
zero offset. The original normal still gives its contact direction.
This follows the same-family offset route described in
[the viability study](research/cp2-grill/sdf-spline-viability.md#3-does-m10-accept-analytic-against-spline-offsets).
The spline-plane control checks the centre against the exact offset
equation `15x + 8y = 17/10`. The control's tube radius is
`1.224362e-14 mm` at 347,890 work. The field proof and work agree
at 1, 2 and 4 workers. The median query took 4.762 us over 99
independently generated points on the work PC.

**Step 7 and 4.2a remain open.** The B-rep spline-plane loft
fillet refuses with `fillet.unsupported-contact`. The comparison used
896,264 work at `e30b732a`, before the basis-table and cap-normal charges.
Forward display coverage for that row waits for the B-rep build.

## Open for the done-when line

- Samples on rolling faces, for the reverse coverage findings. Rational samples exist on rational tori, cylinders and spheres; irrational faces need another route.
- B-rep faces with no equation. `tangent-boss` and `narrow-cone` have horn-torus or fitted B-rep blend faces, which compare refuses before it reaches the field (`surface comparison needs a finite parametric correspondence solver`). The saddle and ellipse rows of steps 2 and 3 are fitted on the B-rep side too. For them the step tests run the same check directly against the pieces.
- The original and mirrored sphere-ring cases from `rolling_quadrics.rs` have membership checks. Both display comparisons remain open because the B-rep cannot display the ring bands. The frozen `sphere-cylinder` recipe has a separate B-rep build refusal.

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
| 3 | Quadric pairs (closed forms per Dupont, Lazard, Lazard, Petitjean 2008) | [#729](https://github.com/BohnBawerick/catallax/pull/729) |
| 4 | Ends at an edge or vertex, P5 at radius 2 | part 1, rim contact and hole rings: [#733](https://github.com/BohnBawerick/catallax/pull/733); part 2, the slot file: [#738](https://github.com/BohnBawerick/catallax/pull/738); part 3, `l-mixed-all` and `vertex-4-one`: [#742](https://github.com/BohnBawerick/catallax/pull/742) |
| 5 | Fillet over fillet, nested path | part 1, chains over an earlier straight round with the ball on its axis side: [#744](https://github.com/BohnBawerick/catallax/pull/744); part 2, sphere corners and cliffs; part 3, an upright under an earlier mitre; the face registry for rolling pieces and compare on rolling faces; part 4, a ball outside the round's cylinder (`mixed-ordered`); part 5, a ball narrower than a mitre's rounds (`narrow-mitres`); part 6, every edge at once and a third round over rolling pieces: this change |
| 6 | Torus pairs | `torus-slice`, `torus-intersection` and a cut below the centre plane: this change |
| 7 | Spline faces (7.0 to 7.7 of the decisions) | Whole patches, independent spine and sweep queries, and forward coverage for the three M11 spline-pair rows. Spline-plane display remains open with B-rep refusal `fillet.unsupported-contact`. |
