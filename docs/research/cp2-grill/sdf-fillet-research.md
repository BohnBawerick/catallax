# Rolling-ball fillets in the Catallax distance field: research and a milestone

Worktree: Catallax repository at bb2f8375 (CX below).
Labels: [read] = I read the file or page text this session; [snippet] = search-result summary only;
[unread/memory] = my own knowledge, not checked this session. Lines marked (inference) are my reasoning.

## 1. What the literature offers

### 1a. Blend operators that are not rolling balls
- Smooth minimum (Quilez). His own article: "all smooth-minimum will always produce non-exact SDF in some regions"; some families under-estimate distance in infinitely large regions. https://iquilezles.org/articles/smin/ [snippet]
- hg_sdf `fOpUnionRound` (Mercury). Source comment: the operators "produce correct distances or distance bounds ... on the condition that the gradients of the two SDFs are at right angles"; off by more than about 30 degrees the Lipschitz bound fails; worst case is close parallel surfaces. https://github.com/danielscherzer/SHADER/blob/master/libs/hg_sdf.glsl , https://in4k.github.io/html_articles/hg_sdf.html [snippet]. Body `max(r,min(a,b)) - length(max(vec2(r-a,r-b),0))` [unread/memory]. It is a quarter-circle in the (a,b) chart, so it is a true rolling ball only when a, b are exact SDFs with orthogonal gradients, which means perpendicular planes (or locally perpendicular faces near the edge).
- FRep / R-function blends (Pasko, Adzhiev, Sourin, Savchenko 1995; Pasko and Savchenko blend `g1+g2±sqrt(g1²+g2²-2αg1g2) + a0/(1+(g1/a1)²+(g2/a2)²)`). The displacement coefficients are "proportional to the algebraic distance", not a radius. Bounded blending (G. Pasko, A. Pasko, Kunii 2002, 2005) localises the blend with a bounding solid and supports a single selected edge or a partial edge blend. A 2024 critique (Synchronized Tracing, arXiv 2304.09673) notes that when operands are exact SDFs the result is "only a conservative SDF", and that standard bounded blends alter the field in an infinite volume. https://www.researchgate.net/publication/7942067 , https://arxiv.org/pdf/2304.09673 [snippet]
- Hoffmann and Hopcroft, potential method (Cornell TR 1985; SIAM "Geometric Modeling: Algorithms and New Trends" 1987 pp. 347-366). Parameters a, b, λ; a=b=r, λ=0 gives an elliptic blend (Curv docs). Algebraic, not constant radius. [snippet] https://github.com/curv3d/curv/blob/master/docs/shapes/Distance_Field_Operations.rst
- Rockwood, displacement method (ACM TOG 1989), Rockwood and Owen 1987 (same SIAM volume). Displacement blends; close to circular only for distance-function operands meeting at right angles [unread/memory].
- Middleditch and Sears, "Blend surfaces for set theoretic volume modelling systems", SIGGRAPH 1985: range-controlled blends [snippet, citation only]. Not constant radius [unread/memory].
- Barthe, Dodgson, Sabin, Wyvill, Gaildrat 2003 (2D potential fields), Gourmel et al. 2013 (gradient-based blend), BlobTree fillets (SMI 2010, UVic) which relate circular-section blends to canal surfaces [snippet]. Modelling operators for graphics, no certified radius.
- Bloomenthal: convolution surfaces and the 1997 implicit-surfaces book [unread/memory]. Not rolling ball.
Conclusion for this group (inference): none gives a constant-radius ball between general surfaces, and none gives better than a Lipschitz bound. They cannot reach B-rep parity. Only hg_sdf round matches, and only for perpendicular planes, a case the field already does exactly.

### 1b. Offsets: the one exact route in the literature
- Rossignac and Requicha. "Constant radius blending in solid modelling", Computers in Mechanical Engineering 3(1) 65-73, **1984** (PAP TM-42). "Offsetting operations in solid modelling", CAGD 3(2) 129-148, **1986** (PAP TM-53). Thesis "Blending and Offsetting Solid Models", PAP TM-54, 1985. The 1987 IBM J. R&D paper: "Constant-radius blends can be produced by combining growing and shrinking operations (collectively called offsetting) with Boolean operations." [snippet] https://www.semanticscholar.org/paper/ace482153ca4f053783340ac9ffbbe6884088079 . Note: the owner's notes say 1986; the blending paper is 1984, the offsetting paper 1986.
- Mathematically (well known, [unread/memory]): closing (grow r, shrink r) gives exact rolling-ball fillets on every concave edge; opening (shrink, grow) gives rounds on every convex edge. Both act on the whole solid, not on chosen edges.
- In a field the composition is not free (inference, standard): if d is an exact SDF, d - r is exact for the grown solid only outside the original; shrinking back by computing (d - r) + r returns d. A true closing needs the exact distance to the grown set's boundary, a new minimisation.
- Tools that do it on voxel grids: PicoGK (LEAP 71, OpenVDB) `DoubleOffset`, `TripleOffset`/`Smoothen`, `Fillet` ("creates a fillet-like effect") https://github.com/leap71/PicoGK/blob/main/PicoGK_Voxels.cs [snippet]; Altair Inspire "Fillet Implicit Geometry" with fillets/rounds/both https://help.altair.com/inspire/en_us/topics/implicit/fillet_implicit_t.htm [snippet]. Accuracy is the voxel size; all edges of the chosen type round.
- nTop: Offset Body, Boolean blocks with a "Blend radius", Smoothen Body; I found no rolling-ball edge fillet block. https://support.ntop.com/hc/en-us/articles/360041442674 , https://support.ntop.com/hc/en-us/articles/360056533394 [snippet]. Whether nTop's Blend radius is circular is unverified.
- Hyperganic: voxel kernel, assumed offset-based smoothing [unread/memory]. Shapr3D: Parasolid B-rep fillets, rolling-ball corners [snippet] https://support.shapr3d.com/hc/en-us/articles/7874402399900 . Not an implicit fillet.
- libfive: `blend_expt` (exponential), `blend_rough`, `blend_difference` [snippet] https://github.com/libfive/libfive/blob/master/libfive/stdlib/stdlib_impl.cpp . Fidget: an evaluator (interval arithmetic, tape simplification, JIT); I found no fillet operator. https://www.mattkeeter.com/projects/fidget/ [snippet]. Neither does true rolling-ball fillets that I could find.
- I found no implicit or field system that builds selected-edge, exact rolling-ball fillets between general surfaces with a certified distance. Absence is from about 12 searches, not proof.

### 1c. Ball-centre path (spine) machinery
- Natural quadrics (plane, sphere, cylinder, cone) are closed under offset. A torus offset is a torus. So between any two analytic CAD faces the spine is an intersection of two surfaces of the same families [unread/memory, standard].
- Quadric-quadric intersections: Dupont, Lazard, Lazard, Petitjean, "Near-optimal parameterization of the intersection of quadrics", SoCG 2003 and JSC 43(3) 2008 (three parts). Exact parametrisation, rational where the curve is not a smooth quartic, at most one extra square root otherwise. C++ implementation (QI; Lazard, Penaranda, Petitjean, CGTA 2006). https://inria.hal.science/inria-00186089 [snippet]. Torus-quadric analysis by projection and lifting: arXiv 2503.08924 [snippet, title only].
- Spline offsets: Elber and Kim 2001, multivariate rational spline constraint solver (IRIT). Zero-dimensional solutions are certified by subdivision with normal cones; one-dimensional sets (a spine is one) are "approximated by fitting ... to discrete points sampled", which is not a certificate. https://cs.nyu.edu/~exact/pap/solvingSys/myungsoo-elber-solve01.pdf [snippet]
- Distance to a rolling-ball face (inference, standard canal-surface fact): the face lies on the pipe of radius r about the spine C. For a point whose nearest spine point is interior and whose direction from it lies inside the contact sector, distance to the face is |dist(x, C) - r|, valid while r is below the spine's reach (no local self-overlap). So a certified point-to-spine distance gives a certified face distance. The field's rolling_region already uses this shape of test (interval panel cover over a 1-parameter centre curve, `d² > r²` or witness per panel) [read].
- Fillet over fillet (inference): an earlier round is a pipe of radius r1 about spine C1. Its offset by r2 is a pipe of radius r1 ± r2 about the same C1. The new spine is {c : dist(c, C1) = r1 ± r2, offset2(c)}. Taking C1's parameter t as an extra unknown, for each t the candidates lie on a circle in C1's normal plane, and meeting a plane or quadric offset gives at most 4 roots of a quartic in the circle angle. Nested, but still one-dimensional.
- Ends (inference): where the ball also touches a third face, edge or vertex, the spine switches piece: face offset ∩ edge tube (radius-r tube about an edge curve), tube ∩ tube, vertex sphere. P5's rim contact is exactly the cylinder offset ∩ mouth-circle tube, which CX rolling_region.rs implements by hand [read]. A localised closing/opening (1b) produces these pieces automatically.

## 2. Which approach can give certified or bounded distance

| Case | Spine | Certified route | Quality |
|---|---|---|---|
| plane-plane, plane-perp-cyl, 3-plane corner | line, circle | exists in CX | exact closed form |
| plane-cylinder any pose | ellipse | closed form; 1-param interval minimise | TrueDistance enclosure |
| plane-cone, plane-sphere | conic, circle | closed form | TrueDistance |
| plane-torus | spiric quartic, per torus angle linear in cos/sin | angle param with sqrt branch | TrueDistance |
| cyl-cyl, sphere-cyl, cone-cyl, coaxial | quartic | angle param with one sqrt (DLLP gives near-optimal radicals) | TrueDistance, more panels near branch points |
| torus-quadric, torus-torus | degree up to 8 | implicit curve, 2-param interval search | TrueDistance, slower |
| fillet over fillet | nested pipe | extra variable t, quartic per t | TrueDistance, 2-param search |
| spline-spline | not rational | precomputed certified spine tube (Hausdorff ε) from the field's own solve | bounded: enclosure width 2ε |
| end regions (P5) | piecewise, edge tubes and vertex spheres | per-piece as above, or localised closing | TrueDistance per piece; region after Boolean is Residual as today |

Localised closing/opening alternative: centres set E = {c : d_loc(c) ≥ r}, d_loc = exact distance to the trimmed support faces near the selected edges (CX distance_features.rs has trimmed face distances [read header]). Membership of x = branch-and-bound over c in B(x, r), 3 unknowns. Covers every pair, ends and corners with one routine. Cost is a 3D interval search per query (inference: 10³ to 10⁵ interval evaluations). Risk: closing rounds every concave feature in the window, and B-rep end treatments that are not ball envelopes (setbacks, face limits) will disagree by design.

## 3. Cost, difficulty, what CX has

Per query (inference, not measured): closed-form spine pieces need a 1-parameter interval cover, tens to a few hundred panels, so tens of microseconds in 96-bit outward arithmetic. Two-parameter searches (torus pairs, fillet over fillet) roughly square that. Spline spines cost a one-off certified fit per fillet, then a cheap polyline query.

Difficulty order (inference): plane-anything < quadric-quadric < ends < fillet over fillet < torus pairs < splines.

What CX already has [read]:
- crates/catallax-field/src/fillet.rs (883 lines): supports resolved from recipe origins, exact material-sector classification, finite runs, clearance proof `finite_round_clear`, refusal `field.fillet-contact` (line 24).
- crates/catallax-field/src/rolling_region.rs + rolling_region/sectors.rs (491 lines): a finite, interval-certified, 1-parameter centre curve region with whole-interval panel covers, 96 bits, ≤ 4096 panels, depth 64, typed refusals (docs/M18-ROLLING-REGION.md). This is the template for every spine piece.
- fillet_network.rs, fillet_smooth.rs: planar corner fans; CornerDefinition (FIELD.md:352).
- offset_spline.rs: reuses catallax-shapes `certify_elliptic_offset_sweep`; precedent for shared proof code.
- catallax-intersect M10 Bernstein cell search (docs/M27-SPLINE-FILLET-SCOPE.md): the direct six-variable spline offset search refused at depth 24 after 2.9M of 20M units.
- B-rep side catallax-shapes has rolling_pipe (plane-oblique cylinder ellipse), rolling_saddle (cyl-cyl), rolling_sphere, rolling_mitre (round over round), rolling_envelope, spline_fillet. The field must derive its spines itself to stay independent; sharing generic arithmetic is fine, sharing these modules would make the check common-mode.

Test coverage today (my grep of crates/catallax-features/tests): 23 fillet test files; 7 mention catallax_field (coaxial, fillet_over_fillet, generator, oblique, saddle, chamfer_after_fillet, trimmed_fillet_clearance), mostly asserting refusal.

## 4. Proposed milestone: a numbered sub-milestone of M18 "The field catches up" (roadmap 4.2); M36 is the next free top-level number per ROADMAP.md

Scope: constant-radius rolling-ball fillets compiled from recipe supports, any pair of analytic faces, earlier fillets, and ball contacts at fillet ends. Out of scope: variable radius, setback vertex blends, conic or G2 sections (non-ball by definition, keep typed refusals).

Steps, easiest first:
1. Generalise rolling_region into one certified point-to-spine-piece routine (parametrised piece, panel cover, sector test, reach check). Re-derive plane-plane and plane-perp-cylinder through it and require bit-identical answers to today's closed forms.
2. Plane with cylinder any pose, cone, sphere (conic spines). Target: oblique_fillets, oblique_cone_fillets, generator_fillets refusals cleared.
3. Natural quadric pairs (quartic spines): saddle_fillets, coaxial_fillets, fillet_sphere_ring, equal_crossing_fillets.
4. Ends: edge-tube and vertex-sphere pieces, generalising the P5 rim work: P5 at r = 2, rim_contact_fillets, hole_ring_fillets, slot_fillets.
5. Fillet over fillet (nested pipe): fillet_over_fillet, chamfer_after_fillet.
6. Torus pairs.
7. Spline-spline via certified spine tube: separate decision, after the M27 route change.

Done when: P5 builds in the field at radius 2 with no fillet refusal (the existing 4.2 line), and for every constant-radius rolling-ball case in the fillet test files (list frozen at milestone start), the field builds the part and catallax-compare finds every B-rep display vertex inside the field's certified distance bound at the part tolerance; each remaining refusal is named with its reason.

## Not verified
- nTop's blend radius profile; Hyperganic internals; Fidget fillet support.
- Exact contents of Rossignac-Requicha on selective (per-edge) blends.
- Per-query costs above are estimates, not measurements.
- Which of the 23 test files are pure rolling-ball cases (fillet_setbacks, variable_fillets, fillet_face_limit, numerical_fillet_admission likely not).
