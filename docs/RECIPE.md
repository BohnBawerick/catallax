# Catallax recipe language

Version 3, format `catallax-recipe/3`. [PRD.md](PRD.md) sections 2, 4 and 10 step 3 ask for this document.

A recipe is the feature tree, the master model of a part. The B-rep kernel and the field kernel both build the part from it. This document gives each operation its parameters, one mathematical definition both kernels can be checked against, the names of the faces it makes, and the inputs that are impossible. It defines how names survive upstream edits, how a lost name fails, the file format, and how SolidDraft's present references map onto Catallax names.

The crate `catallax-recipe` holds the data model: the types, their canonical JSON text and the checks that need no geometry. `catallax-field` evaluates stock solids and booleans, with the certified queries described in [FIELD.md](FIELD.md). `catallax-features::evaluate` builds named B-rep stock solids, sketches, extrusions, revolutions, through-all cuts and booleans. Its current rational analytic representation leaves the legal-input bugs in section 7.14 open. The two evaluators do not yet cover the whole language.

The exchange contract in [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md) is unchanged. Its question format stays the way the independent tester talks to the B-rep half. The recipe language sits beside it and does not replace it.

Where a definition is open, section 7 says so and lists the options. Every design choice in this document comes with its reason and one alternative with that alternative's cost.

## 1. The model

### 1.1 Steps and values

A recipe is a list of steps. Each step has an id and runs one operation. An operation reads the values of earlier steps by their ids and makes one new value. It never changes a value it reads.

A value is one of three kinds:

| Kind | Made by | Read by |
|---|---|---|
| profile | `sketch` | `extrude`, `extrude_through`, `revolve`, `sweep`, `loft` |
| path | `path` | `sweep` |
| body | every other operation | the operations that act on solids |

The recipe names one step as its `result`, which must make a body.

- **Operations make new values and never edit their inputs.** The reason is that one definition per operation then reads as a function from inputs to output, and any step can be rebuilt or compared on its own. SolidDraft's kernel slot works the same way. The alternative is features that edit "the part" in place, such as an extrude that merges into whatever came before. Recipes get shorter, but every such feature hides a boolean whose naming must also be defined.
- **A step reads only earlier steps.** The recipe is then a sequence a rebuild walks once, with no cycles to detect. The alternative is a graph whose order the reader works out. Steps could then be listed in any order, at the cost of a cycle check and a second order to keep stable.

### 1.2 Ids

Step ids, and element ids inside a step (sketch segments, path segments, hole placements), use the grammar of the exchange contract's step names: 1 to 64 characters from `a-z 0-9 _ -`.

- **Faces are named through ids that the user's edits leave alone.** Changing a parameter never changes an id. Inserting or deleting another step or segment never changes one either. The alternative is positional indices, such as "segment 3 of loop 0". Those need no ids in the file, but inserting one segment renames the side face of every segment after it.

### 1.3 Numbers, units and placements

- Lengths are millimetres. Angles are degrees.
- Every geometric number is a JSON string in the exact rational grammar of piece 1, the grammar of the exchange contract (ND2): a decimal with an optional exponent, such as `"1.1"` or `"-2.5e-3"`, or a fraction, such as `"1/3"`. A JSON number where a geometric value belongs is an error. Counts and indices, which are not geometry, are JSON integers.
- A placement is a pose as in the exchange contract (ND1): a quaternion `[w, x, y, z]` of exact numbers, not all zero and of any length, then a translation. A sketch plane is the image of the local plane `z = 0`, with sketch coordinates `(u, v)` along local `x` and `y`. Its normal is local `+z`. That normal is a unit vector with rational coordinates, because a rational quaternion gives a rational orthogonal matrix.
- A direction, such as an extrude direction or a hole axis, is a vector of any non-zero length.

Choices:

- **Exact rationals, as in the exchange contract.** Both kernels then read the same number, and no value passes through a double. The alternative is JSON numbers. They are shorter, but every reader rounds `0.1` to a different binary value before the kernel sees it.
- **Angles in rational degrees.** A rational number of degrees is a rational multiple of pi, so its sine and cosine are algebraic numbers and the B-rep can hold them exactly. Their degree can be high: `cos 1°` has degree 48. The alternative is giving each angle as a rational point on the unit circle, such as `(3/5, 4/5)`. Every value then stays rational, but no user can type "30 degrees".
- **Poses by quaternion, not by origin, normal and x-direction.** A quaternion gives an orthonormal frame by construction, and every rational rotation matrix comes from an integer quaternion. The alternative is three vectors, as SolidDraft's slot uses. Users read those easily, but rational vectors are rarely exactly unit length and orthogonal, so the frame needs a tolerance check.

### 1.3.1 Algebraic evaluation of rational inputs

The 2026-10-22 architecture decision changes evaluated geometry, not the rational number input grammar. A rational degree angle is converted to certified algebraic sine/cosine values. Normalising a direction or computing an arc radius may also introduce algebraic values. Both evaluators keep these expressions exactly in the shared shape definitions; they must not reject a legal input merely because the result is not rational. An arc centred at the origin through `(1,1)` has radius `sqrt(2)` without requiring a new number syntax.

An oblique circular extrusion reduces to an elliptic-cylinder support with the original cap planes and finite sweep domain. A supported meridional revolve reduces to analytic supports with exact angular trims. The recipe and segment origins remain the definition of the operation after reduction. Unsupported trimmed revolution carriers remain implementation gaps, not additional impossible inputs. Generic derived-surface storage for future profiles follows the architecture decision when its full proof/query path is supplied.

### 1.4 Solids, faces and agreement

A body is a solid: a bounded closed set `S` in space that is the closure of its interior. Set operations are regularised: `A ∪* B` is the closure of the interior of `A ∪ B`, and likewise for `−` and `∩`. The empty set is a valid body. Its result is the empty body, not a refusal.

Each definition below defines the point set `S` the step makes. The face structure is the maximal one on its boundary, as in the B-rep half's contract (ND6): no edge separates two faces on the same surface with the same sense. Each face carries the names of section 3.

Agreement between the kernels is agreement on the set and on the names. The B-rep's boundary must be the boundary of `S`. The field must be negative inside `S`, positive outside and zero on its boundary. It need not be the true signed distance, except where a definition reads distance values: offset, shell, thicken, smooth blend and graded features. At every boundary point away from edges, both kernels must report the same origin set.

### 1.5 Refusals

Each operation lists its impossible inputs. Per PRD section 3.2 these are the only refusals a finished kernel may give. Both kernels must refuse exactly these inputs, with the same stable code. The code is assigned when the operation lands in a kernel.

Every operation also refuses a result that touches itself along a curve or at a point where edges meet, because no manifold shell can bound it. This is today's `boolean.self-touching-shell`. Each operation's list leaves it out.

Recipe errors are a separate class. They are mistakes in the recipe, not in the geometry. Section 3.6 lists them.

## 2. Profiles and paths

### 2.1 Sketch

| Parameter | Value |
|---|---|
| `plane` | pose of the sketch plane |
| `loops` | list of loops; each loop is a list of segments |

Segments, in sketch coordinates:

| `kind` | Fields | Meaning |
|---|---|---|
| `line` | `id`, `to` | straight segment from the previous segment's end to `to` |
| `arc` | `id`, `centre`, `to`, `turn` | arc of the circle about `centre` from the previous segment's end to `to`, counterclockwise (`"ccw"`) or clockwise (`"cw"`) seen from the side the normal points to |
| `circle` | `id`, `centre`, `radius` | a full circle, alone in its loop |

Each segment starts where the previous segment in its loop ends. The first segment starts where the last ends. Every loop is closed by construction.

**Definition.** Each loop is a closed curve `L_i` in the sketch plane. The profile's region is

`R = closure { q in the plane, not on any loop : q is enclosed by an odd number of loops }`.

**Names.** A sketch makes no faces. Its segment ids name the faces that sweeps make from them, as `side` roles of the sweeping step.

**Impossible.**

- an arc whose two ends are not at exactly the same distance from its centre;
- a line or arc whose end equals its start;
- a circle with a radius of zero or less;
- a loop that crosses or touches itself;
- two loops that cross or touch.

A loop that touches another loop makes a region whose sweep touches itself along a curve. A sketch with no loops is allowed and has an empty region.

Choices:

- **The recipe stores solved geometry, not constraints.** The kernels need positions, and constraint solving is the application's job. The alternative is storing the constraint system and letting both kernels solve it. Tangencies would then stay exact, see section 7.1, but both kernels would need the same exact solver.
- **Each segment gives only its end point.** Closure needs no tolerance, and two segments can never disagree about a shared point. The alternative is both ends per segment, as SolidDraft does. That reads better on its own, but it needs a join tolerance, and a gap below it is a hidden approximation.
- **The region is "enclosed by an odd number of loops".** Islands inside holes need no extra structure, and loop order and direction do not matter. The alternative is an explicit outer loop with a list of holes. It is simpler to read, but an island inside a hole then needs a second profile.

### 2.2 Path

| Parameter | Value |
|---|---|
| `start` | point where the first segment starts |
| `segments` | list of path segments |

| `kind` | Fields | Meaning |
|---|---|---|
| `line` | `id`, `to` | straight segment to `to` |
| `arc` | `id`, `centre`, `to`, `axis` | arc about `centre`, counterclockwise about `axis` by the right-hand rule, to `to` |

**Definition.** A path is the curve `γ(s)`, `s` in `[0, L]`, the segments joined end to end and parametrised by arc length.

**Names.** A path makes no faces. Its segment ids name sweep faces.

**Impossible.**

- a line whose end equals its start;
- an arc whose ends are not at exactly the same distance from its centre, or whose ends do not both lie in the plane through the centre perpendicular to `axis`;
- an arc whose end equals its start;
- a junction where the two tangents are not exactly parallel and pointing the same way.

All of these checks are exact on rational input. A path with a corner is not defined yet, see section 7.6.

## 3. Names

### 3.1 Origins

An **origin** says which step made a face and which part of the step's operation made it:

```json
{ "step": "base", "role": { "side": { "segment": "top" } } }
```

Each operation's section lists the roles it makes. All of them:

| Role | Made by | Fields |
|---|---|---|
| `start`, `end` | block, cylinder, cone, extrude, revolve, sweep, loft | none |
| `side` | stock solids, extrude, extrude_through, revolve, sweep, loft | `segment`; a fixed segment id for stock solids; for a sweep also `path` |
| `hole` | hole | `at` (placement id), `part` |
| `blend` | fillet | `edge` (edge name) |
| `bevel` | chamfer | `edge` |
| `corner` | fillet, chamfer | `vertex` (vertex name) |
| `inner` | shell | `of` (origin) |
| `round` | shell, offset | `edge` |
| `round_corner` | shell, offset | `vertex` |
| `offset` | thicken | `of` |
| `rim` | thicken | `edge` |
| `copy` | pattern, mirror | `instance`, `of` |
| `beams` | lattice | none |
| `surface` | gyroid | none |
| `smooth` | smooth blend | none |
| `imported` | STEP import | `key` |
| `mesh` | mesh import | none |

Sketches, paths, booleans and drafts make no origins.

Origins have a canonical order: the order of their compact canonical JSON text, compared byte by byte. Any implementation can compute it from the text alone.

### 3.2 Faces carry origin sets

Every face of every body carries a non-empty set of origins. The sets follow these rules:

1. **Made.** A face an operation makes from its own geometry gets that operation's origin.
2. **Kept.** A face of an input that survives into the result, perhaps trimmed, keeps its set.
3. **Split.** When an operation divides a face, every piece keeps the whole set.
4. **Merged.** When faces merge into one, because they lie on the same surface and share an edge, the merged face takes the union of their sets. A union of two blocks with coplanar tops has one top face with both tops' origins.
5. **Replaced.** When an operation replaces a face's surface one for one, as draft and offset do, the new face keeps the old face's set.
6. **Copied.** A copy made by a pattern or mirror wraps each origin: `{"step": "<pattern>", "role": {"copy": {"instance": k, "of": <origin>}}}`.
7. **Gone.** A face with no part left on the boundary is gone. Its origins may live on in faces that merged with it.

The set of a face depends only on the recipe and on the point set each step makes. It never depends on the order in which a kernel lists faces, or on its handles.

- **Origin sets rather than one origin per face.** The B-rep half merges faces to keep topology maximal (ND6), and a merged face really has two makers. The alternative is keeping one origin per face and never merging faces of different origins. Names would stay one-to-one, but the B-rep would have to keep edges that separate nothing, against its contract, and STEP files would carry them.
- **A boolean makes no faces.** Every face of a regularised set operation lies on a face of an input, so rules 2 to 4 name it. The alternative is giving a boolean its own `split` origins for pieces. Pieces would get unique names, but every boolean would rename the faces a later fillet refers to.

### 3.3 Edge and vertex names

An edge lies between two faces. Its **edge name** is the pair of the first origin of each face, in canonical order:

```json
{ "between": [ { "step": "base", "role": "end" },
               { "step": "base", "role": { "side": { "segment": "top" } } } ] }
```

A **vertex name** lists the first origin of each face that meets at the vertex, in canonical order. Ordinary vertices have at least two origins. A singular vertex such as a cone apex has one, the wall's origin. The list is never empty and never repeats an origin.

Edge and vertex names appear inside roles: a fillet names its blend face by the edge it rounds. The edge or vertex they name is on the step's input body.

- **Edges are named by their faces, never by their own maker.** On a manifold solid every edge lies between two faces, so face names are enough. The PRD's own example is "the edge between the side face of extrude 1 and the floor of cut 2". The alternative is an origin per edge as well. Some edges, such as intersection curves, would get more precise names, but every operation would need edge rules as well as face rules, and the two would have to agree.

### 3.4 References

A step that acts on chosen faces or edges names them by **references**, resolved against the body it reads.

A reference term is an **origin pattern**: an origin whose role may be left out. `{"step": "boss"}` matches every origin step `boss` made. A pattern matches a face when it matches any origin in the face's set.

| Reference | Fields | Chooses |
|---|---|---|
| face reference | `face`, optional `touching`, optional `near` | faces matching `face` |
| edge reference | `between` (two patterns), optional `touching`, optional `near` | edges between a face matching one pattern and a face matching the other |
| edge selection | `"all"`, or a list of edge references | every edge, or the union of what the references choose |

Resolution has three steps:

1. Take every face or edge that matches.
2. If `touching` is present, keep a face only when it shares an edge with a face matching each listed pattern. Keep an edge only when one of its end vertices lies on a face matching each listed pattern. A ring edge has no vertices, so `touching` never keeps it.
3. If `near` is absent, the reference chooses all that remain. If `near` is present, it chooses the one at strictly the smallest distance from that point.

A reference without `near` deliberately chooses every piece. "Fillet the edge between the top and the side" still means every piece of that edge after an upstream cut divides it.

- **Topology first, a point last.** `touching` survives dimension changes, which move every point. `near` breaks ties that topology cannot, such as two circles where a sphere meets a coaxial cylinder. The alternative is a point only, as SolidDraft's positional selectors use. It needs no names at all, but any edit that moves the part changes what it chooses.
- **An exact tie under `near` fails.** A rebuild must give the same answer on every machine, and both kernels decide distances exactly. The alternative is a tie-breaking rule, such as the lower name. Rebuilds would always succeed, but a near-symmetric part would choose a side the user never saw.

### 3.5 How names survive upstream edits

| Edit upstream | What happens to names |
|---|---|
| A parameter changes | Nothing. No name holds a dimension. |
| A step is inserted, or a step no name mentions is deleted | Nothing. |
| A step is renamed | The editor rewrites every name that mentions it. Names are structured, so the rewrite is exact. |
| A step is deleted that names mention | The structural check fails each such name (`recipe.unknown-origin`). |
| A sketch segment changes shape, for example a line becomes an arc | Nothing. The role holds the id, not the surface type. |
| A sketch segment is deleted | Names of its faces fail the structural check (`recipe.unknown-element`). |
| A sketch segment is split in two | The sketch editor keeps the old id on one piece and gives the other a new id. The rule for which piece keeps it is the editor's, see section 7.13. |
| New geometry divides a face | References without `near` choose every piece. References with `near` choose the nearest. |
| Two faces merge | A face reference to either origin finds the merged face. An edge between them is gone, so references to it fail. |
| A face disappears | References to it fail (`recipe.missing-reference`). |
| A pattern's count shrinks | Names of copies it no longer makes fail the structural check (`recipe.no-such-copy`). |

### 3.6 How a lost name fails

The structural check, in `catallax-recipe`, catches every name that can never exist. A name may mention only an earlier step, with a role that step's operation makes, with element ids that step has, and with a copy instance it makes. These failures are recipe errors with stable codes such as `recipe.unknown-origin`, `recipe.wrong-role` and `recipe.unknown-element`. The crate lists all of them.

A rebuild catches what only geometry can show:

| Code | When |
|---|---|
| `recipe.missing-reference` | a reference chooses nothing |
| `recipe.ambiguous-reference` | `near` ties exactly, or a two-distance chamfer edge matches its two patterns both ways round |
| `recipe.input-failed` | a step reads a step that failed |

A failed step fails every step that reads it, directly or through other steps. A recipe either builds its result completely or reports the first failing step and its code. A kernel never returns part of a result.

- **No skipping past a failure.** The PRD forbids showing the user a half-built state. The alternative is SolidWorks' behaviour: suppress the failed feature and build the rest. The user would see more of the part, but it would be a body the recipe does not describe.
- **The application keeps the last good build.** Showing that while the user fixes the recipe is the application's choice, not the language's.

### 3.7 What the kernels must report

Each kernel reports, for every face of every body it builds, its origin set. In the B-rep this is a label on each face. In the field kernel it is the origin of the input primitive whose field is active at a boundary point, which min and max operations already identify. The agreement checks of PRD section 3.5 compare these sets at sampled boundary points, as well as the surfaces.

## 4. Operations

Notation: `S` is the input body, `∂S` its boundary, `d(p, X)` the distance from a point to a set, `B(c, r)` the closed ball. For a face `F`, `S_F` is its surface extended beyond the face, and `n_F` its unit normal pointing out of the material.

### 4.0 Stock solids

Sketch, extrude and revolve remain the main way users build parts. Stock solids give both kernels the exact inputs the B-rep builders already accept, without requiring those feature evaluators first. Every stock operation makes a body and takes a required `pose`, the rigid placement of section 1.3. The sets below are in local coordinates before that pose. No dimensions pass through floating point.

| `op` | Parameters in canonical order after `op` | One local solid definition |
|---|---|---|
| `block` | `size: [sx, sy, sz]`, `pose` | `0 <= x <= sx`, `0 <= y <= sy`, `0 <= z <= sz` |
| `cylinder` | `radius: r`, `height: h`, `pose` | `x² + y² <= r²`, `0 <= z <= h` |
| `sphere` | `radius: r`, `pose` | `x² + y² + z² <= r²` |
| `cone` | `base_radius: rb`, `top_radius: rt`, `height: h`, `pose` | `0 <= z <= h`, `x² + y² <= (rb + (rt-rb)z/h)²` |
| `torus` | `major: R`, `minor: r`, `pose` | `(sqrt(x²+y²)-R)² + z² <= r²` |

The cone definition includes capped expanding and contracting frusta, the top apex at `rt = 0`, and the cylinder at `rb = rt`. The current B-rep builder still refuses equal radii, even though that geometry is possible. The field evaluator accepts the cylinder definition; this does not change the B-rep builder. Ring tori require `0 < r < R`, so the equivalent quartic has no extra solid branch.

Impossible inputs:

- block: any side length not positive;
- cylinder: radius or height not positive;
- sphere: radius not positive;
- cone: base radius or height not positive, or top radius negative;
- torus: tube radius not positive, or major radius not greater than tube radius;
- every stock solid: the zero rotation quaternion.

The field evaluator reports `field.impossible-primitive` for invalid dimensions and `shapes.zero-quaternion` for the invalid pose. The recipe structural reader does not perform geometric validation. B-rep refusal codes remain unchanged; unifying them belongs to the B-rep recipe adapter, not to this field implementation.

**Face names.** The creating stock step supplies `step` in every origin.

| Solid | Roles and exact trimmed domains |
|---|---|
| block | `start` at `z = 0`, `end` at `z = sz`; `side` with fixed segment ids `x_min`, `x_max`, `y_min`, `y_max` on the corresponding coordinate bounds, trimmed by the other bounds |
| cylinder | `start` and `end` for the two closed disks; `side` with segment `wall` for `x²+y² = r²`, `0 <= z <= h` |
| sphere | `side` with segment `wall`, the whole sphere |
| cone | `start` for the base disk; `end` for the top disk only when `rt > 0`; `side` with segment `wall` for the conical side between the caps, or the cylindrical side for equal radii |
| torus | `side` with segment `wall`, the whole torus |

For example, a sphere's face is `{"step":"ball","role":{"side":{"segment":"wall"}}}`. `path` is absent for every stock `side`. Fixed segment ids are local to their stock operation, not user-created sketch elements. A cone apex has no zero-area `end` face.

**Edge and vertex names.** Section 3.3 applies without arena indices. Each box edge is between two adjacent named planes; each corner lists its three named planes. A capped cylinder has two ring edges, between `start` and `wall`, and between `end` and `wall`, and no vertices. A frustum has the same rings where their radii are positive. A cone with a top apex has one base ring and one singular vertex named by the singleton wall origin. A sphere and a torus have no physical edges or vertices. Parameter seams and sphere poles create none.

**Equivalent features.** A stock solid must give the same set and named result as its canonical equivalent extrude or revolve once those evaluators exist in the B-rep half. The creating stock step id stays the identity of that equivalent operation; internal construction steps must not leak different origins.

- A block is a rectangle extruded by `sz`. Its profile segments use the fixed side ids above, so the generated `start`, `end` and `side` roles agree directly.
- A cylinder is a circle with segment id `wall`, extruded by `h`.
- A sphere is a half-circle meridian `wall` revolved through 360 degrees.
- A torus is a circle meridian `wall`, centred at radial distance `R`, revolved through 360 degrees.
- A cone is its axial trapezoid or triangle revolved through 360 degrees. The sloped meridian generates `side.wall`. The lower and upper radial meridians generate the named `start` and `end` disks. An expansion to generic revolve must map those radial-segment side origins to these cap origins before exposing the stock result. The axis generates no face.

This remains a naming obligation on the eventual canonical feature expansion. The current B-rep feature evaluator preserves explicitly authored sweep origins but does not yet certify stock-to-feature equivalence or perform the cone-cap origin remapping. The field evaluator does not implement sketch, extrude or revolve.


### 4.1 Extrude

| Parameter | Value |
|---|---|
| `profile` | sketch step, region `R`, unit normal `n` |
| `distance` | `d`, not zero |
| `direction` | optional vector `w` |

**Definition.** The sweep vector is `v = d n`, or `v = d w / |w|` when `direction` is present. The solid is `E = { q + t v : q in R, t in [0, 1] }`.

**Names.**

- `start`: the face `R`, at `t = 0`;
- `end`: the face `R + v`;
- `side` with `segment` `s`: the face `{ q + t v : q on segment s, t in [0, 1] }`.

Two consecutive collinear lines make coplanar side faces. These merge (ND6), and the merged face carries both `side` origins.

**Impossible.** `d = 0`. `w = 0`. `w` parallel to the sketch plane, `w · n = 0`, which gives zero thickness.

- **A negative distance is allowed.** It extrudes the other way, and `start` stays in the sketch plane. Flipping the direction is then a parameter edit that keeps every name. The alternative is requiring `d > 0`, as SolidDraft's slot does, and flipping by `direction`. One sign rule is less to explain, but a flip becomes a different kind of edit.

#### 4.1.1 Through-all cut

Version 3 adds `op: "extrude_through"` with fields in this order:

| Parameter | Value |
|---|---|
| `target` | body `S` |
| `profile` | sketch step, region `R`, unit normal `n` |
| `direction` | optional non-zero vector `w`; defaults to `n` |

**Definition.** Let `T = { q + t w : q in R, t in the real numbers }`. The result is `S −* T`. The tool extends in both directions from the sketch plane. A one-sided blind cut is an ordinary signed `extrude` followed by `boolean` subtraction.

**Names.** Target faces keep their origins. Tool walls carry this step's `side` origins, with the sketch segment ids. There are no `start` or `end` origins because the infinite tool has no caps. Collinear or cocircular consecutive segments merge their origin sets.

**Impossible.** `w = 0`, or `w · n = 0`. Empty profiles leave the target unchanged; an empty target stays empty. The self-touching-shell rule still applies.

The B-rep implementation uses finite construction caps whose exclusion is proved from the target's defining bounds, not an arbitrary large length. Let the target lie within the sphere of radius `B` about model-space zero, and let `o` be the sketch origin. Choose

`k = ((B + |o|₁) |n|₁) / |w · n| + 1`.

For every target point, its line parameter `(p − o) · n / (w · n)` lies strictly between `−k` and `k`. Sweeping the shifted profile `R − k w` by `2 k w` therefore gives exactly the same cut inside `S`, with both caps outside it. Stock and swept definitions supply conservative bounds; union takes their maximum, intersection their minimum, and subtraction keeps the target bound. No finite construction cap can survive on the result boundary.

### 4.2 Revolve

| Parameter | Value |
|---|---|
| `profile` | sketch step, region `R` |
| `axis` | `point` `a` and `direction` `u`, a line in the sketch plane |
| `angle` | `θ` in degrees |

**Definition.** `Rot(φ)` turns by `φ` degrees about the axis, counterclockwise about `u` by the right-hand rule. The solid is the union of `Rot(φ)(R)` for `φ` in `[0, θ]`.

**Names.** `start`, the face `R`. `end`, the face `Rot(θ)(R)`. Neither exists when `θ = 360`. `side` with `segment` `s`, the surface `s` sweeps. A segment that lies on the axis sweeps no area and makes no face.

At 180 degrees, caps joined along an axis segment merge into one physical face carrying both `start` and `end` origins. The shared axis segment is not a physical edge. Disjoint caps remain separate faces.

**Impossible.** `θ <= 0` or `θ > 360`. `u = 0`. An axis that does not lie in the sketch plane. An axis that meets the interior of `R`, because the solid would then pass through itself. The region may touch the axis along segments or at points.

### 4.3 Hole and countersink

| Parameter | Value |
|---|---|
| `target` | body `S` |
| `at` | placements, each with `id`, `point` `p` and `direction` `w`, the drilling direction into the material |
| `diameter` | `D` |
| `depth` | `{"blind": h}` or `"through"` |
| `tip` | `"flat"`, or `{"cone": α}` with `α` the full point angle in degrees |
| `head` | optional `{"counterbore": {"diameter": Dc, "depth": hc}}` or `{"countersink": {"diameter": Dk, "angle": β}}` |

**Definition.** For one placement let `ŵ = w / |w|`. A point `x` has depth `z = (x − p) · ŵ` and distance `ρ` from the axis. The tool is `T = { x : ρ <= r(z) }`, where `r(z)` is:

- for `z < 0`, the entry: `r = r0`, the radius at depth 0. That is `Dc/2` with a counterbore, `Dk/2` with a countersink, and `D/2` otherwise;
- with a counterbore: `Dc/2` for `0 <= z <= hc`, then `D/2` down to `h`;
- with a countersink: `max(D/2, Dk/2 − z tan(β/2))` for `0 <= z <= h`;
- with no head: `D/2` for `0 <= z <= h`;
- below a blind hole with a cone tip: `D/2 − (z − h) tan(α/2)` while that is positive. With a flat tip, nothing below `h`;
- for `"through"`, `h` is infinite and there is no tip.

The result is `S −* (T_1 ∪ T_2 ∪ ...)`, one tool per placement.

**Names.** `{"hole": {"at": <placement id>, "part": <part>}}` with part:

| `part` | Surface |
|---|---|
| `entry` | the cylinder of radius `r0` above depth 0 |
| `counterbore_wall`, `counterbore_floor` | the counterbore's cylinder and its flat floor at `hc` |
| `countersink` | the countersink cone |
| `wall` | the cylinder of diameter `D` |
| `floor` | the flat floor at `h`, flat tip only |
| `tip` | the point cone, cone tip only |

A part exists only where it bounds the result. `entry` appears only where material lies above the placement point within `r0` of the axis. Faces of `S` keep their origins. The structural check rejects parts a hole cannot make, such as `floor` on a through hole.

**Impossible.** `D <= 0`, `h <= 0`, `w = 0`. `α` or `β` outside `(0, 180)`. `Dc <= D`, `Dk <= D` or `hc <= 0`, where the head would have no width or no depth.

Choices:

- **The tool continues above the placement point.** A hole started on a slanted or curved face then leaves no lip of material above its start plane. The alternative is a tool that starts exactly at depth 0. That is more local, but every hole on a non-flat face needs its start point placed by hand above the surface.
- **Placements are points and directions, not positions on a named face.** A hole then depends on no face's frame, and moving the face does not move the hole by surprise. The alternative is a face reference plus 2D coordinates. The hole would follow its face, but the language would need a rule for each face's 2D frame.
- **One step holds several placements, each with an id.** A bolt circle is one feature, and deleting one hole keeps the others' names. The alternative is one step per hole. That is simpler, but a bolt circle of twelve holes becomes twelve steps.

The B-rep boolean path now builds coaxial cylindrical countersinks and apex-cut cone cases. This does not yet implement the named hole operation or its full placement and origin-set rules; general cone/torus poses remain unfinished.

### 4.4 Boolean

| Parameter | Value |
|---|---|
| `kind` | `"union"`, `"subtract"` or `"intersect"` |
| `target`, `tool` | bodies `A`, `B` |

**Definition.** `A ∪* B`, `A −* B` or `A ∩* B`, regularised (section 1.4). `A` and `B` may be the same body.

**Names.** None made. Rules 2 to 4 of section 3.2 name every face.

**Impossible.** Only a result that touches itself (section 1.5).

- **One tool per step.** Each step is then exactly one set operation. The alternative is a list of tools. That is shorter for many unions, but the step's naming would need a rule for the order of its tools.

### 4.5 Fillet

| Parameter | Value |
|---|---|
| `body` | `S` |
| `edges` | edge selection |
| `radius` | `r` |

**Definition.** Take one chosen edge `e` between faces `A` and `B`. Call it convex where the material angle across `e` is less than 180 degrees, and concave where it is more. Let `σ` be the **spine**: the centres `c` at distance `r` from both `S_A` and `S_B`, on the material side for a convex edge and on the outside for a concave edge, in the component that runs along `e`. Each centre touches `S_A` at a foot point `a(c)` and `S_B` at `b(c)`. The **blend surface** is the union, over `c` in `σ`, of the shorter arc from `a(c)` to `b(c)` of the circle about `c` of radius `r` in the plane through `c`, `a(c)` and `b(c)`. The **blend region** `W_e` is the union, over the same `c`, of the part of that plane enclosed by the arc and the curves of `S_A` and `S_B` back to `e`.

The rolling ball `B(c, r)` touches `A` and `B` and rolls along `σ`. Where the edge ends at a vertex with a third face `C`, and the blend surface crosses `S_C` transversally, the blend continues along the spine's extension and is cut by `S_C`.

The result is `(S −* ∪ W_e over convex edges) ∪* (∪ W_e over concave edges)`. Every chosen edge uses the same `r`, and all edges are blended at once, so their order does not matter.

Exactness: two planes give a cylinder, a plane and a perpendicular cylinder give a torus. A blend between two curved faces has no closed STEP form. Both kernels hold it exactly by this definition, and only the STEP file approximates it (PRD section 3.1).

**Names.**

- `{"blend": {"edge": <edge name of e>}}` for the blend face;
- `{"corner": {"vertex": <vertex name>}}` for a patch where blends meet at a vertex, see section 7.2;
- `A`, `B` and `C` keep their origins.

**Impossible.** `r <= 0`. A chosen edge whose foot points leave `A` or `B` before the blend reaches the edge's ends, because no ball of radius `r` touches both faces there. SolidDraft's slot refuses a fillet that is too big in the same way. Section 7.2 lists the open alternative, rolling onto the next face.

- **The rolling-ball definition, as the PRD states it.** It is a pure function of the two surfaces and `r`, and the field kernel can evaluate it exactly because the recipe names both faces. The alternative is the morphological opening `(S ⊖ B_r) ⊕ B_r`. It is fully defined at every corner, but it rounds every convex edge of the body and closes narrow slots, not just the chosen edges.

### 4.6 Chamfer

| Parameter | Value |
|---|---|
| `body` | `S` |
| `edges` | edge selection |
| `distance` | `d1`, on the face matching the first pattern of each edge reference |
| `other_distance` | optional `d2`, on the other face; `d1` when absent |

**Definition.** Take one chosen edge `e` between `A` and `B`. The **contact curve** on `A` is the curve of `S_A` at distance `d1` from the curve of `e`, its intersection with the tube of radius `d1` about `e`, in the component next to `e`. The contact curve on `B` is the same with `d2`. The chamfer face is the ruled surface joining the point of each contact curve that lies in the normal plane of `e` at each point of `e`. `W_e` is the region between the ruled face and `e`. The result combines the `W_e` exactly as for a fillet, and ends meet a third face in the same way.

Exactness: two planes give a plane. A plane and a cylinder meeting at a circle give a cone.

**Names.** `{"bevel": {"edge": <edge name>}}` and `{"corner": {"vertex": <vertex name>}}`. A two-distance chamfer must list its edges, because `"all"` gives no first face. The structural check rejects that combination, and a rebuild fails an edge that matches both patterns both ways round.

**Impossible.** `d1 <= 0`, `d2 <= 0`. A contact curve that leaves its face before the chamfer reaches the edge's ends.

- **Setback measured in space from the edge.** It is the distance a user measures on the part, and on planes it equals the usual setback. The alternative is distance measured along each surface. That is the same on planes, but on a curved face it needs geodesics, which neither kernel computes exactly.

### 4.7 Shell

| Parameter | Value |
|---|---|
| `body` | `S` |
| `open` | face references, possibly none |
| `thickness` | `t` |

**Definition.** Let `K` be the part of `∂S` not in the open faces: the closure of `∂S` minus the open faces' interiors. The result is `{ p in S : d(p, K) <= t }`: the material within `t` of a kept face.

Points near an open face but far from every kept face are removed, so the cavity reaches through the open faces. The cavity is the set of points of `S` farther than `t` from `K`. At convex edges the walls' inner faces meet sharply. At concave edges and vertices of `K`, the inner surface is a cylinder or sphere of radius `t`, so the wall is exactly `t` thick everywhere.

**Names.**

- `{"inner": {"of": <origin>}}` for the inner wall at distance `t` from a kept face;
- `{"round": {"edge": <edge name>}}` for the inner surface at distance `t` from a concave edge;
- `{"round_corner": {"vertex": <vertex name>}}` at distance `t` from a concave vertex;
- kept faces keep their origins. What remains of an open face is its rim, the ring of the face within `t` of `K`, and it keeps the open face's origin.

**Impossible.** `t <= 0`.

- **Distance from the kept faces.** This is the PRD's description of the field side: "the field measures distance from the faces themselves". It gives exactly constant thickness. The alternative is offsetting each face and extending neighbours until they meet. That matches most CAD tools, but the walls are thicker at concave corners, and section 7.3 explains why the choice matters after booleans.

Only inward shells are defined. Section 7.3 lists outward shells as open.

### 4.8 Offset

| Parameter | Value |
|---|---|
| `body` | `S` |
| `distance` | `d` |

**Definition.** For `d > 0`, `{ p : d(p, S) <= d }`. For `d < 0`, `{ p in S : d(p, outside of S) >= |d| }`. For `d = 0`, `S`.

**Names.** Each face's offset replaces it and keeps its origins (rule 5). An outward offset adds `{"round": {"edge": ...}}`, a pipe of radius `d` about each convex edge, and `{"round_corner": {"vertex": ...}}`, a sphere about each convex vertex. An inward offset adds them about concave edges and vertices. Faces that vanish are gone.

**Impossible.** Nothing beyond section 1.5. An inward offset that removes everything gives the empty body.

- **True distance, so convex edges round over.** This is the definition the field kernel evaluates natively, and offsets of planes, cylinders, spheres and tori stay in those families. The alternative is the sharp offset, where faces move and their neighbours extend to meet. Edges stay sharp, but no single distance defines the result, and faces that do not meet after extension have no answer.

### 4.9 Thicken

| Parameter | Value |
|---|---|
| `body` | the body whose faces are thickened |
| `faces` | face references, a set `F` whose shared edges are tangent-continuous |
| `thickness` | `t` |
| `side` | `"front"`, along `n_F`, or `"back"`, against it |

**Definition.** The new solid is the union, over `q` in `F`, of the segment from `q` to `q + t n_F(q)` for `"front"`, or to `q − t n_F(q)` for `"back"`. The result is this solid alone. It is not united with the body.

**Names.** The faces of `F` keep their origins as the base. `{"offset": {"of": <origin>}}` names the far face of each. `{"rim": {"edge": <edge name>}}` names the side wall along each boundary edge of `F`, an edge between a chosen face and a face not chosen.

**Impossible.** `t <= 0`. A face set with a sharp edge between two chosen faces is not defined yet, see section 7.4.

- **A normal sweep, so rims are straight.** That is what CAD users expect of thicken, and it equals one-sided distance where `t` is less than the faces' radius of curvature on that side. The alternative is `{ p : d(p, F) <= t }` on one side. It is defined for any face set, but the rims become rounded.

### 4.10 Draft

| Parameter | Value |
|---|---|
| `body` | `S` |
| `faces` | face references |
| `pull` | `w`, the direction the part leaves its mould |
| `neutral` | plane `N`, by `point` and `normal` |
| `angle` | `θ` in degrees |

**Definition.** For each chosen face `F`, the **neutral curve** is `ν = S_F ∩ N`. At a point `x` of `ν`, split `n_F(x)` into its part along `ŵ = w / |w|` and its part `m` perpendicular to `ŵ`, and let `m̂ = m / |m|`. The drafted surface is

`D_F = { x + s (cos θ ŵ − sin θ m̂(x)) : x on ν, s real }`.

A positive angle tilts the face into the material on the `+w` side of `N`. A plane turns about its neutral line, and a cylinder parallel to `w` becomes a cone with half-angle `θ`. The result is the solid bounded like `S`, with each chosen face's surface replaced by `D_F`, and each neighbouring face extended or trimmed along its own surface to meet it. The face adjacency of `S` must stay the same.

**Names.** None made. Drafted faces keep their origins (rule 5).

**Impossible.** `|θ| >= 90`. `w = 0`. A neutral normal of zero. A face whose normal is parallel to `w` somewhere on its neutral curve, or whose surface misses `N`, because it then has no draft direction. A draft that would reverse the solid, so that the replaced boundary bounds no solid. A draft that swallows a neighbouring face is open, see section 7.5.

- **A neutral plane, not a parting line.** One plane defines every neutral curve, and both kernels compute `S_F ∩ N` exactly for analytic faces. The alternative is drafting about a chosen edge chain, the parting line. That follows the part's own edges, but the neutral curve is no longer planar, and the definition needs the chain's order.

### 4.11 Sweep

| Parameter | Value |
|---|---|
| `profile` | sketch step, region `R` |
| `path` | path step, `γ(s)`, `s` in `[0, L]` |
| `orientation` | `"rotation_minimising"` or `"fixed"` |

**Definition.** `M(s)` is a rigid motion with `M(0)` the identity that moves `γ(0)` to `γ(s)`. With `"fixed"`, `M(s)` is the translation by `γ(s) − γ(0)`. With `"rotation_minimising"`, its rotation `Q(s)` solves `Q' = [T × T'] Q`, `Q(0) = I`, where `T` is the unit tangent. The profile turns only as much as the tangent turns. Along a line `Q` is constant, and along an arc it turns about the arc's axis, so a sweep along lines and arcs is a chain of extrudes and revolves. The solid is the union of `M(s)(R)` over `s` in `[0, L]`.

The profile is carried rigidly from where it lies. It need not start in the plane normal to the path.

**Names.** `start`, the face `R`. `end`, the face `M(L)(R)`. `{"side": {"segment": s, "path": g}}`, the surface profile segment `s` sweeps along path segment `g`.

**Impossible.** Nothing beyond the profile's, the path's and section 1.5.

- **Rotation-minimising, not Frenet.** It is defined on straight segments, where the Frenet frame is not, and it does not flip at inflections. The alternative is the Frenet frame. It needs no integration, but it is undefined on lines and jumps where curvature changes sign.

### 4.12 Loft

| Parameter | Value |
|---|---|
| `sections` | at least two, each a sketch step `profile` and a `start` segment for each of its loops |
| `shape` | `"ruled"` or `"smooth"` |

**Definition.** Loops of different sections correspond by their place in the `start` lists. Corresponding loops are walked from their start segments in their own order. They must have the same number of segments, and segment `j` of one corresponds to segment `j` of the next. On each segment, `σ(t)` for `t` in `[0, 1]` is the point at fraction `t` of the segment's length from its start. For each corresponding chain `j` and each `t`:

- `"ruled"`: between sections `i` and `i + 1`, the line from `σ_ij(t)` to `σ_(i+1)j(t)`;
- `"smooth"`: the natural cubic spline through `σ_0j(t), σ_1j(t), ...` with section `i` at parameter `i`.

With two sections both shapes are the same. The solid is the region bounded by the first section's region, the last section's region, and these side surfaces.

**Names.** `start`, the first section's region. `end`, the last. `{"side": {"segment": s}}`, the side surface of the chain that starts at segment `s` of the first section. In a ruled loft through three or more sections, the spans of one chain are separate faces with the same origin.

**Impossible.** Fewer than two sections, or sections that do not name one start segment in each loop: the structural check catches both. Corresponding loops with different numbers of segments. Side surfaces that cross each other, so that they bound no solid.

- **Explicit start segments.** The correspondence is then the user's, and an edit that moves a section does not twist the loft. The alternative is matching by nearest points. That needs no input, but a small move can change the match and every side name with it.

Section 7.7 lists the open parametrisation question.

### 4.13 Pattern

| Parameter | Value |
|---|---|
| `body` | `S`, the seed |
| `layout` | `{"linear": {"vector": v, "count": n}}` or `{"circular": {"axis": ..., "count": n, "angle": φ}}` |

**Definition.** Copy `k`, for `k` from `0` to `n − 1`, is `T_k(S)`. Linear: `T_k` translates by `k v`. Circular: `T_k` turns by `k φ` degrees about the axis. The result is the regularised union of all copies.

**Names.** Copy 0 is the seed itself and keeps its origins. Copy `k >= 1` wraps them: `{"copy": {"instance": k, "of": <origin>}}`. Copies that touch along faces merge as usual.

**Impossible.** `n = 0`. A zero axis direction.

Choices:

- **A pattern copies a body, and a boolean applies it.** To pattern a cut, pattern its tool and subtract once. Every step keeps one definition, and the field kernel can repeat its domain for a pattern of tools almost for free. The alternative is a feature pattern that reruns an earlier step at each copy. That is closer to how SolidWorks users think, but every copy would need its references resolved again, which is a second naming problem.
- **An angle between copies, not a total angle.** No special case for 360 degrees. The alternative is a total angle with equal spacing. Changing the count would then keep the arc, but full and partial circles would need different spacing rules.

### 4.14 Mirror

| Parameter | Value |
|---|---|
| `body` | `S` |
| `plane` | `point` `a`, `normal` `m` |
| `keep_original` | `true` or `false` |

**Definition.** The reflection `M(x) = x − 2 ((x − a) · m / (m · m)) m`. It is rational for rational `m` of any length. The result is `S ∪* M(S)` when `keep_original` is true, else `M(S)`.

**Names.** Faces of `S` keep their origins. Faces of `M(S)` are `{"copy": {"instance": 1, "of": <origin>}}`.

**Impossible.** `m = 0`.

### 4.15 Lattice

| Parameter | Value |
|---|---|
| `region` | body `S` the lattice fills |
| `cell` | `"cubic"`, `"body_centred"`, `"face_centred"`, `"octet"`, or `{"custom": {"nodes": [...], "beams": [[i, j], ...]}}` |
| `size` | `a`, edge of the cubic cell |
| `frame` | pose `P` of the cell grid |
| `radius` | grade `ρ(p)` |

Preset cells, in unit-cell coordinates:

| Cell | Nodes | Beams |
|---|---|---|
| `cubic` | the 8 corners | the 12 cube edges |
| `body_centred` | corners and the centre | centre to each corner, 8 |
| `face_centred` | corners and the 6 face centres | each face centre to the 4 corners of its face, 24 |
| `octet` | corners and the 6 face centres | the `face_centred` beams and each face centre to its 4 neighbouring face centres, 36 |

**Definition.** The beams are the segments from `P(a (k + n_i))` to `P(a (k + n_j))` for every beam `(i, j)` and every integer vector `k`. Let `β(p)` be the distance from `p` to the nearest beam. The lattice is `{ p : β(p) <= ρ(p) }`, and the result is `S ∩* lattice`.

**Names.** `beams` for the beam surfaces, one origin for all of them. Faces of `S` that cut the lattice keep their origins.

**Impossible.** `a <= 0`. A radius of zero or less anywhere in the region. A beam whose two nodes coincide.

- **The radius is read at the query point.** The field is then one minimum of distances minus one grade, which is cheap and exact. The alternative is reading the radius at the nearest point of the beam, which gives true tapered beams. The geometry is nicer, but the field must find that point first.
- **One origin for all beams.** The PRD says freeform field surfaces have nothing to pick, and a part may hold a million beams. The alternative is a name per beam by cell index and beam. Each beam would then be pickable, but names would depend on the frame, and moving the frame would rename them all.

### 4.16 Gyroid and graded density

| Parameter | Value |
|---|---|
| `region` | body `S` |
| `period` | `a` |
| `level` | grade `c(p)` |
| `kind` | `"sheet"` or `"network"` |
| `frame` | pose `P` |

**Definition.** In local coordinates `(x, y, z) = P⁻¹(p)`, with `X = 2πx/a` and likewise `Y`, `Z`, let `g(p) = sin X cos Y + sin Y cos Z + sin Z cos X`. A sheet is `{ p : |g(p)| <= c(p) }`, and a network is `{ p : g(p) <= c(p) }`. The result is `S ∩*` that set.

**Grades.** A grade is a value that varies over space. It is used for a gyroid's level and a lattice's radius:

| Grade | Value at `p` |
|---|---|
| `{"constant": c}` | `c` |
| `{"linear": {"from": f, "to": g, "start": c0, "end": c1}}` | `c0 + τ (c1 − c0)` with `τ = clamp(((p − f) · (g − f)) / |g − f|², 0, 1)` |
| `{"distance": {"faces": [...], "near": c0, "far": c1, "span": D}}` | `c0 + clamp(d(p, faces) / D, 0, 1) (c1 − c0)`, faces of the region |

Graded density, as PRD section 4 means it, is a lattice radius or gyroid level given by a non-constant grade.

**Names.** `surface` for the gyroid surface. Faces of `S` keep their origins.

**Impossible.** `a <= 0`. A sheet with a level of zero or less anywhere, which is zero thickness. A linear grade with `f = g`. A distance grade with `D <= 0`.

- **A level band, not a true thickness.** The function is the definition, and both kernels hold it exactly (PRD section 3.1). The alternative is the set within `t/2` of the surface `g = 0`. Wall thickness would then be exact in millimetres, but the field would need the distance to an implicit surface, which only numerics can find.

### 4.17 Smooth blend

| Parameter | Value |
|---|---|
| `kind` | `"union"`, `"subtract"` or `"intersect"` |
| `target`, `tool` | bodies `A`, `B` |
| `radius` | `k` |

**Definition.** Let `a` and `b` be the exact signed distances to `A` and `B`, negative inside. The result is `{ p : f(p) <= 0 }` with:

- union: `f = max(k, min(a, b)) − |(max(k − a, 0), max(k − b, 0))|`;
- intersect: `f = min(−k, max(a, b)) + |(max(k + a, 0), max(k + b, 0))|`;
- subtract: the intersect formula with `b` replaced by `−b`.

For two planes meeting at a right angle this is exactly a rolling-ball fillet of radius `k`.

**Names.** `smooth` for the part of the boundary where the blend term differs from plain `min` or `max`. Elsewhere faces keep the origins of `A` and `B`.

**Impossible.** `k <= 0`.

The definition reads true distances. A field kernel that holds only a bound on distance, as after booleans, must compute true distance for `A` and `B` here. The B-rep can only fit splines (PRD section 4).

- **The circular blend.** It matches a fillet in the right-angle case, which gives the agreement checks one exact case to test. The alternative is the quadratic polynomial blend. It is cheaper and has no square root, but it matches no fillet anywhere.

### 4.18 STEP import

| Parameter | Value |
|---|---|
| `file` | `name` and `sha256` of the file's bytes |
| `pose` | placement |
| `decisions` | saved answers to import questions, by question id |

**Definition.** The solid the file's B-rep bounds, after applying the saved decisions, placed by `pose`. The file's units convert to millimetres exactly, for example one inch is exactly 25.4 mm. Per PRD section 3.6 the decisions are Jev's saved answers, such as whether a small gap is intended. Saving them in the recipe makes a rebuild give the same body offline.

**Names.** `{"imported": {"key": "#1234"}}`, the entity name of the face's `ADVANCED_FACE` in that file.

**Impossible.** A file that, after the decisions, bounds no solid: an open shell, or faces that cross. A file whose bytes do not match `sha256` is a recipe error, not geometry.

- **Names are the file's entity names.** The hash fixes the file, so its entity names are fixed too, and they are exact. The alternative is keys built from each face's geometry. Those survive a re-export from the source CAD tool, but two faces on the same surface collide. Section 7.10 covers replacing the file.

### 4.19 Mesh import

| Parameter | Value |
|---|---|
| `file` | `name` and `sha256` |
| `pose` | placement |
| `scale` | millimetres per file unit |

**Definition.** The closed, oriented triangle mesh bounds a polyhedron, and the solid is the region it encloses, placed by `pose` after `scale`. Every coordinate in the file is a binary or decimal number, which is an exact rational. PRD section 7 describes the reconstruction route: an AI model writes ordinary recipe steps from the mesh, and the app measures the deviation. The import step stays as the reference that deviation is measured from.

**Names.** `mesh`, one origin for the whole surface.

**Impossible.** A mesh that is not closed or not consistently oriented. Triangles that cross. `scale <= 0`.

- **A closed mesh is required.** Inside and outside are then decided exactly, by an integer winding number. The alternative is the generalised winding number, which accepts open scans. Real scans have holes, but the inside test would then involve solid angles and a chosen threshold.

## 5. On-disk form

A recipe is a JSON text, RFC 8259, in UTF-8:

| Key | Value |
|---|---|
| `format` | `"catallax-recipe/3"` |
| `name` | the part's name |
| `description` | optional |
| `steps` | list of steps; each is an object with `id`, `op` and the operation's parameters |
| `result` | id of the step whose body is the part |

A small recipe in canonical text:

```json
{
  "format": "catallax-recipe/3",
  "name": "plate",
  "steps": [
    {
      "id": "outline",
      "op": "sketch",
      "plane": {
        "rotation": ["1", "0", "0", "0"],
        "translation": ["0", "0", "0"]
      },
      "loops": [
        [
          {
            "kind": "line",
            "id": "bottom",
            "to": ["40", "0"]
          },
          {
            "kind": "line",
            "id": "right",
            "to": ["40", "30"]
          },
          {
            "kind": "line",
            "id": "top",
            "to": ["0", "30"]
          },
          {
            "kind": "line",
            "id": "left",
            "to": ["0", "0"]
          }
        ],
        [
          {
            "kind": "circle",
            "id": "bore",
            "centre": ["20", "15"],
            "radius": "4"
          }
        ]
      ]
    },
    {
      "id": "plate",
      "op": "extrude",
      "profile": "outline",
      "distance": "5"
    },
    {
      "id": "round",
      "op": "fillet",
      "body": "plate",
      "edges": [
        {
          "between": [
            {
              "step": "plate",
              "role": "end"
            },
            {
              "step": "plate",
              "role": {
                "side": {
                  "segment": "top"
                }
              }
            }
          ]
        }
      ],
      "radius": "1"
    }
  ],
  "result": "round"
}
```

`crates/catallax-recipe/tests/fixtures/every-operation.json` covers the original feature operations and roles. `tests/primitives.rs` covers stock-solid round trips, exact parameters, face references and singleton apex names.

**Reading.** A reader accepts any whitespace, any key order and any spelling of a number in the grammar. It rejects unknown keys, unknown operations and roles, ids outside the grammar, JSON numbers where an exact number belongs, edge names whose two origins are equal, and vertex names with no origins or repeated origins.

**Writing.** The same recipe always gives the same bytes. The canonical text:

- lists keys in the order this document's tables list them, which is the order of the fields in `catallax-recipe`, and leaves out optional keys that are absent or empty;
- writes every number in canonical form: a plain decimal when the value has a terminating one, `"0.0025"`, otherwise a fraction in lowest terms, `"1/3"`. `"-0"` is `"0"`;
- writes edge and vertex names in canonical order (section 3.3), and saved decisions in key order;
- indents by two spaces and puts each object member, and each element of an array of objects or arrays, on its own line. An array of strings or numbers stays on one line, as `["1", "0"]`. A colon and one space follow each key;
- ends with one newline.

**Versions.** The format name carries the version. Any change of meaning or schema makes a new version and an entry in [INTERFACE-CHANGES.md](INTERFACE-CHANGES.md). A reader reads only the versions it knows, and older ones go through an explicit conversion.

Version 2 added the five stock-solid operations and singleton singular-vertex names. Version 3 adds `extrude_through`; older operation meanings do not change. Reviewed version 1 or 2 files migrate by explicitly changing `format` to `catallax-recipe/3`, then parsing, checking and writing canonical text. The reader does not silently relabel older files or accept a version alias. The field evaluator accepts the current schema but still refuses sketches and feature operations with `field.unsupported-operation`.


Choices:

- **Unknown keys are errors.** A reader that dropped a key it did not know would build a different part without saying so. The alternative is ignoring unknown keys for forward compatibility. Older readers could then open newer files, but they would build them wrong.
- **An indented canonical text.** The recipe is the part file, and people diff it in Git. The alternative is one-line canonical JSON, such as RFC 8785. It is smaller and standard for hashing, but a one-number change becomes a one-line diff of the whole part.
- **Canonical numbers on writing, any spelling on reading.** Equal recipes then have equal bytes, and a hash of the text identifies the part. The alternative is keeping each number as the user typed it. `"1.10"` would survive a round trip, but equal parts could hash differently.

## 6. SolidDraft's kernel slot

SolidDraft is the application Catallax serves (PRD section 9.4). This section describes its present kernel slot and recipe in Catallax's own terms, and how each maps onto this language.

### 6.1 The slot

The slot is a set of operations on bodies the kernel owns. The application holds only opaque handles. Every operation returns a new body and leaves its inputs as they were:

- a planar sheet body from a closed profile, given as an outer loop and hole loops, each segment with both its ends;
- extrude by a positive distance along the profile normal or a given unit direction;
- revolve about an axis by an angle in degrees;
- union, subtract and intersect;
- fillet chosen edges with one radius, chamfer them with one width;
- shell with chosen faces left open, to one thickness;
- display, drawing and mass queries, which are outside this language.

The application tags each call with its own name for the operation. For each call the kernel reports which faces, edges and vertices it created, modified, deleted or split, with per-session ids. It also reports a role for each face it created: the profile region, the start or end cap of a sweep, a side swept from segment `j` of loop `i`, a blend of a given edge, an inner wall offset from a given face, or "other". Numbers are binary floating point, in millimetres and degrees.

The application's recipe is an ordered list of these calls, each with an id. Calls read earlier calls by id. Faces and edges are chosen by selectors evaluated against the body when the step runs:

- every entity;
- the entities a given step generated, followed through later splits and deletions;
- the one entity whose centre is nearest a point;
- the entities inside a box;
- the entities every one of several selectors chooses.

### 6.2 The mapping

| SolidDraft | Catallax |
|---|---|
| The application's tag for a call | The step id. Tags outside the id grammar, such as mixed case, are lower-cased, and a converter renames any that then collide. |
| A floating-point number | The shortest decimal that reads back as the same double. That is the number the user typed, not the double's long exact binary value. |
| A profile call that makes a sheet body | A `sketch` step, which makes a profile, not a body. The sheet's "region" face role has no counterpart. |
| Loop `i`, segment `j` | Segment id `l{i}s{j}`, for example `l0s3`. The outer-and-holes region equals the odd-nesting region for every valid profile. |
| A segment with both ends, joined within a tolerance | A segment with only its end. A converter checks each gap is within the application's tolerance and then makes the ends equal exactly. It moves an arc's end exactly onto its circle. Such a point always exists close by: the circle passes through the arc's rational start point, so rational points on it are dense, by stereographic projection from that point. |
| A sketch plane by origin, normal and x-direction | A pose. The world planes are exact. XY, with normal `+z` and `u` along `+x`, is the identity. YZ, with normal `+x` and `u` along `+y`, is `[1, 1, 1, 1]`. ZX, with normal `+y` and `u` along `+z`, is `[1, -1, -1, -1]`. A converter replaces any other plane by a rational quaternion whose frame is within tolerance. Rational rotations are dense, so one always exists. |
| Extrude, revolve | `extrude`, `revolve` with the same parameters. Catallax also defines negative distances. |
| Union, subtract, intersect | `boolean`. Catallax also defines a body combined with itself. |
| Fillet, chamfer, shell | `fillet`, `chamfer`, `shell`, choosing edges and faces by reference |
| Start cap, end cap, side from loop `i` segment `j` | `start`, `end`, `{"side": {"segment": "l{i}s{j}"}}` |
| Blend of an edge, where several blends meet "the nearest of them" | `{"blend": {"edge": ...}}`, and `{"corner": {"vertex": ...}}` where blends meet. Catallax has no nearest rule. |
| Inner wall offset from a face | `{"inner": {"of": ...}}` |
| "Other" | No counterpart. Every Catallax face has a defined role, and booleans make none. |
| Created, modified, split, deleted | Made (rule 1); kept or replaced (rules 2, 5); pieces keep the set (rule 3); gone (rule 7). The slot has no merge record. An adapter that reports a merge as one face modified and one deleted loses an origin, so it must keep the union instead (rule 4). |
| Seam edges, reported as a face listed twice | No counterpart. Catallax topology has no seam edges, so there is nothing to name. |
| Every entity | `"all"` |
| The entities a step generated | A face reference `{"face": {"step": "<step>"}}`. For the edges a boolean generated, edge references between the operands' patterns, such as `{"between": [{"step": "base"}, {"step": "boss"}]}`, one per pair of steps whose faces meet. |
| Nearest a point, inside a box, all of several | No direct counterpart, because these depend on where faces are, not what made them. A converter resolves them once against a built part and writes the names they chose, adding `near` with the old point only where names alone tie. |
| A step marked as expected to fail | No counterpart. A test that expects a step to fail says so outside the recipe. |

The slot's digest compares two builds by rounded geometry. Catallax compares builds by agreement on point sets and names (section 1.4) and needs no digest in the recipe.

## 7. Open items

These definitions are not settled. Each lists the options. Until one is chosen, the case is unfinished work in the PRD's sense, not an impossible input.

### 7.1 Sketch constraints with irrational solutions

A line tangent to two circles usually touches them at irrational points. Rational sketch coordinates can only come near the tangency, and near-tangent faces then reach the boolean engine as near-misses. Options:

1. Rational points only, with the application rounding. Simple, but tangency is lost.
2. Points defined by construction, such as "the tangent point from P to circle C", which both kernels solve exactly as algebraic numbers. This follows DESIGN R4, storing by definition.
3. The full constraint system in the recipe, solved by both kernels. The most general, and the most work to keep identical.

### 7.2 Fillet and chamfer ends, corners and overflow

The definition covers the run of an edge and an end where a third face meets the blend transversally. Open:

- **Corners where three or more blends meet at a vertex.** Options: a setback vertex blend with its own parameters; a sphere when equal radii meet at a convex corner of three mutually perpendicular faces, which is exact; the morphological opening restricted to the vertex's neighbourhood; refusing until one is defined.
- **Ends that meet a face tangentially, or end on a concave edge.**
- **Overflow**, where a contact curve leaves its face. Section 4.5 calls this impossible because no ball touches both named faces. Options: keep that; roll the ball onto the neighbouring face, as most CAD tools do; keep the neighbouring edge and trim the blend at it.
- **Different radii on edges that meet**, and variable-radius fillets.

### 7.3 Shell after booleans, and outward shells

Section 4.7 defines a shell by distance from the kept faces. That rounds the inner surface about concave edges and keeps the wall exactly `t` thick. A face-offset shell, extended until neighbours meet, keeps inner edges sharp. Its wall is `t√2` thick across a right-angled concave edge. The two agree on solids with no concave edges. Booleans and concave profiles make concave edges, so the choice shows on almost every machined part. Options: keep distance (constant thickness, and what the field kernel computes natively); adopt face-offset (what users of other CAD tools expect, but undefined where extended faces do not meet); offer both as a parameter.

An outward shell needs a rule at the open faces. Plain distance from the kept faces grows a rounded lip past the rim. Options: cut the lip off by each open face's surface; define the outward shell as an offset minus the original solid.

### 7.4 Thicken across sharp edges, and both sides

A normal sweep leaves a wedge-shaped gap at a convex edge between chosen faces and overlaps at a concave one. Options: fill the gap with a round, which is the distance definition; extend the far faces to meet; refuse face sets with sharp edges. A two-sided thicken needs names for two far faces. Options: an `offset` role with a side field, or two steps.

### 7.5 Draft that swallows a face

When a drafted face would pass a neighbour completely, the face adjacency cannot stay the same. Options: refuse; remove the neighbour, as CAD tools do; split the draft at the neighbour. Drafting both ways from a parting line inside one face is also open.

### 7.6 Sweep paths with corners, helices and splines

A path corner makes the frame jump. Options: a mitre, where profiles meet at the bisecting plane; a round corner, revolving the profile about the corner; refusal. Helices, which threads and springs need, and spline paths are not in version 1. Twist and scale laws along a sweep are not either.

### 7.7 Loft parametrisation and segment counts

Section 4.12 matches points at equal fractions of arc length. A line-to-line ruled surface is then a bilinear patch, an exact spline. Arc-length correspondence on arcs makes surfaces that are not NURBS. Both kernels hold them exactly by definition, and STEP would approximate them. PRD section 4 says "exact spline faces" for lofts, which needs a rational parametrisation of arcs instead. Open: which parametrisation; whether a smooth loft spaces sections uniformly, as now, or by chord length; and what to do when corresponding loops have different segment counts, where the options are refusing or splitting segments.

### 7.8 Up-to-face and one-sided through-all extrudes

Section 4.1.1 settles a two-sided through-all cut against a named target. An unbounded body-valued extrude, a one-sided through-all cut, and end conditions "up to face" or "up to next" remain open. The latter two need face references and a rule for several intersections.

### 7.9 Feature patterns

Section 4.13 patterns bodies. A feature pattern reruns a step at each copy and resolves its references again per copy. It needs a rule for what a reference means in copy `k`.

### 7.10 Replacing an imported file

A new version of a STEP file has a new hash and new entity names. Options: the application maps old names to new ones by geometric matching, with Jev choosing among measured candidates, and rewrites every reference; or the import step carries a rename table from old keys to new ones.

### 7.11 Open scans

See section 4.19: the generalised winding number would accept them.

### 7.12 Refreshing `near` points

After an edit moves a part, the stored `near` point is stale but still chooses. Options: the application rewrites it after each successful rebuild, which changes the recipe file on every edit; or leave it as written.

### 7.13 Which piece keeps a split segment's id

Section 3.5 leaves this to the sketch editor. A shared rule, such as "the piece containing the old start point", would make two editors behave the same.

### 7.14 B-rep feature geometry bugs

These are valid geometry, not additions to the impossible-input lists. The current B-rep feature implementation refuses them explicitly until the shared exact surface and frame representation is widened:

| Bug | Refusal code | Missing exact geometry |
|---|---|---|
| A signed extrusion whose normalized direction is irrational, such as `[1, 0, 1]`; a quarter-turn revolve about an axis with an irrational unit vector | `features.unsupported-algebraic-placement` | Algebraic placements and rotated frames |
| A non-cardinal rational-degree revolve, such as 45 degrees | `features.unsupported-algebraic-revolve-frame` | Algebraic rotated cap and trimming-curve definitions |
| An oblique extrusion of a circle or arc, such as direction `[3, 0, 4]` from the XY plane | `features.unsupported-elliptic-cylinder` | Elliptic cylinders, not circular cylinders |
| A cylinder generated with irrational radius, including an arc with radius squared 2 | `features.unsupported-algebraic-cylinder-radius` | Cylinders with exact algebraic radii or squared-radius definitions |
| A revolved circular arc producing a sphere with irrational radius | `features.unsupported-algebraic-sphere-radius` | Spheres with exact algebraic radii |
| A revolved arc or circle whose torus major or minor radius is irrational | `features.unsupported-algebraic-torus-radius` | Tori with exact algebraic radii |
| A valid trimmed horn or spindle torus sheet generated from an arc | `features.unsupported-spindle-torus` | Non-ring torus sheets and their branch predicates |

Extrusion supports axis-aligned and Pythagorean directions whose unit vectors are rational. Circular sweeps additionally require a normal direction and a rational surface radius. Revolve supports full turns and partial turns of 90, 180 or 270 degrees when the required frames and dimensions are rational. Zero degrees is impossible, not an empty partial revolve.

The language still specifies angles as rational degrees. A non-cardinal Pythagorean unit-circle point does not correspond to a rational number of degrees, so it has no exact input spelling in this schema. Supporting explicit rational cosine/sine pairs needs a separately versioned angle representation. No rounded degree value is treated as the exact Pythagorean angle.

Inherited cone/torus boolean gaps and unresolved exact decisions remain unfinished kernel work too. These feature limits do not justify approximate geometry, extra impossible-input rules, or a claim that the complete PRD is implemented.
