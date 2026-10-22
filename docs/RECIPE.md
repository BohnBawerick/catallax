# Catallax recipe language

Version 4 specification, format `catallax-recipe/4`. The released data model and evaluators still read version 3. Version 4 below is an implementation contract, not a claim that its schema or operations have landed. [PRD.md](PRD.md) sections 2, 4 and 10 step 3 ask for this document.

A recipe is the feature tree, the master model of a part. The B-rep kernel and the field kernel both build the part from it. This document gives each operation its parameters, one mathematical definition both kernels can be checked against, the names of the faces it makes, and the inputs that are impossible. It defines how names survive upstream edits, how a lost name fails, the file format, and how SolidDraft's present references map onto Catallax names.

The crate `catallax-recipe` holds the data model: the types, their canonical JSON text and the checks that need no geometry. `catallax-field` evaluates stock solids and booleans, with the certified queries described in [FIELD.md](FIELD.md). `catallax-features::evaluate` builds named B-rep stock solids, sketches, extrusions, revolutions, through-all cuts and booleans. Its current rational analytic representation leaves the legal-input bugs in section 7.14 open. The two evaluators do not yet cover the whole language.

The exchange contract in [JUDGE-INTERFACE.md](JUDGE-INTERFACE.md) is unchanged. Its question format stays the way the independent tester talks to the B-rep half. The recipe language sits beside it and does not replace it.

Where a definition is open, section 7 says so and lists the options. Every design choice in this document comes with its reason and one alternative with that alternative's cost.

## 1. The model

### 1.1 Steps and values

A recipe is a list of steps. Each step has an id and runs one operation. An operation reads the values of earlier steps by their ids and makes one new value. It never changes a value it reads.

A value is one of four kinds:

| Kind | Made by | Read by |
|---|---|---|
| profile | `sketch` | `extrude`, `extrude_through`, `extrude_to`, `extrude_next`, `revolve`, `sweep`, `loft` |
| path | `path` | `sweep` |
| body | solid operations below | the operations that act on solids |
| scan | `scan_import` | application reconstruction and deviation queries; no solid operation |

The recipe names one step as its `result`, which must make a body. Version 4 also permits a scan reference value, described in section 4.19.1; it is never a body or the recipe result.

- **Operations make new values and never edit their inputs.** The reason is that one definition per operation then reads as a function from inputs to output, and any step can be rebuilt or compared on its own. SolidDraft's kernel slot works the same way. The alternative is features that edit "the part" in place, such as an extrude that merges into whatever came before. Recipes get shorter, but every such feature hides a boolean whose naming must also be defined.
- **A step reads only earlier steps.** The recipe is then a sequence a rebuild walks once, with no cycles to detect. The alternative is a graph whose order the reader works out. Steps could then be listed in any order, at the cost of a cycle check and a second order to keep stable.

### 1.2 Ids

Step ids, and element ids inside a step (sketch segments, path segments, hole placements), use the grammar of the exchange contract's step names: 1 to 64 characters from `a-z 0-9 _ -`.

- **Faces are named through ids that the user's edits leave alone.** Changing a parameter never changes an id. Inserting or deleting another step or segment never changes one either. The alternative is positional indices, such as "segment 3 of loop 0". Those need no ids in the file, but inserting one segment renames the side face of every segment after it.

### 1.3 Numbers, units and placements

- Lengths are millimetres. Angles are degrees.
- Every input number is a JSON string in the exact rational grammar of piece 1, the grammar of the exchange contract (ND2): a decimal with an optional exponent, such as `"1.1"` or `"-2.5e-3"`, or a fraction, such as `"1/3"`. A JSON number where a geometric value belongs is an error. Counts and indices, which are not geometry, are JSON integers. Section 2.1 also admits a solved algebraic construction in a sketch coordinate or radius. Its coefficients use this same rational grammar; it does not add a rounded numeric spelling.
- A placement is a pose as in the exchange contract (ND1): a quaternion `[w, x, y, z]` of exact numbers, not all zero and of any length, then a translation. A sketch plane is the image of the local plane `z = 0`, with sketch coordinates `(u, v)` along local `x` and `y`. Its normal is local `+z`. That normal is a unit vector with rational coordinates, because a rational quaternion gives a rational orthogonal matrix.
- A direction, such as an extrude direction or a hole axis, is a vector of any non-zero length.

Choices:

- **Exact rationals, as in the exchange contract.** Both kernels then read the same number, and no value passes through a double. The alternative is JSON numbers. They are shorter, but every reader rounds `0.1` to a different binary value before the kernel sees it.
- **Angles in rational degrees.** A rational number of degrees is a rational multiple of pi, so its sine and cosine are algebraic numbers and the B-rep can hold them exactly. Their degree can be high: `cos 1°` has degree 48. The alternative is giving each angle as a rational point on the unit circle, such as `(3/5, 4/5)`. Every value then stays rational, but no user can type "30 degrees".
- **Poses by quaternion, not by origin, normal and x-direction.** A quaternion gives an orthonormal frame by construction, and every rational rotation matrix comes from an integer quaternion. The alternative is three vectors, as SolidDraft's slot uses. Users read those easily, but rational vectors are rarely exactly unit length and orthogonal, so the frame needs a tolerance check.

### 1.3.1 Algebraic evaluation of rational inputs

The 2026-10-22 architecture decision changes evaluated geometry, not the rational number input grammar. A rational degree angle is converted to certified algebraic sine/cosine values. Normalising a direction or computing an arc radius may also introduce algebraic values. Both evaluators keep these expressions exactly in the shared shape definitions; they must not reject a legal input merely because the result is not rational. An arc centred at the origin through `(1,1)` has radius `sqrt(2)` without requiring a new number syntax.

An oblique circular extrusion reduces to an elliptic-cylinder support with the original cap planes and finite sweep domain. A supported meridional revolve reduces to analytic supports with exact angular trims. The recipe and segment origins remain the definition of the operation after reduction. Unsupported trimmed revolution carriers remain implementation gaps, not additional impossible inputs. Generic derived-surface storage for future profiles follows the architecture decision when its full proof/query path is supplied.

The definitions in this document do not enlarge `Real` to arbitrary real-valued functions. Algebraic sketch solutions, rational-degree rotations, rational spline charts and algebraic offset/contact equations can use the shared algebraic foundation. A helix or the rotation-minimising frame of a general spline needs certified function and differential-equation definitions beyond `Real`. Their exact definitions remain required even while their implementation is deferred.

Every derived carrier must retain its inputs, branch selectors, parameter domains and orientation. Analytic reduction must preserve names and trims. It is complete only with certified charts, predicates, intersections, integration and export, as required by [ARCHITECTURE.md](ARCHITECTURE.md#exact-geometry-definitions-decision-of-2026-10-22). Missing carriers or proofs are implementation gaps, never impossible geometry.

The [whole-result error policy](ARCHITECTURE.md#composable-geometric-error-decision-of-2026-10-22) applies to every definition below. Model error needs a certificate strictly below `1/1_000_000` mm for the entire recipe result, with a transfer proof through each operation. A root isolating interval, display mesh or sampled comparison is not permission to perturb a definition. STEP approximation has its own bound and never replaces the model. Neither tolerance nor thread order may choose a root, a contact branch, a reference or a topology.

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

- **The recipe stores solved geometry, not constraints.** The kernels need positions, and constraint solving is the application's job. The alternative is storing the constraint system and letting both kernels solve it. That could preserve editing intent, but both kernels would need the same exact solver. The solved constructions below preserve tangency without that solver.
- **Each segment gives only its end point.** Closure needs no tolerance, and two segments can never disagree about a shared point. The alternative is both ends per segment, as SolidDraft does. That reads better on its own, but it needs a join tolerance, and a gap below it is a hidden approximation.
- **The region is "enclosed by an odd number of loops".** Islands inside holes need no extra structure, and loop order and direction do not matter. The alternative is an explicit outer loop with a list of holes. It is simpler to read, but an island inside a hole then needs a second profile.

#### 2.1.1 Solved algebraic sketch values

Version 4 permits a sketch `to` or `centre` coordinate, or a circle `radius`, to be either an existing rational string or a root definition:

```json
{ "root": { "coefficients": ["-2", "0", "1"], "index": 1 } }
```

The coefficients define `p(x) = c0 + c1 x + ... + cn x^n`. `index` is the zero-based position among its distinct real roots in increasing order. The example defines positive `sqrt(2)`. Coefficients are rational strings, the polynomial is nonconstant, and the selected root must exist. A writer removes trailing zeros, divides out repeated polynomial factors, clears denominators, divides the integer coefficients by their gcd and makes the leading coefficient positive. The index refers to the distinct roots before and after this normalisation. A zero or constant polynomial, or an out-of-range index, is a recipe error.

This is a solved construction, not an approximate coordinate or a constraint to solve during rebuild. The application chooses the branch and writes it. Kernels construct a certified algebraic leaf and retain its definition in `Real`; rational isolation intervals are caches. Different polynomial definitions of the same root need not have identical JSON. Canonical writing normalises the supplied definition, and geometric equality still needs proof.

For example, a line tangent to a unit circle at `(sqrt(2)/2, sqrt(2)/2)` can store each coordinate as the positive root of `2x² - 1`. Both incidence and tangency then hold exactly. All sketch validation and `side.segment` naming remain those of section 2.1. An inconsistent solved sketch is impossible under its existing rules. High algebraic degree or an unfinished root solver is not a new refusal. An underdetermined constraint system must be resolved by the application before it writes a sketch.

This keeps constraint editing outside both kernels while allowing every algebraic solved sketch. Rounding to rationals would lose tangency. Storing the whole constraint system would add solver branch choices to rebuilds.

### 2.2 Path

| Parameter | Value |
|---|---|
| `start` | point where the first segment starts |
| `segments` | list of path segments |

| `kind` | Fields | Meaning |
|---|---|---|
| `line` | `id`, `to` | straight segment to `to` |
| `arc` | `id`, `centre`, `to`, `axis` | arc about `centre`, counterclockwise about `axis` by the right-hand rule, to `to` |
| `helix` | `id`, `centre`, `axis`, `angle`, `rise` | rotation about the axis with simultaneous axial translation, as below |
| `spline` | `id`, `degree`, `control_points`, `weights`, `knots` | clamped rational B-spline, as below |

**Definition.** A path is the curve `γ(s)`, `s` in `[0, L]`, the segments joined end to end and parametrised by arc length.

**Names.** A path makes no faces. Its segment ids name sweep faces.

**Impossible.**

- a line whose end equals its start;
- an arc whose ends are not at exactly the same distance from its centre, or whose ends do not both lie in the plane through the centre perpendicular to `axis`;
- an arc whose end equals its start;
- a zero arc or helix axis;
- a helix with zero radius or zero angle.

Version 4 accepts corners. Every segment may have an optional `corner_axis`, last in its field order. It specifies the positive rotation axis at the join before that segment only when the incoming and outgoing unit tangents are opposite. It must be nonzero and perpendicular to those tangents. Section 4.11 defines the turn. A rotation-minimising sweep at a reversal without this axis has an incomplete recipe, not impossible geometry. A fixed-orientation sweep needs no such axis. A supplied axis at a non-reversing join, including the path start, is a recipe error.

For a helix starting at `p`, let `a` be the unit `axis` and `c` its `centre`. Require `(p-c) · a = 0`; otherwise the given centre and start do not define this helix and the input is impossible. For `u` in `[0,1]`, its curve is

`γ(u) = c + Rot_a(u angle)(p-c) + u rise a`.

`angle` is a nonzero signed rational number of degrees; `rise` is a signed length. Zero rise gives a circular path, including full turns. The endpoint is derived, so this segment has no `to` field. At rational endpoint angles it is algebraic. A nonzero-rise helix between those endpoints is not an algebraic curve; the zero-rise circle is.

A spline with degree `d >= 1`, controls `P_i` and positive weights `w_i` is `γ(u) = Σ N_i,d(u) w_i P_i / Σ N_i,d(u) w_i`. Use the Cox-de Boor basis on a nondecreasing rational knot vector over `[0,1]`, clamped with `d+1` copies of each end knot. There are `controls + d + 1` knots and one weight per control; interior knot multiplicities are at most `d`. There must be at least `d+1` controls. The first control equals the preceding endpoint exactly, and the last is this segment's endpoint. A wrong count, knot range, endpoint join or nonpositive weight is a recipe error. Repeated controls are allowed. A wholly constant segment is impossible as a path segment. At isolated stationary parameters use the one-sided tangent limits of the geometric curve, splitting its internal traversal there. If those limits reverse, a rotation-minimising traversal needs an explicit split into segments and a `corner_axis`.

An empty segment list has zero path length; its swept profile regularises to the empty body. Paths may cross themselves or repeat a point. The sweep's regularised result decides validity, not the path alone. Segment ids remain the naming source through internal knot spans and turn patches. No internal subdivision allocates a new user id. Twist and scale laws are not added by this revision; the operation remains the rigid sweep in section 4.11.

## 3. Names

### 3.1 Origins

An **origin** says which step made a face and which part of the step's operation made it:

```json
{ "step": "base", "role": { "side": { "segment": "top" } } }
```

Each operation's section lists the roles it makes. All of them:

| Role | Made by | Fields |
|---|---|---|
| `start`, `end` | block, cylinder, cone, extrude, extrude_to, extrude_next, revolve, sweep, loft; extrude_through makes only `start` when one-sided | none |
| `side` | stock solids, extrude, extrude_through, extrude_to, extrude_next, revolve, sweep, loft | `segment`; a fixed segment id for stock solids; for a sweep also `path` |
| `envelope` | sweep | `path` |
| `hole` | hole | `at` (placement id), `part` |
| `blend` | fillet | `edge` (edge name) |
| `bevel` | chamfer | `edge` |
| `corner` | fillet, chamfer | `vertex` (vertex name) |
| `inner` | shell | `of` (origin) |
| `round` | shell, offset, thicken | `edge` |
| `round_corner` | shell, offset, thicken | `vertex` |
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

#### 3.4.1 Stored near points

`near` is a model-space point and part of the recipe's definition. A rebuild never rewrites it, including after a successful upstream edit. The selected entity minimises Euclidean distance to its closed trimmed set, not to its supporting surface, centre or sampled mesh. Exact equality between the smallest distances is `recipe.ambiguous-reference`.

An editor may offer an explicit retarget edit. It resolves the old recipe, lets the user choose the intended entity in the edited recipe, and saves a new rational point only after proving it uniquely selects that entity. A tie remains a tie until the edit changes it. The edit and its undo record include the new point; rebuild order, caches and the last successful build have no effect. Automatic refresh was rejected because the same saved recipe could otherwise acquire different future meanings on two machines.

### 3.5 How names survive upstream edits

| Edit upstream | What happens to names |
|---|---|
| A parameter changes | Nothing. No name holds a dimension. |
| A step is inserted, or a step no name mentions is deleted | Nothing. |
| A step is renamed | The editor rewrites every name that mentions it. Names are structured, so the rewrite is exact. |
| A step is deleted that names mention | The structural check fails each such name (`recipe.unknown-origin`). |
| A sketch segment changes shape, for example a line becomes an arc | Nothing. The role holds the id, not the surface type. |
| A sketch segment is deleted | Names of its faces fail the structural check (`recipe.unknown-element`). |
| A sketch segment is split in two | The piece starting at the old directed start keeps the old id. Other pieces get fresh saved ids, under section 3.5.1. |
| New geometry divides a face | References without `near` choose every piece. References with `near` choose the nearest. |
| Two faces merge | A face reference to either origin finds the merged face. An edge between them is gone, so references to it fail. |
| A face disappears | References to it fail (`recipe.missing-reference`). |
| A pattern's count shrinks | Names of copies it no longer makes fail the structural check (`recipe.no-such-copy`). |

#### 3.5.1 Segment split identity

For a directed line or arc, use its traversal from its old start to its old end. Split at strictly interior, distinct parameters `0 < t1 < ... < tk < 1`. The interval `[0,t1]` keeps the old id. Each later interval receives a fresh id supplied by the edit and saved in traversal order. Its actual spelling is editorial; it must be unique within the sketch and satisfy section 1.2. Duplicate or endpoint split parameters are recipe-edit errors. Splitting makes no face by itself, and every swept piece gets its own resulting `side.segment` origin.

A circle has the canonical start `centre + (radius,0)` in sketch coordinates and counterclockwise traversal. Split it into at least two arcs. The arc beginning at that start keeps the old id if the start is a split point; otherwise the arc containing it in its relative interior keeps the old id. Store the resulting arcs in cyclic traversal order. The same rule applies to more than two cuts. A reversed or reshaped segment uses its direction immediately before the split, not a coordinate ordering or history from an earlier revision.

References to the old id select only the retained piece after this edit. The editor must explicitly rewrite a downstream selection to include new pieces when that is intended. Kernel-internal subdivision is different: it keeps the original segment origin on every piece under section 3.2. Allowing each editor to choose the retained id would make identical split edits rename different downstream faces.

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

Version 4 fillet, chamfer and draft also report the canonical consumed-face history defined in their sections. This is provenance, not a selectable origin on a surviving face. Empty results retain the history. Scan values report reference measurements under section 4.19.1 and never masquerade as named solid faces.

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
| `side` | version 4 optional `"both"`, `"forward"` or `"backward"`; defaults to `"both"` |

**Definition.** Let `T = { q + t w : q in R, t in the real numbers }`. The result is `S −* T`. With `side: "both"`, the tool extends in both directions from the sketch plane. Version 4 uses `t >= 0` for `"forward"` and `t <= 0` for `"backward"`. A one-sided blind cut is an ordinary signed `extrude` followed by `boolean` subtraction.

**Names.** Target faces keep their origins. Tool walls carry this step's `side` origins, with the sketch segment ids. A two-sided tool has no `start` or `end` origins. A one-sided tool has a `start` cap in the sketch plane, named by this step wherever it bounds the cut result, and no `end`. The surviving tool boundary is oriented outward from the result, as for boolean subtraction. Collinear or cocircular consecutive segments merge their origin sets.

**Impossible.** `w = 0`, or `w · n = 0`. Empty profiles leave the target unchanged; an empty target stays empty. The self-touching-shell rule still applies.

The B-rep implementation uses finite construction caps whose exclusion is proved from the target's defining bounds, not an arbitrary large length. Let the target lie within the sphere of radius `B` about model-space zero, and let `o` be the sketch origin. Choose

`k = ((B + |o|₁) |n|₁) / |w · n| + 1`.

For every target point, its line parameter `(p − o) · n / (w · n)` lies strictly between `−k` and `k`. Sweeping the shifted profile `R − k w` by `2 k w` therefore gives exactly the same cut inside `S`, with both caps outside it. Stock and swept definitions supply conservative bounds; union takes their maximum, intersection their minimum, and subtraction keeps the target bound. For a one-sided cut use `[0,k]` or `[-k,0]`; only its real cap at zero may survive. No artificial construction cap can survive on the result boundary.

#### 4.1.2 Extrude to a face or to the next boundary

Version 4 adds two body-valued operations. They do not implicitly add to or cut their target.

| Parameter | Value |
|---|---|
| `target` | body `S` supplying the limiting boundary |
| `profile` | sketch with region `R` and unit normal `n` |
| `direction` | optional nonzero vector `w`, defaults to `n` |
| `limit` | face reference, required for `extrude_to`, absent for `extrude_next` |

Resolve `limit` once against `S`. For `extrude_to`, let `H(q)` be the set of strictly positive `t` such that `q + t w` is on one of those closed trimmed faces. Choose its least element `h(q)`. Multiple selected faces are deliberate, and equal hits at one point are one hit. A tangent contact counts. For `extrude_next`, use the first strictly positive boundary hit at which the ray changes between inside and outside `S`. Mere tangency does not count. Starting on the target boundary does not stop the extrusion at zero; starting within material uses the first exit.

The result is `closure interior { q + t w : q in R, 0 <= t <= h(q) }`. Require a positive finite hit for every point in the relative interior of `R` and a uniform finite bound on `h` there. Boundary values come from closure. This defines a bounded solid even when the first-hit branch changes across the profile. Empty `R` gives the empty body after structural reference validation.

`start` names surviving pieces at `t=0`. `side.segment` names surviving walls over profile segments. `end` names the limiting cap and any curtains at discontinuities of `h` inside `R`. Those curtains are part of the closure boundary, not missing side faces. If a physical face contains more than one such contribution, it carries their union. The target contributes a limit, not its own origins: this extrusion creates every output face. Subsequent booleans use the usual naming rules.

A zero direction or `w · n = 0` is impossible. Missing reference matches are recipe errors. If positive hits accumulate at zero without a least positive hit, the ray has no first positive cap under this definition. A profile ray without the required hit, or hits with no uniform finite bound, cannot make the requested finite capped extrusion and is impossible. A proof or intersection method that has not yet found an existing hit is an implementation gap. Section 1.5 still applies to the resulting shell.

Using trimmed faces and the first positive hit makes cavities and several intersections deterministic. Using an untrimmed support could cap against an extension the user did not select. An unbounded body-valued extrusion is excluded by the bounded-body type; through-all remains an operation on a bounded target.

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
| `overflow` | optional `"roll_over"`, the default, or `{"keep_neighbours": limits}` under section 4.5.1 |

**Definition.** Take one chosen edge `e` between faces `A` and `B`. Call it convex where the material angle across `e` is less than 180 degrees, and concave where it is more. Let `σ` be the **spine**: the centres `c` at distance `r` from both `S_A` and `S_B`, on the material side for a convex edge and on the outside for a concave edge, in the component that runs along `e`. Each centre touches `S_A` at a foot point `a(c)` and `S_B` at `b(c)`. The **blend surface** is the union, over `c` in `σ`, of the shorter arc from `a(c)` to `b(c)` of the circle about `c` of radius `r` in the plane through `c`, `a(c)` and `b(c)`. The **blend region** `W_e` is the union, over the same `c`, of the part of that plane enclosed by the arc and the curves of `S_A` and `S_B` back to `e`.

The rolling ball `B(c, r)` touches `A` and `B` and rolls along `σ`. Where the edge ends at a vertex with a third face `C`, and the blend surface crosses `S_C` transversally, the blend continues along the spine's extension and is cut by `S_C`.

The result is `(S −* ∪ W_e over convex edges) ∪* (∪ W_e over concave edges)`. Every chosen edge uses the same `r`, and all edges are blended at once, so their order does not matter.

Exactness: two planes give a cylinder, a plane and a perpendicular cylinder give a torus. A blend between two curved faces has no closed STEP form. Both kernels hold it exactly by this definition, and only the STEP file approximates it (PRD section 3.1).

**Names.**

- `{"blend": {"edge": <edge name of e>}}` for the blend face;
- `{"corner": {"vertex": <vertex name>}}` for a patch where blends meet at a vertex, see section 7.2;
- `A`, `B` and `C` keep their origins.

**Impossible.** `r <= 0`, section 1.5, or incompatible explicit keep-neighbours constraints under section 4.5.1. Foot points leaving the trimmed faces are valid roll-over geometry, not proof of impossibility.

- **The rolling-ball definition, as the PRD states it.** It is a pure function of the two surfaces and `r`, and the field kernel can evaluate it exactly because the recipe names both faces. The alternative is the morphological opening `(S ⊖ B_r) ⊕ B_r`. It is fully defined at every corner, but it rounds unselected convex edges too and can remove thin material.

#### 4.5.1 Ends, corners and overflow

Version 4 uses `roll_over` by default. Continue the selected blend across neighbouring input faces when a contact foot reaches a trimming boundary. The neighbouring face can be trimmed away completely. The alternative `keep_neighbours` mode supplies explicit setback limits and preserves those neighbours. Neither mode may be chosen from a size heuristic.

Use a radius homotopy `ρ = λ r`, `0 < λ <= 1`, starting at the selected input edges. At a regular contact, centres satisfy `|c-a| = |c-b| = ρ`, `c-a` and `c-b` are normal to their contacting input supports, with the convex/concave material sides specified in section 4.5. Retain the contact sheets whose limit at zero is the selected edge and follow their oriented continuation. On a contact sheet, stay with the same pair of source faces until a foot reaches a trimmed boundary. At that event replace that foot's face by the adjacent input face, following the same continuous centre/foot limits and material sense. Do not start another edge run merely because it is connected to a selected vertex. If an event has several continuations with those limits, retain all and trim their hidden portions by the material-set rule below. Exact contact equations and domain predicates decide the event, never a sampled nearest face.

At a selected multi-edge junction, join the active contact sheets through their common higher-contact strata. A centre with three or more simultaneous contacts contributes the portion of its sphere whose directions lie in the spherical convex hull of its active contact normals on the appropriate material side. Here spherical convex hull means the unit vectors in their positive cone, restricted to the same local material sector as the incoming sheets. This supplies the vertex patch; it is not a fitted spline through trim points. At three perpendicular convex faces with equal radii, it is the sphere about the three inward offsets. Contact sheets and these spherical sectors are an envelope construction: retain only patches exposed on the boundary after all selected contributions have combined. Tangential contacts use the one-sided limits of the same sheets; coincident contributions merge. A vanished contact branch makes no new envelope beyond its collapse; the material already consumed or added by its affected region remains consumed or added.

For each active convex sheet, `W_e` is the swept region between the input boundary chain and the inward envelope; for a concave sheet it is the region between that chain and the outward envelope. Continue the oriented region as its contact faces change, including the vertex sectors above. Take the union of its affected regions over `0 < λ <= 1`, so material does not reappear when its last contact sheet collapses. The result remains `(S −* union W_convex) ∪* union W_concave`, all selections at once. Overlapping additions/removals therefore have the same order-independent set meaning as section 4.5. The boundary is the boundary of that regularised set, not every patch generated by the parametrisation. A fully consumed solid is empty. A parameter fold, a tangential end or loss of an input face is not itself impossible geometry.

In `keep_neighbours`, `limits` is a list of `{edge, at, setback}` records. `edge` is an edge reference resolving to one selected input edge, `at` is one of its endpoint vertex names, and positive `setback` is distance along that edge from the endpoint. The two ends may have different setbacks; a ring has no ends and cannot have such a record. Stop the run at each specified station and close its affected region against the plane normal to the edge there, keeping the side toward the remaining run. A limit on a nonsmooth station must choose a smooth subedge explicitly. The planes are part of the definition, not approximate construction caps. Every unselected neighbouring input face must survive with its whole trimmed point set and original origins. If the specified radius and limits cannot satisfy that preservation, or the two stations leave no run, the simultaneous constraints have no solution. Report that incompatibility; never reduce the radius, move a station or switch modes. Missing/ambiguous endpoint references are recipe errors. Arc-length stations on general curves may need certified derived definitions beyond algebraic `Real`.

Every exposed edge-run patch keeps `blend.edge` for the originally selected input edge, including after roll-over. Junction sectors and setback closures use `corner.vertex` for their original input vertex. Surviving input faces keep their origins; merged contributions union origins. The operation also reports `consumed`, sorted canonically: records with `source`, the origin set of an input face removed completely, and `by`, the origins of the blend/corner patches or affected regions that consumed it. If several input faces share the same origin set, coalesce identical records and union `by`. This history does not put a vanished source origin on a new face or make a stale face reference succeed. A fully consumed result still reports its history.

This version retains one constant radius per operation. Different radii use separate, explicitly ordered operations; each resolves its own input edges and preserves its own history. Variable-radius intent must be written as a new operation with an explicit radius law and compatible junction values before either kernel can evaluate it. It is not approximated by a sequence of guessed constant radii. This separates the settled constant-radius operation from a future schema extension, without declaring variable-radius geometry impossible.

Roll-over was chosen because reaching a neighbour does not make the requested geometry impossible. Explicit limits retain the local alternative when preserving that neighbour matters. Algebraic contacts fit `Real` and `RealPoly`; the complete continuation, canal, spherical-sector and trimming proof paths remain implementation work on deferred derived carriers. These equations and branch rules must be shared by both evaluators.

### 4.6 Chamfer

| Parameter | Value |
|---|---|
| `body` | `S` |
| `edges` | edge selection |
| `distance` | `d1`, on the face matching the first pattern of each edge reference |
| `other_distance` | optional `d2`, on the other face; `d1` when absent |
| `overflow` | optional `"roll_over"`, the default, or `{"keep_neighbours": limits}` |

**Definition.** Take one chosen edge `e` between `A` and `B`. The **contact curve** on `A` is the curve of `S_A` at distance `d1` from the curve of `e`, its intersection with the tube of radius `d1` about `e`, in the component next to `e`. The contact curve on `B` is the same with `d2`. The chamfer face is the ruled surface joining the point of each contact curve that lies in the normal plane of `e` at each point of `e`. `W_e` is the region between the ruled face and `e`. The result combines the `W_e` exactly as for a fillet, and ends meet a third face in the same way.

Exactness: two planes give a plane. A plane and a cylinder meeting at a circle give a cone.

**Names.** `{"bevel": {"edge": <edge name>}}` and `{"corner": {"vertex": <vertex name>}}`. A two-distance chamfer must list its edges, because `"all"` gives no first face. The structural check rejects that combination, and a rebuild fails an edge that matches both patterns both ways round.

**Impossible.** `d1 <= 0`, `d2 <= 0`, section 1.5, or incompatible explicit keep-neighbours constraints. Leaving a trimmed face triggers continuation onto the neighbouring face.

- **Setback measured in space from the edge.** It is the distance a user measures on the part, and on planes it equals the usual setback. The alternative is distance measured along each surface. That is the same on planes, but on a curved face it needs geodesics, which neither kernel computes exactly.

Use the continuation and regularised region rules of section 4.5.1 with distances `λ d1, λ d2` in place of radius. Each contact remains on its tube about the original selected edge and in that edge's normal section, continuing across a neighbour when it reaches a trim. At a vertex, close the active bevels by the planar facets of the convex hull of their endpoint contact positions within each local material sector. For a nonconvex fan, form these hulls separately in its incident convex sectors and combine their oriented affected regions with the same convex-removal/concave-addition rule. The sectors are the cells cut out by the input face tangent planes at that vertex, with the input material classification; keep all incident active sectors. A collapsed facet makes no face. This defines a faceted corner rather than an unspecified smooth interpolant. Tangent-plane limits are taken on the incident smooth face branches at a singular vertex.

`bevel.edge` continues to name an edge run across neighbours. `corner.vertex` names the junction facets and setback closures. Report consumed input faces exactly as for fillet. The explicit keep-neighbours mode has the same station planes and preservation constraint. For unequal distances, overlapping references must agree on which input face receives each distance; conflicting assignments are a recipe error, never resolved by list order. Missing contact, hull, continuation or trimming proofs remain implementation gaps.

### 4.7 Shell

| Parameter | Value |
|---|---|
| `body` | `S` |
| `open` | face references, possibly none |
| `thickness` | `t` |

**Definition.** Resolve `open` against the final input body, including boolean splits and merges. Let `K` be the union of the closures of the faces not selected as open. The result is `closure interior { p in S : d(p, K) <= t }`: the material within `t` of a kept face. Distance to the empty set is infinite, so opening every face gives the empty body. Regularisation removes zero-thickness remnants.

Points near an open face but far from every kept face are removed, so the cavity reaches through the open faces. The cavity is the set of points of `S` farther than `t` from `K`. At convex edges the walls' inner faces meet sharply. At concave straight edges and vertices of `K`, the inner surface is a cylinder or sphere of radius `t`. A curved edge gives a pipe envelope. The new inner boundary is the distance level `t` wherever a cavity remains. Opposing walls can consume the cavity entirely; this is a valid result, not a thickness failure.

**Names.**

- `{"inner": {"of": <origin>}}` for the inner wall at distance `t` from a kept face;
- `{"round": {"edge": <edge name>}}` for the inner surface at distance `t` from a concave edge;
- `{"round_corner": {"vertex": <vertex name>}}` at distance `t` from a concave vertex;
- kept faces keep their origins. What remains of an open face is its rim, the ring of the face within `t` of `K`, and it keeps the open face's origin.

**Impossible.** `t <= 0`.

- **Distance from the kept faces.** This is the PRD's description of the field side: "the field measures distance from the faces themselves". It gives exactly constant thickness. The alternative is offsetting each face and extending neighbours until they meet. That matches most CAD tools, but the walls are thicker at concave corners, and section 7.3 explains why the choice matters after booleans.

Distance semantics also apply after booleans. There is no face-offset variant hidden behind a topology change. A nearest witness in a kept face's relative interior contributes `inner.of` for each of its origins. Closest edge or vertex witnesses instead generate `round.edge` or `round_corner.vertex`. All generating witnesses contribute, and maximal merged faces union their origins under section 3.2. Naming uses the input body, never the order of distance queries.

An outward closed shell is already the composition `offset(S,t) −* S` for `t > 0`. For an open outward shell, first form that composition `C`, then apply `extrude_through` to `C`, or model a cutter `T` explicitly and form `C −* T`. The cutter defines the lip and opening, including which remote wall it reaches. The recipe stores these ordinary steps and their origins; there is no automatic clipping against an infinite open-face support. The outer offset inherits source names and adds its rounds, the inner boundary inherits the subtracted body's names, and the cutter supplies its wall names. References may therefore need `near` to distinguish inner and outer copies of one origin.

This settles outward-shell construction without changing `shell` to accept negative thickness. Extending `shell.open` to infer a cutter from arbitrary curved faces would hide another geometric choice. The explicit cutter is longer to write but gives both kernels the same lip, and works after booleans. Distance computations require trimmed-face metrics; a field residual is not a substitute. Curved offsets and pipes beyond analytic reductions need derived carriers.

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
| `body` | body `S` whose faces are thickened |
| `faces` | face references, their closed union `F` |
| `thickness` | `t` |
| `side` | `"front"`, toward the outside of `S`, or `"back"`, toward its interior |

**Definition.** Version 4 extends normal sweeps through sharp selected joins. In the relative interior of each smooth selected face, take the outward unit normal for `front` and its negative for `back`. At an interior edge or vertex of `F`, also take the unit directions `v` for which `q` is a local nearest point of `F` to `q + ε v` for every sufficiently small positive `ε`, with that point outside `S` for `front` and inside `S` for `back`. These are the local projection normal sectors. Take the closure of this relation in `(q,v)` to include the boundary of `F`; do not add new sectors at that boundary.

Let this closed relation be `G`. The result is

`closure interior { q + λ v : (q,v) in G, 0 <= λ <= t }`.

At a convex selected edge, the front sectors fill the gap by a circular round; at a concave edge the separate front sweeps overlap and their union trims the hidden part. Back reverses these cases. A boundary rim gets only limiting directions from the selected interior, so it remains the normal-sweep wall, not a rounded distance lip. Curvature folds and overlapping sweeps are regularised as sets, never rejected merely because two parameter sheets overlap. Section 1.5 applies to a non-manifold result. The result is this solid alone, without an implicit union with `S`.

**Names.** Surviving base patches keep the selected face origins. A surviving face-normal patch at `λ=t` has `offset.of` for each contributing origin. Exposed patches from sharp-edge sectors have `round.edge`; vertex sectors have `round_corner.vertex`, using input names. The sweep of a boundary edge of `F` has `rim.edge`. At a boundary vertex, a limiting sharp-edge sector can sweep an end fan that closes the gap between rim sheets. Give that fan the union of `rim.edge` origins of its incident boundary edges. It adds no spherical lip beyond the limiting sector. Singular envelopes formed before `λ=t` by a face sweep also use that face's `offset.of`; only exposed boundary portions receive names. Coincident contributions union their origins, and buried contributions have no face.

**Impossible.** `t <= 0`, plus section 1.5. A missing selection is a reference error. Sharp selected edges, excessive overlap and disappearance of a base face are not impossible-input rules.

A two-sided thicken is two steps on the same face selection, one `front` and one `back`, followed by boolean union. Their distinct step ids distinguish both far sides without a side field in `offset.of`. Thicknesses may differ. The zero-level base is internal and vanishes wherever both sides cover it. Each half uses the definition above; the boolean handles shared boundaries and origin unions.

Projection sectors preserve the existing normal sweep on smooth faces and straight rims while closing sharp-edge gaps. Extending offset supports until they meet would produce mitres with different thickness. Rounding the entire boundary of `F` would change the rim. Algebraic input faces give algebraic normal-sector and envelope equations, but general envelopes need certified derived carriers. A singular parameter chart alone is not evidence of impossible geometry.

### 4.10 Draft

| Parameter | Value |
|---|---|
| `body` | `S` |
| `faces` | face references |
| `pull` | `w`, the direction the part leaves its mould |
| `neutral` | plane `N`, by `point` and `normal` |
| `angle` | `θ` in degrees |
| `keep_regions` | optional list of exact model-space points overriding neutral-plane retention |

**Definition.** For each chosen face `F`, the **neutral curve** is `ν = S_F ∩ N`. At a point `x` of `ν`, split `n_F(x)` into its part along `ŵ = w / |w|` and its part `m` perpendicular to `ŵ`, and let `m̂ = m / |m|`. The drafted surface is

`D_F = { x + s (cos θ ŵ − sin θ m̂(x)) : x on ν, s real }`.

A positive angle tilts the face into the material on the `+w` side of `N`. A plane turns about its neutral line, and a cylinder parallel to `w` becomes a cone with half-angle `θ`. The result is the solid bounded like `S`, with each chosen face's surface replaced by `D_F`, and each neighbouring face extended or trimmed along its own surface to meet it. Section 4.10.1 defines material selection when the arrangement changes. Swallowed neighbours disappear; old adjacency is not a constraint.

**Names.** None made. Drafted faces keep their origins (rule 5).

**Impossible.** `|θ| >= 90`. `w = 0`. A neutral normal of zero. A face whose normal is parallel to `w` somewhere on its neutral curve, or whose surface misses `N`, because it then has no draft direction. A draft that would reverse the solid, so that the replaced boundary bounds no solid. Swallowing a neighbour is allowed under section 4.10.1.

- **A neutral plane, not a parting line.** One plane defines every neutral curve, and both kernels compute `S_F ∩ N` exactly for analytic faces. The alternative is drafting about a chosen edge chain, the parting line. That follows the part's own edges, but the neutral curve is no longer planar, and the definition needs the chain's order.

#### 4.10.1 Swallowed faces and retained regions

Build the arrangement of the drafted supports and the unchanged input supports, retaining their exact branch/domain predicates. Extend or trim their candidate patches at exact intersections. The candidate material cells are those bounded with the input material sense: crossing a surviving outward-oriented patch goes from material to exterior, and crossing it backwards goes from exterior to material. Cancel coincident opposite patches and union coincident same-sense origins. A vanished neighbour contributes no separating patch. Use winding number one for the oriented closed candidate boundary; cells of winding zero are exterior. A contradictory orientation, nonclosing boundary or self-touching shell is impossible. Failure to construct an existing closed arrangement is an implementation gap.

Let `C` be the resulting regular closed candidate material and let `C_j` be the closures of its connected interior components. Default retention keeps exactly those `C_j` whose intersection with the neutral plane contains a nonempty relatively open two-dimensional subset. Thus

`draft(S) = union* { C_j : interior_N(C_j ∩ N) is nonempty }`.

A point or edge contact with the neutral plane does not anchor a volume. When a drafted face consumes a neighbour, remove that neighbour and recompute the arrangement and components. All components that still reach the neutral plane survive, with no largest-piece or nearest-point tie-break. A component detached from the neutral plane disappears. If no component reaches it, the default result is empty; an off-body neutral plane therefore needs explicit kept-region picks when material is intended to survive.

When `keep_regions` is present, replace the neutral-plane test by `C_j` containing at least one listed point in its three-dimensional interior. The list is saved intent, not a point cache. An empty list explicitly chooses the empty result. A point on a boundary or outside all candidate material cells is a recipe reference error, and several points in one component select it once. Points may be retargeted only by an explicit edit. The override can retain remote components that automatic neutral-plane retention drops.

Drafted patches keep the origin sets of their source faces. Split neighbours keep their origins, consumed neighbours lose them, and merged patches union names. Report `consumed` with the same `source`/`by` convention as fillet; use the replacing drafted-face origins in `by`. The report preserves the removed face's history without making a lost reference resolve. Edges come from the final adjacent origin sets. The operation adds no undeclared cap.

A parting-plane draft with different angles on the two sides is an explicit composition. Intersect `S` with each half-space of `N` using finite cutters proved to cover `S`, draft each half with its signed angle and its own retention picks, then regularised-union the two outputs. The artificial cut faces disappear where the halves meet; any exposed remainder keeps its cutter origins. The recipe stores all these steps and their ordinary names. A nonplanar parting-line law is a different future operation, not inferred from this neutral-plane definition.

Automatic neutral-plane retention handles consumed faces without a new user pick. Explicit retained-region points provide the override when that default removes wanted material. Planar and axial cylindrical draft reduce to planes and cones on the algebraic foundation; general neutral curves, ruled supports, arrangement orientation and cell selection require the deferred derived-surface proof path.

### 4.11 Sweep

| Parameter | Value |
|---|---|
| `profile` | sketch step, region `R` |
| `path` | path step, `γ(s)`, `s` in `[0, L]` |
| `orientation` | `"rotation_minimising"` or `"fixed"` |

**Definition.** `M(s)` is a rigid motion with `M(0)` the identity that moves `γ(0)` to `γ(s)`. With `"fixed"`, `M(s)` is the translation by `γ(s) − γ(0)`. With `"rotation_minimising"`, its rotation `Q(s)` solves `Q' = [T × T'] Q`, `Q(0) = I`, where `T` is the unit tangent. The profile turns only as much as the tangent turns. For this rotation-minimising mode, `Q` is constant along a line and turns about an arc's axis along an arc, giving a chain of extrusion and revolution motions. Fixed-orientation translation along an arc is a different motion and may need a derived carrier. Write the motion explicitly as `M(s)(x) = γ(s) + Q(s)(x - γ(0))`. The solid is the regularised union of these profile positions and the corner positions defined below. Internal parameter overlaps are set union, not automatic self-intersection refusals.

The profile is carried rigidly from where it lies. It need not start in the plane normal to the path.

**Names.** `start`, the face `R`. `end`, the face `M(L)(R)`. `{"side": {"segment": s, "path": g}}`, the surface profile segment `s` sweeps along path segment `g`.

**Impossible.** Nothing beyond the profile's, the path's and section 1.5.

- **Rotation-minimising, not Frenet.** It is defined on straight segments, where the Frenet frame is not, and it does not flip at inflections. The alternative is the Frenet frame. It needs no integration, but it is undefined on lines and jumps where curvature changes sign.

#### 4.11.1 Corners and general paths

With fixed orientation, `Q=I` throughout, including corners, so the translated profile moves continuously. With rotation-minimising orientation at a corner `c`, let `T-` and `T+` be the one-sided tangents. If they agree, do nothing. Otherwise turn through `α = arccos(T- · T+)`, measured in degrees in `[0,180]`, about `a = (T- × T+)/|T- × T+|`. At a reversal use the outgoing segment's `corner_axis` and `α=180` instead. Include every position

`c + Rot_a(β) Q- (R - γ(0)), 0 <= β <= α`.

Continue the next segment with `Q+ = Rot_a(α) Q-`. This is a stationary round turn of the entire profile, not an inferred fillet radius on the path. The corner patch uses `side` with the generating profile segment and the outgoing path segment id. It makes no ordinary cap at the join. An internal spline knot uses the spline's id. Turning the profile can add material or overlap earlier material; the regularised union defines both cases. A zero-volume remnant disappears, and only the resulting shell is subject to section 1.5.

At the two path endpoints, keep `start` and `end` only where their profiles contribute actual boundary. A general oblique sweep can also expose an envelope generated by the profile interior, where the swept-volume map loses rank. Such a patch gets `{envelope: {path: g}}`, with `g` the generating path segment, or the outgoing segment at a corner. This role also covers an exposed interior profile at an internal join. Keep all generating origins on maximal merged faces. Internal folds that are buried in the union create no face. A closed path does not force the final frame to the initial frame. Any frame holonomy follows the differential equation, and any surviving endpoint patches keep those roles. Silently distributing an extra closing twist would define another part.

The stationary round turn gives one answer without a mitre plane that might fail to meet a remote or oblique profile. Mitred sweeps are not an alias for this operation. They would need a separate mode and clipping definition. Missing axes at reversals are explicit recipe errors, so coordinate axes or thread order never choose a half-turn.

Lines, circular arcs and their round turns can reduce to exact algebraic extrude/revolve supports where the profile permits; general oblique placements may still require derived carriers. Fixed sweeps of rational spline paths have algebraic rational charts. Helices need transcendental curve definitions. Rotation-minimising spline or helix sweeps need the specified frame ODE and certified transport, not sampled frames stored as geometry. A point evaluator or a polyline approximation is not that proof path. The algebraic foundation alone does not implement these general sweeps.

### 4.12 Loft

| Parameter | Value |
|---|---|
| `sections` | at least two, each a sketch step `profile` and a `start` segment for each of its loops |
| `shape` | `"ruled"` or `"smooth"` |

**Definition.** Loops correspond by their place in each section's `start` list, which must cover every loop once. All sections have the same number of loops. Walk each loop from its designated segment in stored order. Version 4 uses the following rational charts instead of fraction of arc length.

A line has its affine chart over `[0,1]`. Split a circular arc of directed angle `φ` into `ceil(|φ|/90°)` equal-angle spans, each at most 90 degrees. Each span uses a rational quadratic Bezier chart with endpoint controls `P0,P2`, middle control `P1` at the intersection of the endpoint tangents, and weights `1, cos(span_angle/2), 1`. Subdivision endpoints and weights are algebraic and retain exact branch selectors. A full circle starts at local `centre + (radius,0)`, goes counterclockwise and has four quarter-circle spans. Assign equal parameter intervals to the spans of each original segment.

For a loop with `m` original segments, assign its segment `j` the interval `[j/m,(j+1)/m]` of a common loop parameter `u`. Within that interval use the segment's chart above. Take the union of all segment and chart breakpoints from corresponding loops in all sections and refine each chart at those exact rational values. This permits unequal segment counts without changing any authored id. It makes the correspondence explicit and reproducible; it is not nearest-point matching or equal arc length.

Let `C_i(u)` be the resulting position in section `i`. On each common parameter interval:

- `"ruled"` uses `(1-v) C_i(u) + v C_(i+1)(u)` for `0 <= v <= 1` between adjacent sections;
- `"smooth"` uses the natural interpolating cubic spline in station parameter `v`, with value `C_i(u)` at integer station `i`, zero second derivative at the first and last stations, and continuous first and second derivatives at interior stations.

With two sections the two shapes agree. Uniform station spacing is part of the definition even when sections are unevenly placed. Chord-length spacing would change the solid when a section moved and would introduce another parametrisation choice. Rational section charts and this linear interpolation system give piecewise rational surfaces with algebraic coefficients. They can be represented as exact rational spline patches after a common-denominator refinement; no arc-length inversion or transcendental parameter is required.

Orient the first cap inward against increasing station, the last outward, and the side patches consistently with their loop's material side under the profile nesting rule. The body is the regular closed region bounded by those caps and side patches. A collapsed side patch can disappear. A self-crossing or inconsistent shell that cannot bound that region is impossible; an internal chart seam is not a physical edge. Zero-volume output regularises to the empty body. Inability to prove patch intersections is an implementation gap.

**Names.** `start` and `end` name surviving first and last caps. On an open common interval, `side.segment` uses the original segment of the first section whose parameter interval contains it. Later sections' extra breaks subdivide that origin but do not create names. Internal chart refinement preserves it. If several first-section segments contribute to a maximal face, that face carries all their side origins. An interval boundary alone has no area and invents no face. Ruled spans at different stations may be separate faces with the same origin.

**Recipe errors.** Fewer than two sections, unequal numbers of loops, a missing or repeated start segment, or a start not in its assigned loop. Unequal segment counts are valid. Profile failures and section 1.5 remain the geometric refusal rules.

Explicit start segments keep correspondence independent of geometric nearest-point guesses. Automatic chart refinement accepts different counts while making every split derivable from the input. The trade-off is that adding an authored segment changes its loop's uniform parameter allocation and can change the loft; editors must treat that as a geometry edit. Exact NURBS and general rational patch predicates, intersections, integration and export are still missing carriers, not part of the analytic foundation merely because their coefficients are algebraic.

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

#### 4.13.1 Feature pattern expansion

A feature pattern is an application edit that writes ordinary recipe steps, not another meaning of `pattern`. Given an ordered list of transformations `T_k`, the editor instantiates the feature's operation graph at each placement and saves every instantiated step, input id and reference. For a single target-changing feature, its mathematical result is `B_0 = S`, `B_(k+1) = F_k(B_k)`, where `F_k` is the saved operation for instance `k`. More complex feature graphs use the same earlier-step rule. The explicit step order matters whenever features interact.

Transform geometric inputs in model space by `T_k`, including points used in `near`; directions use its rotation. Rebind profile/path/body ids to the saved instance inputs. Each reference must explicitly name either an unchanged external origin or the corresponding creating step in that instance. Resolve it against the body that instance actually reads, using section 3.4. No implicit search for a geometrically similar face is allowed. `touching` patterns and nested origins are rewritten by the same map. If a transformed parameter has no spelling in that format's schema, expansion requires a later exact construction form; it must not round it to fit.

The editor allocates and saves a unique step id for each instance operation under section 1.2. Increasing the count preserves earlier ids; deleting instances follows the ordinary lost-name rules. Output faces carry the creating expanded steps' origins. They do not receive `copy.instance`, which belongs to copying a completed body. Thus a hole re-evaluated on a differently shaped target still has its own `hole.at` parts, and a fillet refers to the input edge it actually found.

An unresolved or ambiguous reference fails the expanded recipe as usual. Any geometry refusal is the underlying operation's refusal. Expansion has no extra impossible-geometry rule. A compact kernel-level feature-pattern opcode would need a new schema and the same explicit binding rules. Expansion costs file size; it avoids an undocumented per-copy binding language. It is buildable for exactly representable instances as those operations land. General circular placements can produce algebraic parameters outside the current non-sketch input schema; an exact construction extension is required before claiming all feature-pattern placements.

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

#### 4.18.1 Replacing an imported file

Replacement is an explicit recipe edit, never a rebuild side effect. The editor writes the new file name and content hash, resolves import decisions for those bytes, and rewrites every affected structured reference to the new `ADVANCED_FACE` keys in the same transaction. Decisions whose question ids depend on old bytes are not carried forward without validation. The old recipe and old bytes remain the undoable previous revision. A kernel evaluates only the new saved recipe and never guesses a correspondence.

A replacement map is an editorial relation from old face keys to sets of new keys. One-to-one matches rewrite directly. Several old faces may map to one new face; an old edge between them then vanishes and its reference must be edited or fail. A split maps one old face to several new keys. An editor may expand a list selection to those keys only after verifying the intended selection. A single reference or nested `of`/edge name that cannot express the split needs an explicit user selection or fails. Unmatched old keys never bind to a numerically identical key by accident: each match must be confirmed against the old and new file hashes.

The new import's origin is still `{step: import_id, role: {imported: {key: new_key}}}`. No alias is carried invisibly by the kernel. Candidate matching may use measured geometry, but closeness is evidence for an edit, not proof of identity or permission to alter geometry within 1 nm. The editor must not silently save an arbitrary member of a tied candidate set. A rename table inside the import was rejected because it would preserve stale keys indefinitely and need its own split/merge semantics. Ordinary rewritten references keep the saved recipe self-contained.

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

#### 4.19.1 Open scans as reference data

Version 4 adds `scan_import` with `file`, `pose` and `scale` in the same order as mesh import. Its value is the placed union `M` of the file's closed triangles, a reference surface with no inside or outside. Boundary edges, inconsistent triangle orientations, crossings and isolated samples are allowed. Degenerate triangles contribute their segments or points. Empty data is an empty reference. Coordinates remain exact file rationals under exact scale and placement; `scale <= 0` or a hash mismatch is a recipe error for this non-solid value.

A scan has no body faces, solid face origins or volume. Source items may be selected by their index in the hash-fixed file for reconstruction measurements, never as a solid face reference. `scan_import` cannot be the recipe result or an input to a solid boolean, shell, offset, fillet or other body operation. Those type mismatches are recipe errors, not geometric refusals. Replacing a scan uses the same explicit file/hash edit as STEP replacement; sample indices are not assumed to survive a replacement.

The application reconstructs a solid by writing ordinary recipe steps against the reference and evaluating those steps under their existing definitions. It reports measured deviation against the saved scan hash and recipe revision. For scan surface `M` and reconstructed body boundary `B`, report separate directed bounds

`d_scan_to_body = sup_(x in M) inf_(y in B) |x-y|`,
`d_body_to_scan = sup_(y in B) inf_(x in M) |x-y|`.

Report certified intervals only when the whole surfaces were covered with a proof. A finite sample report must state the sampled domain and maximum observed distances, not claim a Hausdorff bound. Direction matters: an intended closure over a hole in the scan can have a large body-to-scan deviation even when scan-to-body is small. Distance to an empty comparison set is infinite; two empty sets have zero directed deviation. A nonempty scan versus an empty reconstruction has infinite scan-to-body deviation.

The reconstruction's model error remains relative to its own exact recipe under ARCHITECTURE.md. Scan deviation is a separate measurement, not a claim that the unknown scanned object is correct within 1 nm. There is no winding threshold, automatic hole closure or inferred solid in `scan_import`. This preserves the observed data and makes every solid assumption visible in the reconstruction recipe. Reference storage needs no new surface carrier beyond exact triangles; certified deviation to reconstructed curved faces needs the corresponding closest-point proof paths.

## 5. On-disk form

A recipe is a JSON text, RFC 8259, in UTF-8:

| Key | Value |
|---|---|
| `format` | `"catallax-recipe/4"` |
| `name` | the part's name |
| `description` | optional |
| `steps` | list of steps; each is an object with `id`, `op` and the operation's parameters |
| `result` | id of the step whose body is the part |

A small recipe in canonical text:

```json
{
  "format": "catallax-recipe/4",
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

The example shows the version 4 specification. Until the schema implementation lands, version 4 files are not accepted by the crate. `crates/catallax-recipe/tests/fixtures/every-operation.json` covers the implemented version 3 feature operations and roles. `tests/primitives.rs` covers stock-solid round trips, exact parameters, face references and singleton apex names.

**Reading.** A reader accepts any whitespace, any key order and any spelling of a number in the grammar. It rejects unknown keys, unknown operations and roles, ids outside the grammar, JSON numbers where an exact number belongs, edge names whose two origins are equal, and vertex names with no origins or repeated origins.

**Writing.** The same recipe always gives the same bytes. The canonical text:

- lists keys in the order this document's tables list them, which the version 4 data model must follow, and leaves out optional keys that are absent or empty;
- writes every number in canonical form: a plain decimal when the value has a terminating one, `"0.0025"`, otherwise a fraction in lowest terms, `"1/3"`. `"-0"` is `"0"`;
- writes edge and vertex names in canonical order (section 3.3), and saved decisions in key order. Root definitions write `root`, then `coefficients` and `index`, with coefficients normalised under section 2.1.1;
- indents by two spaces and puts each object member, and each element of an array of objects or arrays, on its own line. An array of strings or numbers stays on one line, as `["1", "0"]`. A colon and one space follow each key;
- ends with one newline.

**Versions.** The format name carries the version. Any change of meaning or schema makes a new version and an entry in [INTERFACE-CHANGES.md](INTERFACE-CHANGES.md). A reader reads only the versions it knows, and older ones go through an explicit conversion.

Version 2 added the five stock-solid operations and singleton singular-vertex names. Version 3 adds `extrude_through`; older operation meanings do not change. Reviewed version 1 or 2 files migrate by explicitly changing `format` to `catallax-recipe/3`, then parsing, checking and writing canonical text. The reader does not silently relabel older files or accept a version alias. Both kernels evaluate the version 3 sketch profiles, signed extrusions, through-all cuts and full or rational-frame partial revolutions. Section 7.14 lists the specific unfinished algebraic feature cases.

Version 4 changes loft correspondence and segment-count handling, sharp-join thicken semantics, path corner acceptance and fillet/chamfer overflow semantics and draft region retention. It adds solved sketch roots, helix/spline paths, one-sided through-all, bounded face/next extrusions, reference-only scan values, blend overflow modes and draft kept-region picks. `shell` now states its regularisation and empty-kept-set case explicitly. The version bump prevents an old reader from silently producing the old arc-length loft or rejecting the newly defined geometry.

Do not relabel every version 3 recipe mechanically. Unchanged operations can migrate after reference validation. A loft containing arcs or unequal authored segmentation needs explicit review of correspondence; preserve the old version when its old solid is intended. Thicken edits involving sharp joins need review of the sectors and new `round` roles. The no-refresh `near` rule and split-ID policy govern new editing actions, not a rewrite of old saved ids. The interpreter and schema update will be a separate implementation change; version 3 fixtures, code and exchange-contract formats remain unchanged in this documentation decision.


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

### 6.3 Application entry

`catallax-app` owns the converter and the B-rep application entry. It has no build-time dependency on a SolidDraft checkout. The SolidDraft backend belongs in SolidDraft's repository.

- `soliddraft::convert(json)` reads `soliddraft-recipe/1` and returns a checked `catallax-recipe/3` plus an original-tag to step-id table. It retains the description, normalises tags without collisions, uses shortest round-trip decimal numbers, and gives profile segments their `l{i}s{j}` names. Rectangle arithmetic happens after conversion to exact numbers.
- Named and offset XY, YZ and ZX planes convert exactly. Full frames convert to rational quaternions, with exact checks that their axes stay within the source frame tolerance of `1e-9`. Loop joins use SolidDraft's `1e-7` mm tolerance before making shared endpoints exact. Arcs use rational stereographic projection onto their circles, with an exact endpoint displacement check.
- `build(&recipe, &mesh_options, &mut decider)` returns `ApplicationBody`. Its `built` contains the exact model, result body, intermediate bodies and immutable recipe definition. Its `display` retains the mesh answer's marks. Its `names` maps every result face to its complete origin set and every edge or vertex to its optional persistent name.
- `build_soliddraft(json, &mesh_options, &mut decider)` converts and builds in one call, returning the body and the tag table. Neither entry computes mass properties or writes STEP as part of rebuilding.
- `ApplicationBody::mass_properties(&accuracy, &mut decider)` requests proven volume, area, centroid and unit-density centroidal inertia for that result. `ApplicationBody::inertia` requests only the tensor. Both borrow the body; a measurement refusal leaves the build and display available. The tensor is a symmetric 3-by-3 enclosure matrix in model axes, with negative products of inertia off the diagonal and length units to the fifth power.

For picking, a display triangle's `face` indexes `names.faces` directly. An `EdgeMesh` identifies its B-rep edge and indexes the shared `Mesh::positions` to form that edge's polyline. Its closed polylines repeat their first index. Handle indices belong only to this returned model. A backend may allocate its own opaque ids, but must retain an explicit correspondence rather than cast indices or guess face order. Persist the recipe origins and edge/vertex names, not handles. Merged faces keep every creating origin, not just the first.

Failures distinguish conversion errors, feature-build errors, invalid topology and display refusals. `Error::code()` exposes their stable codes without parsing prose. A feature-build error retains its failing step and the original kernel refusal. Fillet, chamfer and shell return `ConversionError::NotYet` with `features.unsupported-operation`; their positional selectors are not guessed or discarded. An all-arc loop whose rounded endpoints cannot yet be closed on the same exact circles returns `soliddraft.arc-closure`. Uncertified frame or arc conversion returns a corresponding `soliddraft.*-conversion` not-yet code. These are implementation gaps, not new impossible-input rules. Expected-failure annotations are checked and removed, never treated as instructions to skip an operation or return an earlier body.

The example runs the same library entry and writes face-grouped display OBJ:

```sh
cargo run --release -p catallax-app --example build -- recipe.json display.obj
```

The bracket fixture exercises a flush plate/flange union, a boss and bore, four holes and a cross-plane slot. Its tests edit the flange height, pick the resulting faces through their triangles, retain merged origins, and check edge names and certified polylines. The application still needs its backend and UI integration before this repository's entry becomes a SolidDraft editing session.


## 7. Decision record and remaining open items

Items 7.1 to 7.13 are recorded here so their original questions remain traceable. A settled definition is not an implementation-coverage claim. The product choices for blend overflow, draft retention and scans are settled in their operation sections; remaining carrier and schema work is identified separately. Unfinished carriers, algorithms and proofs do not add impossible-input rules.

### 7.1 Sketch constraints with irrational solutions

Settled in section 2.1.1. The application saves solved algebraic roots with exact selectors. Both kernels read those definitions through the shared algebraic foundation; neither solves the application's constraint system or rounds tangencies.

### 7.2 Fillet and chamfer ends, corners and overflow

Settled for the constant-radius/distance operations in sections 4.5.1 and 4.6. Roll-over consumes neighbours by default and records consumed origins. Explicit setback limits provide keep-neighbours mode. Continuation, spherical fillet junctions and faceted chamfer junctions define ends and corners; hidden sheets trim by regularised material sets. Different radii use ordered steps. A variable-radius operation needs an explicit law and compatible junction schema, rather than changing this operation silently. General canal/corner carriers need certified derived surfaces.

### 7.3 Shell after booleans, and outward shells

Settled in section 4.7. Keep distance from the closed union of kept faces, including boolean splits and merges. Outward shells use offset minus the original and an explicit cutter for openings. Analytic reductions use the algebraic foundation; general trimmed-face metrics, offsets and pipes need further proof and carriers.

### 7.4 Thicken across sharp edges, and both sides

Settled in section 4.9. Local projection sectors fill sharp gaps, overlaps regularise, and boundary rims stay normal-swept. Two-sided material is the union of separately named front and back steps. Algebraic normal equations fit the foundation; general envelopes need derived carriers.

### 7.5 Draft that swallows a face

Settled in section 4.10.1. Retain candidate material components connected to the neutral plane by area, with explicit kept-region points as an override. Swallowed faces disappear with recorded history. Parting-plane draft composes two half-body drafts and a union. Nonplanar parting-line laws require a separate operation. Planar and axial cylindrical cases reduce analytically; general draft needs ruled carriers and certified arrangements.

### 7.6 Sweep paths with corners, helices and splines

Settled for rigid sweeps in sections 2.2 and 4.11.1. Corners turn the profile through the shortest stationary rotation, with a saved axis for a reversal. Helices and rational splines have exact path definitions. General rotation-minimising frames require certified transport; helices also require transcendental functions. These exceed the algebraic foundation. Twist and scale are not part of this revision's rigid operation.

### 7.7 Loft parametrisation and segment counts

Settled in section 4.12. Use rational arc charts, uniform integer section stations and exact common refinement of loop parameters. Unequal segment counts are valid. Names follow original first-section segments. Exact rational spline carriers and their proof paths must land after the analytic foundation. Version 4 replaces the version 3 arc-length meaning.

### 7.8 Up-to-face and one-sided through-all extrudes

Settled in sections 4.1.1 and 4.1.2. Through-all cuts may be one-sided. Bounded extrusions use the first positive trimmed-face hit or first material crossing for up-to-next. Missing finite caps are impossible end conditions; failure to compute existing hits is unfinished coverage. No unbounded body value is introduced. Analytic limit faces build on the foundation; arbitrary limits and discontinuity curtains may need derived carriers.

### 7.9 Feature patterns

Settled in section 4.13.1 as saved expansion into ordinary steps, with explicit per-instance bindings and order. Expanded operations create their own names. `pattern` keeps its body-copy meaning. No new surface kind is required beyond those needed by the expanded operations. Exact parameter serialization still needs extending for algebraic placements not expressible by existing non-sketch fields; rounding them is forbidden.

### 7.10 Replacing an imported file

Settled in section 4.18.1 as an explicit transactional edit of file hash, decisions and every affected reference. Geometric candidates inform an edit but do not establish identity. No kernel rename table or new surface kind is required.

### 7.11 Open scans

Settled in section 4.19.1. `scan_import` produces reference data, never an inferred solid. Reconstruction writes ordinary recipe steps and reports separate directed deviations against the scan. Closed mesh import is unchanged. Triangle reference storage is algebraic; certified deviation to curved reconstructions needs their distance-query proofs.

### 7.12 Refreshing `near` points

Settled in section 3.4.1. Rebuild never changes the saved point. Retargeting is an explicit recipe edit with the usual unique-nearest proof and tie error. Exact trimmed-face distance queries are required, not a new surface kind.

### 7.13 Which piece keeps a split segment's id

Settled in section 3.5.1. The old directed-start piece retains the id; circles use the canonical local positive-u point and counterclockwise traversal. Other ids are saved by the edit. Kernel subdivision keeps the original origin on all pieces and is not an authored split.

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
