# Piece 5: the boolean engine (`catallax-boolean`)

Status: complete for blind test writing and blind building.

This is the input the test writer and the creator get, together with `docs/DESIGN.md`, `docs/ARCHITECTURE.md`, `docs/JUDGE-INTERFACE.md` (the answer's topology) and the descriptions of pieces 1 to 4 in `docs/pieces/`. This piece builds on piece 3's `Model`, decides through pieces 1, 2 and 4, and uses their outcomes, facts and marks exactly as those descriptions define them. Nothing else is needed.

## Approved algebraic-definition migration

The 2026-10-22 architecture decision keeps the five boolean stages and their public result types. Widen internal coefficients, projection frames, probe geometry and germ calculations to piece 1 `Real`/`RealPoly` where definitions can be algebraic. Retain rational fast paths and construction facts. A failed conversion to `Rational` is not a geometric predicate or proof of impossibility.

Use piece 4's certified closure query when building rings, choosing cyclic order or requesting a full trace. Open conics use finite endpoints and their natural linear order. Keep cone nappes and hyperbola branches through intersection, imprint, classification and sewing. Polynomial carrier coincidence does not establish branch or trimmed-face coincidence.

Fix batch 5 joins the common foundation after its merge. Its contact-order and probe regressions must also cover algebraic inputs, without increasing global budgets or using sampled signs. Unsupported cases remain explicit implementation gaps; the wider types do not turn an unresolved branch into a bounded coincidence.

## What this piece is for

Piece 5 is the one boolean engine of the kernel (DESIGN.md, R7). It unites, subtracts and intersects two bodies of one `Model`, through five named stages that later operations (fillet, shell, offset) reuse:

1. **Intersect**: find where the two bodies' boundaries meet;
2. **Imprint**: split the faces and edges of both bodies along those meetings;
3. **Classify**: label every face piece inside, outside or on the other body;
4. **Keep**: choose the pieces the operation keeps;
5. **Sew**: join the kept pieces into a new body with maximal faces and edges.

It also answers whether a rational point is inside, outside or on the boundary of a body (`classify_point`) and gives a rational point strictly inside a face (`face_point`).

The operands never change. The result is a new body in the same model.

## Number types

| What | Type | Meaning |
| --- | --- | --- |
| Geometry of every entity | `SurfaceId`, `CurveId`, `PointId` in the model's `GeomStore` | Exact definitions (pieces 2 and 4) |
| Query points of `classify_point`, results of `face_point` | `RPoint3` | Exact rationals |
| Points used inside the stages | `PointId` and `Point3` | Exact real numbers |
| Marks carried by the result | `MarkLog` | The record of every not-exact decision behind the body |

This piece adds no number type and holds no tolerance.

## The bodies it works on

An operand is a `Solid` or `Empty` body (piece 3) that passes `Model::check`, with faces and boundary curves covered by piece 4. Its faces and edges need not be maximal, and it may have several lumps and voids. An operand of another body kind is refused with `boolean.not-solid`.

Cone faces use an apex-safe radial chart and exact parametric probes with certified chart paths. Probe frames and dimensions may be algebraic; rational definitions retain their native arithmetic path. Point loops survive imprint and sew when they remain isolated apex boundaries. Torus regions bounded by complete parallel or meridian circles use exact two-angle component labels, independent of arbitrary parameter seams. General noncircular torus boundaries still require a certified periodic chart. Missing intersection and chart paths are temporary coverage failures, not invalid operands.

Cylindrical winding uses a periodic diamond chart. Its dyadic rectangle corners keep phase as an integer numerator and positive denominator through quadrant selection, seam lifting and extrema comparisons. Only the two extrema become binary balls, through `Ball::from_ratio`. This avoids repeated general GCD normalization without changing the phase bounds, pole rejection, precision or winding decisions.

The result is a `Solid` or `Empty` body that:

1. passes `Model::check` with no fault;
2. holds exactly the material of the operation (below);
3. has maximal faces and edges (below);
4. follows the touching rules of `JUDGE-INTERFACE.md` (below);
5. carries the marks of the operation followed by the marks of the first operand, then of the second.

### The material of the operation

The material of a solid body is the closure of the open set it bounds. For operands `A` and `B` the result's material is the regularised set operation: the closure of the interior of `A ∪ B`, `A ∩ B` or `A − B`.
So material of zero thickness never survives: two blocks sharing a face intersect in an empty body, and a block minus a body that only touches it is the block.
When the material is empty, the result has kind `Empty` and no lumps.

### Maximal faces and edges

After sew, in the result:

- no edge separates two faces on the same surface (`same_surface`, piece 2) with the same sense;
- no edge has both its coedges in one face (a slit, such as a tangent segment or line left inside a face), and when such an edge is removed its vertices are removed too unless another edge uses them;
- no vertex joins exactly two edges that lie on the same curve (`same_curve`, piece 4) and bound the same two faces;
- an edge that covers a whole closed curve and whose vertex no other edge uses is a ring edge, with no vertex, unless that vertex is a singular point of the curve (a point piece 4 listed in the component's `singular`), which stays.

Two faces on the same surface that share no edge stay separate faces. When a decision behind a merge ends `NotExact`, the operation refuses with `boolean.merge-ambiguous`, stage `sew`: merging would need incidence facts that no proof gives (DESIGN.md, R5), and not merging would leave the topology not maximal (`JUDGE-INTERFACE.md`, ND6).

A face boundary that passes through a vertex twice (both branches at a node, for example) is one coedge loop, walked keeping the face on the left, never two loops that share the vertex. Piece 3's Euler count depends on it.

### Touching material

Where the result's material touches itself, the two sides of the contact share no face, edge or vertex (`JUDGE-INTERFACE.md`, ND4):

- pieces of material that touch each other only along curves or at points are separate lumps, each with its own faces, edges and vertices, which may use the same curve and point definitions;
- one connected piece of material that touches itself only at isolated points that are not where its edges meet is one lump, and those points add no entity;
- one connected piece of material that touches itself along a curve, or at a point where its edges meet, cannot be built as a manifold shell, and sew refuses with `boolean.self-touching-shell`.

Inside one shell, every edge has exactly two coedges of opposite sense, and the faces around every vertex form one cycle. Sew decides both from topology alone: an edge with four kept coedges, or a vertex whose kept faces form two or more cycles, is a contact. When the material on the two sides belongs to different lumps, sew gives each side its own copy of the edge or vertex; when it belongs to one connected piece, sew refuses.

A void is a shell inside another lump's outer shell. A void that touches the outer shell at points is a second shell of the same lump.

## The five stages

The stage functions are public so that later operations reuse them. Their outputs are opaque types read through accessors. Each stage adds to the model (new entities, new geometry) and changes nothing already there.

### Intersect

`intersect(model, a, b, dec)` returns an `IntersectionGraph`: the meeting of the boundaries of `a` and `b`, as section edges and isolated points.

- A **section edge** is an arc of a curve that lies on the boundary of both bodies. `SectionEdge { curve, start, end, faces_of_first, faces_of_second }` gives the curve, the end points (`None` for both on a whole closed curve, which is a ring), the faces of `a` the arc lies on (one face when the arc crosses a face's interior, two when it runs along an edge of `a`, more at a non-manifold operand edge) and the faces of `b`, each list in `faces_of_body` order.
- An **isolated point** is a point of both boundaries that is on no section edge: a tangent point of two faces, or a point where an edge of one body touches a face or an edge of the other.

Where two faces lie on the same surface, their common region is two-dimensional; the graph holds its boundary as section edges, which are arcs of each body's edges that lie inside the other body's face.

Section edges are split at every point where the set of faces they lie on changes, and at every vertex of either body that lies on them, so that each section edge lies on one fixed set of faces. They are otherwise as long as possible. Every curve is one piece 4 made (`intersect_surfaces`) or an edge curve of an operand; every end point is a point of the store.

`edges()` lists the section edges and `points()` the isolated points, in a deterministic order that is not part of the contract.

### Imprint

`imprint(model, graph, dec)` returns `Imprinted`: the faces and edges of both bodies split along the graph, as new topology in the model.

- Every face of each body is replaced by its **pieces**: faces on the same surface, with the same sense, whose union is the original face and which meet only along edges. A face the graph does not touch is its own single piece (the original `FaceId` may be reused).
- Every section edge is covered by edges of the imprinted topology, used by coedges of pieces of both bodies. These are the **common edges**. A section can be subdivided at other canonical curve cuts; its covering edges use the same partition as coincident operand edges.
- Every edge of each body is split at every section point on it.

The pieces satisfy the classification rule: **the interior of each piece, apart from the isolated points of the graph, lies entirely inside the other body, entirely outside it, or entirely on its boundary**. This is what imprint is for, and tests check it through `classify` and `face_point`.

Accessors: `pieces_of_first()` and `pieces_of_second()` list the pieces of each body, and `common_edges()` the common edges, in deterministic orders that are not part of the contract.

### Classify

`classify(model, imprinted, dec)` returns a `Classification`. For each piece `P` of one body, `class_of(P)` gives its class against the other body `X`:

| Class | Meaning |
| --- | --- |
| `Inside` | the interior of `P`, apart from isolated points of the graph, lies in the interior of `X` |
| `Outside` | the interior of `P`, apart from isolated points of the graph, lies outside `X` |
| `OnSame` | `P` lies on the boundary of `X`, in a face of `X` whose outward normal agrees with `P`'s |
| `OnOpposite` | `P` lies on the boundary of `X`, in a face of `X` whose outward normal is opposite to `P`'s |

`class_of` returns `None` for a face that is not a piece of either body.
The method is the creator's choice. A method that fits: classify the pieces next to each common edge by the radial order of the four faces around it, decided exactly from their normals and tangents (with second-order terms where faces are tangent), carry each class across edges that are not common edges, and classify any group of pieces that touches no common edge by `classify_point` of one of its rational points.
A rational probe that is an isolated point of the graph gives `ProvenCoincident` and decides no class. A method using such probes must choose another point rather than assigning an on-boundary class from that point.

### Keep

`keep(classification, imprinted, op)` returns a `Selection` by this fixed table. It makes no decision.

| Operation | Pieces of the first body | Pieces of the second body |
| --- | --- | --- |
| `Unite` | `Outside`, `OnSame` | `Outside` |
| `Intersect` | `Inside`, `OnSame` | `Inside` |
| `Subtract` | `Outside`, `OnOpposite` | `Inside`, with sense `Reversed` |

`faces()` lists `(piece, sense)`: the pieces of the first body in `pieces_of_first` order, then those of the second in `pieces_of_second` order, each with `Sense::Same` except the second body's pieces in a subtraction.

### Sew

`sew(model, imprinted, selection, dec)` builds the result body from the selected pieces:

1. a piece with sense `Reversed` becomes a new face on the same surface with the opposite sense and its loops travelled the other way;
2. the pieces are joined along their shared edges into shells, applying the touching rules above;
3. faces and edges are merged until they are maximal, and redundant vertices removed;
4. shells are grouped into lumps: each outer shell starts a lump, and each void shell joins the lump whose outer shell encloses it;
5. the body is added with kind `Solid`, or `Empty` when no piece was selected.

Sew records an `Incident` fact, with origin `FactOrigin::Construction("boolean.sew")`, for every vertex point it puts on an edge's curve and for every edge curve it puts on a face's surface, unless the fact already holds, so that the result passes `check`'s `topology.missing-incidence`. Each such fact is guaranteed by the definitions or by a `Proven` or `ProvenCoincident` decision; sew never records one on a `NotExact` (DESIGN.md, R5).
The order of lumps, shells and faces in the result is deterministic but not part of the contract.

### The whole operation

`boolean(model, a, b, op, dec)` runs the five stages in order and returns `BooleanAnswer { body, marks, report }`, or `BooleanRefusal { stage, refusal }` naming the first stage that broke.

- `marks` holds the marks of the five stages in order, then the marks of `a`, then those of `b`. The result body's `marks` field holds the same log.
- `report.stages` holds one `StageRecord` per executed stage, in order, with the statistics of the decisions the stage made and the number of marks it added. Each counter of a record's `stats` other than `max_bits` is what the stage added to the `Decider`'s statistics: that counter of `dec.stats()` when the stage ended minus that counter when it began. A record's `max_bits` is `dec.stats().max_bits` when the stage ended. Piece 1 keeps only the running maximum, so this bounds the stage's own highest precision from above and equals it when no earlier decision went higher.
- With an `Empty` operand, `boolean` may skip the stages and follow the set rules: `A ∪ ∅ = A`, `∅ ∪ B = B`, `A ∩ ∅ = ∅ ∩ B = ∅`, `A − ∅ = A`, `∅ − B = ∅`. Its report then holds no stage records. A result equal to an operand may reuse that operand's lumps.

When both operands are the same `BodyId`, `A − A` produces a new empty body without geometric decisions. `A ∪ A` and `A ∩ A` retain the material and run maximal normalisation, so a nonmaximal operand does not bypass the result invariant. Their reports contain only the executed sew stage. All three rules retain the solid-or-empty input check and the normal operand mark order. Recipe evaluation uses these engine rules; the legacy question exchange retains its existing repeated-operand mark and statistics contract.

## Points and bodies

`classify_point(model, body, point, dec)` decides where a rational point lies:

- `Proven(Containment::Inside)` when it is in the interior of the body's material;
- `Proven(Containment::Outside)` when it is outside the material;
- `ProvenCoincident` when it is on the boundary: on a face, an edge or a vertex;
- `NotExact` when the point could not be separated from a face of the body; its bound is a `Distance` bound: the point lies within that distance of the face;
- `Refused` with `boolean.classify-ambiguous` when no decision path could settle it within the budget, or with a refusal from pieces 1, 2 or 4.

An `Empty` body gives `Proven(Outside)` for every point. The method is the creator's choice; a method that fits casts rational rays from the point and counts crossings, and decides whether a ray's hit lies in a face by probes on the face's surface. Because the point is rational, every ray and probe can be chosen with rational data, so piece 4 can intersect them exactly. That is why the query point is an `RPoint3`. It takes `&mut Model` because it may add auxiliary planes, curves and points to the store; it never changes topology.

`face_point(model, face, dec)` returns a rational point strictly inside a face: on the face's surface, in the face, and on none of its edges or vertices. Rational points are dense on every surface the builders make, because their frames are rational and orthonormal. A face on a hand-built surface with no rational point near it is refused with `boolean.no-rational-point`. It never carries marks: a decision that ends `NotExact` while proving the point inside the face makes it try another candidate, and when none is left it refuses with `boolean.no-rational-point`.
The point may also lie on another face, including at an isolated touching point or on a tangent line between touching lumps, when that contact is not an edge or vertex of this face. A caller needing a positive distance from every other face checks the point and finds another when that check fails.

`PreparedFaceQueries::new(&mut model)` provides repeated exact trimmed-face
membership within one `Decider`, numeric policy and operation-budget scope.
It reuses certified boundary traces, projected interval boxes and spherical
reference windings. Every query still proves exact support and closed-trim
membership. A homotopy proof must separate the query from the enclosed boundary;
points on a closed rim use exact incidence checks. Exact surface probes remain
available when the chart proof is unresolved. Numeric refusals do not reset the
allowance.

The projected-box key includes the exact model surface and chart choice, edge,
trace piece, parameter interval and precision. The cache retains at most 4096
entries; a full cache computes further boxes without storing them. Cancellation
is checked on hits. Drop the batch before leaving its budget scope or changing
existing geometry or topology. Auxiliary query geometry may be appended.
Ordinary construction and single-point queries retain their previous proof order.


## Decisions, outcomes and marks

- Every geometric decision goes through piece 2, piece 4 or the `Decider`. A rounded display coordinate is not an exact predicate. No angle enters a decision (DESIGN.md, R8).
- A `NotExact` outcome in intersect, imprint or sew is accepted with a mark built with `Mark::from_not_exact`, with the bound kind of the piece that returned it, only when its degenerate branch needs no fact that is not already true: a contact treated as tangent, and a hit that piece 4 merged or made tangent (its facts hold). Every other degenerate branch (two points treated as one point, a point treated as on a curve or surface, two curves or surfaces treated as the same, and a hit piece 4 kept on a branch, which carries no incidence on its curve) would need sew to record an `Incident` or `Coincident` fact that no proof gives. DESIGN.md R5 forbids that, so the operation refuses: with `boolean.imprint-ambiguous` in intersect and imprint, and with `boolean.merge-ambiguous` in sew. It never picks a side from cached coordinates.
- Probe constructions retain their exact point/curve-to-plane incidences, and a `ProvenCoincident` side query may supply a point-to-surface fact with that proof's provenance. An unmarked, complete surface-pair answer containing one curve, together with proven incidence on both supports, proves incidence on that curve. With multiple components, membership of the supports does not choose a branch. Auxiliary probe planes and their derived sections remain valid predicate facts, but do not seed subsequent probes; curves used as inputs to construct a probe plane remain eligible seeds.
- In classify, a containment that cannot be decided is a refusal with `boolean.classify-ambiguous`: there is no degenerate class to take.
- Wherever two edges meet at a vertex with the same tangent, imprint first compares exact curvature. If curvature also agrees, regular germs may be compared by exact Taylor coefficients through order eight, using one transverse linear parameter for both curves. Two quadric sections of a common quadric have contact order at most eight unless they share a component. No nearby sample point decides the order. If the available branch derivatives do not prove a strict order, imprint refuses with `boolean.imprint-ambiguous`.
- A `Proven` or `ProvenCoincident` outcome is never wrong (DESIGN.md, R2).
- The same model, operands and `Decider` state give the same result, the same new handles in the same order, the same marks and the same statistics on every machine and thread (DESIGN.md, R10).

## How this piece uses the earlier pieces

- **Outcomes.** Stages return `Answered`; `boolean` wraps the first refusal in a `BooleanRefusal` with its stage. `classify_point` returns an `Outcome`.
- **Budgets.** Every decision runs on the caller's `Decider`. This piece never creates a `Decider` or a policy, and it reads statistics only through `stats()`, never `take_stats()`, so the caller's running totals are kept. The deadline is honoured through the decisions it makes; a decision after the deadline refuses with `numbers.deadline`, and the boolean stops in that stage.
- **Facts.** Every point and curve made by piece 4 comes with its incidence facts. This piece records the `Incident` facts listed under sew, and may record `Coincident` facts with origin `Proof` for points and curves it proves the same. Piece 2's and piece 4's decisions read them first.
- **Marks.** This piece builds marks from the `NotExact` outcomes it accepts, keeps their codes, and passes them to the result in the order stated.
- **Refusals.** This piece adds the codes below and passes every refusal of pieces 1, 2 and 4 through unchanged inside a `BooleanRefusal`.

## Refusal codes

| Code | Stage | When |
| --- | --- | --- |
| `boolean.not-solid` | intersect | an operand's kind is neither `Solid` nor `Empty` |
| `boolean.self-touching-shell` | sew | one connected piece of the result touches itself along a curve, or at a point where its edges meet |
| `boolean.classify-ambiguous` | where `classify_point` was called (usually classify) | a piece's class, or a point's containment, could not be decided |
| `boolean.imprint-ambiguous` | intersect or imprint | edges meeting at a vertex could not be ordered around it, or a decision ended `NotExact` where the degenerate branch would need a fact no proof gives |
| `boolean.merge-ambiguous` | sew | a decision behind a merge of faces, edges or vertices ended `NotExact` |
| `boolean.no-rational-point` | where `face_point` was called (usually classify) | no rational point was found on a face's surface where one is needed |
| `boolean.torus-chart-required` | where face membership or projection was called | a torus boundary is not covered by the complete-circle periodic region path |

Refusals of pieces 1, 2 and 4 keep their codes. The stage of a refusal from `classify_point` or `face_point` called inside a stage is that stage. Called directly, they return the code alone.

## Programming errors

These may panic. None can come from a question file:

- a `BodyId`, `FaceId` or other handle from another model, or out of range;
- an operand that fails `Model::check`;
- the same `BodyId` passed as both operands;
- passing to `imprint`, `classify`, `keep` or `sew` an output of an earlier stage made on another model or on bodies changed since (bodies never change, so this means a value from another model).

## Exact public interface

The crate re-exports these types at its root. Everything listed is `pub`; fields not listed are private. The compiled crate in `crates/catallax-boolean` is authoritative if a signature listing becomes stale; update this contract with any API change.

```rust
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum BoolOp { Unite, Subtract, Intersect }
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Stage { Intersect, Imprint, Classify, Keep, Sew }
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct StageRecord { pub stage: Stage, pub stats: DecisionStats, pub marks: usize }
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct StageReport { pub stages: Vec<StageRecord> }
#[derive(Clone, Debug)]
pub struct BooleanAnswer { pub body: BodyId, pub marks: MarkLog, pub report: StageReport }
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BooleanRefusal { pub stage: Stage, pub refusal: Refusal }
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Containment { Inside, Outside }
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum FaceClass { Inside, Outside, OnSame, OnOpposite }

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SectionEdge {
    pub curve: CurveId,
    pub start: Option<PointId>,
    pub end: Option<PointId>,
    pub faces_of_first: Vec<FaceId>,
    pub faces_of_second: Vec<FaceId>,
}
#[derive(Clone, Debug)]
pub struct IntersectionGraph;
impl IntersectionGraph {
    pub fn first(&self) -> BodyId;
    pub fn second(&self) -> BodyId;
    pub fn edges(&self) -> &[SectionEdge];
    pub fn points(&self) -> &[PointId];
}
#[derive(Clone, Debug)]
pub struct Imprinted;
impl Imprinted {
    pub fn pieces_of_first(&self) -> Vec<FaceId>;
    pub fn pieces_of_second(&self) -> Vec<FaceId>;
    pub fn common_edges(&self) -> Vec<EdgeId>;
}
#[derive(Clone, Debug)]
pub struct Classification;
impl Classification { pub fn class_of(&self, face: FaceId) -> Option<FaceClass>; }
#[derive(Clone, Debug)]
pub struct Selection;
impl Selection { pub fn faces(&self) -> Vec<(FaceId, Sense)>; }

pub fn boolean(model: &mut Model, a: BodyId, b: BodyId, op: BoolOp, dec: &mut Decider) -> Result<BooleanAnswer, BooleanRefusal>;
pub fn intersect(model: &mut Model, a: BodyId, b: BodyId, dec: &mut Decider) -> Answered<IntersectionGraph>;
pub fn imprint(model: &mut Model, graph: &IntersectionGraph, dec: &mut Decider) -> Answered<Imprinted>;
pub fn classify(model: &mut Model, imprinted: &Imprinted, dec: &mut Decider) -> Answered<Classification>;
pub fn keep(classification: &Classification, imprinted: &Imprinted, op: BoolOp) -> Selection;
pub fn sew(model: &mut Model, imprinted: &Imprinted, selection: &Selection, dec: &mut Decider) -> Answered<BodyId>;
pub fn classify_point(model: &mut Model, body: BodyId, point: &RPoint3, dec: &mut Decider) -> Outcome<Containment>;
pub fn face_point(model: &mut Model, face: FaceId, dec: &mut Decider) -> Answered<RPoint3>;
```

## What defines correct behaviour

The test writer builds tests from the cases below. Every expected answer follows from this description and maths, so tests compare against it, never against a run of the code.

### How the cases are written

Operands are built with piece 3's builders:

- `block(s, t)`: the block of size `s` moved by the translation `t`;
- `cyl(r, h, t)`: the cylinder of radius `r` and height `h` (local axis `z` from 0 to `h`) moved by `t`; `xcyl(r, h, t)` is the same with the rotation of quaternion `(1, 0, 1, 0)`, whose matrix `[[0, 0, 1], [0, 1, 0], [-1, 0, 0]]` takes the local `z` axis to `+x`, then the translation `t`;
- `ball(R, t)`: the sphere of radius `R` moved by `t`.

A result is checked by:

1. `Model::check` gives no fault;
2. its counts: lumps, shells, faces, edges, ring edges, vertices, faces by surface kind, and edges by kind, where an edge's kind is the kind of piece 4's `curve_form` of its curve (`line`, `circle`, `ellipse`, `branch`), as the manifest counts them (the manifest's key for `branch` is `intersection`). The test computes edge kinds with `curve_form` of each edge of `edges_of_body`; piece 3's `EntityCounts::edges_by_kind` counts definitions by `CurveKind` and differs for intersection curves with a line or circle form;
3. `classify_point` of the listed rational points against the result;
4. every face of the result is on one of the operands' surfaces (`same_surface`), with the sense the material requires: at a rational point `q` strictly inside the face, the point moved by a small rational step along the outward normal (the surface's normal times the face's sense) is outside the result by the test's own rational oracle for the operation's material, and the point moved the other way is inside. The step is chosen smaller than the distance from `q` to every other face, which the test bounds with its oracle. `q` is `face_point` of the face when that distance is positive. The distance is zero when `face_point` lands where the face touches another face of the result with no edge or vertex there (ND4): an isolated touching point, such as `(5, 5, 10)` in cases 19 and 20, or a point of a tangent line between touching lumps, as in cases 12 and 26. No step works there, so the test uses its own `q` instead: a rational point it proves on the face's surface (implicit exactly 0), inside the face and off every other face, such as `(5, 5, 4)` on the void face of case 20 `A − B`.

`F / E (R) / V` below means faces / edges (ring edges among them) / vertices. Unless a case says otherwise the result has 1 lump and 1 shell.

The counts follow from the material and the rules above. Each was checked against the Euler characteristic of piece 3 (`V + R - E + sum of (2 - loops)` over faces, `= 2 - 2g` per shell), given in brackets where the shell is not a sphere.

### Two blocks

1. General: `A = block([2, 2, 2], 0)`, `B = block([2, 2, 2], (1, 1, 1))`.
   - `A ∪ B`: 12 / 30 / 20, all planes and lines. Inside: `(1/2, 1/2, 1/2)`, `(5/2, 5/2, 5/2)`, `(3/2, 3/2, 3/2)`. Outside: `(5/2, 1/2, 1/2)`. On: `(2, 1/2, 1/2)`.
   - `A ∩ B`: the cube `[1, 2]^3`: 6 / 12 / 8. Inside: `(3/2, 3/2, 3/2)`. Outside: `(1/2, 1/2, 1/2)`.
   - `A − B`: 9 / 21 / 14. Inside: `(1/2, 1/2, 1/2)`. Outside: `(3/2, 3/2, 3/2)`. On: `(3/2, 3/2, 1)`.
   - `B − A`: 9 / 21 / 14.
2. Face contact: `A = block([1, 1, 1], 0)`, `B = block([1, 1, 1], (1, 0, 0))`. `A ∪ B`: the box `[0, 2] x [0, 1]^2`, 6 / 12 / 8 (the shared face is gone and the coplanar faces merged). On: `(1, 1/2, 0)`. `A ∩ B`: `Empty`. `A − B`: 6 / 12 / 8, the faces on `A`'s surfaces.
3. Partial face contact: `A = block([2, 2, 2], 0)`, `B = block([2, 2, 2], (2, 1, 0))`. `A ∪ B`: the prism over the 8-gon `(0, 0), (2, 0), (2, 1), (4, 1), (4, 3), (2, 3), (2, 2), (0, 2)`: 10 / 24 / 16. Inside: `(2, 3/2, 1)`. `A ∩ B`: `Empty`. `A − B`: 6 / 12 / 8.
4. Edge contact: `A = block([1, 1, 1], 0)`, `B = block([1, 1, 1], (1, 1, 0))`. `A ∪ B`: 2 lumps, 2 shells, 12 / 24 / 16. On: `(1, 1, 1/2)`. `A ∩ B`: `Empty`. `A − B`: 6 / 12 / 8.
5. Corner contact: `B = block([1, 1, 1], (1, 1, 1))`, `A` as in case 4. `A ∪ B`: 2 lumps, 2 shells, 12 / 24 / 16. `A ∩ B`: `Empty`.
6. Through, the coaxial pose: `A = block([4, 4, 4], 0)`, `B = block([2, 2, 6], (1, 1, -1))`.
   - `A − B`: 10 / 24 / 16 [genus 1]. Outside: `(2, 2, 2)`. Inside: `(1/2, 2, 2)`.
   - `A ∪ B`: 16 / 36 / 24.
   - `A ∩ B`: the box `[1, 3]^2 x [0, 4]`: 6 / 12 / 8.
7. Equal size: `A = block([2, 2, 2], 0)`, `B = block([2, 2, 2], (1, 0, 0))`. `A ∪ B`, `A ∩ B` and `A − B` are the boxes `[0, 3] x [0, 2]^2`, `[1, 2] x [0, 2]^2` and `[0, 1] x [0, 2]^2`: each 6 / 12 / 8.
8. Identical: `B` a second block built with the same call as `A` (another body; passing `A` as both operands is a programming error). `A ∪ B` and `A ∩ B` give `A`'s counts, `A − B` is `Empty`.
9. Near miss: `A = block([1, 1, 1], 0)`, `B = block([1, 1, 1], (1 + 10^-30, 0, 0))`: `A ∪ B` has 2 lumps, 12 / 24 / 16; `A ∩ B` is `Empty`; `A − B` is 6 / 12 / 8. With `(1 - 10^-30, 0, 0)`: `A ∪ B`, `A ∩ B` and `A − B` are each 6 / 12 / 8, and `(1 - 10^-30/2, 1/2, 1/2)` is inside `A ∩ B`.

### A block and a cylinder

In cases 10 to 15, `A = block([10, 10, 10], 0)`.

10. Through hole: `B = cyl(3, 12, (5, 5, -1))`.
    - `A − B`: 7 / 14 (2) / 8, faces `{plane: 6, cylinder: 1}`, edges `{line: 12, circle: 2}` [genus 1]. Outside: `(5, 5, 5)`, `(5, 5, 10)`. Inside: `(1, 1, 1)`, `(5, 5 + 7/2, 5)`. On: `(5, 8, 5)`, `(5, 5 + 3, 0)`.
    - `A ∪ B`: 10 / 16 (4) / 8, `{plane: 8, cylinder: 2}`, `{line: 12, circle: 4}`. Inside: `(5, 5, 11 - 1/2)`, `(5, 5, 5)`.
    - `A ∩ B`: 3 / 2 (2) / 0, `{plane: 2, cylinder: 1}`, `{circle: 2}`. On: `(5, 5, 10)`.
11. Blind hole: `B = cyl(3, 5, (5, 5, 6))`. `A − B`: 8 / 14 (2) / 8, `{plane: 7, cylinder: 1}`. Outside: `(5, 5, 7)`. Inside: `(5, 5, 5)`.
12. Tangent line, from outside: `B = xcyl(2, 20, (-5, 5, 12))`, whose axis is the line `y = 5, z = 12`, at distance 2 above the top face. `A ∪ B`: 2 lumps, 2 shells, 9 / 14 (2) / 8. On: `(5, 5, 10)`. `A ∩ B`: `Empty`. `A − B`: 6 / 12 / 8.
13. Near miss across that tangent: with `(-5, 5, 12 + 10^-30)` the answers are those of case 12. With `(-5, 5, 12 - 10^-30)` the cylinder dips `10^-30` into the block, along a strip of the top face of width `2 sqrt(4 * 10^-30 - 10^-60)`:
    - `A ∪ B`: 10 / 20 (2) / 12, `{plane: 9, cylinder: 1}`, `{line: 16, circle: 4}`;
    - `A ∩ B`: 4 / 6 / 4, `{plane: 3, cylinder: 1}`, `{line: 4, circle: 2}`. Inside: `(5, 5, 10 - 10^-31)`;
    - `A − B`: 8 / 18 / 12, `{plane: 7, cylinder: 1}`, `{line: 16, circle: 2}`.
    The two section lines of the strip are intersection curves with a line form (piece 4, configuration PC3, with `h2 = 4 * 10^-30 - 10^-60`, not a rational square).
14. Equal radius, the cylinder inscribed in the block's section: `B = cyl(5, 12, (5, 5, -1))`, tangent to the four side faces along lines.
    - `A ∩ B`: 3 / 2 (2) / 0. The tangent lines leave no edge, and no vertex remains on the circles.
    - `A − B`: 4 lumps, 4 shells, 20 / 36 (0) / 24, `{plane: 16, cylinder: 4}`, `{line: 28, circle: 8}`. Each lump is one corner: 5 faces, 9 edges, 6 vertices, and the lumps that touch along a tangent line each have their own edge on it. Inside: `(1, 1, 5)`. On: `(0, 5, 5)`.
    - `A ∪ B`: 16 / 30 (2) / 16, `{plane: 14, cylinder: 2}`, `{line: 20, circle: 10}`. The top face of the block becomes four corner faces that touch at the four points `(5, 0, 10)`, `(10, 5, 10)`, `(5, 10, 10)`, `(0, 5, 10)`, which are vertices. This is not self-touching: at each point the faces around the vertex form one cycle (at `(5, 0, 10)`: the side face `y = 0`, a corner face, the part of the cylinder's side above the block, and the other corner face), so the boundary is a manifold there.
15. Along the seam: `B = cyl(3, 12, (0, 5, -1))`, whose axis lies in the face `x = 0`. The groove `B` cuts in `A` crosses `B`'s seam, the line `(3, 5, z)`.
    - `A − B`: 8 / 18 (0) / 12, `{plane: 7, cylinder: 1}`, `{line: 16, circle: 2}`. Outside: `(1, 5, 5)`. Inside: `(4, 5, 5)`.
    - `A ∩ B`: the half cylinder: 4 / 6 / 4, `{plane: 3, cylinder: 1}`, `{line: 4, circle: 2}`.
16. Inside, no contact: `B = cyl(2, 4, (5, 5, 3))`. `A − B`: 1 lump, 2 shells, 9 / 14 (2) / 8. `A ∪ B`: `A`'s counts, 6 / 12 / 8. `A ∩ B`: 3 / 2 (2) / 0.

### A block and a sphere

In cases 17 to 21, `A = block([10, 10, 10], 0)`.

17. At a corner: `B = ball(3, (10, 10, 10))`. `A − B`: 7 / 15 / 10, `{plane: 6, sphere: 1}`, `{line: 12, circle: 3}`. `A ∩ B`: 4 / 6 / 4, `{plane: 3, sphere: 1}`, `{line: 3, circle: 3}`. `A ∪ B`: 7 / 15 / 10.
18. Through a face: `B = ball(3, (5, 5, 10))`. `A − B`: 7 / 13 (1) / 8. `A ∩ B`: the half ball, 2 / 1 (1) / 0. `A ∪ B`: 7 / 13 (1) / 8.
19. Tangent from outside: `B = ball(3, (5, 5, 13))`. `A ∪ B`: 2 lumps, 2 shells, 7 / 12 / 8, with no vertex at the touching point `(5, 5, 10)`. `A ∩ B`: `Empty`. `A − B`: 6 / 12 / 8. Near miss: with `(5, 5, 13 + 10^-30)` the same; with `(5, 5, 13 - 10^-30)`, `A ∪ B` and `A − B` are 7 / 13 (1) / 8 and `A ∩ B` is 2 / 1 (1) / 0.
20. Tangent from inside: `B = ball(3, (5, 5, 7))`, touching the top face at `(5, 5, 10)`. `A − B`: 1 lump, 2 shells, 7 / 12 / 8 (the void touches the outer shell at a point, which adds no entity). `A ∪ B`: 6 / 12 / 8. `A ∩ B`: 1 / 0 / 0.
21. Equal radius, the inscribed ball `B = ball(5, (5, 5, 5))`: `A − B`: 1 lump, 2 shells, 7 / 12 / 8. `A ∪ B`: 6 / 12 / 8. `A ∩ B`: 1 / 0 / 0.
22. Along the seam and through the poles: `A = block([10, 10, 10], 0)`, `B = ball(3, (10, 5, 5))`, centred on the face `x = 10`; the section is the great circle through both poles of `B`. `A − B`: 7 / 13 (1) / 8. `A ∩ B`: 2 / 1 (1) / 0.

### Two cylinders

23. Crossing, the smaller through the larger: `A = cyl(2, 10, (0, 0, -5))`, `B = xcyl(1, 10, (-5, 0, 0))` (piece 4, configuration CC5).
    - `A ∪ B`: 7 / 6 (6) / 0, `{plane: 4, cylinder: 3}`, `{circle: 4, branch: 2}`.
    - `A ∩ B`: 3 / 2 (2) / 0, `{cylinder: 3}`, `{branch: 2}`.
    - `A − B`: 4 / 4 (4) / 0, `{plane: 2, cylinder: 2}` [genus 1]. Outside: `(0, 0, 0)`. Inside: `(0, 3/2, 0)`.
    - `B − A`: 2 lumps, 2 shells, 6 / 4 (4) / 0.
24. Steinmetz, equal radii and crossing axes: `A = cyl(1, 4, (0, 0, -2))`, `B = xcyl(1, 4, (-2, 0, 0))` (CC6).
    - `A ∩ B`: 4 / 4 (0) / 2, `{cylinder: 4}`, `{ellipse: 4}`. The vertices are `(0, 1, 0)` and `(0, -1, 0)`.
    - `A ∪ B`: 8 / 8 (4) / 2, `{plane: 4, cylinder: 4}`, `{circle: 4, ellipse: 4}`.
    - `A − B`: 2 lumps, 2 shells, 6 / 6 (2) / 4. The two lumps (`z > 0` and `z < 0`) touch at `(0, ±1, 0)`, and each has its own vertices there.
25. Parallel, overlapping: `A = cyl(5, 10, 0)`, `B = cyl(5, 10, (8, 0, 0))` (CP3, lines through `(4, 3, z)` and `(4, -3, z)`). `A ∪ B`: 4 / 6 (0) / 4, `{plane: 2, cylinder: 2}`, `{line: 2, circle: 4}`. `A ∩ B`: 4 / 6 / 4. `A − B`: 4 / 6 / 4.
26. Parallel, tangent: `B = cyl(5, 10, (10, 0, 0))`. `A ∪ B`: 2 lumps, 2 shells, 6 / 4 (4) / 0. `A ∩ B`: `Empty`.
27. Coaxial: `A = cyl(2, 5, 0)`, `B = cyl(2, 5, (0, 0, 5))`: `A ∪ B` is one cylinder, 3 / 2 (2) / 0; `A ∩ B` is `Empty`. `A = cyl(3, 10, 0)`, `B = cyl(1, 12, (0, 0, -1))`: `A − B` is the tube, 4 / 4 (4) / 0 [genus 1]; `A ∪ B` is 7 / 6 (6) / 0; `A ∩ B` is 3 / 2 (2) / 0.
28. Inside tangent, two branches through a node: `A = cyl(2, 10, (0, 0, -5))`, `B = xcyl(1, 10, (-5, 1, 0))` (CC4, node `(0, 2, 0)`).
    - `A ∪ B`: 7 / 6 (4) / 1, `{plane: 4, cylinder: 3}`, `{circle: 4, branch: 2}`.
    - `A ∩ B`: 3 / 2 (0) / 1, `{cylinder: 3}`, `{branch: 2}`.
    - `B − A`: 2 lumps, 2 shells, 6 / 4 (2) / 2; the lumps touch at the node and each has its own vertex there.
    - `A − B`: refused with `boolean.self-touching-shell`, stage `sew`: the material touches itself at the node, where edges meet.
29. Tangent at a point: `A = cyl(2, 10, (0, 0, -5))`, `B = xcyl(1, 10, (-5, 3, 0))` (CC2). `A ∪ B`: 2 lumps, 2 shells, 6 / 4 (4) / 0, with no vertex at `(0, 2, 0)`. With `(-5, 3 - 10^-30, 0)` (CC3, a tiny bite): `A ∪ B` is 6 / 5 (5) / 0.

### A cylinder and a sphere

30. Coaxial, through: `A = cyl(3, 20, (0, 0, -10))`, `B = ball(5, 0)` (CS1, circles at `z = ±4`).
    - `A ∩ B`: 3 / 2 (2) / 0, `{cylinder: 1, sphere: 2}`.
    - `A ∪ B`: 5 / 4 (4) / 0, `{plane: 2, cylinder: 2, sphere: 1}`.
    - `B − A`: 2 / 2 (2) / 0 [genus 1].
    - `A − B`: 2 lumps, 2 shells, 6 / 4 (4) / 0.
    The same with `B = ball(6, 0)` (`k2 = 27`, intersection curves with a circle form) gives the same counts, the section edges counting as circles.
31. Equal radius, the ball inscribed in the cylinder: `A = cyl(5, 20, (0, 0, -10))`, `B = ball(5, 0)`, touching along the equator (CS2). `A ∪ B`: `A`'s counts, 3 / 2 (2) / 0. `A ∩ B`: 1 / 0 / 0. `B − A`: `Empty`. `A − B`: 2 lumps, 2 shells, 6 / 4 (4) / 0; the upper and lower parts touch along the equator, and each has its own ring edge on it.
32. Viviani, a node: `A = cyl(1, 6, (0, 0, -3))`, `B = ball(2, (1, 0, 0))` (CS7, node `(-1, 0, 0)`).
    - `A ∩ B`: 3 / 2 (0) / 1, `{cylinder: 1, sphere: 2}`, `{branch: 2}`.
    - `A ∪ B`: 5 / 4 (2) / 1, `{plane: 2, cylinder: 2, sphere: 1}`, `{circle: 2, branch: 2}`.
    - `A − B`: 2 lumps, 2 shells, 6 / 4 (2) / 2.
    - `B − A`: refused with `boolean.self-touching-shell`, stage `sew`.
33. One branch, a bite: `A = cyl(1, 6, (0, 0, -3))`, `B = ball(2, (2, 0, 0))` (CS6). `A ∪ B`: 4 / 3 (3) / 0. `A ∩ B`: 2 / 1 (1) / 0. `A − B`: 4 / 3 (3) / 0. `B − A`: 2 / 1 (1) / 0.
34. Tangent at a point: `A = cyl(1, 6, (0, 0, -3))`, `B = ball(2, (3, 0, 0))` (CS5). `A ∪ B`: 2 lumps, 2 shells, 4 / 2 (2) / 0.
35. Along the seam: `A = cyl(1, 6, (0, 0, -3))`, `B = ball(1, (1, 0, 0))`, centred on `A`'s seam line `(1, 0, z)` (CS6). The counts are those of case 33.

### Two spheres

36. Crossing: `A = ball(5, 0)`, `B = ball(5, (6, 0, 0))`. `A ∪ B`, `A ∩ B` and `A − B`: each 2 / 1 (1) / 0. Inside `A ∩ B`: `(3, 0, 0)`. On `A ∪ B`: `(3, 0, 4)`.
37. Tangent from outside: `B = ball(5, (10, 0, 0))`. `A ∪ B`: 2 lumps, 2 shells, 2 / 0 / 0, with no vertex at `(5, 0, 0)`. `A ∩ B`: `Empty`. With `(10 + 10^-30, 0, 0)`: the same. With `(10 - 10^-30, 0, 0)`: `A ∪ B` is 2 / 1 (1) / 0.
38. Tangent from inside: `B = ball(2, (3, 0, 0))`. `A − B`: 1 lump, 2 shells, 2 / 0 / 0. `A ∪ B`: 1 / 0 / 0. `A ∩ B`: 1 / 0 / 0. `B − A`: `Empty`.
39. Concentric, the coaxial pose: `B = ball(3, 0)`. `A − B`: 1 lump, 2 shells, 2 / 0 / 0. `A ∪ B` and `A ∩ B`: 1 / 0 / 0.
40. Equal radius and the same centre, `B` another body with `A`'s fields: `A ∪ B` and `A ∩ B` are 1 / 0 / 0; `A − B` is `Empty`.

### Chains

41. Two crossing holes: `A = block([10, 10, 10], 0)`, `B = cyl(2, 12, (5, 5, -1))`, `C = xcyl(2, 12, (-1, 5, 5))`. `(A − B) − C`: 10 / 20 (4) / 10, `{plane: 6, cylinder: 4}`, `{line: 12, circle: 4, ellipse: 4}` [genus 3]. Outside: `(5, 5, 5)`. Inside: `(1, 1, 1)`.
42. A frame from four blocks: `B1 = block([3, 1, 1], 0)`, `B2 = block([1, 3, 1], (2, 0, 0))`, `B3 = block([3, 1, 1], (0, 2, 0))`, `B4 = block([1, 3, 1], 0)`. `B1 ∪ B2`: 8 / 18 / 12. `(B1 ∪ B2) ∪ B3`: 10 / 24 / 16. `((B1 ∪ B2) ∪ B3) ∪ B4`: 10 / 24 / 16 [genus 1]. Outside: `(3/2, 3/2, 1/2)`.
43. Two lens cuts: `A = ball(5, 0)`, `B = ball(5, (6, 0, 0))`, `C = ball(5, (-6, 0, 0))`. `(A − B) − C`: 3 / 2 (2) / 0.
44. A chain of five unions of blocks `block([2, 1, 1], (k, 0, 0))` for `k = 0` to 5, folded left: the result is the box `[0, 7] x [0, 1]^2`, 6 / 12 / 8, after every step.

### Edge cases

45. Empty operands: with `E` from `empty(model)` and `A` any operand above, `A ∪ E` and `E ∪ A` have `A`'s counts, `A ∩ E`, `E ∩ A` and `E − A` are `Empty`, `A − E` has `A`'s counts, and `E ∪ E` is `Empty`.
46. A body of another kind would be refused with `boolean.not-solid`, stage `intersect`; piece 3 has no such kind yet, so no test builds one. Valid cone and torus cases are checked against exact mass references, closed-shell incidence, Euler characteristic and STEP export, rather than an expectation that their kind is unsupported.
47. A connected piece touching itself along a curve: `A = block([10, 10, 10], 0)`, `B = cyl(5, 6, (5, 5, 2))`. `A − B` is the block minus an inscribed cylinder that stops short of the top and bottom faces; the corner regions are joined through the slabs `z < 2` and `z > 8` and touch each other along four segments. Refused with `boolean.self-touching-shell`, stage `sew`.

### Stages

48. The graph of case 10 (`A`, `B`): two section edges, each a ring on a circle of radius 3 about the line `x = y = 5` at `z = 0` and `z = 10`, each with one face of `A` (the bottom or top plane) and one of `B` (its side); no isolated point.
49. The graph of case 37 at `(10, 0, 0)`: no section edge; one isolated point, at `(5, 0, 0)`.
50. The graph of case 2: the boundary of the common square `x = 1`, as four section edges along `A`'s and `B`'s edges, each lying on two faces of each body.
51. Classification, for every case 1 to 40: for every piece `P` of either body, `face_point(P)` lies on `P`'s surface (its implicit is exactly 0) and is classified by the test's own maths against the other operand; `class_of(P)` must be `Inside` or `Outside` accordingly. A probe on the other operand's boundary gives `OnSame` or `OnOpposite` by the two outward normals unless it is an isolated point of the graph, in which case the test chooses another probe. This checks imprint's classification rule too.
52. Keep: for every case, `keep` returns exactly the pieces the table names, with the senses it names, in the stated order.
53. Reports: for operands that are not `Empty`, `report.stages` has five records in stage order. With `before` and `after` the `Decider`'s statistics just before and after the call, each counter other than `max_bits` in `after` equals that counter in `before` plus its sum over the five records. The records' `max_bits` never decrease from one stage to the next, the first is at least `before.max_bits`, and the last equals `after.max_bits`. The test runs this on a fresh `Decider` and on one that has already made decisions (for example a `classify_point` call before the boolean), so a build that clears the caller's totals fails.

### Points

54. `classify_point` on primitives, with exact answers from their implicit formulas: for `block([2, 3, 5], 0)`, `(1, 1, 1)` inside, `(2, 1, 1)` on, `(3, 1, 1)` outside, `(2, 3, 5)` on (a corner); for `cyl(2, 5, 0)`, `(0, 0, 1)` inside, `(2, 0, 1)` on, `(6/5, 8/5, 5)` on (the top rim), `(2, 2, 1)` outside; for `ball(5, 0)`, `(3, 4, 0)` on, `(3, 4, 10^-30)` outside, `(3, 4 - 10^-30, 0)` inside.
55. A generated family on random rational points, near and far from the boundaries, against the primitives under random poses and against the results of cases 1 to 40, with the answer computed by the test from the operands' implicit formulas and the operation's set rule. Points that are on an operand boundary are skipped unless the case lists them.

### Metamorphic properties

56. For every case 1 to 40, and for random rational points not on either operand's boundary: `classify_point(A ∪ B, p)` is `Inside` exactly when the point is inside `A` or inside `B`; `A ∩ B` when inside both; `A − B` when inside `A` and outside `B`. `A ∪ B` and `B ∪ A` have equal counts, and so do `A ∩ B` and `B ∩ A`.
57. Scale and placement (DESIGN.md, R9): a variant of a case multiplies every size and translation of its operands by `10^k`, then moves them by a pose. The 20 poses are random rational poses (integer quaternion components up to 100, translations up to `10^6`) drawn once from a fixed seed. Each variant gives the case's counts, and the case's listed points, scaled and moved the same way, classify the same way. The default run gives every operation of cases 1 to 44 three variants:
    - `k = -30` at the identity pose;
    - `k = 30` at the identity pose;
    - `k` drawn from `-29..=29` by a fixed seed, with pose number `i mod 20` for the `i`-th operation in case order (counting from 0), so every pose is used.

    The steps of a chain share the variant of its first operation, and `i` counts every step.

    A long-run qualification test, ignored in the default run and run on demand, gives every operation the variants with every `k` in `-30..=30` at the identity pose and with every pose at `k` = -30, 0 and 30: 121 variants per operation. It prints its wall time, and no time limit is part of the contract.
58. Determinism: each case run twice on two fresh models gives equal counts, equal marks and equal statistics; run on eight threads, one model each, identical results.
59. Budgets: under `BudgetPolicy { max_precision_bits: 64, max_work: 10_000 }` (this budget admits piece 1's `f64` filter, which costs one work unit), each case gives its default answer, a marked answer, or a refusal; never a wrong `Proven` containment, and never a result that fails `check`. A `Decider` whose deadline has passed: each case either refuses with `numbers.deadline` or gives its default-budget answer (ARCHITECTURE.md, "The time budget"); the stage of the refusal is not part of the contract.

### Volume and running time of the tests

The generated families (cases 51, 55, 56 and 57) must together run at least 30,000 maths-answered questions per test run, from fixed seeds; case 57 counts only its default-run variants, and its long-run test is not part of the default run. Booleans are expensive; spend most of the count on `classify_point` against a few hundred results, and keep the default test run under about 15 minutes in release mode on one core.

## Values checked for this description

The architect derived every count above from the material of the operation and the touching and maximality rules, and checked each against piece 3's Euler characteristic for its genus. The points listed as inside, outside or on were checked against the operands' implicit formulas in exact rationals. Volumes for the cross-checks of piece 6 (for example `1000 - 90 pi` for case 10 `A − B`, `16/3` for case 24 `A ∩ B`, and `1000 - 80 pi + 128/3` for case 41) are listed in piece 6's description.

## Speed: measured, never gated

No test fails for speed. A benchmark target measures and prints:

1. time per stage and per operation, by shape pair and pose (cases 1 to 40);
2. decisions per stage and the share settled by facts and by piece 1's filter;
3. growth of number size along chains: the largest precision used and the size of the largest observable stored rational after each step of cases 41 to 44 and of a chain of five alternating unions and subtractions of blocks, cylinders and spheres under random rational poses. Observable rationals include public rational fields, implicit coefficients, root boxes and root-system coefficients reachable from the result. Opaque `Real` coordinates are counted separately; no interface is added to inspect them;
4. time of `classify_point` against face count.

The numbers are written to the benchmark output so later runs can be compared.

## Out of scope for this piece

- Fillet, shell and offset, which reuse these stages later.
- Sheets, wires, meshes and lattices.
- Any tolerance. This piece has none.
