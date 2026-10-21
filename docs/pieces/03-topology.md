# Piece 3: faces and edges (`catallax-topology`)

Status: complete for blind test writing and blind building.

This is the input the test writer and the creator get, together with `docs/DESIGN.md`, `docs/ARCHITECTURE.md`, piece 1's description `docs/pieces/01-numbers-and-decisions.md` and piece 2's description `docs/pieces/02-shapes.md`. This piece stores its geometry in piece 2's `GeomStore` and uses piece 1's facts, refusals and marks exactly as those descriptions define them. Nothing else is needed.

## Approved algebraic-definition migration

Under the 2026-10-22 architecture decision, geometry handles and append-only topology stay unchanged. Piece 2's poses and stored definitions become algebraic; the rational builder size/radius inputs below remain valid convenience APIs for the existing exchange. Builders convert exact rational inputs to the shared definitions without rounding. A transformed primitive may contain algebraic fields even when its local dimensions are rational. The rational-only pose statements below describe the retained baseline.

Recipe builders retain each operation's exact profile, sweep law, finite trim and origin set when reducing a side to an analytic carrier. Reusing a carrier does not merge operation identities. Oblique cylinder caps remain the translated sketch planes. Open conics require bounded edges with endpoints; an edge without endpoints requires certified closure from piece 4. Topology checks may use recorded closure facts rather than depending on the intersection crate. They must not infer closure from a curve being non-linear.

The rational subset of sketch/extrude/revolve may land before the shared foundation. Temporary failures for irrational dimensions, frames or unsupported carriers must name the missing coverage and remain repair cases. Only the recipe's proven impossible inputs justify an impossible-geometry refusal. Coverage of a legal feature requires its topology, measurement and export paths, not just successful side-surface construction.

## What this piece is for

Piece 3 holds the boundary representation: bodies, lumps, shells, faces, loops, coedges, edges and vertices, all over piece 2's exact geometry.
It builds the primitive solids that questions name (block, cylinder, sphere, cone, torus) and checks a body's topology.

The vocabulary carries curved geometry from the first commit (DESIGN.md, R6). A face lies on any surface, an edge on any curve, and a vertex at any point definition, including the meeting point of three curved surfaces.
Nothing in this piece assumes planar faces, straight edges, vertices where three planes meet, or a loop that is the "outer" boundary.

Coordinates are a cache (DESIGN.md, R4). A vertex holds a `PointId`, an edge a `CurveId`, a face a `SurfaceId`. No type of this piece has a coordinate field, and this piece makes no geometric decision: its builders work in exact rationals and its check reads recorded facts.

## Number types

| What | Type | Meaning |
| --- | --- | --- |
| Builder sizes and radii | `Rational` | Exact, given by the caller |
| Poses | piece 2 `Pose` | Exact rational rotation and translation |
| Geometry of every entity | `SurfaceId`, `CurveId`, `PointId` into the model's `GeomStore` | Exact definitions (piece 2) |
| Marks carried by a body | `MarkLog` (piece 1) | The record of every not-exact decision behind the body |

This piece holds no `Real`, `Ball`, `Interval` or `f64` value of its own.

## The model

`Model` owns one `GeomStore` and one write-once arena per entity kind. Entities never change once added; operations create new entities and return new handles.

| Entity | Fields | Rules |
| --- | --- | --- |
| `Vertex` | `point: PointId` | Any point definition. |
| `Edge` | `curve: CurveId`, `sense: Sense`, `start`, `end: Option<VertexId>` | `sense` relates the edge direction to the curve's natural direction. `start` and `end` are given along the edge direction. Both `None` is a ring edge covering a whole closed curve; otherwise both are `Some`, and they are equal only on a closed curve with one vertex. |
| `Coedge` | `edge: EdgeId`, `sense: Sense` | One use of an edge by one loop. `sense` relates the loop's direction of travel to the edge direction. An edge may have any number of coedges. |
| `Loop` | `Coedges(Vec<CoedgeId>)` or `Vertex(VertexId)` | A cycle of coedges, each ending where the next starts; or a single point boundary, such as a cone's apex. |
| `Face` | `surface: SurfaceId`, `sense: Sense`, `loops: Vec<LoopId>` | The face normal is the surface normal (`Same`) or its reverse (`Reversed`). No loop is outer. Zero loops means the whole closed surface (a sphere or a torus). |
| `Shell` | `faces: Vec<FaceId>` | A connected set of faces. |
| `Lump` | `shells: Vec<ShellId>` | The first shell is the outer boundary; the rest bound voids. |
| `Body` | `kind: BodyKind`, `lumps: Vec<LumpId>`, `marks: MarkLog` | `marks` are the marks of the operation that produced the body, merged with those of its operands. |

`BodyKind` is `Solid` or `Empty`. It is non-exhaustive: sheets, wires, meshes and lattices come later (DESIGN.md, scope).

A coedge's start vertex is its edge's `start` when its sense is `Same`, and its edge's `end` when `Reversed`; its end vertex is the other one. A ring edge's coedges have neither.

### Orientation

- A face's normal points out of the material.
- Travelling along a coedge loop with the face normal pointing towards the viewer, the face lies on the left. In vectors: at a point of the loop where the travel direction is `d` and the face normal is `n`, the face lies in the direction `n x d`.
- A void shell's faces point into the void, which is out of the material.

### Sharing

- Different bodies may share entities. A result may list an operand's unchanged faces, loops and edges, since nothing ever changes once added.
- Within one body, no face, edge or vertex is used by two different shells. Where material touches itself, each side keeps its own entities (`JUDGE-INTERFACE.md`, the answer's topology). The geometry definitions under them (curves, points, surfaces) may be shared.
- Within one shell, every edge is used by exactly two coedges, and their senses are opposite.

### Adding entities and back-links

`add_vertex`, `add_edge`, `add_coedge`, `add_loop`, `add_face`, `add_shell`, `add_lump` and `add_body` append and return the new handle, from index 0 in each arena.
They check nothing beyond handle range: a model is built bottom-up, so a half-built body cannot be judged until it is complete. `check` judges a complete body.
A handle out of range, or from another model, is a programming error and may panic. So are the getters (`vertex`, `edge`, `coedge`, `loop_`, `face`, `shell`, `lump`, `body`) on such a handle.

The model keeps two back-links as entities are added:

- `coedges_of_edge(e)`: every coedge of the model whose edge is `e`, in creation order, whether or not it is in a face yet.
- `owner_of_coedge(c)`: the loop and face that own `c`, or `None` before a face lists a loop that holds it. When `add_face` adds a face, each coedge in each of its loops that has no owner yet gets that loop and that face, in the face's loop order. The first face added wins; a coedge listed twice is a fault that `check` reports.

## Queries

Bodies are walked in one fixed order: lumps in order, then each lump's shells in order, each shell's faces in order, each face's loops in order, and each loop's coedges in order.

- `faces_of_body(b)`: faces in that order, each listed once, at first sight.
- `edges_of_body(b)`: the edge of each coedge in that order, each listed once, at first sight.
- `vertices_of_body(b)`: in that order, for each coedge its start vertex then its end vertex (none for a ring edge's coedge), and for each vertex loop its vertex; each listed once, at first sight.
- `counts(b)`: `lumps` and `shells` are the lengths of the lump and shell lists walked; `faces`, `edges` and `vertices` are the lengths of the three lists above; `loops` is the number of loops of those faces; `ring_edges` is the number of those edges with neither end; `faces_by_kind` and `edges_by_kind` count those faces by their surface's `SurfaceKind` and those edges by their curve's `CurveKind`. A kind with a count of zero has no entry in the map. An empty body has every count zero and empty maps.

## Topology check

`check(b)` returns every topological fault of a body, or an empty list for a sound body.
It reads the entities and the fact book only, and makes no numeric decision. Geometric facts that topology needs (an edge's curve lies on its face's surface) are checked as recorded facts. A geometric fault with no topological sign (a vertex moved off its surfaces with its facts still recorded) is out of its reach.

| Code | Fault |
| --- | --- |
| `topology.body-kind` | A `Solid` body with no lump, or an `Empty` body with a lump. |
| `topology.empty-part` | A lump with no shell, a shell with no face, or a coedge loop with no coedge. |
| `topology.shared-entity` | Within the body: a lump, shell, face, loop or coedge listed twice anywhere in the walk (including a face in two shells, a loop in two faces, a coedge in two loops), or an edge or vertex used in two different shells. |
| `topology.loop-not-closed` | In a coedge loop, a coedge's end vertex is not the next coedge's start vertex (the last coedge's end against the first's start); a ring edge's coedge that is not alone in its loop. Vertices are compared as handles. |
| `topology.edge-ends` | An edge with exactly one of `start` and `end`; or a ring edge, or an edge whose `start` equals its `end`, on a curve with `is_closed() == Some(false)`. |
| `topology.edge-use` | Within one shell, an edge used by a number of coedges other than two, or by two coedges of the same sense. |
| `topology.unbounded-face` | A face with no loops whose surface is a plane, a cylinder or a cone. |
| `topology.vertex-loop` | A face with a vertex loop and no coedge loop, or a face where a vertex loop's vertex also appears in another loop of the same face. |
| `topology.disconnected-shell` | A shell whose faces are not all connected through edges they share; two faces are neighbours when coedges of both use one edge. A shell with one face is connected. |
| `topology.euler` | A shell whose Euler characteristic, counted as below, is not of the form `2 - 2g` for a whole number `g >= 0`. |
| `topology.missing-incidence` | A fact the body needs is not recorded: an edge's curve `Incident` on the surface of a face that uses the edge; an edge's start or end point `Incident` on the edge's curve; a vertex loop's point `Incident` on its face's surface. Each is asked of the fact book with `holds(Incident, ...)`. |

The Euler characteristic of a shell is `V + R - E + sum over faces of chi(face)`, where `V` counts the shell's distinct vertices (vertex-loop vertices included), `R` its distinct ring edges (each ring edge counts as one edge with one virtual vertex), `E` its distinct edges (ring edges included), and `chi(face)` is `2 - k` for a face with `k` loops (vertex loops included), except that a loopless face on a torus counts 0.
This reads each face as a sphere with `k` holes and the closed torus as a torus. The count is a necessary condition only.

Values for the builders: block `8 + 0 - 12 + 6 = 2`; cylinder and frustum `0 + 2 - 2 + (1 + 1 + 0) = 2`; cone to an apex `1 + 1 - 1 + (1 + 0) = 2`; sphere `2`; torus `0` (genus 1).

The order of the returned faults is deterministic but not part of the contract. One planted fault may cause several faults. Tests check that the expected code is among them and that a sound body gives an empty list.

## Primitive solids

Each builder validates its sizes, adds exact geometry to the store, records the construction facts below, adds the topology, and returns a new `Solid` body with one lump, one shell and an empty mark log.
It makes no decision: every check is exact rational arithmetic.

A refused builder returns the refusal and adds nothing: every arena of the model and of its store, and the fact book, keep their lengths. No builder ever panics on a size.

### Validation

Checks run in this order, and the first that fails gives the refusal:

| Builder | `topology.non-positive-size` when | Then |
| --- | --- | --- |
| `block(size, pose)` | any of `sx`, `sy`, `sz` is `<= 0` | |
| `cylinder(radius, height, pose)` | `radius <= 0` or `height <= 0` | |
| `sphere(radius, pose)` | `radius <= 0` | |
| `cone(base_radius, top_radius, height, pose)` | `base_radius <= 0`, `height <= 0` or `top_radius < 0` | `topology.cone-is-cylinder` when `top_radius == base_radius` |
| `torus(major, minor, pose)` | `major <= 0` or `minor <= 0` | `topology.torus-self-intersecting` when `minor >= major` |

Because a builder validates first, it never meets a refusal from piece 2's constructors.

### Geometry

Each builder makes its geometry in the local frame below and then moves it: every surface is `local.transformed(pose)`, every curve `local.transformed(pose)`, and every vertex `PointDef::Explicit(pose.apply_point(local))`. Tests may compare each field with `==`.
Every builder uses `ref_dir = (1, 0, 0)` for its circles and its curved surfaces, so seams lie on the local half-plane `y = 0, x > 0`, where the judge contract (`JUDGE-INTERFACE.md`, seams and poles) promises them.

**Block** `[0, sx] x [0, sy] x [0, sz]`. Name the corner `(i sx, j sy, k sz)` by `ijk`, for `i, j, k` in `{0, 1}`.

| Face | Plane origin | Plane normal | Plane `ref_dir` | Loop, as a cycle of corners |
| --- | --- | --- | --- | --- |
| `x = 0` | `(0, 0, 0)` | `(-1, 0, 0)` | `(0, 1, 0)` | `000, 001, 011, 010` |
| `x = sx` | `(sx, 0, 0)` | `(1, 0, 0)` | `(0, 1, 0)` | `100, 110, 111, 101` |
| `y = 0` | `(0, 0, 0)` | `(0, -1, 0)` | `(0, 0, 1)` | `000, 100, 101, 001` |
| `y = sy` | `(0, sy, 0)` | `(0, 1, 0)` | `(0, 0, 1)` | `010, 011, 111, 110` |
| `z = 0` | `(0, 0, 0)` | `(0, 0, -1)` | `(1, 0, 0)` | `000, 010, 110, 100` |
| `z = sz` | `(0, 0, sz)` | `(0, 0, 1)` | `(1, 0, 0)` | `001, 101, 111, 011` |

- Faces are in the shell in this order, each with sense `Same`.
- 8 vertices, one per corner. 12 edges, one per side of the box. Each edge lies on the line whose `origin` is its corner nearer the local origin and whose `dir` is `(1, 0, 0)`, `(0, 1, 0)` or `(0, 0, 1)` along that side; its sense is `Same`, its `start` that nearer corner and its `end` the other.
- Each loop visits its corners in the order listed, up to its starting corner. A coedge's sense is `Same` when the loop travels its edge from `start` to `end`, else `Reversed`.

**Cylinder**, axis along `z` from 0 to `height`.

| Face | Surface | Sense | Loops |
| --- | --- | --- | --- |
| bottom | plane: origin `(0, 0, 0)`, normal `(0, 0, -1)`, `ref_dir` `(1, 0, 0)` | `Same` | bottom ring, coedge `Reversed` |
| top | plane: origin `(0, 0, height)`, normal `(0, 0, 1)`, `ref_dir` `(1, 0, 0)` | `Same` | top ring, coedge `Same` |
| side | cylinder: origin `(0, 0, 0)`, axis `(0, 0, 1)`, `radius`, `ref_dir` `(1, 0, 0)` | `Same` | bottom ring, coedge `Same`; top ring, coedge `Reversed` |

The bottom ring is a ring edge (sense `Same`) on the circle with centre `(0, 0, 0)`, normal `(0, 0, 1)`, `radius_sq = radius^2`, `ref_dir` `(1, 0, 0)`; the top ring the same with centre `(0, 0, height)`. No vertices.

**Cone**, axis along `z`, `base_radius` at `z = 0` and `top_radius` at `z = height`.
The cone surface has `ref_dir` `(1, 0, 0)` and

- when `top_radius < base_radius`: apex `(0, 0, height base_radius / (base_radius - top_radius))`, axis `(0, 0, -1)`, slope `(base_radius - top_radius) / height`;
- when `top_radius > base_radius`: apex `(0, 0, -height base_radius / (top_radius - base_radius))`, axis `(0, 0, 1)`, slope `(top_radius - base_radius) / height`.

Faces, in order: the bottom plane as for the cylinder; the top plane as for the cylinder when `top_radius > 0`; the side, on the cone surface. All have sense `Same`, since the cone's normal points away from its axis.
The rings are as for the cylinder, with `radius_sq` `base_radius^2` at the bottom and `top_radius^2` at the top, and the coedge senses are the cylinder's.
When `top_radius = 0` there is no top plane and no top ring. The apex `(0, 0, height)` is one vertex, and the side has two loops: the bottom ring (coedge `Same`) and a vertex loop on the apex.

**Sphere**, centred at the origin: one face with no loops, sense `Same`, on the sphere with centre `(0, 0, 0)`, `radius`, axis `(0, 0, 1)` and `ref_dir` `(1, 0, 0)`.

**Torus**, centred at the origin about the `z` axis: one face with no loops, sense `Same`, on the torus with centre `(0, 0, 0)`, axis `(0, 0, 1)`, `major`, `minor` and `ref_dir` `(1, 0, 0)`.

`empty(model)` adds a body of kind `Empty` with no lumps and an empty mark log, and adds nothing else.

### What each builder adds

| Builder | Surfaces | Curves | Points | Vertices | Edges | Coedges | Loops | Faces |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| block | 6 | 12 | 8 | 8 | 12 | 24 | 6 | 6 |
| cylinder | 3 | 2 | 0 | 0 | 2 | 4 | 4 | 3 |
| cone, `top_radius > 0` | 3 | 2 | 0 | 0 | 2 | 4 | 4 | 3 |
| cone, `top_radius = 0` | 2 | 1 | 1 | 1 | 1 | 2 | 3 | 2 |
| sphere | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| torus | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

Each also adds one shell, one lump and one body. No builder reuses geometry or topology that was in the model before the call.

### Construction facts

Each builder records at least these facts, with origin `FactOrigin::Construction("topology.block")` (`"topology.cylinder"`, `"topology.sphere"`, `"topology.cone"`, `"topology.torus"`):

- `Incident(curve, surface)` for each edge's curve and the surface of each face that uses the edge;
- `Incident(point, curve)` for each edge's start and end point and the edge's curve;
- `Incident(point, surface)` for each vertex loop's point and its face's surface;
- `Parallel(a, b)` for the planes of each pair of opposite block faces, and for the bottom and top planes of a cylinder or of a cone with `top_radius > 0`.

A builder may record more facts, but only facts its construction makes true. So these must not hold after a block: `Parallel` between the planes `x = 0` and `y = 0`; `Incident` of the corner `000` on the plane `x = sx`.

## How this piece uses piece 1

- **Outcomes.** This piece returns no `Outcome`: builders return `Result<BodyId, Refusal>` and `check` returns faults. It decides nothing that could end not exact.
- **Budgets.** No function here takes a `Decider`. Builders compare exact rationals, whose signs are known without a decision, and `check` reads facts.
- **Facts.** Builders record facts in the store's `FactBook<GeomKey>`; `check` reads them with `holds`. Pieces 4 and 5 record the facts of the geometry they add, so the bodies they build pass the same check.
- **Marks.** `Body::marks` carries a `MarkLog`. Builders and `empty` give an empty log. Piece 5 fills it for its results.
- **Refusals.** Builders refuse with the codes below. Faults are not refusals: `check` reports them as data.

## Refusal and fault codes

Refusals: `topology.non-positive-size`, `topology.cone-is-cylinder`, `topology.torus-self-intersecting`.

Faults: `topology.body-kind`, `topology.empty-part`, `topology.shared-entity`, `topology.loop-not-closed`, `topology.edge-ends`, `topology.edge-use`, `topology.unbounded-face`, `topology.vertex-loop`, `topology.disconnected-shell`, `topology.euler`, `topology.missing-incidence`.

## Exact public interface

The crate re-exports these types at its root. Everything listed is `pub`; fields not listed are private. The compiled crate in `crates/catallax-topology` is authoritative if a signature listing becomes stale; update this contract with any API change.

```rust
pub type VertexId = Handle<Vertex>;
pub type EdgeId = Handle<Edge>;
pub type CoedgeId = Handle<Coedge>;
pub type LoopId = Handle<Loop>;
pub type FaceId = Handle<Face>;
pub type ShellId = Handle<Shell>;
pub type LumpId = Handle<Lump>;
pub type BodyId = Handle<Body>;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Sense { Same, Reversed }
// Each struct below: #[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Vertex { pub point: PointId }
pub struct Edge { pub curve: CurveId, pub sense: Sense, pub start: Option<VertexId>, pub end: Option<VertexId> }
pub struct Coedge { pub edge: EdgeId, pub sense: Sense }
pub enum Loop { Coedges(Vec<CoedgeId>), Vertex(VertexId) }
pub struct Face { pub surface: SurfaceId, pub sense: Sense, pub loops: Vec<LoopId> }
pub struct Shell { pub faces: Vec<FaceId> }
pub struct Lump { pub shells: Vec<ShellId> }
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum BodyKind { Solid, Empty }
#[derive(Clone, PartialEq, Debug)]
pub struct Body { pub kind: BodyKind, pub lumps: Vec<LumpId>, pub marks: MarkLog }

#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct EntityCounts {
    pub lumps: usize,
    pub shells: usize,
    pub faces: usize,
    pub loops: usize,
    pub edges: usize,
    pub ring_edges: usize,
    pub vertices: usize,
    pub faces_by_kind: BTreeMap<SurfaceKind, usize>,
    pub edges_by_kind: BTreeMap<CurveKind, usize>,
}
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct TopologyFault { pub code: &'static str, pub detail: String }

#[derive(Clone, Debug, Default)]
pub struct Model;
impl Model {
    pub fn new() -> Self;
    pub fn geom(&self) -> &GeomStore;
    pub fn geom_mut(&mut self) -> &mut GeomStore;
    pub fn add_vertex(&mut self, v: Vertex) -> VertexId;
    pub fn add_edge(&mut self, e: Edge) -> EdgeId;
    pub fn add_coedge(&mut self, c: Coedge) -> CoedgeId;
    pub fn add_loop(&mut self, l: Loop) -> LoopId;
    pub fn add_face(&mut self, f: Face) -> FaceId;
    pub fn add_shell(&mut self, s: Shell) -> ShellId;
    pub fn add_lump(&mut self, l: Lump) -> LumpId;
    pub fn add_body(&mut self, b: Body) -> BodyId;
    pub fn vertex(&self, id: VertexId) -> &Vertex;
    pub fn edge(&self, id: EdgeId) -> &Edge;
    pub fn coedge(&self, id: CoedgeId) -> &Coedge;
    pub fn loop_(&self, id: LoopId) -> &Loop;
    pub fn face(&self, id: FaceId) -> &Face;
    pub fn shell(&self, id: ShellId) -> &Shell;
    pub fn lump(&self, id: LumpId) -> &Lump;
    pub fn body(&self, id: BodyId) -> &Body;
    pub fn coedges_of_edge(&self, e: EdgeId) -> Vec<CoedgeId>;
    pub fn owner_of_coedge(&self, c: CoedgeId) -> Option<(LoopId, FaceId)>;
    pub fn faces_of_body(&self, b: BodyId) -> Vec<FaceId>;
    pub fn edges_of_body(&self, b: BodyId) -> Vec<EdgeId>;
    pub fn vertices_of_body(&self, b: BodyId) -> Vec<VertexId>;
    pub fn counts(&self, b: BodyId) -> EntityCounts;
    pub fn check(&self, b: BodyId) -> Vec<TopologyFault>;
}

pub fn block(model: &mut Model, size: [Rational; 3], pose: &Pose) -> Result<BodyId, Refusal>;
pub fn cylinder(model: &mut Model, radius: Rational, height: Rational, pose: &Pose) -> Result<BodyId, Refusal>;
pub fn sphere(model: &mut Model, radius: Rational, pose: &Pose) -> Result<BodyId, Refusal>;
pub fn cone(model: &mut Model, base_radius: Rational, top_radius: Rational, height: Rational, pose: &Pose) -> Result<BodyId, Refusal>;
pub fn torus(model: &mut Model, major: Rational, minor: Rational, pose: &Pose) -> Result<BodyId, Refusal>;
pub fn empty(model: &mut Model) -> BodyId;
```

## What defines correct behaviour

The test writer builds tests from the cases below. Every expected answer follows from this description and maths, so tests compare against it, never against a run of the code.
Tests may use piece 2's types to build hand-made bodies, and piece 2's `==` on surfaces, curves and points to compare fields. Where a case needs a side of a surface, the test evaluates the surface's implicit formula (piece 2's table) itself in exact rationals.

### Builders

1. `block([2, 3, 5], identity)`: `counts` is 1 lump, 1 shell, 6 faces, 6 loops, 12 edges, 0 ring edges, 8 vertices, `{Plane: 6}`, `{Line: 12}`. `check` is empty. The six planes are `==` to the table's; the 8 vertices are the corners; every coedge cycle matches the table up to its starting corner.
2. `cylinder(2, 5, identity)`: 3 faces, 4 loops, 2 edges, both ring edges, 0 vertices, `{Plane: 2, Cylinder: 1}`, `{Circle: 2}`; the side face has two loops; the coedge senses are those of the table; `check` is empty.
3. `sphere(3, identity)`: 1 face with 0 loops, 0 edges, 0 vertices, `{Sphere: 1}`, no edge kinds; `check` is empty.
4. Cones: `cone(2, 1, 2, identity)` has apex `(0, 0, 4)`, axis `(0, 0, -1)`, slope `1/2`; `cone(1, 3, 2, identity)` has apex `(0, 0, -1)`, axis `(0, 0, 1)`, slope 1. Both have 3 faces, 4 loops and 2 ring edges, `{Plane: 2, Cone: 1}`, `{Circle: 2}`. `cone(3, 0, 4, identity)` has apex `(0, 0, 4)`, axis `(0, 0, -1)`, slope `3/4`: 2 faces, 3 loops, 1 ring edge, 1 vertex at `(0, 0, 4)`, `{Plane: 1, Cone: 1}`, `{Circle: 1}`. Each `check` is empty.
5. `torus(3, 1, identity)`: 1 face with no loops, `{Torus: 1}`; `check` is empty (Euler characteristic 0, genus 1).
6. `empty(model)`: kind `Empty`, no lumps, an empty mark log, every count zero, `check` empty.
7. Arena growth: each builder adds exactly the numbers of the table "What each builder adds" to the model and store arenas.
8. Posed builders: `block([1, 1, 1], pose)` with rotation `(1, 0, 0, 1)` and translation `(10, 0, 0)` has the corner `100` at `(10, 1, 0)`, and its face `x = 0` lies on the plane with origin `(10, 0, 0)`, normal `(0, -1, 0)` and `ref_dir` `(-1, 0, 0)`. For random rational poses and every builder, every surface and curve is `==` to the local one transformed by the pose, every vertex point is the moved corner, and `counts` and `check` are those of the identity pose.
9. Outward normals and loop direction, by maths: for each builder at identity and at random poses, a point inside the material (the block's centre, a point on the cylinder's or cone's axis at mid-height, the sphere's centre, the torus's point `(major, 0, 0)` moved by the pose) has `f < 0` for every face whose sense is `Same` (evaluated from the implicit formula). For the block, every loop satisfies the left-hand rule: at each edge, `n x d` points into the face, with `n` the outward normal and `d` the coedge's direction.
10. Facts: after each builder, every fact listed under "Construction facts" holds, and the two false relations listed there do not.
11. Refusals: `block` with a size of 0 or -1 in any of the three places; `cylinder` with radius 0 or height -1; `sphere` with radius 0; `cone` with `base_radius` 0, `height` 0 or `top_radius` -1 give `topology.non-positive-size`. `cone(2, 2, 1)` gives `topology.cone-is-cylinder`; `cone(-1, -1, 1)` gives `topology.non-positive-size` (the order of the checks). `torus` with `minor` 0 gives `topology.non-positive-size`; `torus(1, 1)` and `torus(1, 2)` give `topology.torus-self-intersecting`. After each refusal every arena length and the fact count are unchanged. None panics.
12. Scale: every builder with its sizes multiplied by `10^k` for `k` in `-30..=30` keeps its counts and an empty `check`, and every geometric field scales by `10^k` (squared radii by `10^(2k)`; directions and slopes unchanged).

### Queries

13. On a block: `coedges_of_edge(e)` holds exactly two coedges for each edge, in creation order. `owner_of_coedge(c)` gives the loop and face that hold `c`. A coedge added but not yet in a face has no owner. For a body built by hand that lists one loop in two faces, the owner is the first face added.
14. `faces_of_body`, `edges_of_body` and `vertices_of_body` follow the walk order of "Queries", which the test computes from the body's structure; each entity appears once. `vertices_of_body` of `cone(3, 0, 4)` is the apex alone.
15. A body built by hand from two blocks as two lumps gives 2 lumps, 2 shells, 12 faces, 24 edges and 16 vertices.

### Sound bodies built by hand

These must give an empty `check`. Each needs its incidence facts recorded by the test. Every edge in them has sense `Same`, and coedge senses are given against that.

16. Two lumps touching along a line: a block `[0, 1]^3` and a block `[1, 2] x [1, 2] x [0, 1]` as two lumps of one body. They meet along the segment from `(1, 1, 0)` to `(1, 1, 1)`; each lump has its own edge there, and both edges may use the same line. Counts: 2 lumps, 24 edges, 16 vertices.
17. A lump with a void: the outer shell of the block `[0, 10]^3` and an inner shell on the block `[4, 6]^3`, with every inner face of sense `Reversed` and every inner loop travelled the other way (coedge order reversed and coedge senses swapped). Counts: 1 lump, 2 shells, 12 faces; each shell's Euler characteristic is 2.
18. A sphere of radius 2 split at its equator: the ring edge on the circle with centre `0`, normal `(0, 0, 1)`, `radius_sq` 4, `ref_dir` `(1, 0, 0)`; the north face on the sphere, sense `Same`, with one loop using the ring with coedge `Same`; the south face on the same sphere, sense `Same`, with one loop using the ring with coedge `Reversed`. Two faces on one surface are allowed; merging them is piece 5's concern. Euler characteristic `0 + 1 - 1 + (1 + 1) = 2`.
19. A torus with major 3 and minor 1 cut along one meridian: the ring edge on the circle with centre `(3, 0, 0)`, normal `(0, 1, 0)`, `radius_sq` 1, `ref_dir` `(1, 0, 0)`; one face on the torus with two loops, each holding one coedge of that ring edge, of opposite senses. The one face uses the edge twice. Euler characteristic `0 + 1 - 1 + 0 = 0`.
20. A block with a hole: the block `[0, 4]^3` minus the cylinder of radius 1 about the line `x = 2, y = 2`. Faces: the four side planes as in the block table; the planes `z = 0` and `z = 4` each with two loops, the square and a ring; the hole face on the cylinder with origin `(2, 2, 0)`, axis `(0, 0, 1)`, radius 1, `ref_dir` `(1, 0, 0)`, with sense `Reversed` (its normal points into the hole, out of the material) and two loops. The rings are the circles with centre `(2, 2, 0)` and `(2, 2, 4)`, normal `(0, 0, 1)`, `radius_sq` 1, `ref_dir` `(1, 0, 0)`. Coedge senses on the rings: top plane `Reversed`, hole face at the top `Same`, bottom plane `Same`, hole face at the bottom `Reversed`. Counts: 7 faces, 10 loops, 14 edges with 2 ring edges, 8 vertices, `{Plane: 6, Cylinder: 1}`, `{Line: 12, Circle: 2}`. Euler characteristic `8 + 2 - 14 + (4 + 0 + 0 + 0) = 0`, genus 1.

### Faulty bodies built by hand

Each is built from new entities, mostly reusing a builder's. Each must give the listed code among its faults.

21. A block's shell without one face: `topology.edge-use` (four edges with one coedge each), and `topology.euler` (characteristic 1).
22. A block face rebuilt with one coedge's sense flipped: `topology.loop-not-closed` and `topology.edge-use`.
23. A block face whose loop lacks one coedge: `topology.loop-not-closed`.
24. A block's shell listing one face twice: `topology.shared-entity`. A seventh face that reuses an existing face's loop: `topology.shared-entity`.
25. A ring edge on a line, and an edge with a `start` and no `end`: `topology.edge-ends`.
26. A plane face with no loops: `topology.unbounded-face`.
27. A face whose only loop is a vertex loop: `topology.vertex-loop`.
28. A cylinder rebuilt with its top ring on a new circle for which no fact was recorded: `topology.missing-incidence`.
29. A shell of two loopless faces on one sphere: `topology.euler` (characteristic 4) and `topology.disconnected-shell`.
30. A shell holding all twelve faces of two separate blocks: `topology.disconnected-shell`.
31. A body whose two lumps list the same shell, and a body listing one lump twice: `topology.shared-entity`.
32. A `Solid` body with no lumps, and an `Empty` body with one lump: `topology.body-kind`. A lump with no shells: `topology.empty-part`.

### Determinism

33. The same sequence of builder calls on two new models gives the same handles and `==` entities, geometry and facts. Running it on eight threads, one model each, gives identical models.

### Volume and running time of the tests

The generated families (cases 8, 9, 12 and 33) must together run at least 10,000 builder-and-check questions per test run, from fixed seeds.
Keep the default test run of this piece under about 5 minutes in release mode on one core. Everything here is exact rational arithmetic and handle bookkeeping, so this is a generous limit.

## Speed: measured, never gated

No test fails for speed. A benchmark target measures and prints:

1. the time to build each primitive, at identity and at a random pose;
2. the time of `check`, `counts`, `faces_of_body`, `edges_of_body` and `vertices_of_body` against face count, on prisms built by hand: the prism of height 1 over the convex polygon with vertices `(i, i^2, 0)` for `i = 0, ..., n - 1`, with `n + 2` faces, `3n` edges and `2n` vertices, for `n` up to 10,000;
3. the time of `coedges_of_edge` and `owner_of_coedge` on the largest prism.

The numbers are written to the benchmark output so later runs can be compared.

## Out of scope for this piece

- Every geometric decision: sides, intersections, containment (pieces 2, 4 and 5).
- Splitting and merging faces and edges, and building bodies from booleans (piece 5).
- Checking that a vertex lies on its surfaces by numbers; `check` trusts the recorded facts.
- Sheets, wires, meshes and lattices: later `BodyKind` variants.
