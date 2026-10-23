# Piece 6: measurement (`catallax-measure`)

Status: complete for blind test writing and blind building.

This is the input the test writer and the creator get, together with `docs/DESIGN.md`, `docs/ARCHITECTURE.md` and the descriptions of pieces 1 to 5 in `docs/pieces/`. This piece reads piece 3's `Model`, evaluates piece 2's surface parametrisations and piece 4's traces, and uses piece 1's `Ball` arithmetic, outcomes and marks. Bounds reuse piece 5's exact trimmed-face classifier on a private copy of the model.

## Approved algebraic-definition migration

The 2026-10-22 architecture decision widens surface constants, frame transformations and curve forms to `Real`. Compile exact expressions into the existing certified integration path; coefficient enclosures contribute to every remainder and property bound. Keep the rational path and public `Accuracy`, `Enclosure` and mass-property interfaces unchanged. No conversion to rational approximations is permitted inside a proof.

Add elliptic-cylinder support using piece 2's chart and derivatives. Its area density is `sqrt(major_sq sin(u)^2 + minor_sq cos(u)^2)`; compute its perimeter as a certified integral, without adding transcendental nodes to `Real`. Integrate over the actual oblique cap/edge trims. For an extruded planar region of area `A`, use `A * abs(v . n)` as an independent volume check, where `v` is the exact sweep vector and `n` the unit profile normal. A full right elliptic cylinder provides the independent volume check `pi * sqrt(major_sq * minor_sq) * h`.

Use certified closure and finite conic traces from piece 4. A parabola or hyperbola support is unbounded but a solid's face/edge domain is bounded; integrate that domain. Preserve seam and singularity treatment. Add algebraic placement invariance, open-conic trims and oblique circular extrusions to the independent tests before claiming feature coverage. Existing budget and accuracy failure policies remain unchanged and remain separate from PRD completion.

## Original-target properties of certified offset prisms

`elliptic_offset_prism_mass_properties(definition, axial, accuracy, dec)` and
`elliptic_offset_prism_volume_area_centroid(...)` integrate the retained
`EllipticOffsetDefinition` and exact axial interval. They measure the nominal
product solid, not its fitted B-spline carrier or display mesh. The three-
property path does not compute inertia. An empty interval denotes empty
material; nonregular inward profiles refuse.

The existing validated Taylor integrator encloses the exact principal-quarter
normal-offset trace. Reflection symmetry supplies the whole profile; Green's
formula supplies its area and second moments, and the trace speed supplies
perimeter. Product-solid formulas give volume, total cap/side area, centroid
and centroidal unit-density inertia in the original exact frame, then transport
the tensor to model axes. No sampled fit or floating-point quadrature result
becomes a property enclosure.

These definition-level APIs do not themselves validate a body certificate.
`CertifiedApplicationBody` binds them to a verified immutable body and target.
Requested property widths stay independent of the carrier's local tolerance;
the carrier error is not added to make a nominal-property request pass.
Library accuracy marks retain their existing meaning; the opt-in manifest/4
consumer reports `measure.accuracy-not-reached` rather than claiming a requested
width it did not achieve.


## What this piece is for

Piece 6 computes body mass properties, face areas and centroids, edge lengths and half-arc-length points, and body and face axis-aligned bounds. Results are proven enclosures with exact rational ends for the geometry as built.
It works from exact definitions only: the surfaces' parametrisations and their derivatives (piece 2), the certified traces of edge arcs (piece 4), and `Ball` enclosures of `pi` and elementary functions (piece 1). It never measures a mesh or a sample.

## Number types

| What | Type | Meaning |
| --- | --- | --- |
| Results | `Enclosure { lo, hi }`, both `Rational` | `lo <= true value <= hi`, proven |
| Requested accuracy | `Accuracy { relative, absolute }`, both `Rational` | The width asked for |
| Work inside the computation | `Ball` | Proven enclosures at a working precision this piece chooses |
| Parameter values, angles | `Ball` | Caches for integration; never a decision |

## What it measures

For a body `b` and its material `M` (the closure of the open set its shells bound, piece 5):

- **volume**: the volume of `M`;
- **area**: the sum of the areas of the faces of `b`, each counted once, in `faces_of_body` order;
- **centre of mass**: `(1 / volume) * integral of p over M`, the centroid of uniform density;
- **inertia**: the unit-density tensor about the centroid `c`, with diagonal entry `integral ((p_j-c_j)^2 + (p_k-c_k)^2)` and off-diagonal entry `-integral (p_i-c_i)(p_j-c_j)`, expressed in model axes and length units to the fifth power;
- **face area and centroid**: area and first area moments of the part of the surface inside its loops;
- **edge length and midpoint**: trimmed arc length and its half-arc-length point;
- **bounds**: inclusive axis-aligned limits for the exact body or face, without display padding.

The model-based APIs below measure the body as built. When that body carries marks (it was built through not-exact decisions), the enclosures are proven for that body, which may differ from the true answer within those marks (`JUDGE-INTERFACE.md`, ND10). This piece does not add the body's marks to its answer; they stay on the body, and the front end reports both. The certified-prism APIs above instead measure the explicitly retained original target and keep its represented-model error separate.

An `Empty` body has volume and area exactly 0 (`lo = hi = 0`), centroid `None` and inertia `None`.
A `Solid` body whose faces lie on planes, circular or elliptic cylinders, cones, spheres or tori, and whose edges lie on curves that piece 4 can trace, is measured. That covers every supported primitive and boolean result; uncovered surfaces retain `measure.unsupported-surface`.

## Rules

1. **Containment.** Every result encloses the true value: `lo <= value <= hi`. A true value outside its enclosure is a wrong proven answer (DESIGN.md, R2).
2. **Accuracy.** `Accuracy { relative, absolute }` asks for a width `hi - lo` of at most `max(relative * m, absolute)`, where `m` is the smaller of `|lo|` and `|hi|` when `lo` and `hi` have the same sign, and 0 otherwise. The default is `relative = 10^-12`, `absolute = 10^-30`. The rule applies to each centroid coordinate and each tensor entry.
3. **When the budget does not allow it.** The answer is still a proven enclosure, only wider, and it carries one mark per property that missed its width: code `measure.accuracy-not-reached`, kind `BoundKind::Value`, bound the width reached. `mass_properties` adds at most four marks, in the order volume, area, centroid, inertia. Centroid and inertia each contribute one mark with the largest entry width if any entry missed its target.
4. **Exact results.** The enclosure has `lo == hi` in exactly these cases, where the property is rational and the computation can be done in rationals:
   - the volume, centroid and inertia of a body whose faces are all planes and whose vertices are all `Explicit` rational points;
   - the area of such a body when every plane's `normal` has a rational length (true for every posed block, since rotations keep the unit normals unit).
   Every other enclosure, including those of piece 5's planar results whose vertices are meeting points, has `lo < hi` or, by chance, `lo == hi`; tests check containment, not `lo < hi`.
5. **Centroid.** `None` exactly when the body is `Empty`. For a `Solid` body the volume enclosure must exclude zero before dividing; this piece refines it until it does, ignoring the `absolute` floor for that purpose. When the budget runs out first, the answer is refused with `measure.centroid-not-separated`.
6. **Determinism.** The same body, accuracy and policy give the same enclosure ends on every machine and every run (DESIGN.md, R10).

## Budgets and the deadline

This piece's arithmetic is `Ball` arithmetic, which piece 1 does not count. It counts its own work the same way on every machine: one `Ball` operation at working precision `p` costs `ceil(p / 64)^2` units. It never uses more than `policy().max_precision_bits` bits or more than `policy().max_work` units on one property, where `policy()` is the `Decider`'s. Exact rational arithmetic is not counted. Whether the `Ball` operations inside piece 4's `TracePiece` methods count is the creator's choice; tests do not measure the count. When the next refinement would pass either limit, it stops and returns the enclosure reached with the mark of rule 3. When no enclosure was reached at all (for example under `max_work = 0` for a property that needs `Ball` operations), that property is refused with `measure.budget-exhausted`.

The exact quantities it needs from other pieces (vertex positions through `GeomStore::enclose_point`, arcs through piece 4's `trace`) go through the `Decider`, which counts them and refuses them after the deadline. So after the deadline a measurement refuses with `numbers.deadline` at its first such call, and one that needs none (a block with `Explicit` vertices, a builder sphere, whose properties need only `Ball::pi`) may answer as usual (ARCHITECTURE.md, "The time budget").

Inertia uses the radial divergence field `F = p p_i p_j / 5` for each raw second moment, since `div F = p_i p_j`. The six scalar refinements each receive one sixth of the tensor's Ball-work allowance. Exact planar moments use rational triangle integrals. Curved moments use the same analytic forms, whole-surface corrections and certified quadrature as volume and centroid. The centroid enclosure's midpoint supplies a translated integration frame; interval parallel-axis subtraction retains its uncertainty. A missed final tensor width produces the inertia accuracy mark, never an unproved point value.

## How the integrals may be computed

The method is the creator's choice, within the rules above. A method that fits the definitions:

1. By the divergence theorem, the volume is `(1/3)` times the sum over faces of the integral of `p . n` (with `n` the outward unit normal: the surface normal times the face's sense), and the first moments are `(1/2)` times the sums of the integrals of `p_i^2 n_i`.
2. Each face integral is an integral over the face's region of the surface's `(u, v)` domain, with the area element `|dS/du x dS/dv|` from piece 2's `derivatives`. By Stokes's theorem on the surface, it equals the integral, around the face's loops, of a 1-form whose exterior derivative is the integrand. For the analytic surfaces such 1-forms have closed forms.
3. Each loop integral is a sum over coedges of integrals along piece 4's traces of the edge arcs, in the coedge's direction. Validated quadrature over each trace piece uses `taylor` enclosures: a Taylor expansion of the integrand on a sub-range, integrated term by term, plus a proven remainder from the highest-order enclosure over the whole sub-range.

Two traps decide whether such a method is right:

- **Seams.** The angle parameter `u` of a cylinder, cone, sphere or torus jumps by `2 pi` at the seam (piece 2, `parameters_of`). A 1-form must be single-valued on the surface: use `du`, which is defined everywhere off the axis or poles, and never `u` itself, or add the seam's jump explicitly. A cylinder band whose loops are two rings has no edge on the seam at all, and its area must still come out as `2 pi r h`.
- **Poles and apexes.** At a sphere's poles, a cone's apex and on a cylinder's axis, `u` is undefined, so a 1-form built on `du` is singular there. A face that contains a pole (a whole sphere, a hemisphere, a sphere cap) needs the limit term of a small circle around the pole; a cone face with a vertex loop at the apex needs the same at the apex.

"Along the seam" regression poses test exactly these traps: faces whose edges lie on a seam, or cross it, or pass through a pole.

Complete torus parallel circles use the exact Fourier mean of the boundary form over the azimuth, rather than traced quadrature around a full period. On a certified cone generator through the apex, `q × dq = 0` makes all five boundary forms identically zero; that proof removes the apparent `1/z` endpoint singularity. Neither path substitutes a sampled integral or drops a nonzero boundary term.

## How this piece uses the earlier pieces

- **Outcomes.** Every function returns `Answered`. Refusals from geometry, tracing and exact trim membership pass through unchanged.
- **Budgets.** As stated above. This piece never creates a `Decider` or a policy.
- **Facts.** Integral measurements read no facts and record none. Bounds delegate exact trim predicates to piece 5 on a private model copy; the caller's facts remain unchanged.
- **Marks.** Only `measure.accuracy-not-reached`, built with `Mark::new`. A `trace` from piece 4 never carries marks.
- **Geometry.** Surface points and derivatives from piece 2's `Surface::evaluate` and `derivatives`; parameters of edge points from `parameters_of` when a method needs them; arcs from piece 4's `trace` with the coedge's edge and its `start` and `end` vertices' points (a ring edge is traced with `None, None`); vertex-loop points from `GeomStore::enclose_point`.

## Refusal codes

| Code | When |
| --- | --- |
| `measure.not-solid` | the body's kind is neither `Solid` nor `Empty` (piece 3 has no such kind yet, so no test reaches it) |
| `measure.unsupported-surface` | a face's surface kind this piece does not measure (none today; the enum is non-exhaustive) |
| `measure.centroid-not-separated` | the normalizing volume or face area enclosure could not be separated from zero |
| `measure.budget-exhausted` | no enclosure of a property was reached within the budget |
| `measure.bounds-not-separated` | an exact extrema ordering or membership predicate is unresolved |
| `measure.unbounded-face` | a face or edge has no finite trimmed extrema |
| `measure.unsupported-curve` | a future curve kind has no extrema definition |

Refusals of pieces 1, 2, 4 and the reused piece 5 trim classifier pass through unchanged.

## Programming errors

These may panic. None can come from a question file:

- a `BodyId` or `FaceId` from another model, or out of range;
- a body that fails `Model::check`;
- an `Accuracy` with a negative part.

## Exact public interface

The crate re-exports these types at its root. Everything listed is `pub`; fields not listed are private. The compiled crate in `crates/catallax-measure` is authoritative if a signature listing becomes stale; update this contract with any API change.

```rust
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Enclosure { pub lo: Rational, pub hi: Rational }
impl Enclosure {
    pub fn width(&self) -> Rational;                  // hi - lo
    pub fn contains(&self, value: &Rational) -> bool; // lo <= value <= hi
    pub fn is_exact(&self) -> bool;                   // lo == hi
}
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Accuracy { pub relative: Rational, pub absolute: Rational }
// Accuracy: Default (relative 10^-12, absolute 10^-30)
pub type InertiaTensor = [[Enclosure; 3]; 3];
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct MassProperties { pub volume: Enclosure, pub area: Enclosure, pub centroid: Option<[Enclosure; 3]>, pub inertia: Option<InertiaTensor> }
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct VolumeAreaCentroid { pub volume: Enclosure, pub area: Enclosure, pub centroid: Option<[Enclosure; 3]> }

pub fn mass_properties(model: &Model, body: BodyId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<MassProperties>;
pub fn volume_area_centroid(model: &Model, body: BodyId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<VolumeAreaCentroid>;
pub fn centroid(model: &Model, body: BodyId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<Option<[Enclosure; 3]>>;
pub fn inertia(model: &Model, body: BodyId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<Option<InertiaTensor>>;
pub fn volume(model: &Model, body: BodyId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<Enclosure>;
pub fn area(model: &Model, body: BodyId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<Enclosure>;
pub fn face_area(model: &Model, face: FaceId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<Enclosure>;
pub fn face_centroid(model: &Model, face: FaceId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<[Enclosure; 3]>;
pub fn edge_length(model: &Model, edge: EdgeId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<Enclosure>;
pub fn edge_midpoint(model: &Model, edge: EdgeId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<[Enclosure; 3]>;
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundingBox { pub min: [Rational; 3], pub max: [Rational; 3] }
pub fn bounding_box(model: &Model, body: BodyId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<Option<BoundingBox>>;
pub fn face_bounding_box(model: &Model, face: FaceId, accuracy: &Accuracy, dec: &mut Decider) -> Answered<BoundingBox>;
```

`volume`, `area`, `volume_area_centroid` and `mass_properties` give the same volume and area enclosures for the same body, accuracy and fresh `Decider`. `volume_area_centroid` retains the original three-property computation without requesting inertia. `mass_properties` requests all four properties. Application builds do not measure automatically; `ApplicationBody::mass_properties` and `ApplicationBody::inertia` borrow an already built result.

### Face, edge and bounds queries

`face_centroid` divides the surface's first area moments by its area, independent
of material orientation. Rational polygonal faces use exact triangle moments.
Curved faces use boundary forms, including the sphere pole and torus holonomy
corrections used by body measurement. Area must exclude zero before division.

`edge_length` measures the trimmed arc. `edge_midpoint` is its half-arc-length
point, not the average of endpoints or a parameter midpoint. Closed edges start
at their vertex, or at the carrier seam when vertexless. Lines and whole circles
have direct formulas. Other arcs use certified speed quadrature; midpoint panels
enclose every parameter interval whose cumulative length can contain half the
total. Budget exhaustion retains the last proven length or midpoint enclosure
with its accuracy mark.

Bounds have inclusive rational limits and no display chord padding. An empty
body has no box. Bounds select vertices, exact axis extrema of conic edges and
interior points whose surface normal is parallel to an axis. Compact carriers
with trims use piece 5's exact `face_contains` classifier on a private model copy;
the query does not change the caller's geometry or topology. That classifier
accepts exact `Point3` as well as the existing rational points. Measurement's
dependency on the Boolean crate is acyclic.

General intersection-branch extrema use endpoint ranges padded by a certified
second-derivative interpolation error. Only panels that can improve an unmet
extremum refine. Each request shares one deterministic arithmetic allowance,
including body requests over multiple faces. If this fallback exhausts it,
`measure.budget-exhausted` refuses the query rather than returning an unfinished
refinement. An unsettled exact candidate predicate returns
`measure.bounds-not-separated`; an unbounded face returns
`measure.unbounded-face`. Numeric extrema enclosures that miss the requested
width carry the ordinary accuracy marks. No query measures display triangles.

## What defines correct behaviour

The test writer builds tests from the cases below. Every expected value is a closed form, so tests check that the enclosure contains it, using their own independent high-precision evaluation of the closed form (an interval or ball library used only as a dev-dependency, or a table from a separate program). A test of containment against a transcendental value encloses that value first and checks that the enclosures overlap and that the closed form's enclosure, narrowed until it is much narrower than the result, lies inside it.

Bodies come from piece 3's builders, or are built by hand with piece 3's `add_*` functions as described, with every coedge sense following piece 3's orientation rule and the incidence facts recorded. Results of piece 5 are listed for the cross-check once piece 5 exists; they are not needed to test this piece.

### Primitives

1. `block([a, b, c], identity)` for `(a, b, c) = (2, 3, 5)` and `(1/3, 7, 10^-3)`: volume `abc`, area `2 (ab + bc + ca)` and centroid `(a/2, b/2, c/2)`, all exact (`lo == hi`).
2. `cylinder(2, 5, identity)`: volume `20 pi`, area `28 pi`, centroid `(0, 0, 5/2)`.
3. `sphere(5, identity)`: volume `500 pi / 3`, area `100 pi`, centroid `0`. Its face has no loops, and its poles and seam are inside the face.
4. `cone(2, 1, 2, identity)`: volume `14 pi / 3`, area `(5 + 3 sqrt(5)) pi`, centroid `(0, 0, 11/14)`. `cone(3, 0, 4, identity)`: volume `12 pi`, area `24 pi`, centroid `(0, 0, 1)`; the apex is a vertex loop. In general, for base radius `R`, top radius `r` and height `h`: volume `pi h (R^2 + R r + r^2) / 3`, area `pi (R^2 + r^2) + pi (R + r) sqrt((R - r)^2 + h^2)`, centroid height `h (R^2 + 2 R r + 3 r^2) / (4 (R^2 + R r + r^2))`.
5. `torus(3, 1, identity)`: volume `6 pi^2`, area `12 pi^2`, centroid `0`. Both of its seams are inside its one loopless face.
6. `empty(model)`: volume and area exactly 0, centroid `None`.

### Bodies built by hand

7. The block with a through hole of piece 3's case 20 (`[0, 4]^3` minus the cylinder of radius 1 about `x = 2, y = 2`): volume `64 - 4 pi`, area `96 + 6 pi`, centroid `(2, 2, 2)`.
8. The block with a void of piece 3's case 17 (`[0, 10]^3` with the void `[4, 6]^3`): volume 992, area 624, centroid `(5, 5, 5)`, all exact.
9. The upper half ball of radius `R = 2`: the sphere face (centre `0`, radius 2, axis `(0, 0, 1)`, `ref_dir (1, 0, 0)`, sense `Same`) with one loop, the equator ring edge (circle with centre `0`, normal `(0, 0, 1)`, `radius_sq` 4, `ref_dir (1, 0, 0)`, edge sense `Same`) with coedge `Same`; and the disc face (plane with origin `0`, normal `(0, 0, -1)`, `ref_dir (1, 0, 0)`, sense `Same`) with one loop, the same ring with coedge `Reversed`. Volume `2 pi R^3 / 3 = 16 pi / 3`, area `3 pi R^2 = 12 pi`, centroid `(0, 0, 3R/8) = (0, 0, 3/4)`. The north pole is inside the curved face.
10. A sphere cap: the ball of radius `R = 5` cut by the plane `z = 3`, keeping `z >= 3` (cap height `h = 2`): the sphere face (centre `0`, axis `(0, 0, 1)`) with one loop on the circle with centre `(0, 0, 3)`, normal `(0, 0, 1)`, `radius_sq` 16, and the disc on the plane with origin `(0, 0, 3)`, normal `(0, 0, -1)`. As in case 9, every surface and circle has `ref_dir (1, 0, 0)`, the ring edge has sense `Same`, its coedge is `Same` on the sphere face and `Reversed` on the disc, and both faces have sense `Same`. Volume `pi h^2 (3R - h) / 3 = 52 pi / 3`, area `2 pi R h + pi (2 R h - h^2) = 36 pi`, centroid height `3 (2R - h)^2 / (4 (3R - h)) = 48/13`.
11. The octant of the ball of radius `R = 3` (`x, y, z >= 0`): three quarter-disc faces on the planes `x = 0`, `y = 0`, `z = 0` with outward normals `(-1, 0, 0)`, `(0, -1, 0)`, `(0, 0, -1)`; the sphere face; three edges on the coordinate axes from `0` to `(3, 0, 0)`, `(0, 3, 0)`, `(0, 0, 3)`; three quarter-circle edges, each on the circle of `radius_sq` 9 centred at `0` in one coordinate plane; vertices `0`, `(3, 0, 0)`, `(0, 3, 0)`, `(0, 0, 3)`. Volume `pi R^3 / 6 = 9 pi / 2`, area `5 pi R^2 / 4 = 45 pi / 4`, centroid `(3R/8, 3R/8, 3R/8) = (9/8, 9/8, 9/8)`.
12. Along the seam, the half cylinder `y >= 0` of radius `r = 2` and height `h = 3` on the `z` axis: the cylinder face (origin `0`, axis `(0, 0, 1)`, radius 2, `ref_dir (1, 0, 0)`), which runs from the seam line `(2, 0, z)` to the line `(-2, 0, z)`; the rectangle on the plane `y = 0` (normal `(0, -1, 0)`); the half discs on `z = 0` and `z = 3`. Volume `pi r^2 h / 2 = 6 pi`, area `pi r h + 2 r h + pi r^2 = 10 pi + 12`, centroid `(0, 4r / (3 pi), h/2) = (0, 8 / (3 pi), 3/2)`. The half cylinder `y <= 0`, whose face runs from `(-2, 0, z)` to the seam: the same volume and area, centroid `(0, -8 / (3 pi), 3/2)`.
13. The Steinmetz solid of two cylinders of radius `r = 1` with axes `z` and `x` through `0`: four faces (two on each cylinder), four edges on the two ellipses of piece 4's case 21 (each ellipse split at `(0, 1, 0)` and `(0, -1, 0)`), two vertices. Volume `16 r^3 / 3 = 16/3`, area `16 r^2 = 16`, centroid `0`. The edges are ellipse arcs. With `ref_dir (1, 0, 0)` on the `z` cylinder and `(0, 0, -1)` on the `x` cylinder (the frame `xcyl` of piece 5 gives), the face `x > 0` of the first and the face `z < 0` of the second each contain a stretch of their cylinder's seam.
14. `face_area` of each face of cases 1 to 13, against the closed form of that face (for example `pi r^2` for a cylinder cap, `2 pi r h` for its side, `4 pi R^2` for a whole sphere, `pi R^2 / 4` for a quarter disc), and the sum of the face areas against the body's area.

### Planted faults

15. No test plants a fault. A body that fails `Model::check` is a programming error here, and a face whose sense disagrees with its loops passes `check` but bounds no solid, so its measurement is not specified. That the measurement reads the body is checked instead by the hand-built bodies of cases 9 to 13, which no builder makes.

### Accuracy, budgets and scale

16. For the primitives of cases 2 to 5 and the bodies of cases 9 to 13, with `relative` from `10^-6` to `10^-30` and `absolute = 10^-40`: each enclosure contains the closed form, and its width meets rule 2. The enclosures for tighter accuracies need not lie inside those for looser ones.
17. Under `BudgetPolicy { max_precision_bits: 64, max_work: 10_000 }`, with the default accuracy: each property either contains the closed form, carrying `measure.accuracy-not-reached` with the width reached as the bound when it misses the width, or is refused by piece 1, 2 or 4 (for example `trace` when locating a vertex ends `NotExact`). Under `BudgetPolicy { max_precision_bits: 64, max_work: 0 }`: the bodies of rule 4 answer exactly, since rational arithmetic is not counted; every other property is refused with `measure.budget-exhausted` or by piece 1, 2 or 4. Never a false enclosure.
18. Scale (DESIGN.md, R9): cases 1 to 13 with every length multiplied by `10^k` for `k` in `-30..=30`: volumes scale by `10^(3k)`, areas by `10^(2k)`, centroids by `10^k`, and every enclosure contains the scaled closed form. Placement: the same bodies moved by 20 random rational poses (integer quaternion components up to 100, translations up to `10^6`): volume and area unchanged, centroid moved by the pose.
19. Determinism: the same body and accuracy give identical enclosure ends with a fresh `Decider`, on repeat and on eight threads.
20. A `Decider` whose deadline has passed: each measurement of cases 1 to 13 either refuses with `numbers.deadline` or gives its default-budget answer, never another (ARCHITECTURE.md, "The time budget").

### Results of piece 5, for the cross-check

These are the closed forms of piece 5's cases, for the integration tests once piece 5 is built:

| Piece 5 case | Body | Volume | Area | Centroid |
| --- | --- | --- | --- | --- |
| 1 | `A ∪ B` | 15 | 42 | `(3/2, 3/2, 3/2)` |
| 1 | `A ∩ B` | 1 | 6 | `(3/2, 3/2, 3/2)` |
| 1 | `A − B` | 7 | 24 | `(13/14, 13/14, 13/14)` |
| 10 | `A − B` | `1000 - 90 pi` | `600 + 42 pi` | `(5, 5, 5)` |
| 10 | `A ∩ B` | `90 pi` | `78 pi` | `(5, 5, 5)` |
| 14 | `A − B`, four lumps | `1000 - 250 pi` | `600 + 50 pi` | `(5, 5, 5)` |
| 24 | `A ∩ B` | `16/3` | 16 | `0` |
| 30 | `A ∩ B` | `244 pi / 3` | `68 pi` | `0` |
| 41 | `(A − B) − C` | `1000 - 80 pi + 128/3` | | |

21. How the table's values follow. Case 1: the cubes have volume 8 and overlap in a unit cube, so the union is `8 + 8 - 1`; each cube loses three unit squares of its surface to the other, so the union's area is `24 + 24 - 6` and the difference's is `24 - 3 + 3`; the centroids are volume-weighted means of `(1, 1, 1)`, `(2, 2, 2)` and the overlap's `(3/2, 3/2, 3/2)`. Case 10: the hole removes `9 pi * 10` of volume, `2 * 9 pi` of top and bottom, and adds the wall `2 pi * 3 * 10`. Case 14: the four side faces keep their full area 400, the top and bottom keep `2 (100 - 25 pi)`, and the inner wall adds `2 pi * 5 * 10`. Case 30: a cylinder of radius 3 between `z = -4` and `z = 4` with two spherical caps of height 1 on the ball of radius 5, so `72 pi + 2 pi (1^2)(15 - 1) / 3` and `48 pi + 2 * (2 pi * 5 * 1)`. Case 41: the two holes of radius 2 and length 10 overlap in a Steinmetz solid of volume `16 * 2^3 / 3`.

### Volume and running time of the tests

Cases 16, 18 and 19 are the generated families. They must together run at least 3,000 maths-answered enclosures per test run, from fixed seeds. Measurements at `10^-30` are expensive; keep the default test run under about 10 minutes in release mode on one core, and run the tightest accuracies on a few bodies only.

## Values checked for this description

The architect checked these closed forms with `sympy` and 30-digit decimals in an independent program: the cone frustum's volume, area and centroid height (by integrating the cross-section), `14 pi / 3`, `(5 + 3 sqrt(5)) pi = 36.78240746107114537178...` and `11/14` for `cone(2, 1, 2)`, `12 pi`, `24 pi` and 1 for `cone(3, 0, 4)`; the half ball's `3R/8`, the cap's `3 (2R - h)^2 / (4 (3R - h))`, the half cylinder's `4r / (3 pi)`, the octant's `3R/8`, the Steinmetz volume 16/3 (by integrating square cross-sections); `1000 - 90 pi = 717.2566611769186085383620955...`, `600 + 42 pi = 731.9468914507713160154310220...` and `1000 - 80 pi + 128/3 = 791.3392543794832075896551960...`.

## Speed: measured, never gated

No test fails for speed. A benchmark target measures and prints:

1. time of each property against requested accuracy (`10^-6` to `10^-30`) for each primitive and the bodies of cases 7 to 13;
2. the `Decider` counters each call adds (`DecisionStats`, including its `max_bits`), the wall time, and the width of every enclosure returned. The precision and `Ball` work of measurement's own integration are not reached through the `Decider`, so no public interface can observe them; the benchmark records them as not observable until an observation interface exists;
3. time for the results of piece 5's first-milestone cases, once they exist.

The numbers are written to the benchmark output so later runs can be compared.

## Out of scope for this piece

- Moments of inertia and other integrals (later, on the same machinery).
- Measuring meshes or sampled geometry, ever.
- Any tolerance. This piece has none.
