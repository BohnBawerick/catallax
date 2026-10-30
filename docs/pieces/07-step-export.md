# Piece 7: STEP export (`catallax-step`)

Status: complete for blind test writing and blind building.

This is the input the test writer and the creator get, together with `docs/DESIGN.md`, `docs/ARCHITECTURE.md`, `docs/JUDGE-INTERFACE.md` (the STEP file) and the descriptions of pieces 1 to 4 in `docs/pieces/`. This piece reads piece 3's `Model`, piece 2's definitions and piece 4's closed forms and traces, and uses piece 1's numbers, outcomes and marks exactly as those descriptions define them. Nothing else is needed.

## Approved algebraic-definition migration

The 2026-10-22 architecture decision preserves the exchange text, topology conventions, entity ordering and report fields below. It widens the exact data supplied to the exporter. Algebraic positions, directions, radii and conic dimensions pass through the existing certified decimal enclosure path. Frame normalisation and every coefficient's rounding contribute to the deviation over the finite face/edge extent. No decimal reconstruction becomes the kernel's definition.

Extend the geometry mapping as follows:

| Exact support | Export representation |
| --- | --- |
| Elliptic cylinder | `SURFACE_OF_LINEAR_EXTRUSION` of the transverse `ELLIPSE` in piece 2's frame, with extrusion vector along its unit axis; retain the actual face trims |
| Parabola | `PARABOLA` at the vertex, reference along the opening axis, with focal distance `focal` |
| One hyperbola branch | `HYPERBOLA` with semi-axis `sqrt(transverse_sq)` and semi-imaginary-axis `sqrt(conjugate_sq)`; orient the placement to select the stored branch and preserve the edge's natural direction |

These STEP entity definitions are documented by STEP Tools for [linear extrusion](https://steptools.com/docs/stp_aim/html/t_surface_of_linear_extrusion.html), [parabola](https://steptools.com/docs/stp_aim/html/t_parabola.html) and [hyperbola](https://www.steptools.com/docs/stp_aim/html/t_hyperbola.html). Verify the chosen placements and bounded trims by reading back independent fixtures. If an entity is unsuitable for the selected AP242 representation, use a certified rational spline representation of that finite patch/arc; an unproved sampled fit is not an alternative. An exact analytic representation still has decimal-rounding error.

Open curves must have finite edge endpoints. Existing ring handling applies only after piece 4 proves closure. Keep the existing entity sequence and bytes for the successful rational corpus during conic preparation. The new geometry must not change manifest kind tags or `JUDGE-INTERFACE.md`. The existing tolerance/report policy below remains the baseline; this migration does not settle the separate PRD-wide 1 nm enforcement gap.

## Certified spline-prism export

`write_spline_prism_step(model, body, &CertifiedSplineStepOptions, dec, out)`
is the separate exporter for the checked first spline-prism family. It uses
`SplinePrism::check` before constructing the file and writes a genuine
degree-(3,1) `B_SPLINE_SURFACE_WITH_KNOTS`, two shared degree-3 boundary curves
and planar cap faces. Surface rows and edge curves share the same rounded
controls. Model face names are `face_<handle>`; export-only ring vertices are
reported separately from the model's vertex-free topology.

The finite decimal construction uses an exact rational integer-direction frame,
common rounded Hermite endpoints and tangents, and an exact sweep vector.
Shared knot lengths retain exact incidence and C1 joins, including periodic
closure. The exporter rechecks regular convex embedding and cap compatibility
after rounding. Its whole-domain two-way file-to-model bound includes every
control displacement and covers the complete product solid; it is not a
sample maximum. Other frames/knot formats remain explicit coverage gaps.

`file_tolerance` is a nonnegative required bound in mm.
`decimal_places` controls the common origin grid, not a promise that every
coefficient has that many decimals: exact directions and knot lengths can
require more. A coarse grid or an unproved bound refuses before writing bytes,
including `step.tolerance-not-reached`; the exporter never enlarges the request.
`StepReport.max_deviation` concerns the represented model only. The certified
application exchange separately composes it with the body's proved
model-to-definition bound and binds both to physical regions.

The existing `write_step` contract, analytic entity ordering and retained
legacy bytes below are unchanged. It must not relabel a spline as an ellipse
or silently use this new export contract. Independent native STEP import checks
topology separately from the exact finite-decimal proof. An external kernel's
nonadaptive volume estimate is not a certified property enclosure.


## Rolling-envelope junction export

`write_rolling_junction_step(model, body, &CertifiedSplineStepOptions, dec, out)`
writes a certified body with rolling-envelope faces through the ordinary exporter
and refuses with `step.tolerance-not-reached` before any byte leaves if
`max_deviation` exceeds `file_tolerance`. Only a body with such a face reaches
this plan, so every other file keeps its bytes.

Each envelope is one complex rational degree-(8,2) B-spline surface. Its four
boundary edges are complex rational B-spline curves over that surface's own
rounded control rows, so a shared row cannot drift from the surface. A reader
without pcurves cannot place a curve on a curved face, so each row also has a straight
line on the envelope and a cubic pcurve on its analytic neighbour.
The neighbour's parameter is an arctangent of the row's parameter, which no
polynomial matches. `RollingEnvelope::parameter_table` fits it with cubic Hermite
pieces from interval jets at 10^-11 rad, and the bound times the chart's radius
joins `max_deviation`.

Cylinder and cone faces beside a junction end in chains of arcs and rows, so
their seams start at a chain vertex and need no boundary split. Tori take the
cylinder's cut, which keeps an arc's circle parameter equal on both faces.

## What this piece is for

Piece 7 writes bodies as STEP AP242 files (ISO 10303-21 text) that other kernels read as valid solids. It is the outbound half of conversion at the border (DESIGN.md): exact geometry becomes ordinary STEP entities, and the export proves and reports the largest distance between what it wrote and the exact model.

## Number types

| What | Type | Meaning |
| --- | --- | --- |
| Geometry read | the model's definitions; `Point3`, `RVector3`, `Real` from piece 4's `curve_form` | Exact |
| Numbers written | decimal text | Rounded to `significant_digits` significant digits |
| `max_deviation`, `curve_tolerance` | `Rational` | Millimetres, exact |
| Work inside the export | `Ball` | Proven enclosures |

## The file

### Text

The file is ASCII text with `\n` line ends and no other white space outside strings. It is exactly:

```
ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('catallax export'),'2;1');
FILE_NAME('<name>.step','1970-01-01T00:00:00',(''),(''),'catallax','catallax','');
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
ENDSEC;
DATA;
<one entity per line>
ENDSEC;
END-ISO-10303-21;
```

`<name>` is `StepOptions::product_name` written as a STEP string (below). Each entity line is `#n=ENTITY(...);`, with `n` from 1 up and no gap, in increasing order, and no space anywhere outside strings.

**Strings.** A string is written between apostrophes. An apostrophe inside it is doubled, a backslash is doubled, and a character outside printable ASCII (`0x20` to `0x7E`) is written as `\X2\` followed by its UTF-16 code units as four upper-case hexadecimal digits each, then `\X0\`.

**Numbers.** A real number is written with `significant_digits` significant digits (a value below 1 is treated as 1), as a decimal within one unit in its last written place of the exact value. Its text is piece 1's `to_decimal_bounds` notation with two changes that STEP needs: a number without a decimal point or exponent gets a trailing `.` (`12000.`), and a number in scientific notation has `E` for `e` and a `.` after its first digit when it has no other (`1.E+30`, `2.5E-7`). Zero is `0.`. Which of the two nearest decimals is written is the creator's choice; it is deterministic. When the budget allows no enclosure narrow enough to prove a decimal within that unit, the export refuses with `step.budget-exhausted` (Budgets and the deadline). Integers that STEP needs as integers (a B-spline degree, knot multiplicities, the year 2014) are written as plain digits.

An enclosure that contains zero does not establish an exact zero. If decimal refinement stalls, piece 1 may prove the number coincident with zero; only that proof permits `0.` with zero rounding error. A nonzero or unresolved value keeps the existing certification refusal. Successful decimal certificates retain their existing bounds.

### The fixed entities

Entities `#1` to `#20` are always these, with `<name>` the product name, `<u>` the uncertainty (below) and `<items>` the list `#17` followed by the solids:

```
#1=APPLICATION_CONTEXT('managed model based 3d engineering');
#2=APPLICATION_PROTOCOL_DEFINITION('international standard','ap242_managed_model_based_3d_engineering',2014,#1);
#3=PRODUCT_CONTEXT('',#1,'mechanical');
#4=PRODUCT('<name>','<name>','',(#3));
#5=PRODUCT_DEFINITION_FORMATION('','',#4);
#6=PRODUCT_DEFINITION_CONTEXT('part definition',#1,'design');
#7=PRODUCT_DEFINITION('design','',#5,#6);
#8=PRODUCT_DEFINITION_SHAPE('','',#7);
#9=SHAPE_DEFINITION_REPRESENTATION(#8,#16);
#10=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));
#11=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));
#12=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());
#13=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(<u>),#10,'distance_accuracy_value','maximum deviation from the exact model');
#14=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#13))GLOBAL_UNIT_ASSIGNED_CONTEXT((#10,#11,#12))REPRESENTATION_CONTEXT('',''));
#15=PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#4));
#16=ADVANCED_BREP_SHAPE_REPRESENTATION('<name>',(<items>),#14);
#17=AXIS2_PLACEMENT_3D('',#18,#19,#20);
#18=CARTESIAN_POINT('',(0.,0.,0.));
#19=DIRECTION('',(0.,0.,1.));
#20=DIRECTION('',(1.,0.,0.));
```

`<u>` is the larger of `curve_tolerance` and `max_deviation`, written rounded up to `significant_digits` digits.
An empty body list, or bodies that are all `Empty`, give a file whose `<items>` is `#17` alone.

### Solids

For each body in `bodies` order and each of its lumps in order, one item of `<items>`:

- a lump with one shell: `MANIFOLD_SOLID_BREP('',#shell)`;
- a lump with voids: `BREP_WITH_VOIDS('',#outer,(#void1,...))`, each void an `ORIENTED_CLOSED_SHELL('',*,#shell,.F.)`.

A shell is `CLOSED_SHELL('',(#face1,...))` with its faces in shell order. A void shell is written with its faces reversed (each `ADVANCED_FACE`'s `same_sense` negated and each `FACE_BOUND`'s orientation `.F.`), so that the oriented shell's `.F.` gives back the model's faces, which point into the void.
Lumps share no faces, edges or vertices (piece 3, `JUDGE-INTERFACE.md`, ND4), so each solid's topology entities are its own. Geometry entities (curves and surfaces) are written once per definition and may be shared between solids. Points are not: each `VERTEX_POINT` has its own `CARTESIAN_POINT`, even when two vertices of touching lumps use one `PointId`.

### Topology

| Model entity | STEP entity |
| --- | --- |
| face | `ADVANCED_FACE('',(#bound1,...),#surface,<same_sense>)`, `same_sense` `.T.` when the face's sense is `Same` |
| coedge loop | `FACE_BOUND('',#loop,.T.)` with `EDGE_LOOP('',(#oe1,...))` |
| vertex loop | `FACE_BOUND('',#loop,.T.)` with `VERTEX_LOOP('',#vertex)` |
| coedge | `ORIENTED_EDGE('',*,*,#edge,<orientation>)`, `.T.` when the coedge's sense is `Same` |
| edge | `EDGE_CURVE('',#start,#end,#curve,<same_sense>)`, `.T.` when the edge's sense is `Same` |
| vertex | `VERTEX_POINT('',#point)` with `CARTESIAN_POINT('',(x,y,z))` |

No `FACE_OUTER_BOUND` is written, because the model has no outer loop. Loops, coedges and faces keep the model's orientation, which is STEP's: the face lies to the left of each bound, seen from the side its normal points to.
No `PCURVE`, `SURFACE_CURVE` or `SEAM_CURVE` is written (ND5).

**Entities STEP needs that the model lacks.** They are counted in `added_vertices` and `added_edges`, and the manifest reports them apart from the answer's topology.

- A ring edge gets one vertex, at its curve's parameter origin (piece 2's `parameter_origin`: a circle's or ellipse's `t = 0` point, an intersection curve's witness); its `EDGE_CURVE` starts and ends there.
- A loopless sphere face gets two vertices at its poles, `center - radius e3` and `center + radius e3`, and one seam edge from the first to the second along the half meridian through `center + radius e1`, on a `CIRCLE` with centre `center`, axis `-e2`, reference direction `e1` and the sphere's radius. Its bound is `EDGE_LOOP('',(#oe+,#oe-))`: the seam edge `.T.`, then `.F.`.
- A loopless torus face gets one vertex at `center + (major + minor) e1`, and two edges, each a closed curve starting and ending there: the tube circle with centre `center + major e1`, axis `-e2`, reference direction `e1` and radius `minor`; and the outer equator with centre `center`, axis `e3`, reference direction `e1` and radius `major + minor`. Its bound is `EDGE_LOOP` of the tube `.T.`, the equator `.T.`, the tube `.F.`, the equator `.F.`.

Here `e1`, `e2` and `e3` are the surface's frame (piece 2). These choices put the added entities on the seams of `JUDGE-INTERFACE.md`.

#### Decimal-collapsed planar wire gaps

An exact planar face can have two disjoint wires whose straight edges become
overlapping segments after decimal rounding. Keeping those separate wires would
write invalid topology. The exporter may reconnect them when one written segment
is strictly contained in the other, their travels oppose, and the two neighbouring
faces have parallel exact planar supports and coincident written supports.
The actual rounded line carriers, not just their endpoints, must coincide.

The gap is a constant translation between parallel exact lines. Its exact rational
L1 width bounds the Euclidean displacement over the entire removed ribbon and the
neighbouring face's change of support. That bound must not exceed either the
requested `curve_tolerance` or the unchanged 1 nm modelling tolerance. The two
wires become one and the parallel neighbouring faces become one. Split edges and
merged-face rounding errors include the displacement bound. No geometry checker,
sampling, increased tolerance or model-side merge participates in export.

`StepReport.regularizations` records each reconnection's face, edge pair,
neighbouring face pair and proven `max_deviation`. The report's topology counts
describe the file, which can now have fewer faces than the unchanged exact model.
The seam/ring `added_edges` and `added_vertices` counters retain their meanings.
Files with separated written wires follow the original writer and keep their
bytes.

#### Planar steps too narrow to read

OpenCASCADE cannot translate an edge of at most 1e-7 mm: its raw translator
rejects the shape, and its standard reader removes the edge and loses area on the
walls beside it. That floor is fixed and does not follow the file's declared
uncertainty. An exact planar face can be such a ribbon: a four-edge wire with two
straight end edges at or below the floor, written, between two longer straight
edges at least ten times their length. The exporter removes that ribbon when its
two long edges border planar walls facing the same way, so the ribbon is a step
between them, and one wall can move onto the other's plane.

A wall can move when it has one wire of straight edges, every neighbour across
those edges is planar and still contains the moved edge, and each of its vertices
has one edge outside the wall or already lies on the kept plane. A ribbon vertex
then lands on the far end of its end edge. Any other vertex with one outside edge
slides along it to the kept plane, and a vertex with more stays where it is. A
vertex that lands on a vertex the neighbour already has needs no other check
against that neighbour, whose earlier repair may already carry a deviation. Of the
walls that can move, the one with the smaller area moves, which changes the file's
volume least. All of this is exact and spends no decision work until a short-ended
ribbon exists. A shell holding a vertex loop keeps its topology.

The moved wall maps piecewise linearly, so no point of it, the ribbon or a
neighbour's trim moves further than a vertex does. The largest L1 vertex movement
is the repair's displacement bound, which must not exceed the requested
`curve_tolerance` or the 1 nm modelling tolerance. Changed faces and edges carry
that bound in their error, so `max_deviation` includes it. A larger bound, or any
failed check, keeps the original topology and bytes.

Each repair removes the ribbon and the moved wall, two faces, and the ribbon's
four edges. `StepRegularization::removed_faces` and `removed_edges` record the
count for both repair classes, so a consumer can reconcile the file with the exact
body: exact faces equal `StepReport.faces` plus every record's `removed_faces`, and
exact edges plus `added_edges` equal `StepReport.edges` plus every record's
`removed_edges`. The exact model, its names, its measured properties and its marks
do not change. The written solid differs from the exact one by at most the reported
deviation; a 1e-7 mm step on a 900 mm2 wall changes the file's volume by 9e-5 mm3.

A step wider than the floor keeps its bytes, however close. So does a whole body
scaled down until every edge is short, because its long edges fall below the floor
or its width stays large against its length.

#### Planar wedges too narrow to read

A planar face's wire can also turn back on itself across one straight edge at or
below the floor, written, between two straight edges at least ten times its
length. When those two edges border planar walls facing the same way, the face
holds a wedge closing to that short edge, and the exporter closes it. One wall
moves onto the other's plane by the rules above: its vertex on the short edge
lands on the edge's other end, so the two long edges fold onto one line. The short
edge and both long edges leave the file, the moved wall becomes part of the kept
one, and one straight edge joins the long edges' far ends unless they are one
vertex. A three-sided wedge face is then empty and leaves too. The new edge's
error is the larger of the two long edges' errors plus the displacement bound.

Each repair removes the moved wall and, for a three-sided wedge, the wedge face:
one or two faces. It removes the three edges and adds the joining edge when there
is one: two or three edges. The new edge can itself be too short to read, as when
a plate's front leans 2e-17 mm proud of a flange and the flange starts 2e-28 mm in
from the plate's side; the next repair then closes the wedge it bounds. A wedge
wider than the floor, or scaled down with its whole body, keeps its bytes as a
step does.

Other narrow shapes, such as a ribbon with more than four edges, a narrowing with
no short edge, or a curved neighbour, keep the original topology and remain
uncorrected. Other classes of export regularization are not implied by these
rules.


### Geometry

Directions are written as unit vectors: the exact direction divided by its exact length, rounded.

| Model geometry | STEP entity |
| --- | --- |
| plane | `PLANE('',#a)`, `#a` placed at `origin` with axis `normal` and reference direction `ref_dir` |
| cylinder | `CYLINDRICAL_SURFACE('',#a,radius)`, placed at `origin`, axis `axis`, reference `ref_dir` |
| cone | `CONICAL_SURFACE('',#a,0.,semi_angle)`, placed at `apex`, axis `axis`, reference `ref_dir`, with `semi_angle = atan(slope)` in radians |
| sphere | `SPHERICAL_SURFACE('',#a,radius)`, placed at `center`, axis `axis`, reference `ref_dir` |
| torus | `TOROIDAL_SURFACE('',#a,major,minor)`, placed at `center`, axis `axis`, reference `ref_dir` |
| a curve whose `curve_form` is a line | `LINE('',#p,#v)` with `#p` at `origin` and `VECTOR('',#d,1.)`, `#d` along `dir` |
| a circle form | `CIRCLE('',#a,sqrt(radius_sq))`, placed at `center`, axis `normal`, reference `ref_dir` |
| an ellipse form | `ELLIPSE('',#a,sqrt(major_sq),sqrt(minor_sq))`, placed at `center`, axis `normal`, reference `major_dir` |
| a branch form | one `B_SPLINE_CURVE_WITH_KNOTS` per edge, below |

In the retained analytic `write_step` path, a placement is `AXIS2_PLACEMENT_3D('',#location,#axis,#ref)` with a `CARTESIAN_POINT` and two `DIRECTION` entities of its own. Each `LINE`, `VECTOR`, placement and fitted branch B-spline has its own points and directions; nothing below a curve or surface entity is shared. The separate certified spline-prism exporter above deliberately shares its surface-boundary controls.
One curve entity is written per `CurveId` with a closed form, and shared by all edges on it. One surface entity is written per `SurfaceId`.
Thus analytic supports keep their analytic form (ND7), their seams stay where the model has them, and analytic curves stay analytic even when they are intersection curves with irrational data. The separate spline path exports the actual spline rather than an analytic proxy.

**Branches.** An edge on a curve with no closed form is written as one non-rational `B_SPLINE_CURVE_WITH_KNOTS('',degree,(#cp1,...),.UNSPECIFIED.,<closed>,.F.,(multiplicities),(knots),.UNSPECIFIED.)`, fitted to piece 4's `trace` of the edge's arc in the edge's direction. It is clamped (its end knots have multiplicity `degree + 1`), its first and last control points are the written points of the edge's start and end vertices (the same decimals), `<closed>` is `.T.` exactly when they are the same vertex, and the `EDGE_CURVE`'s `same_sense` is `.T.`. Degree, knots and control points are the creator's choice, within the deviation rule. Each B-spline counts in `approximated_curves`.

### Order and numbering

Entities from `#21` on are numbered in the order a depth-first walk first reaches them: start at `#16`'s items after `#17`; at each entity, visit its references in the order they appear in its own text, and number an entity when it is first reached. An entity reached again keeps its number. So the text of every entity is fixed by the model, and the same body and options give the same bytes on every machine (DESIGN.md, R10).

## The deviation

`StepReport::max_deviation` is a proven rational upper bound, in millimetres, on all of:

1. the distance between each written vertex point and the exact vertex;
2. for each edge written, the edges this piece adds included, the distance from every point of the exact arc to the written curve, and from every point of the written curve between the written vertex points (the part the edge uses) to the exact arc;
3. for each face, the distance from every point of the exact face to the written surface.

It covers decimal rounding and spline fitting. It is 0 exactly when every written number equals its exact value and no curve is approximated (for example a block with integer sizes at an axis-aligned pose).

`max_deviation` is at most `curve_tolerance`, unless the answer carries a mark `step.tolerance-not-reached` of kind `BoundKind::Distance` whose bound is `max_deviation` (ND5). This piece adds that mark, once, when:

- rounding at `significant_digits` alone makes the deviation larger than `curve_tolerance`, or
- the budget did not allow a spline fit within `curve_tolerance`; it then writes the best fit it proved.

How the bounds are proven is the creator's choice. A method that fits: bound the effect of rounding a point, a direction or a radius on each entity over the extent of the face or edge that uses it (a plane's tilt error times the face's radius about its placement point, a sphere's centre and radius errors added), and bound a spline's distance from its arc on each knot span from Taylor enclosures of the trace and exact evaluation of the polynomial span.

## Budgets and the deadline

This piece counts its `Ball` work as piece 6 does: one operation at working precision `p` costs `ceil(p / 64)^2` units, and it never uses more than the `Decider` policy's `max_precision_bits` or `max_work` on one curve or surface. Exact values it rounds (vertex positions, irrational centres and radii) are enclosed through `Decider::enclose` and `GeomStore::enclose_point`, which refuse after the deadline. So after the deadline an export refuses with `numbers.deadline` at its first such call and writes nothing, and one that needs none may finish (ARCHITECTURE.md, "The time budget").

The export refuses with `step.budget-exhausted` and writes nothing when the budget allows no proven bound on something the file needs. Two examples are `atan(slope)` for a cone under `max_work = 0`, and a spline fit to a branch whose first bound needs more work than is left. It also refuses when an enclosure is too wide to write a number within one unit in its last written place. A spline fit whose proven bound misses `curve_tolerance` is not refused: it is written with the mark `step.tolerance-not-reached` (The deviation). The export never builds a `numbers.*` refusal of its own; those come only from piece 1's calls.

## Writing

`write_step(model, bodies, options, dec, out)` builds the whole file first and writes it to `out` only when everything succeeded, so a refusal writes nothing. An error from `out` is refused with `step.io-error`, whose detail is the error's text; the bytes already written are then undefined.
The returned `StepReport` counts what the file holds: `solids` (items of `#16` other than `#17`), `faces`, `edges` (`EDGE_CURVE` entities) and `vertices` (`VERTEX_POINT` entities), each including what was added, and `added_vertices`, `added_edges`, `approximated_curves`, `max_deviation` and the per-gap `regularizations` records.

## How this piece uses the earlier pieces

- **Outcomes.** `write_step` returns `Answered<StepReport>`, whose marks hold only `step.tolerance-not-reached`. A refusal of pieces 1, 2 or 4 is returned unchanged.
- **Budgets.** As stated above; this piece never creates a `Decider` or a policy.
- **Facts.** This piece reads no facts and records none.
- **Geometry.** Piece 2's definitions and `parameter_origin`; piece 4's `curve_form` for every edge curve, and `trace` for every branch edge.

## Refusal codes

| Code | When |
| --- | --- |
| `step.not-solid` | a body's kind is neither `Solid` nor `Empty` (piece 3 has no such kind yet, so no test reaches it) |
| `step.unsupported-geometry` | a surface or curve form this piece cannot write (none today; the enums are non-exhaustive) |
| `step.io-error` | writing to `out` failed |
| `step.budget-exhausted` | the budget allowed no proven bound on a value or spline fit the file needs, or no enclosure narrow enough to write a number within one unit in its last written place |

Refusals of pieces 1, 2 and 4 pass through unchanged.

## Programming errors

These may panic. None can come from a question file:

- a `BodyId` from another model, or out of range;
- a body that fails `Model::check`.

## Exact public interface

The crate re-exports these types at its root. Everything listed is `pub`. The compiled crate in `crates/catallax-step` is authoritative if a signature listing becomes stale; update this contract with any API change.

```rust
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct StepOptions { pub product_name: String, pub curve_tolerance: Rational, pub significant_digits: u32 }
// StepOptions: Default ("catallax", 1/10^7 mm, 17 digits)
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct StepReport {
    pub max_deviation: Rational,
    pub approximated_curves: usize,
    pub added_vertices: usize,
    pub added_edges: usize,
    pub solids: usize,
    pub faces: usize,
    pub edges: usize,
    pub vertices: usize,
}
pub fn write_step<W: std::io::Write>(model: &Model, bodies: &[BodyId], options: &StepOptions, dec: &mut Decider, out: &mut W) -> Answered<StepReport>;
```

## What defines correct behaviour

The test writer builds tests from the cases below, parsing the file with a small test-only reader of the syntax above (or a third-party STEP parser used only as a dev-dependency). Every expected answer follows from this description and maths.

### The text

1. `block([2, 3, 5], identity)` with default options: the file starts with the header above with `catallax` as the name and `<u> = 1.E-7`, holds entities `#1` to `#20` exactly as listed, and has these entity counts in all (the fixed entities included): `MANIFOLD_SOLID_BREP` 1, `CLOSED_SHELL` 1, `ADVANCED_FACE` 6, `FACE_BOUND` 6, `EDGE_LOOP` 6, `ORIENTED_EDGE` 24, `EDGE_CURVE` 12, `VERTEX_POINT` 8, `CARTESIAN_POINT` 27, `LINE` 12, `VECTOR` 12, `DIRECTION` 26, `AXIS2_PLACEMENT_3D` 7, `PLANE` 6; 170 entities. The report: `solids 1`, `faces 6`, `edges 12`, `vertices 8`, nothing added or approximated, `max_deviation 0`, no mark. The eight vertex points are the corners, written `(0.,0.,0.)`, `(2.,0.,0.)` and so on.
2. Numbering: every reference `#n` in the file names an entity that exists, every entity from `#21` on is reached by the walk, and renumbering by the walk rule reproduces the file.
3. `cylinder(2, 5, identity)`: `ADVANCED_FACE` 3, `FACE_BOUND` 4, `EDGE_LOOP` 4, `ORIENTED_EDGE` 4, `EDGE_CURVE` 2, `VERTEX_POINT` 2 at `(2.,0.,0.)` and `(2.,0.,5.)` (the seam), `CIRCLE` 2 of radius `2.`, `PLANE` 2, `CYLINDRICAL_SURFACE` 1; `added_vertices 2`, `added_edges 0`, `max_deviation 0`.
4. `sphere(5, identity)`: `ADVANCED_FACE` 1, `EDGE_LOOP` 1 with two `ORIENTED_EDGE`s of one `EDGE_CURVE` (`.T.` then `.F.`), `VERTEX_POINT` 2 at `(0.,0.,-5.)` and `(0.,0.,5.)`, one `CIRCLE` through `(5.,0.,0.)`, `SPHERICAL_SURFACE` 1 of radius `5.`; `added_vertices 2`, `added_edges 1`.
5. `torus(3, 1, identity)`: one `VERTEX_POINT` at `(4.,0.,0.)`, two `EDGE_CURVE`s on `CIRCLE`s of radius `1.` (centre `(3.,0.,0.)`) and `4.` (centre `0`), one `EDGE_LOOP` of four `ORIENTED_EDGE`s, `TOROIDAL_SURFACE` 1 with `3.` and `1.`; `added_vertices 1`, `added_edges 2`.
6. `cone(3, 0, 4, identity)`: no `VERTEX_LOOP`; the cone face's one `EDGE_LOOP` is the base ring, then a seam edge from the apex vertex `(0.,0.,4.)` to the ring vertex `(3.,0.,0.)` `.F.`, then the same seam `.T.`. A `CONICAL_SURFACE` placed at `(0.,0.,4.)` with axis `(0.,0.,-1.)`, radius `0.` and `semi_angle` within one unit in the last place of `atan(3/4) = 0.6435011087932843868...`; `added_vertices 1` (the base ring), `added_edges 1` (the seam). Its `max_deviation` is positive and at most `curve_tolerance`, since the angle is rounded.
7. Strings: the product name `a'b\c` is written `'a''b\\c'`; `é` is written `'\X2\00E9\X0\'`.
8. `empty(model)` alone, and an empty `bodies` list: the twenty fixed entities only, with `<items>` `#17`.

### Several solids and voids

9. A body of two lumps (piece 3's case 15, two blocks): two `MANIFOLD_SOLID_BREP`s, 12 faces, 24 edges, 16 vertices.
10. Piece 3's case 17, a block with a void: one `BREP_WITH_VOIDS` with one `ORIENTED_CLOSED_SHELL` (`.F.`) whose `CLOSED_SHELL`'s six faces have `same_sense .T.` and `FACE_BOUND`s with `.F.`.
11. Two bodies passed together: their lumps in order, one file.

### Positions and the deviation

12. `block([1/3, 1/3, 1/3], identity)` with 5 digits: the corner `(1/3, 1/3, 1/3)` is written `(0.33333,0.33333,0.33333)` or with `0.33334`; `max_deviation` is at least `sqrt(3) / 300000` (about `5.774e-6`, the distance between that corner and `(0.33333, 0.33333, 0.33333)` when that is what was written, larger otherwise), and the answer carries `step.tolerance-not-reached` with that bound. With 17 digits: `max_deviation` at most `10^-7` and no mark.
13. Posed: `block([2, 3, 5], pose)` with rotation `(2, 1, 0, 0)` (matrix with entries `0, ±3/5, ±4/5, 1`) and translation `(10^6, -7, 1/3)`: every written point is within one unit in its last place of the exact posed corner, every direction within the same of the exact unit vector, and `max_deviation` is at least the largest true distance between a written and an exact vertex.
14. Scale: `block([2, 3, 5] * 10^30, identity)` with 17 digits: the sizes are exact decimals, so `max_deviation` is 0. `block([10^30 / 3, 1, 1], identity)`: rounding at 17 digits moves the corner by about `3.3 * 10^12`, so `max_deviation` is at least that and the mark is present.
15. The written geometry holds the exact geometry within `max_deviation`: for every case here, the test parses each `CARTESIAN_POINT`, `DIRECTION` and radius, and checks that each exact vertex is within `max_deviation` of its written point (exact rationals), and that each exact vertex is within `max_deviation` of the written surface of every face that uses it (the distance to a written plane, cylinder or sphere in exact arithmetic, or in a proven enclosure).

### Intersection curves

These need piece 4 to make the curves.

16. A body with edges on intersection curves with a line form or a circle form (piece 5's case 13 or the irrational case of 30, or a hand-built body on piece 4's case 6 lines): those edges are written as `LINE` and `CIRCLE`, never as B-splines, and `approximated_curves` is 0.
17. A body with branch edges (piece 5's case 23 `A ∩ B`, or a hand-built body on piece 4's case 25 loop): each branch edge is one `B_SPLINE_CURVE_WITH_KNOTS`, clamped, with its first and last control points equal (as text) to its vertices' points; `approximated_curves` counts them. The test evaluates each B-spline exactly from its written decimals at 1,000 parameters per span and checks that every point is within `max_deviation` of both of the branch's surfaces, by the exact distance to each surface enclosed with `Ball` arithmetic (piece 2's distance conversions are upper bounds that can be up to about twice the distance, so they may be used only against `2 max_deviation`), and that `max_deviation` is at most `10^-7`. With `curve_tolerance` `10^-12`: the same check at `10^-12`, or the mark with the bound reached.
18. Under `BudgetPolicy { max_precision_bits: 64, max_work: 10_000 }` (this budget admits piece 1's `f64` filter), case 17 either meets the tolerance, carries the mark with its proven bound, or is refused by piece 1, 2 or 4 with that piece's code unchanged or with `step.budget-exhausted`, writing nothing to `out`. Piece 4's `trace` may refuse at this budget, for example with `numbers.budget-exhausted`, since it refuses whenever a decision locating the arc's ends ends `NotExact` or piece 1 refuses. When the export answers, the bound always holds. Under `BudgetPolicy { max_precision_bits: 64, max_work: 0 }` every export of cases 1 to 17 gives the same three outcomes: an answer within the tolerance, an answer with the mark whose bound holds, or one of those refusals.

### Determinism

19. Exporting the same body twice gives the same bytes. Exporting it on eight threads, from eight clones of the model, gives the same bytes.
20. A `Decider` whose deadline has passed: each export either refuses with `numbers.deadline`, writing nothing to `out`, or writes the same bytes as with no deadline, never others (ARCHITECTURE.md, "The time budget").

### Reading back

Tests in this repository check the STEP structure above. Interoperability with other CAD readers needs separate measurement. A rejected form is an exchange-contract question, not permission for a silent local fallback.

### Volume and running time of the tests

Cases 13, 14 and 15 are the generated families: run them on at least 1,000 random rational poses and sizes per test run, from fixed seeds. Keep the default test run under about 5 minutes in release mode on one core.

## Values checked for this description

The architect counted the entities of cases 1, 3, 4, 5 and 6 from the mapping above, and checked `atan(3/4) = 0.643501108793284386802809228717...`, `sqrt(3) / 300000 = 5.7735026918962576...e-6` and `10^30 / 3` rounded at 17 digits with `sympy`.

## Speed: measured, never gated

No test fails for speed. A benchmark target measures and prints:

1. time and file size against face count, for the primitives and for the results of piece 5's first-milestone cases once they exist;
2. B-spline control-point counts and fitting time against `curve_tolerance` (`10^-6` to `10^-12`) for the branches of piece 4's cases 17 to 26.

The numbers are written to the benchmark output so later runs can be compared.

## Out of scope for this piece

- STEP import (conversion on the way in comes later).
- Assemblies, colours, product metadata beyond one part.
- Free-form surfaces, which the model does not have yet.
- Any tolerance inside the kernel. The deviation here is a reported bound at the border.
