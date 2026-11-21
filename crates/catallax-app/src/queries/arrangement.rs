//! Clash or clear for two solids whose overlap Boolean refused.
//!
//! OpenCASCADE's `BRepExtrema_DistShapeShape::SolidTreatment` classifies the
//! vertices of one shape against the other solid and reports an inner
//! solution once one lies inside. SMLib's `SmApiBrepRelationship`
//! (`SM_API/src/SmApiQueries.cpp`) first takes the minimum distance between
//! the two boundaries. When it is positive, every point of one boundary has
//! one class against the other solid, so one classified point tells
//! separate solids from one inside the other. When the boundaries meet,
//! vertices, one point per face and edge pieces strictly inside the other
//! solid show that the solids overlap.
//!
//! Every step here is proven, and the steps run in that order. Vertices are
//! classified first, then one point strictly inside each face; a sample
//! outside the other solid's proven box is outside without a ray. A
//! boundary shell with one sample strictly inside the other solid and one
//! strictly outside crosses its boundary. Otherwise every pair of faces
//! whose boxes meet is proven apart by the certified interval search, and
//! the samples decide per shell, as SMLib decides per solid, since a solid
//! may have voids: with every shell outside the other solid the pair is
//! clear, and its gap comes from [`catallax_measure::solid_gap`], whose
//! noncrossing precondition now holds; with every shell of one solid inside
//! the other and none of the other's inside it, the overlap is the inner
//! solid itself. A point strictly inside both solids proves a clash on its
//! own. When every shell is sampled outside the other solid, no shell lies
//! wholly inside it, so no solid's volume can be the overlap: points of a
//! dyadic grid over the two proven boxes' overlap are then classified first,
//! coarsest level first, and one inside both settles the clash before any
//! faces are proven apart.
//!
//! A clash proven from samples needs a positive lower end for its overlap's
//! enclosure: a ball inside both solids. Take p on the boundary of X
//! strictly inside Y, with a proven lower bound r of its distance to Y's
//! boundary; the ball of radius r about p lies in Y. A step of r/2 along
//! X's inward normal at p gives q. Once q is classified inside X with a
//! proven lower bound rho of its distance to X's boundary, the ball of
//! radius min(rho, r/2) about q lies in both solids, so the overlap is at
//! least 4/3 pi min(rho, r/2)^3, more than 4 min(rho, r/2)^3. A grid point
//! inside both solids grows its ball the same way. The two proven boxes'
//! overlap bounds the overlap above. Such an enclosure is the bounded
//! answer kind and carries the `measure.accuracy-not-reached` mark when it
//! is wider than requested.

use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, btree_map::Entry};

use catallax_boolean::{Containment, classify_point, face_point};
use catallax_measure::{Accuracy, Enclosure, distance::PreparedFace};
use catallax_numbers::{
    Answer, Answered, Ball, Decider, MarkLog, Outcome, Rational, Real, Refusal, Sign,
};
use catallax_shapes::{
    Point3, RationalPatch, Side as Carrier, SplineSupport, Surface, Vector3, distance::Options,
};
use catallax_topology::{BodyId, FaceId, Loop, Model, Sense, VertexId};

use super::{BodyRelation, box_accuracy, gap_relation, mark_width};

/// Bits of the outward enclosures that compare sample coordinates with
/// boxes, as `PreparedFace` encloses face boxes and boundary points.
const BITS: u32 = 128;

/// Dyadic levels of the grid probed over the two boxes' overlap: at most
/// 73 points. Source: the three levels of the parallel planar witness search
/// of `catallax_measure::distance::PreparedFace` (`PLANAR_LEVELS`).
const GRID_LEVELS: u32 = 3;

/// Exact lower and upper corners.
type Corners = [[Rational; 2]; 3];

/// Where a boundary sample lies against the other solid.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Class {
    Inside,
    Outside,
    Boundary,
    Unknown,
}

/// A point of one solid's boundary and the faces of that solid it lies on.
struct Sample {
    point: Point3,
    faces: Vec<FaceId>,
}

/// One boundary shell and the classes its samples have had so far.
struct Shell {
    faces: Vec<FaceId>,
    inside: bool,
    outside: bool,
}

/// One operand: its proven box, its faces with their proven boxes, and its
/// boundary shells.
struct Side {
    body: BodyId,
    corners: Option<Corners>,
    faces: Vec<(FaceId, Option<Corners>)>,
    shells: Vec<Shell>,
}

fn box_corners(b: &catallax_measure::BoundingBox) -> Corners {
    std::array::from_fn(|i| [b.min[i].clone(), b.max[i].clone()])
}

impl Side {
    fn new(model: &Model, body: BodyId, dec: &mut Decider) -> Result<Self, Refusal> {
        let accuracy = box_accuracy();
        let corners = catallax_measure::bounding_box(model, body, &accuracy, dec)
            .ok()
            .filter(|b| b.marks.is_empty())
            .and_then(|b| b.value)
            .map(|b| box_corners(&b));
        dec.check_deadline()?;
        let mut faces = Vec::new();
        for face in model.faces_of_body(body) {
            let bounds = catallax_measure::face_bounding_box(model, face, &accuracy, dec)
                .ok()
                .filter(|b| b.marks.is_empty())
                .map(|b| box_corners(&b.value));
            dec.check_deadline()?;
            faces.push((face, bounds));
        }
        let mut shells = Vec::new();
        for &lump in &model.body(body).lumps {
            for &shell in &model.lump(lump).shells {
                shells.push(Shell {
                    faces: model.shell(shell).faces.clone(),
                    inside: false,
                    outside: false,
                });
            }
        }
        Ok(Self {
            body,
            corners,
            faces,
            shells,
        })
    }
}

/// Distinct vertices of some faces in order, with the faces each lies on.
fn shell_vertices(model: &Model, faces: &[FaceId]) -> Vec<(VertexId, Vec<FaceId>)> {
    let mut out: Vec<(VertexId, Vec<FaceId>)> = Vec::new();
    let mut add = |vertex: VertexId, face: FaceId| match out.iter_mut().find(|v| v.0 == vertex) {
        Some(v) if !v.1.contains(&face) => v.1.push(face),
        Some(_) => {}
        None => out.push((vertex, vec![face])),
    };
    for &face in faces {
        for lp in &model.face(face).loops {
            match model.loop_(*lp) {
                Loop::Vertex(v) => add(*v, face),
                Loop::Coedges(coedges) => {
                    for c in coedges {
                        let edge = model.edge(model.coedge(*c).edge);
                        for v in [edge.start, edge.end].into_iter().flatten() {
                            add(v, face);
                        }
                    }
                }
            }
        }
    }
    out
}

fn enclose(point: &Point3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_budgeted_at(&point.x, BITS)?,
        dec.enclose_budgeted_at(&point.y, BITS)?,
        dec.enclose_budgeted_at(&point.z, BITS)?,
    ])
}

/// The largest coordinate gap between a point's enclosure and a box: a
/// lower bound of their Euclidean distance.
fn point_gap(c: &Corners, p: &[Ball; 3]) -> Rational {
    (0..3)
        .map(|i| (&c[i][0] - &p[i].hi()).max(&p[i].lo() - &c[i][1]))
        .fold(Rational::zero(), Rational::max)
}

/// The largest coordinate gap between two boxes.
fn box_gap(a: &Corners, b: &Corners) -> Rational {
    (0..3)
        .map(|i| (&a[i][0] - &b[i][1]).max(&b[i][0] - &a[i][1]))
        .fold(Rational::zero(), Rational::max)
}

/// Certified distance searches and chart control boxes, one per face.
#[derive(Default)]
struct Faces {
    prepared: BTreeMap<FaceId, Option<PreparedFace>>,
    charts: BTreeMap<FaceId, Option<Vec<Corners>>>,
}

impl Faces {
    fn prepare(&mut self, model: &Model, face: FaceId, dec: &mut Decider) -> Result<(), Refusal> {
        if let Entry::Vacant(slot) = self.prepared.entry(face) {
            let prepared = PreparedFace::new(model, face, dec).ok();
            dec.check_deadline()?;
            slot.insert(prepared);
        }
        Ok(())
    }

    /// The control boxes of a rational spline face's charts. A rational
    /// patch with positive weights lies in the convex hull of its control
    /// points (Piegl and Tiller, The NURBS Book, section 4.4), so the union
    /// of the boxes holds the face; an atlas face is the union of its charts.
    fn charts(
        &mut self,
        model: &Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<&Vec<Corners>>, Refusal> {
        if let Entry::Vacant(slot) = self.charts.entry(face) {
            let patches = match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::RationalPatch(patch)) => vec![patch.clone()],
                Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => {
                    atlas.charts().iter().map(|c| c.patch().clone()).collect()
                }
                _ => Vec::new(),
            };
            let mut boxes = Some(Vec::new());
            for patch in &patches {
                let carrier = patch.carrier();
                let mut positive = true;
                for w in carrier.weights() {
                    positive &= dec.enclose_budgeted_at(w, BITS)?.lo().sign() == Sign::Positive;
                }
                let mut hull: Option<Corners> = None;
                for c in carrier.controls() {
                    let p = enclose(c, dec)?;
                    let b: Corners = std::array::from_fn(|i| [p[i].lo(), p[i].hi()]);
                    hull = Some(match hull {
                        None => b,
                        Some(h) => std::array::from_fn(|i| {
                            [
                                h[i][0].clone().min(b[i][0].clone()),
                                h[i][1].clone().max(b[i][1].clone()),
                            ]
                        }),
                    });
                }
                match (positive, hull, boxes.as_mut()) {
                    (true, Some(h), Some(list)) => list.push(h),
                    _ => boxes = None,
                }
            }
            slot.insert(boxes.filter(|list| !list.is_empty()));
        }
        Ok(self.charts[&face].as_ref())
    }

    /// A proven lower bound of the distance from a point, enclosed by `p`,
    /// to one face: its whole analytic carrier's ([`carrier_floor`]), or the
    /// least over the patches of a rational spline face's charts
    /// ([`patch_floor`]). A ball needs a positive radius, not the distance to
    /// the certified search's default width.
    fn floor(
        &mut self,
        model: &Model,
        face: FaceId,
        p: &[Ball; 3],
        ceiling: Option<&Rational>,
        dec: &mut Decider,
    ) -> Result<Rational, Refusal> {
        let surface = model.geom().surface(model.face(face).surface);
        if let Some(bound) = carrier_floor(surface, p, dec)? {
            return Ok(bound);
        }
        let patches = match surface {
            Surface::BSpline(SplineSupport::RationalPatch(patch)) => vec![patch.clone()],
            Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => {
                atlas.charts().iter().map(|c| c.patch().clone()).collect()
            }
            _ => return Ok(Rational::zero()),
        };
        let mut least: Option<Rational> = None;
        for patch in &patches {
            let bound = patch_floor(patch, p, ceiling, dec)?;
            least = Some(least.map_or(bound.clone(), |l| l.min(bound)));
        }
        Ok(least.unwrap_or_else(Rational::zero))
    }

    /// The least gap between the chart control boxes of one face and those
    /// of another, or its proven box, when every pair is apart.
    fn charts_apart(
        &mut self,
        model: &Model,
        (f, fb): (FaceId, &Option<Corners>),
        (g, gb): (FaceId, &Option<Corners>),
        dec: &mut Decider,
    ) -> Result<Option<Rational>, Refusal> {
        let a = self
            .charts(model, f, dec)?
            .cloned()
            .or_else(|| fb.clone().map(|b| vec![b]));
        let b = self
            .charts(model, g, dec)?
            .cloned()
            .or_else(|| gb.clone().map(|b| vec![b]));
        let (Some(a), Some(b)) = (a, b) else {
            return Ok(None);
        };
        let least = a
            .iter()
            .flat_map(|x| b.iter().map(|y| box_gap(x, y)))
            .min()
            .expect("nonempty boxes");
        Ok((least.sign() == Sign::Positive).then_some(least))
    }
}

/// The lower end of an enclosure's absolute value.
fn abs_lo(b: &Ball) -> Rational {
    if b.lo().sign() == Sign::Positive {
        b.lo()
    } else if b.hi().sign() == Sign::Negative {
        -b.hi()
    } else {
        Rational::zero()
    }
}

/// A proven lower bound of the distance from a point, enclosed by `p`, to an
/// analytic face's whole carrier, from the closed forms of OpenCASCADE's
/// `Extrema_ExtPElS` for a point and a plane, cylinder, cone, sphere or
/// torus, evaluated in ball arithmetic. A torus, or a sheet of one, is the
/// set at its minor radius from its core circle, and the distance to that
/// circle is 1-Lipschitz, so the distance to the circle less the minor radius
/// bounds the distance to the torus. A cone's nappe lies on the line through
/// its apex in each meridian half-plane. `None` for other carriers.
fn carrier_floor(
    surface: &Surface,
    p: &[Ball; 3],
    dec: &mut Decider,
) -> Result<Option<Rational>, Refusal> {
    fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
        &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
    }
    fn minus(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
        std::array::from_fn(|i| &a[i] - &b[i])
    }
    // Axial coordinate and distance from the axis through the origin along `a`.
    fn axial(w: &[Ball; 3], a: &[Ball; 3]) -> Option<(Ball, Ball)> {
        let h = dot(w, a).checked_div(&dot(a, a).sqrt()?)?;
        let rho = (&dot(w, w) - &h.square()).sqrt()?;
        Some((h, rho))
    }
    fn vector(v: [&Real; 3], dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        Ok([
            dec.enclose_budgeted_at(v[0], BITS)?,
            dec.enclose_budgeted_at(v[1], BITS)?,
            dec.enclose_budgeted_at(v[2], BITS)?,
        ])
    }
    fn point(v: &Point3) -> [&Real; 3] {
        [&v.x, &v.y, &v.z]
    }
    fn direction(v: &Vector3) -> [&Real; 3] {
        [&v.x, &v.y, &v.z]
    }
    let value = match surface {
        Surface::Plane(s) => {
            let n = vector(direction(&s.normal), dec)?;
            let w = minus(p, &vector(point(&s.origin), dec)?);
            dot(&n, &n)
                .sqrt()
                .and_then(|length| dot(&w, &n).checked_div(&length))
        }
        Surface::Cylinder(s) => {
            let w = minus(p, &vector(point(&s.origin), dec)?);
            let a = vector(direction(&s.axis), dec)?;
            let radius = dec.enclose_budgeted_at(&s.radius_sq, BITS)?.sqrt();
            axial(&w, &a).zip(radius).map(|((_, rho), r)| &rho - &r)
        }
        Surface::Sphere(s) => {
            let w = minus(p, &vector(point(&s.center), dec)?);
            let radius = dec.enclose_budgeted_at(&s.radius, BITS)?;
            dot(&w, &w).sqrt().map(|d| &d - &radius)
        }
        Surface::Torus(catallax_shapes::Torus {
            center,
            axis,
            major,
            minor,
            ..
        })
        | Surface::TorusSheet(catallax_shapes::TorusSheet {
            center,
            axis,
            major,
            minor,
            ..
        }) => {
            let w = minus(p, &vector(point(center), dec)?);
            let a = vector(direction(axis), dec)?;
            let [major, minor] = [
                dec.enclose_budgeted_at(major, BITS)?,
                dec.enclose_budgeted_at(minor, BITS)?,
            ];
            axial(&w, &a)
                .and_then(|(h, rho)| (&h.square() + &(&rho - &major).square()).sqrt())
                .map(|core| &core - &minor)
        }
        Surface::Cone(s) => {
            let w = minus(p, &vector(point(&s.apex), dec)?);
            let a = vector(direction(&s.axis), dec)?;
            let slope = dec.enclose_budgeted_at(&s.slope, BITS)?;
            let one = Ball::from_rational(&Rational::from_i64(1), BITS);
            axial(&w, &a).and_then(|(h, rho)| {
                (&rho - &(&slope * &h)).checked_div(&(&one + &slope.square()).sqrt()?)
            })
        }
        _ => return Ok(None),
    };
    Ok(Some(value.as_ref().map_or_else(Rational::zero, abs_lo)))
}

/// A proven lower bound of the distance from a point, enclosed by `p`, to a
/// rational patch, capped at `ceiling`: best-first subdivision of its
/// parameter domain into quarters, each cell's image enclosed by
/// `BSplineSurface::evaluate_box`. Once the nearest remaining cell's box is
/// apart from the point, its gap bounds every cell. Cells at least `ceiling`
/// away are dropped. Zero when the certified search's default cell count
/// runs out first.
fn patch_floor(
    patch: &RationalPatch,
    p: &[Ball; 3],
    ceiling: Option<&Rational>,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    type Cell = [[Rational; 2]; 2];
    let carrier = patch.carrier();
    let gap = |cell: &Cell, dec: &mut Decider| -> Result<Rational, Refusal> {
        let image = carrier.evaluate_box(cell[0].each_ref(), cell[1].each_ref(), BITS, dec)?;
        Ok(point_gap(
            &std::array::from_fn(|i| [image[i].lo(), image[i].hi()]),
            p,
        ))
    };
    let mut queue: BinaryHeap<Reverse<(Rational, usize, Cell)>> = BinaryHeap::new();
    let root = patch.domain().clone();
    let root_gap = gap(&root, dec)?;
    if ceiling.is_some_and(|c| root_gap >= *c) {
        return Ok(ceiling.cloned().expect("a ceiling"));
    }
    queue.push(Reverse((root_gap, 0, root)));
    let mut cells = 1;
    while let Some(Reverse((least, _, cell))) = queue.pop() {
        dec.check_deadline()?;
        if least.sign() == Sign::Positive {
            return Ok(least);
        }
        if cells + 4 > Options::default().max_cells {
            return Ok(Rational::zero());
        }
        let mid = [0, 1].map(|k| (&cell[k][0] + &cell[k][1]).mul_pow2(-1));
        for (i, j) in [(0, 0), (0, 1), (1, 0), (1, 1)] {
            let u = if i == 0 {
                [cell[0][0].clone(), mid[0].clone()]
            } else {
                [mid[0].clone(), cell[0][1].clone()]
            };
            let v = if j == 0 {
                [cell[1][0].clone(), mid[1].clone()]
            } else {
                [mid[1].clone(), cell[1][1].clone()]
            };
            let child = [u, v];
            let g = gap(&child, dec)?;
            cells += 1;
            if ceiling.is_none_or(|c| g < *c) {
                queue.push(Reverse((g, cells, child)));
            }
        }
    }
    Ok(ceiling.cloned().unwrap_or_else(Rational::zero))
}

/// A proven lower bound of the distance from `point` to the side's boundary,
/// capped at `ceiling`. Faces are visited in order of their boxes' distance
/// from the point, and the rest are skipped once that distance reaches the
/// bound found, as `solid_gap` prunes face pairs.
fn clearance(
    model: &Model,
    faces: &mut Faces,
    side: &Side,
    point: &Point3,
    ceiling: Option<&Rational>,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let p = enclose(point, dec)?;
    let mut order: Vec<(Rational, FaceId)> = side
        .faces
        .iter()
        .map(|(face, bounds)| {
            let gap = bounds
                .as_ref()
                .map_or_else(Rational::zero, |b| point_gap(b, &p));
            (gap, *face)
        })
        .collect();
    order.sort();
    let mut best = ceiling.cloned();
    for (gap, face) in order {
        dec.check_deadline()?;
        if best.as_ref().is_some_and(|b| gap >= *b) {
            break;
        }
        let floor = faces.floor(model, face, &p, best.as_ref(), dec)?.max(gap);
        best = Some(best.map_or(floor.clone(), |b| b.min(floor)));
    }
    Ok(best.unwrap_or_else(Rational::zero))
}

fn classify(
    scratch: &mut Model,
    side: &Side,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Class, Refusal> {
    if let Some(c) = &side.corners
        && point_gap(c, &enclose(point, dec)?).sign() == Sign::Positive
    {
        return Ok(Class::Outside);
    }
    // A solid bounded by one whole closed analytic surface, a sphere or a
    // torus, is one side of it: the exact sign of its equation, with no ray.
    if let [(face, _)] = side.faces.as_slice()
        && scratch.face(*face).loops.is_empty()
        && !matches!(
            scratch.geom().surface(scratch.face(*face).surface),
            Surface::BSpline(_)
        )
    {
        let surface = scratch.face(*face).surface;
        // The face normal points out of the material.
        let outward = match scratch.face(*face).sense {
            Sense::Same => Carrier::Front,
            Sense::Reversed => Carrier::Back,
        };
        match scratch.geom().side_of(surface, point, dec) {
            Outcome::Proven(s) if s == outward => return Ok(Class::Outside),
            Outcome::Proven(_) => return Ok(Class::Inside),
            Outcome::ProvenCoincident(_) => return Ok(Class::Boundary),
            Outcome::NotExact(_) | Outcome::Refused(_) => dec.check_deadline()?,
        }
    }
    Ok(match classify_point(scratch, side.body, point, dec) {
        Outcome::Proven(Containment::Inside) => Class::Inside,
        Outcome::Proven(Containment::Outside) => Class::Outside,
        Outcome::ProvenCoincident(_) => Class::Boundary,
        Outcome::NotExact(_) | Outcome::Refused(_) => {
            dec.check_deadline()?;
            Class::Unknown
        }
    })
}

/// The search state shared by the steps.
struct Arrangement<'m> {
    model: &'m Model,
    scratch: Model,
    sides: [Side; 2],
    faces: Faces,
    /// Boundary points of `sides[i]` strictly inside the other solid.
    witnesses: Vec<(usize, Sample)>,
    /// Some boundary point lies on the other solid's boundary.
    contact: bool,
}

impl Arrangement<'_> {
    /// Classifies each shell's vertices against the other solid until one
    /// shell has vertices on both sides, which proves the boundaries cross.
    fn vertices(&mut self, dec: &mut Decider) -> Result<bool, Refusal> {
        for x in 0..2 {
            for s in 0..self.sides[x].shells.len() {
                let samples = shell_vertices(self.model, &self.sides[x].shells[s].faces);
                for (vertex, faces) in samples {
                    dec.check_deadline()?;
                    let point = self.model.geom().position(self.model.vertex(vertex).point);
                    let class = classify(&mut self.scratch, &self.sides[1 - x], &point, dec)?;
                    if self.record(x, s, class, Sample { point, faces }) {
                        return Ok(true);
                    }
                }
            }
        }
        Ok(false)
    }

    /// Classifies one point strictly inside each face against the other
    /// solid, as SMLib samples faces, until a shell has crossed.
    fn face_points(&mut self, dec: &mut Decider) -> Result<bool, Refusal> {
        for x in 0..2 {
            for s in 0..self.sides[x].shells.len() {
                for face in self.sides[x].shells[s].faces.clone() {
                    let Ok(point) = face_point(&mut self.scratch, face, dec) else {
                        dec.check_deadline()?;
                        continue;
                    };
                    let point = point.value;
                    let class = classify(&mut self.scratch, &self.sides[1 - x], &point, dec)?;
                    let faces = vec![face];
                    if self.record(x, s, class, Sample { point, faces }) {
                        return Ok(true);
                    }
                }
            }
        }
        Ok(false)
    }

    /// Whether every shell of both solids has a sample outside the other
    /// solid and none inside it. Boundaries that do not cross then leave the
    /// two solids separate.
    fn outside(&self) -> bool {
        self.sides
            .iter()
            .all(|side| side.shells.iter().all(|s| s.outside && !s.inside))
    }

    /// Records a sample's class; true once its shell has crossed.
    fn record(&mut self, x: usize, s: usize, class: Class, sample: Sample) -> bool {
        let shell = &mut self.sides[x].shells[s];
        match class {
            Class::Inside => {
                shell.inside = true;
                self.witnesses.push((x, sample));
            }
            Class::Outside => shell.outside = true,
            Class::Boundary => self.contact = true,
            Class::Unknown => {}
        }
        shell.inside && shell.outside
    }

    /// Proves that no face of one solid meets a face of the other, or names
    /// the first pair it cannot. Pairs whose proven boxes, or the control
    /// boxes of their rational charts, are apart need nothing more. The rest
    /// go to the certified interval search
    /// ([`PreparedFace::face_to_face_below`]) with the search's default width
    /// as ceiling, all sharing its default cell count, as `solid_gap` settles
    /// its open pairs. A pair is apart once the proven lower bound is positive.
    /// Returns the least proven lower bound over all pairs, a lower bound of
    /// the distance between the two boundaries.
    fn apart(&mut self, dec: &mut Decider) -> Result<Result<Rational, String>, Refusal> {
        let options = Options::default();
        let mut cells = options.max_cells;
        let mut least: Option<Rational> = None;
        let mut bound = |lower: Rational| {
            least = Some(least.take().map_or(lower.clone(), |l| l.min(lower)));
        };
        let [a, b] = &self.sides;
        for (f, fb) in &a.faces {
            for (g, gb) in &b.faces {
                dec.check_deadline()?;
                if let (Some(x), Some(y)) = (fb, gb)
                    && box_gap(x, y).sign() == Sign::Positive
                {
                    bound(box_gap(x, y));
                    continue;
                }
                if let Some(gap) = self
                    .faces
                    .charts_apart(self.model, (*f, fb), (*g, gb), dec)?
                {
                    bound(gap);
                    continue;
                }
                self.faces.prepare(self.model, *f, dec)?;
                self.faces.prepare(self.model, *g, dec)?;
                let (Some(Some(x)), Some(Some(y))) =
                    (self.faces.prepared.get(f), self.faces.prepared.get(g))
                else {
                    return Ok(Err(format!(
                        "{f:?} and {g:?} have no certified distance search"
                    )));
                };
                if cells == 0 {
                    return Ok(Err(format!(
                        "{f:?} and {g:?} are not proven apart within the cell count"
                    )));
                }
                let bounds =
                    match x.face_to_face_below(y, &mut self.scratch, &options.width, cells, dec) {
                        Ok(bounds) => bounds,
                        Err(refusal) => {
                            dec.check_deadline()?;
                            return Ok(Err(format!(
                                "{f:?} and {g:?} are not proven apart: {}",
                                refusal.code
                            )));
                        }
                    };
                cells = cells.saturating_sub(bounds.cells);
                if bounds.lower.sign() != Sign::Positive {
                    return Ok(Err(format!("{f:?} and {g:?} are not proven apart")));
                }
                bound(bounds.lower);
            }
        }
        Ok(Ok(least.expect("faces on both sides")))
    }

    /// With the boundaries apart, each shell lies wholly inside or wholly
    /// outside the other solid, so any one sample decides it. Returns, per
    /// side, whether each shell is inside, or `None` when a shell has no
    /// sample strictly inside or outside.
    fn classes(&self) -> Option<[Vec<bool>; 2]> {
        let decided = |side: &Side| -> Option<Vec<bool>> {
            side.shells
                .iter()
                .map(|s| (s.inside || s.outside).then_some(s.inside))
                .collect()
        };
        Some([decided(&self.sides[0])?, decided(&self.sides[1])?])
    }

    /// The outward normal of `face` at a point on it, scaled to at most unit
    /// length: divided by an upper bound of its 1-norm, which is at least
    /// its length.
    fn outward(
        &self,
        face: FaceId,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<Vector3>, Refusal> {
        let surface = self.model.geom().surface(self.model.face(face).surface);
        let Ok(normal) = surface.normal_at(point, dec) else {
            dec.check_deadline()?;
            return Ok(None);
        };
        // The face normal points out of the material.
        let normal = match self.model.face(face).sense {
            Sense::Same => normal,
            Sense::Reversed => normal.scaled(&Real::from_i64(-1)),
        };
        let mut length = Rational::zero();
        for c in [&normal.x, &normal.y, &normal.z] {
            let c = dec.enclose_budgeted_at(c, BITS)?;
            length = length + c.lo().abs().max(c.hi().abs());
        }
        if length.sign() != Sign::Positive {
            return Ok(None);
        }
        Ok(Some(
            normal.scaled(&Real::from(Rational::from_i64(1) / length)),
        ))
    }

    /// A proven positive radius of a ball inside both solids, grown from a
    /// boundary point of `sides[x]` strictly inside the other solid. The
    /// steps inward try each face's normal at the point, then, on an edge or
    /// at a vertex, the mean of those normals.
    fn ball(
        &mut self,
        x: usize,
        sample: &Sample,
        dec: &mut Decider,
    ) -> Result<Option<Rational>, Refusal> {
        let reach = clearance(
            self.model,
            &mut self.faces,
            &self.sides[1 - x],
            &sample.point,
            None,
            dec,
        )?;
        if reach.sign() != Sign::Positive {
            return Ok(None);
        }
        let step = reach.mul_pow2(-1);
        let mut directions = Vec::new();
        for &face in &sample.faces {
            if let Some(n) = self.outward(face, &sample.point, dec)? {
                directions.push(n);
            }
        }
        if directions.len() > 1 {
            let count =
                Real::from(Rational::from_i64(1) / Rational::from_i64(directions.len() as i64));
            let mean = directions
                .iter()
                .skip(1)
                .fold(directions[0].clone(), |s, n| s.plus(n))
                .scaled(&count);
            directions.push(mean);
        }
        for normal in directions {
            dec.check_deadline()?;
            // The direction is at most unit length, so q is at most `step` from p.
            let q = sample
                .point
                .plus(&normal.scaled(&Real::from(-step.clone())));
            if classify(&mut self.scratch, &self.sides[x], &q, dec)? != Class::Inside {
                continue;
            }
            let rho = clearance(
                self.model,
                &mut self.faces,
                &self.sides[x],
                &q,
                Some(&step),
                dec,
            )?;
            if rho.sign() == Sign::Positive {
                return Ok(Some(rho));
            }
        }
        Ok(None)
    }

    /// A ball inside both solids about a point of a dyadic grid over the two
    /// proven boxes' overlap, coarsest level first, as the parallel planar
    /// witness search of `PreparedFace` probes dyadic grid centres.
    fn grid(
        &mut self,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<Option<Answer<BodyRelation>>, Refusal> {
        let (Some(a), Some(b)) = (&self.sides[0].corners, &self.sides[1].corners) else {
            return Ok(None);
        };
        let span: [[Rational; 2]; 3] = std::array::from_fn(|i| {
            [
                a[i][0].clone().max(b[i][0].clone()),
                a[i][1].clone().min(b[i][1].clone()),
            ]
        });
        if span.iter().any(|s| s[1] <= s[0]) {
            return Ok(None);
        }
        for level in 0..GRID_LEVELS {
            let n = 1_i64 << level;
            for cell in 0..n * n * n {
                dec.check_deadline()?;
                let index = [cell % n, cell / n % n, cell / (n * n)];
                let c: [Real; 3] = std::array::from_fn(|i| {
                    let t = Rational::from_ratio(2 * index[i] + 1, 2 * n).expect("positive count");
                    Real::from(&span[i][0] + &(&t * &(&span[i][1] - &span[i][0])))
                });
                let point = Point3::new(c[0].clone(), c[1].clone(), c[2].clone());
                if classify(&mut self.scratch, &self.sides[0], &point, dec)? != Class::Inside
                    || classify(&mut self.scratch, &self.sides[1], &point, dec)? != Class::Inside
                {
                    continue;
                }
                let mut r = None;
                for x in 0..2 {
                    r = Some(clearance(
                        self.model,
                        &mut self.faces,
                        &self.sides[x],
                        &point,
                        r.as_ref(),
                        dec,
                    )?);
                }
                if let Some(r) = r.filter(|r| r.sign() == Sign::Positive) {
                    return Ok(self.bounded(&r, accuracy));
                }
            }
        }
        Ok(None)
    }

    /// The overlap enclosure from a ball of radius `r` inside both solids:
    /// at least 4 r^3, at most the two proven boxes' overlap.
    fn bounded(&self, r: &Rational, accuracy: &Accuracy) -> Option<Answer<BodyRelation>> {
        let (Some(a), Some(b)) = (&self.sides[0].corners, &self.sides[1].corners) else {
            return None;
        };
        let hi = (0..3).fold(Rational::from_i64(1), |v, i| {
            let side = (&a[i][1]).min(&b[i][1]) - (&a[i][0]).max(&b[i][0]);
            v * side.max(Rational::zero())
        });
        let lo = Rational::from_i64(4) * r.clone() * r.clone() * r.clone();
        let overlap = Enclosure { lo, hi };
        let mut marks = self.model.body(self.sides[0].body).marks.clone();
        marks.extend(&self.model.body(self.sides[1].body).marks);
        mark_width(
            &mut marks,
            &overlap,
            accuracy,
            "overlap volume is enclosed by an inscribed ball and the boxes' overlap",
        );
        Some(Answer::with_marks(BodyRelation::Clash { overlap }, marks))
    }
}

/// Clash or clear where the overlap Boolean refused with `refused`. Pairs
/// that stay undecided keep that refusal, with what this search proved
/// added to its detail.
pub(super) fn relation(
    model: &Model,
    first: BodyId,
    second: BodyId,
    accuracy: &Accuracy,
    mut refused: Refusal,
    dec: &mut Decider,
) -> Answered<BodyRelation> {
    // A Boolean that spent the whole allowance leaves nothing to search with.
    if dec.check_deadline().is_err() {
        return Err(refused);
    }
    let mut search = Arrangement {
        model,
        scratch: model.clone(),
        sides: [
            Side::new(model, first, dec)?,
            Side::new(model, second, dec)?,
        ],
        faces: Faces::default(),
        witnesses: Vec::new(),
        contact: false,
    };
    // Boundary samples: vertices, as OpenCASCADE, then one point inside each
    // face, as SMLib.
    let crossing = search.vertices(dec)? || search.face_points(dec)?;
    let mut unproven = None;
    let mut probed = false;
    if !crossing {
        // With every shell sampled outside the other solid, neither solid
        // can be the overlap, so a point inside both settles the clash
        // without proving the faces apart.
        if search.outside() {
            probed = true;
            if let Some(answer) = search.grid(accuracy, dec)? {
                return Ok(answer);
            }
        }
        let apart = search.apart(dec)?;
        unproven = apart.as_ref().err().cloned();
        if let Ok(floor) = &apart
            && let Some(inside) = search.classes()
        {
            let marks = || {
                let mut marks = model.body(first).marks.clone();
                marks.extend(&model.body(second).marks);
                marks
            };
            let none = |v: &[bool]| v.iter().all(|s| !s);
            let all = |v: &[bool]| v.iter().all(|s| *s);
            if none(&inside[0]) && none(&inside[1]) {
                return gap_relation(model, first, second, marks(), Some(floor), accuracy, dec);
            }
            // The intersection's boundary is the inner solid's whole boundary.
            for x in 0..2 {
                if all(&inside[x]) && none(&inside[1 - x]) {
                    match catallax_measure::volume(model, search.sides[x].body, accuracy, dec) {
                        Ok(volume) if volume.value.lo.sign() == Sign::Positive => {
                            let mut marks: MarkLog = marks();
                            marks.extend(&volume.marks);
                            return Ok(Answer::with_marks(
                                BodyRelation::Clash {
                                    overlap: volume.value,
                                },
                                marks,
                            ));
                        }
                        _ => dec.check_deadline()?,
                    }
                }
            }
        }
    }
    // A boundary point strictly inside the other solid proves a clash; the
    // ball grows from it, or from a grid point inside both.
    let witnesses = std::mem::take(&mut search.witnesses);
    for (x, sample) in &witnesses {
        if let Some(r) = search.ball(*x, sample, dec)?
            && let Some(answer) = search.bounded(&r, accuracy)
        {
            return Ok(answer);
        }
    }
    let overlapping = !witnesses.is_empty();
    if (overlapping || crossing)
        && !probed
        && let Some(answer) = search.grid(accuracy, dec)?
    {
        return Ok(answer);
    }
    let why = if overlapping {
        "a boundary point lies strictly inside the other solid, so the overlap is positive, but \
         no ball inside both solids was proven"
            .to_string()
    } else if crossing {
        "the boundaries cross".to_string()
    } else {
        let mut why = unproven.unwrap_or_else(|| "a shell's class is undecided".into());
        if search.contact {
            why.push_str("; a boundary point lies on the other solid's boundary");
        }
        why
    };
    refused.detail = format!("{}; arrangement: {why}", refused.detail);
    Err(refused)
}
