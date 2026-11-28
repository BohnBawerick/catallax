//! Sub-tolerance gaps between the operands' planes, merged once before the
//! boolean decides anything.
//!
//! A tool plane parallel to the plane of a nearby straight-edged target face,
//! with all its faces' corners closer than the 1 nm modelling tolerance to it,
//! is that plane. A tilted plane is not merged: turning a face leaves a later
//! feature sketched on its old plane a new, tilted gap, such as a skin over a
//! slot cut through it. The exporter's wedge repair handles a tilted gap. The tool is copied once with the target's plane in its place: each
//! corner on it is again the meeting of its planes, and each straight edge on
//! it the meeting of its two faces' planes. The exact boolean then decides an
//! exact coincidence, so the gap never becomes a sliver face, and later
//! operations and exports never see it.
//!
//! The answer is no longer exact. Its mark is the largest corner movement,
//! measured as an L1 length, which bounds every boundary point's movement:
//! the faces that move are planar and their edges straight. Nothing changes
//! unless every moved corner and edge is proved on its planes, the movement is
//! strictly under the modelling tolerance and under a tenth of the length of
//! every edge it touches.

use std::collections::{BTreeMap, BTreeSet, btree_map::Entry};

use catallax_numbers::{
    BoundKind, Decider, EvidenceBound, EvidenceUnit, Mark, ModellingTolerance, Rational, Refusal,
    Sign,
};
use catallax_shapes::{
    Curve, CurveId, Line, Meet, MeetOf, Plane, PointDef, PointId, Position, RPoint3, RVector3,
    Surface, SurfaceId,
};
use catallax_topology::{BodyId, FaceId, Loop, Model, VertexId};

/// The mark of a boolean whose tool had a plane merged with a target plane.
pub const SUB_TOLERANCE_MERGE: &str = "boolean.sub-tolerance-merge";

/// `normal · x = offset`.
#[derive(Clone)]
struct Equation {
    normal: RVector3,
    offset: Rational,
}

fn vector(p: &RPoint3) -> RVector3 {
    p.minus(&RPoint3::origin())
}

impl Equation {
    fn of(model: &Model, surface: SurfaceId) -> Option<Self> {
        let Surface::Plane(plane) = model.geom().surface(surface) else {
            return None;
        };
        let normal = plane.normal.as_rational()?;
        let origin = plane.origin.as_rational()?;
        Some(Self {
            offset: normal.dot(&vector(&origin)),
            normal,
        })
    }

    fn value(&self, p: &RPoint3) -> Rational {
        self.normal.dot(&vector(p)) - self.offset.clone()
    }

    /// Whether `p` is strictly closer than `reach` to the plane.
    fn near(&self, p: &RPoint3, reach: &Rational) -> bool {
        let value = self.value(p);
        &value * &value < &(reach * reach) * &self.normal.norm_squared()
    }

    fn same(&self, other: &Self) -> bool {
        self.normal.cross(&other.normal).is_zero()
            && &other.offset * &self.normal.norm_squared()
                == &self.offset * &self.normal.dot(&other.normal)
    }

    /// The same plane, with its normal on the side of `normal`.
    fn facing(&self, normal: &RVector3) -> Self {
        if self.normal.dot(normal).sign() == Sign::Negative {
            Self {
                normal: self.normal.scaled(&-Rational::one()),
                offset: -self.offset.clone(),
            }
        } else {
            self.clone()
        }
    }
}

/// The one point on every plane, from the first three independent ones.
fn meet(planes: &[Equation]) -> Option<RPoint3> {
    let n = planes.len();
    let point = (0..n).find_map(|i| {
        (i + 1..n).find_map(|j| {
            (j + 1..n).find_map(|k| {
                let [p, q, r] = [&planes[i], &planes[j], &planes[k]];
                let det = p.normal.dot(&q.normal.cross(&r.normal));
                if det.is_zero() {
                    return None;
                }
                let sum = q
                    .normal
                    .cross(&r.normal)
                    .scaled(&p.offset)
                    .plus(&r.normal.cross(&p.normal).scaled(&q.offset))
                    .plus(&p.normal.cross(&q.normal).scaled(&r.offset));
                Some(RPoint3::origin().plus(&sum.scaled(&(Rational::one() / det))))
            })
        })
    })?;
    planes
        .iter()
        .all(|plane| plane.value(&point).is_zero())
        .then_some(point)
}

fn l1(v: &RVector3) -> Rational {
    v.x.abs() + v.y.abs() + v.z.abs()
}

fn largest_component(v: &RVector3) -> Rational {
    v.x.abs().max(v.y.abs()).max(v.z.abs())
}

/// Whether `p` lies on the line through `origin` along `dir`.
fn on_line(origin: &RPoint3, dir: &RVector3, p: &RPoint3) -> bool {
    p.minus(origin).cross(dir).is_zero()
}

/// An axis-aligned box.
struct Bounds {
    lo: [Rational; 3],
    hi: [Rational; 3],
}

impl Bounds {
    fn of<'a>(points: impl IntoIterator<Item = &'a RPoint3>) -> Option<Self> {
        let mut bounds: Option<Self> = None;
        for p in points {
            let p = [&p.x, &p.y, &p.z];
            match &mut bounds {
                None => {
                    bounds = Some(Self {
                        lo: p.map(Clone::clone),
                        hi: p.map(Clone::clone),
                    })
                }
                Some(b) => {
                    for axis in 0..3 {
                        if p[axis] < &b.lo[axis] {
                            b.lo[axis] = p[axis].clone();
                        }
                        if p[axis] > &b.hi[axis] {
                            b.hi[axis] = p[axis].clone();
                        }
                    }
                }
            }
        }
        bounds
    }

    /// Whether the boxes come within `reach` of each other on every axis.
    fn meets(&self, other: &Self, reach: &Rational) -> bool {
        (0..3).all(|axis| {
            &self.lo[axis] - reach <= other.hi[axis] && &other.lo[axis] - reach <= self.hi[axis]
        })
    }
}

fn position(model: &Model, vertex: VertexId) -> Option<RPoint3> {
    model
        .geom()
        .position(model.vertex(vertex).point)
        .as_rational()
}

/// The corners of a face bounded by straight edges only.
fn corners(model: &Model, face: FaceId) -> Option<Vec<RPoint3>> {
    let mut vertices = BTreeSet::new();
    for lp in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return None;
        };
        for coedge in coedges {
            let edge = model.edge(model.coedge(*coedge).edge);
            if !matches!(model.geom().curve(edge.curve), Curve::Line(_)) {
                return None;
            }
            vertices.extend([edge.start?, edge.end?]);
        }
    }
    vertices.into_iter().map(|v| position(model, v)).collect()
}

fn vertices_of(model: &Model, face: FaceId) -> BTreeSet<VertexId> {
    let mut vertices = BTreeSet::new();
    for lp in &model.face(face).loops {
        match model.loop_(*lp) {
            Loop::Vertex(v) => {
                vertices.insert(*v);
            }
            Loop::Coedges(coedges) => {
                for coedge in coedges {
                    let edge = model.edge(model.coedge(*coedge).edge);
                    vertices.extend(edge.start.into_iter().chain(edge.end));
                }
            }
        }
    }
    vertices
}

/// For each tool plane within the tolerance of a target plane, that plane,
/// and the work units the search took: one per target face, and per pair of
/// a tool plane and a target face plus one per corner tested.
fn targets(
    model: &Model,
    target: BodyId,
    tool: &[FaceId],
    tolerance: &Rational,
) -> (BTreeMap<SurfaceId, Equation>, u64) {
    let mut work = 0;
    let mut planes = Vec::new();
    for face in model.faces_of_body(target) {
        work += 1;
        let Some(equation) = Equation::of(model, model.face(face).surface) else {
            continue;
        };
        if let Some(bounds) = corners(model, face).and_then(|c| Bounds::of(&c)) {
            planes.push((equation, bounds));
        }
    }
    let mut surfaces: BTreeMap<SurfaceId, Vec<FaceId>> = BTreeMap::new();
    for face in tool {
        surfaces
            .entry(model.face(*face).surface)
            .or_default()
            .push(*face);
    }
    let mut targets = BTreeMap::new();
    'surfaces: for (surface, faces) in surfaces {
        let Some(own) = Equation::of(model, surface) else {
            continue;
        };
        let mut points = Vec::new();
        for face in faces {
            let Some(corners) = corners(model, face) else {
                continue 'surfaces;
            };
            points.extend(corners);
        }
        let Some(bounds) = Bounds::of(&points) else {
            continue;
        };
        let mut found: Option<&Equation> = None;
        for (plane, near) in &planes {
            work += 1 + points.len() as u64;
            if !near.meets(&bounds, tolerance) {
                continue;
            }
            // An exact coincidence needs no merge, and leaves this plane as it is.
            if own.same(plane) {
                continue 'surfaces;
            }
            if !plane.normal.cross(&own.normal).is_zero()
                || !points.iter().all(|p| plane.near(p, tolerance))
            {
                continue;
            }
            match found {
                // Two different target planes within reach: neither is chosen.
                Some(earlier) if !earlier.same(plane) => continue 'surfaces,
                _ => found = Some(plane),
            }
        }
        if let Some(plane) = found {
            targets.insert(surface, plane.facing(&own.normal));
        }
    }
    (targets, work)
}

/// Copies `tool` with each plane within the modelling tolerance of a target
/// plane replaced by that plane, or `None` when no plane merges.
pub(crate) fn merge(
    model: &mut Model,
    target: BodyId,
    tool: BodyId,
    dec: &mut Decider,
) -> Result<Option<BodyId>, Refusal> {
    if model.blend_base(tool).is_some() || model.chart_expansion(tool).is_some() {
        return Ok(None);
    }
    let tolerance = ModellingTolerance::default();
    let reach = tolerance.millimetres().clone();
    let faces = model.faces_of_body(tool);
    // Exact rational tests over the operands' corners alone. A part with no
    // such gap is charged nothing, so its answer and work stay as they were.
    let (targets, searched) = targets(model, target, &faces, &reach);
    if targets.is_empty() {
        return Ok(None);
    }
    dec.charge_traversal(searched)?;
    let plane = |model: &Model, surface: &SurfaceId| {
        targets
            .get(surface)
            .cloned()
            .or_else(|| Equation::of(model, *surface))
    };

    // Every tool vertex's planes; those on a merged plane move.
    let mut incident: BTreeMap<VertexId, BTreeSet<SurfaceId>> = BTreeMap::new();
    let mut moving = BTreeSet::new();
    for face in &faces {
        let surface = model.face(*face).surface;
        for vertex in vertices_of(model, *face) {
            incident.entry(vertex).or_default().insert(surface);
            if targets.contains_key(&surface) {
                moving.insert(vertex);
            }
        }
    }
    let mut positions = BTreeMap::new();
    let mut bound = Rational::zero();
    for vertex in &moving {
        dec.charge_traversal(1)?;
        let mut planes: Vec<Equation> = Vec::new();
        for surface in &incident[vertex] {
            let Some(equation) = plane(model, surface) else {
                return Ok(None);
            };
            if !planes.iter().any(|p| p.same(&equation)) {
                planes.push(equation);
            }
        }
        let (Some(new), Some(old)) = (meet(&planes), position(model, *vertex)) else {
            return Ok(None);
        };
        bound = bound.max(l1(&new.minus(&old)));
        positions.insert(*vertex, new);
    }
    if bound.is_zero()
        || tolerance
            .check(&EvidenceBound {
                upper: bound.clone(),
                unit: EvidenceUnit::Millimetres,
            })
            .is_err()
    {
        return Ok(None);
    }

    // Every straight edge at a moved corner keeps its direction and a tenth
    // of its length; one on a merged plane is the meeting of its new planes.
    let tool_faces: BTreeSet<FaceId> = faces.iter().copied().collect();
    let at = |model: &Model, vertex: VertexId| {
        positions
            .get(&vertex)
            .cloned()
            .or_else(|| position(model, vertex))
    };
    let edges = model.edges_of_body(tool);
    let mut lines: BTreeMap<CurveId, (RPoint3, RVector3)> = BTreeMap::new();
    for id in &edges {
        let edge = model.edge(*id).clone();
        let owners: BTreeSet<FaceId> = model
            .coedges_of_edge(*id)
            .into_iter()
            .filter_map(|coedge| model.owner_of_coedge(coedge))
            .map(|(_, face)| face)
            .filter(|face| tool_faces.contains(face))
            .collect();
        let on_merged = owners
            .iter()
            .any(|face| targets.contains_key(&model.face(*face).surface));
        let touched = [edge.start, edge.end]
            .into_iter()
            .flatten()
            .any(|vertex| positions.contains_key(&vertex));
        if !on_merged && !touched {
            continue;
        }
        dec.charge_traversal(1)?;
        let Curve::Line(line) = model.geom().curve(edge.curve) else {
            return Ok(None);
        };
        let (Some(origin), Some(dir)) = (line.origin.as_rational(), line.dir.as_rational()) else {
            return Ok(None);
        };
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Ok(None);
        };
        let (Some(old_start), Some(old_end), Some(start), Some(end)) = (
            position(model, start),
            position(model, end),
            at(model, start),
            at(model, end),
        ) else {
            return Ok(None);
        };
        let span = end.minus(&start);
        if span.dot(&old_end.minus(&old_start)).sign() != Sign::Positive
            || &bound * &Rational::from(10) >= largest_component(&span)
        {
            return Ok(None);
        }
        let (origin, dir) = if on_merged {
            let surfaces: Vec<SurfaceId> = owners
                .iter()
                .map(|face| model.face(*face).surface)
                .collect();
            let [p, q] = surfaces.as_slice() else {
                return Ok(None);
            };
            let (Some(p), Some(q)) = (plane(model, p), plane(model, q)) else {
                return Ok(None);
            };
            let meeting = p.normal.cross(&q.normal);
            let meeting = match meeting.dot(&dir).sign() {
                Sign::Zero => return Ok(None),
                Sign::Negative => meeting.scaled(&-Rational::one()),
                Sign::Positive => meeting,
            };
            (start.clone(), meeting)
        } else {
            (origin, dir)
        };
        if !on_line(&origin, &dir, &start) || !on_line(&origin, &dir, &end) {
            return Ok(None);
        }
        if on_merged {
            match lines.entry(edge.curve) {
                Entry::Vacant(entry) => {
                    entry.insert((origin, dir));
                }
                Entry::Occupied(entry) => {
                    let (o, d) = entry.get();
                    if !on_line(o, d, &start) || !on_line(o, d, &end) {
                        return Ok(None);
                    }
                }
            }
        }
    }
    // A replaced line or point that another edge or vertex shares must still
    // hold that edge's corners and that vertex.
    for id in &edges {
        let edge = model.edge(*id);
        if let Some((origin, dir)) = lines.get(&edge.curve) {
            for vertex in [edge.start, edge.end].into_iter().flatten() {
                match at(model, vertex) {
                    Some(p) if on_line(origin, dir, &p) => {}
                    _ => return Ok(None),
                }
            }
        }
    }
    let moved_points: BTreeSet<PointId> = positions
        .keys()
        .map(|vertex| model.vertex(*vertex).point)
        .collect();
    if model.vertices_of_body(tool).into_iter().any(|vertex| {
        !positions.contains_key(&vertex) && moved_points.contains(&model.vertex(vertex).point)
    }) {
        return Ok(None);
    }

    // Each merged plane keeps its own frame, moved onto the target plane.
    let mut surfaces = BTreeMap::new();
    for (surface, equation) in &targets {
        let Surface::Plane(plane) = model.geom().surface(*surface).clone() else {
            return Ok(None);
        };
        let (Some(origin), Some(ref_dir)) =
            (plane.origin.as_rational(), plane.ref_dir.as_rational())
        else {
            return Ok(None);
        };
        let norm = equation.normal.norm_squared();
        let origin = origin.plus(
            &equation
                .normal
                .scaled(&(-equation.value(&origin) / norm.clone())),
        );
        let ref_dir = ref_dir.minus(
            &equation
                .normal
                .scaled(&(equation.normal.dot(&ref_dir) / norm)),
        );
        if ref_dir.is_zero() {
            return Ok(None);
        }
        let id = model.geom_mut().add_surface(Surface::Plane(Plane {
            origin: origin.to_real(),
            normal: equation.normal.to_real(),
            ref_dir: ref_dir.to_real(),
        }));
        surfaces.insert(*surface, id);
    }
    let mut curves = BTreeMap::new();
    for (curve, (origin, dir)) in lines {
        let id = model.geom_mut().add_curve(Curve::Line(Line {
            origin: origin.to_real(),
            dir: dir.to_real(),
        }));
        curves.insert(curve, id);
    }
    let mut points = BTreeMap::new();
    for (vertex, new) in &positions {
        let source = model.vertex(*vertex).point;
        if position(model, *vertex).as_ref() == Some(new) {
            continue;
        }
        let of = incident[vertex]
            .iter()
            .map(|surface| surfaces.get(surface).copied().unwrap_or(*surface))
            .collect();
        let id = model.geom_mut().add_point(PointDef::Meet(Meet {
            of: MeetOf::Surfaces(of),
            position: Position::Closed(new.to_real()),
        }));
        points.insert(source, id);
    }

    let copy = catallax_topology::substituted_body(model, tool, surfaces, curves, points);
    let mut body = model.body(copy).clone();
    body.marks.push(Mark::new(
        SUB_TOLERANCE_MERGE,
        BoundKind::Distance,
        bound.clone(),
        format!(
            "{} tool plane(s) within the {} mm modelling tolerance of the target's merged with them; corners moved at most {bound} mm",
            targets.len(),
            tolerance.millimetres()
        ),
    ));
    Ok(Some(model.add_body(body)))
}
