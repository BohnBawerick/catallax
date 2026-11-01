//! Rail-first elementary fillet construction over an immutable input skin.

use std::collections::{BTreeMap, BTreeSet};

use catallax_intersect::{Component, CurveSurface, Cyclic, SurfacePair};
use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign};
use catallax_recipe::{EdgeSelection, Origin, Role, StepId};
use catallax_shapes::{
    Circle, Curve, Line, Point3, PointDef, PointId, Surface, Torus, TorusSheet, Vector3,
};
use catallax_topology::{BodyId, Edge, EdgeId, FaceId, Loop, Model, Sense, Vertex, VertexId};

use crate::{
    edge_blend::{self, Replacement},
    extrude::reverse,
    fillet::{FilletDefinition, separated},
    fillet_corner::CornerDefinition,
    profile::sign,
};

#[path = "fillet_envelope.rs"]
mod fillet_envelope;
#[path = "fillet_tangent.rs"]
mod fillet_tangent;

use fillet_envelope::Clearance;

type Uses = BTreeMap<EdgeId, Vec<(FaceId, Sense)>>;

fn bad(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-topology", detail)
}
fn vertex(model: &mut Model, p: Point3) -> VertexId {
    let point = model.geom_mut().add_point(PointDef::Explicit(p));
    model.add_vertex(Vertex { point })
}
fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}
fn uses(model: &Model, body: BodyId) -> Uses {
    let mut uses: Uses = BTreeMap::new();
    for face in model.faces_of_body(body) {
        for lp in &model.face(face).loops {
            if let Loop::Coedges(cs) = model.loop_(*lp) {
                for c in cs {
                    let c = model.coedge(*c);
                    uses.entry(c.edge).or_default().push((face, c.sense));
                }
            }
        }
    }
    uses
}
fn matching(model: &Model, face: FaceId, pattern: &catallax_recipe::OriginPattern) -> bool {
    model.face_origins(face).iter().any(|o| pattern.matches(o))
}
/// Chooses the edge at strictly the smallest distance from a stored `near`
/// point, measured to its closed trimmed set. Only a straight segment has that
/// exact distance here, so curved candidates keep the refusal.
fn nearest(
    model: &Model,
    edges: &[EdgeId],
    near: &catallax_recipe::Point3,
    dec: &mut Decider,
) -> Result<EdgeId, Refusal> {
    let [x, y, z] = near.0.each_ref().map(|n| n.value().clone().into());
    let point = Point3::new(x, y, z);
    let mut distances = Vec::with_capacity(edges.len());
    for edge in edges {
        let e = model.edge(*edge);
        let (Curve::Line(_), Some(start), Some(end)) =
            (model.geom().curve(e.curve), e.start, e.end)
        else {
            return Err(bad(
                "nearest-edge selection among curved edges needs an exact trimmed-edge distance query",
            ));
        };
        let [a, b] = [start, end].map(|v| model.geom().position(model.vertex(v).point));
        let along = b.minus(&a);
        let offset = point.minus(&a);
        let foot = offset.dot(&along);
        let length_sq = along.norm_squared();
        distances.push(if sign(&foot, dec)? != Sign::Positive {
            offset.norm_squared()
        } else if sign(&(&foot - &length_sq), dec)? != Sign::Negative {
            point.minus(&b).norm_squared()
        } else {
            offset.norm_squared() - foot.square() / length_sq
        });
    }
    Ok(edges[crate::draft::near::unique_nearest(&distances, dec)?])
}
fn resolve(
    model: &Model,
    selection: &EdgeSelection,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Vec<EdgeId>, Refusal> {
    match selection {
        EdgeSelection::All => Ok(uses.keys().copied().collect()),
        EdgeSelection::Refs(refs) => {
            let mut selected = BTreeSet::new();
            for reference in refs {
                let mut found = Vec::new();
                for (edge, incident) in uses {
                    let [(a, _), (b, _)] = incident.as_slice() else {
                        continue;
                    };
                    if !((matching(model, *a, &reference.between[0])
                        && matching(model, *b, &reference.between[1]))
                        || (matching(model, *b, &reference.between[0])
                            && matching(model, *a, &reference.between[1])))
                    {
                        continue;
                    }
                    let e = model.edge(*edge);
                    let touches = reference.touching.iter().all(|pattern| {
                        uses.iter().any(|(other, incident)| {
                            let other = model.edge(*other);
                            [e.start, e.end]
                                .into_iter()
                                .flatten()
                                .any(|v| other.start == Some(v) || other.end == Some(v))
                                && incident.iter().any(|(f, _)| matching(model, *f, pattern))
                        })
                    });
                    if touches {
                        found.push(*edge);
                    }
                }
                if found.is_empty() {
                    return Err(Refusal::new(
                        "fillet.lost-edge",
                        format!("no input edge resolves {reference:?}"),
                    ));
                }
                if let Some(near) = &reference.near
                    && found.len() > 1
                {
                    found = vec![nearest(model, &found, near, dec)?];
                }
                selected.extend(found);
            }
            Ok(selected.into_iter().collect())
        }
    }
}

struct Run {
    edge: EdgeId,
    faces: [FaceId; 2],
    senses: [Sense; 2],
    definition: FilletDefinition,
    ends: BTreeMap<(VertexId, FaceId), VertexId>,
    caps: BTreeMap<VertexId, EdgeId>,
}
fn outward(
    model: &Model,
    face: FaceId,
    at: &Point3,
    dec: &mut Decider,
) -> Result<catallax_shapes::Vector3, Refusal> {
    let f = model.face(face);
    let surface = model.geom().surface(f.surface);
    let n = surface.normal_at(at, dec)?;
    if sign(&n.norm_squared(), dec)? == Sign::Zero {
        return Err(bad("singular support has no normal"));
    }
    Ok(n.scaled(&Real::from_i64(if f.sense == Sense::Same { 1 } else { -1 })))
}
fn convex(
    model: &Model,
    edge: EdgeId,
    incident: &[(FaceId, Sense)],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [(a, sense), (b, _)] = incident else {
        return Err(bad("a fillet edge must have two incident faces"));
    };
    let e = model.edge(edge);
    let orientation = if *sense == e.sense { 1 } else { -1 };
    match model.geom().curve(e.curve) {
        Curve::Circle(c) => {
            let plane = [*a, *b]
                .into_iter()
                .find(|f| {
                    matches!(
                        model.geom().surface(model.face(*f).surface),
                        Surface::Plane(_)
                    )
                })
                .ok_or_else(|| bad("circular elementary fillets need one plane"))?;
            let coedge = incident.iter().find(|(f, _)| *f == plane).unwrap().1;
            let orientation = if coedge == e.sense { 1 } else { -1 };
            let n = outward(model, plane, &c.center, dec)?;
            Ok(sign(&(c.normal.dot(&n) * Real::from_i64(orientation)), dec)? == Sign::Positive)
        }
        Curve::Line(l) => {
            let na = outward(model, *a, &l.origin, dec)?;
            let nb = outward(model, *b, &l.origin, dec)?;
            let value = nb.dot(&na.cross(&l.dir)) * Real::from_i64(orientation);
            match sign(&value, dec)? {
                Sign::Negative => Ok(true),
                Sign::Positive => Ok(false),
                Sign::Zero => Err(bad("tangent faces do not define a regular sharp edge")),
            }
        }
        _ => Err(bad("elementary fillets need a line or circle")),
    }
}
fn line_plane(
    line: &Line,
    plane: &catallax_shapes::Plane,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let divisor = line.dir.dot(&plane.normal);
    if sign(&divisor, dec)? == Sign::Zero {
        return Err(bad("end face is not transversal to the contact rail"));
    }
    Ok(line.origin.plus(
        &line
            .dir
            .scaled(&(plane.origin.minus(&line.origin).dot(&plane.normal) / divisor)),
    ))
}
fn line_line(a: &Line, b: &Line, dec: &mut Decider) -> Result<Point3, Refusal> {
    let cross = a.dir.cross(&b.dir);
    let divisor = cross.norm_squared();
    if sign(&divisor, dec)? == Sign::Zero {
        return Err(bad("parallel contact rails need chain continuation"));
    }
    let t = b.origin.minus(&a.origin).cross(&b.dir).dot(&cross) / divisor;
    let p = a.origin.plus(&a.dir.scaled(&t));
    if sign(&p.minus(&b.origin).cross(&b.dir).norm_squared(), dec)? != Sign::Zero {
        return Err(bad("contact rails do not meet"));
    }
    Ok(p)
}
fn rail(run: &Run, face: FaceId) -> Result<&Line, Refusal> {
    let index = run
        .faces
        .iter()
        .position(|f| *f == face)
        .ok_or_else(|| bad("rail has no input support"))?;
    let Curve::Line(line) = &run.definition.contacts[index] else {
        return Err(bad("vertex trimming needs straight rails"));
    };
    Ok(line)
}
fn decided<T>(outcome: Outcome<T>) -> Result<Option<T>, Refusal> {
    match outcome {
        Outcome::Proven(value) => Ok(Some(value)),
        Outcome::ProvenCoincident(_) => Ok(None),
        Outcome::Refused(e) => Err(e),
        Outcome::NotExact(e) => Err(Refusal::new(e.code, e.detail)),
    }
}
/// Whether a point of the edge's carrier lies on its trimmed arc, ends included.
fn on_trimmed(model: &Model, edge: EdgeId, p: PointId, dec: &mut Decider) -> Result<bool, Refusal> {
    let e = model.edge(edge);
    let (Some(start), Some(end)) = (e.start, e.end) else {
        return Err(bad("a curved end wall needs a bounded input edge"));
    };
    let geom = model.geom();
    if decided(catallax_intersect::on_curve(geom, e.curve, p, dec))?.is_some() {
        return Ok(false);
    }
    let [start, end] = [start, end].map(|v| model.vertex(v).point);
    let point = geom.position(p);
    for v in [start, end] {
        if same(&point, &geom.position(v), dec)? {
            return Ok(true);
        }
    }
    let (from, to) = if e.sense == Sense::Same {
        (start, end)
    } else {
        (end, start)
    };
    if catallax_intersect::curve_is_closed(geom, e.curve, dec)? {
        let order = catallax_intersect::cyclic_order(geom, e.curve, [from, p, to], dec);
        return Ok(decided(order)? == Some(Cyclic::Forward));
    }
    let first = decided(catallax_intersect::order_on_curve(
        geom, e.curve, from, p, dec,
    ))?;
    let second = decided(catallax_intersect::order_on_curve(
        geom, e.curve, p, to, dec,
    ))?;
    Ok(first.is_some() && first == second)
}
/// Where a straight rail leaves its support across a curved end wall: the
/// carrier hit on that support's input edge to the wall at this vertex. The
/// carrier's other hit lies beyond the edge and is never a strip end.
fn wall_foot(
    model: &mut Model,
    uses: &Uses,
    run: &Run,
    face: FaceId,
    wall: FaceId,
    old: VertexId,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let boundary = uses
        .iter()
        .find(|(edge, incident)| {
            let e = model.edge(**edge);
            (e.start == Some(old) || e.end == Some(old))
                && incident.iter().any(|(f, _)| *f == face)
                && incident.iter().any(|(f, _)| *f == wall)
        })
        .map(|(edge, _)| *edge)
        .ok_or_else(|| bad("a strip support does not meet its end wall at the vertex"))?;
    let line = model
        .geom_mut()
        .add_curve(Curve::Line(rail(run, face)?.clone()));
    let surface = model.face(wall).surface;
    let meeting =
        catallax_intersect::intersect_curve_surface(model.geom_mut(), line, surface, dec)?.value;
    let CurveSurface::Meets(hits) = meeting else {
        return Err(bad("a contact rail does not cross its curved end wall"));
    };
    let mut foot = None;
    for hit in hits {
        if on_trimmed(model, boundary, hit.point, dec)? && foot.replace(hit.point).is_some() {
            return Err(bad("a contact rail meets its end-wall edge twice"));
        }
    }
    let foot = foot.ok_or_else(|| {
        Refusal::new(
            "fillet.unsupported-overflow",
            "a contact rail leaves its end-wall edge and needs neighbour continuation",
        )
    })?;
    Ok(model.geom().position(foot))
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &v.norm_squared().sqrt()))
}
/// The ball centre of a strip in the normal plane of its edge through `at`.
fn spine(definition: &FilletDefinition, at: &Point3) -> Result<Point3, Refusal> {
    match &definition.surface {
        Surface::Cylinder(c) => {
            let d = at.minus(&c.origin);
            Ok(c.origin
                .plus(&c.axis.scaled(&(d.dot(&c.axis) / c.axis.norm_squared()))))
        }
        Surface::Torus(Torus {
            center,
            axis,
            major,
            ..
        })
        | Surface::TorusSheet(TorusSheet {
            center,
            axis,
            major,
            ..
        }) => {
            let d = at.minus(center);
            let radial = d.minus(&axis.scaled(&(d.dot(axis) / axis.norm_squared())));
            Ok(center.plus(&unit(&radial).scaled(major)))
        }
        _ => Err(bad("tangent strip joins need cylinder or torus strips")),
    }
}
/// Two strips meeting where their edges are tangent and their outer supports
/// are tangent share one ball and one cross-section circle: the meridian
/// through the vertex. No corner face has finite area there. Returns that
/// circle and the common contacts on the shared and the outer supports.
fn tangent_join(
    model: &mut Model,
    first: &Run,
    second: &Run,
    common: FaceId,
    outer: [FaceId; 2],
    old: VertexId,
    dec: &mut Decider,
) -> Result<Option<(Circle, Point3, Point3)>, Refusal> {
    let point = model.vertex(old).point;
    let at = model.geom().position(point);
    let [t, u] = [first, second].map(|run| {
        catallax_intersect::tangent_at(model.geom(), model.edge(run.edge).curve, point, dec)
            .map(|[_, leaving]| leaving)
    });
    let (t, u) = (t?, u?);
    if sign(&t.cross(&u).norm_squared(), dec)? != Sign::Zero {
        return Ok(None);
    }
    if first.definition.convex != second.definition.convex {
        return Err(bad(
            "a tangent join of convex and concave strips needs a transition",
        ));
    }
    let centre = spine(&first.definition, &at)?;
    if !same(&centre, &spine(&second.definition, &at)?, dec)? {
        return Err(bad(
            "tangent strips have different ball centres at their join",
        ));
    }
    let reach = if first.definition.convex {
        first.definition.radius.clone()
    } else {
        -&first.definition.radius
    };
    let mut foot = |face: FaceId| -> Result<Point3, Refusal> {
        Ok(centre.plus(&unit(&outward(model, face, &at, dec)?).scaled(&reach)))
    };
    let p = foot(common)?;
    let [qa, qb] = [foot(outer[0])?, foot(outer[1])?];
    if !same(&qa, &qb, dec)? {
        return Err(bad(
            "tangent strips need tangent outer supports at their join",
        ));
    }
    // Each contact must lie on its own strip's contact curve.
    for (run, face, at) in [
        (first, common, &p),
        (second, common, &p),
        (first, outer[0], &qa),
        (second, outer[1], &qa),
    ] {
        let index = run.faces.iter().position(|f| *f == face).unwrap();
        let curve = model
            .geom_mut()
            .add_curve(run.definition.contacts[index].clone());
        let id = model.geom_mut().add_point(PointDef::Explicit(at.clone()));
        if decided(catallax_intersect::on_curve(model.geom(), curve, id, dec))?.is_some() {
            return Err(bad(
                "a tangent join contact leaves its strip's contact curve",
            ));
        }
    }
    let circle = Circle {
        center: centre.clone(),
        normal: t,
        radius_sq: first.definition.radius.square(),
        ref_dir: p.minus(&centre),
    };
    Ok(Some((circle, p, qa)))
}
/// The shorter meridian arc from the shared contact to the outer contact.
fn meridian_edge(
    model: &mut Model,
    meridian: Circle,
    from: VertexId,
    to: VertexId,
    dec: &mut Decider,
) -> Result<EdgeId, Refusal> {
    let [p, q] = [from, to].map(|v| model.geom().position(model.vertex(v).point));
    let turn = p
        .minus(&meridian.center)
        .cross(&q.minus(&meridian.center))
        .dot(&meridian.normal);
    let sense = match sign(&turn, dec)? {
        Sign::Positive => Sense::Same,
        Sign::Negative => Sense::Reversed,
        Sign::Zero => return Err(bad("a tangent join meridian has no shorter arc")),
    };
    let curve = model.geom_mut().add_curve(Curve::Circle(meridian));
    Ok(model.add_edge(Edge {
        curve,
        sense,
        start: Some(from),
        end: Some(to),
    }))
}
/// Whether two analytic faces are tangent with agreeing normals at a point.
/// A parametric face has no event normal here and never continues a run.
fn tangent_faces(
    model: &Model,
    a: FaceId,
    b: FaceId,
    at: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if [a, b].into_iter().any(|f| {
        matches!(
            model.geom().surface(model.face(f).surface),
            Surface::BSpline(_)
        )
    }) {
        return Ok(false);
    }
    let na = outward(model, a, at, dec)?;
    let nb = outward(model, b, at, dec)?;
    Ok(sign(&na.cross(&nb).norm_squared(), dec)? == Sign::Zero
        && sign(&na.dot(&nb), dec)? == Sign::Positive)
}
/// Version 4 roll-over at a smooth input boundary. Where a contact foot
/// reaches an edge across which its face continues tangentially, the run
/// continues along the edge between the strip's other face and the next face.
/// A continued edge keeps the blend name of the selected edge it continues;
/// explicitly selected edges keep their own. Connection alone never starts a
/// run: the faces must be tangent at the event vertex.
fn continue_tangent(
    model: &Model,
    uses: &Uses,
    selected: Vec<EdgeId>,
    dec: &mut Decider,
) -> Result<Vec<(EdgeId, EdgeId)>, Refusal> {
    let touches = |edge: EdgeId, v: VertexId| {
        let e = model.edge(edge);
        e.start == Some(v) || e.end == Some(v)
    };
    let between = |edge: &EdgeId, a: FaceId, b: FaceId| {
        uses[edge].iter().any(|(f, _)| *f == a) && uses[edge].iter().any(|(f, _)| *f == b)
    };
    let mut sources: BTreeMap<EdgeId, EdgeId> = selected.iter().map(|e| (*e, *e)).collect();
    let mut pending = selected;
    while let Some(edge) = pending.pop() {
        let source = sources[&edge];
        let [(first, _), (second, _)] = uses[&edge].as_slice() else {
            continue;
        };
        let e = model.edge(edge);
        for v in [e.start, e.end].into_iter().flatten() {
            let at = model.geom().position(model.vertex(v).point);
            for (stay, leave) in [(*first, *second), (*second, *first)] {
                for (boundary, incident) in uses {
                    if *boundary == edge || !touches(*boundary, v) {
                        continue;
                    }
                    let [(a, _), (b, _)] = incident.as_slice() else {
                        continue;
                    };
                    let next = match (*a == leave, *b == leave) {
                        (true, false) => *b,
                        (false, true) => *a,
                        _ => continue,
                    };
                    if next == stay || !tangent_faces(model, leave, next, &at, dec)? {
                        continue;
                    }
                    if let Some(continued) = uses
                        .keys()
                        .find(|h| **h != edge && touches(**h, v) && between(h, stay, next))
                        && !sources.contains_key(continued)
                    {
                        sources.insert(*continued, source);
                        pending.push(*continued);
                    }
                }
            }
        }
    }
    Ok(sources.into_iter().collect())
}
fn plane(model: &Model, face: FaceId) -> Result<&catallax_shapes::Plane, Refusal> {
    let Surface::Plane(p) = model.geom().surface(model.face(face).surface) else {
        return Err(bad(
            "transversal line-strip ends currently require planar neighbours",
        ));
    };
    Ok(p)
}

fn connector(
    model: &mut Model,
    a: &Surface,
    b: &Surface,
    from: VertexId,
    to: VertexId,
    sectors: &[&FilletDefinition],
    dec: &mut Decider,
) -> Result<EdgeId, Refusal> {
    let a = model.geom_mut().add_surface(a.clone());
    let b = model.geom_mut().add_surface(b.clone());
    let meeting = catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value;
    let SurfacePair::Meets(components) = meeting else {
        return Err(bad("blend end supports have no regular common curve"));
    };
    let p = model.vertex(from).point;
    let q = model.vertex(to).point;
    for component in components {
        let Component::Curve { curve, .. } = component else {
            continue;
        };
        let on = |model: &Model, p, dec: &mut Decider| -> Result<bool, Refusal> {
            match catallax_intersect::on_curve(model.geom(), curve, p, dec) {
                Outcome::ProvenCoincident(_) => Ok(true),
                Outcome::Proven(_) => Ok(false),
                Outcome::Refused(e) => Err(e),
                Outcome::NotExact(e) => Err(Refusal::new(e.code, e.detail)),
            }
        };
        if !on(model, p, dec)? || !on(model, q, dec)? {
            continue;
        }
        for (sense, start, end) in [(Sense::Same, p, q), (Sense::Reversed, q, p)] {
            let sample =
                catallax_intersect::arc_point(model.geom(), curve, Some(start), Some(end), dec)?
                    .value;
            let mut inside = true;
            for sector in sectors {
                inside &= sector.in_sector(&sample, dec)?;
            }
            if inside {
                return Ok(model.add_edge(Edge {
                    curve,
                    sense,
                    start: Some(from),
                    end: Some(to),
                }));
            }
        }
    }
    Err(bad("no end-curve arc lies in every active contact sector"))
}

fn guard(
    model: &mut Model,
    body: BodyId,
    run: &Run,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let edge = model.edge(run.edge);
    let ends = [edge.start, edge.end];
    let neighbours: BTreeSet<_> = uses
        .iter()
        .filter(|(e, _)| {
            let e = model.edge(**e);
            ends.into_iter()
                .flatten()
                .any(|v| e.start == Some(v) || e.end == Some(v))
        })
        .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
        .collect();
    for face in model.faces_of_body(body) {
        if run.faces.contains(&face) || neighbours.contains(&face) {
            continue;
        }
        let other = model.geom().surface(model.face(face).surface).clone();
        if separated(&run.definition.surface, &other, dec)? {
            continue;
        }
        let a = model.geom_mut().add_surface(run.definition.surface.clone());
        let b = model.geom_mut().add_surface(other);
        if !matches!(
            catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value,
            SurfacePair::Disjoint
        ) {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a non-end neighbour intersects the blend carrier; its trimmed contact continuation must be resolved before publishing",
            ));
        }
    }
    Ok(())
}

fn corner_guard(
    model: &mut Model,
    body: BodyId,
    corner: &CornerDefinition,
    incident: &BTreeSet<FaceId>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for face in model.faces_of_body(body) {
        if incident.contains(&face) {
            continue;
        }
        let other = model.geom().surface(model.face(face).surface).clone();
        if separated(&corner.surface, &other, dec)? {
            continue;
        }
        let a = model.geom_mut().add_surface(corner.surface.clone());
        let b = model.geom_mut().add_surface(other);
        if !matches!(
            catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value,
            SurfacePair::Disjoint
        ) {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a non-incident face meets the corner sphere; its trimmed continuation must be resolved before publishing",
            ));
        }
    }
    Ok(())
}

pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    selection: &EdgeSelection,
    radius: &Real,
    roll_over: bool,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<FilletDefinition>), Refusal> {
    if sign(radius, dec)? != Sign::Positive {
        return Err(Refusal::new(
            "fillet.non-positive-radius",
            "fillet radius must be positive",
        ));
    }
    let uses = uses(model, body);
    let selected = resolve(model, selection, &uses, dec)?;
    if selected.is_empty() {
        return Ok((body, Vec::new()));
    }
    let selected = if roll_over {
        continue_tangent(model, &uses, selected, dec)?
    } else {
        selected.into_iter().map(|edge| (edge, edge)).collect()
    };
    let mut runs = Vec::with_capacity(selected.len());
    for (edge, source) in selected {
        let incident = &uses[&edge];
        let [(a, sa), (b, sb)] = incident.as_slice() else {
            return Err(bad("a selected edge must have two incident faces"));
        };
        let name = model
            .edge_name(body, source)
            .ok_or_else(|| bad("a selected edge needs an unambiguous origin-pair name"))?;
        let faces = [*a, *b];
        let definition = FilletDefinition::new(
            Origin {
                step: step.clone(),
                role: Role::Blend { edge: name },
            },
            faces.map(|f| model.geom().surface(model.face(f).surface).clone()),
            faces.map(|f| model.face(f).sense == Sense::Same),
            model.geom().curve(model.edge(edge).curve).clone(),
            radius.clone(),
            convex(model, edge, incident, dec)?,
            dec,
        )?;
        runs.push(Run {
            edge,
            faces,
            senses: [*sa, *sb],
            definition,
            ends: BTreeMap::new(),
            caps: BTreeMap::new(),
        });
    }
    let trimming = runs
        .iter()
        .any(|run| {
            matches!(
                run.definition.supports,
                [Surface::Cylinder(_), Surface::Cylinder(_)]
            )
        })
        .then(|| catallax_numbers::OperationSpan::begin("brep", step, "fillet.trimming", dec));
    for run in &runs {
        if let Err(refusal) = guard(model, body, run, &uses, dec) {
            // A lone closed ring may clear every countersink mouth, touch them
            // at single points, or continue over them; every other overflow,
            // and any ring outside that stratum, still refuses.
            let continuable = refusal.code == "fillet.unsupported-overflow"
                || refusal.code.starts_with("intersect.");
            if !continuable || runs.len() != 1 {
                return Err(refusal);
            }
            match fillet_envelope::strip_clear(
                model,
                body,
                run.edge,
                run.faces,
                &run.definition,
                dec,
            )? {
                Clearance::Clear => continue,
                Clearance::Tangent(mouths) => {
                    if let Some(result) = fillet_tangent::build(
                        model,
                        body,
                        run.edge,
                        run.faces,
                        &run.definition,
                        &mouths,
                        &uses,
                        dec,
                    )? {
                        return Ok((result, vec![run.definition.clone()]));
                    }
                }
                Clearance::Blocked => {}
            }
            if let Some(result) = fillet_envelope::build(
                model,
                body,
                run.edge,
                run.faces,
                &run.definition,
                &uses,
                dec,
            )? {
                return Ok((result, vec![run.definition.clone()]));
            }
            return Err(refusal);
        }
    }
    for (i, a) in runs.iter().enumerate() {
        for b in &runs[i + 1..] {
            let a_edge = model.edge(a.edge);
            let b_edge = model.edge(b.edge);
            if [a_edge.start, a_edge.end]
                .into_iter()
                .flatten()
                .any(|v| b_edge.start == Some(v) || b_edge.end == Some(v))
            {
                continue;
            }
            for (index, face) in a.faces.iter().enumerate() {
                if let Some(other) = b.faces.iter().position(|f| f == face)
                    && a.definition
                        .wall_arc_meets(index, &b.definition, other, dec)?
                {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "two strips consume one arc of their shared cylinder wall and require simultaneous contact-region trimming",
                    ));
                }
            }
            if separated(&a.definition.surface, &b.definition.surface, dec)?
                || a.definition.sectors_apart(&b.definition, dec)?
            {
                continue;
            }
            let a = model.geom_mut().add_surface(a.definition.surface.clone());
            let b = model.geom_mut().add_surface(b.definition.surface.clone());
            if !matches!(
                catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value,
                SurfacePair::Disjoint
            ) {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "nonincident blend carriers meet and require simultaneous contact-region trimming",
                ));
            }
        }
    }
    let mut edit = Replacement::default();
    let mut end_changes = BTreeMap::new();
    let mut corners = Vec::new();
    for old in model.vertices_of_body(body) {
        let active: Vec<_> = runs
            .iter()
            .enumerate()
            .filter_map(|(i, run)| {
                let e = model.edge(run.edge);
                (e.start == Some(old) || e.end == Some(old)).then_some(i)
            })
            .collect();
        if active.is_empty() {
            continue;
        }
        let incident: BTreeSet<_> = uses
            .iter()
            .filter(|(edge, _)| {
                let e = model.edge(**edge);
                e.start == Some(old) || e.end == Some(old)
            })
            .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
            .collect();
        if incident.len() != 3 {
            return Err(bad("line-strip end needs a regular three-face vertex"));
        }
        match active.as_slice() {
            [i] => {
                let faces = runs[*i].faces;
                let third = *incident.iter().find(|f| !faces.contains(f)).unwrap();
                // A cylindrical end wall cuts the extended strip; its end
                // curve is the blend/wall intersection, not an invented cap.
                let (a, b) = if matches!(
                    model.geom().surface(model.face(third).surface),
                    Surface::Cylinder(_)
                ) {
                    (
                        wall_foot(model, &uses, &runs[*i], faces[0], third, old, dec)?,
                        wall_foot(model, &uses, &runs[*i], faces[1], third, old, dec)?,
                    )
                } else {
                    (
                        line_plane(rail(&runs[*i], faces[0])?, plane(model, third)?, dec)?,
                        line_plane(rail(&runs[*i], faces[1])?, plane(model, third)?, dec)?,
                    )
                };
                let [a, b] = [vertex(model, a), vertex(model, b)];
                for (face, v) in faces.into_iter().zip([a, b]) {
                    runs[*i].ends.insert((old, face), v);
                }
                let surface = runs[*i].definition.surface.clone();
                let third_surface = model.geom().surface(model.face(third).surface).clone();
                let cap = connector(
                    model,
                    &surface,
                    &third_surface,
                    a,
                    b,
                    &[&runs[*i].definition],
                    dec,
                )?;
                runs[*i].caps.insert(old, cap);
                edit.gaps.insert((third, old), cap);
                for (edge, adjacent) in &uses {
                    if *edge == runs[*i].edge {
                        continue;
                    }
                    let e = model.edge(*edge);
                    if e.start != Some(old) && e.end != Some(old) {
                        continue;
                    }
                    let v = if adjacent.iter().any(|(f, _)| *f == faces[0]) {
                        a
                    } else {
                        b
                    };
                    end_changes.insert((*edge, old), v);
                }
            }
            [i, j] => {
                let common = *runs[*i]
                    .faces
                    .iter()
                    .find(|f| runs[*j].faces.contains(f))
                    .ok_or_else(|| bad("selected strips do not share a support at their end"))?;
                let a = *runs[*i].faces.iter().find(|f| **f != common).unwrap();
                let b = *runs[*j].faces.iter().find(|f| **f != common).unwrap();
                let join = tangent_join(model, &runs[*i], &runs[*j], common, [a, b], old, dec)?;
                let (p, qa) = match &join {
                    Some((_, p, q)) => (p.clone(), q.clone()),
                    None => {
                        let p = line_line(rail(&runs[*i], common)?, rail(&runs[*j], common)?, dec)?;
                        let qa = line_plane(rail(&runs[*i], a)?, plane(model, b)?, dec)?;
                        let qb = line_plane(rail(&runs[*j], b)?, plane(model, a)?, dec)?;
                        if !same(&qa, &qb, dec)? {
                            return Err(bad(
                                "unequal endpoint contacts require a corner continuation patch",
                            ));
                        }
                        (p, qa)
                    }
                };
                let p = vertex(model, p);
                let q = vertex(model, qa);
                for index in [*i, *j] {
                    for face in runs[index].faces {
                        runs[index]
                            .ends
                            .insert((old, face), if face == common { p } else { q });
                    }
                }
                let cap = match join {
                    Some((meridian, ..)) => meridian_edge(model, meridian, p, q, dec)?,
                    None => connector(
                        model,
                        &runs[*i].definition.surface,
                        &runs[*j].definition.surface,
                        p,
                        q,
                        &[&runs[*i].definition, &runs[*j].definition],
                        dec,
                    )?,
                };
                for index in [*i, *j] {
                    runs[index].caps.insert(old, cap);
                }
                for edge in uses.keys() {
                    if *edge == runs[*i].edge || *edge == runs[*j].edge {
                        continue;
                    }
                    let e = model.edge(*edge);
                    if e.start == Some(old) || e.end == Some(old) {
                        end_changes.insert((*edge, old), q);
                    }
                }
            }
            [i, j, k] => {
                let indices = [*i, *j, *k];
                let name = model
                    .vertex_name(body, old)
                    .ok_or_else(|| bad("a corner vertex needs an unambiguous origin name"))?;
                let corner = CornerDefinition::new(
                    Origin {
                        step: step.clone(),
                        role: Role::Corner { vertex: name },
                    },
                    indices.map(|index| &runs[index].definition),
                    dec,
                )?;
                let mut faces = [None; 3];
                for (strip, index) in indices.iter().enumerate() {
                    for (support, face) in corner.strips[strip].into_iter().zip(runs[*index].faces)
                    {
                        if faces[support].replace(face).is_some_and(|f| f != face) {
                            return Err(bad("one corner support carries two input faces"));
                        }
                    }
                }
                if uses.keys().any(|edge| {
                    let e = model.edge(*edge);
                    (e.start == Some(old) || e.end == Some(old))
                        && indices.iter().all(|index| runs[*index].edge != *edge)
                }) {
                    return Err(bad("an unselected edge meets the three-strip corner"));
                }
                let feet = corner.feet.clone().map(|p| vertex(model, p));
                for (strip, index) in indices.iter().enumerate() {
                    let [a, b] = corner.strips[strip];
                    let run = &mut runs[*index];
                    run.ends.insert((old, run.faces[0]), feet[a]);
                    run.ends.insert((old, run.faces[1]), feet[b]);
                    let curve = model.geom_mut().add_curve(corner.boundary(strip));
                    let cap = model.add_edge(Edge {
                        curve,
                        sense: Sense::Same,
                        start: Some(feet[a]),
                        end: Some(feet[b]),
                    });
                    runs[*index].caps.insert(old, cap);
                }
                corner_guard(model, body, &corner, &incident, dec)?;
                corners.push((old, indices, corner));
            }
            _ => {
                return Err(bad(
                    "four or more selected strips at a vertex need a general corner patch",
                ));
            }
        }
    }
    for (old, incident) in &uses {
        let e = model.edge(*old).clone();
        let start = e
            .start
            .map(|v| end_changes.get(&(*old, v)).copied().unwrap_or(v));
        let end = e
            .end
            .map(|v| end_changes.get(&(*old, v)).copied().unwrap_or(v));
        if start == e.start && end == e.end {
            continue;
        }
        let edge = model.add_edge(Edge { start, end, ..e });
        for (face, _) in incident {
            edit.edges.insert((*face, *old), edge);
        }
    }
    let mut cap_uses = BTreeMap::new();
    for run in &runs {
        let input = model.edge(run.edge).clone();
        let mut rails = Vec::with_capacity(2);
        for (i, face) in run.faces.iter().enumerate() {
            for old in [input.start, input.end].into_iter().flatten() {
                let v = run.ends[&(old, *face)];
                let p = model.geom().position(model.vertex(v).point);
                if !catallax_boolean::face_contains(model, *face, &p, dec)? {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "a contact foot leaves the input face and needs neighbour continuation",
                    ));
                }
            }
            let curve = model
                .geom_mut()
                .add_curve(run.definition.contacts[i].clone());
            let edge = model.add_edge(Edge {
                curve,
                sense: input.sense,
                start: input.start.map(|v| run.ends[&(v, *face)]),
                end: input.end.map(|v| run.ends[&(v, *face)]),
            });
            edit.edges.insert((*face, run.edge), edge);
            rails.push((edge, reverse(run.senses[i])));
        }
        let paths = if input.start.is_none() && input.end.is_none() {
            vec![vec![rails[0]], vec![rails[1]]]
        } else {
            let mut path = vec![rails[0]];
            let end = edge_blend::ends(model, rails[0].0, rails[0].1).1.unwrap();
            let old = if rails[0].1 == Sense::Same {
                input.end.unwrap()
            } else {
                input.start.unwrap()
            };
            let cap = run.caps[&old];
            let sense = if model.edge(cap).start == Some(end) {
                Sense::Same
            } else {
                Sense::Reversed
            };
            cap_uses.insert(cap, sense);
            path.push((cap, sense));
            path.push(rails[1]);
            let end = edge_blend::ends(model, rails[1].0, rails[1].1).1.unwrap();
            let old = if rails[1].1 == Sense::Same {
                input.end.unwrap()
            } else {
                input.start.unwrap()
            };
            let cap = run.caps[&old];
            let sense = if model.edge(cap).start == Some(end) {
                Sense::Same
            } else {
                Sense::Reversed
            };
            cap_uses.insert(cap, sense);
            path.push((cap, sense));
            vec![path]
        };
        let face = crate::extrude::face(
            model,
            run.definition.surface.clone(),
            if run.definition.convex {
                Sense::Same
            } else {
                Sense::Reversed
            },
            paths,
            None,
            vec![run.definition.origin.clone()],
        );
        edit.faces.push((run.faces[0], face));
    }
    for (old, indices, corner) in &corners {
        // Each quarter circle is used once by its strip and once, reversed,
        // by the sphere sector.
        let mut caps: Vec<_> = indices
            .iter()
            .map(|index| {
                let cap = runs[*index].caps[old];
                (cap, reverse(cap_uses[&cap]))
            })
            .collect();
        let mut path = vec![caps.remove(0)];
        while !caps.is_empty() {
            let end = edge_blend::ends(model, path.last().unwrap().0, path.last().unwrap().1).1;
            let next = caps
                .iter()
                .position(|(cap, sense)| edge_blend::ends(model, *cap, *sense).0 == end)
                .ok_or_else(|| bad("corner boundary arcs do not close"))?;
            path.push(caps.remove(next));
        }
        let face = crate::extrude::face(
            model,
            corner.surface.clone(),
            Sense::Same,
            vec![path],
            None,
            vec![corner.origin.clone()],
        );
        edit.faces.push((runs[indices[0]].faces[0], face));
    }
    let result = edge_blend::apply(model, body, &edit)?;
    if let Some(trimming) = trimming {
        trimming.finish(dec, None, None);
    }
    Ok((result, runs.into_iter().map(|run| run.definition).collect()))
}
