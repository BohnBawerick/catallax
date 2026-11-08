//! Rail-first elementary fillet construction over an immutable input skin.

use std::collections::{BTreeMap, BTreeSet};

use catallax_intersect::{Component, CurveSurface, Cyclic, SurfacePair};
use catallax_numbers::{
    BoundKind, Decider, Fact, FactKind, FactOrigin, Mark, Outcome, Rational, Real, Refusal, Sign,
};
use catallax_recipe::{EdgeSelection, Origin, Role, StepId};
use catallax_shapes::{
    Circle, Curve, GeomKey, Line, Point3, PointDef, PointId, RollingMitre, SplineSupport, Surface,
    Torus, TorusSheet, Vector3,
};
use catallax_topology::{BodyId, Edge, EdgeId, FaceId, Loop, Model, Sense, Vertex, VertexId};

use crate::{
    ConsumedFace,
    edge_blend::{self, Gap, Replacement},
    extrude::reverse,
    fillet::{
        FilletDefinition, UNSUPPORTED_CORNER, UNSUPPORTED_END, hull_clear, hulls_apart,
        hulls_touch, separated,
    },
    fillet_corner::CornerDefinition,
    profile::sign,
};

#[path = "fillet_collapse.rs"]
mod fillet_collapse;
#[path = "fillet_consume.rs"]
mod fillet_consume;
#[path = "fillet_envelope.rs"]
mod fillet_envelope;
#[path = "fillet_equal_crossing.rs"]
mod fillet_equal_crossing;
#[path = "fillet_inscribed.rs"]
mod fillet_inscribed;
#[path = "fillet_limit.rs"]
mod fillet_limit;
#[path = "fillet_mixed.rs"]
mod fillet_mixed;
#[path = "fillet_setback.rs"]
pub(crate) mod fillet_setback;
#[path = "fillet_tangent.rs"]
mod fillet_tangent;
#[path = "fillet_vertex.rs"]
mod fillet_vertex;
#[path = "fillet_vertex_end.rs"]
mod fillet_vertex_end;
#[path = "fillet_vertex_network.rs"]
mod fillet_vertex_network;
#[path = "fillet_vertex_smooth.rs"]
mod fillet_vertex_smooth;

#[path = "fillet_variable.rs"]
mod variable;

pub(crate) use variable::build as build_variable;

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
        EdgeSelection::Refs(refs) => resolve_refs(model, refs, uses, dec),
    }
}
fn resolve_refs(
    model: &Model,
    refs: &[catallax_recipe::EdgeRef],
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Vec<EdgeId>, Refusal> {
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

struct Run {
    edge: EdgeId,
    faces: [FaceId; 2],
    senses: [Sense; 2],
    definition: FilletDefinition,
    ends: BTreeMap<(VertexId, FaceId), VertexId>,
    /// The edges closing the strip face at each input end vertex.
    caps: BTreeMap<VertexId, Vec<EdgeId>>,
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
/// The point of a circle at its reference direction.
fn circle_point(c: &Circle) -> Point3 {
    c.center.plus(
        &c.ref_dir
            .scaled(&(&c.radius_sq / &c.ref_dir.norm_squared()).sqrt()),
    )
}
/// Proves each contact of a ring strip leaves the edge into its own face, the
/// side where the face continues away from the edge. A wrong material side
/// puts a rail beyond the face it replaces: the loops still close and the
/// shell stays oriented, so no structural check sees it.
fn rails_on_faces(model: &Model, run: &Run, dec: &mut Decider) -> Result<(), Refusal> {
    let e = model.edge(run.edge);
    let Curve::Circle(edge) = model.geom().curve(e.curve) else {
        return Ok(());
    };
    let at = circle_point(edge);
    let along = edge.normal.cross(&at.minus(&edge.center));
    for (index, face) in run.faces.into_iter().enumerate() {
        let Curve::Circle(contact) = &run.definition.contacts[index] else {
            return Err(bad("a ring strip needs circular contacts"));
        };
        let foot = circle_point(contact);
        // A face lies to the left of its coedges about its outward normal.
        let forward = if run.senses[index] == e.sense { 1 } else { -1 };
        let into = outward(model, face, &at, dec)?.cross(&along.scaled(&Real::from_i64(forward)));
        if sign(&foot.minus(&at).dot(&into), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "fillet.contact-off-face",
                "a contact rail leaves the edge away from its face, so the ball is on the wrong material side",
            ));
        }
    }
    Ok(())
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
            // The dihedral sign at one point of the ring. Both faces count: a
            // plane inside its rim is a convex boss top over an outward wall
            // and a concave hole floor over an inward one.
            let at = circle_point(c);
            let along = c.normal.cross(&at.minus(&c.center));
            let na = outward(model, *a, &at, dec)?;
            let nb = outward(model, *b, &at, dec)?;
            let value = nb.dot(&na.cross(&along)) * Real::from_i64(orientation);
            match sign(&value, dec)? {
                Sign::Negative => Ok(true),
                Sign::Positive => Ok(false),
                Sign::Zero => Err(bad("tangent faces do not define a regular sharp edge")),
            }
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
        Curve::Ellipse(c) => {
            // The end of the first axis, where the curve runs along the second.
            let at = c.center.plus(
                &c.major_dir
                    .scaled(&(&c.major_sq / &c.major_dir.norm_squared()).sqrt()),
            );
            let along = c.normal.cross(&c.major_dir);
            let na = outward(model, *a, &at, dec)?;
            let nb = outward(model, *b, &at, dec)?;
            let value = nb.dot(&na.cross(&along)) * Real::from_i64(orientation);
            match sign(&value, dec)? {
                Sign::Negative => Ok(true),
                Sign::Positive => Ok(false),
                Sign::Zero => Err(bad("tangent faces do not define a regular sharp edge")),
            }
        }
        Curve::Intersection(ring) => {
            let point = ring.witness;
            let at = model.geom().position(point);
            let [_, along] = catallax_intersect::tangent_at(model.geom(), e.curve, point, dec)?;
            let na = outward(model, *a, &at, dec)?;
            let nb = outward(model, *b, &at, dec)?;
            let value = nb.dot(&na.cross(&along)) * Real::from_i64(orientation);
            match sign(&value, dec)? {
                Sign::Negative => Ok(true),
                Sign::Positive => Ok(false),
                Sign::Zero => Err(bad("tangent faces do not define a regular sharp edge")),
            }
        }
        _ => Err(bad("elementary fillets need a line or a conic")),
    }
}
/// The exact circle of an intersection edge whose pair plan proves it
/// circular, such as a plane square to a torus's axis cutting it at an
/// irrational radius. It runs the way the edge runs at its witness, so its
/// rails replace the edge as a stored circle's would. `None` for every other
/// edge, and for cylinder pairs, whose rings keep their crossing strips.
fn exact_ring_circle(
    model: &Model,
    edge: EdgeId,
    supports: &[Surface; 2],
    dec: &mut Decider,
) -> Result<Option<Circle>, Refusal> {
    let id = model.edge(edge).curve;
    let Curve::Intersection(ring) = model.geom().curve(id) else {
        return Ok(None);
    };
    if matches!(supports, [Surface::Cylinder(_), Surface::Cylinder(_)]) {
        return Ok(None);
    }
    let catallax_intersect::CurveForm::Circle {
        center,
        normal,
        radius_sq,
        ref_dir,
    } = catallax_intersect::curve_form(model.geom(), id, dec)?
    else {
        return Ok(None);
    };
    let at = model.geom().position(ring.witness);
    // The witness must lie on the plan's circle: the edge is this branch.
    let from = at.minus(&center);
    if sign(&from.dot(&normal), dec)? != Sign::Zero
        || sign(&(&from.norm_squared() - &radius_sq), dec)? != Sign::Zero
    {
        return Ok(None);
    }
    let [_, along] = catallax_intersect::tangent_at(model.geom(), id, ring.witness, dec)?;
    let normal = match sign(&normal.cross(&from).dot(&along), dec)? {
        Sign::Positive => normal,
        Sign::Negative => normal.scaled(&Real::from_i64(-1)),
        Sign::Zero => return Ok(None),
    };
    Ok(Some(Circle {
        center,
        normal,
        radius_sq,
        ref_dir,
    }))
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
        catallax_intersect::intersect_curve_surface_framed(model.geom_mut(), line, surface, dec)?
            .value;
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
        // The ball of a sphere strip does not move.
        Surface::Sphere(s) => Ok(s.center.clone()),
        _ => Err(bad(
            "tangent strip joins need cylinder, torus or sphere strips",
        )),
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
    // Each contact must lie on its own strip's contact curve, or be the
    // point a sphere strip's contact collapses to.
    for (run, face, at) in [
        (first, common, &p),
        (second, common, &p),
        (first, outer[0], &qa),
        (second, outer[1], &qa),
    ] {
        let index = run.faces.iter().position(|f| *f == face).unwrap();
        if run.definition.point_contact() == Some(index) {
            let Curve::Circle(point) = &run.definition.contacts[index] else {
                unreachable!("a sphere strip's contacts are circles");
            };
            if !same(&point.center, at, dec)? {
                return Err(bad(
                    "a tangent join contact misses its sphere strip's point contact",
                ));
            }
            continue;
        }
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
/// Where a ring strip meets a planar end face through its axis: the contact
/// on each support in that plane, in strip face order, and the meridian
/// through them about the ball centre there. Any other end plane cuts the
/// ring's carrier in a spiric curve, which no edge carries yet.
fn ring_end(
    model: &Model,
    run: &Run,
    third: FaceId,
    old: VertexId,
    dec: &mut Decider,
) -> Result<(Point3, Point3, Circle), Refusal> {
    let end = |detail: &str| Refusal::new(UNSUPPORTED_END, detail);
    let edge = model.edge(run.edge);
    if edge.start == edge.end {
        return Err(end("an arc strip whose two ends meet at one vertex"));
    }
    let Surface::Plane(wall) = model.geom().surface(model.face(third).surface) else {
        return Err(end("an arc strip end needs a planar end face"));
    };
    let (Surface::Torus(Torus {
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
    })) = &run.definition.surface
    else {
        return Err(end("an arc strip end needs a toroidal strip"));
    };
    if sign(&wall.normal.dot(axis), dec)? != Sign::Zero
        || sign(&center.minus(&wall.origin).dot(&wall.normal), dec)? != Sign::Zero
    {
        return Err(end(
            "an arc strip ends at a meridian only on an end plane through its ring axis",
        ));
    }
    let at = model.geom().position(model.vertex(old).point);
    let d = at.minus(center);
    let out = unit(&d.minus(&axis.scaled(&(d.dot(axis) / axis.norm_squared()))));
    let centre = center.plus(&out.scaled(major));
    let [a, b] = run
        .definition
        .contacts
        .each_ref()
        .map(|contact| match contact {
            Curve::Circle(c) => Ok(c.center.plus(&out.scaled(&c.radius_sq.sqrt()))),
            _ => Err(end("an arc strip end needs circular contacts")),
        });
    let (a, b) = (a?, b?);
    let meridian = Circle {
        center: centre.clone(),
        normal: wall.normal.clone(),
        radius_sq: run.definition.radius.square(),
        ref_dir: a.minus(&centre),
    };
    Ok((a, b, meridian))
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
/// Whether a support's normal along a curve on it turns rigidly with the
/// curve: along a line, on a plane or on a cylinder along its axis; round a
/// circle, on a surface of revolution about the circle's axis, on any sphere,
/// or on a torus whose tube it is. Two such faces tangent with agreeing
/// normals at one point of the curve are tangent along all of it.
fn steady(surface: &Surface, curve: &Curve, dec: &mut Decider) -> Result<bool, Refusal> {
    let parallel = |a: &Vector3, b: &Vector3, dec: &mut Decider| {
        Ok::<_, Refusal>(sign(&a.cross(b).norm_squared(), dec)? == Sign::Zero)
    };
    match (curve, surface) {
        (Curve::Line(_), Surface::Plane(_)) | (Curve::Circle(_), Surface::Sphere(_)) => Ok(true),
        (Curve::Line(line), Surface::Cylinder(c)) => parallel(&line.dir, &c.axis, dec),
        (Curve::Circle(k), Surface::Plane(p)) => parallel(&k.normal, &p.normal, dec),
        (Curve::Circle(k), Surface::Cylinder(c)) => Ok(parallel(&k.normal, &c.axis, dec)?
            && parallel(&k.center.minus(&c.origin), &c.axis, dec)?),
        (
            Curve::Circle(k),
            Surface::Torus(Torus {
                center,
                axis,
                major,
                minor,
                ..
            })
            | Surface::TorusSheet(TorusSheet {
                center,
                axis,
                major,
                minor,
                ..
            }),
        ) => {
            let radial = k.center.minus(center);
            if parallel(&k.normal, axis, dec)? && parallel(&radial, axis, dec)? {
                return Ok(true);
            }
            // A tube circle: centred on the core circle, in the plane through
            // the axis, with the tube's radius.
            Ok(sign(&radial.dot(axis), dec)? == Sign::Zero
                && sign(&(radial.norm_squared() - major.square()), dec)? == Sign::Zero
                && sign(&k.normal.dot(axis), dec)? == Sign::Zero
                && sign(&k.normal.dot(&radial), dec)? == Sign::Zero
                && sign(&(&k.radius_sq - &minor.square()), dec)? == Sign::Zero)
        }
        _ => Ok(false),
    }
}
/// Whether an input edge is smooth along all of it, so no ball rounds it. An
/// earlier strip's contact rail, between the strip and one of its own
/// supports, is smooth by the strip's retained definition whatever its
/// carrier; other edges are smooth when their faces are tangent with
/// agreeing normals at one point and both turn rigidly along the edge.
fn smooth(
    model: &Model,
    edge: EdgeId,
    incident: &[(FaceId, Sense)],
    earlier: &[FilletDefinition],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [(a, _), (b, _)] = incident else {
        return Ok(false);
    };
    let surface = |f: &FaceId| model.geom().surface(model.face(*f).surface);
    let e = model.edge(edge);
    let curve = model.geom().curve(e.curve);
    for (blend, other) in [(a, b), (b, a)] {
        let origins = model.face_origins(*blend);
        if earlier.iter().any(|strip| {
            origins.contains(&strip.origin)
                && strip
                    .supports
                    .iter()
                    .zip(&strip.contacts)
                    .any(|(support, rail)| {
                        support.same_definition(surface(other)) && rail.same_definition(curve)
                    })
        }) {
            return Ok(true);
        }
    }
    let at = match (curve, e.start) {
        (_, Some(start)) => model.geom().position(model.vertex(start).point),
        (Curve::Line(line), None) => line.origin.clone(),
        (Curve::Circle(circle), None) => circle_point(circle),
        _ => return Ok(false),
    };
    Ok(steady(surface(a), curve, dec)?
        && steady(surface(b), curve, dec)?
        && tangent_faces(model, *a, *b, &at, dec)?)
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
    radius: &Real,
    earlier: &[FilletDefinition],
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
        let convex = convex(model, edge, &uses[&edge], dec)?;
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
                    // A ball wider than an earlier round of its convexity
                    // cannot roll onto it: the strip takes the cliff there.
                    if narrower_round(model, uses, *boundary, leave, radius, convex, earlier, dec)?
                        .is_some()
                    {
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

/// The arc of the curve shared by two supports from one vertex to another
/// whose sample lies in every sector. An arc that starts inside a sector can
/// leave it and return, so one through `avoid`, where a closed curve crosses
/// that sector's other side, is skipped. With a `face`, the arc must lie on
/// it: every arc that passes those tests is collected, and when more than one
/// does, the one whose sample the face contains is taken.
#[allow(clippy::too_many_arguments)]
fn connector(
    model: &mut Model,
    a: &Surface,
    b: &Surface,
    from: VertexId,
    to: VertexId,
    sectors: &[&FilletDefinition],
    avoid: Option<&Point3>,
    face: Option<FaceId>,
    dec: &mut Decider,
) -> Result<EdgeId, Refusal> {
    let avoid = avoid.map(|p| model.geom_mut().add_point(PointDef::Explicit(p.clone())));
    let a = model.geom_mut().add_surface(a.clone());
    let b = model.geom_mut().add_surface(b.clone());
    let meeting = catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value;
    let SurfacePair::Meets(components) = meeting else {
        return Err(bad("blend end supports have no regular common curve"));
    };
    let p = model.vertex(from).point;
    let q = model.vertex(to).point;
    let mut found = Vec::new();
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
            if let Some(avoid) = avoid
                && catallax_intersect::curve_is_closed(model.geom(), curve, dec)?
                && decided(catallax_intersect::cyclic_order(
                    model.geom(),
                    curve,
                    [start, avoid, end],
                    dec,
                ))? == Some(Cyclic::Forward)
            {
                continue;
            }
            let sample =
                catallax_intersect::arc_point(model.geom(), curve, Some(start), Some(end), dec)?
                    .value;
            let mut inside = true;
            for sector in sectors {
                inside &= sector.in_sector(&sample, dec)?;
            }
            if inside {
                if face.is_none() {
                    return Ok(model.add_edge(Edge {
                        curve,
                        sense,
                        start: Some(from),
                        end: Some(to),
                    }));
                }
                found.push((curve, sense, sample));
            }
        }
    }
    if let Some(face) = face
        && found.len() > 1
    {
        let mut kept = Vec::new();
        for (curve, sense, sample) in found {
            if catallax_boolean::face_contains(model, face, &sample, dec)? {
                kept.push((curve, sense, sample));
            }
        }
        found = kept;
    }
    match found.as_slice() {
        [] => Err(bad("no end-curve arc lies in every active contact sector")),
        [(curve, sense, _)] => Ok(model.add_edge(Edge {
            curve: *curve,
            sense: *sense,
            start: Some(from),
            end: Some(to),
        })),
        _ => Err(bad("two end-curve arcs lie on one cut face")),
    }
}

/// The planes that cut a strip at its two end vertices, each as a point and
/// a normal: the normal plane of its edge through a corner's ball centre in
/// `centres`, the one face there besides the strip's own two when that face
/// is planar, or the normal plane of its edge where another selected run
/// continues it tangentially, the meridian of their join.
fn end_cuts(
    model: &Model,
    run: &Run,
    runs: &[Run],
    uses: &Uses,
    centres: &BTreeMap<(EdgeId, VertexId), Point3>,
    dec: &mut Decider,
) -> Result<Option<[(Point3, Vector3); 2]>, Refusal> {
    let e = model.edge(run.edge);
    let mut cut_at = |v: VertexId| -> Result<Option<(Point3, Vector3)>, Refusal> {
        if let (Some(centre), Curve::Line(line)) =
            (centres.get(&(run.edge, v)), model.geom().curve(e.curve))
        {
            return Ok(Some((centre.clone(), line.dir.clone())));
        }
        let others: BTreeSet<_> = uses
            .iter()
            .filter(|(edge, _)| {
                let e = model.edge(**edge);
                e.start == Some(v) || e.end == Some(v)
            })
            .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
            .filter(|f| !run.faces.contains(f))
            .collect();
        if let [face] = others.into_iter().collect::<Vec<_>>()[..]
            && let Surface::Plane(p) = model.geom().surface(model.face(face).surface)
        {
            return Ok(Some((p.origin.clone(), p.normal.clone())));
        }
        let Curve::Line(line) = model.geom().curve(e.curve) else {
            return Ok(None);
        };
        let point = model.vertex(v).point;
        for other in runs {
            let o = model.edge(other.edge);
            if other.edge == run.edge || (o.start != Some(v) && o.end != Some(v)) {
                continue;
            }
            let [_, along] = catallax_intersect::tangent_at(model.geom(), o.curve, point, dec)?;
            if sign(&along.cross(&line.dir).norm_squared(), dec)? == Sign::Zero {
                return Ok(Some((model.geom().position(point), line.dir.clone())));
            }
        }
        Ok(None)
    };
    let (Some(start), Some(end)) = (e.start, e.end) else {
        return Ok(None);
    };
    Ok(cut_at(start)?.zip(cut_at(end)?).map(|(a, b)| [a, b]))
}
/// Points whose convex hull holds a strip's region: a straight strip between
/// its end cuts, or an arc strip over its input arc. A strip ends on its end
/// faces, tangent joins and meridians there, so only that part can meet
/// another face.
fn hull(
    model: &Model,
    run: &Run,
    runs: &[Run],
    uses: &Uses,
    centres: &BTreeMap<(EdgeId, VertexId), Point3>,
    dec: &mut Decider,
) -> Result<Option<Vec<Point3>>, Refusal> {
    let e = model.edge(run.edge);
    if let (Curve::Circle(_), Some(start), Some(end)) =
        (model.geom().curve(e.curve), e.start, e.end)
    {
        let [from, to] = if e.sense == Sense::Same {
            [start, end]
        } else {
            [end, start]
        }
        .map(|v| model.geom().position(model.vertex(v).point));
        return run.definition.arc_hull([&from, &to], dec);
    }
    Ok(match end_cuts(model, run, runs, uses, centres, dec)? {
        Some([(a, m), (b, n)]) => run
            .definition
            .cut_hull([(&a, &m), (&b, &n)], dec)?
            .map(Vec::from),
        None => None,
    })
}

/// Orders and orients the edges closing one strip end into a path from `from`.
fn chain(model: &Model, caps: &[EdgeId], from: VertexId) -> Result<Vec<(EdgeId, Sense)>, Refusal> {
    let mut left = caps.to_vec();
    let mut at = from;
    let mut path = Vec::with_capacity(left.len());
    while !left.is_empty() {
        let next = left
            .iter()
            .position(|cap| {
                let e = model.edge(*cap);
                e.start == Some(at) || e.end == Some(at)
            })
            .ok_or_else(|| bad("the edges closing a strip end do not join"))?;
        let cap = left.remove(next);
        let e = model.edge(cap);
        let (sense, to) = if e.start == Some(at) {
            (Sense::Same, e.end)
        } else {
            (Sense::Reversed, e.start)
        };
        at = to.ok_or_else(|| bad("a strip end edge has no end"))?;
        path.push((cap, sense));
    }
    Ok(path)
}

/// Two straight strips of one convexity sharing one face meet at a vertex where their outer
/// rails cross the unselected third edge at different points `hits`. Each
/// strip's region is cut by the other's outer face, and their union leaves
/// the deeper strip running past the shallower one's outer rail: that rail
/// ends where it enters the deeper carrier, from where the deeper strip is cut
/// by the shallower one's outer face down to its own hit on the third edge.
/// The shallower strip's cut by the deeper one's outer face lies inside the
/// deeper strip's band there, so nothing else shows. Returns the deeper hit,
/// the shallower strip's index in `runs`, and the point where its outer rail
/// enters the deeper carrier.
fn unequal_mitre(
    model: &mut Model,
    uses: &Uses,
    runs: [&Run; 2],
    outer: [FaceId; 2],
    hits: [Point3; 2],
    old: VertexId,
    dec: &mut Decider,
) -> Result<(Point3, usize, Point3), Refusal> {
    let corner = |detail: &str| Refusal::new(UNSUPPORTED_CORNER, detail);
    let third: Vec<_> = uses
        .iter()
        .filter(|(edge, _)| {
            let e = model.edge(**edge);
            (e.start == Some(old) || e.end == Some(old))
                && runs.iter().all(|run| run.edge != **edge)
        })
        .collect();
    let [(third, incident)] = third[..] else {
        return Err(corner(
            "two strips with unequal reach need one unselected edge between them",
        ));
    };
    let (third, incident) = (*third, incident.clone());
    if convex(model, third, &incident, dec)? != runs[0].definition.convex {
        return Err(corner(
            "strips with unequal reach across an edge of the other convexity need a general corner patch",
        ));
    }
    let at = model.geom().position(model.vertex(old).point);
    let [first, second] = hits.each_ref().map(|hit| hit.minus(&at));
    // Both rails reach the third edge on the side where it leaves the vertex.
    if sign(&first.dot(&second), dec)? != Sign::Positive {
        return Err(corner(
            "outer rails meeting the unselected edge on both sides of the vertex need a general corner patch",
        ));
    }
    let deep = match sign(&(second.norm_squared() - first.norm_squared()), dec)? {
        Sign::Positive => 1,
        Sign::Negative => 0,
        Sign::Zero => return Err(corner("unequal contacts at one distance from the vertex")),
    };
    let shallow = 1 - deep;
    let q = model
        .geom_mut()
        .add_point(PointDef::Explicit(hits[deep].clone()));
    if !on_trimmed(model, third, q, dec)? {
        return Err(Refusal::new(
            "fillet.unsupported-overflow",
            "the deeper strip at a vertex runs past the unselected edge and needs neighbour continuation",
        ));
    }
    let line = model
        .geom_mut()
        .add_curve(Curve::Line(rail(runs[shallow], outer[shallow])?.clone()));
    let carrier = model
        .geom_mut()
        .add_surface(runs[deep].definition.surface.clone());
    let CurveSurface::Meets(entries) =
        catallax_intersect::intersect_curve_surface_framed(model.geom_mut(), line, carrier, dec)?
            .value
    else {
        return Err(corner(
            "the shallower outer rail does not enter the deeper strip",
        ));
    };
    let mut entry = None;
    for hit in entries {
        let p = model.geom().position(hit.point);
        if runs[deep].definition.in_sector(&p, dec)? && entry.replace(p).is_some() {
            return Err(corner(
                "the shallower outer rail enters the deeper strip twice",
            ));
        }
    }
    let entry =
        entry.ok_or_else(|| corner("the shallower outer rail does not enter the deeper strip"))?;
    Ok((hits[deep].clone(), shallow, entry))
}

/// A strip ending where two earlier strips meet at a mitre, each strip
/// tangent to one of its faces. When their retained definitions and the
/// strip form the trihedral corner the three would have closed together,
/// that corner's sphere sector closes the end and swallows the mitre: its
/// boundary arcs trim both earlier strips, and its contact on their common
/// support is the mitre's other end.
struct Mitre {
    corner: CornerDefinition,
    /// The earlier strips' faces, in the corner's strip order after the run.
    faces: [FaceId; 2],
    /// The mitre between them, and its end away from the strip.
    edge: EdgeId,
    apex: VertexId,
    /// The index of their common support in the corner, whose contact is
    /// the apex.
    shared: usize,
}
/// The earlier strip whose face this is, by its origin and carrier. A
/// continued chain's strips share their origin.
fn retained<'a>(
    model: &Model,
    face: FaceId,
    earlier: &'a [FilletDefinition],
) -> Option<&'a FilletDefinition> {
    let origins = model.face_origins(face);
    let surface = model.geom().surface(model.face(face).surface);
    earlier
        .iter()
        .find(|strip| origins.contains(&strip.origin) && strip.surface.same_definition(surface))
}
#[allow(clippy::too_many_arguments)]
fn mitre(
    model: &Model,
    body: BodyId,
    uses: &Uses,
    run: &Run,
    old: VertexId,
    earlier: &[FilletDefinition],
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<Mitre>, Refusal> {
    let at = |edge: &EdgeId| {
        let e = model.edge(*edge);
        e.start == Some(old) || e.end == Some(old)
    };
    let others: BTreeSet<_> = uses
        .iter()
        .filter(|(edge, _)| at(edge))
        .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
        .filter(|f| !run.faces.contains(f))
        .collect();
    let [c, d] = others.into_iter().collect::<Vec<_>>()[..] else {
        return Ok(None);
    };
    let (Some(first), Some(second)) = (retained(model, c, earlier), retained(model, d, earlier))
    else {
        return Ok(None);
    };
    let Some(edge) = uses
        .iter()
        .find(|(edge, incident)| {
            at(edge)
                && incident.iter().any(|(f, _)| *f == c)
                && incident.iter().any(|(f, _)| *f == d)
        })
        .map(|(edge, _)| *edge)
    else {
        return Ok(None);
    };
    // A smooth join of two earlier strips, such as a chain's meridian, is
    // no mitre.
    if smooth(model, edge, &uses[&edge], earlier, dec)? {
        return Ok(None);
    }
    let corner_end = |detail: &str| Refusal::new(UNSUPPORTED_CORNER, detail);
    let e = model.edge(edge);
    let apex = if e.start == Some(old) { e.end } else { e.start }
        .ok_or_else(|| corner_end("a mitre between earlier strips needs two ends"))?;
    let name = model
        .vertex_name(body, old)
        .ok_or_else(|| bad("a corner vertex needs an unambiguous origin name"))?;
    let corner = CornerDefinition::new(
        Origin {
            step: step.clone(),
            role: Role::Corner { vertex: name },
        },
        [&run.definition, first, second],
        dec,
    )?;
    let shared = corner.strips[1]
        .into_iter()
        .find(|i| corner.strips[2].contains(i))
        .ok_or_else(|| corner_end("the earlier strips share no support"))?;
    let tip = model.geom().position(model.vertex(apex).point);
    if !same(&corner.sector.feet[shared], &tip, dec)? {
        return Err(corner_end(
            "the earlier strips' mitre does not end at the corner's contact on their common support",
        ));
    }
    Ok(Some(Mitre {
        corner,
        faces: [c, d],
        edge,
        apex,
        shared,
    }))
}
/// A straight strip ending where two earlier rounds of one radius, both
/// wider than its ball and of its convexity, meet at a mitre. Section 4.5.1
/// rolls the ball over the mitre to its apex along a [`RollingMitre`] band.
struct RollingCorner {
    band: RollingMitre,
    // The rounds, rail 0 on the first.
    faces: [FaceId; 2],
    edge: EdgeId,
    apex: VertexId,
    // The feet in the run's face order, and the band rail each foot starts.
    feet: [Point3; 2],
    rails: [usize; 2],
    origin: Origin,
}

/// The rational with the smallest denominator in `[low, high]`, by continued
/// fractions, or `None` past a denominator of 64 bits.
fn simplest_between(low: &Rational, high: &Rational) -> Option<Rational> {
    if low > high {
        return None;
    }
    let zero = Rational::zero();
    if low <= &zero && &zero <= high {
        return Some(zero);
    }
    if high < &zero {
        return simplest_between(&-high, &-low).map(|q| -q);
    }
    let (mut low, mut high) = (low.clone(), high.clone());
    // Each step keeps `x = whole + 1 / rest` with `rest` in the reciprocal range.
    let mut wholes = Vec::new();
    let tail = loop {
        if wholes.len() > 128 {
            return None;
        }
        let whole = Rational::from(low.floor());
        if whole == low {
            break whole;
        }
        let next = &whole + &Rational::one();
        if next <= high {
            break next;
        }
        wholes.push(whole.clone());
        let (a, b) = ((&high - &whole).recip()?, (&low - &whole).recip()?);
        (low, high) = (a, b);
    };
    let value = wholes
        .into_iter()
        .rev()
        .try_fold(tail, |rest, whole| Some(&whole + &rest.recip()?))?;
    (value.denominator().bit_length() <= 64).then_some(value)
}

/// A point with each coordinate replaced by the simple rational it provably
/// equals, where one exists. A rolling mitre's feet often lie on rational
/// faces while the band's law reaches them through square roots; a short
/// definition keeps every later proof about them short.
fn simplest_point(point: &Point3, dec: &mut Decider) -> Result<Point3, Refusal> {
    let mut simplify = |value: &Real| -> Result<Real, Refusal> {
        if value.as_rational().is_some() {
            return Ok(value.clone());
        }
        let ball = dec.enclose_at(value, 128)?;
        if let Some(candidate) = simplest_between(&ball.lo(), &ball.hi())
            && sign(&(value - &Real::from(candidate.clone())), dec)? == Sign::Zero
        {
            return Ok(Real::from(candidate));
        }
        Ok(value.clone())
    };
    Ok(Point3::new(
        simplify(&point.x)?,
        simplify(&point.y)?,
        simplify(&point.z)?,
    ))
}

#[allow(clippy::too_many_arguments)]
fn rolling_corner(
    model: &Model,
    body: BodyId,
    uses: &Uses,
    run: &Run,
    old: VertexId,
    earlier: &[FilletDefinition],
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<RollingCorner>, Refusal> {
    let ([Surface::Plane(_), Surface::Plane(_)], Curve::Line(_), Surface::Cylinder(strip)) = (
        &run.definition.supports,
        &run.definition.edge,
        &run.definition.surface,
    ) else {
        return Ok(None);
    };
    let at = |edge: &EdgeId| {
        let e = model.edge(*edge);
        e.start == Some(old) || e.end == Some(old)
    };
    let others: BTreeSet<_> = uses
        .iter()
        .filter(|(edge, _)| at(edge))
        .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
        .filter(|f| !run.faces.contains(f))
        .collect();
    let [c, d] = others.into_iter().collect::<Vec<_>>()[..] else {
        return Ok(None);
    };
    let (Some(first), Some(second)) = (retained(model, c, earlier), retained(model, d, earlier))
    else {
        return Ok(None);
    };
    let (Surface::Cylinder(one), Surface::Cylinder(two)) = (&first.surface, &second.surface) else {
        return Ok(None);
    };
    let convex = run.definition.convex;
    if first.convex != convex
        || second.convex != convex
        || sign(&(&first.radius - &second.radius), dec)? != Sign::Zero
        || sign(&(&first.radius - &run.definition.radius), dec)? != Sign::Positive
    {
        return Ok(None);
    }
    let Some(edge) = uses
        .iter()
        .find(|(edge, incident)| {
            at(edge)
                && incident.iter().any(|(f, _)| *f == c)
                && incident.iter().any(|(f, _)| *f == d)
        })
        .map(|(edge, _)| *edge)
    else {
        return Ok(None);
    };
    if smooth(model, edge, &uses[&edge], earlier, dec)? {
        return Ok(None);
    }
    let corner_end = |detail: &str| Refusal::new(UNSUPPORTED_CORNER, detail);
    let e = model.edge(edge);
    let apex = if e.start == Some(old) { e.end } else { e.start }
        .ok_or_else(|| corner_end("a mitre between earlier strips needs two ends"))?;
    let axes = [one, two].map(|c| Line {
        origin: c.origin.clone(),
        dir: c.axis.clone(),
    });
    if sign(&axes[0].dir.cross(&axes[1].dir).norm_squared(), dec)? == Sign::Zero {
        return Ok(None);
    }
    let Ok(meet) = line_line(&axes[0], &axes[1], dec) else {
        return Err(corner_end("the rounds' axes do not meet at their mitre"));
    };
    // Each axis points from the corner along its round.
    let corner = model.geom().position(model.vertex(old).point);
    let mut toward = Vec::with_capacity(2);
    for axis in &axes {
        toward.push(match sign(&meet.minus(&corner).dot(&axis.dir), dec)? {
            Sign::Positive => axis.dir.clone(),
            Sign::Negative => axis.dir.scaled(&Real::from_i64(-1)),
            Sign::Zero => return Err(corner_end("a round's axis is square to its corner")),
        });
    }
    // `y = a1 x a2` points to the apex for one order of the rounds.
    let tip = model.geom().position(model.vertex(apex).point);
    let mut found = None;
    for (order, faces) in [([0, 1], [c, d]), ([1, 0], [d, c])] {
        let band = RollingMitre::new(
            meet.clone(),
            order.map(|k| toward[k].clone()),
            first.radius.clone(),
            run.definition.radius.clone(),
            convex,
            dec,
        )?;
        if same(&band.apex(), &tip, dec)? {
            found = Some((band, faces));
            break;
        }
    }
    let (band, faces) = found.ok_or_else(|| {
        corner_end("the rolling ball's band does not close at the mitre's far end")
    })?;
    // The band's end section is the strip's own section through its axis.
    let centre = band.centre(&Rational::zero());
    let gap = centre.minus(&strip.origin);
    let off_axis = &gap.cross(&strip.axis).norm_squared() / &strip.axis.norm_squared();
    if sign(&off_axis, dec)? != Sign::Zero
        || sign(&strip.axis.cross(&band.frame()[1]).norm_squared(), dec)? != Sign::Zero
    {
        return Err(corner_end(
            "the rolling ball's band does not start on the strip's end section",
        ));
    }
    let mut starts = Vec::with_capacity(2);
    for v in [Rational::zero(), Rational::one()] {
        starts.push(simplest_point(&band.point(&Rational::zero(), &v), dec)?);
    }
    let mut feet = Vec::with_capacity(2);
    let mut rails = Vec::with_capacity(2);
    for face in run.faces {
        let plane = plane(model, face)?;
        let mut found = None;
        for (rail, p) in starts.iter().enumerate() {
            if sign(&p.minus(&plane.origin).dot(&plane.normal), dec)? == Sign::Zero {
                found = Some((rail, p.clone()));
            }
        }
        let (rail, p) = found.ok_or_else(|| corner_end("a strip face holds no band rail"))?;
        rails.push(rail);
        feet.push(p);
    }
    if rails[0] == rails[1] {
        return Err(corner_end("both strip faces hold one band rail"));
    }
    let name = model
        .vertex_name(body, old)
        .ok_or_else(|| bad("a corner vertex needs an unambiguous origin name"))?;
    Ok(Some(RollingCorner {
        band,
        faces,
        edge,
        apex,
        feet: [feet[0].clone(), feet[1].clone()],
        rails: [rails[0], rails[1]],
        origin: Origin {
            step: step.clone(),
            role: Role::Corner { vertex: name },
        },
    }))
}

/// Whether an edge is a mitre between two earlier strips at whose end
/// another selected edge closes with a [`Mitre`] corner. That corner's ball
/// is as wide as both strips and sits on both their axes, so the mitre's own
/// ball has that one centre and rounds nothing beyond the corner; if the
/// corner does not close, the selected edge refuses there.
fn swallowed_mitre(
    model: &Model,
    edge: EdgeId,
    selected: &[EdgeId],
    uses: &Uses,
    earlier: &[FilletDefinition],
) -> bool {
    let [(c, _), (d, _)] = uses[&edge][..] else {
        return false;
    };
    if retained(model, c, earlier).is_none() || retained(model, d, earlier).is_none() {
        return false;
    }
    let e = model.edge(edge);
    [e.start, e.end].into_iter().flatten().any(|v| {
        let at = |other: &EdgeId| {
            let o = model.edge(*other);
            o.start == Some(v) || o.end == Some(v)
        };
        let faces: BTreeSet<_> = uses
            .iter()
            .filter(|(other, _)| at(other))
            .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
            .collect();
        faces.len() == 4
            && selected.iter().any(|other| {
                *other != edge
                    && at(other)
                    && uses[other]
                        .iter()
                        .all(|(f, _)| *f != c && *f != d && faces.contains(f))
            })
    })
}
/// A straight strip between planes whose foot reaches, at one of its input
/// end vertices, an earlier round of its convexity narrower than its ball and
/// tangent to the foot's face. The ball cannot roll onto the round, so the
/// strip takes the cliff rule of recipe section 4.5.1: it keeps both planar
/// supports, runs on along its edge and the faces it meets cut it. The cut
/// walks from that foot across those faces along the carrier until it reaches
/// the strip's other face, at the other rail's end.
struct Cliff {
    /// The rail ends, in the run's face order.
    feet: [VertexId; 2],
    /// The cut curves in order from one foot to the other.
    caps: Vec<EdgeId>,
    /// Each cut face's loop takes its cut curve where its removed vertices were.
    gaps: Vec<(FaceId, VertexId, EdgeId)>,
    /// Edges wholly inside the removed region, with each face using them.
    removed: Vec<(FaceId, EdgeId)>,
    /// Edges leaving the removed region, from their removed end to the cut.
    ends: Vec<((EdgeId, VertexId), VertexId)>,
    /// Faces the cut crosses, which the guard leaves to the cut.
    faces: BTreeSet<FaceId>,
    /// A point on the edge's line past every cut curve along the edge, where
    /// a plane normal to the edge bounds the strip's region.
    cut: Point3,
}

/// The earlier round across a strip face's edge at a vertex that the ball
/// cannot roll onto: a retained strip of the run's convexity, narrower than
/// `radius`, tangent to that face along the edge.
#[allow(clippy::too_many_arguments)]
fn narrower_round(
    model: &Model,
    uses: &Uses,
    edge: EdgeId,
    face: FaceId,
    radius: &Real,
    convex: bool,
    earlier: &[FilletDefinition],
    dec: &mut Decider,
) -> Result<Option<FaceId>, Refusal> {
    let Some(round) = uses[&edge].iter().map(|(f, _)| *f).find(|f| *f != face) else {
        return Ok(None);
    };
    let Some(strip) = retained(model, round, earlier) else {
        return Ok(None);
    };
    if strip.convex != convex || sign(&(radius - &strip.radius), dec)? != Sign::Positive {
        return Ok(None);
    }
    Ok(smooth(model, edge, &uses[&edge], earlier, dec)?.then_some(round))
}

/// The cliff end of `run` at `old`, or `None` when neither foot reaches a
/// narrower round there.
fn cliff(
    model: &mut Model,
    uses: &Uses,
    run: &Run,
    old: VertexId,
    earlier: &[FilletDefinition],
    dec: &mut Decider,
) -> Result<Option<Cliff>, Refusal> {
    let ([Surface::Plane(_), Surface::Plane(_)], Curve::Line(line), Surface::Cylinder(cylinder)) = (
        &run.definition.supports,
        &run.definition.edge,
        &run.definition.surface,
    ) else {
        return Ok(None);
    };
    let (line, cylinder) = (line.clone(), cylinder.clone());
    let at = |model: &Model, edge: &EdgeId| {
        let e = model.edge(*edge);
        e.start == Some(old) || e.end == Some(old)
    };
    // Each run face's one other edge at the vertex, and whether the round
    // across it is a cliff.
    let mut sides = Vec::with_capacity(2);
    for face in run.faces {
        let edges: Vec<_> = uses
            .iter()
            .filter(|(edge, incident)| {
                **edge != run.edge && at(model, edge) && incident.iter().any(|(f, _)| *f == face)
            })
            .map(|(edge, _)| *edge)
            .collect();
        let [edge] = edges[..] else {
            return Ok(None);
        };
        let round = narrower_round(
            model,
            uses,
            edge,
            face,
            &run.definition.radius,
            run.definition.convex,
            earlier,
            dec,
        )?;
        sides.push((edge, round));
    }
    let start = match (sides[0].1, sides[1].1) {
        (Some(_), _) => 0,
        (None, Some(_)) => 1,
        (None, None) => return Ok(None),
    };
    let other = 1 - start;
    let cliff_end = |detail: &str| Refusal::new(UNSUPPORTED_CORNER, detail);
    // A cliff foot is where its rail meets the round's straight rail.
    let mut foot = |model: &mut Model, k: usize| -> Result<(VertexId, Point3), Refusal> {
        let edge = sides[k].0;
        let Curve::Line(boundary) = model.geom().curve(model.edge(edge).curve).clone() else {
            return Err(cliff_end(
                "a round's rail on a planar strip face must be straight",
            ));
        };
        let p = line_line(rail(run, run.faces[k])?, &boundary, dec)?;
        let id = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
        // line_line already proved incidence on its second line.
        let curve = model.edge(edge).curve;
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(id),
            b: GeomKey::Curve(curve),
            at: None,
            origin: FactOrigin::Construction("fillet.cliff-foot"),
        });
        if !on_trimmed(model, edge, id, dec)? {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a contact rail passes the round's rail outside its edge and needs neighbour continuation",
            ));
        }
        Ok((model.add_vertex(Vertex { point: id }), p))
    };
    let first = foot(model, start)?;
    let known = match sides[other].1 {
        Some(_) => Some(foot(model, other)?),
        None => None,
    };
    let carrier = model.geom_mut().add_surface(run.definition.surface.clone());
    let [rail_start, rail_other] = [start, other].map(|k| rail(run, run.faces[k]).cloned());
    let (_rail_start, rail_other) = (rail_start?, rail_other?);
    // Walk the cut from the first foot across the faces it meets.
    let mut face = sides[start].1.expect("the starting side is a cliff");
    let mut entry = sides[start].0;
    let (mut from, mut from_p) = first.clone();
    let mut crossings = vec![(sides[start].0, first.0)];
    if let Some((v, _)) = &known {
        crossings.push((sides[other].0, *v));
    }
    let mut caps = Vec::new();
    let mut visited: Vec<FaceId> = Vec::new();
    let last = loop {
        if visited.contains(&face) || visited.len() == 8 {
            return Err(cliff_end("a cliff cut returns to a face it crossed"));
        }
        let mut found: Option<(Point3, EdgeId, bool)> = None;
        let boundary: Vec<_> = uses
            .iter()
            .filter(|(edge, incident)| **edge != entry && incident.iter().any(|(f, _)| *f == face))
            .map(|(edge, incident)| (*edge, incident.iter().any(|(f, _)| *f == run.faces[other])))
            .collect();
        for (edge, closes) in boundary {
            let curve = model.edge(edge).curve;
            let points = match (closes, model.geom().curve(curve)) {
                (true, Curve::Line(boundary)) => {
                    let boundary = boundary.clone();
                    if sign(&rail_other.dir.cross(&boundary.dir).norm_squared(), dec)? == Sign::Zero
                    {
                        continue;
                    }
                    vec![line_line(&rail_other, &boundary, dec)?]
                }
                _ => {
                    match catallax_intersect::intersect_curve_surface_framed(
                        model.geom_mut(),
                        curve,
                        carrier,
                        dec,
                    )?
                    .value
                    {
                        CurveSurface::Meets(hits) => hits
                            .into_iter()
                            .map(|hit| model.geom().position(hit.point))
                            .collect(),
                        _ => Vec::new(),
                    }
                }
            };
            for p in points {
                if same(&p, &from_p, dec)? || !run.definition.in_sector(&p, dec)? {
                    continue;
                }
                let id = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
                if !on_trimmed(model, edge, id, dec)? {
                    continue;
                }
                if found.replace((p, edge, closes)).is_some() {
                    return Err(cliff_end("a cliff cut leaves a face at two points"));
                }
            }
        }
        let (p, edge, closes) =
            found.ok_or_else(|| cliff_end("a cliff cut does not leave a face it enters"))?;
        let to = match (&known, closes) {
            (Some((v, q)), true) => {
                if !same(q, &p, dec)? {
                    return Err(cliff_end("a cliff cut misses the other cliff foot"));
                }
                *v
            }
            _ => vertex(model, p.clone()),
        };
        let surface = model.geom().surface(model.face(face).surface).clone();
        // A carrier tangent to a round meets it on a curve with a node at
        // the foot, whose other arc to the exit runs through the other foot.
        let avoid = match (&known, caps.is_empty()) {
            (Some((_, q)), true) => Some(q.clone()),
            (Some(_), false) if closes => Some(first.1.clone()),
            _ => None,
        };
        caps.push(connector(
            model,
            &run.definition.surface,
            &surface,
            from,
            to,
            &[&run.definition],
            avoid.as_ref(),
            Some(face),
            dec,
        )?);
        visited.push(face);
        if !crossings.iter().any(|(e, _)| *e == edge) {
            crossings.push((edge, to));
        }
        if closes {
            break (to, p);
        }
        face = uses[&edge]
            .iter()
            .map(|(f, _)| *f)
            .find(|f| *f != face)
            .ok_or_else(|| bad("a cut edge has one face"))?;
        entry = edge;
        from = to;
        from_p = p;
    };
    let mut feet = [first.0; 2];
    feet[other] = last.0;
    // The removed region: past the edge's end, inside the strip's sector and
    // outside its carrier.
    let toward = {
        let o = model.edge(run.edge);
        let (Some(a), Some(b)) = (o.start, o.end) else {
            return Ok(None);
        };
        let [a, b] = [a, b].map(|v| model.geom().position(model.vertex(v).point));
        let own = model.geom().position(model.vertex(old).point);
        let far = if same(&a, &own, dec)? { b } else { a };
        own.minus(&far)
    };
    let own = model.geom().position(model.vertex(old).point);
    let inside = |p: &Point3, dec: &mut Decider| -> Result<bool, Refusal> {
        if sign(&p.minus(&own).dot(&toward), dec)? == Sign::Negative
            || !run.definition.in_sector(p, dec)?
        {
            return Ok(false);
        }
        let d = p.minus(&cylinder.origin);
        let radial = d.minus(
            &cylinder
                .axis
                .scaled(&(d.dot(&cylinder.axis) / cylinder.axis.norm_squared())),
        );
        Ok(sign(&(&radial.norm_squared() - &cylinder.radius_sq), dec)? == Sign::Positive)
    };
    let mut affected: BTreeSet<FaceId> = visited.iter().copied().collect();
    affected.extend(run.faces);
    affected.extend(
        uses.iter()
            .filter(|(edge, _)| at(model, edge))
            .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f)),
    );
    let mut gone = BTreeSet::new();
    let mut removed = Vec::new();
    let mut ends = Vec::new();
    for (edge, incident) in uses {
        if *edge == run.edge || !incident.iter().any(|(f, _)| affected.contains(f)) {
            continue;
        }
        let e = model.edge(*edge).clone();
        let mut out = Vec::with_capacity(2);
        for v in [e.start, e.end].into_iter().flatten() {
            let p = model.geom().position(model.vertex(v).point);
            if inside(&p, dec)? {
                out.push(v);
            }
        }
        let crossing = crossings.iter().find(|(c, _)| c == edge).map(|(_, v)| *v);
        match (out.len(), crossing) {
            (0, None) => {}
            (2, None) => {
                for (f, _) in incident {
                    if !affected.contains(f) {
                        return Err(cliff_end(
                            "a cliff removes an edge of a face it does not cut",
                        ));
                    }
                    removed.push((*f, *edge));
                }
                gone.extend(out);
            }
            (1, Some(cut)) => {
                ends.push(((*edge, out[0]), cut));
                gone.insert(out[0]);
            }
            _ => return Err(cliff_end("a cliff cut crosses an edge it does not end")),
        }
    }
    let mut gaps = Vec::new();
    for (face, cap) in visited.iter().zip(&caps) {
        for (edge, incident) in uses {
            if !incident.iter().any(|(f, _)| f == face) {
                continue;
            }
            let e = model.edge(*edge);
            for v in [e.start, e.end].into_iter().flatten() {
                if gone.contains(&v) && !gaps.iter().any(|(f, w, _)| f == face && *w == v) {
                    gaps.push((*face, v, *cap));
                }
            }
        }
    }
    // Bound how far the cut curves reach along the edge from certified
    // enclosures of their arcs, and cut the strip's region past that. The
    // ends' own reach is a floor; boxes reaching no more than a hair past it
    // are accepted, and others are split until they do.
    let width = Rational::parse("1e-30").expect("width");
    let dir = [
        dec.enclose(&toward.x, &width)?,
        dec.enclose(&toward.y, &width)?,
        dec.enclose(&toward.z, &width)?,
    ];
    let mut floor: Option<Rational> = None;
    for (_, v) in &crossings {
        let p = model.geom().position(model.vertex(*v).point);
        let value = dec
            .enclose(&p.minus(&Point3::origin()).dot(&toward), &width)?
            .hi();
        if floor.as_ref().is_none_or(|f| value > *f) {
            floor = Some(value);
        }
    }
    let floor = floor.expect("a cliff has crossings");
    let slack = &floor + &Rational::parse("1/1000").expect("slack");
    let precision = dir[0].precision();
    let mut reach = floor.clone();
    for cap in &caps {
        let e = model.edge(*cap).clone();
        let [a, b] = [e.start, e.end].map(|v| v.map(|v| model.vertex(v).point));
        let (a, b) = if e.sense == Sense::Same {
            (a, b)
        } else {
            (b, a)
        };
        let traced = catallax_intersect::trace(model.geom(), e.curve, a, b, dec)?.value;
        for piece in traced.pieces() {
            let mut open = vec![(Rational::zero(), Rational::one(), 0u32)];
            while let Some((low, high, depth)) = open.pop() {
                let t = catallax_numbers::Ball::from_rational(&low, precision)
                    .hull(&catallax_numbers::Ball::from_rational(&high, precision));
                let boxed = piece.enclose(&t, dec)?;
                let hi = (&(&(&boxed[0] * &dir[0]) + &(&boxed[1] * &dir[1]))
                    + &(&boxed[2] * &dir[2]))
                    .hi();
                if hi <= slack || depth == 16 {
                    if hi > reach {
                        reach = hi;
                    }
                    continue;
                }
                let middle = &(&low + &high) / &Rational::from_i64(2);
                open.push((low, middle.clone(), depth + 1));
                open.push((middle, high, depth + 1));
            }
        }
    }
    // The edge point whose projection along `toward` is `reach`.
    let base = line.origin.minus(&Point3::origin()).dot(&toward);
    let cut = line
        .origin
        .plus(&toward.scaled(&(&(&Real::from(reach) - &base) / &toward.norm_squared())));
    let mut faces: BTreeSet<FaceId> = visited.into_iter().collect();
    faces.extend(affected);
    Ok(Some(Cliff {
        feet,
        caps,
        gaps,
        removed,
        ends,
        faces,
        cut,
    }))
}

/// Chains oriented boundary arcs into one closed loop.
fn close(model: &Model, mut arcs: Vec<(EdgeId, Sense)>) -> Result<Vec<(EdgeId, Sense)>, Refusal> {
    let mut path = vec![arcs.remove(0)];
    while !arcs.is_empty() {
        let end = edge_blend::ends(model, path.last().unwrap().0, path.last().unwrap().1).1;
        let next = arcs
            .iter()
            .position(|(arc, sense)| edge_blend::ends(model, *arc, *sense).0 == end)
            .ok_or_else(|| bad("corner boundary arcs do not close"))?;
        path.push(arcs.remove(next));
    }
    Ok(path)
}

fn guard(
    model: &mut Model,
    body: BodyId,
    run: &Run,
    hull: Option<&[Point3]>,
    uses: &Uses,
    cut: &BTreeSet<FaceId>,
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
        if run.faces.contains(&face) || neighbours.contains(&face) || cut.contains(&face) {
            continue;
        }
        let other = model.geom().surface(model.face(face).surface).clone();
        if separated(&run.definition.surface, &other, dec)? {
            continue;
        }
        if let Some(hull) = hull {
            if hull_clear(hull, &other, dec)? {
                continue;
            }
            if let Some(corners) = polygon(model, face)
                && hulls_apart(hull, &corners, dec)?
            {
                continue;
            }
        }
        // Another piece of a curved strip's support touches the carrier only
        // along the contact rail, and that rail's feet lie in the run's face.
        let curved = run
            .definition
            .supports
            .iter()
            .any(|support| matches!(support, Surface::Cylinder(_) | Surface::Cone(_)));
        if curved
            && matches!(
                run.definition.edge,
                Curve::Line(_) | Curve::Ellipse(_) | Curve::Intersection(_)
            )
            && run
                .definition
                .supports
                .iter()
                .any(|support| support.same_definition(&other))
        {
            continue;
        }
        // A band has no surface pair solver. Its hull proves it off a plane,
        // and nothing yet proves it off any other surface.
        if let Surface::BSpline(band) = &run.definition.surface {
            let clear = match (band, &other) {
                (SplineSupport::RollingPipe(pipe), Surface::Plane(plane)) => {
                    Some(pipe.clear_of(plane, dec)?)
                }
                (SplineSupport::RollingSaddle(saddle), Surface::Plane(plane)) => {
                    Some(saddle.clear_of(plane, dec)?)
                }
                (SplineSupport::RollingSphere(band), Surface::Plane(plane)) => {
                    Some(band.clear_of(plane, dec)?)
                }
                (
                    SplineSupport::RollingPipe(_)
                    | SplineSupport::RollingSaddle(_)
                    | SplineSupport::RollingSphere(_),
                    _,
                ) => Some(false),
                _ => None,
            };
            match clear {
                Some(true) => continue,
                Some(false) => {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "a face reaches the hull of the blend band and needs trimmed continuation",
                    ));
                }
                None => {}
            }
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

/// The corners of a plane face bounded only by straight edges, whose convex
/// hull therefore holds the whole face. Corners may repeat.
fn polygon(model: &Model, face: FaceId) -> Option<Vec<Point3>> {
    if !matches!(
        model.geom().surface(model.face(face).surface),
        Surface::Plane(_)
    ) {
        return None;
    }
    let mut corners = Vec::new();
    for lp in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return None;
        };
        for coedge in coedges {
            let edge = model.edge(model.coedge(*coedge).edge);
            if !matches!(model.geom().curve(edge.curve), Curve::Line(_)) {
                return None;
            }
            // Coedge senses vary, so take both ends: every corner then appears.
            for end in [edge.start?, edge.end?] {
                corners.push(model.geom().position(model.vertex(end).point));
            }
        }
    }
    Some(corners)
}

fn corner_guard(
    model: &mut Model,
    body: BodyId,
    corner: &CornerDefinition,
    incident: &BTreeSet<FaceId>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    surface_guard(model, body, &corner.sector.surface, incident, dec)
}

/// Proves every face but the incident ones clear of a corner's support.
fn surface_guard(
    model: &mut Model,
    body: BodyId,
    surface: &Surface,
    incident: &BTreeSet<FaceId>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for face in model.faces_of_body(body) {
        if incident.contains(&face) {
            continue;
        }
        let other = model.geom().surface(model.face(face).surface).clone();
        if separated(surface, &other, dec)? {
            continue;
        }
        let a = model.geom_mut().add_surface(surface.clone());
        let b = model.geom_mut().add_surface(other);
        if !matches!(
            catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value,
            SurfacePair::Disjoint
        ) {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a non-incident face meets the corner's support; its trimmed continuation must be resolved before publishing",
            ));
        }
    }
    Ok(())
}

/// The sense in which a face's loop uses an edge.
fn edge_sense(model: &Model, face: FaceId, edge: EdgeId) -> Result<Sense, Refusal> {
    model
        .face(face)
        .loops
        .iter()
        .find_map(|lp| match model.loop_(*lp) {
            Loop::Coedges(cs) => cs
                .iter()
                .map(|c| model.coedge(*c))
                .find(|c| c.edge == edge)
                .map(|c| c.sense),
            Loop::Vertex(_) => None,
        })
        .ok_or_else(|| bad("an earlier strip does not use its mitre"))
}

type Strip = (EdgeId, [FaceId; 2], FilletDefinition);

/// Whether two runs' input edges share a vertex.
fn adjacent(model: &Model, a: &Run, b: &Run) -> bool {
    let [a, b] = [a, b].map(|run| model.edge(run.edge));
    [a.start, a.end]
        .into_iter()
        .flatten()
        .any(|v| b.start == Some(v) || b.end == Some(v))
}

/// In the active positive cone, write d=sum(alpha_i*n_i), alpha_i>=0.
/// Input halfspaces give n_i.d<=r=n_i.h, with h=vertex-centre, hence
/// d.(h-d)>=0. The whole cone cell lies in the ball with diameter [0,h],
/// not just on the radius-r sphere. Its enclosing box also bounds removal.
fn corner_cell_box(centre: &Point3, vertex: &Point3) -> [Point3; 8] {
    let delta = vertex.minus(centre);
    let half = Real::one() / Real::from_i64(2);
    let middle = centre.plus(&delta.scaled(&half));
    let reach = &delta.norm_squared().sqrt() * &half;
    std::array::from_fn(|index| {
        let offset = |axis: usize| {
            if index & (1_usize << axis) == 0_usize {
                -&reach
            } else {
                reach.clone()
            }
        };
        Point3::new(
            &middle.x + &offset(0),
            &middle.y + &offset(1),
            &middle.z + &offset(2),
        )
    })
}

enum RegionBound {
    Sheet {
        /// [virtual edge, A rail, C rail] on each cut; hull[3] is the input vertex.
        hull: [Point3; 6],
        centre: Point3,
        points: std::cell::OnceCell<[Point3; 14]>,
    },
    Network {
        vertex: Point3,
        centres: [Point3; 2],
        points: std::cell::OnceCell<[[Point3; 8]; 2]>,
    },
}

struct Continuation {
    definition: FilletDefinition,
    region: RegionBound,
}

impl Continuation {
    fn region_hull(&self) -> &[Point3] {
        match &self.region {
            RegionBound::Sheet {
                hull,
                centre,
                points,
            } => points.get_or_init(|| {
                let cell = corner_cell_box(centre, &hull[3]);
                let mut points = hull.iter().cloned().chain(cell);
                std::array::from_fn(|_| points.next().unwrap())
            }),
            RegionBound::Network {
                vertex,
                centres,
                points,
            } => points
                .get_or_init(|| {
                    centres
                        .each_ref()
                        .map(|centre| corner_cell_box(centre, vertex))
                })
                .as_flattened(),
        }
    }

    fn clear_of(
        &self,
        own: &[usize],
        hulls: &[Option<Vec<Point3>>],
        others: &[Self],
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        for (index, hull) in hulls.iter().enumerate() {
            if own.contains(&index) {
                continue;
            }
            if matches!(hull, Some(other) if hulls_apart(self.region_hull(), other, dec)?) {
                continue;
            }
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a junction's separation from another blend is unproved",
            ));
        }
        for other in others {
            if !hulls_apart(self.region_hull(), other.region_hull(), dec)? {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "the junction affected regions are not proved apart",
                ));
            }
        }
        Ok(())
    }
}

fn corner_cap_use(
    model: &Model,
    run: &Run,
    old: VertexId,
    cap: EdgeId,
) -> Result<(EdgeId, Sense), Refusal> {
    let input = model.edge(run.edge);
    let at_first = if reverse(run.senses[0]) == Sense::Same {
        input.end
    } else {
        input.start
    };
    let from = run.ends[&(old, run.faces[if at_first == Some(old) { 0 } else { 1 }])];
    let edge = model.edge(cap);
    let sense = if edge.start == Some(from) {
        Sense::Reversed
    } else if edge.end == Some(from) {
        Sense::Same
    } else {
        return Err(bad("a corner cap misses its strip's oriented endpoint"));
    };
    Ok((cap, sense))
}

/// Rounds the selected edges of `body`. `earlier` holds the retained
/// definitions of earlier strips, which tell their smooth rails apart from
/// sharp edges whatever their carriers. Input faces the blends remove
/// completely are added to `consumed`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    selection: &EdgeSelection,
    radius: &Real,
    roll_over: bool,
    closure: catallax_recipe::FilletCorner,
    step: &StepId,
    earlier: &[FilletDefinition],
    consumed: &mut Vec<ConsumedFace>,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<FilletDefinition>), Refusal> {
    if sign(radius, dec)? != Sign::Positive {
        return Err(Refusal::new(
            "fillet.non-positive-radius",
            "fillet radius must be positive",
        ));
    }
    catallax_numbers::prefer_reduced_angle_series();
    let selected = select(model, body, selection, roll_over, earlier, dec)?;
    if roll_over
        && let Some(result) = fillet_mixed::build(
            model, body, &selected, radius, closure, step, earlier, consumed, dec,
        )?
    {
        return Ok(result);
    }
    strips_of(
        model, body, selected, radius, roll_over, closure, step, earlier, consumed, dec,
    )
}

/// The input edges a selection rounds: every sharp edge for `"all"`, less
/// the mitres an earlier roll-over already swallowed, or the named edges.
fn select(
    model: &Model,
    body: BodyId,
    selection: &EdgeSelection,
    roll_over: bool,
    earlier: &[FilletDefinition],
    dec: &mut Decider,
) -> Result<Vec<EdgeId>, Refusal> {
    let uses = uses(model, body);
    let mut selected = resolve(model, selection, &uses, dec)?;
    // Every edge leaves out the smooth ones, which no ball rounds; a smooth
    // edge chosen by name still refuses below.
    if matches!(selection, EdgeSelection::All) {
        let mut sharp = Vec::with_capacity(selected.len());
        for edge in selected {
            if !smooth(model, edge, &uses[&edge], earlier, dec)? {
                sharp.push(edge);
            }
        }
        selected = sharp;
        if roll_over {
            let mitres: BTreeSet<_> = selected
                .iter()
                .copied()
                .filter(|edge| swallowed_mitre(model, *edge, &selected, &uses, earlier))
                .collect();
            selected.retain(|edge| !mitres.contains(edge));
        }
    }
    Ok(selected)
}

/// Rounds the given input edges of `body` in one pass, explaining a refusal
/// by the radius where a proof applies.
#[allow(clippy::too_many_arguments)]
fn strips_of(
    model: &mut Model,
    body: BodyId,
    selected: Vec<EdgeId>,
    radius: &Real,
    roll_over: bool,
    closure: catallax_recipe::FilletCorner,
    step: &StepId,
    earlier: &[FilletDefinition],
    consumed: &mut Vec<ConsumedFace>,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<FilletDefinition>), Refusal> {
    let mut strips = Vec::new();
    let mut unfinished = None;
    build_strips(
        model,
        body,
        selected,
        radius,
        roll_over,
        closure,
        step,
        earlier,
        &mut strips,
        &mut unfinished,
        consumed,
        dec,
    )
    .map_err(|refusal| {
        let refusal = match unfinished {
            Some(code) if refusal.code == "fillet.unsupported-topology" => {
                Refusal::new(code, refusal.detail)
            }
            _ => refusal,
        };
        if strips.is_empty() {
            return refusal;
        }
        let uses = uses(model, body);
        fillet_limit::width(model, &uses, &strips, refusal, dec)
    })
}

/// Builds the strips, recording each one whose definition it reached so a
/// refusal can be explained by the radius, and the refusal code that names
/// the strip end or junction it is closing, if any.
#[allow(clippy::too_many_arguments)]
fn build_strips(
    model: &mut Model,
    body: BodyId,
    selected: Vec<EdgeId>,
    radius: &Real,
    roll_over: bool,
    closure: catallax_recipe::FilletCorner,
    step: &StepId,
    earlier: &[FilletDefinition],
    strips: &mut Vec<Strip>,
    unfinished: &mut Option<&'static str>,
    consumed: &mut Vec<ConsumedFace>,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<FilletDefinition>), Refusal> {
    let uses = uses(model, body);
    if selected.is_empty() {
        return Ok((body, Vec::new()));
    }
    if let Some(result) =
        fillet_equal_crossing::build(model, body, &selected, &uses, radius, step, dec)?
    {
        return Ok(result);
    }
    if roll_over
        && let Some((result, definitions, history)) =
            fillet_inscribed::build(model, body, &selected, &uses, radius, step, dec)?
    {
        consumed.extend(history);
        return Ok((result, definitions));
    }
    let selected = if roll_over {
        continue_tangent(model, &uses, selected, radius, earlier, dec)?
    } else {
        selected.into_iter().map(|edge| (edge, edge)).collect()
    };
    let single_run = selected.len() == 1;
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
        let origin = Origin {
            step: step.clone(),
            role: Role::Blend { edge: name },
        };
        let supports = faces.map(|f| model.geom().surface(model.face(f).surface).clone());
        let same_sense = faces.map(|f| model.face(f).sense == Sense::Same);
        let curve = match exact_ring_circle(model, edge, &supports, dec)? {
            Some(circle) => Curve::Circle(circle),
            None => model.geom().curve(model.edge(edge).curve).clone(),
        };
        let convex = convex(model, edge, incident, dec)?;
        let given = (supports.clone(), curve.clone());
        let terminal = if roll_over
            && single_run
            && model.edge(edge).start.is_none()
            && model.edge(edge).end.is_none()
        {
            fillet_collapse::terminal(model, body, edge, faces, radius, convex, dec)?
        } else {
            None
        };
        let (contact_radius, region) = match terminal {
            Some(fillet_collapse::Limit {
                radius: physical,
                region,
            }) => (Some(physical), Some(region)),
            None => (None, None),
        };
        let mut definition = if let Curve::Intersection(ring) = &curve {
            let input = model.edge(edge);
            if input.start.is_some() || input.end.is_some() {
                return Err(bad("a curved intersection fillet needs one closed ring"));
            }
            let witness = model.geom().position(ring.witness);
            FilletDefinition::ring(
                origin,
                supports,
                same_sense,
                curve,
                &witness,
                radius.clone(),
                convex,
                dec,
            )
        } else {
            FilletDefinition::new(
                origin,
                supports,
                same_sense,
                curve,
                contact_radius.unwrap_or_else(|| radius.clone()),
                convex,
                dec,
            )
        }
        .map_err(|refusal| {
            let refusal = fillet_equal_crossing::corner_refusal(
                model, edge, faces, radius, convex, refusal, dec,
            );
            fillet_limit::contact(refusal, &given.0, same_sense, &given.1, radius, convex, dec)
        })?;
        if let Some(region) = region {
            if region.beyond {
                definition.requested_radius = Some(radius.clone());
            }
            let fillet_collapse::Collapse {
                body: result,
                consumed: removed,
            } = fillet_collapse::build(model, body, region, &definition, dec)?;
            let mut history = removed.map(|face| {
                face.map(|face| crate::ConsumedFace {
                    step: step.clone(),
                    source: model.face_origins(face).to_vec(),
                    by: vec![definition.origin.clone()],
                })
            });
            if let [Some(first), Some(second)] = &history {
                match first.source.cmp(&second.source) {
                    std::cmp::Ordering::Less => {}
                    std::cmp::Ordering::Equal => history[1] = None,
                    std::cmp::Ordering::Greater => history.swap(0, 1),
                }
            }
            consumed.extend(history.into_iter().flatten());
            return Ok((result, vec![definition]));
        }
        if definition.point_contact().is_some() && model.edge(edge).start.is_none() {
            let refusal = Refusal::new(
                "fillet.unsupported-contact",
                "a ball as wide as the ring's cylinder consumes the face inside the ring",
            );
            return Err(fillet_limit::contact(
                refusal, &given.0, same_sense, &given.1, radius, convex, dec,
            ));
        }
        strips.push((edge, faces, definition.clone()));
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
                run.definition.edge,
                Curve::Line(_) | Curve::Ellipse(_) | Curve::Intersection(_)
            ) && run
                .definition
                .supports
                .iter()
                .any(|support| matches!(support, Surface::Cylinder(_) | Surface::Cone(_)))
        })
        .then(|| catallax_numbers::OperationSpan::begin("brep", step, "fillet.trimming", dec));
    // A single straight run whose foot reaches a convex boundary edge before
    // the requested radius ends at that event.
    if roll_over
        && let [run] = runs.as_slice()
        && let Some((result, history, definition)) =
            fillet_consume::terminal(model, body, run, &uses, dec)?
    {
        consumed.extend(history);
        return Ok((result, vec![definition]));
    }
    for run in &runs {
        rails_on_faces(model, run, dec)?;
    }
    // Roll-over closes a strip ending at an earlier mitre with the corner the
    // three strips share, a higher-contact stratum of version 4.
    // A strip whose ball is wider than both rounds of such a mitre cannot
    // roll onto them and takes the cliff rule instead.
    let mut mitres = BTreeMap::new();
    let mut rollings = BTreeMap::new();
    let mut cliffs = BTreeMap::new();
    if roll_over {
        for (i, run) in runs.iter().enumerate() {
            let e = model.edge(run.edge);
            for old in [e.start, e.end].into_iter().flatten() {
                if let Some(found) = cliff(model, &uses, run, old, earlier, dec)? {
                    cliffs.insert((i, old), found);
                } else if let Some(found) =
                    rolling_corner(model, body, &uses, run, old, earlier, step, dec)?
                {
                    rollings.insert((i, old), found);
                } else if let Some(found) = mitre(model, body, &uses, run, old, earlier, step, dec)?
                {
                    mitres.insert((i, old), found);
                }
            }
        }
    }
    let centres = mitres
        .iter()
        .map(|((i, old), m)| ((runs[*i].edge, *old), m.corner.centre().clone()))
        .chain(
            cliffs
                .iter()
                .map(|((i, old), c)| ((runs[*i].edge, *old), c.cut.clone())),
        )
        .chain(
            rollings
                .iter()
                .map(|((i, old), r): (&(usize, VertexId), &RollingCorner)| {
                    ((runs[*i].edge, *old), r.band.centre(&Rational::zero()))
                }),
        )
        .collect();
    let hulls = runs
        .iter()
        .map(|run| hull(model, run, &runs, &uses, &centres, dec))
        .collect::<Result<Vec<_>, _>>()?;
    for (index, (run, hull)) in runs.iter().zip(&hulls).enumerate() {
        // A cliff end's cut faces are trimmed by its cut.
        let cut: BTreeSet<_> = cliffs
            .iter()
            .filter(|((i, _), _)| *i == index)
            .flat_map(|(_, c)| c.faces.iter().copied())
            .collect();
        if let Err(refusal) = guard(model, body, run, hull.as_deref(), &uses, &cut, dec) {
            // A lone closed ring may clear every countersink mouth, touch them
            // at single points, or continue over them; every other overflow,
            // and any ring outside that stratum, still refuses.
            let continuable = refusal.code == "fillet.unsupported-overflow"
                || refusal.code.starts_with("intersect.");
            let closed_ring = matches!(run.definition.edge, Curve::Intersection(_));
            if !continuable || runs.len() != 1 || closed_ring {
                return Err(refusal);
            }
            if roll_over
                && let Some((result, history)) =
                    fillet_consume::build(model, body, run, &uses, dec)?
            {
                consumed.extend(history);
                return Ok((result, vec![run.definition.clone()]));
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
        for (j, b) in runs.iter().enumerate().skip(i + 1) {
            if adjacent(model, a, b) {
                continue;
            }
            // Two strips joined at both ends of a sphere strip meet at its
            // point contact, where all three share a vertex.
            let touch = runs.iter().find_map(|run| {
                let point = run.definition.point_contact()?;
                let Curve::Circle(circle) = &run.definition.contacts[point] else {
                    return None;
                };
                (adjacent(model, run, a) && adjacent(model, run, b)).then_some(&circle.center)
            });
            for (index, face) in a.faces.iter().enumerate() {
                if let Some(other) = b.faces.iter().position(|f| f == face)
                    && a.definition.bands_meet(index, &b.definition, other, dec)?
                {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "two strips consume one band of their shared support and require simultaneous contact-region trimming",
                    ));
                }
            }
            if separated(&a.definition.surface, &b.definition.surface, dec)?
                || a.definition.sectors_apart(&b.definition, dec)?
            {
                continue;
            }
            if let (Some(first), Some(second)) = (&hulls[i], &hulls[j])
                && match touch {
                    Some(touch) => hulls_touch(first, second, touch, dec)?,
                    None => hulls_apart(first, second, dec)?,
                }
            {
                continue;
            }
            // Nothing solves a fitted band against another surface.
            if [a, b].iter().any(|run| {
                matches!(
                    run.definition.surface,
                    Surface::BSpline(
                        SplineSupport::RollingSaddle(_) | SplineSupport::RollingSphere(_)
                    )
                )
            }) {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "two blend bands are not proved apart and need trimmed continuation",
                ));
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
    let mut swallowed = Vec::new();
    let mut rolled = Vec::new();
    let mut bounded_runs = BTreeSet::new();
    let mut continuations: Vec<Continuation> = Vec::new();
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
        // Every end or junction case not built yet says which it is.
        *unfinished = Some(if active.len() == 1 {
            UNSUPPORTED_END
        } else {
            UNSUPPORTED_CORNER
        });
        if let [i] = active[..]
            && let Some(c) = cliffs.get(&(i, old))
        {
            // The rails end at the feet, the cut curves close the strip and
            // fill the cut faces' loops, and the edges leaving the removed
            // region end on the cut.
            let faces = runs[i].faces;
            for (face, v) in faces.into_iter().zip(c.feet) {
                runs[i].ends.insert((old, face), v);
            }
            runs[i].caps.insert(old, c.caps.clone());
            for (face, v, cap) in &c.gaps {
                if edit.gaps.insert((*face, *v), Gap::Edge(*cap)).is_some() {
                    return Err(Refusal::new(
                        UNSUPPORTED_CORNER,
                        "two cliff ends cut one face at one vertex",
                    ));
                }
            }
            edit.removed.extend(c.removed.iter().copied());
            for (key, v) in &c.ends {
                end_changes.insert(*key, *v);
            }
            continue;
        }
        if let [i] = active[..]
            && let Some(r) = rollings.get(&(i, old))
        {
            // The rails end at the band's feet, the band's end section caps
            // the strip, and each band rail replaces the mitre in its round.
            let faces = runs[i].faces;
            let feet = r.feet.clone().map(|p| vertex(model, p));
            runs[i].ends.insert((old, faces[0]), feet[0]);
            runs[i].ends.insert((old, faces[1]), feet[1]);
            let foot_of = |rail: usize| feet[usize::from(r.rails[1] == rail)];
            let curve = model
                .geom_mut()
                .add_curve(Curve::Circle(r.band.section(dec)?));
            let cap = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(foot_of(0)),
                end: Some(foot_of(1)),
            });
            runs[i].caps.insert(old, vec![cap]);
            let mitre_edge = model.edge(r.edge).clone();
            let mut arcs = Vec::with_capacity(2);
            for (rail, face) in r.faces.into_iter().enumerate() {
                let own = foot_of(rail);
                let curve = model
                    .geom_mut()
                    .add_curve(Curve::Ellipse(r.band.rail(rail, dec)?));
                let (start, end, sense) = if mitre_edge.start == Some(old) {
                    (own, r.apex, Sense::Same)
                } else {
                    (r.apex, own, Sense::Reversed)
                };
                let arc = model.add_edge(Edge {
                    curve,
                    sense,
                    start: Some(start),
                    end: Some(end),
                });
                edit.edges.insert((face, r.edge), arc);
                arcs.push((face, arc));
            }
            for (edge, adjacent) in &uses {
                let e = model.edge(*edge);
                if *edge == runs[i].edge
                    || *edge == r.edge
                    || (e.start != Some(old) && e.end != Some(old))
                {
                    continue;
                }
                let v = if adjacent.iter().any(|(f, _)| *f == faces[0]) {
                    feet[0]
                } else if adjacent.iter().any(|(f, _)| *f == faces[1]) {
                    feet[1]
                } else {
                    return Err(bad("an edge at a rolling mitre meets neither strip face"));
                };
                end_changes.insert((*edge, old), v);
            }
            let mut near = incident.clone();
            near.extend(
                uses.iter()
                    .filter_map(|(edge, adjacent)| {
                        let e = model.edge(*edge);
                        (e.start == Some(r.apex) || e.end == Some(r.apex)).then_some(adjacent)
                    })
                    .flat_map(|adjacent| adjacent.iter().map(|(f, _)| *f)),
            );
            let surface = Surface::BSpline(SplineSupport::RollingMitre(r.band.clone()));
            surface_guard(model, body, &surface, &near, dec)?;
            rolled.push((old, i, arcs));
            continue;
        }
        if let [i] = active[..]
            && let Some(m) = mitres.get(&(i, old))
        {
            let (corner, shared) = (&m.corner, m.shared);
            let feet: Vec<_> = (0..3)
                .map(|s| {
                    if s == shared {
                        m.apex
                    } else {
                        vertex(model, corner.sector.feet[s].clone())
                    }
                })
                .collect();
            let faces = runs[i].faces;
            let [a, b] = corner.strips[0];
            runs[i].ends.insert((old, faces[0]), feet[a]);
            runs[i].ends.insert((old, faces[1]), feet[b]);
            let curve = model.geom_mut().add_curve(corner.boundary(0));
            let cap = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(feet[a]),
                end: Some(feet[b]),
            });
            runs[i].caps.insert(old, vec![cap]);
            // Each earlier strip's boundary arc replaces the mitre in its
            // loop, from the mitre's far end to its own foot.
            let mitre_edge = model.edge(m.edge).clone();
            let mut arcs = Vec::with_capacity(2);
            for (k, face) in [(1, m.faces[0]), (2, m.faces[1])] {
                let [p, q] = corner.strips[k];
                let own = feet[if p == shared { q } else { p }];
                let (start, end) = if mitre_edge.start == Some(old) {
                    (own, m.apex)
                } else {
                    (m.apex, own)
                };
                let curve = model.geom_mut().add_curve(corner.boundary(k));
                let arc = model.add_edge(Edge {
                    curve,
                    sense: if start == feet[p] {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    },
                    start: Some(start),
                    end: Some(end),
                });
                edit.edges.insert((face, m.edge), arc);
                arcs.push((face, arc));
            }
            for (edge, adjacent) in &uses {
                let e = model.edge(*edge);
                if *edge == runs[i].edge
                    || *edge == m.edge
                    || (e.start != Some(old) && e.end != Some(old))
                {
                    continue;
                }
                let v = if adjacent.iter().any(|(f, _)| *f == faces[0]) {
                    feet[a]
                } else if adjacent.iter().any(|(f, _)| *f == faces[1]) {
                    feet[b]
                } else {
                    return Err(bad("an edge at a mitre corner meets neither strip face"));
                };
                end_changes.insert((*edge, old), v);
            }
            let mut near = incident.clone();
            near.extend(
                uses.iter()
                    .filter_map(|(edge, adjacent)| {
                        let e = model.edge(*edge);
                        (e.start == Some(m.apex) || e.end == Some(m.apex)).then_some(adjacent)
                    })
                    .flat_map(|adjacent| adjacent.iter().map(|(f, _)| *f)),
            );
            corner_guard(model, body, corner, &near, dec)?;
            swallowed.push((old, i, arcs));
            continue;
        }
        if incident.len() > 3 && active.len() == incident.len() {
            let mut network = None;
            let sphere = fillet_vertex::build(
                model,
                body,
                old,
                &active,
                &mut runs,
                &incident,
                &uses,
                step,
                &mut edit,
                &mut network,
                dec,
            )?;
            bounded_runs.extend(active.iter().copied());
            if let Some((origin, surface, convex)) = sphere {
                corners.push((old, active, origin, surface, convex));
            } else {
                let junction = network
                    .ok_or_else(|| bad("a split corner did not emit its contact network"))?;
                junction.clear_of(&active, &hulls, &continuations, dec)?;
                continuations.push(junction);
            }
            continue;
        }
        if let [i] = active[..]
            && incident.len() == 4
        {
            fillet_vertex_end::build(
                model,
                body,
                old,
                &mut runs[i],
                &incident,
                &uses,
                &mut edit,
                &mut end_changes,
                dec,
            )?;
            bounded_runs.insert(i);
            continue;
        }
        if let [i, j] = active[..]
            && incident.len() == 4
        {
            let (first, second) = runs.split_at_mut(j);
            let continuation = fillet_vertex_end::build_pair(
                model,
                body,
                old,
                [&mut first[i], &mut second[0]],
                &incident,
                &uses,
                &mut edit,
                &mut end_changes,
                step,
                closure,
                dec,
            )?;
            if let Some(continuation) = continuation {
                continuation.clear_of(&[i, j], &hulls, &continuations, dec)?;
                continuations.push(continuation);
            }
            bounded_runs.extend([i, j]);
            continue;
        }
        if incident.len() != 3 {
            return Err(bad("a strip end needs a regular three-face vertex"));
        }
        match active.as_slice() {
            [i] => {
                let faces = runs[*i].faces;
                let third = *incident.iter().find(|f| !faces.contains(f)).unwrap();
                // A cylindrical or spherical end wall cuts the extended strip;
                // its end curve is the blend/wall intersection, not an
                // invented cap. A ring strip ends on a plane through its axis
                // at a meridian.
                let ring = matches!(
                    &runs[*i].definition.contacts,
                    [Curve::Circle(_), Curve::Circle(_)]
                );
                let curved = matches!(
                    model.geom().surface(model.face(third).surface),
                    Surface::Cylinder(_) | Surface::Sphere(_)
                );
                let (a, b, meridian) = if ring {
                    let (a, b, meridian) = ring_end(model, &runs[*i], third, old, dec)?;
                    (a, b, Some(meridian))
                } else if curved {
                    (
                        wall_foot(model, &uses, &runs[*i], faces[0], third, old, dec)?,
                        wall_foot(model, &uses, &runs[*i], faces[1], third, old, dec)?,
                        None,
                    )
                } else {
                    (
                        line_plane(rail(&runs[*i], faces[0])?, plane(model, third)?, dec)?,
                        line_plane(rail(&runs[*i], faces[1])?, plane(model, third)?, dec)?,
                        None,
                    )
                };
                let [a, b] = [vertex(model, a), vertex(model, b)];
                for (face, v) in faces.into_iter().zip([a, b]) {
                    runs[*i].ends.insert((old, face), v);
                }
                let cap = match meridian {
                    Some(meridian) => meridian_edge(model, meridian, a, b, dec)?,
                    None => {
                        let surface = runs[*i].definition.surface.clone();
                        let third_surface = model.geom().surface(model.face(third).surface).clone();
                        connector(
                            model,
                            &surface,
                            &third_surface,
                            a,
                            b,
                            &[&runs[*i].definition],
                            None,
                            None,
                            dec,
                        )?
                    }
                };
                runs[*i].caps.insert(old, vec![cap]);
                edit.gaps.insert((third, old), Gap::Edge(cap));
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
                let (p, qa, unequal) = match &join {
                    Some((_, p, q)) => (p.clone(), q.clone(), None),
                    None => {
                        if closure == catallax_recipe::FilletCorner::Smooth {
                            return Err(Refusal::new(
                                UNSUPPORTED_CORNER,
                                "the smooth junction lacks a certified terminating contact sheet",
                            ));
                        }
                        if runs[*i].definition.convex != runs[*j].definition.convex {
                            return Err(bad(
                                "a convex and a concave strip meeting at a vertex need a general corner patch",
                            ));
                        }
                        let p = line_line(rail(&runs[*i], common)?, rail(&runs[*j], common)?, dec)?;
                        let qa = line_plane(rail(&runs[*i], a)?, plane(model, b)?, dec)?;
                        let qb = line_plane(rail(&runs[*j], b)?, plane(model, a)?, dec)?;
                        if same(&qa, &qb, dec)? {
                            (p, qa, None)
                        } else {
                            let (q, shallow, m) = unequal_mitre(
                                model,
                                &uses,
                                [&runs[*i], &runs[*j]],
                                [a, b],
                                [qa, qb],
                                old,
                                dec,
                            )?;
                            (p, q, Some((shallow, m)))
                        }
                    }
                };
                // A sphere strip's contact on its plane is one point, so its
                // two ends share that vertex there.
                let point = |face: FaceId| {
                    [*i, *j].into_iter().find_map(|index| {
                        let run = &runs[index];
                        let at = run.faces.iter().position(|f| *f == face)?;
                        if run.definition.point_contact() != Some(at) {
                            return None;
                        }
                        run.ends
                            .iter()
                            .find(|((_, f), _)| *f == face)
                            .map(|(_, v)| *v)
                    })
                };
                let p = match point(common) {
                    Some(v) => v,
                    None => vertex(model, p),
                };
                let q = match point(a).or_else(|| point(b)) {
                    Some(v) => v,
                    None => vertex(model, qa),
                };
                let unequal = unequal.map(|(shallow, m)| (shallow, vertex(model, m)));
                for (index, run) in [*i, *j].into_iter().enumerate() {
                    for face in runs[run].faces {
                        let v = match unequal {
                            _ if face == common => p,
                            Some((shallow, m)) if shallow == index => m,
                            _ => q,
                        };
                        runs[run].ends.insert((old, face), v);
                    }
                }
                let (from, to) = match unequal {
                    Some((_, m)) => (p, m),
                    None => (p, q),
                };
                let cap = match join {
                    Some((meridian, ..)) => meridian_edge(model, meridian, from, to, dec)?,
                    None => connector(
                        model,
                        &runs[*i].definition.surface,
                        &runs[*j].definition.surface,
                        from,
                        to,
                        &[&runs[*i].definition, &runs[*j].definition],
                        None,
                        None,
                        dec,
                    )?,
                };
                for index in [*i, *j] {
                    runs[index].caps.insert(old, vec![cap]);
                }
                // The deeper strip runs on past the shallower one's outer rail
                // and is cut by that rail's face down to the unselected edge.
                if let Some((shallow, m)) = unequal {
                    let [deep, shallow] = if shallow == 0 { [*j, *i] } else { [*i, *j] };
                    let face = if shallow == *i { a } else { b };
                    let surface = model.geom().surface(model.face(face).surface).clone();
                    // The cut stays in the deeper sector: it never reaches
                    // the deeper strip's shared rail.
                    let side = line_plane(rail(&runs[deep], common)?, plane(model, face)?, dec)?;
                    let cut = connector(
                        model,
                        &runs[deep].definition.surface,
                        &surface,
                        m,
                        q,
                        &[&runs[deep].definition],
                        Some(&side),
                        None,
                        dec,
                    )?;
                    runs[deep].caps.get_mut(&old).unwrap().push(cut);
                    edit.gaps.insert((face, old), Gap::Edge(cut));
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
                let feet = corner.sector.feet.clone().map(|p| vertex(model, p));
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
                    runs[*index].caps.insert(old, vec![cap]);
                }
                corner_guard(model, body, &corner, &incident, dec)?;
                corners.push((
                    old,
                    active,
                    corner.origin,
                    corner.sector.surface,
                    corner.sector.convex,
                ));
            }
            _ => {
                return Err(bad(
                    "four or more selected strips at a vertex need a general corner patch",
                ));
            }
        }
    }
    if !continuations.is_empty() {
        for (old, _, _, surface, _) in &corners {
            let Surface::Sphere(sphere) = surface else {
                return Err(Refusal::new(
                    UNSUPPORTED_CORNER,
                    "a neighbouring corner needs a complete affected-region bound",
                ));
            };
            let at = model.geom().position(model.vertex(*old).point);
            let bounds = corner_cell_box(&sphere.center, &at);
            for continuation in &continuations {
                if !hulls_apart(continuation.region_hull(), &bounds, dec)? {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "a continued junction meets another corner and needs region trimming",
                    ));
                }
            }
        }
    }
    for i in bounded_runs {
        let run = &runs[i];
        let input = model.edge(run.edge);
        let [Some(start), Some(end)] = [input.start, input.end] else {
            return Err(bad("a bounded fillet run needs finite strip ends"));
        };
        let original = model
            .geom()
            .position(model.vertex(end).point)
            .minus(&model.geom().position(model.vertex(start).point));
        for face in run.faces {
            let a = model
                .geom()
                .position(model.vertex(run.ends[&(start, face)]).point);
            let b = model
                .geom()
                .position(model.vertex(run.ends[&(end, face)]).point);
            if sign(&b.minus(&a).dot(&original), dec)? != Sign::Positive {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "fillet ends consume a strip and need contact-region trimming",
                ));
            }
        }
    }
    *unfinished = None;
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
        let point = run.definition.point_contact();
        let mut rails = Vec::with_capacity(2);
        for (i, face) in run.faces.iter().enumerate() {
            // Every other face is clear of the band, so a closed rail stays in
            // its face if one of its points does.
            if let (
                Surface::BSpline(SplineSupport::RollingPipe(pipe)),
                Curve::Ellipse(_) | Curve::BSpline(_),
            ) = (&run.definition.surface, &run.definition.contacts[i])
            {
                let row = if matches!(run.definition.supports[i], Surface::Plane(_)) {
                    Rational::zero()
                } else {
                    Rational::one()
                };
                let p = pipe.point(&Rational::zero(), &row);
                if !catallax_boolean::face_contains(model, *face, &p, dec)? {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "a contact rail leaves the input face and needs neighbour continuation",
                    ));
                }
            }
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
            // A rail runs up the carrier's chart, which need not be the input
            // edge's way round the ring.
            let mut sense = input.sense;
            if let (
                Surface::BSpline(SplineSupport::RollingPipe(pipe)),
                Curve::BSpline(_),
                Curve::Ellipse(ring),
            ) = (
                &run.definition.surface,
                &run.definition.contacts[i],
                &run.definition.edge,
            ) && sign(&ring.normal.dot(&pipe.chart_normal()), dec)? == Sign::Negative
            {
                sense = reverse(sense);
            }
            if let (Surface::BSpline(band), Curve::Intersection(ring)) =
                (&run.definition.surface, &run.definition.edge)
                && matches!(
                    band,
                    SplineSupport::RollingSaddle(_) | SplineSupport::RollingSphere(_)
                )
            {
                let contact = &run.definition.contacts[i];
                let Curve::BSpline(rail) = contact else {
                    return Err(bad("a saddle contact is not its rail"));
                };
                let p = rail.control(0).clone();
                if !catallax_boolean::face_contains(model, *face, &p, dec)? {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "a contact rail leaves the input face and needs neighbour continuation",
                    ));
                }
                let at = model.geom().position(ring.witness);
                let [_, along] =
                    catallax_intersect::tangent_at(model.geom(), input.curve, ring.witness, dec)?;
                let with_chart = match band {
                    SplineSupport::RollingSaddle(saddle) => {
                        saddle.definition().runs_with_chart(&at, &along, dec)?
                    }
                    SplineSupport::RollingSphere(band) => {
                        band.definition().runs_with_chart(&at, &along, dec)?
                    }
                    _ => unreachable!("matched above"),
                };
                if !with_chart {
                    sense = reverse(sense);
                }
            }
            // A point contact leaves no rail: its face loses the input edge.
            if point == Some(i) {
                edit.removed.insert((*face, run.edge));
                rails.push((None, reverse(run.senses[i]), *face));
                continue;
            }
            let curve = model
                .geom_mut()
                .add_curve(run.definition.contacts[i].clone());
            let edge = model.add_edge(Edge {
                curve,
                sense,
                start: input.start.map(|v| run.ends[&(v, *face)]),
                end: input.end.map(|v| run.ends[&(v, *face)]),
            });
            edit.edges.insert((*face, run.edge), edge);
            rails.push((Some(edge), reverse(run.senses[i]), *face));
        }
        let paths = if input.start.is_none() && input.end.is_none() {
            rails
                .iter()
                .map(|(edge, sense, _)| {
                    edge.map(|edge| vec![(edge, *sense)])
                        .ok_or_else(|| bad("a ring strip keeps both rails"))
                })
                .collect::<Result<_, _>>()?
        } else {
            let mut path = Vec::with_capacity(4);
            for (edge, sense, face) in rails {
                if let Some(edge) = edge {
                    path.push((edge, sense));
                }
                let old = if sense == Sense::Same {
                    input.end.unwrap()
                } else {
                    input.start.unwrap()
                };
                for (cap, sense) in chain(model, &run.caps[&old], run.ends[&(old, face)])? {
                    cap_uses.insert(cap, sense);
                    path.push((cap, sense));
                }
            }
            vec![path]
        };
        // A fitted ring band's normal points to the ball, where the other
        // bands' points away.
        let toward_ball = matches!(
            run.definition.surface,
            Surface::BSpline(SplineSupport::RollingSaddle(_) | SplineSupport::RollingSphere(_))
        );
        let face = crate::extrude::face(
            model,
            run.definition.surface.clone(),
            if run.definition.convex != toward_ball {
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
    for (old, indices, origin, surface, convex) in &corners {
        // Each great-circle arc is shared by its strip and the sphere polygon.
        let caps: Vec<_> = indices
            .iter()
            .map(|index| {
                let [cap] = runs[*index].caps[old][..] else {
                    unreachable!("a corner strip end has one cap");
                };
                (cap, reverse(cap_uses[&cap]))
            })
            .collect();
        let path = close(model, caps)?;
        // A concave corner's material lies outside its ball.
        let face = crate::extrude::face(
            model,
            surface.clone(),
            if *convex {
                Sense::Same
            } else {
                Sense::Reversed
            },
            vec![path],
            None,
            vec![origin.clone()],
        );
        edit.faces.push((runs[indices[0]].faces[0], face));
    }
    for (old, i, arcs) in &swallowed {
        let corner = &mitres[&(*i, *old)].corner;
        let [cap] = runs[*i].caps[old][..] else {
            unreachable!("a mitre corner's strip end has one cap");
        };
        let mut sides = vec![(cap, reverse(cap_uses[&cap]))];
        // Each arc is used once by its earlier strip, as the mitre was, and
        // once, reversed, by the sphere sector.
        let mitre = mitres[&(*i, *old)].edge;
        for (face, arc) in arcs {
            let sense = model
                .face(*face)
                .loops
                .iter()
                .find_map(|lp| match model.loop_(*lp) {
                    Loop::Coedges(cs) => cs
                        .iter()
                        .map(|c| model.coedge(*c))
                        .find(|c| c.edge == mitre)
                        .map(|c| c.sense),
                    Loop::Vertex(_) => None,
                })
                .ok_or_else(|| bad("an earlier strip does not use its mitre"))?;
            sides.push((*arc, reverse(sense)));
        }
        let path = close(model, sides)?;
        let face = crate::extrude::face(
            model,
            corner.sector.surface.clone(),
            if corner.sector.convex {
                Sense::Same
            } else {
                Sense::Reversed
            },
            vec![path],
            None,
            vec![corner.origin.clone()],
        );
        edit.faces.push((runs[*i].faces[0], face));
    }
    for (old, i, arcs) in &rolled {
        let r = &rollings[&(*i, *old)];
        let [cap] = runs[*i].caps[old][..] else {
            unreachable!("a rolling mitre's strip end has one cap");
        };
        let mut sides = vec![(cap, reverse(cap_uses[&cap]))];
        // Each rail is used once by its round, as the mitre was, and once,
        // reversed, by the band.
        for (face, arc) in arcs {
            sides.push((*arc, reverse(edge_sense(model, *face, r.edge)?)));
        }
        let path = close(model, sides)?;
        let face = crate::extrude::face(
            model,
            Surface::BSpline(SplineSupport::RollingMitre(r.band.clone())),
            if r.band.convex() {
                Sense::Same
            } else {
                Sense::Reversed
            },
            vec![path],
            None,
            vec![r.origin.clone()],
        );
        edit.faces.push((runs[*i].faces[0], face));
    }
    let result = edge_blend::apply(model, body, &edit)?;
    if let Some(trimming) = trimming {
        trimming.finish(dec, None, None);
    }
    // Measurement reads rings as the body they rounded plus their exact change,
    // which holds only when nothing else changed it in this step.
    let result = if runs.iter().all(|run| {
        matches!(
            run.definition.surface,
            Surface::BSpline(SplineSupport::RollingSaddle(_) | SplineSupport::RollingSphere(_))
        )
    }) {
        let fit = runs.iter().try_fold(Rational::zero(), |worst, run| {
            let bound = match &run.definition.surface {
                Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => saddle.bound(),
                Surface::BSpline(SplineSupport::RollingSphere(band)) => band.bound(),
                _ => return Err(bad("a ring run has no fitted carrier")),
            };
            Ok(worst.max(bound.upper.clone()))
        })?;
        let mut marked = model.body(result).clone();
        marked.marks.push(Mark::new(
            crate::certified::FIT_MARK,
            BoundKind::Distance,
            fit,
            "fitted rolling-ball ring faces; the mark alone is not authority",
        ));
        let marked = model.add_body(marked);
        model.record_blend_base(marked, body);
        marked
    } else {
        result
    };
    let mut definitions = Vec::with_capacity(runs.len() + continuations.len());
    definitions.extend(runs.into_iter().map(|run| run.definition));
    definitions.extend(
        continuations
            .into_iter()
            .map(|junction| junction.definition),
    );
    Ok((result, definitions))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::BudgetPolicy;
    use catallax_recipe::Recipe;

    fn q(text: &str) -> Real {
        Real::from(Rational::parse(text).unwrap())
    }

    /// A straight strip hull from x = -1 to x = 1 whose cross-section has its
    /// edge point and two contacts at the given (y, z) points.
    fn strip(section: [[&str; 2]; 3]) -> Vec<Point3> {
        ["-1", "1"]
            .into_iter()
            .flat_map(|x| section.map(|[y, z]| Point3::new(q(x), q(y), q(z))))
            .collect()
    }

    #[test]
    fn a_block_face_polygon_keeps_every_corner_whatever_the_coedge_senses() {
        let mut model = Model::default();
        let body = catallax_topology::block(
            &mut model,
            [1, 1, 1].map(Rational::from_i64),
            &catallax_shapes::Pose::identity(),
        )
        .unwrap();
        let mut dec = Decider::new(BudgetPolicy::default());
        let face = model
            .faces_of_body(body)
            .into_iter()
            .find(|face| {
                matches!(model.geom().surface(model.face(*face).surface),
                    Surface::Plane(p) if p.normal.x.as_rational() == Some(Rational::from_i64(-1)))
            })
            .unwrap();
        let corners = polygon(&model, face).unwrap();
        for [y, z] in [["0", "0"], ["1", "0"], ["0", "1"], ["1", "1"]] {
            let corner = Point3::new(q("0"), q(y), q(z));
            assert!(
                corners.iter().any(|c| same(c, &corner, &mut dec).unwrap()),
                "missing ({y}, {z})"
            );
        }
        // This strip lies on the far side of its chord plane y + z = 17/10
        // from three corners, but crosses the square near (1, 1).
        let crossing = strip([["11/10", "11/10"], ["3/5", "11/10"], ["11/10", "3/5"]]);
        assert!(!hulls_apart(&crossing, &corners, &mut dec).unwrap());
        let apart = strip([["31/10", "31/10"], ["13/5", "31/10"], ["31/10", "13/5"]]);
        assert!(hulls_apart(&apart, &corners, &mut dec).unwrap());
    }

    // The ring at height `z` between a plane and a perpendicular cylinder, its
    // strip built with the given material side, and the rail check's verdict.
    fn rails(steps: &str, body: &str, z: i64, convex: bool) -> Result<(), Refusal> {
        let recipe = Recipe::parse(&format!(
            r#"{{"format":"catallax-recipe/3","name":"ring","steps":[{steps}],"result":"{body}"}}"#
        ))
        .unwrap();
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 20_000_000,
            ..Default::default()
        });
        let mut built = crate::evaluate(&recipe, &mut dec).unwrap();
        let model = &mut built.model;
        let uses = uses(model, built.result);
        let (edge, incident) = uses
            .iter()
            .find(|(edge, incident)| {
                matches!(model.geom().curve(model.edge(**edge).curve),
                    Curve::Circle(c) if c.center.z.as_rational() == Some(Rational::from_i64(z)))
                    && incident.iter().any(|(face, _)| {
                        matches!(
                            model.geom().surface(model.face(*face).surface),
                            Surface::Cylinder(_)
                        )
                    })
            })
            .unwrap();
        let [(a, sa), (b, sb)] = incident.as_slice() else {
            panic!("two faces");
        };
        let faces = [*a, *b];
        let definition = FilletDefinition::new(
            Origin {
                step: StepId::new("round").unwrap(),
                role: Role::Blend {
                    edge: model.edge_name(built.result, *edge).unwrap(),
                },
            },
            faces.map(|f| model.geom().surface(model.face(f).surface).clone()),
            faces.map(|f| model.face(f).sense == Sense::Same),
            model.geom().curve(model.edge(*edge).curve).clone(),
            Real::one(),
            convex,
            &mut dec,
        )?;
        let run = Run {
            edge: *edge,
            faces,
            senses: [*sa, *sb],
            definition,
            ends: BTreeMap::new(),
            caps: BTreeMap::new(),
        };
        rails_on_faces(model, &run, &mut dec)
    }

    const BLOCK: &str = r#"{"id":"block-section","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"line","id":"s","to":["20","-20"]},{"kind":"line","id":"e","to":["20","20"]},{"kind":"line","id":"n","to":["-20","20"]},{"kind":"line","id":"w","to":["-20","-20"]}]]},{"id":"block","op":"extrude","profile":"block-section","distance":"20"}"#;

    // A strip on the wrong material side puts a rail beyond its face. The
    // check refuses it whatever the convexity test said, on every kind of
    // ring, and passes the strip on the right side.
    #[test]
    fn a_rail_off_its_face_is_refused() {
        let boss = r#"{"id":"boss","op":"cylinder","radius":"3","height":"10","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}"#;
        let mouth = format!(
            r#"{BLOCK},{{"id":"hole","op":"cylinder","radius":"3","height":"30","pose":{{"rotation":["1","0","0","0"],"translation":["0","0","-5"]}}}},{{"id":"cut","op":"boolean","kind":"subtract","target":"block","tool":"hole"}}"#
        );
        let blind = format!(
            r#"{BLOCK},{{"id":"hole","op":"cylinder","radius":"3","height":"20","pose":{{"rotation":["1","0","0","0"],"translation":["0","0","10"]}}}},{{"id":"cut","op":"boolean","kind":"subtract","target":"block","tool":"hole"}}"#
        );
        let base = format!(
            r#"{BLOCK},{{"id":"boss","op":"cylinder","radius":"3","height":"20","pose":{{"rotation":["1","0","0","0"],"translation":["0","0","10"]}}}},{{"id":"joined","op":"boolean","kind":"union","target":"block","tool":"boss"}}"#
        );
        for (steps, body, z, convex) in [
            (boss.to_string(), "boss", 10, true),
            (mouth, "cut", 20, true),
            (blind, "cut", 10, false),
            (base, "joined", 20, false),
        ] {
            assert_eq!(rails(&steps, body, z, convex), Ok(()), "{body} right side");
            let refusal = rails(&steps, body, z, !convex).unwrap_err();
            assert_eq!(refusal.code, "fillet.contact-off-face", "{body} wrong side");
        }
    }
}
