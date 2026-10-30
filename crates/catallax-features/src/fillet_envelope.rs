//! Rim-contact continuation of a closed plane/cylinder fillet ring.
//!
//! When countersink mouths cut the plane close to the junction, the rolling
//! ball rests on each mouth edge between two transitions instead of on the
//! plane. The ring then splits into elementary torus pieces between the holes
//! and one checked rolling-envelope face per hole. The plane's junction loop
//! and the touched mouth loops merge into one loop, and each mouth edge splits
//! into the arc under the envelope and the arc the plane keeps.
//!
//! Only this stratum is covered: one closed concave circular ring between a
//! plane and a perpendicular cylinder, with every intruding face either a
//! touched coaxial-parallel countersink cone below the plane or a face whose
//! finite bounds lie strictly below it. Anything else keeps refusing.

use std::collections::BTreeMap;

use catallax_numbers::{BoundKind, Decider, Mark, Rational, Real, Refusal, Sign};
use catallax_shapes::{
    BSplineCurve, Circle, Cone, Curve, Point3, PointDef, RollingEnvelope,
    RollingEnvelopeDefinition, SplineBoundary, Surface, Vector3,
};
use catallax_topology::{
    BodyId, Edge, EdgeId, FaceId, Loop, Lump, Model, Sense, Shell, Vertex, VertexId,
};

use crate::{
    extrude::reverse,
    fillet::{FilletDefinition, separated},
    profile::sign,
};

type Uses = BTreeMap<EdgeId, Vec<(FaceId, Sense)>>;

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("fillet.unsupported-overflow", detail)
}
fn real_sqrt(value: &Real) -> Real {
    value
        .as_rational()
        .and_then(|r| crate::extrude::sqrt(&r))
        .map(Real::from)
        .unwrap_or_else(|| value.sqrt())
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &real_sqrt(&v.norm_squared())))
}
fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Zero)
}
fn positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Positive)
}
fn sense_sign(sense: Sense) -> i8 {
    if sense == Sense::Same { 1 } else { -1 }
}
fn vertex(model: &mut Model, p: Point3) -> VertexId {
    let point = model.geom_mut().add_point(PointDef::Explicit(p));
    model.add_vertex(Vertex { point })
}

// A mouth that the ball reaches, with its exact envelope definition.
struct Hole {
    face: FaceId,
    edge: EdgeId,
    circle: Circle,
    definition: RollingEnvelopeDefinition,
}

// Per-hole vertices: [cylinder foot, plane foot] at the early and late
// transition in counter-clockwise order about the axis.
struct Feet {
    early: [VertexId; 2],
    late: [VertexId; 2],
    // Whether the envelope's u direction runs from early to late.
    forward: bool,
}

/// Builds the continued ring, or returns `None` when the overflow is not this
/// stratum so the caller keeps its original refusal.
pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    edge: EdgeId,
    faces: [FaceId; 2],
    definition: &FilletDefinition,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<Option<BodyId>, Refusal> {
    let input = model.edge(edge).clone();
    let (Curve::Circle(_), None, None, false) = (
        model.geom().curve(input.curve),
        input.start,
        input.end,
        definition.convex,
    ) else {
        return Ok(None);
    };
    let Surface::Torus(torus) = &definition.surface else {
        return Ok(None);
    };
    let surfaces = faces.map(|f| model.geom().surface(model.face(f).surface).clone());
    let (plane_index, plane, cylinder) = match (&surfaces[0], &surfaces[1]) {
        (Surface::Plane(p), Surface::Cylinder(c)) => (0, p.clone(), c.clone()),
        (Surface::Cylinder(c), Surface::Plane(p)) => (1, p.clone(), c.clone()),
        _ => return Ok(None),
    };
    let plane_face = faces[plane_index];
    let cylinder_face = faces[1 - plane_index];
    // The ball sits on the plane face's outward side.
    let outward = if model.face(plane_face).sense == Sense::Same {
        1
    } else {
        -1
    };
    let axis = unit(&plane.normal).scaled(&Real::from_i64(outward));
    if !zero(&axis.cross(&cylinder.axis).norm_squared(), dec)? {
        return Ok(None);
    }
    let lift = plane.origin.minus(&cylinder.origin).dot(&axis);
    let origin = cylinder.origin.plus(&axis.scaled(&lift));
    let body_radius = real_sqrt(&cylinder.radius_sq);
    let radius = torus.minor.clone();
    let centre = &body_radius + &radius;

    // Touched mouths: single-edge hole loops of the plane face whose other face
    // is a coaxial-parallel cone that narrows below the plane.
    let mut holes = Vec::new();
    let mut touched_faces = Vec::new();
    for lp in model.face(plane_face).loops.clone() {
        let Loop::Coedges(cs) = model.loop_(lp) else {
            continue;
        };
        let [c] = cs.as_slice() else { continue };
        let hole_edge = model.coedge(*c).edge;
        if hole_edge == edge {
            continue;
        }
        let Curve::Circle(circle) = model.geom().curve(model.edge(hole_edge).curve).clone() else {
            continue;
        };
        let Some(cone_face) = uses
            .get(&hole_edge)
            .and_then(|incident| incident.iter().find(|(f, _)| *f != plane_face))
            .map(|(f, _)| *f)
        else {
            continue;
        };
        let Surface::Cone(cone) = model.geom().surface(model.face(cone_face).surface).clone()
        else {
            continue;
        };
        let offset_vector = circle.center.minus(&origin);
        let offset = real_sqrt(&offset_vector.norm_squared());
        let mouth = real_sqrt(&circle.radius_sq);
        let reach = &(&offset - &mouth) - &centre;
        if positive(&reach, dec)? {
            // Beyond the plane rail: the ball never reaches this mouth.
            continue;
        }
        if zero(&reach, dec)? {
            return Err(bad(
                "a mouth tangent to the plane rail needs a tangent continuation",
            ));
        }
        let Some(slope) = narrowing(&cone, &circle, &axis, dec)? else {
            return Ok(None);
        };
        let toward = offset_vector.scaled(&(&Real::one() / &offset));
        let hole_definition = RollingEnvelopeDefinition::new(
            origin.clone(),
            axis.clone(),
            toward,
            body_radius.clone(),
            mouth,
            offset,
            radius.clone(),
            [Real::one(), slope],
            dec,
        )?;
        touched_faces.push(cone_face);
        holes.push(Hole {
            face: cone_face,
            edge: hole_edge,
            circle,
            definition: hole_definition,
        });
    }
    if holes.is_empty() {
        return Ok(None);
    }
    // Disjoint mouths keep each ring arc between holes outside every mouth.
    for (i, a) in holes.iter().enumerate() {
        for c in &holes[i + 1..] {
            let gap = &(a.definition.mouth() + c.definition.mouth()).square()
                - &a.circle.center.minus(&c.circle.center).norm_squared();
            if !positive(&-&gap, dec)? {
                return Ok(None);
            }
        }
    }

    // The strip lies on the ball side of the plane, so a touched cone face
    // must lie on the material side, meeting the plane only at its mouth.
    for hole in &holes {
        if !below_plane(model, hole.face, &origin, &axis, false, dec)? {
            return Ok(None);
        }
    }
    // Every other face must be clear of the strip: separated whole supports,
    // or finite bounds strictly on the material side of the plane.
    for face in model.faces_of_body(body) {
        if face == plane_face || face == cylinder_face || touched_faces.contains(&face) {
            continue;
        }
        let other = model.geom().surface(model.face(face).surface).clone();
        if separated(&definition.surface, &other, dec)? {
            continue;
        }
        // The strip stays within radial distance body + radius of the axis
        // and outside the body cylinder, so a coaxial cylinder wholly beyond
        // either is clear even when its finite face reaches the plane.
        if let Surface::Cylinder(c) = &other
            && zero(&c.axis.cross(&axis).norm_squared(), dec)?
            && zero(&c.origin.minus(&origin).cross(&axis).norm_squared(), dec)?
            && (positive(&(&c.radius_sq - &(&centre + &radius).square()), dec)?
                || positive(&(&cylinder.radius_sq - &c.radius_sq), dec)?)
        {
            continue;
        }
        if !below_plane(model, face, &origin, &axis, true, dec)? {
            return Ok(None);
        }
    }

    // Order the holes counter-clockwise about the axis and require their
    // transition spans to alternate without overlap.
    let reference = holes[0].definition.toward().clone();
    let second = axis.cross(&reference);
    let angle_key = |p: &Point3, dec: &mut Decider| -> Result<(u8, Real, Real), Refusal> {
        let d = p.minus(&origin);
        let (x, y) = (d.dot(&reference), d.dot(&second));
        let upper = positive(&y, dec)? || (zero(&y, dec)? && positive(&x, dec)?);
        Ok((u8::from(!upper), x, y))
    };
    let before = |a: &Point3, b: &Point3, dec: &mut Decider| -> Result<bool, Refusal> {
        let (ha, xa, ya) = angle_key(a, dec)?;
        let (hb, xb, yb) = angle_key(b, dec)?;
        if ha != hb {
            return Ok(ha < hb);
        }
        positive(&(&(&xa * &yb) - &(&ya * &xb)), dec)
    };
    let mut feet_points = Vec::with_capacity(holes.len());
    for hole in &holes {
        let negative = hole.definition.corners(Sign::Negative);
        let positive_corner = hole.definition.corners(Sign::Positive);
        // Within one hole the transitions are less than a half turn apart.
        let d0 = negative[1].minus(&origin);
        let d1 = positive_corner[1].minus(&origin);
        let forward = positive(&d0.cross(&d1).dot(&axis), dec)?;
        let (early, late) = if forward {
            (negative, positive_corner)
        } else {
            (positive_corner, negative)
        };
        feet_points.push((early, late, forward));
    }
    let mut order: Vec<usize> = (0..holes.len()).collect();
    // Insertion sort with exact comparisons on the early plane feet.
    for i in 1..order.len() {
        let mut j = i;
        while j > 0
            && before(
                &feet_points[order[j]].0[1],
                &feet_points[order[j - 1]].0[1],
                dec,
            )?
        {
            order.swap(j, j - 1);
            j -= 1;
        }
    }
    // Walking early, late, early, ... in that order must turn once: exactly
    // one step goes back in the angular order, at the wrap.
    let sequence: Vec<&Point3> = order
        .iter()
        .flat_map(|&i| [&feet_points[i].0[1], &feet_points[i].1[1]])
        .collect();
    let mut wraps = 0;
    for k in 0..sequence.len() {
        if !before(sequence[k], sequence[(k + 1) % sequence.len()], dec)? {
            wraps += 1;
        }
    }
    if wraps != 1 {
        return Err(bad("neighbouring mouth contacts overlap along the ring"));
    }

    // Fit one envelope and reuse it for patterned holes with the same law.
    let mut envelopes: Vec<Option<RollingEnvelope>> = (0..holes.len()).map(|_| None).collect();
    let mut fitted: Vec<RollingEnvelope> = Vec::new();
    for (i, hole) in holes.iter().enumerate() {
        let mut found = None;
        for candidate in &fitted {
            if let Ok(copy) = candidate.reframed(hole.definition.clone(), dec) {
                found = Some(copy);
                break;
            }
        }
        let envelope = match found {
            Some(copy) => copy,
            None => {
                let envelope = RollingEnvelope::fit(hole.definition.clone(), dec)?;
                fitted.push(envelope.clone());
                envelope
            }
        };
        envelopes[i] = Some(envelope);
    }
    let envelopes: Vec<RollingEnvelope> = envelopes.into_iter().map(Option::unwrap).collect();

    // Vertices and the shared transition meridians.
    let mut feet = Vec::with_capacity(holes.len());
    for (early, late, forward) in &feet_points {
        feet.push(Feet {
            early: early.clone().map(|p| vertex(model, p)),
            late: late.clone().map(|p| vertex(model, p)),
            forward: *forward,
        });
    }
    let [cylinder_contact, plane_contact] = if plane_index == 0 {
        [
            definition.contacts[1].clone(),
            definition.contacts[0].clone(),
        ]
    } else {
        [
            definition.contacts[0].clone(),
            definition.contacts[1].clone(),
        ]
    };
    let ccw_sense = |curve: &Curve, dec: &mut Decider| -> Result<Sense, Refusal> {
        let Curve::Circle(c) = curve else {
            return Err(bad("rail contact is not a circle"));
        };
        Ok(if positive(&c.normal.dot(&axis), dec)? {
            Sense::Same
        } else {
            Sense::Reversed
        })
    };
    let cylinder_sense = ccw_sense(&cylinder_contact, dec)?;
    let plane_sense = ccw_sense(&plane_contact, dec)?;
    let cylinder_rail = model.geom_mut().add_curve(cylinder_contact);
    let plane_rail = model.geom_mut().add_curve(plane_contact);
    let meridian = |model: &mut Model, foot: [VertexId; 2]| -> EdgeId {
        let [cylinder_foot, plane_foot] =
            foot.map(|v| model.geom().position(model.vertex(v).point));
        // The ball centre sits radially outward of the cylinder foot.
        let radial = plane_foot.minus(&origin);
        let outward = radial
            .minus(&axis.scaled(&radial.dot(&axis)))
            .scaled(&(&Real::one() / &centre));
        let centre = cylinder_foot.plus(&outward.scaled(&radius));
        // With normal ref x axis, counter-clockwise runs from the cylinder foot
        // (-ref) to the plane foot (-axis) through the corner-facing quarter.
        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: centre,
            normal: outward.cross(&axis),
            radius_sq: radius.square(),
            ref_dir: outward,
        }));
        model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(foot[0]),
            end: Some(foot[1]),
        })
    };
    let mut meridians = Vec::with_capacity(holes.len());
    for f in &feet {
        meridians.push([meridian(model, f.early), meridian(model, f.late)]);
    }

    // Ring arcs between consecutive holes, counter-clockwise about the axis.
    let mut cylinder_arcs = Vec::with_capacity(holes.len());
    let mut plane_arcs = Vec::with_capacity(holes.len());
    for (k, &i) in order.iter().enumerate() {
        let next = order[(k + 1) % order.len()];
        // The sense makes every arc run counter-clockwise from `from` to `to`.
        let arc = |model: &mut Model, curve, sense, from: VertexId, to: VertexId| {
            model.add_edge(Edge {
                curve,
                sense,
                start: Some(from),
                end: Some(to),
            })
        };
        cylinder_arcs.push(arc(
            model,
            cylinder_rail,
            cylinder_sense,
            feet[i].late[0],
            feet[next].early[0],
        ));
        plane_arcs.push(arc(
            model,
            plane_rail,
            plane_sense,
            feet[i].late[1],
            feet[next].early[1],
        ));
    }

    // Each envelope's cylinder rail runs in its chart's u direction, and each
    // mouth circle splits at the two plane feet into near and far arcs.
    let mut rails = Vec::with_capacity(holes.len());
    let mut mouths = Vec::with_capacity(holes.len());
    for (i, hole) in holes.iter().enumerate() {
        let curve = model
            .geom_mut()
            .add_curve(Curve::BSpline(BSplineCurve::boundary(
                envelopes[i].carrier(),
                SplineBoundary::VMin,
            )));
        let (start, end) = if feet[i].forward {
            (feet[i].early[0], feet[i].late[0])
        } else {
            (feet[i].late[0], feet[i].early[0])
        };
        rails.push(model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(start),
            end: Some(end),
        }));
        let mouth_curve = model.edge(hole.edge).curve;
        let sense = if positive(&hole.circle.normal.dot(&axis), dec)? {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let near_point = hole.circle.center.plus(
            &hole
                .definition
                .toward()
                .scaled(&(-&real_sqrt(&hole.circle.radius_sq))),
        );
        // Counter-clockwise about the axis around the mouth centre.
        let orient = |p: &Point3, q: &Point3, dec: &mut Decider| -> Result<bool, Refusal> {
            let c = &hole.circle.center;
            positive(&p.minus(c).cross(&q.minus(c)).dot(&axis), dec)
        };
        let position = |v: VertexId| model.geom().position(model.vertex(v).point);
        let (e, l) = (position(feet[i].early[1]), position(feet[i].late[1]));
        let minor = orient(&e, &l, dec)?;
        let contains = if minor {
            orient(&e, &near_point, dec)? && orient(&near_point, &l, dec)?
        } else {
            orient(&e, &near_point, dec)? || orient(&near_point, &l, dec)?
        };
        // Near arc counter-clockwise from early to late, or from late to early.
        let (near_from, near_to) = if contains {
            (feet[i].early[1], feet[i].late[1])
        } else {
            (feet[i].late[1], feet[i].early[1])
        };
        let ccw = |model: &mut Model, from: VertexId, to: VertexId| {
            model.add_edge(Edge {
                curve: mouth_curve,
                sense,
                start: Some(from),
                end: Some(to),
            })
        };
        let near = ccw(model, near_from, near_to);
        let far = ccw(model, near_to, near_from);
        mouths.push([near, far]);
    }

    // Traversal direction of each replaced closed edge in its old faces.
    let turn =
        |model: &Model, face: FaceId, closed: EdgeId, dec: &mut Decider| -> Result<i8, Refusal> {
            let used = uses[&closed]
                .iter()
                .find(|(f, _)| *f == face)
                .map(|(_, s)| *s)
                .ok_or_else(|| bad("replaced boundary is not on its face"))?;
            let e = model.edge(closed);
            let Curve::Circle(c) = model.geom().curve(e.curve) else {
                return Err(bad("replaced boundary is not a circle"));
            };
            let normal = if positive(&c.normal.dot(&axis), dec)? {
                1
            } else {
                -1
            };
            Ok(normal * sense_sign(e.sense) * sense_sign(used))
        };
    let with = |direction: i8| {
        if direction > 0 {
            Sense::Same
        } else {
            Sense::Reversed
        }
    };
    let plane_turn = turn(model, plane_face, edge, dec)?;
    let cylinder_turn = turn(model, cylinder_face, edge, dec)?;

    // Rebuilt plane loop: plane arcs plus far mouth arcs.
    let mut plane_loop = Vec::new();
    for (k, _) in order.iter().enumerate() {
        plane_loop.push((plane_arcs[k], Some(with(plane_turn))));
    }
    let mut cone_loops = Vec::new();
    for (i, hole) in holes.iter().enumerate() {
        let hole_turn = turn(model, plane_face, hole.edge, dec)?;
        plane_loop.push((mouths[i][1], Some(with(hole_turn))));
        let cone_turn = turn(model, hole.face, hole.edge, dec)?;
        cone_loops.push(chain(
            model,
            vec![
                (mouths[i][0], Some(with(cone_turn))),
                (mouths[i][1], Some(with(cone_turn))),
            ],
        )?);
    }
    let plane_loop = chain(model, plane_loop)?;

    // Rebuilt cylinder loop: cylinder arcs plus envelope rails.
    let mut cylinder_loop = Vec::new();
    for (k, _) in order.iter().enumerate() {
        cylinder_loop.push((cylinder_arcs[k], Some(with(cylinder_turn))));
    }
    for (i, f) in feet.iter().enumerate() {
        let ccw = f.forward;
        cylinder_loop.push((
            rails[i],
            Some(if ccw == (cylinder_turn > 0) {
                Sense::Same
            } else {
                Sense::Reversed
            }),
        ));
    }
    let cylinder_loop = chain(model, cylinder_loop)?;
    let used = |path: &[(EdgeId, Sense)], edge: EdgeId| {
        path.iter().find(|(e, _)| *e == edge).map(|(_, s)| *s)
    };

    // Torus pieces between consecutive holes.
    let mut new_faces = Vec::new();
    let mut meridian_use = BTreeMap::new();
    for (k, &i) in order.iter().enumerate() {
        let next = order[(k + 1) % order.len()];
        let path = chain(
            model,
            vec![
                (
                    cylinder_arcs[k],
                    Some(reverse(used(&cylinder_loop, cylinder_arcs[k]).unwrap())),
                ),
                (meridians[next][0], None),
                (
                    plane_arcs[k],
                    Some(reverse(used(&plane_loop, plane_arcs[k]).unwrap())),
                ),
                (meridians[i][1], None),
            ],
        )?;
        for (e, s) in &path {
            meridian_use.insert(*e, *s);
        }
        new_faces.push(crate::extrude::face(
            model,
            definition.surface.clone(),
            Sense::Reversed,
            vec![path],
            None,
            vec![definition.origin.clone()],
        ));
    }
    for (i, _) in holes.iter().enumerate() {
        let path = chain(
            model,
            vec![
                (
                    rails[i],
                    Some(reverse(used(&cylinder_loop, rails[i]).unwrap())),
                ),
                (
                    meridians[i][1],
                    Some(reverse(meridian_use[&meridians[i][1]])),
                ),
                (
                    mouths[i][0],
                    Some(reverse(used(&cone_loops[i], mouths[i][0]).unwrap())),
                ),
                (
                    meridians[i][0],
                    Some(reverse(meridian_use[&meridians[i][0]])),
                ),
            ],
        )?;
        new_faces.push(crate::extrude::face(
            model,
            Surface::BSpline(envelopes[i].clone().into()),
            Sense::Same,
            vec![path],
            None,
            vec![definition.origin.clone()],
        ));
    }

    // Replace the plane, cylinder and touched cone faces.
    let mut replaced = BTreeMap::new();
    let rebuild =
        |model: &mut Model, face: FaceId, drop: &[EdgeId], add: Vec<(EdgeId, Sense)>| -> FaceId {
            let old = model.face(face).clone();
            let mut paths = Vec::with_capacity(old.loops.len());
            for lp in &old.loops {
                let Loop::Coedges(cs) = model.loop_(*lp) else {
                    continue;
                };
                let path: Vec<_> = cs
                    .iter()
                    .map(|c| {
                        let c = model.coedge(*c);
                        (c.edge, c.sense)
                    })
                    .collect();
                if path.iter().any(|(e, _)| drop.contains(e)) {
                    continue;
                }
                paths.push(path);
            }
            paths.push(add);
            let surface = model.geom().surface(old.surface).clone();
            let origins = model.face_origins(face).to_vec();
            crate::extrude::face(model, surface, old.sense, paths, None, origins)
        };
    let mut plane_drop = vec![edge];
    plane_drop.extend(holes.iter().map(|h| h.edge));
    replaced.insert(
        plane_face,
        rebuild(model, plane_face, &plane_drop, plane_loop),
    );
    replaced.insert(
        cylinder_face,
        rebuild(model, cylinder_face, &[edge], cylinder_loop),
    );
    for (i, hole) in holes.iter().enumerate() {
        replaced.insert(
            hole.face,
            rebuild(model, hole.face, &[hole.edge], cone_loops[i].clone()),
        );
    }

    let input = model.body(body).clone();
    let mut lumps = Vec::with_capacity(input.lumps.len());
    for lump in &input.lumps {
        let mut shells = Vec::new();
        for shell in model.lump(*lump).shells.clone() {
            let old = model.shell(shell).faces.clone();
            let owner = old.contains(&plane_face);
            let mut out: Vec<FaceId> = old
                .iter()
                .map(|f| replaced.get(f).copied().unwrap_or(*f))
                .collect();
            if owner {
                out.extend(new_faces.iter().copied());
            }
            shells.push(model.add_shell(Shell { faces: out }));
        }
        lumps.push(model.add_lump(Lump { shells }));
    }
    let mut output = input;
    output.lumps = lumps;
    let fit = envelopes.iter().fold(Rational::zero(), |worst, envelope| {
        worst.max(envelope.bound().upper.clone())
    });
    output.marks.push(Mark::new(
        crate::certified::FIT_MARK,
        BoundKind::Distance,
        fit,
        "fitted rolling-envelope faces; the mark alone is not authority",
    ));
    let output = model.add_body(output);
    model.record_blend_base(output, body);
    let faults = model.check(output);
    if !faults.is_empty() {
        return Err(Refusal::new(
            "edge-blend.invalid-topology",
            format!("{faults:?}"),
        ));
    }
    Ok(Some(output))
}

// Whether a face's finite bounds lie on the material side of the plane,
// strictly or touching it.
fn below_plane(
    model: &Model,
    face: FaceId,
    origin: &Point3,
    axis: &Vector3,
    strict: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let bounds = catallax_measure::face_bounding_box(
        model,
        face,
        &catallax_measure::Accuracy::default(),
        dec,
    )?
    .value;
    for corner in 0..8 {
        let pick = |axis: usize| {
            Real::from_rational(if corner & (1 << axis) == 0 {
                bounds.min[axis].clone()
            } else {
                bounds.max[axis].clone()
            })
        };
        let p = Point3::new(pick(0), pick(1), pick(2));
        let side = sign(&p.minus(origin).dot(axis), dec)?;
        if side == Sign::Positive || (strict && side == Sign::Zero) {
            return Ok(false);
        }
    }
    Ok(true)
}

// The cone's narrowing slope below the plane when its axis is the mouth axis
// and its section in the plane is the mouth circle; otherwise None.
fn narrowing(
    cone: &Cone,
    circle: &Circle,
    axis: &Vector3,
    dec: &mut Decider,
) -> Result<Option<Real>, Refusal> {
    if !zero(&cone.axis.cross(axis).norm_squared(), dec)? {
        return Ok(None);
    }
    let apex = cone.apex.minus(&circle.center);
    if !zero(&apex.cross(axis).norm_squared(), dec)? {
        return Ok(None);
    }
    // The apex lies below the plane, so the cone narrows into the material.
    let depth = -&apex.dot(axis);
    if !positive(&depth, dec)? {
        return Ok(None);
    }
    let slope = &real_sqrt(&circle.radius_sq) / &depth;
    // The stored slope must give the mouth radius at the plane.
    let stored = &cone.slope / &real_sqrt(&cone.axis.norm_squared());
    if !zero(&(&stored - &slope), dec)? && !zero(&(&cone.slope - &slope), dec)? {
        return Ok(None);
    }
    Ok(Some(slope))
}

// Orders coedges into one closed loop. A fixed sense is kept; a free one is
// chosen so the loop continues.
fn chain(
    model: &Model,
    mut items: Vec<(EdgeId, Option<Sense>)>,
) -> Result<Vec<(EdgeId, Sense)>, Refusal> {
    let ends = |edge: EdgeId, sense: Sense| {
        let e = model.edge(edge);
        if sense == Sense::Same {
            (e.start, e.end)
        } else {
            (e.end, e.start)
        }
    };
    let first = items
        .iter()
        .position(|(_, s)| s.is_some())
        .ok_or_else(|| bad("a continued loop has no oriented edge"))?;
    let (edge, sense) = items.remove(first);
    let sense = sense.expect("fixed sense");
    let mut path = vec![(edge, sense)];
    let start = ends(edge, sense).0;
    let mut current = ends(edge, sense).1;
    while !items.is_empty() {
        let next = items.iter().position(|(e, s)| match s {
            Some(s) => ends(*e, *s).0 == current,
            None => ends(*e, Sense::Same).0 == current || ends(*e, Sense::Reversed).0 == current,
        });
        let Some(next) = next else {
            return Err(bad("continued boundary pieces do not close a loop"));
        };
        let (e, s) = items.remove(next);
        let s = s.unwrap_or(if ends(e, Sense::Same).0 == current {
            Sense::Same
        } else {
            Sense::Reversed
        });
        current = ends(e, s).1;
        path.push((e, s));
    }
    if current != start {
        return Err(bad("continued boundary loop is not closed"));
    }
    Ok(path)
}
