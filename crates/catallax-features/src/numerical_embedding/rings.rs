//! Explicit numerical vertices and seam cells for retained circular ring bands.

use super::*;
use catallax_shapes::{Circle, Line};

pub(super) fn is_ring(model: &Model, id: EdgeId) -> bool {
    let edge = model.edge(id);
    edge.start.is_none()
        && edge.end.is_none()
        && matches!(model.geom().curve(edge.curve), Curve::Circle(_))
}

fn root(parents: &mut [usize], mut index: usize) -> usize {
    while parents[index] != index {
        parents[index] = parents[parents[index]];
        index = parents[index];
    }
    index
}

fn meridians(model: &Model, body: BodyId, dec: &mut Decider) -> Result<BTreeMap<EdgeId, [f64; 3]>, Refusal> {
    let ids: Vec<_> = model.edges_of_body(body).into_iter().filter(|id| is_ring(model, *id)).collect();
    if ids.is_empty() {
        return Ok(BTreeMap::new());
    }
    let index: BTreeMap<_, _> = ids.iter().enumerate().map(|(i, id)| (*id, i)).collect();
    let mut parents: Vec<_> = (0..ids.len()).collect();
    let mut anchors = Vec::new();
    for id in model.faces_of_body(body) {
        let face = model.face(id);
        if !matches!(model.geom().surface(face.surface), Surface::Cylinder(_) | Surface::Cone(_) | Surface::Torus(_)) {
            continue;
        }
        let mut members = Vec::new();
        let mut anchor = None;
        for lp in &face.loops {
            let Loop::Coedges(uses) = model.loop_(*lp) else { continue };
            for coedge in uses {
                dec.charge_work(1)?;
                let coedge = model.coedge(*coedge);
                if let Some(member) = index.get(&coedge.edge) {
                    members.push(*member);
                } else if anchor.is_none() {
                    let edge = model.edge(coedge.edge);
                    let vertex = if coedge.sense == Sense::Same { edge.start } else { edge.end };
                    if let Some(vertex) = vertex {
                        anchor = Some(f64_point(&model.geom().position(model.vertex(vertex).point), dec)?);
                    }
                }
            }
        }
        let Some(first) = members.first().copied() else { continue };
        for member in &members[1..] {
            let a = root(&mut parents, first);
            let b = root(&mut parents, *member);
            parents[b] = a;
        }
        if let Some(anchor) = anchor { anchors.push((first, anchor)); }
    }
    let mut groups = BTreeMap::new();
    for (member, anchor) in anchors {
        groups.entry(root(&mut parents, member)).or_insert(anchor);
    }
    Ok(ids.into_iter().enumerate().filter_map(|(member, id)| {
        groups.get(&root(&mut parents, member)).map(|anchor| (id, *anchor))
    }).collect())
}

pub(super) fn vertices(
    model: &Model,
    body: BodyId,
    vertices: &mut Vec<NumericalVertex>,
    approx: &mut Vec<[f64; 3]>,
    dec: &mut Decider,
) -> Result<BTreeMap<EdgeId, usize>, Refusal> {
    let mut rings = BTreeMap::new();
    let meridians = meridians(model, body, dec)?;
    for id in model.edges_of_body(body) {
        if !is_ring(model, id) {
            continue;
        }
        dec.charge_work(1)?;
        let Curve::Circle(circle) = model.geom().curve(model.edge(id).curve) else {
            unreachable!("circular ring")
        };
        let frame = Circle64::new(circle, dec)?;
        let n = f64_point(
            &Point3::new(circle.normal.x.clone(), circle.normal.y.clone(), circle.normal.z.clone()),
            dec,
        )?;
        let n = n.map(|v| v / dot(n, n).sqrt());
        // A finite rim chooses the band's seam. Otherwise coaxial rings use
        // one common world meridian independently of their stored frames.
        let delta = match meridians.get(&id) {
            Some(anchor) => sub(*anchor, frame.center),
            None => {
                let axis = usize::from(n[0].abs() >= 0.9);
                std::array::from_fn(|i| if i == axis { 1.0 } else { 0.0 })
            }
        };
        let axial = dot(delta, n);
        let radial = std::array::from_fn(|i| delta[i] - axial * n[i]);
        let phase = dot(radial, frame.y).atan2(dot(radial, frame.x));
        let point = frame.point(phase);
        let position = [exact(point[0])?, exact(point[1])?, exact(point[2])?];
        rings.insert(id, vertices.len());
        vertices.push(NumericalVertex { key: format!("{id:?}:ring-vertex"), position });
        approx.push(point);
    }
    Ok(rings)
}

fn position(vertex: &NumericalVertex) -> Point3 {
    let [x, y, z] = vertex.position.clone();
    Point3::new(x.into(), y.into(), z.into())
}

fn angular_direction(
    surface: &Surface,
    trim: &NumericalTrim,
    edge: &NumericalEdge,
    dec: &mut Decider,
) -> Result<f64, Refusal> {
    let Curve::Circle(circle) = &edge.curve else {
        return Err(unsupported("a periodic band boundary must be a circular ring"));
    };
    let frame = Circle64::new(circle, dec)?;
    let parameter = if trim.sense == Sense::Same { &edge.parameters[0] } else { &edge.parameters[1] };
    let (sin, cos) = parameter.to_f64_nearest().sin_cos();
    let tangent = std::array::from_fn(|i| -sin * frame.x[i] + cos * frame.y[i]);
    let at = trim.uv[0].each_ref().map(|q| Ball::from_rational(q, BITS));
    let du = surface.derivatives(&at[0], &at[1], dec)?[0]
        .each_ref().map(|v| v.mid().to_f64_nearest());
    let travel = if trim.sense == Sense::Same { 1.0 } else { -1.0 };
    Ok(dot(tangent, du) * (&edge.parameters[1] - &edge.parameters[0]).to_f64_nearest() * travel)
}

pub(super) fn close_band(
    face: &mut NumericalFace,
    vertices: &[NumericalVertex],
    edges: &mut Vec<NumericalEdge>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if matches!(face.surface, Surface::Plane(_) | Surface::BSpline(_)) {
        return Ok(());
    }
    let closed = |trim: &NumericalTrim| {
        let edge = &edges[trim.edge];
        edge.vertices[0] == edge.vertices[1] && matches!(edge.curve, Curve::Circle(_))
    };
    if !face.loops.iter().flatten().any(closed) {
        return Ok(());
    }
    let [first, second] = face.loops.as_slice() else {
        return Err(unsupported("a complete native band needs two winding boundaries"));
    };
    let (a, mut rim) = match (first.as_slice(), second.as_slice()) {
        ([a], other) if closed(a) => (a.clone(), other.to_vec()),
        (other, [a]) if closed(a) => (a.clone(), other.to_vec()),
        _ => return Err(unsupported("a periodic band needs one retained circular ring at its seam")),
    };
    let mut a = a;
    let direction = angular_direction(&face.surface, &a, &edges[a.edge], dec)?;
    let period = Ball::pi(BITS).mid().mul_pow2(1);
    if let [b] = rim.as_mut_slice() && closed(b) {
        let other = angular_direction(&face.surface, b, &edges[b.edge], dec)?;
        if !direction.is_finite() || !other.is_finite() || direction * other >= 0.0 {
            return Err(fit_failure("complete band rings have inconsistent angular traversal"));
        }
        let turn = if direction > 0.0 { period.clone() } else { -&period };
        a.uv[1][0] = &a.uv[0][0] + &turn;
        b.uv[0][0] = a.uv[1][0].clone();
        b.uv[1][0] = a.uv[0][0].clone();
    } else {
        let first = rim.first().ok_or_else(|| unsupported("a periodic rim is empty"))?;
        let last = rim.last().unwrap();
        let turn = &last.uv[1][0] - &first.uv[0][0];
        if turn.abs() != period || direction * turn.to_f64_nearest() >= 0.0 {
            return Err(fit_failure("a finite periodic rim must wind once opposite its ring"));
        }
        a.uv[0][0] = last.uv[1][0].clone();
        a.uv[1][0] = first.uv[0][0].clone();
    }
    let first = rim.first().unwrap();
    let last = rim.last().unwrap();
    let rim_vertex = edges[first.edge].vertices[usize::from(first.sense == Sense::Reversed)];
    let rim_end = edges[last.edge].vertices[usize::from(last.sense == Sense::Same)];
    if rim_vertex != rim_end {
        return Err(fit_failure("a periodic rim does not close at its seam vertex"));
    }
    let ends = [edges[a.edge].vertices[0], rim_vertex];
    let from = position(&vertices[ends[0]]);
    let to = position(&vertices[ends[1]]);
    let (curve, parameters) = match &face.surface {
        Surface::Cylinder(_) | Surface::Cone(_) => (
            Curve::Line(Line::new(from.clone(), to.minus(&from), dec)?),
            [Rational::zero(), Rational::one()],
        ),
        Surface::Torus(torus) => {
            let delta = from.minus(&torus.center);
            let axial = torus.axis.scaled(&(delta.dot(&torus.axis) / torus.axis.norm_squared()));
            let radial = delta.minus(&axial);
            let radial = radial.scaled(&(Real::one() / radial.norm_squared().sqrt()));
            let center = torus.center.plus(&radial.scaled(&torus.major));
            let normal = radial.cross(&torus.axis);
            let circle = Circle::new(center, normal, torus.minor.square(), radial, dec)?;
            (Curve::Circle(circle), [a.uv[1][1].clone(), first.uv[0][1].clone()])
        }
        _ => return Err(unsupported("complete rings need a regular cylinder, cone or torus band")),
    };
    let seam = edges.len();
    edges.push(NumericalEdge {
        key: format!("{}:periodic-seam", face.key), curve, vertices: ends, parameters,
    });
    let forward = NumericalTrim {
        edge: seam, sense: Sense::Same, uv: [a.uv[1].clone(), first.uv[0].clone()],
    };
    let reverse = NumericalTrim {
        edge: seam, sense: Sense::Reversed, uv: [last.uv[1].clone(), a.uv[0].clone()],
    };
    let mut boundary = Vec::with_capacity(rim.len() + 3);
    boundary.push(a);
    boundary.push(forward);
    boundary.append(&mut rim);
    boundary.push(reverse);
    face.loops = vec![boundary];
    Ok(())
}
