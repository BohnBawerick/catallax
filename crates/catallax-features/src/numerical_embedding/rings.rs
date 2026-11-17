//! Explicit numerical vertices and seam cells for retained circular ring bands.

use super::*;
use catallax_shapes::{Circle, Line};

pub(super) fn is_ring(model: &Model, id: EdgeId) -> bool {
    let edge = model.edge(id);
    edge.start.is_none()
        && edge.end.is_none()
        && matches!(model.geom().curve(edge.curve), Curve::Circle(_))
}

pub(super) fn vertices(
    model: &Model,
    body: BodyId,
    vertices: &mut Vec<NumericalVertex>,
    approx: &mut Vec<[f64; 3]>,
    dec: &mut Decider,
) -> Result<BTreeMap<EdgeId, usize>, Refusal> {
    let mut rings = BTreeMap::new();
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
        // Coaxial rings choose the same world meridian, independently of their
        // stored circle frames and normal signs.
        let axis = usize::from(n[0].abs() >= 0.9);
        let radial = std::array::from_fn(|i| (if i == axis { 1.0 } else { 0.0 }) - n[axis] * n[i]);
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
        return Err(unsupported("a complete native ring band needs two circular boundaries"));
    };
    let ([a], [b]) = (first.as_slice(), second.as_slice()) else {
        return Err(unsupported("a complete native ring band needs single-ring boundaries"));
    };
    if !closed(a) || !closed(b) {
        return Err(unsupported("both complete band boundaries must be circular rings"));
    }
    let mut a = a.clone();
    let mut b = b.clone();
    let direction = angular_direction(&face.surface, &a, &edges[a.edge], dec)?;
    let other = angular_direction(&face.surface, &b, &edges[b.edge], dec)?;
    if !direction.is_finite() || !other.is_finite() || direction * other >= 0.0 {
        return Err(fit_failure("complete band rings have inconsistent angular traversal"));
    }
    let period = Ball::pi(BITS).mid().mul_pow2(1);
    let turn = if direction > 0.0 { period } else { -period };
    a.uv[1][0] = &a.uv[0][0] + &turn;
    b.uv[0][0] = a.uv[1][0].clone();
    b.uv[1][0] = a.uv[0][0].clone();
    let ends = [edges[a.edge].vertices[0], edges[b.edge].vertices[0]];
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
            (Curve::Circle(circle), [a.uv[1][1].clone(), b.uv[0][1].clone()])
        }
        _ => return Err(unsupported("complete rings need a regular cylinder, cone or torus band")),
    };
    let seam = edges.len();
    edges.push(NumericalEdge {
        key: format!("{}:periodic-seam", face.key), curve, vertices: ends, parameters,
    });
    let forward = NumericalTrim {
        edge: seam, sense: Sense::Same, uv: [a.uv[1].clone(), b.uv[0].clone()],
    };
    let reverse = NumericalTrim {
        edge: seam, sense: Sense::Reversed, uv: [b.uv[1].clone(), a.uv[0].clone()],
    };
    face.loops = vec![vec![a, forward, b, reverse]];
    Ok(())
}
