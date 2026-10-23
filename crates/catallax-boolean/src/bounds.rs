//! Conservative finite face bounds used only to prove irrelevant support pairs.
//! A strictly nonzero support enclosure excludes an intersection; a box overlap
//! never claims a meeting or changes the analytic intersection method.

use crate::{
    geometry::sign,
    region::{Location, locate},
};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Curve, Point3, PointDef, Surface};
use catallax_topology::{FaceId, Loop, Model};

fn enclose(x: &Real, dec: &mut Decider) -> Result<Ball, Refusal> {
    if let Some(q) = x.rational_ref() {
        Ok(Ball::from_rational(q, 128))
    } else {
        dec.enclose_budgeted_at(x, 128)
    }
}
fn sphere_box(center: &Point3, radius: &Real, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    let radius = enclose(radius, dec)?;
    let mut result = Vec::with_capacity(3);
    for x in [&center.x, &center.y, &center.z] {
        let c = enclose(x, dec)?;
        result.push((&c - &radius).hull(&(&c + &radius)));
    }
    Ok(result.try_into().expect("three coordinates"))
}
fn zero(x: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(x, dec, "boolean.bound-ambiguous")? == Sign::Zero)
}
fn face_box(
    model: &mut Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    let surface = model.geom().surface(model.face(face).surface).clone();
    let (origin, axis, reference) = match &surface {
        Surface::Sphere(s) => return sphere_box(&s.center, &s.radius, dec).map(Some),
        Surface::Torus(s) => return sphere_box(&s.center, &(&s.major + &s.minor), dec).map(Some),
        Surface::Cylinder(s) => (&s.origin, &s.axis, &s.ref_dir),
        Surface::Cone(s) => (&s.apex, &s.axis, &s.ref_dir),
        _ => return Ok(None),
    };
    let mut ends = Vec::with_capacity(2);
    for l in &model.face(face).loops {
        match model.loop_(*l) {
            Loop::Vertex(v) => {
                let PointDef::Explicit(p) = model.geom().point(model.vertex(*v).point) else {
                    return Ok(None);
                };
                if !matches!(&surface, Surface::Cone(_))
                    || !zero(&p.minus(origin).norm_squared(), dec)?
                {
                    return Ok(None);
                }
                ends.push((p.clone(), Real::zero()));
            }
            Loop::Coedges(cs) if cs.len() == 1 => {
                let edge = model.edge(model.coedge(cs[0]).edge);
                if edge.start.is_some() || edge.end.is_some() {
                    return Ok(None);
                }
                let Curve::Circle(c) = model.geom().curve(edge.curve) else {
                    return Ok(None);
                };
                if !zero(&c.center.minus(origin).cross(axis).norm_squared(), dec)?
                    || !zero(&c.normal.cross(axis).norm_squared(), dec)?
                {
                    return Ok(None);
                }
                let expected = match &surface {
                    Surface::Cylinder(s) => s.radius_sq.clone(),
                    Surface::Cone(s) => {
                        let axial = c.center.minus(&s.apex).dot(&s.axis);
                        if sign(&axial, dec, "boolean.bound-ambiguous")? == Sign::Negative {
                            return Ok(None);
                        }
                        &s.slope.square() * &(&axial.square() / &s.axis.norm_squared())
                    }
                    _ => unreachable!(),
                };
                if !zero(&(&expected - &c.radius_sq), dec)? {
                    return Ok(None);
                }
                let radius = if let Some(q) = c.radius_sq.rational_ref() {
                    if *q == Rational::one() {
                        Real::one()
                    } else {
                        c.radius_sq.sqrt()
                    }
                } else {
                    c.radius_sq.sqrt()
                };
                ends.push((c.center.clone(), radius));
            }
            _ => return Ok(None),
        }
    }
    if ends.len() != 2 {
        return Ok(None);
    }
    // Two complete coaxial circles cut a cylinder/cone nappe into axial bands.
    // Proving the middle band is the trimmed face excludes the unbounded cells.
    let half = Real::from(Rational::from_ratio(1, 2).expect("fixed denominator"));
    let center = ends[0].0.plus(&ends[1].0.minus(&ends[0].0).scaled(&half));
    let radius = &(&ends[0].1 + &ends[1].1) * &half;
    let norm = reference.norm_squared();
    let length = if norm.rational_ref() == Some(&Rational::one()) {
        Real::one()
    } else {
        norm.sqrt()
    };
    let sample = center.plus(&reference.scaled(&(&radius / &length)));
    let point = model.geom_mut().add_point(PointDef::Explicit(sample));
    if locate(model, face, point, dec)? != Location::Inside {
        return Ok(None);
    }
    let a = sphere_box(&ends[0].0, &ends[0].1, dec)?;
    let b = sphere_box(&ends[1].0, &ends[1].1, dec)?;
    Ok(Some(std::array::from_fn(|i| a[i].hull(&b[i]))))
}

pub(crate) fn disjoint(
    model: &mut Model,
    a: FaceId,
    b: FaceId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for (bounded, other) in [(a, b), (b, a)] {
        if let Some(bounds) = face_box(model, bounded, dec)? {
            let polynomial = model.geom().surface(model.face(other).surface).implicit();
            let value = polynomial.eval_ball(&bounds, dec)?;
            if !value.contains_zero() {
                return Ok(true);
            }
        }
    }
    Ok(false)
}
