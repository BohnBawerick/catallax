//! Conservative finite face bounds used only to prove irrelevant support pairs.
//! A strictly nonzero support enclosure excludes an intersection; a box overlap
//! never claims a meeting or changes the analytic intersection method.

use crate::{
    geometry::sign,
    region::{Location, locate},
};
use catallax_numbers::{Ball, Decider, FactKind, Rational, Real, Refusal, Sign};
use catallax_shapes::{Circle, Curve, GeomKey, Plane, Point3, PointDef, Surface};
use catallax_topology::{FaceId, Loop, Model};
use std::collections::BTreeMap;

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

fn circle_box(circle: &Circle, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    let normal_squared = circle.normal.norm_squared();
    let coordinates = [&circle.center.x, &circle.center.y, &circle.center.z];
    let normal = [&circle.normal.x, &circle.normal.y, &circle.normal.z];
    let mut bounds = Vec::with_capacity(3);
    for axis in 0..3 {
        let factor = &Real::one() - &(&normal[axis].square() / &normal_squared);
        let width = if factor.rational_ref().is_some_and(Rational::is_zero) {
            Real::zero()
        } else {
            (&circle.radius_sq * &factor).sqrt()
        };
        let width = enclose(&width, dec)?;
        let center = enclose(coordinates[axis], dec)?;
        bounds.push((&center - &width).hull(&(&center + &width)));
    }
    Ok(bounds.try_into().expect("three coordinates"))
}

fn planar_box(
    model: &mut Model,
    face: FaceId,
    plane: &Plane,
    dec: &mut Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    let mut bounds: Option<[Ball; 3]> = None;
    for lp in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Ok(None);
        };
        for coedge in coedges {
            let edge = model.edge(model.coedge(*coedge).edge);
            let next = match model.geom().curve(edge.curve) {
                Curve::Circle(circle) => circle_box(circle, dec)?,
                Curve::Line(_) => {
                    let (Some(start), Some(end)) = (edge.start, edge.end) else {
                        return Ok(None);
                    };
                    let a = model.geom().position(model.vertex(start).point);
                    let b = model.geom().position(model.vertex(end).point);
                    let mut result = Vec::with_capacity(3);
                    for (a, b) in [(&a.x, &b.x), (&a.y, &b.y), (&a.z, &b.z)] {
                        result.push(enclose(a, dec)?.hull(&enclose(b, dec)?));
                    }
                    result.try_into().expect("three coordinates")
                }
                _ => return Ok(None),
            };
            if let Some(bounds) = &mut bounds {
                for i in 0..3 {
                    bounds[i] = bounds[i].hull(&next[i]);
                }
            } else {
                bounds = Some(next);
            }
        }
    }
    let Some(bounds) = bounds else {
        return Ok(None);
    };
    // All boundaries lie in this box. Prove the unbounded planar component is
    // outside using a point farther from the origin than the entire box. Every
    // exterior point is in that component, so every retained cell is bounded.
    let mut radius = Rational::zero();
    for (axis, coordinate) in [&plane.origin.x, &plane.origin.y, &plane.origin.z]
        .into_iter()
        .enumerate()
    {
        let origin = enclose(coordinate, dec)?;
        radius = radius
            + (&bounds[axis].lo() - &origin.hi())
                .abs()
                .max((&bounds[axis].hi() - &origin.lo()).abs());
    }
    let norm = plane.ref_dir.norm_squared();
    let length = if norm.rational_ref() == Some(&Rational::one()) {
        Real::one()
    } else {
        norm.sqrt()
    };
    let sample = plane.origin.plus(
        &plane
            .ref_dir
            .scaled(&(&Real::from(radius + Rational::one()) / &length)),
    );
    let point = model.geom_mut().add_point(PointDef::Explicit(sample));
    if locate(model, face, point, dec)? != Location::Outside {
        return Ok(None);
    }
    Ok(Some(bounds))
}

fn face_box(
    model: &mut Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    let surface = model.geom().surface(model.face(face).surface).clone();
    let (origin, axis, reference) = match &surface {
        Surface::Plane(plane) => return planar_box(model, face, plane, dec),
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
                let same_apex = !surface.has_rational_definition()
                    && p.x.same_definition(&origin.x)
                    && p.y.same_definition(&origin.y)
                    && p.z.same_definition(&origin.z);
                if !matches!(&surface, Surface::Cone(_))
                    || (!same_apex && !zero(&p.minus(origin).norm_squared(), dec)?)
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
                // Incidence on a circular carrier and coaxiality already prove
                // this entire ring is its section, including the retained nappe.
                // Keep the cheap rational path and avoid eliminating algebraic
                // construction identities that the producer has recorded.
                let constructed_section = !surface.has_rational_definition()
                    && [FactKind::Incident, FactKind::Coaxial]
                        .into_iter()
                        .all(|kind| {
                            model.geom().facts().holds(
                                kind,
                                GeomKey::Curve(edge.curve),
                                GeomKey::Surface(model.face(face).surface),
                            )
                        });
                if constructed_section {
                    dec.note_fact_decision();
                } else {
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

#[derive(Default)]
pub(crate) struct FaceBounds(BTreeMap<FaceId, Option<[Ball; 3]>>);

impl FaceBounds {
    pub(crate) fn disjoint(
        &mut self,
        model: &mut Model,
        a: FaceId,
        b: FaceId,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for (bounded, other) in [(a, b), (b, a)] {
            if let std::collections::btree_map::Entry::Vacant(entry) = self.0.entry(bounded) {
                entry.insert(face_box(model, bounded, dec)?);
            }
            if let Some(bounds) = self.0.get(&bounded).expect("cached face bounds") {
                let Some(polynomial) = model.geom().surface(model.face(other).surface).implicit()
                else {
                    continue;
                };
                let value = polynomial.eval_ball(bounds, dec)?;
                if !value.contains_zero() {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Pose, Vector3};
    use catallax_topology::{Face, Sense, cylinder};

    #[test]
    fn planar_bounds_exclude_only_supports_outside_the_finite_disk() {
        let r = Rational::from_i64;
        let mut model = Model::new();
        let body = cylinder(&mut model, r(2), r(3), &Pose::identity()).unwrap();
        let cap = model
            .faces_of_body(body)
            .into_iter()
            .find(|f| {
                matches!(
                    model.geom().surface(model.face(*f).surface),
                    Surface::Plane(_)
                )
            })
            .unwrap();
        let mut dec = Decider::default();
        let mut bounds = FaceBounds::default();
        for (x, excluded) in [(3, true), (2, false), (1, false)] {
            let surface = model.geom_mut().add_surface(Surface::Plane(Plane {
                origin: Point3::new(Real::from_i64(x), Real::zero(), Real::zero()),
                normal: Vector3::new(Real::one(), Real::zero(), Real::zero()),
                ref_dir: Vector3::new(Real::zero(), Real::one(), Real::zero()),
            }));
            let other = model.add_face(Face {
                surface,
                sense: Sense::Same,
                loops: vec![],
            });
            assert_eq!(
                bounds.disjoint(&mut model, cap, other, &mut dec).unwrap(),
                excluded
            );
        }
    }
}
