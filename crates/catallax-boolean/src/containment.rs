//! Exact rays against trimmed faces. Boundary-edge hits never choose a side.

use crate::Containment;
use crate::geometry::{r, refusal, sign};
use crate::region::{Location, locate};
use catallax_intersect::{Contact, CurveSurface, intersect_curve_surface};
use catallax_numbers::{Coincidence, Decider, Outcome, Refusal, Sign, ZeroProof};
use catallax_shapes::{Curve, Line, Point3, PointDef, RVector3, Surface, Vector3};
use catallax_topology::{BodyId, BodyKind, FaceId, Model, SplinePrism};

pub(crate) fn faces(
    model: &mut Model,
    faces: &[FaceId],
    point: &Point3,
    dec: &mut Decider,
) -> Outcome<Containment> {
    match classify(model, faces, point, dec) {
        Ok(Some(value)) => Outcome::Proven(value),
        Ok(None) => Outcome::ProvenCoincident(Coincidence {
            proof: ZeroProof::ExactArithmetic,
        }),
        Err(error) => Outcome::Refused(error),
    }
}
fn classify(
    model: &mut Model,
    faces: &[FaceId],
    point: &Point3,
    dec: &mut Decider,
) -> Result<Option<Containment>, Refusal> {
    if faces.is_empty() {
        return Ok(Some(Containment::Outside));
    }
    let query = model
        .geom_mut()
        .add_point(PointDef::Explicit(point.clone()));
    for &face in faces {
        let surface = model.face(face).surface;
        // A proved miss of a finite support owes no global side of it.
        if let Surface::BSpline(support) = model.geom().surface(surface)
            && support.misses_point(point, dec)?
        {
            continue;
        }
        match model.geom().side(surface, query, dec) {
            Outcome::Proven(_) => {}
            Outcome::ProvenCoincident(_) => {
                if locate(model, face, query, dec)? != Location::Outside {
                    return Ok(None);
                }
            }
            Outcome::NotExact(n) => return Err(refusal("boolean.classify-ambiguous", n.detail)),
            Outcome::Refused(error) => return Err(error),
        }
    }
    // An exact cap normal avoids unnecessary high-height side intersections.
    let normal = faces.iter().find_map(|&face| {
        let surface = model.face(face).surface;
        match model.geom().surface(surface) {
            Surface::Plane(plane) => Some((surface, plane.normal.clone())),
            _ => None,
        }
    });
    let mut rejected_direction = None;
    if let Some((surface, direction)) = normal {
        if let Some(inside) = ray_parity_owned(model, faces, point, direction, dec)? {
            return Ok(Some(if inside {
                Containment::Inside
            } else {
                Containment::Outside
            }));
        }
        let Surface::Plane(plane) = model.geom().surface(surface) else {
            unreachable!("selected plane does not change");
        };
        // This is duplicate-work suppression, not a geometry predicate. An
        // algebraic direction with no cheap rational key simply remains eligible.
        if let Some([x, y, z]) = plane.normal.rational_components()
            && x.sign() == Sign::Positive
        {
            rejected_direction = Some((y / x, z / x));
        }
    }
    for choice in 0..crate::geometry::DIRECTION_ATTEMPTS {
        let (a, b) = crate::geometry::probe_parameters(choice);
        if rejected_direction
            .as_ref()
            .is_some_and(|(y, z)| &a == y && &b == z)
        {
            continue;
        }
        let direction = RVector3::new(r(1), a, b).to_real();
        if let Some(inside) = ray_parity_owned(model, faces, point, direction, dec)? {
            return Ok(Some(if inside {
                Containment::Inside
            } else {
                Containment::Outside
            }));
        }
    }
    Err(refusal(
        "boolean.classify-ambiguous",
        "no rational ray avoided the face boundaries",
    ))
}
pub(crate) fn ray_parity(
    model: &mut Model,
    faces: &[FaceId],
    point: &Point3,
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    ray_parity_owned(model, faces, point, direction.clone(), dec)
}
fn ray_parity_owned(
    model: &mut Model,
    faces: &[FaceId],
    point: &Point3,
    direction: Vector3,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let curve = model
        .geom_mut()
        .add_curve(Curve::Line(Line::new(point.clone(), direction, dec)?));
    let mut crossings = 0usize;
    let mut ray_frame = None;
    for &face in faces {
        let surface = model.face(face).surface;
        let answer = intersect_curve_surface(model.geom_mut(), curve, surface, dec)?;
        if !answer.marks.is_empty() {
            return Err(refusal(
                "boolean.classify-ambiguous",
                "a ray intersection required an unresolved decision",
            ));
        }
        let hits = match answer.value {
            CurveSurface::Disjoint => continue,
            CurveSurface::Contained => return Ok(None),
            CurveSurface::Meets(hits) => hits,
        };
        for hit in hits {
            let (origin, direction) = ray_frame.get_or_insert_with(|| {
                let Curve::Line(line) = model.geom().curve(curve) else {
                    unreachable!("parity ray is stored as a line");
                };
                (line.origin.clone(), line.dir.clone())
            });
            let displacement = model.geom().position(hit.point).minus(origin);
            if sign(
                &displacement.dot(direction),
                dec,
                "boolean.classify-ambiguous",
            )? != Sign::Positive
            {
                continue;
            }
            match locate(model, face, hit.point, dec)? {
                Location::Outside => {}
                Location::Boundary => return Ok(None),
                Location::Inside => {
                    if hit.contact == Contact::Crossing {
                        crossings += 1;
                    }
                }
            }
        }
    }
    Ok(Some(crossings % 2 == 1))
}
// A complete checked prism has an exact closed form, with no ray to choose.
fn prism(
    model: &Model,
    body: BodyId,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Option<Outcome<Containment>>, Refusal> {
    let Some(prism) = SplinePrism::of_body(model, body, dec)? else {
        return Ok(None);
    };
    Ok(Some(match prism.locate(point, dec)? {
        Sign::Negative => Outcome::Proven(Containment::Inside),
        Sign::Positive => Outcome::Proven(Containment::Outside),
        Sign::Zero => Outcome::ProvenCoincident(Coincidence {
            proof: ZeroProof::ExactArithmetic,
        }),
    }))
}
pub(crate) fn classify_point(
    model: &mut Model,
    body: BodyId,
    point: &Point3,
    dec: &mut Decider,
) -> Outcome<Containment> {
    match model.body(body).kind {
        BodyKind::Empty => Outcome::Proven(Containment::Outside),
        BodyKind::Solid => match prism(model, body, point, dec) {
            Ok(Some(answer)) => answer,
            Ok(None) => faces(model, &model.faces_of_body(body), point, dec),
            Err(error) => Outcome::Refused(error),
        },
        _ => Outcome::Refused(refusal(
            "boolean.not-solid",
            "point containment requires a solid or empty body",
        )),
    }
}
