//! Rational rays against trimmed faces. Rays through boundary edges are
//! discarded, so neither vertex multiplicity nor a seam chooses a side.

use catallax_intersect::{Contact, CurveSurface, intersect_curve_surface};
use catallax_numbers::{Coincidence, Decider, Outcome, Refusal, Sign, ZeroProof};
use catallax_shapes::{Curve, Line, PointDef, RPoint3, RVector3, Surface};
use catallax_topology::{BodyId, BodyKind, FaceId, Model};

use crate::Containment;
use crate::geometry::{r, refusal, sign};
use crate::region::{Location, locate};

pub(crate) fn faces(
    model: &mut Model,
    faces: &[FaceId],
    point: &RPoint3,
    dec: &mut Decider,
) -> Outcome<Containment> {
    match classify(model, faces, point, dec) {
        Ok(Some(value)) => Outcome::Proven(value),
        Ok(None) => Outcome::ProvenCoincident(Coincidence {
            proof: ZeroProof::ExactArithmetic,
        }),
        Err(e) => Outcome::Refused(e),
    }
}

fn classify(
    model: &mut Model,
    faces: &[FaceId],
    point: &RPoint3,
    dec: &mut Decider,
) -> Result<Option<Containment>, Refusal> {
    if faces.is_empty() {
        return Ok(Some(Containment::Outside));
    }
    let query = model
        .geom_mut()
        .add_point(PointDef::Explicit(point.clone()));
    for &face in faces {
        match model.geom().side(model.face(face).surface, query, dec) {
            Outcome::Proven(_) => {}
            Outcome::ProvenCoincident(_) => {
                if locate(model, face, query, dec)? != Location::Outside {
                    return Ok(None);
                }
            }
            Outcome::NotExact(n) => return Err(refusal("boolean.classify-ambiguous", n.detail)),
            Outcome::Refused(e) => return Err(e),
        }
    }
    // Start in the body's rational frame. A cap normal is axial for its
    // cylinder, avoiding an unnecessary high-height quadratic on the side.
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
            unreachable!("the selected plane definition does not change");
        };
        // World directions have X=1; only positive multiples repeat this half-ray.
        if plane.normal.x.sign() == Sign::Positive {
            rejected_direction = Some((
                &plane.normal.y / &plane.normal.x,
                &plane.normal.z / &plane.normal.x,
            ));
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
        let direction = RVector3::new(r(1), a, b);
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
    point: &RPoint3,
    direction: &RVector3,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    ray_parity_owned(model, faces, point, direction.clone(), dec)
}

fn ray_parity_owned(
    model: &mut Model,
    faces: &[FaceId],
    point: &RPoint3,
    direction: RVector3,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let curve = model
        .geom_mut()
        .add_curve(Curve::Line(Line::new(point.clone(), direction)?));
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
                    unreachable!("a parity ray is inserted as a line");
                };
                (line.origin.to_real(), line.dir.to_real())
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

pub(crate) fn classify_point(
    model: &mut Model,
    body: BodyId,
    point: &RPoint3,
    dec: &mut Decider,
) -> Outcome<Containment> {
    match model.body(body).kind {
        BodyKind::Empty => Outcome::Proven(Containment::Outside),
        BodyKind::Solid => faces(model, &model.faces_of_body(body), point, dec),
        _ => Outcome::Refused(refusal(
            "boolean.not-solid",
            "point containment requires a solid or empty body",
        )),
    }
}
