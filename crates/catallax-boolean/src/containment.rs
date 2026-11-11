//! Exact rays against trimmed faces. Boundary-edge hits never choose a side.

use crate::Containment;
use crate::geometry::{r, refusal, sign, undecided};
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
/// Whether the query point lies on `face`, boundary included.
fn on_face(
    model: &mut Model,
    face: FaceId,
    query: catallax_shapes::PointId,
    point: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let surface = model.face(face).surface;
    if matches!(
        model.geom().surface(surface),
        Surface::BSpline(catallax_shapes::SplineSupport::BezierExtrusion(_))
    ) {
        return Ok(crate::geometry::on_surface(
            model,
            surface,
            query,
            dec,
            "boolean.classify-ambiguous",
        )? && locate(model, face, query, dec)? != Location::Outside);
    }
    // A proved miss of a finite support owes no global side of it.
    if let Surface::BSpline(support) = model.geom().surface(surface)
        && support.misses_point(point, dec)?
    {
        return Ok(false);
    }
    match model.geom().side(surface, query, dec) {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(_) => Ok(locate(model, face, query, dec)? != Location::Outside),
        Outcome::NotExact(n) => Err(undecided("boolean.classify-ambiguous", n)),
        Outcome::Refused(error) => Err(error),
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
    // Each face is tested on its own; workers test later faces ahead.
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("point sides");
    let mut ahead = catallax_numbers::run_ahead!(dec, faces.len());
    for (index, &face) in faces.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        ahead.at(index, || {
            let (model, faces, point) = (model.clone(), faces.to_vec(), point.clone());
            move |index: usize, dec: &mut Decider| {
                let _ = on_face(&mut model.clone(), faces[index], query, &point, dec);
            }
        });
        if on_face(model, face, query, point, dec)? {
            return Ok(None);
        }
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
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
/// How many times the ray on line `curve` crosses the inside of `face`, or
/// `None` when the ray lies in the face's surface or meets its boundary.
fn ray_crossings(
    model: &mut Model,
    curve: catallax_shapes::CurveId,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<usize>, Refusal> {
    let surface = model.face(face).surface;
    let answer = intersect_curve_surface(model.geom_mut(), curve, surface, dec)?;
    if !answer.marks.is_empty() {
        return Err(refusal(
            "boolean.classify-ambiguous",
            "a ray intersection required an unresolved decision",
        ));
    }
    let hits = match answer.value {
        CurveSurface::Disjoint => return Ok(Some(0)),
        CurveSurface::Contained => return Ok(None),
        CurveSurface::Meets(hits) => hits,
    };
    let Curve::Line(line) = model.geom().curve(curve) else {
        unreachable!("parity ray is stored as a line");
    };
    let (origin, direction) = (line.origin.clone(), line.dir.clone());
    let mut crossings = 0usize;
    for hit in hits {
        let displacement = model.geom().position(hit.point).minus(&origin);
        if sign(
            &displacement.dot(&direction),
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
    Ok(Some(crossings))
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
    // Each face's crossings are counted on their own; workers count later
    // faces ahead.
    #[cfg(feature = "work-span")]
    let span_tasks = catallax_numbers::work_span::parallel("ray crossings");
    let mut ahead = catallax_numbers::run_ahead!(dec, faces.len());
    for (index, &face) in faces.iter().enumerate() {
        #[cfg(feature = "work-span")]
        let _span_task = span_tasks.task();
        ahead.at(index, || {
            let (model, faces) = (model.clone(), faces.to_vec());
            move |index: usize, dec: &mut Decider| {
                let _ = ray_crossings(&mut model.clone(), curve, faces[index], dec);
            }
        });
        let Some(found) = ray_crossings(model, curve, face, dec)? else {
            return Ok(None);
        };
        crossings += found;
    }
    drop(ahead);
    #[cfg(feature = "work-span")]
    drop(span_tasks);
    Ok(Some(crossings % 2 == 1))
}
// A complete checked prism has an exact closed form, with no ray to choose.
fn prism(
    model: &Model,
    body: BodyId,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Option<Outcome<Containment>>, Refusal> {
    if let Some(cell) = catallax_topology::DraftedCylinder::of_body(model, body, dec)? {
        return Ok(Some(match cell.locate(point, dec)? {
            Sign::Negative => Outcome::Proven(Containment::Inside),
            Sign::Positive => Outcome::Proven(Containment::Outside),
            Sign::Zero => Outcome::ProvenCoincident(Coincidence {
                proof: ZeroProof::ExactArithmetic,
            }),
        }));
    }
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
