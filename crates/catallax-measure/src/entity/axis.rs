//! Centre distances and angles from exact axes, centres and normals.

use super::{EdgeCurve, PairEntity, edge_curve, enclose, face_surface, unit, zero};
use crate::{Accuracy, Enclosure};
use catallax_numbers::{Answer, Answered, Ball, Decider, MarkLog, Rational, Real, Refusal, Sign};
use catallax_shapes::{Point3, Surface, Vector3};
use catallax_topology::Model;

/// The centre feature of an entity.
enum Centre {
    Point(Point3),
    Line(Point3, Vector3),
    Plane(Point3, Vector3),
}

/// A direction and whether it is the normal of a plane, circle or torus core
/// rather than an axis.
struct Direction {
    vector: Vector3,
    normal: bool,
}

fn centre(model: &Model, entity: PairEntity, dec: &mut Decider) -> Result<Centre, Refusal> {
    Ok(match entity {
        PairEntity::Face(face) => match face_surface(model, face)? {
            Surface::Plane(s) => Centre::Plane(s.origin.clone(), s.normal.clone()),
            Surface::Cylinder(s) => Centre::Line(s.origin.clone(), s.axis.clone()),
            Surface::Cone(s) => Centre::Line(s.apex.clone(), s.axis.clone()),
            Surface::Sphere(s) => Centre::Point(s.center.clone()),
            Surface::Torus(s) => Centre::Point(s.center.clone()),
            _ => unreachable!("face_surface admits only analytic carriers"),
        },
        PairEntity::Edge(edge) => match edge_curve(model, edge, dec)? {
            EdgeCurve::Segment(a, b) => {
                let d = b.minus(&a);
                Centre::Line(a, d)
            }
            EdgeCurve::Circle { center, .. } => Centre::Point(center),
        },
    })
}

fn direction(model: &Model, entity: PairEntity, dec: &mut Decider) -> Result<Direction, Refusal> {
    let axis = |vector: &Vector3| Direction {
        vector: vector.clone(),
        normal: false,
    };
    Ok(match entity {
        PairEntity::Face(face) => match face_surface(model, face)? {
            Surface::Plane(s) => Direction {
                vector: s.normal.clone(),
                normal: true,
            },
            Surface::Cylinder(s) => axis(&s.axis),
            Surface::Cone(s) => axis(&s.axis),
            Surface::Torus(s) => Direction {
                vector: s.axis.clone(),
                normal: true,
            },
            Surface::Sphere(_) => {
                return Err(Refusal::new(
                    "measure.angle-domain",
                    format!("face {face:?} is a sphere, which has no direction"),
                ));
            }
            _ => unreachable!("face_surface admits only analytic carriers"),
        },
        PairEntity::Edge(edge) => match edge_curve(model, edge, dec)? {
            EdgeCurve::Segment(a, b) => axis(&b.minus(&a)),
            EdgeCurve::Circle { normal, .. } => Direction {
                vector: normal,
                normal: true,
            },
        },
    })
}

fn not_parallel() -> Refusal {
    Refusal::new(
        "measure.centre-domain",
        "the centre line or plane meets the other plane, so no centre distance is defined",
    )
}

/// Squared distance between two centre features.
fn centre_square(a: &Centre, b: &Centre, dec: &mut Decider) -> Result<Real, Refusal> {
    let line_point =
        |p: &Point3, o: &Point3, d: &Vector3| p.minus(o).cross(d).norm_squared() / d.norm_squared();
    let plane_point =
        |p: &Point3, o: &Point3, n: &Vector3| p.minus(o).dot(n).square() / n.norm_squared();
    Ok(match (a, b) {
        (Centre::Point(p), Centre::Point(q)) => p.minus(q).norm_squared(),
        (Centre::Point(p), Centre::Line(o, d)) | (Centre::Line(o, d), Centre::Point(p)) => {
            line_point(p, o, d)
        }
        (Centre::Point(p), Centre::Plane(o, n)) | (Centre::Plane(o, n), Centre::Point(p)) => {
            plane_point(p, o, n)
        }
        (Centre::Line(o1, d1), Centre::Line(o2, d2)) => {
            let c = d1.cross(d2);
            let c2 = c.norm_squared();
            if zero(&c2, dec)? {
                line_point(o2, o1, d1)
            } else {
                o2.minus(o1).dot(&c).square() / c2
            }
        }
        (Centre::Line(o, d), Centre::Plane(p, n)) | (Centre::Plane(p, n), Centre::Line(o, d)) => {
            if !zero(&d.dot(n), dec)? {
                return Err(not_parallel());
            }
            plane_point(o, p, n)
        }
        (Centre::Plane(o1, n1), Centre::Plane(o2, n2)) => {
            if !zero(&n1.cross(n2).norm_squared(), dec)? {
                return Err(not_parallel());
            }
            plane_point(o2, o1, n1)
        }
    })
}

/// Distance between the entities' centres. Lines and planes that meet a plane
/// refuse `measure.centre-domain`.
pub(super) fn centre_distance(
    model: &Model,
    a: PairEntity,
    b: PairEntity,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    let (ca, cb) = (centre(model, a, dec)?, centre(model, b, dec)?);
    let square = centre_square(&ca, &cb, dec)?;
    enclose(&super::root(&square), accuracy, "centre distance", dec)
}

/// Unoriented acute angle in degrees. Squared cosines 0, 1/4, 1/2, 3/4 and 1
/// answer 90, 60, 45, 30 and 0 exactly; others are outward enclosures of the
/// principal arctangent.
pub(super) fn angle(
    model: &Model,
    a: PairEntity,
    b: PairEntity,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    let (da, db) = (direction(model, a, dec)?, direction(model, b, dec)?);
    let (u, v) = (unit(&da.vector), unit(&db.vector));
    let cos2 = u.dot(&v).square();
    // A line against a plane measures the complement of its angle to the normal.
    let complement = da.normal != db.normal;
    let exact = match cos2.rational_ref() {
        Some(q) => Some(q.clone()),
        None => dec.rational_value(&cos2)?,
    };
    if let Some(q) = exact {
        let table = [(0, 1, 90), (1, 4, 60), (1, 2, 45), (3, 4, 30), (1, 1, 0)];
        if let Some((_, _, degrees)) = table
            .iter()
            .find(|(n, d, _)| q == Rational::from_ratio(*n, *d).unwrap())
        {
            let degrees = if complement { 90 - degrees } else { *degrees };
            return Ok(Answer::exact(Enclosure::exact(Rational::from_i64(degrees))));
        }
    }
    let sin2 = Real::one() - cos2.clone();
    let mut enclosure = None;
    for bits in [128u32, 256] {
        let bits = bits.min(dec.policy().max_precision_bits);
        let c = dec.enclose_budgeted_at(&cos2, bits)?;
        let s = dec.enclose_budgeted_at(&sin2, bits)?;
        let clamp = |x: Ball| {
            if x.lo().sign() == Sign::Negative {
                Ball::from_rational(&Rational::zero(), bits).hull(&x)
            } else {
                x
            }
        };
        let (Some(c), Some(s)) = (clamp(c).sqrt(), clamp(s).sqrt()) else {
            return Err(Refusal::new(
                "measure.angle-domain",
                "the angle's sine and cosine could not be separated from zero",
            ));
        };
        let (y, x) = if complement { (c, s) } else { (s, c) };
        let Some(radians) = Ball::atan2(&y, &x) else {
            return Err(Refusal::new(
                "measure.angle-domain",
                "the angle's sine and cosine could not be separated from zero",
            ));
        };
        let pi = Ball::pi(bits);
        let degrees = (&radians * &Ball::from_rational(&Rational::from_i64(180), bits))
            .checked_div(&pi)
            .expect("pi is nonzero");
        let e = Enclosure::from_ball(&degrees);
        let reached = accuracy.reached(&e);
        enclosure = Some(e);
        if reached {
            break;
        }
    }
    let enclosure = enclosure.expect("at least one refinement");
    // Exact ends never leave [0, 90].
    let enclosure = Enclosure {
        lo: enclosure.lo.max(Rational::zero()),
        hi: enclosure.hi.min(Rational::from_i64(90)),
    };
    let mut marks = MarkLog::new();
    accuracy.mark(&mut marks, &[&enclosure], "angle");
    Ok(Answer::with_marks(enclosure, marks))
}
