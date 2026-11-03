//! Certified readouts between two faces or edges with plane, cylinder, cone,
//! sphere or torus carriers, or line or circle curves.

use crate::{Accuracy, Enclosure};
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{
    Answer, Answered, Decider, Integer, MarkLog, Outcome, Rational, Real, Refusal, Sign, StrictSign,
};
use catallax_shapes::{Point3, Surface, Vector3};
use catallax_topology::{EdgeId, FaceId, Model};

mod axis;
mod critical;
mod gap;

/// One end of an entity-pair readout.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum PairEntity {
    /// The whole trimmed face.
    Face(FaceId),
    /// The whole trimmed edge.
    Edge(EdgeId),
}

/// The scalar requested from two entities.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum PairRequest {
    /// Minimum distance over both closed trimmed entities, in model units.
    MinimumDistance,
    /// Distance between the entities' centres, in model units: the axis of a
    /// cylinder or cone, the centre of a sphere, torus or circle, the carrier
    /// line of a straight edge and the supporting plane of a planar face.
    CentreDistance,
    /// Unoriented acute angle in degrees between the entities' directions: the
    /// axis of a cylinder or cone, a straight edge's direction, and the normal
    /// of a planar face or of the plane of a circle or a torus's core circle.
    /// An axis and a normal meet at 90 degrees less the angle between them.
    AngleDegrees,
}

pub(crate) fn measure(
    model: &Model,
    a: PairEntity,
    b: PairEntity,
    request: PairRequest,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    accuracy.validate();
    dec.check_deadline()?;
    match request {
        PairRequest::MinimumDistance => gap::minimum(model, a, b, accuracy, dec),
        PairRequest::CentreDistance => axis::centre_distance(model, a, b, accuracy, dec),
        PairRequest::AngleDegrees => axis::angle(model, a, b, accuracy, dec),
    }
}

/// Exact sign, refusing when it cannot be decided.
fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}

fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Zero)
}

/// The nonnegative square root, kept rational for rational perfect squares.
fn root(value: &Real) -> Real {
    fn integer_root(n: &Integer) -> Option<Integer> {
        if n.sign() == Sign::Negative {
            return None;
        }
        if n.is_zero() {
            return Some(Integer::zero());
        }
        let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
        loop {
            let y = (&x + &n.div_rem_floor(&x)?.0)
                .div_rem_floor(&Integer::from_i64(2))?
                .0;
            if y >= x {
                return (&x * &x == *n).then_some(x);
            }
            x = y;
        }
    }
    value
        .rational_ref()
        .and_then(|q| Rational::new(integer_root(q.numerator())?, integer_root(q.denominator())?))
        .map_or_else(|| value.sqrt(), Real::from)
}

/// `v / |v|` for a nonzero vector.
fn unit(v: &Vector3) -> Vector3 {
    let length = root(&v.norm_squared());
    let inverse = match length.rational_ref() {
        Some(q) => Real::from(q.recip().expect("nonzero length")),
        None => Real::one() / length,
    };
    v.scaled(&inverse)
}

/// An enclosure of a nonnegative exact value: exact when the value reduces to a
/// rational, otherwise refined up to 256 bits and marked when still too wide.
fn enclose(
    value: &Real,
    accuracy: &Accuracy,
    property: &str,
    dec: &mut Decider,
) -> Answered<Enclosure> {
    let exact = match value.rational_ref() {
        Some(q) => Some(q.clone()),
        None => dec.rational_value(value)?,
    };
    let enclosure = match exact {
        Some(q) => Enclosure::exact(q),
        None => {
            let mut enclosure = None;
            for bits in [128, 256] {
                let bits = bits.min(dec.policy().max_precision_bits);
                let e = Enclosure::from_ball(&dec.enclose_budgeted_at(value, bits)?);
                let reached = accuracy.reached(&e);
                enclosure = Some(e);
                if reached {
                    break;
                }
            }
            enclosure.expect("at least one refinement")
        }
    };
    let mut marks = MarkLog::new();
    accuracy.mark(&mut marks, &[&enclosure], property);
    Ok(Answer::with_marks(enclosure, marks))
}

fn unsupported_surface(face: FaceId) -> Refusal {
    Refusal::new(
        "measure.unsupported-surface",
        format!("face {face:?} is not a plane, cylinder, cone, sphere or ring torus"),
    )
}

fn unsupported_curve(edge: EdgeId) -> Refusal {
    Refusal::new(
        "measure.unsupported-curve",
        format!("edge {edge:?} is not a finite straight edge or a circle"),
    )
}

/// The analytic carrier of a supported face.
fn face_surface(model: &Model, face: FaceId) -> Result<&Surface, Refusal> {
    let surface = model.geom().surface(model.face(face).surface);
    match surface {
        Surface::Plane(_)
        | Surface::Cylinder(_)
        | Surface::Cone(_)
        | Surface::Sphere(_)
        | Surface::Torus(_) => Ok(surface),
        _ => Err(unsupported_surface(face)),
    }
}

/// A finite straight edge's exact ends, or a circle edge's exact circle. Curves
/// stored as intersections count when their exact form is a line or a circle.
enum EdgeCurve {
    Segment(Point3, Point3),
    Circle {
        center: Point3,
        normal: Vector3,
        radius_sq: Real,
        ref_dir: Vector3,
    },
}

fn edge_curve(model: &Model, edge: EdgeId, dec: &mut Decider) -> Result<EdgeCurve, Refusal> {
    let e = model.edge(edge);
    match curve_form(model.geom(), e.curve, dec)? {
        CurveForm::Line { .. } => {
            let (Some(start), Some(end)) = (e.start, e.end) else {
                return Err(unsupported_curve(edge));
            };
            let at = |v| model.geom().position(model.vertex(v).point);
            Ok(EdgeCurve::Segment(at(start), at(end)))
        }
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        } => Ok(EdgeCurve::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        }),
        _ => Err(unsupported_curve(edge)),
    }
}

#[cfg(test)]
mod tests;
