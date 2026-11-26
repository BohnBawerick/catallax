//! Nominal values of dimensions, measured from the exact geometry.
//!
//! A size is read off a cylinder, a distance between planes and cylinder
//! axes. Every geometric condition a dimension rests on is decided exactly:
//! the faces of a round feature share one axis and radius, the two ends of a
//! distance are parallel. A condition that does not hold refuses
//! `step.pmi-not-measurable` with the item's name; nothing is estimated.

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, StrictSign};
use catallax_shapes::{Point3, Surface, Vector3};
use catallax_topology::{FaceId, Model};

use super::{Anchor, At, Feature, HoleHead, ItemKind, Pmi, Size, SizeKind};
use crate::number::sqrt_real;

/// What a feature is for measuring.
enum Shape {
    Plane {
        point: Point3,
        normal: Vector3,
    },
    /// A cylinder, or several faces of one.
    Axis {
        point: Point3,
        direction: Vector3,
        radius_sq: Real,
    },
}

fn refuse(item: &str, what: impl std::fmt::Display) -> Refusal {
    Refusal::new("step.pmi-not-measurable", format!("{item}: {what}"))
}

fn is_zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::Refused(refusal) => Err(refusal),
        Outcome::NotExact(inexact) => Err(Refusal::from(inexact)),
    }
}

fn is_positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(sign) => Ok(sign == StrictSign::Positive),
        Outcome::ProvenCoincident(_) => Ok(false),
        Outcome::Refused(refusal) => Err(refusal),
        Outcome::NotExact(inexact) => Err(Refusal::from(inexact)),
    }
}

fn is_zero_vector(value: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(is_zero(&value.x, dec)? && is_zero(&value.y, dec)? && is_zero(&value.z, dec)?)
}

fn magnitude(value: &Real) -> Real {
    sqrt_real(&value.square())
}

fn face_shape(model: &Model, face: FaceId, item: &str) -> Result<Shape, Refusal> {
    match model.geom().surface(model.face(face).surface) {
        Surface::Plane(plane) => Ok(Shape::Plane {
            point: plane.origin.clone(),
            normal: plane.normal.clone(),
        }),
        Surface::Cylinder(cylinder) => Ok(Shape::Axis {
            point: cylinder.origin.clone(),
            direction: cylinder.axis.clone(),
            radius_sq: cylinder.radius_sq.clone(),
        }),
        _ => Err(refuse(item, "only planes and cylinders can be dimensioned")),
    }
}

/// True when two shapes are the same plane, or cylinders of one axis and radius.
fn coincident(a: &Shape, b: &Shape, dec: &mut Decider) -> Result<bool, Refusal> {
    match (a, b) {
        (
            Shape::Plane { point, normal },
            Shape::Plane {
                point: other_point,
                normal: other_normal,
            },
        ) => Ok(is_zero_vector(&normal.cross(other_normal), dec)?
            && is_zero(&other_point.minus(point).dot(normal), dec)?),
        (
            Shape::Axis {
                point,
                direction,
                radius_sq,
            },
            Shape::Axis {
                point: other_point,
                direction: other_direction,
                radius_sq: other_radius_sq,
            },
        ) => Ok(is_zero_vector(&direction.cross(other_direction), dec)?
            && is_zero_vector(&other_point.minus(point).cross(direction), dec)?
            && is_zero(&(radius_sq - other_radius_sq), dec)?),
        _ => Ok(false),
    }
}

/// The one shape all faces of a feature lie on.
fn feature_shape(
    model: &Model,
    feature: &[FaceId],
    item: &str,
    dec: &mut Decider,
) -> Result<Shape, Refusal> {
    let mut faces = feature.iter();
    let first = faces
        .next()
        .ok_or_else(|| refuse(item, "a feature has no face"))?;
    let shape = face_shape(model, *first, item)?;
    for face in faces {
        if !coincident(&shape, &face_shape(model, *face, item)?, dec)? {
            return Err(refuse(
                item,
                "the feature's faces are not one plane or one cylinder",
            ));
        }
    }
    Ok(shape)
}

fn midpoint(a: &Point3, b: &Point3) -> Point3 {
    let half = Real::from_rational(Rational::from_ratio(1, 2).expect("a nonzero denominator"));
    a.plus(&b.minus(a).scaled(&half))
}

/// The shape a distance ends on.
fn anchor_shape(
    model: &Model,
    anchor: &Anchor<FaceId>,
    item: &str,
    dec: &mut Decider,
) -> Result<Shape, Refusal> {
    if anchor.at != At::Middle {
        let shape = feature_shape(model, &anchor.feature, item, dec)?;
        if matches!(anchor.at, At::Inner | At::Outer) && matches!(shape, Shape::Plane { .. }) {
            return Err(refuse(
                item,
                "only a round feature has an inner and an outer side",
            ));
        }
        return Ok(shape);
    }
    let [first, second] = anchor.feature.as_slice() else {
        return Err(refuse(item, "the middle is taken between two faces"));
    };
    match (
        face_shape(model, *first, item)?,
        face_shape(model, *second, item)?,
    ) {
        (
            Shape::Plane { point, normal },
            Shape::Plane {
                point: other_point,
                normal: other_normal,
            },
        ) => {
            if !is_zero_vector(&normal.cross(&other_normal), dec)? {
                return Err(refuse(item, "the two faces are not parallel"));
            }
            Ok(Shape::Plane {
                point: midpoint(&point, &other_point),
                normal,
            })
        }
        (
            Shape::Axis {
                point,
                direction,
                radius_sq,
            },
            Shape::Axis {
                point: other_point,
                direction: other_direction,
                radius_sq: other_radius_sq,
            },
        ) => {
            if !is_zero_vector(&direction.cross(&other_direction), dec)? {
                return Err(refuse(item, "the two axes are not parallel"));
            }
            if !is_zero(&(&radius_sq - &other_radius_sq), dec)? {
                return Err(refuse(item, "the two round faces differ in radius"));
            }
            Ok(Shape::Axis {
                point: midpoint(&point, &other_point),
                direction,
                radius_sq,
            })
        }
        _ => Err(refuse(item, "the two faces are not of one kind")),
    }
}

/// The distance between the planes and axes two shapes stand for.
fn separation(a: &Shape, b: &Shape, item: &str, dec: &mut Decider) -> Result<Real, Refusal> {
    match (a, b) {
        (
            Shape::Plane { point, normal },
            Shape::Plane {
                point: other_point,
                normal: other_normal,
            },
        ) => {
            if !is_zero_vector(&normal.cross(other_normal), dec)? {
                return Err(refuse(item, "the two planes are not parallel"));
            }
            Ok(
                magnitude(&other_point.minus(point).dot(normal))
                    / sqrt_real(&normal.norm_squared()),
            )
        }
        (
            Shape::Plane { point, normal },
            Shape::Axis {
                point: axis_point,
                direction,
                ..
            },
        )
        | (
            Shape::Axis {
                point: axis_point,
                direction,
                ..
            },
            Shape::Plane { point, normal },
        ) => {
            if !is_zero(&direction.dot(normal), dec)? {
                return Err(refuse(item, "the axis is not parallel to the plane"));
            }
            Ok(magnitude(&axis_point.minus(point).dot(normal)) / sqrt_real(&normal.norm_squared()))
        }
        (
            Shape::Axis {
                point, direction, ..
            },
            Shape::Axis {
                point: other_point,
                direction: other_direction,
                ..
            },
        ) => {
            if !is_zero_vector(&direction.cross(other_direction), dec)? {
                return Err(refuse(item, "the two axes are not parallel"));
            }
            let across = other_point.minus(point).cross(direction);
            Ok(sqrt_real(&across.norm_squared()) / sqrt_real(&direction.norm_squared()))
        }
    }
}

/// What a round end adds to a distance between axes: its radius outwards,
/// minus its radius inwards.
fn side(shape: &Shape, at: At) -> Option<Real> {
    let Shape::Axis { radius_sq, .. } = shape else {
        return None;
    };
    match at {
        At::Outer => Some(sqrt_real(radius_sq)),
        At::Inner => Some(Real::zero() - sqrt_real(radius_sq)),
        At::Feature | At::Middle => None,
    }
}

fn size(
    model: &Model,
    kind: &SizeKind,
    features: &[Feature<FaceId>],
    item: &str,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let mut common: Option<Real> = None;
    for feature in features {
        let Shape::Axis { radius_sq, .. } = feature_shape(model, feature, item, dec)? else {
            return Err(refuse(item, "a diameter or radius belongs to a cylinder"));
        };
        match &common {
            None => common = Some(radius_sq),
            Some(first) => {
                if !is_zero(&(first - &radius_sq), dec)? {
                    return Err(refuse(item, "the pattern's features differ in size"));
                }
            }
        }
    }
    let radius = sqrt_real(&common.ok_or_else(|| refuse(item, "the size names no feature"))?);
    Ok(match kind {
        SizeKind::Diameter | SizeKind::Hole(_) => &radius + &radius,
        SizeKind::Radius => radius,
    })
}

fn hole_size(
    model: &Model,
    size: &Size<FaceId>,
    head: Option<&HoleHead>,
    item: &str,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let [pilot, heads @ ..] = size.features.as_slice() else {
        return Err(refuse(item, "the hole names no pilot"));
    };
    let Shape::Axis {
        point,
        direction,
        radius_sq,
    } = feature_shape(model, pilot, item, dec)?
    else {
        return Err(refuse(item, "a hole's pilot belongs to a cylinder"));
    };
    let radius = sqrt_real(&radius_sq);
    let diameter = &radius + &radius;
    let (head, faces) = match (head, heads) {
        (None, []) => return Ok(diameter),
        (Some(head), [faces]) if !faces.is_empty() => (head, faces),
        _ => return Err(refuse(item, "the hole's head has no unique face set")),
    };
    let nominal = match head {
        HoleHead::Counterbore { diameter, .. } | HoleHead::Countersink { diameter, .. } => {
            Real::from_rational(diameter.clone())
        }
    };
    match dec.sign(&(&nominal - &diameter)) {
        Outcome::Proven(StrictSign::Positive) => {}
        Outcome::Refused(refusal) => return Err(refusal),
        Outcome::NotExact(inexact) => return Err(Refusal::from(inexact)),
        _ => {
            return Err(refuse(
                item,
                "the head diameter must be larger than the pilot",
            ));
        }
    }
    if matches!(head, HoleHead::Counterbore { .. }) {
        let Shape::Axis {
            point: origin,
            direction: axis,
            radius_sq: head_radius,
        } = feature_shape(model, faces, item, dec)?
        else {
            return Err(refuse(item, "the counterbore head is not cylindrical"));
        };
        if !is_zero_vector(&axis.cross(&direction), dec)?
            || !is_zero_vector(&origin.minus(&point).cross(&direction), dec)?
        {
            return Err(refuse(
                item,
                "the counterbore head does not share the pilot axis",
            ));
        }
        match dec.sign(&(&head_radius - &radius_sq)) {
            Outcome::Proven(StrictSign::Positive) => {}
            Outcome::Refused(refusal) => return Err(refusal),
            Outcome::NotExact(inexact) => return Err(Refusal::from(inexact)),
            _ => {
                return Err(refuse(
                    item,
                    "the counterbore head is not wider than the pilot",
                ));
            }
        }
    } else {
        let mut first = None;
        for face in faces {
            let Surface::Cone(cone) = model.geom().surface(model.face(*face).surface) else {
                return Err(refuse(item, "the countersink head is not conical"));
            };
            if !is_zero_vector(&cone.axis.cross(&direction), dec)?
                || !is_zero_vector(&cone.apex.minus(&point).cross(&direction), dec)?
            {
                return Err(refuse(
                    item,
                    "the countersink head does not share the pilot axis",
                ));
            }
            if let Some(other) = first {
                let other: &catallax_shapes::Cone = other;
                if !is_zero_vector(&cone.apex.minus(&other.apex), dec)?
                    || !is_positive(&cone.axis.dot(&other.axis), dec)?
                    || !is_zero(&(&cone.slope - &other.slope), dec)?
                {
                    return Err(refuse(item, "the countersink head faces are not one cone"));
                }
            } else {
                first = Some(cone);
            }
        }
    }
    Ok(diameter)
}

/// The nominal value of every dimension of `pmi` in millimetres, measured
/// on `model`. The result has one entry per item, `None` for an item that is
/// not a dimension. Any nominal an item states is ignored.
///
/// Refuses `step.pmi-not-measurable` when a dimension does not fit its
/// faces, such as a diameter on a plane or a distance between faces that are
/// not parallel.
pub fn nominals(
    model: &Model,
    pmi: &Pmi<FaceId>,
    dec: &mut Decider,
) -> Result<Vec<Option<Real>>, Refusal> {
    let mut values = Vec::with_capacity(pmi.items.len());
    for item in &pmi.items {
        let id = item.id.as_str();
        values.push(match &item.kind {
            ItemKind::Size(dimension) => Some(match &dimension.kind {
                SizeKind::Hole(hole) => hole_size(model, dimension, hole.head.as_ref(), id, dec)?,
                kind => size(model, kind, &dimension.features, id, dec)?,
            }),
            ItemKind::Location(location) => {
                let from = anchor_shape(model, &location.from, id, dec)?;
                let to = anchor_shape(model, &location.to, id, dec)?;
                let mut distance = separation(&from, &to, id, dec)?;
                for (shape, at) in [(&from, location.from.at), (&to, location.to.at)] {
                    if let Some(radius) = side(shape, at) {
                        distance = &distance + &radius;
                    }
                }
                match dec.sign(&distance) {
                    Outcome::Proven(StrictSign::Positive) => {}
                    Outcome::Refused(refusal) => return Err(refusal),
                    Outcome::NotExact(inexact) => return Err(Refusal::from(inexact)),
                    Outcome::Proven(StrictSign::Negative) | Outcome::ProvenCoincident(_) => {
                        return Err(refuse(id, "the two ends are not apart"));
                    }
                }
                Some(distance)
            }
            ItemKind::Datum(_) | ItemKind::Tolerance(_) | ItemKind::Texture(_) => None,
        });
    }
    Ok(values)
}
