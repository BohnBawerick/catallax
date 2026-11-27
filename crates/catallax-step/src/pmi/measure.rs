//! Nominal values of dimensions, measured from the exact geometry.
//!
//! Predicates on the dimension's faces are exact. Angular values are
//! converted to degrees with a certified decimal rounding bound.

use std::collections::BTreeSet;

use catallax_numbers::{Ball, Decider, Outcome, Rational, Real, Refusal, StrictSign};
use catallax_shapes::{Point3, Surface, Vector3};
use catallax_topology::{EdgeId, FaceId, Loop, Model, Sense};

use super::{Anchor, At, Feature, HoleHead, Item, ItemKind, Pmi, Size, SizeKind};
use crate::number::{self, sqrt_real};

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

struct OrientedPlane {
    point: Point3,
    normal: Vector3,
}

fn oriented_plane(
    model: &Model,
    feature: &[FaceId],
    item: &str,
    dec: &mut Decider,
) -> Result<OrientedPlane, Refusal> {
    let plane = |id| {
        let face = model.face(id);
        let Surface::Plane(plane) = model.geom().surface(face.surface) else {
            return Err(refuse(item, "the dimension needs planar faces"));
        };
        let normal = if face.sense == Sense::Same {
            plane.normal.clone()
        } else {
            plane.normal.scaled(&(-Real::one()))
        };
        Ok(OrientedPlane {
            point: plane.origin.clone(),
            normal,
        })
    };
    let (first, rest) = feature
        .split_first()
        .ok_or_else(|| refuse(item, "a feature has no face"))?;
    let first = plane(*first)?;
    for id in rest {
        let other = plane(*id)?;
        if !is_zero_vector(&first.normal.cross(&other.normal), dec)?
            || !is_positive(&first.normal.dot(&other.normal), dec)?
            || !is_zero(&other.point.minus(&first.point).dot(&first.normal), dec)?
        {
            return Err(refuse(
                item,
                "the feature's faces are not one oriented plane",
            ));
        }
    }
    Ok(first)
}

fn face_edges(model: &Model, face: FaceId) -> impl Iterator<Item = EdgeId> + '_ {
    model
        .face(face)
        .loops
        .iter()
        .filter_map(|id| match model.loop_(*id) {
            Loop::Coedges(edges) => Some(edges),
            Loop::Vertex(_) => None,
        })
        .flatten()
        .map(|id| model.coedge(*id).edge)
}

/// The setback measured along `side`, facing away from the theoretical
/// corner. `across` is the other side's outward normal projected into `side`,
/// scaled by the squared norm of this side's normal.
fn setback(
    side: &OrientedPlane,
    other: &OrientedPlane,
    bevel: &OrientedPlane,
    item: &str,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let n = &side.normal;
    let m = &other.normal;
    let across = m.scaled(&n.norm_squared()).minus(&n.scaled(&m.dot(n)));
    let slope = across.dot(&bevel.normal);
    if !is_positive(&slope, dec)? {
        return Err(refuse(
            item,
            "the bevel must cut back the theoretical corner",
        ));
    }
    let at_corner = other.point.minus(&side.point).dot(m) / across.dot(m);
    let corner = side.point.plus(&across.scaled(&at_corner));
    let distance = corner.minus(&bevel.point).dot(&bevel.normal) / slope;
    if !is_positive(&distance, dec)? {
        return Err(refuse(
            item,
            "the bevel must cut back the theoretical corner",
        ));
    }
    Ok(distance * sqrt_real(&across.norm_squared()))
}

fn bevel_planes(
    model: &Model,
    bevel: &Feature<FaceId>,
    sides: &[Feature<FaceId>; 2],
    item: &str,
    dec: &mut Decider,
) -> Result<[OrientedPlane; 3], Refusal> {
    let a = oriented_plane(model, &sides[0], item, dec)?;
    let b = oriented_plane(model, &sides[1], item, dec)?;
    let h = oriented_plane(model, bevel, item, dec)?;
    let along = a.normal.cross(&b.normal);
    let edges: BTreeSet<_> = bevel.iter().flat_map(|id| face_edges(model, *id)).collect();
    let adjacent = |side: &Feature<FaceId>| {
        side.iter()
            .flat_map(|id| face_edges(model, *id))
            .any(|edge| edges.contains(&edge))
    };
    if is_zero_vector(&along, dec)?
        || !is_zero(&along.dot(&h.normal), dec)?
        || !sides.iter().all(adjacent)
    {
        return Err(refuse(
            item,
            "the bevel must meet both sides along their corner direction",
        ));
    }
    Ok([a, b, h])
}

fn angle_deg(rise: &Real, run: &Real, item: &str, dec: &mut Decider) -> Result<Rational, Refusal> {
    if is_zero(run, dec)? {
        return Ok(Rational::from_i64(90));
    }
    const DIGITS: u32 = 15;
    let precision = number::work_precision(dec.policy(), DIGITS * 4 + 48, 2)
        .ok_or_else(|| number::budget("no work available for the PMI angle"))?;
    let rise = dec.enclose_at(rise, precision)?;
    let run = dec.enclose_at(run, precision)?;
    let pi = Ball::pi_checked(precision, dec)?;
    let radians = crate::uv::angle(&rise, &run, &pi, Ball::atan_checked, dec)?;
    let scale = Ball::from_rational(&Rational::from_i64(180), precision);
    let degrees = (&radians * &scale)
        .checked_div(&pi)
        .expect("pi is positive");
    let number = number::certified(number::ball_number(&degrees, DIGITS), DIGITS)?;
    if number.value <= Rational::zero() || number.value >= Rational::from_i64(180) {
        return Err(refuse(
            item,
            "the rounded angle is not between 0 and 180 degrees",
        ));
    }
    Ok(number.value)
}

pub(super) struct Measured {
    pub(super) value: Real,
    pub(super) bevel_angle: Option<Rational>,
}

pub(super) fn dimension(
    model: &Model,
    item: &Item<FaceId>,
    dec: &mut Decider,
) -> Result<Option<Measured>, Refusal> {
    let id = item.id.as_str();
    let mut bevel_angle = None;
    let value = match &item.kind {
        ItemKind::Size(dimension) => match &dimension.kind {
            SizeKind::Hole(hole) => hole_size(model, dimension, hole.head.as_ref(), id, dec)?,
            kind => size(model, kind, &dimension.features, id, dec)?,
        },
        ItemKind::Location(location) => {
            let from = anchor_shape(model, &location.from, id, dec)?;
            let to = anchor_shape(model, &location.to, id, dec)?;
            let mut distance = separation(&from, &to, id, dec)?;
            for (shape, at) in [(&from, location.from.at), (&to, location.to.at)] {
                if let Some(radius) = side(shape, at) {
                    distance = &distance + &radius;
                }
            }
            if !is_positive(&distance, dec)? {
                return Err(refuse(id, "the two ends are not apart"));
            }
            distance
        }
        ItemKind::Width { flats: [a, b], .. } => {
            let a = oriented_plane(model, a, id, dec)?;
            let b = oriented_plane(model, b, id, dec)?;
            let gap = b.point.minus(&a.point).dot(&a.normal);
            if !is_zero_vector(&a.normal.cross(&b.normal), dec)?
                || !is_positive(&(-a.normal.dot(&b.normal)), dec)?
                || is_zero(&gap, dec)?
            {
                return Err(refuse(
                    id,
                    "width needs two distinct opposed parallel planes",
                ));
            }
            magnitude(&gap) / sqrt_real(&a.normal.norm_squared())
        }
        ItemKind::Chamfer { bevel, sides, .. } | ItemKind::EdgeBreak { bevel, sides, .. } => {
            let [a, b, h] = bevel_planes(model, bevel, sides, id, dec)?;
            let first = setback(&a, &b, &h, id, dec)?;
            let second = setback(&b, &a, &h, id, dec)?;
            if matches!(&item.kind, ItemKind::Chamfer { .. }) {
                let rise = sqrt_real(&a.normal.cross(&h.normal).norm_squared());
                bevel_angle = Some(angle_deg(
                    &rise,
                    &magnitude(&a.normal.dot(&h.normal)),
                    id,
                    dec,
                )?);
                first
            } else {
                match dec.sign(&(&second - &first)) {
                    Outcome::Proven(StrictSign::Positive) => second,
                    Outcome::Proven(StrictSign::Negative) | Outcome::ProvenCoincident(_) => first,
                    Outcome::Refused(refusal) => return Err(refusal),
                    Outcome::NotExact(inexact) => return Err(Refusal::from(inexact)),
                }
            }
        }
        ItemKind::BasicAngle { planes: [a, b], .. } => {
            let a = oriented_plane(model, a, id, dec)?;
            let b = oriented_plane(model, b, id, dec)?;
            let cross = a.normal.cross(&b.normal);
            if is_zero_vector(&cross, dec)? {
                return Err(refuse(id, "basic angle needs two nonparallel planes"));
            }
            angle_deg(
                &sqrt_real(&cross.norm_squared()),
                &(-a.normal.dot(&b.normal)),
                id,
                dec,
            )?
            .into()
        }
        ItemKind::Datum(_) | ItemKind::Tolerance(_) | ItemKind::Texture(_) => return Ok(None),
    };
    Ok(Some(Measured { value, bevel_angle }))
}

/// Each dimension's primary nominal on `model`, in millimetres (degrees
/// for a basic angle). For a chamfer this is the reference-side setback.
/// There is one entry per item, `None` for non-dimensions, and imported
/// nominals are ignored. Refuses a dimension that does not fit its faces.
pub fn nominals(
    model: &Model,
    pmi: &Pmi<FaceId>,
    dec: &mut Decider,
) -> Result<Vec<Option<Real>>, Refusal> {
    pmi.items
        .iter()
        .map(|item| Ok(dimension(model, item, dec)?.map(|measured| measured.value)))
        .collect()
}
