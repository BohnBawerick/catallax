//! Analytic curves and certified intersection components stored by definition.

use crate::{GeomStore, Point3, PointId, Pose, SurfaceId, Vector3, math};
use catallax_numbers::{Ball, Decider, Real, RealPoly, Refusal, Sign};

/// The kind of a curve.
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum CurveKind {
    /// Straight line.
    Line,
    /// Circle.
    Circle,
    /// Ellipse.
    Ellipse,
    /// Intersection component.
    Intersection,
}

/// A line `origin + t dir`.
#[derive(Clone, Debug)]
pub struct Line {
    /// Parameter origin.
    pub origin: Point3,
    /// Nonzero natural direction.
    pub dir: Vector3,
}

/// A circle, counter-clockwise when viewed from its normal tip.
#[derive(Clone, Debug)]
pub struct Circle {
    /// Centre.
    pub center: Point3,
    /// Nonzero plane normal.
    pub normal: Vector3,
    /// Positive squared radius.
    pub radius_sq: Real,
    /// Nonzero perpendicular seam direction.
    pub ref_dir: Vector3,
}

/// An ellipse, counter-clockwise when viewed from its normal tip.
#[derive(Clone, Debug)]
pub struct Ellipse {
    /// Centre.
    pub center: Point3,
    /// Nonzero plane normal.
    pub normal: Vector3,
    /// Nonzero perpendicular major-axis and seam direction.
    pub major_dir: Vector3,
    /// Squared major semi-axis.
    pub major_sq: Real,
    /// Squared minor semi-axis, with `0 < minor_sq < major_sq`.
    pub minor_sq: Real,
}

/// A connected regular component of the named surfaces through its witness.
/// Natural direction is `grad f0 x grad f1`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct IntersectionCurve {
    /// Defining surface pair.
    pub surfaces: [SurfaceId; 2],
    /// Regular witness point.
    pub witness: PointId,
}

impl Line {
    /// Validates all expression domains before proving the direction nonzero.
    pub fn new(origin: Point3, dir: Vector3, dec: &mut Decider) -> Result<Self, Refusal> {
        math::validate_point(&origin, dec)?;
        math::validate_vector(&dir, dec)?;
        if math::sign(&dir.norm_squared(), dec)? == Sign::Zero {
            return Err(Refusal::new("shapes.zero-vector", "line direction is zero"));
        }
        Ok(Self { origin, dir })
    }
}
impl Circle {
    /// Validates expression domains, directions, radius and perpendicularity.
    pub fn new(
        center: Point3,
        normal: Vector3,
        radius_sq: Real,
        ref_dir: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&center, dec)?;
        math::validate_vector(&normal, dec)?;
        math::validate(&radius_sq, dec)?;
        math::validate_vector(&ref_dir, dec)?;
        math::directions(&normal, &ref_dir, dec)?;
        math::positive(&radius_sq, "shapes.non-positive-radius", dec)?;
        math::perpendicular(&normal, &ref_dir, dec)?;
        Ok(Self {
            center,
            normal,
            radius_sq,
            ref_dir,
        })
    }
}
impl Ellipse {
    /// Validates expression domains before proving the strict ellipse axes.
    pub fn new(
        center: Point3,
        normal: Vector3,
        major_dir: Vector3,
        major_sq: Real,
        minor_sq: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&center, dec)?;
        math::validate_vector(&normal, dec)?;
        math::validate_vector(&major_dir, dec)?;
        math::validate(&major_sq, dec)?;
        math::validate(&minor_sq, dec)?;
        math::directions(&normal, &major_dir, dec)?;
        if math::sign(&minor_sq, dec)? != Sign::Positive
            || math::sign(&(&major_sq - &minor_sq), dec)? != Sign::Positive
        {
            return Err(Refusal::new(
                "shapes.ellipse-axes",
                "ellipse requires 0 < minor_sq < major_sq",
            ));
        }
        math::perpendicular(&normal, &major_dir, dec)?;
        Ok(Self {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        })
    }
}

/// An exact analytic curve or a certified intersection component.
#[non_exhaustive]
#[derive(Clone, Debug)]
pub enum Curve {
    /// Line.
    Line(Line),
    /// Circle.
    Circle(Circle),
    /// Ellipse.
    Ellipse(Ellipse),
    /// Intersection component.
    Intersection(IntersectionCurve),
}

impl Curve {
    /// Curve kind.
    pub fn kind(&self) -> CurveKind {
        match self {
            Self::Line(_) => CurveKind::Line,
            Self::Circle(_) => CurveKind::Circle,
            Self::Ellipse(_) => CurveKind::Ellipse,
            Self::Intersection(_) => CurveKind::Intersection,
        }
    }

    fn definition_fields(&self) -> Option<[Option<&Real>; 11]> {
        let (origin, axis, reference, major, minor) = match self {
            Self::Line(c) => (&c.origin, &c.dir, None, None, None),
            Self::Circle(c) => (
                &c.center,
                &c.normal,
                Some(&c.ref_dir),
                Some(&c.radius_sq),
                None,
            ),
            Self::Ellipse(c) => (
                &c.center,
                &c.normal,
                Some(&c.major_dir),
                Some(&c.major_sq),
                Some(&c.minor_sq),
            ),
            Self::Intersection(_) => return None,
        };
        Some([
            Some(&origin.x),
            Some(&origin.y),
            Some(&origin.z),
            Some(&axis.x),
            Some(&axis.y),
            Some(&axis.z),
            reference.map(|v| &v.x),
            reference.map(|v| &v.y),
            reference.map(|v| &v.z),
            major,
            minor,
        ])
    }

    /// Whether all analytic defining scalars have cheap rational values.
    /// Intersection components require their carrier and selector context.
    pub fn has_rational_definition(&self) -> bool {
        self.definition_fields().is_some_and(|fields| {
            fields
                .into_iter()
                .flatten()
                .all(|value| value.rational_ref().is_some())
        })
    }

    /// Safe definition identity, not a test for different geometric point sets.
    /// A false result requires a geometric predicate before concluding apartness.
    pub fn same_definition(&self, other: &Self) -> bool {
        if self.kind() != other.kind() {
            return false;
        }
        match (self.definition_fields(), other.definition_fields()) {
            (Some(a), Some(b)) => a.into_iter().zip(b).all(|(a, b)| match (a, b) {
                (Some(a), Some(b)) => a.same_definition(b),
                (None, None) => true,
                _ => false,
            }),
            (None, None) => {
                matches!((self, other), (Self::Intersection(a), Self::Intersection(b)) if a == b)
            }
            _ => false,
        }
    }

    /// Defining equations in the documented order. The equations contain an
    /// intersection component's carrier; its witness still selects the component.
    pub fn equations(&self, store: &GeomStore) -> Vec<RealPoly> {
        match self {
            Self::Line(c) => math::cross(&math::displacement(&c.origin), &c.dir)
                .into_iter()
                .filter(|p| !p.rational().is_some_and(|r| r.is_zero()))
                .collect(),
            Self::Circle(c) => {
                let d = math::displacement(&c.center);
                vec![
                    math::dot(&d, &c.normal),
                    math::norm(&d) - RealPoly::constant(3, c.radius_sq.clone()),
                ]
            }
            Self::Ellipse(c) => {
                let d = math::displacement(&c.center);
                let second = c.normal.cross(&c.major_dir);
                let a = math::dot(&d, &c.major_dir);
                let b = math::dot(&d, &second);
                vec![
                    math::dot(&d, &c.normal),
                    (&a * &a).scale(&(&Real::one() / &(&c.major_dir.norm_squared() * &c.major_sq)))
                        + (&b * &b)
                            .scale(&(&Real::one() / &(&second.norm_squared() * &c.minor_sq)))
                        - RealPoly::constant(3, Real::one()),
                ]
            }
            Self::Intersection(c) => c.surfaces.map(|s| store.implicit(s)).into(),
        }
    }

    /// Cheap closure query; component closure requires an intersection proof.
    pub fn is_closed(&self) -> Option<bool> {
        match self {
            Self::Line(_) => Some(false),
            Self::Circle(_) | Self::Ellipse(_) => Some(true),
            Self::Intersection(_) => None,
        }
    }

    /// Applies a validated rigid pose; intersection definitions retain their objects.
    pub fn transformed(&self, pose: &Pose) -> Option<Curve> {
        Some(match self {
            Self::Line(c) => Self::Line(Line {
                origin: pose.apply_point(&c.origin),
                dir: pose.apply_vector(&c.dir),
            }),
            Self::Circle(c) => Self::Circle(Circle {
                center: pose.apply_point(&c.center),
                normal: pose.apply_vector(&c.normal),
                radius_sq: c.radius_sq.clone(),
                ref_dir: pose.apply_vector(&c.ref_dir),
            }),
            Self::Ellipse(c) => Self::Ellipse(Ellipse {
                center: pose.apply_point(&c.center),
                normal: pose.apply_vector(&c.normal),
                major_dir: pose.apply_vector(&c.major_dir),
                major_sq: c.major_sq.clone(),
                minor_sq: c.minor_sq.clone(),
            }),
            Self::Intersection(_) => return None,
        })
    }

    /// Certified chart enclosure; `None` for intersection components.
    pub fn evaluate(&self, t: &Ball, dec: &Decider) -> Result<Option<[Ball; 3]>, Refusal> {
        let precision = t.precision();
        let (center, normal, reference, major, minor) = match self {
            Self::Line(c) => {
                let origin = math::point(&c.origin, precision, dec)?;
                let direction = [
                    dec.enclose_at(&c.dir.x, precision)?,
                    dec.enclose_at(&c.dir.y, precision)?,
                    dec.enclose_at(&c.dir.z, precision)?,
                ];
                return Ok(Some(math::add(&origin, &math::scale(&direction, t))));
            }
            Self::Circle(c) => (&c.center, &c.normal, &c.ref_dir, &c.radius_sq, &c.radius_sq),
            Self::Ellipse(c) => (&c.center, &c.normal, &c.major_dir, &c.major_sq, &c.minor_sq),
            Self::Intersection(_) => return Ok(None),
        };
        let [e1, e2, _] = math::frame(normal, reference, precision, dec)?;
        // Keep the rational chart's rounding order: enclose squared dimension,
        // then take its square root at the input precision.
        let a = dec
            .enclose_at(major, precision)?
            .sqrt()
            .expect("positive squared radius")
            * t.cos();
        let b = dec
            .enclose_at(minor, precision)?
            .sqrt()
            .expect("positive squared radius")
            * t.sin();
        Ok(Some(math::add(
            &math::point(center, precision, dec)?,
            &math::add(&math::scale(&e1, &a), &math::scale(&e2, &b)),
        )))
    }

    /// Exact parameter-origin point, retaining expressions and meeting definitions.
    pub fn parameter_origin(&self, store: &GeomStore) -> Point3 {
        match self {
            Self::Line(c) => c.origin.clone(),
            Self::Circle(c) => c.center.plus(
                &c.ref_dir
                    .scaled(&(&c.radius_sq / &c.ref_dir.norm_squared()).sqrt()),
            ),
            Self::Ellipse(c) => c.center.plus(
                &c.major_dir
                    .scaled(&(&c.major_sq / &c.major_dir.norm_squared()).sqrt()),
            ),
            Self::Intersection(c) => store.position(c.witness),
        }
    }
}
