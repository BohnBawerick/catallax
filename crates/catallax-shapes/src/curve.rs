//! Analytic curves and intersection curves stored by definition.

use catallax_numbers::{Ball, MPoly, Rational, Refusal};

use crate::math;
use crate::{GeomStore, Point3, PointId, Pose, RPoint3, RVector3, SurfaceId};

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
    /// Branch of the meeting of two surfaces.
    Intersection,
}

/// A line `origin + t dir`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Line {
    /// Point at `t = 0`.
    pub origin: RPoint3,
    /// Direction; non-zero; the natural direction.
    pub dir: RVector3,
}

/// A circle. Natural direction is counter-clockwise seen from the tip of `normal`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Circle {
    /// Centre.
    pub center: RPoint3,
    /// Normal of the circle's plane; non-zero.
    pub normal: RVector3,
    /// Squared radius; positive.
    pub radius_sq: Rational,
    /// Direction of `t = 0`; perpendicular to `normal`.
    pub ref_dir: RVector3,
}

/// An ellipse. Natural direction is counter-clockwise seen from the tip of `normal`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Ellipse {
    /// Centre.
    pub center: RPoint3,
    /// Normal of the ellipse's plane; non-zero.
    pub normal: RVector3,
    /// Direction of the major axis and of `t = 0`; perpendicular to `normal`.
    pub major_dir: RVector3,
    /// Squared major semi-axis.
    pub major_sq: Rational,
    /// Squared minor semi-axis; `0 < minor_sq < major_sq` (equal axes make a circle).
    pub minor_sq: Rational,
}

/// The connected component, of the meeting set of two surfaces with its
/// singular points removed, that contains `witness`. Natural direction is
/// `grad f0 x grad f1`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct IntersectionCurve {
    /// The two surfaces.
    pub surfaces: [SurfaceId; 2],
    /// A regular point of the branch.
    pub witness: PointId,
}

impl Line {
    /// A validated line.
    pub fn new(origin: RPoint3, dir: RVector3) -> Result<Self, Refusal> {
        if dir.is_zero() {
            return Err(Refusal::new("shapes.zero-vector", "line direction is zero"));
        }
        Ok(Self { origin, dir })
    }
}

impl Circle {
    /// A validated circle.
    pub fn new(
        center: RPoint3,
        normal: RVector3,
        radius_sq: Rational,
        ref_dir: RVector3,
    ) -> Result<Self, Refusal> {
        math::directions(&normal, &ref_dir)?;
        math::positive(&radius_sq, "shapes.non-positive-radius")?;
        math::perpendicular(&normal, &ref_dir)?;
        Ok(Self {
            center,
            normal,
            radius_sq,
            ref_dir,
        })
    }
}

impl Ellipse {
    /// A validated ellipse.
    pub fn new(
        center: RPoint3,
        normal: RVector3,
        major_dir: RVector3,
        major_sq: Rational,
        minor_sq: Rational,
    ) -> Result<Self, Refusal> {
        math::directions(&normal, &major_dir)?;
        if minor_sq <= Rational::zero() || minor_sq >= major_sq {
            return Err(Refusal::new(
                "shapes.ellipse-axes",
                "ellipse requires 0 < minor_sq < major_sq",
            ));
        }
        math::perpendicular(&normal, &major_dir)?;
        Ok(Self {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        })
    }
}

/// A curve.
#[non_exhaustive]
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Curve {
    /// Line.
    Line(Line),
    /// Circle.
    Circle(Circle),
    /// Ellipse.
    Ellipse(Ellipse),
    /// Intersection branch.
    Intersection(IntersectionCurve),
}

impl Curve {
    /// The kind.
    pub fn kind(&self) -> CurveKind {
        match self {
            Self::Line(_) => CurveKind::Line,
            Self::Circle(_) => CurveKind::Circle,
            Self::Ellipse(_) => CurveKind::Ellipse,
            Self::Intersection(_) => CurveKind::Intersection,
        }
    }

    /// Polynomials in `x, y, z` whose common zero set contains the curve, in
    /// the exact form and order of the piece description; for a line, a
    /// circle and an ellipse the common zero set is exactly the curve.
    pub fn equations(&self, store: &GeomStore) -> Vec<MPoly> {
        match self {
            Self::Line(c) => math::cross(&math::displacement(&c.origin), &c.dir)
                .into_iter()
                .filter(|p| !p.is_zero())
                .collect(),
            Self::Circle(c) => {
                let d = math::displacement(&c.center);
                vec![
                    math::dot(&d, &c.normal),
                    math::norm(&d) - MPoly::constant(3, c.radius_sq.clone()),
                ]
            }
            Self::Ellipse(c) => {
                let d = math::displacement(&c.center);
                let d2 = c.normal.cross(&c.major_dir);
                let a = math::dot(&d, &c.major_dir);
                let b = math::dot(&d, &d2);
                vec![
                    math::dot(&d, &c.normal),
                    (&a * &a).scale(&(&c.major_dir.norm_squared() * &c.major_sq).recip().unwrap())
                        + (&b * &b).scale(&(&d2.norm_squared() * &c.minor_sq).recip().unwrap())
                        - MPoly::constant(3, Rational::one()),
                ]
            }
            Self::Intersection(c) => c.surfaces.map(|s| store.implicit(s)).into(),
        }
    }

    /// Whether the curve is closed; `None` for intersection curves (piece 4 decides).
    pub fn is_closed(&self) -> Option<bool> {
        match self {
            Self::Line(_) => Some(false),
            Self::Circle(_) | Self::Ellipse(_) => Some(true),
            Self::Intersection(_) => None,
        }
    }

    /// The same curve moved by `pose`; `None` for intersection curves.
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

    /// Enclosure of the point at parameter `t`; `None` for intersection curves.
    pub fn evaluate(&self, t: &Ball) -> Option<[Ball; 3]> {
        let precision = t.precision();
        let (center, normal, reference, major, minor) = match self {
            Self::Line(c) => {
                let origin = math::point(&c.origin, precision);
                let direction = [&c.dir.x, &c.dir.y, &c.dir.z].map(|x| math::ball(x, precision));
                return Some(math::add(&origin, &math::scale(&direction, t)));
            }
            Self::Circle(c) => (&c.center, &c.normal, &c.ref_dir, &c.radius_sq, &c.radius_sq),
            Self::Ellipse(c) => (&c.center, &c.normal, &c.major_dir, &c.major_sq, &c.minor_sq),
            Self::Intersection(_) => return None,
        };
        let [e1, e2, _] = math::frame(normal, reference, precision);
        let a = math::ball(major, precision)
            .sqrt()
            .expect("positive squared radius")
            * t.cos();
        let b = math::ball(minor, precision)
            .sqrt()
            .expect("positive squared radius")
            * t.sin();
        Some(math::add(
            &math::point(center, precision),
            &math::add(&math::scale(&e1, &a), &math::scale(&e2, &b)),
        ))
    }

    /// The exact point at `t = 0`: a line's origin, the point of a circle or
    /// ellipse in its reference direction, an intersection curve's witness.
    pub fn parameter_origin(&self, store: &GeomStore) -> Point3 {
        match self {
            Self::Line(c) => c.origin.to_real(),
            Self::Circle(c) => c.center.to_real().plus(&c.ref_dir.to_real().scaled(
                &catallax_numbers::Real::from(&c.radius_sq / &c.ref_dir.norm_squared()).sqrt(),
            )),
            Self::Ellipse(c) => c.center.to_real().plus(&c.major_dir.to_real().scaled(
                &catallax_numbers::Real::from(&c.major_sq / &c.major_dir.norm_squared()).sqrt(),
            )),
            Self::Intersection(c) => store.position(c.witness),
        }
    }
}
