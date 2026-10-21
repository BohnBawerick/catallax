//! Analytic surfaces with exact rational parameters.

use catallax_numbers::{Ball, MPoly, Rational, Real, Refusal, Sign, StrictSign};

use crate::math;
use crate::{Point3, Pose, RPoint3, RVector3, Vector3};

/// The kind of a surface.
#[non_exhaustive]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum SurfaceKind {
    /// Plane.
    Plane,
    /// Circular cylinder.
    Cylinder,
    /// Circular cone.
    Cone,
    /// Sphere.
    Sphere,
    /// Torus.
    Torus,
}

/// What values a surface parameter takes.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum ParamKind {
    /// All reals.
    Line,
    /// Period `2 pi`, seam at 0.
    Periodic,
    /// `[-pi/2, pi/2]`.
    Latitude,
    /// `[0, infinity)`.
    Ray,
}

/// A plane. Normal is `normal`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Plane {
    /// A point of the plane, the parameter origin.
    pub origin: RPoint3,
    /// Normal direction; non-zero.
    pub normal: RVector3,
    /// In-plane direction of the `u` parameter; perpendicular to `normal`.
    pub ref_dir: RVector3,
}

/// A circular cylinder. Normal points away from the axis.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Cylinder {
    /// A point of the axis, where `v = 0`.
    pub origin: RPoint3,
    /// Axis direction; non-zero.
    pub axis: RVector3,
    /// Radius; positive.
    pub radius: Rational,
    /// Seam direction; perpendicular to `axis`.
    pub ref_dir: RVector3,
}

/// A circular cone, the nappe on the side `axis` points to. Normal points away from the axis.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Cone {
    /// Apex.
    pub apex: RPoint3,
    /// Axis direction from the apex into the nappe; non-zero.
    pub axis: RVector3,
    /// Radius gained per unit distance along the axis; positive.
    pub slope: Rational,
    /// Seam direction; perpendicular to `axis`.
    pub ref_dir: RVector3,
}

/// A sphere. Normal points away from the centre.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Sphere {
    /// Centre.
    pub center: RPoint3,
    /// Radius; positive.
    pub radius: Rational,
    /// Pole direction; non-zero.
    pub axis: RVector3,
    /// Seam direction; perpendicular to `axis`.
    pub ref_dir: RVector3,
}

/// A torus with `0 < minor < major`. Normal points away from the tube's core circle.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Torus {
    /// Centre.
    pub center: RPoint3,
    /// Axis direction; non-zero.
    pub axis: RVector3,
    /// Radius of the core circle.
    pub major: Rational,
    /// Radius of the tube.
    pub minor: Rational,
    /// Seam direction; perpendicular to `axis`.
    pub ref_dir: RVector3,
}

impl Plane {
    /// A validated plane.
    pub fn new(origin: RPoint3, normal: RVector3, ref_dir: RVector3) -> Result<Self, Refusal> {
        math::directions(&normal, &ref_dir)?;
        math::perpendicular(&normal, &ref_dir)?;
        Ok(Self {
            origin,
            normal,
            ref_dir,
        })
    }
}

impl Cylinder {
    /// A validated cylinder.
    pub fn new(
        origin: RPoint3,
        axis: RVector3,
        radius: Rational,
        ref_dir: RVector3,
    ) -> Result<Self, Refusal> {
        math::directions(&axis, &ref_dir)?;
        math::positive(&radius, "shapes.non-positive-radius")?;
        math::perpendicular(&axis, &ref_dir)?;
        Ok(Self {
            origin,
            axis,
            radius,
            ref_dir,
        })
    }
}

impl Cone {
    /// A validated cone.
    pub fn new(
        apex: RPoint3,
        axis: RVector3,
        slope: Rational,
        ref_dir: RVector3,
    ) -> Result<Self, Refusal> {
        math::directions(&axis, &ref_dir)?;
        math::positive(&slope, "shapes.non-positive-slope")?;
        math::perpendicular(&axis, &ref_dir)?;
        Ok(Self {
            apex,
            axis,
            slope,
            ref_dir,
        })
    }
}

impl Sphere {
    /// A validated sphere.
    pub fn new(
        center: RPoint3,
        radius: Rational,
        axis: RVector3,
        ref_dir: RVector3,
    ) -> Result<Self, Refusal> {
        math::directions(&axis, &ref_dir)?;
        math::positive(&radius, "shapes.non-positive-radius")?;
        math::perpendicular(&axis, &ref_dir)?;
        Ok(Self {
            center,
            radius,
            axis,
            ref_dir,
        })
    }
}

impl Torus {
    /// A validated torus.
    pub fn new(
        center: RPoint3,
        axis: RVector3,
        major: Rational,
        minor: Rational,
        ref_dir: RVector3,
    ) -> Result<Self, Refusal> {
        math::directions(&axis, &ref_dir)?;
        math::positive(&major, "shapes.non-positive-radius")?;
        math::positive(&minor, "shapes.non-positive-radius")?;
        if minor >= major {
            return Err(Refusal::new(
                "shapes.torus-self-intersecting",
                "minor radius must be less than major radius",
            ));
        }
        math::perpendicular(&axis, &ref_dir)?;
        Ok(Self {
            center,
            axis,
            major,
            minor,
            ref_dir,
        })
    }
}

/// An analytic surface.
#[non_exhaustive]
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Surface {
    /// Plane.
    Plane(Plane),
    /// Cylinder.
    Cylinder(Cylinder),
    /// Cone.
    Cone(Cone),
    /// Sphere.
    Sphere(Sphere),
    /// Torus.
    Torus(Torus),
}

impl Surface {
    /// The kind.
    pub fn kind(&self) -> SurfaceKind {
        match self {
            Self::Plane(_) => SurfaceKind::Plane,
            Self::Cylinder(_) => SurfaceKind::Cylinder,
            Self::Cone(_) => SurfaceKind::Cone,
            Self::Sphere(_) => SurfaceKind::Sphere,
            Self::Torus(_) => SurfaceKind::Torus,
        }
    }

    /// Implicit polynomial in `x, y, z`, positive on the normal side, exactly
    /// as the piece description's table gives it, with no scaling.
    pub fn implicit(&self) -> MPoly {
        let (origin, axis, _) = self.definition_frame();
        let d = math::displacement(origin);
        let constant = |r| MPoly::constant(3, r);
        match self {
            Self::Plane(_) => math::dot(&d, axis),
            Self::Cylinder(s) => {
                math::norm(&math::cross(&d, axis))
                    - constant(&s.radius * &s.radius * axis.norm_squared())
            }
            Self::Cone(s) => {
                let h = math::dot(&d, axis);
                math::norm(&math::cross(&d, axis)) - (&h * &h).scale(&(&s.slope * &s.slope))
            }
            Self::Sphere(s) => math::norm(&d) - constant(&s.radius * &s.radius),
            Self::Torus(s) => {
                let q = math::norm(&d);
                let h = math::dot(&d, axis);
                let radial = &q - &(&h * &h).scale(&axis.norm_squared().recip().unwrap());
                let a = q + constant(&s.major * &s.major - &s.minor * &s.minor);
                &a * &a - radial.scale(&(Rational::from_i64(4) * (&s.major * &s.major)))
            }
        }
    }

    /// Gradient of the implicit polynomial at `p`.
    pub fn gradient(&self, p: &Point3) -> Vector3 {
        let (origin, axis, _) = self.definition_frame();
        let d = p.minus(&origin.to_real());
        let a = axis.to_real();
        let n = Real::from(axis.norm_squared());
        let h = d.dot(&a);
        let two = Real::from_i64(2);
        match self {
            Self::Plane(_) => a,
            Self::Cylinder(_) => d.scaled(&n).minus(&a.scaled(&h)).scaled(&two),
            Self::Cone(s) => {
                let k = Real::from(Rational::one() + &s.slope * &s.slope);
                d.scaled(&n).minus(&a.scaled(&(h * k))).scaled(&two)
            }
            Self::Sphere(_) => d.scaled(&two),
            Self::Torus(s) => {
                let major_sq = Real::from(&s.major * &s.major);
                let q = &d.norm_squared() + &major_sq - Real::from(&s.minor * &s.minor);
                let radial = d.minus(&a.scaled(&(h / n)));
                d.scaled(&(Real::from_i64(4) * q))
                    .minus(&radial.scaled(&(Real::from_i64(8) * major_sq)))
            }
        }
    }

    /// The same surface moved by `pose`.
    pub fn transformed(&self, pose: &Pose) -> Surface {
        match self {
            Self::Plane(s) => Self::Plane(Plane {
                origin: pose.apply_point(&s.origin),
                normal: pose.apply_vector(&s.normal),
                ref_dir: pose.apply_vector(&s.ref_dir),
            }),
            Self::Cylinder(s) => Self::Cylinder(Cylinder {
                origin: pose.apply_point(&s.origin),
                axis: pose.apply_vector(&s.axis),
                radius: s.radius.clone(),
                ref_dir: pose.apply_vector(&s.ref_dir),
            }),
            Self::Cone(s) => Self::Cone(Cone {
                apex: pose.apply_point(&s.apex),
                axis: pose.apply_vector(&s.axis),
                slope: s.slope.clone(),
                ref_dir: pose.apply_vector(&s.ref_dir),
            }),
            Self::Sphere(s) => Self::Sphere(Sphere {
                center: pose.apply_point(&s.center),
                radius: s.radius.clone(),
                axis: pose.apply_vector(&s.axis),
                ref_dir: pose.apply_vector(&s.ref_dir),
            }),
            Self::Torus(s) => Self::Torus(Torus {
                center: pose.apply_point(&s.center),
                axis: pose.apply_vector(&s.axis),
                major: s.major.clone(),
                minor: s.minor.clone(),
                ref_dir: pose.apply_vector(&s.ref_dir),
            }),
        }
    }

    /// Kinds of the `u` and `v` parameters.
    pub fn domain(&self) -> [ParamKind; 2] {
        match self {
            Self::Plane(_) => [ParamKind::Line, ParamKind::Line],
            Self::Cylinder(_) => [ParamKind::Periodic, ParamKind::Line],
            Self::Cone(_) => [ParamKind::Periodic, ParamKind::Ray],
            Self::Sphere(_) => [ParamKind::Periodic, ParamKind::Latitude],
            Self::Torus(_) => [ParamKind::Periodic, ParamKind::Periodic],
        }
    }

    /// Enclosure of the point at parameters `(u, v)`.
    pub fn evaluate(&self, u: &Ball, v: &Ball) -> [Ball; 3] {
        let precision = u.precision().max(v.precision());
        let (origin, axis, reference) = self.definition_frame();
        let frame = math::frame(axis, reference, precision);
        let r = |x: &Rational| math::ball(x, precision);
        let zero = || r(&Rational::zero());
        let coefficients = match self {
            Self::Plane(_) => [u.clone(), v.clone(), zero()],
            Self::Cylinder(s) => [r(&s.radius) * u.cos(), r(&s.radius) * u.sin(), v.clone()],
            Self::Cone(s) => {
                let radial = v * &r(&s.slope);
                [&radial * &u.cos(), &radial * &u.sin(), v.clone()]
            }
            Self::Sphere(s) => {
                let radius = r(&s.radius);
                let radial = &radius * &v.cos();
                [&radial * &u.cos(), &radial * &u.sin(), radius * v.sin()]
            }
            Self::Torus(s) => {
                let radial = r(&s.major) + r(&s.minor) * v.cos();
                [
                    &radial * &u.cos(),
                    &radial * &u.sin(),
                    r(&s.minor) * v.sin(),
                ]
            }
        };
        let offset = math::add(
            &math::add(
                &math::scale(&frame[0], &coefficients[0]),
                &math::scale(&frame[1], &coefficients[1]),
            ),
            &math::scale(&frame[2], &coefficients[2]),
        );
        math::add(&math::point(origin, precision), &offset)
    }

    /// Enclosures of the derivatives `[dS/du, dS/dv]` of the parametrisation
    /// at `(u, v)`.
    pub fn derivatives(&self, u: &Ball, v: &Ball) -> [[Ball; 3]; 2] {
        let precision = u.precision().max(v.precision());
        let (_, axis, reference) = self.definition_frame();
        let [e1, e2, e3] = math::frame(axis, reference, precision);
        if matches!(self, Self::Plane(_)) {
            return [e1, e2];
        }
        let (sin_u, cos_u) = (u.sin(), u.cos());
        let tangent = math::add(&math::scale(&e1, &(-&sin_u)), &math::scale(&e2, &cos_u));
        let r = |x: &Rational| math::ball(x, precision);
        if let Self::Cylinder(s) = self {
            return [math::scale(&tangent, &r(&s.radius)), e3];
        }
        let radial = math::add(&math::scale(&e1, &cos_u), &math::scale(&e2, &sin_u));
        match self {
            Self::Plane(_) | Self::Cylinder(_) => unreachable!(),
            Self::Cone(s) => {
                let slope = r(&s.slope);
                [
                    math::scale(&tangent, &(v * &slope)),
                    math::add(&e3, &math::scale(&radial, &slope)),
                ]
            }
            Self::Sphere(s) => {
                let radius = r(&s.radius);
                let cos_v = v.cos();
                [
                    math::scale(&tangent, &(&radius * &cos_v)),
                    math::scale(
                        &math::add(
                            &math::scale(&radial, &(-v.sin())),
                            &math::scale(&e3, &cos_v),
                        ),
                        &radius,
                    ),
                ]
            }
            Self::Torus(s) => {
                let minor = r(&s.minor);
                let cos_v = v.cos();
                [
                    math::scale(&tangent, &(r(&s.major) + &minor * &cos_v)),
                    math::scale(
                        &math::add(
                            &math::scale(&radial, &(-v.sin())),
                            &math::scale(&e3, &cos_v),
                        ),
                        &minor,
                    ),
                ]
            }
        }
    }

    /// Enclosure of the parameters of a point enclosed by `p`, by the formulas
    /// of the piece description; `None` when the box meets a seam, the axis,
    /// a pole or another place where those formulas jump or are undefined.
    pub fn parameters_of(&self, p: &[Ball; 3]) -> Option<[Ball; 2]> {
        let precision = p.iter().map(Ball::precision).max().unwrap();
        let work = precision.max(64).saturating_add(32);
        let (origin, axis, reference) = self.definition_frame();
        let x = math::project(p, origin, reference, work);
        let second = axis.cross(reference);
        let y = math::project(p, origin, &second, work);
        if matches!(self, Self::Plane(_)) {
            return Some([x.with_precision(precision), y.with_precision(precision)]);
        }
        let z = math::project(p, origin, axis, work);
        // Refuse only when a known point of the box lies on the closed seam.
        // Uncertain interval projections instead receive a full angle enclosure.
        let mid = RPoint3::new(p[0].mid(), p[1].mid(), p[2].mid()).minus(origin);
        if mid.dot(&second).is_zero() && mid.dot(reference).sign() != Sign::Negative {
            return None;
        }
        if let Self::Torus(s) = self
            && mid.dot(axis).is_zero()
            && mid.norm_squared() >= &s.major * &s.major
        {
            return None;
        }
        let u = angle(&x, &y);
        if matches!(self, Self::Cylinder(_) | Self::Cone(_)) {
            return Some([u.with_precision(precision), z.with_precision(precision)]);
        }
        let rho = (x.square() + y.square()).sqrt().expect("sum of squares");
        let v = match self {
            Self::Plane(_) | Self::Cylinder(_) | Self::Cone(_) => unreachable!(),
            Self::Sphere(_) => Ball::atan2(&z, &rho).unwrap_or_else(|| {
                let half = Ball::pi(work) * math::ball(&Rational::from_ratio(1, 2).unwrap(), work);
                half.hull(&(-&half))
            }),
            Self::Torus(s) => angle(&(rho - math::ball(&s.major, work)), &z),
        };
        Some([u.with_precision(precision), v.with_precision(precision)])
    }

    pub(crate) fn definition_frame(&self) -> (&RPoint3, &RVector3, &RVector3) {
        match self {
            Self::Plane(s) => (&s.origin, &s.normal, &s.ref_dir),
            Self::Cylinder(s) => (&s.origin, &s.axis, &s.ref_dir),
            Self::Cone(s) => (&s.apex, &s.axis, &s.ref_dir),
            Self::Sphere(s) => (&s.center, &s.axis, &s.ref_dir),
            Self::Torus(s) => (&s.center, &s.axis, &s.ref_dir),
        }
    }

    pub(crate) fn value(&self, p: &Point3) -> Real {
        let (origin, axis, _) = self.definition_frame();
        let d = p.minus(&origin.to_real());
        let a = axis.to_real();
        match self {
            Self::Plane(_) => d.dot(&a),
            Self::Cylinder(s) => {
                d.cross(&a).norm_squared() - Real::from(&s.radius * &s.radius * axis.norm_squared())
            }
            Self::Cone(s) => {
                d.cross(&a).norm_squared() - Real::from(&s.slope * &s.slope) * d.dot(&a).square()
            }
            Self::Sphere(s) => d.norm_squared() - Real::from(&s.radius * &s.radius),
            Self::Torus(s) => {
                let q = d.norm_squared();
                let major_sq = Real::from(&s.major * &s.major);
                let shifted = &q + &major_sq - Real::from(&s.minor * &s.minor);
                shifted.square()
                    - Real::from_i64(4)
                        * major_sq
                        * (q - d.dot(&a).square() / Real::from(axis.norm_squared()))
            }
        }
    }

    pub(crate) fn distance_bound(&self, b: &Rational) -> Rational {
        match self {
            Self::Plane(s) => b * &math::inverse_length_upper(&s.normal),
            Self::Cylinder(s) => b / &(&s.axis.norm_squared() * &s.radius),
            Self::Cone(s) => math::sqrt_upper(&(b / &s.axis.norm_squared())),
            Self::Sphere(s) => b / &s.radius,
            Self::Torus(s) => b / &(&s.minor * &(&s.major * &s.major - &s.minor * &s.minor)),
        }
    }
}

fn angle(x: &Ball, y: &Ball) -> Ball {
    let precision = x.precision().max(y.precision());
    if x.sign() == Some(StrictSign::Negative) {
        // This branch avoids atan2's negative-x jump before shifting to [0, 2pi).
        return Ball::pi(precision) - Ball::atan2(y, &(-x)).expect("positive x coordinate");
    }
    match y.sign() {
        Some(StrictSign::Positive) => Ball::atan2(y, x).expect("nonzero y coordinate"),
        sign => {
            let period = Ball::pi(precision) * math::ball(&Rational::from_i64(2), precision);
            if sign == Some(StrictSign::Negative) {
                Ball::atan2(y, x).expect("nonzero y coordinate") + period
            } else {
                math::ball(&Rational::zero(), precision).hull(&period)
            }
        }
    }
}
