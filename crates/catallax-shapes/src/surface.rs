//! Analytic and certified finite parametric surface definitions.

use crate::{Point3, Pose, RPoint3, SplineSupport, Vector3, math};
use catallax_numbers::{Ball, Decider, Rational, Real, RealPoly, Refusal, Sign, StrictSign};

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
    /// Outer sheet of a horn or spindle torus.
    TorusSheet,
    /// Elliptic cylinder.
    EllipticCylinder,
    /// Certified finite B-spline support.
    BSpline,
}

/// The domain of a surface parameter.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum ParamKind {
    /// All reals.
    Line,
    /// Period `2 pi`, seam at zero.
    Periodic,
    /// `[-pi/2, pi/2]`.
    Latitude,
    /// `[0, infinity)`.
    Ray,
    /// A finite closed interval, given by the spline's exact knot domain.
    Finite,
    /// A finite interval whose endpoints are identified by a checked transition.
    PeriodicFinite,
}

/// A plane, oriented by its normal.
#[derive(Clone, Debug)]
pub struct Plane {
    /// Parameter origin.
    pub origin: Point3,
    /// Nonzero normal.
    pub normal: Vector3,
    /// Nonzero in-plane reference direction.
    pub ref_dir: Vector3,
}

/// A circular cylinder, oriented away from its axis.
#[derive(Clone, Debug)]
pub struct Cylinder {
    /// Axis point where `v = 0`.
    pub origin: Point3,
    /// Nonzero axis direction.
    pub axis: Vector3,
    /// Positive squared radius, retained without extracting a root.
    pub radius_sq: Real,
    /// Nonzero perpendicular seam direction.
    pub ref_dir: Vector3,
}

/// The nappe of a circular cone in its axis direction.
#[derive(Clone, Debug)]
pub struct Cone {
    /// Apex.
    pub apex: Point3,
    /// Nonzero direction into the nappe.
    pub axis: Vector3,
    /// Positive radius gained per unit axial distance.
    pub slope: Real,
    /// Nonzero perpendicular seam direction.
    pub ref_dir: Vector3,
}

/// A sphere, oriented outwards.
#[derive(Clone, Debug)]
pub struct Sphere {
    /// Centre.
    pub center: Point3,
    /// Positive radius.
    pub radius: Real,
    /// Nonzero pole direction.
    pub axis: Vector3,
    /// Nonzero perpendicular seam direction.
    pub ref_dir: Vector3,
}

/// A ring torus with `0 < minor < major`.
#[derive(Clone, Debug)]
pub struct Torus {
    /// Centre.
    pub center: Point3,
    /// Nonzero axis direction.
    pub axis: Vector3,
    /// Positive core-circle radius.
    pub major: Real,
    /// Positive tube radius, less than the major radius.
    pub minor: Real,
    /// Nonzero perpendicular seam direction.
    pub ref_dir: Vector3,
}

/// The outer sheet of a horn or spindle torus, `0 < major <= minor`: the
/// points `center + (major + minor cos v) e(u) + minor sin v axis` with
/// `major + minor cos v > 0`, oriented away from the tube centre circle. The
/// self-intersecting inner sheet is not part of this support, so its side is
/// that of the outer tube, `(rho - major)^2 + z^2 - minor^2`, not the sign of
/// the complete quartic carrier. It has no global polynomial equation.
#[derive(Clone, Debug)]
pub struct TorusSheet {
    /// Centre.
    pub center: Point3,
    /// Nonzero axis direction.
    pub axis: Vector3,
    /// Positive core-circle radius, at most the tube radius.
    pub major: Real,
    /// Positive tube radius.
    pub minor: Real,
    /// Nonzero perpendicular seam direction.
    pub ref_dir: Vector3,
}

/// An elliptic cylinder, oriented away from its axis.
#[derive(Clone, Debug)]
pub struct EllipticCylinder {
    /// Axis point where `v = 0`.
    pub origin: Point3,
    /// Nonzero axis direction.
    pub axis: Vector3,
    /// Nonzero perpendicular major-axis and seam direction.
    pub major_dir: Vector3,
    /// Squared major semi-axis.
    pub major_sq: Real,
    /// Squared minor semi-axis, with `0 < minor_sq < major_sq`.
    pub minor_sq: Real,
}

impl Plane {
    /// Validates field domains, direction lengths and perpendicularity, in order.
    pub fn new(
        origin: Point3,
        normal: Vector3,
        ref_dir: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&origin, dec)?;
        math::validate_vector(&normal, dec)?;
        math::validate_vector(&ref_dir, dec)?;
        math::directions(&normal, &ref_dir, dec)?;
        math::perpendicular(&normal, &ref_dir, dec)?;
        Ok(Self {
            origin,
            normal,
            ref_dir,
        })
    }
}
impl Cylinder {
    /// Validates field domains before proving the cylinder invariants.
    pub fn new(
        origin: Point3,
        axis: Vector3,
        radius: Real,
        ref_dir: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&origin, dec)?;
        math::validate_vector(&axis, dec)?;
        math::validate(&radius, dec)?;
        math::validate_vector(&ref_dir, dec)?;
        math::directions(&axis, &ref_dir, dec)?;
        math::positive(&radius, "shapes.non-positive-radius", dec)?;
        math::perpendicular(&axis, &ref_dir, dec)?;
        Ok(Self {
            origin,
            axis,
            radius_sq: radius.square(),
            ref_dir,
        })
    }

    /// Validates a cylinder directly from its positive squared radius.
    pub fn from_radius_sq(
        origin: Point3,
        axis: Vector3,
        radius_sq: Real,
        ref_dir: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&origin, dec)?;
        math::validate_vector(&axis, dec)?;
        math::validate(&radius_sq, dec)?;
        math::validate_vector(&ref_dir, dec)?;
        math::directions(&axis, &ref_dir, dec)?;
        math::positive(&radius_sq, "shapes.non-positive-radius", dec)?;
        math::perpendicular(&axis, &ref_dir, dec)?;
        Ok(Self {
            origin,
            axis,
            radius_sq,
            ref_dir,
        })
    }

    /// Positive radius for metric calculations and exchange formats.
    pub fn radius(&self) -> Real {
        math::exact_sqrt(&self.radius_sq)
    }

    /// A point-on-circle direction, without roots when the seam already has
    /// the defining radius, as it does for an extruded circular arc.
    pub fn radial_reference(&self) -> Vector3 {
        self.ref_dir.scaled(&math::exact_sqrt(
            &(&self.radius_sq / &self.ref_dir.norm_squared()),
        ))
    }
}
impl Cone {
    /// Validates field domains before proving the cone invariants.
    pub fn new(
        apex: Point3,
        axis: Vector3,
        slope: Real,
        ref_dir: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&apex, dec)?;
        math::validate_vector(&axis, dec)?;
        math::validate(&slope, dec)?;
        math::validate_vector(&ref_dir, dec)?;
        math::directions(&axis, &ref_dir, dec)?;
        math::positive(&slope, "shapes.non-positive-slope", dec)?;
        math::perpendicular(&axis, &ref_dir, dec)?;
        Ok(Self {
            apex,
            axis,
            slope,
            ref_dir,
        })
    }
}
impl Sphere {
    /// Validates field domains before proving the sphere invariants.
    pub fn new(
        center: Point3,
        radius: Real,
        axis: Vector3,
        ref_dir: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&center, dec)?;
        math::validate(&radius, dec)?;
        math::validate_vector(&axis, dec)?;
        math::validate_vector(&ref_dir, dec)?;
        math::directions(&axis, &ref_dir, dec)?;
        math::positive(&radius, "shapes.non-positive-radius", dec)?;
        math::perpendicular(&axis, &ref_dir, dec)?;
        Ok(Self {
            center,
            radius,
            axis,
            ref_dir,
        })
    }
}
impl Torus {
    /// Validates field domains before proving the ring-torus invariants.
    pub fn new(
        center: Point3,
        axis: Vector3,
        major: Real,
        minor: Real,
        ref_dir: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&center, dec)?;
        math::validate_vector(&axis, dec)?;
        math::validate(&major, dec)?;
        math::validate(&minor, dec)?;
        math::validate_vector(&ref_dir, dec)?;
        math::directions(&axis, &ref_dir, dec)?;
        math::positive(&major, "shapes.non-positive-radius", dec)?;
        math::positive(&minor, "shapes.non-positive-radius", dec)?;
        if math::sign(&(&major - &minor), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "shapes.torus-self-intersecting",
                "minor radius must be less than major radius",
            ));
        }
        math::perpendicular(&axis, &ref_dir, dec)?;
        Ok(Self {
            center,
            axis,
            major,
            minor,
            ref_dir,
        })
    }
}
impl TorusSheet {
    /// Validates field domains before proving `0 < major <= minor`.
    pub fn new(
        center: Point3,
        axis: Vector3,
        major: Real,
        minor: Real,
        ref_dir: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&center, dec)?;
        math::validate_vector(&axis, dec)?;
        math::validate(&major, dec)?;
        math::validate(&minor, dec)?;
        math::validate_vector(&ref_dir, dec)?;
        math::directions(&axis, &ref_dir, dec)?;
        math::positive(&major, "shapes.non-positive-radius", dec)?;
        math::positive(&minor, "shapes.non-positive-radius", dec)?;
        if math::sign(&(&minor - &major), dec)? == Sign::Negative {
            return Err(Refusal::new(
                "shapes.torus-sheet-ring",
                "a ring torus is a complete torus, not an outer sheet",
            ));
        }
        math::perpendicular(&axis, &ref_dir, dec)?;
        Ok(Self {
            center,
            axis,
            major,
            minor,
            ref_dir,
        })
    }

    /// Axial height and distance from the axis, both relative to the centre.
    fn cylindrical(&self, p: &Point3) -> (Vector3, Real, Real) {
        let d = p.minus(&self.center);
        let norm = self.axis.norm_squared();
        let height = d.dot(&self.axis);
        let radial = d.minus(&self.axis.scaled(&(&height / &norm)));
        let rho = radial.norm_squared().sqrt();
        (radial, rho, &height / &norm.sqrt())
    }

    /// Outward normal covector away from the tube centre on this sheet's side.
    pub fn normal(&self, p: &Point3) -> Vector3 {
        let (radial, rho, _) = self.cylindrical(p);
        let d = p.minus(&self.center);
        let norm = self.axis.norm_squared();
        let axial = self.axis.scaled(&(&d.dot(&self.axis) / &norm));
        radial.scaled(&(&(&rho - &self.major) / &rho)).plus(&axial)
    }
}

impl EllipticCylinder {
    /// Validates the strict ellipse axes. Equal axes must be reduced to a
    /// circular cylinder by the defining operation before calling this constructor.
    pub fn new(
        origin: Point3,
        axis: Vector3,
        major_dir: Vector3,
        major_sq: Real,
        minor_sq: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&origin, dec)?;
        math::validate_vector(&axis, dec)?;
        math::validate_vector(&major_dir, dec)?;
        math::validate(&major_sq, dec)?;
        math::validate(&minor_sq, dec)?;
        math::directions(&axis, &major_dir, dec)?;
        if math::sign(&minor_sq, dec)? != Sign::Positive
            || math::sign(&(&major_sq - &minor_sq), dec)? != Sign::Positive
        {
            return Err(Refusal::new(
                "shapes.ellipse-axes",
                "elliptic cylinder requires 0 < minor_sq < major_sq",
            ));
        }
        math::perpendicular(&axis, &major_dir, dec)?;
        Ok(Self {
            origin,
            axis,
            major_dir,
            major_sq,
            minor_sq,
        })
    }
}

/// An analytic support or certified finite parametric support. Face trims belong to topology.
#[non_exhaustive]
#[derive(Clone, Debug)]
pub enum Surface {
    /// Plane.
    Plane(Plane),
    /// Circular cylinder.
    Cylinder(Cylinder),
    /// Circular cone.
    Cone(Cone),
    /// Sphere.
    Sphere(Sphere),
    /// Ring torus.
    Torus(Torus),
    /// Outer sheet of a horn or spindle torus.
    TorusSheet(TorusSheet),
    /// Elliptic cylinder.
    EllipticCylinder(EllipticCylinder),
    /// Genuine finite B-spline carrier with its producer's checked chart family.
    BSpline(SplineSupport),
}

impl Surface {
    /// Surface kind.
    pub fn kind(&self) -> SurfaceKind {
        match self {
            Self::Plane(_) => SurfaceKind::Plane,
            Self::Cylinder(_) => SurfaceKind::Cylinder,
            Self::Cone(_) => SurfaceKind::Cone,
            Self::Sphere(_) => SurfaceKind::Sphere,
            Self::Torus(_) => SurfaceKind::Torus,
            Self::TorusSheet(_) => SurfaceKind::TorusSheet,
            Self::EllipticCylinder(_) => SurfaceKind::EllipticCylinder,
            Self::BSpline(_) => SurfaceKind::BSpline,
        }
    }

    fn definition_fields(&self) -> Option<[Option<&Real>; 11]> {
        let (origin, axis, reference) = self.definition_frame()?;
        let dimensions = match self {
            Self::Plane(_) => [None, None],
            Self::Cylinder(s) => [Some(&s.radius_sq), None],
            Self::Cone(s) => [Some(&s.slope), None],
            Self::Sphere(s) => [Some(&s.radius), None],
            Self::Torus(s) => [Some(&s.major), Some(&s.minor)],
            Self::TorusSheet(s) => [Some(&s.major), Some(&s.minor)],
            Self::EllipticCylinder(s) => [Some(&s.major_sq), Some(&s.minor_sq)],
            Self::BSpline(_) => return None,
        };
        Some([
            Some(&origin.x),
            Some(&origin.y),
            Some(&origin.z),
            Some(&axis.x),
            Some(&axis.y),
            Some(&axis.z),
            Some(&reference.x),
            Some(&reference.y),
            Some(&reference.z),
            dimensions[0],
            dimensions[1],
        ])
    }

    /// Whether every defining scalar has a cheap exact rational value.
    pub fn has_rational_definition(&self) -> bool {
        if let Self::BSpline(spline) = self {
            return spline.carrier().has_rational_definition();
        }
        self.definition_fields().is_some_and(|fields| {
            fields
                .into_iter()
                .flatten()
                .all(|value| value.rational_ref().is_some())
        })
    }

    /// Safe construction-key agreement, not a geometry-equality predicate.
    /// Different definitions may still describe the same point set.
    pub fn same_definition(&self, other: &Self) -> bool {
        if self.kind() != other.kind() {
            return false;
        }
        if let (Self::BSpline(a), Self::BSpline(b)) = (self, other) {
            return a.carrier().same_definition(b.carrier());
        }
        match (self.definition_fields(), other.definition_fields()) {
            (Some(a), Some(b)) => a.into_iter().zip(b).all(|(a, b)| match (a, b) {
                (Some(a), Some(b)) => a.same_definition(b),
                (None, None) => true,
                _ => false,
            }),
            _ => false,
        }
    }

    /// Whether a cached result retains the same input expression-DAG nodes.
    /// Value-equal but separately constructed inputs must not acquire aliases
    /// through a cache: those aliases can change deterministic proof work.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        if self.kind() != other.kind() {
            return false;
        }
        if let (Self::BSpline(a), Self::BSpline(b)) = (self, other) {
            return a.shares_definition_nodes(b);
        }
        match (self.definition_fields(), other.definition_fields()) {
            (Some(a), Some(b)) => a.into_iter().zip(b).all(|(a, b)| match (a, b) {
                (Some(a), Some(b)) => a.same_node(b),
                (None, None) => true,
                _ => false,
            }),
            _ => false,
        }
    }

    /// Global analytic equation at its specified scaling. A finite multi-span
    /// spline has no such equation; absence is not the zero polynomial.
    pub fn implicit(&self) -> Option<RealPoly> {
        if let Self::TorusSheet(_) = self {
            return None;
        }
        let (origin, axis, _) = self.definition_frame()?;
        let d = math::displacement(origin);
        let constant = |r| RealPoly::constant(3, r);
        Some(match self {
            Self::Plane(_) => math::dot(&d, axis),
            Self::Cylinder(s) => {
                math::norm(&math::cross(&d, axis)) - constant(&s.radius_sq * &axis.norm_squared())
            }
            Self::Cone(s) => {
                let h = math::dot(&d, axis);
                math::norm(&math::cross(&d, axis)) - (&h * &h).scale(&s.slope.square())
            }
            Self::Sphere(s) => math::norm(&d) - constant(s.radius.square()),
            Self::Torus(s) => {
                let q = math::norm(&d);
                let h = math::dot(&d, axis);
                let radial = &q - &(&h * &h).scale(&(&Real::one() / &axis.norm_squared()));
                let a = q + constant(&s.major.square() - &s.minor.square());
                &a * &a - radial.scale(&(&Real::from_i64(4) * &s.major.square()))
            }
            Self::TorusSheet(_) => unreachable!("an outer sheet has no polynomial"),
            Self::EllipticCylinder(s) => {
                let second = s.axis.cross(&s.major_dir);
                let a = math::dot(&d, &s.major_dir);
                let b = math::dot(&d, &second);
                (&a * &a).scale(&(&Real::one() / &(&s.major_sq * &s.major_dir.norm_squared())))
                    + (&b * &b).scale(&(&Real::one() / &(&s.minor_sq * &second.norm_squared())))
                    - constant(Real::one())
            }
            Self::BSpline(_) => return None,
        })
    }

    /// Gradient of the analytic equation. Parametric normals instead require
    /// a finite chart preimage and are queried with `normal_at`.
    pub fn gradient(&self, p: &Point3) -> Option<Vector3> {
        if let Self::TorusSheet(s) = self {
            return Some(s.normal(p).scaled(&Real::from_i64(2)));
        }
        let (origin, axis, _) = self.definition_frame()?;
        let d = p.minus(origin);
        let n = axis.norm_squared();
        let h = d.dot(axis);
        let two = Real::from_i64(2);
        Some(match self {
            Self::Plane(_) => axis.clone(),
            Self::Cylinder(_) => d.scaled(&n).minus(&axis.scaled(&h)).scaled(&two),
            Self::Cone(s) => d
                .scaled(&n)
                .minus(&axis.scaled(&(&h * &(&Real::one() + &s.slope.square()))))
                .scaled(&two),
            Self::Sphere(_) => d.scaled(&two),
            Self::Torus(s) => {
                let major_sq = s.major.square();
                let q = &d.norm_squared() + &major_sq - s.minor.square();
                let radial = d.minus(&axis.scaled(&(&h / &n)));
                d.scaled(&(&Real::from_i64(4) * &q))
                    .minus(&radial.scaled(&(&Real::from_i64(8) * &major_sq)))
            }
            Self::TorusSheet(_) => unreachable!("handled above"),
            Self::EllipticCylinder(s) => {
                let second = s.axis.cross(&s.major_dir);
                s.major_dir
                    .scaled(&(&d.dot(&s.major_dir) / &(&s.major_sq * &s.major_dir.norm_squared())))
                    .plus(
                        &second
                            .scaled(&(&d.dot(&second) / &(&s.minor_sq * &second.norm_squared()))),
                    )
                    .scaled(&two)
            }
            Self::BSpline(_) => return None,
        })
    }

    /// Outward normal covector at a surface point, not necessarily a unit vector.
    /// Parametric supports prove finite incidence before returning their normal.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        match self {
            Self::BSpline(spline) => spline.normal_at(point, dec),
            _ => self.gradient(point).ok_or_else(no_global_equation),
        }
    }

    /// Applies a validated rigid pose, retaining dimensions exactly.
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
                radius_sq: s.radius_sq.clone(),
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
            Self::TorusSheet(s) => Self::TorusSheet(TorusSheet {
                center: pose.apply_point(&s.center),
                axis: pose.apply_vector(&s.axis),
                major: s.major.clone(),
                minor: s.minor.clone(),
                ref_dir: pose.apply_vector(&s.ref_dir),
            }),
            Self::EllipticCylinder(s) => Self::EllipticCylinder(EllipticCylinder {
                origin: pose.apply_point(&s.origin),
                axis: pose.apply_vector(&s.axis),
                major_dir: pose.apply_vector(&s.major_dir),
                major_sq: s.major_sq.clone(),
                minor_sq: s.minor_sq.clone(),
            }),
            Self::BSpline(spline) => Self::BSpline(spline.transformed(pose)),
        }
    }

    /// Parameter domains.
    pub fn domain(&self) -> [ParamKind; 2] {
        match self {
            Self::Plane(_) => [ParamKind::Line, ParamKind::Line],
            Self::Cylinder(_) | Self::EllipticCylinder(_) => [ParamKind::Periodic, ParamKind::Line],
            Self::Cone(_) => [ParamKind::Periodic, ParamKind::Ray],
            Self::Sphere(_) => [ParamKind::Periodic, ParamKind::Latitude],
            Self::Torus(_) | Self::TorusSheet(_) => [ParamKind::Periodic, ParamKind::Periodic],
            Self::BSpline(spline) => spline.param_kinds(),
        }
    }

    /// Certified chart enclosure at the input precision.
    pub fn evaluate(&self, u: &Ball, v: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        if let Self::BSpline(spline) = self {
            return spline.carrier().evaluate(u, v, dec);
        }
        let precision = u.precision().max(v.precision());
        let (origin, axis, reference) = self.definition_frame().ok_or_else(no_global_equation)?;
        let frame = math::frame(axis, reference, precision, dec)?;
        let r = |x: &Real| dec.enclose_at(x, precision);
        let zero = || math::ball(&Rational::zero(), precision);
        let coefficients = match self {
            Self::Plane(_) => [u.clone(), v.clone(), zero()],
            Self::Cylinder(s) => {
                let radius = r(&s.radius())?;
                [&radius * &u.cos(), &radius * &u.sin(), v.clone()]
            }
            Self::Cone(s) => {
                let radial = v * &r(&s.slope)?;
                [&radial * &u.cos(), &radial * &u.sin(), v.clone()]
            }
            Self::Sphere(s) => {
                let radius = r(&s.radius)?;
                let radial = &radius * &v.cos();
                [&radial * &u.cos(), &radial * &u.sin(), radius * v.sin()]
            }
            Self::Torus(Torus { major, minor, .. })
            | Self::TorusSheet(TorusSheet { major, minor, .. }) => {
                let radial = r(major)? + r(minor)? * v.cos();
                [&radial * &u.cos(), &radial * &u.sin(), r(minor)? * v.sin()]
            }
            Self::EllipticCylinder(s) => [
                r(&s.major_sq.sqrt())? * u.cos(),
                r(&s.minor_sq.sqrt())? * u.sin(),
                v.clone(),
            ],
            Self::BSpline(spline) => return spline.carrier().evaluate(u, v, dec),
        };
        let offset = math::add(
            &math::add(
                &math::scale(&frame[0], &coefficients[0]),
                &math::scale(&frame[1], &coefficients[1]),
            ),
            &math::scale(&frame[2], &coefficients[2]),
        );
        Ok(math::add(&math::point(origin, precision, dec)?, &offset))
    }

    /// Certified chart derivatives `[dS/du, dS/dv]`.
    pub fn derivatives(
        &self,
        u: &Ball,
        v: &Ball,
        dec: &Decider,
    ) -> Result<[[Ball; 3]; 2], Refusal> {
        if let Self::BSpline(spline) = self {
            return spline.carrier().derivatives(u, v, dec);
        }
        let precision = u.precision().max(v.precision());
        let (_, axis, reference) = self.definition_frame().ok_or_else(no_global_equation)?;
        let [e1, e2, e3] = math::frame(axis, reference, precision, dec)?;
        if matches!(self, Self::Plane(_)) {
            return Ok([e1, e2]);
        }
        let (sin_u, cos_u) = (u.sin(), u.cos());
        let tangent = math::add(&math::scale(&e1, &(-&sin_u)), &math::scale(&e2, &cos_u));
        let r = |x: &Real| dec.enclose_at(x, precision);
        if let Self::Cylinder(s) = self {
            return Ok([math::scale(&tangent, &r(&s.radius())?), e3]);
        }
        if let Self::EllipticCylinder(s) = self {
            return Ok([
                math::add(
                    &math::scale(&e1, &(-r(&s.major_sq.sqrt())? * sin_u)),
                    &math::scale(&e2, &(r(&s.minor_sq.sqrt())? * cos_u)),
                ),
                e3,
            ]);
        }
        let radial = math::add(&math::scale(&e1, &cos_u), &math::scale(&e2, &sin_u));
        Ok(match self {
            Self::Plane(_) | Self::Cylinder(_) | Self::EllipticCylinder(_) => unreachable!(),
            Self::Cone(s) => {
                let slope = r(&s.slope)?;
                [
                    math::scale(&tangent, &(v * &slope)),
                    math::add(&e3, &math::scale(&radial, &slope)),
                ]
            }
            Self::Sphere(s) => {
                let radius = r(&s.radius)?;
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
            Self::Torus(Torus { major, minor, .. })
            | Self::TorusSheet(TorusSheet { major, minor, .. }) => {
                let minor = r(minor)?;
                let cos_v = v.cos();
                [
                    math::scale(&tangent, &(r(major)? + &minor * &cos_v)),
                    math::scale(
                        &math::add(
                            &math::scale(&radial, &(-v.sin())),
                            &math::scale(&e3, &cos_v),
                        ),
                        &minor,
                    ),
                ]
            }
            Self::BSpline(spline) => return spline.carrier().derivatives(u, v, dec),
        })
    }

    /// Chart inverse. Seam/domain ambiguity is distinct from an enclosure failure.
    pub fn parameters_of(
        &self,
        p: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        if let Self::BSpline(spline) = self {
            return spline.parameters_of(p, dec);
        }
        let precision = p.iter().map(Ball::precision).max().unwrap();
        let work = precision.max(64).saturating_add(32);
        let (origin, axis, reference) = self.definition_frame().ok_or_else(no_global_equation)?;
        let mut x = math::project(p, origin, reference, work, dec)?;
        let second = axis.cross(reference);
        let mut y = math::project(p, origin, &second, work, dec)?;
        if matches!(self, Self::Plane(_)) {
            return Ok(Some([
                x.with_precision(precision),
                y.with_precision(precision),
            ]));
        }
        let z = math::project(p, origin, axis, work, dec)?;
        if let (Some(origin), Some(reference), Some(axis), Some(second)) = (
            origin.as_rational(),
            reference.as_rational(),
            axis.as_rational(),
            second.as_rational(),
        ) {
            let mid = RPoint3::new(p[0].mid(), p[1].mid(), p[2].mid()).minus(&origin);
            if mid.dot(&second).is_zero() && mid.dot(&reference).sign() != Sign::Negative {
                return Ok(None);
            }
            if let Self::Torus(Torus { major, .. }) | Self::TorusSheet(TorusSheet { major, .. }) =
                self
                && let Some(major) = major.rational_ref()
                && mid.dot(&axis).is_zero()
                && mid.norm_squared() >= major * major
            {
                return Ok(None);
            }
        } else if y.lo() <= Rational::zero()
            && y.hi() >= Rational::zero()
            && x.hi() >= Rational::zero()
        {
            return Ok(None);
        }
        if let Self::EllipticCylinder(s) = self {
            x = x
                .checked_div(&dec.enclose_at(&s.major_sq.sqrt(), work)?)
                .ok_or_else(|| {
                    Refusal::new(
                        "numbers.divisor-not-separated",
                        "major semi-axis enclosure contains zero",
                    )
                })?;
            y = y
                .checked_div(&dec.enclose_at(&s.minor_sq.sqrt(), work)?)
                .ok_or_else(|| {
                    Refusal::new(
                        "numbers.divisor-not-separated",
                        "minor semi-axis enclosure contains zero",
                    )
                })?;
        }
        let u = angle(&x, &y);
        if matches!(
            self,
            Self::Cylinder(_) | Self::Cone(_) | Self::EllipticCylinder(_)
        ) {
            return Ok(Some([
                u.with_precision(precision),
                z.with_precision(precision),
            ]));
        }
        let rho = (x.square() + y.square()).sqrt().expect("sum of squares");
        let v = match self {
            Self::Sphere(_) => Ball::atan2(&z, &rho).unwrap_or_else(|| {
                let half = Ball::pi(work) * math::ball(&Rational::from_ratio(1, 2).unwrap(), work);
                half.hull(&(-&half))
            }),
            Self::Torus(Torus { major, .. }) | Self::TorusSheet(TorusSheet { major, .. }) => {
                angle(&(rho - dec.enclose_at(major, work)?), &z)
            }
            _ => unreachable!(),
        };
        Ok(Some([
            u.with_precision(precision),
            v.with_precision(precision),
        ]))
    }

    pub(crate) fn definition_frame(&self) -> Option<(&Point3, &Vector3, &Vector3)> {
        Some(match self {
            Self::Plane(s) => (&s.origin, &s.normal, &s.ref_dir),
            Self::Cylinder(s) => (&s.origin, &s.axis, &s.ref_dir),
            Self::Cone(s) => (&s.apex, &s.axis, &s.ref_dir),
            Self::Sphere(s) => (&s.center, &s.axis, &s.ref_dir),
            Self::Torus(s) => (&s.center, &s.axis, &s.ref_dir),
            Self::TorusSheet(s) => (&s.center, &s.axis, &s.ref_dir),
            Self::EllipticCylinder(s) => (&s.origin, &s.axis, &s.major_dir),
            Self::BSpline(_) => return None,
        })
    }

    pub(crate) fn value(&self, p: &Point3) -> Option<Real> {
        let (origin, axis, _) = self.definition_frame()?;
        let d = p.minus(origin);
        Some(match self {
            Self::Plane(_) => d.dot(axis),
            Self::Cylinder(s) => d.cross(axis).norm_squared() - &s.radius_sq * &axis.norm_squared(),
            Self::Cone(s) => d.cross(axis).norm_squared() - s.slope.square() * d.dot(axis).square(),
            Self::Sphere(s) => d.norm_squared() - s.radius.square(),
            Self::Torus(s) => {
                let q = d.norm_squared();
                let major_sq = s.major.square();
                let shifted = &q + &major_sq - s.minor.square();
                shifted.square()
                    - Real::from_i64(4)
                        * major_sq
                        * (q - d.dot(axis).square() / axis.norm_squared())
            }
            Self::TorusSheet(s) => {
                let (_, rho, z) = s.cylindrical(p);
                (&rho - &s.major).square() + z.square() - s.minor.square()
            }
            Self::EllipticCylinder(s) => {
                let second = s.axis.cross(&s.major_dir);
                &(&d.dot(&s.major_dir).square() / &(&s.major_sq * &s.major_dir.norm_squared()))
                    + &(&d.dot(&second).square() / &(&s.minor_sq * &second.norm_squared()))
                    - Real::one()
            }
            Self::BSpline(_) => return None,
        })
    }

    pub(crate) fn distance_bound(&self, b: &Rational, dec: &Decider) -> Result<Rational, Refusal> {
        match self {
            Self::Plane(s) => {
                if let Some(norm) = math::rational_norm_squared(&s.normal) {
                    return Ok(b * &math::sqrt_upper(&norm.recip().expect("nonzero normal")));
                }
            }
            Self::Cylinder(s) => {
                if let Some(norm) = math::rational_norm_squared(&s.axis)
                    && let Some(radius) = s.radius().rational_ref()
                {
                    return Ok(b / &(&norm * radius));
                }
            }
            Self::Cone(s) => {
                if let Some(norm) = math::rational_norm_squared(&s.axis) {
                    return Ok(math::sqrt_upper(&(b / &norm)));
                }
            }
            Self::Sphere(s) => {
                if let Some(radius) = s.radius.rational_ref() {
                    return Ok(b / radius);
                }
            }
            Self::Torus(s) => {
                if let (Some(major), Some(minor)) = (s.major.rational_ref(), s.minor.rational_ref())
                {
                    return Ok(b / &(minor * &(major * major - minor * minor)));
                }
            }
            // In a meridian half-plane the residual is `D^2 - minor^2` for the
            // distance `D` to the tube centre. Where the nearest tube point is
            // retained, `|D - minor| <= b/minor`. Elsewhere the point lies in
            // the triangle between the tube centre and the axis points
            // `(0, +-h)`, `h^2 = minor^2 - major^2`, and there `b >= h|p - e|`
            // for the nearer axis point `e`. A horn sheet has no such triangle.
            Self::TorusSheet(s) => {
                if let (Some(major), Some(minor)) = (s.major.rational_ref(), s.minor.rational_ref())
                {
                    let h_sq = minor * minor - major * major;
                    return Ok(match h_sq.recip() {
                        Some(inverse) => b * &math::sqrt_upper(&inverse),
                        None => b / minor,
                    });
                }
            }
            Self::EllipticCylinder(_) => {}
            Self::BSpline(_) => return Err(no_global_equation()),
        }
        let bound = Real::from(b.clone());
        let value = match self {
            Self::Plane(s) => &bound / &s.normal.norm_squared().sqrt(),
            Self::Cylinder(s) => &bound / &(&s.axis.norm_squared() * &s.radius()),
            Self::Cone(s) => (&bound / &s.axis.norm_squared()).sqrt(),
            Self::Sphere(s) => &bound / &s.radius,
            Self::Torus(s) => &bound / &(&s.minor * &(&s.major.square() - &s.minor.square())),
            Self::TorusSheet(s) => &bound / &(s.minor.square() - s.major.square()).sqrt(),
            Self::EllipticCylinder(s) => &bound * &s.major_sq.sqrt(),
            Self::BSpline(_) => return Err(no_global_equation()),
        };
        if let Some(value) = value.as_rational() {
            return Ok(value);
        }
        Ok(dec.enclose_at(&value, 64)?.hi())
    }
}

fn angle(x: &Ball, y: &Ball) -> Ball {
    let precision = x.precision().max(y.precision());
    if x.sign() == Some(StrictSign::Negative) {
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

fn no_global_equation() -> Refusal {
    Refusal::new(
        "shapes.parametric-surface",
        "this query requires a global analytic support equation",
    )
}
