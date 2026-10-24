//! Exact distance-preserving maps, including orientation-reversing reflections.

use crate::{
    Circle, Cone, Curve, Cylinder, Ellipse, EllipticCylinder, Line, Plane, Point3, Pose, RPoint3,
    RVector3, Rotation, Sphere, Surface, Torus, Vector3,
};
use catallax_numbers::{Decider, Rational, Real, Refusal, sin_cos_degrees};

/// An orthogonal affine map. Constructors establish orthogonality without fitting.
#[derive(Clone, Debug)]
pub struct Isometry {
    matrix: [[Real; 3]; 3],
    /// Translation applied after the orthogonal linear map.
    pub translation: Vector3,
    reversing: bool,
}

impl From<&Pose> for Isometry {
    fn from(pose: &Pose) -> Self {
        Self {
            matrix: pose.rotation.matrix().clone(),
            translation: pose.translation.clone(),
            reversing: false,
        }
    }
}

impl Isometry {
    /// Identity map.
    pub fn identity() -> Self {
        Self::from(&Pose::identity())
    }

    /// Exact translation.
    pub fn translation(vector: RVector3) -> Self {
        Self::from(&Pose::from_rational(Rotation::identity(), vector))
    }

    /// Rotation by a rational number of degrees about a rational axis.
    /// Trigonometric values and axis normalisation remain exact algebraic values.
    pub fn around_axis(
        point: &RPoint3,
        direction: &RVector3,
        angle: &Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if direction.is_zero() {
            return Err(Refusal::new(
                "pattern.zero-axis",
                "circular pattern axis must not be zero",
            ));
        }
        let [s, c] = sin_cos_degrees(angle, dec)?;
        let norm = direction.norm_squared();
        let root = Real::from(norm.clone()).sqrt();
        let root = dec.rational_value(&root)?.map(Real::from).unwrap_or(root);
        let n = [&direction.x, &direction.y, &direction.z];
        let cross = [
            [Rational::zero(), -&direction.z, direction.y.clone()],
            [direction.z.clone(), Rational::zero(), -&direction.x],
            [-&direction.y, direction.x.clone(), Rational::zero()],
        ];
        let mut matrix = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let diagonal = if i == j { c.clone() } else { Real::zero() };
                diagonal
                    + (&Real::one() - &c) * Real::from(&(n[i] * n[j]) / &norm)
                    + &(&s * &Real::from(cross[i][j].clone())) / &root
            })
        });
        for value in matrix.iter_mut().flatten() {
            if value.rational_ref().is_none()
                && let Some(q) = dec.rational_value(value)?
            {
                *value = q.into();
            }
        }
        let mut map = Self {
            matrix,
            translation: Vector3::zero(),
            reversing: false,
        };
        let p = Vector3::new(
            point.x.clone().into(),
            point.y.clone().into(),
            point.z.clone().into(),
        );
        map.translation = p.minus(&map.apply_vector(&p));
        Ok(map)
    }

    /// Reflection in a rational plane. The normal need not have unit length.
    pub fn reflection(point: &RPoint3, normal: &RVector3) -> Result<Self, Refusal> {
        if normal.is_zero() {
            return Err(Refusal::new(
                "mirror.zero-normal",
                "mirror plane normal must not be zero",
            ));
        }
        let n = [&normal.x, &normal.y, &normal.z];
        let two = Rational::from_i64(2);
        let norm = normal.norm_squared();
        let matrix = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                Real::from(Rational::from_i64(i64::from(i == j)) - &(&(&two * n[i]) * n[j]) / &norm)
            })
        });
        let offset = normal.scaled(&((&two * &point.minus(&RPoint3::origin()).dot(normal)) / norm));
        Ok(Self {
            matrix,
            translation: offset.to_real(),
            reversing: true,
        })
    }

    /// Whether this map reverses orientation.
    pub fn reversing(&self) -> bool {
        self.reversing
    }

    /// Exact orthogonal matrix, acting on column vectors.
    pub fn matrix(&self) -> &[[Real; 3]; 3] {
        &self.matrix
    }

    /// Borrowed rational matrix entries, without allocating coefficient copies.
    pub fn rational_matrix(&self) -> Option<[[&Rational; 3]; 3]> {
        let row = |r: usize| {
            Some([
                self.matrix[r][0].rational_ref()?,
                self.matrix[r][1].rational_ref()?,
                self.matrix[r][2].rational_ref()?,
            ])
        };
        Some([row(0)?, row(1)?, row(2)?])
    }

    /// Image of a vector, without translation.
    pub fn apply_vector(&self, v: &Vector3) -> Vector3 {
        if let (Some(m), Some(v)) = (self.rational_matrix(), v.rational_components()) {
            let [x, y, z] = std::array::from_fn(|i| {
                (0..3).fold(Rational::zero(), |sum, j| sum + m[i][j] * v[j])
            });
            return Vector3::new(x.into(), y.into(), z.into());
        }
        let row = |i: usize| {
            &self.matrix[i][0] * &v.x + &self.matrix[i][1] * &v.y + &self.matrix[i][2] * &v.z
        };
        Vector3::new(row(0), row(1), row(2))
    }

    /// Image of a point, retaining exact expressions.
    pub fn apply_point(&self, p: &Point3) -> Point3 {
        if let (Some(m), Some(p), Some(t)) = (
            self.rational_matrix(),
            p.rational_components(),
            self.translation.rational_components(),
        ) {
            let [x, y, z] =
                std::array::from_fn(|i| (0..3).fold(t[i].clone(), |sum, j| sum + m[i][j] * p[j]));
            return Point3::new(x.into(), y.into(), z.into());
        }
        let v = self
            .apply_vector(&Vector3::new(p.x.clone(), p.y.clone(), p.z.clone()))
            .plus(&self.translation);
        Point3::new(v.x, v.y, v.z)
    }

    /// Rational image if every resulting coordinate has a native rational form.
    pub fn apply_rational_point(&self, p: &RPoint3) -> Option<RPoint3> {
        self.apply_point(&p.to_real()).as_rational()
    }

    /// Apply this map first, then `next`.
    pub fn then(&self, next: &Self) -> Self {
        if let (Some(a), Some(b)) = (self.rational_matrix(), next.rational_matrix()) {
            return Self {
                matrix: std::array::from_fn(|i| {
                    std::array::from_fn(|j| {
                        Real::from((0..3).fold(Rational::zero(), |sum, k| sum + b[i][k] * a[k][j]))
                    })
                }),
                translation: next.apply_vector(&self.translation).plus(&next.translation),
                reversing: self.reversing != next.reversing,
            };
        }
        Self {
            matrix: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..3).fold(Real::zero(), |sum, k| {
                        sum + &next.matrix[i][k] * &self.matrix[k][j]
                    })
                })
            }),
            translation: next.apply_vector(&self.translation).plus(&next.translation),
            reversing: self.reversing != next.reversing,
        }
    }

    /// Inverse obtained by transposition, not numeric matrix inversion.
    pub fn inverse(&self) -> Self {
        let mut inverse = Self {
            matrix: std::array::from_fn(|i| std::array::from_fn(|j| self.matrix[j][i].clone())),
            translation: Vector3::zero(),
            reversing: self.reversing,
        };
        inverse.translation = inverse
            .apply_vector(&self.translation)
            .scaled(&Real::from_i64(-1));
        inverse
    }

    /// Re-expresses provably rational coefficients without changing the map.
    /// Inconclusive reductions retain their exact expressions, never approximations.
    pub fn reduced(&self, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut result = self.clone();
        for value in result.matrix.iter_mut().flatten().chain([
            &mut result.translation.x,
            &mut result.translation.y,
            &mut result.translation.z,
        ]) {
            if value.rational_ref().is_none()
                && let Some(q) = dec.rational_value(value)?
            {
                *value = q.into();
            }
        }
        Ok(result)
    }

    /// Proves equality of maps under the caller's decision policy.
    pub fn same_motion(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        if let (Some(a), Some(b), Some(at), Some(bt)) = (
            self.rational_matrix(),
            other.rational_matrix(),
            self.translation.rational_components(),
            other.translation.rational_components(),
        ) {
            return Ok(a == b && at == bt);
        }
        for (a, b) in self
            .matrix
            .iter()
            .flatten()
            .zip(other.matrix.iter().flatten())
        {
            if crate::math::sign(&(a - b), dec)? != catallax_numbers::Sign::Zero {
                return Ok(false);
            }
        }
        let d = self.translation.minus(&other.translation);
        for v in [&d.x, &d.y, &d.z] {
            if crate::math::sign(v, dec)? != catallax_numbers::Sign::Zero {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Maps a supporting surface's definition. Dimensions are unchanged.
    pub fn surface(&self, surface: &Surface) -> Surface {
        let p = |p: &Point3| self.apply_point(p);
        let v = |v: &Vector3| self.apply_vector(v);
        match surface {
            Surface::Plane(s) => Surface::Plane(Plane {
                origin: p(&s.origin),
                normal: v(&s.normal),
                ref_dir: v(&s.ref_dir),
            }),
            Surface::Cylinder(s) => Surface::Cylinder(Cylinder {
                origin: p(&s.origin),
                axis: v(&s.axis),
                radius_sq: s.radius_sq.clone(),
                ref_dir: v(&s.ref_dir),
            }),
            Surface::Cone(s) => Surface::Cone(Cone {
                apex: p(&s.apex),
                axis: v(&s.axis),
                slope: s.slope.clone(),
                ref_dir: v(&s.ref_dir),
            }),
            Surface::Sphere(s) => Surface::Sphere(Sphere {
                center: p(&s.center),
                axis: v(&s.axis),
                radius: s.radius.clone(),
                ref_dir: v(&s.ref_dir),
            }),
            Surface::Torus(s) => Surface::Torus(Torus {
                center: p(&s.center),
                axis: v(&s.axis),
                major: s.major.clone(),
                minor: s.minor.clone(),
                ref_dir: v(&s.ref_dir),
            }),
            Surface::EllipticCylinder(s) => Surface::EllipticCylinder(EllipticCylinder {
                origin: p(&s.origin),
                axis: v(&s.axis),
                major_dir: v(&s.major_dir),
                major_sq: s.major_sq.clone(),
                minor_sq: s.minor_sq.clone(),
            }),
            Surface::BSpline(s) => Surface::BSpline(s.mapped(self)),
        }
    }

    /// Maps an explicit analytic or spline curve, preserving parameter direction.
    /// Intersection curves must instead remap their defining objects in the store.
    pub fn curve(&self, curve: &Curve) -> Option<Curve> {
        let p = |p: &Point3| self.apply_point(p);
        let v = |v: &Vector3| self.apply_vector(v);
        let axial = |n: &Vector3| {
            let n = v(n);
            if self.reversing {
                n.scaled(&Real::from_i64(-1))
            } else {
                n
            }
        };
        Some(match curve {
            Curve::Line(c) => Curve::Line(Line {
                origin: p(&c.origin),
                dir: v(&c.dir),
            }),
            Curve::Circle(c) => Curve::Circle(Circle {
                center: p(&c.center),
                normal: axial(&c.normal),
                radius_sq: c.radius_sq.clone(),
                ref_dir: v(&c.ref_dir),
            }),
            Curve::Ellipse(c) => Curve::Ellipse(Ellipse {
                center: p(&c.center),
                normal: axial(&c.normal),
                major_dir: v(&c.major_dir),
                major_sq: c.major_sq.clone(),
                minor_sq: c.minor_sq.clone(),
            }),
            Curve::BSpline(c) => Curve::BSpline(c.mapped(self)),
            Curve::Intersection(_) => return None,
        })
    }
}
