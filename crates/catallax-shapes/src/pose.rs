//! Exact algebraic rotations and rigid poses.

use crate::{Point3, RPoint3, RVector3, Vector3, math};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};

/// A certified proper orthogonal rotation, retaining its exact matrix definition.
#[derive(Clone, Debug)]
pub struct Rotation {
    matrix: [[Real; 3]; 3],
}

impl Rotation {
    /// Checked rational-input quaternion convenience constructor.
    pub fn from_quaternion(
        w: Rational,
        x: Rational,
        y: Rational,
        z: Rational,
    ) -> Result<Self, Refusal> {
        Self::rational_quaternion(&w, &x, &y, &z)
    }

    fn rational_quaternion(
        w: &Rational,
        x: &Rational,
        y: &Rational,
        z: &Rational,
    ) -> Result<Self, Refusal> {
        let (ww, xx, yy, zz) = (w * w, x * x, y * y, z * z);
        let n = (&ww + &xx) + (&yy + &zz);
        if n.is_zero() {
            return Err(Refusal::new(
                "shapes.zero-quaternion",
                "rotation quaternion is zero",
            ));
        }
        let two = Rational::from_i64(2);
        let (xy, xz, yz, wx, wy, wz) = (x * y, x * z, y * z, w * x, w * y, w * z);
        let diagonal =
            |a: &Rational, b: &Rational, c: &Rational, d: &Rational| &((a + b) - (c + d)) / &n;
        let off = |v: Rational| &(&v * &two) / &n;
        Ok(Self::from_rational_matrix([
            [diagonal(&ww, &xx, &yy, &zz), off(&xy - &wz), off(&xz + &wy)],
            [off(&xy + &wz), diagonal(&ww, &yy, &xx, &zz), off(&yz - &wx)],
            [off(&xz - &wy), off(&yz + &wx), diagonal(&ww, &zz, &xx, &yy)],
        ]))
    }

    /// Validates algebraic quaternion domains and proves its norm nonzero.
    pub fn from_real_quaternion(
        w: Real,
        x: Real,
        y: Real,
        z: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        for value in [&w, &x, &y, &z] {
            math::validate(value, dec)?;
        }
        if let (Some(w), Some(x), Some(y), Some(z)) = (
            w.rational_ref(),
            x.rational_ref(),
            y.rational_ref(),
            z.rational_ref(),
        ) {
            return Self::rational_quaternion(w, x, y, z);
        }
        let (ww, xx, yy, zz) = (w.square(), x.square(), y.square(), z.square());
        let n = &(&ww + &xx) + &(&yy + &zz);
        if math::sign(&n, dec)? == Sign::Zero {
            return Err(Refusal::new(
                "shapes.zero-quaternion",
                "rotation quaternion is zero",
            ));
        }
        let two = Real::from_i64(2);
        let (xy, xz, yz, wx, wy, wz) = (&x * &y, &x * &z, &y * &z, &w * &x, &w * &y, &w * &z);
        let diagonal = |a: &Real, b: &Real, c: &Real, d: &Real| &((a + b) - (c + d)) / &n;
        let off = |v: Real| &(&v * &two) / &n;
        Ok(Self {
            matrix: [
                [diagonal(&ww, &xx, &yy, &zz), off(&xy - &wz), off(&xz + &wy)],
                [off(&xy + &wz), diagonal(&ww, &yy, &xx, &zz), off(&yz - &wx)],
                [off(&xz - &wy), off(&yz + &wx), diagonal(&ww, &zz, &xx, &yy)],
            ],
        })
    }

    fn from_rational_matrix(matrix: [[Rational; 3]; 3]) -> Self {
        Self {
            matrix: matrix.map(|row| row.map(Real::from)),
        }
    }

    /// The identity rotation.
    pub fn identity() -> Self {
        let zero = Real::zero();
        let one = Real::one();
        Self {
            matrix: std::array::from_fn(|i| {
                std::array::from_fn(|j| if i == j { one.clone() } else { zero.clone() })
            }),
        }
    }

    /// Exact matrix, row by row, acting on column vectors.
    pub fn matrix(&self) -> &[[Real; 3]; 3] {
        &self.matrix
    }

    /// Borrowed rational entries, without duplicating matrix storage or coefficients.
    pub fn rational_matrix(&self) -> Option<[[&Rational; 3]; 3]> {
        fn row(values: &[Real; 3]) -> Option<[&Rational; 3]> {
            Some([
                values[0].rational_ref()?,
                values[1].rational_ref()?,
                values[2].rational_ref()?,
            ])
        }
        Some([
            row(&self.matrix[0])?,
            row(&self.matrix[1])?,
            row(&self.matrix[2])?,
        ])
    }

    /// Applies `self` first, then `next`, preserving proper orthogonality by construction.
    pub fn then(&self, next: &Rotation) -> Rotation {
        if let (Some(a), Some(b)) = (self.rational_matrix(), next.rational_matrix()) {
            return Self::from_rational_matrix(std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..3).fold(Rational::zero(), |sum, k| sum + b[i][k] * a[k][j])
                })
            }));
        }
        Self {
            matrix: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..3).fold(Real::zero(), |sum, k| {
                        sum + &next.matrix[i][k] * &self.matrix[k][j]
                    })
                })
            }),
        }
    }

    /// Proves equality of rotation matrices under the caller's decision policy.
    pub fn same_rotation(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        if let (Some(a), Some(b)) = (self.rational_matrix(), other.rational_matrix()) {
            return Ok(a == b);
        }
        for (a, b) in self
            .matrix
            .iter()
            .flatten()
            .zip(other.matrix.iter().flatten())
        {
            if math::sign(&(a - b), dec)? != Sign::Zero {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// The inverse, obtained by transposition without a new proof.
    pub fn inverse(&self) -> Rotation {
        Self {
            matrix: std::array::from_fn(|i| std::array::from_fn(|j| self.matrix[j][i].clone())),
        }
    }

    /// Image of an exact vector.
    pub fn apply(&self, vector: &Vector3) -> Vector3 {
        if let (Some(matrix), Some(v)) = (self.rational_matrix(), vector.rational_components()) {
            let [x, y, z] = rational_apply(&matrix, &v);
            return Vector3::new(x.into(), y.into(), z.into());
        }
        let row = |i: usize| {
            &self.matrix[i][0] * &vector.x
                + &self.matrix[i][1] * &vector.y
                + &self.matrix[i][2] * &vector.z
        };
        Vector3::new(row(0), row(1), row(2))
    }

    /// Image of an exact vector.
    pub fn apply_real(&self, vector: &Vector3) -> Vector3 {
        self.apply(vector)
    }
}

fn rational_apply(matrix: &[[&Rational; 3]; 3], vector: &[&Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| (0..3).fold(Rational::zero(), |sum, j| sum + matrix[i][j] * vector[j]))
}

/// A validated rigid motion, with an exact algebraic translation.
#[derive(Clone, Debug)]
pub struct Pose {
    /// Proper orthogonal rotation.
    pub rotation: Rotation,
    /// Translation applied after the rotation.
    pub translation: Vector3,
}

impl Pose {
    /// The identity pose.
    pub fn identity() -> Self {
        Self {
            rotation: Rotation::identity(),
            translation: Vector3::zero(),
        }
    }

    /// A validated rotation and an exact rational translation need no numeric proof.
    pub fn from_rational(rotation: Rotation, translation: RVector3) -> Self {
        Self {
            rotation,
            translation: translation.into(),
        }
    }

    /// Validates translation domains before publishing the pose.
    pub fn new(
        rotation: Rotation,
        translation: Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_vector(&translation, dec)?;
        Ok(Self {
            rotation,
            translation,
        })
    }

    /// Applies `self` first, then `next`, preserving validated domains.
    pub fn then(&self, next: &Pose) -> Pose {
        if let (Some(matrix), Some(a), Some(b)) = (
            next.rotation.rational_matrix(),
            self.translation.rational_components(),
            next.translation.rational_components(),
        ) {
            let [x, y, z] = rational_apply(&matrix, &a);
            return Self {
                rotation: self.rotation.then(&next.rotation),
                translation: Vector3::new(
                    (&x + b[0]).into(),
                    (&y + b[1]).into(),
                    (&z + b[2]).into(),
                ),
            };
        }
        Self {
            rotation: self.rotation.then(&next.rotation),
            translation: next
                .rotation
                .apply(&self.translation)
                .plus(&next.translation),
        }
    }

    /// Proves equality of rigid motions, without comparing expression keys.
    pub fn same_motion(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        if !self.rotation.same_rotation(&other.rotation, dec)? {
            return Ok(false);
        }
        if let (Some(a), Some(b)) = (
            self.translation.rational_components(),
            other.translation.rational_components(),
        ) {
            return Ok(a == b);
        }
        let difference = self.translation.minus(&other.translation);
        for component in [&difference.x, &difference.y, &difference.z] {
            if math::sign(component, dec)? != Sign::Zero {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Inverse of a validated pose.
    pub fn inverse(&self) -> Pose {
        let rotation = self.rotation.inverse();
        if let (Some(matrix), Some(translation)) = (
            rotation.rational_matrix(),
            self.translation.rational_components(),
        ) {
            let [x, y, z] = rational_apply(&matrix, &translation);
            return Self {
                rotation,
                translation: Vector3::new((-x).into(), (-y).into(), (-z).into()),
            };
        }
        let translation = rotation
            .apply(&self.translation)
            .scaled(&Real::from_i64(-1));
        Self {
            rotation,
            translation,
        }
    }

    /// Image of an exact point.
    pub fn apply_point(&self, point: &Point3) -> Point3 {
        if let (Some(matrix), Some(point), Some(translation)) = (
            self.rotation.rational_matrix(),
            point.rational_components(),
            self.translation.rational_components(),
        ) {
            let [x, y, z] = rational_apply(&matrix, &point);
            return Point3::new(
                (&x + translation[0]).into(),
                (&y + translation[1]).into(),
                (&z + translation[2]).into(),
            );
        }
        let v = self
            .rotation
            .apply(&Vector3::new(
                point.x.clone(),
                point.y.clone(),
                point.z.clone(),
            ))
            .plus(&self.translation);
        Point3::new(v.x, v.y, v.z)
    }

    /// Image of an exact direction, without translation.
    pub fn apply_vector(&self, vector: &Vector3) -> Vector3 {
        self.rotation.apply(vector)
    }

    /// Native rational image, when the stored frame and translation are rational.
    /// `None` does not round or attempt an algebraic-to-rational proof.
    pub fn apply_rational_point(&self, point: &RPoint3) -> Option<RPoint3> {
        let matrix = self.rotation.rational_matrix()?;
        let translation = self.translation.rational_components()?;
        let [x, y, z] = rational_apply(&matrix, &[&point.x, &point.y, &point.z]);
        Some(RPoint3::new(
            &x + translation[0],
            &y + translation[1],
            &z + translation[2],
        ))
    }

    /// Native rational direction image, without constructing temporary `Real` leaves.
    pub fn apply_rational_vector(&self, vector: &RVector3) -> Option<RVector3> {
        let matrix = self.rotation.rational_matrix()?;
        let [x, y, z] = rational_apply(&matrix, &[&vector.x, &vector.y, &vector.z]);
        Some(RVector3::new(x, y, z))
    }

    /// Image of an exact point.
    pub fn apply_real_point(&self, point: &Point3) -> Point3 {
        self.apply_point(point)
    }
}
