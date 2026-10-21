//! Exact rational rotations and rigid poses.

use catallax_numbers::{Rational, Refusal};

use crate::{Point3, RPoint3, RVector3, Vector3};

/// A rotation with an exact rational matrix, built from a rational quaternion.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Rotation {
    matrix: [[Rational; 3]; 3],
}

impl Rotation {
    /// The rotation of the non-zero quaternion `(w, x, y, z)`, with the matrix
    /// given in the piece description. Refuses the zero quaternion with
    /// `shapes.zero-quaternion`.
    pub fn from_quaternion(
        w: Rational,
        x: Rational,
        y: Rational,
        z: Rational,
    ) -> Result<Self, Refusal> {
        let (ww, xx, yy, zz) = (&w * &w, &x * &x, &y * &y, &z * &z);
        let n = (&ww + &xx) + (&yy + &zz);
        if n.is_zero() {
            return Err(Refusal::new(
                "shapes.zero-quaternion",
                "rotation quaternion is zero",
            ));
        }
        let two = Rational::from_i64(2);
        let (xy, xz, yz, wx, wy, wz) = (&x * &y, &x * &z, &y * &z, &w * &x, &w * &y, &w * &z);
        let diagonal =
            |a: &Rational, b: &Rational, c: &Rational, d: &Rational| &((a + b) - (c + d)) / &n;
        let off_diagonal = |v: Rational| &(&v * &two) / &n;
        Ok(Self {
            matrix: [
                [
                    diagonal(&ww, &xx, &yy, &zz),
                    off_diagonal(&xy - &wz),
                    off_diagonal(&xz + &wy),
                ],
                [
                    off_diagonal(&xy + &wz),
                    diagonal(&ww, &yy, &xx, &zz),
                    off_diagonal(&yz - &wx),
                ],
                [
                    off_diagonal(&xz - &wy),
                    off_diagonal(&yz + &wx),
                    diagonal(&ww, &zz, &xx, &yy),
                ],
            ],
        })
    }

    /// The identity.
    pub fn identity() -> Self {
        Self {
            matrix: std::array::from_fn(|i| {
                std::array::from_fn(|j| Rational::from_i64(i64::from(i == j)))
            }),
        }
    }

    /// The matrix, row by row; it acts on column vectors.
    pub fn matrix(&self) -> &[[Rational; 3]; 3] {
        &self.matrix
    }

    /// The rotation that applies `self` first, then `next`.
    pub fn then(&self, next: &Rotation) -> Rotation {
        Self {
            matrix: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..3).fold(Rational::zero(), |sum, k| {
                        sum + &next.matrix[i][k] * &self.matrix[k][j]
                    })
                })
            }),
        }
    }

    /// The inverse (the transpose).
    pub fn inverse(&self) -> Rotation {
        Self {
            matrix: std::array::from_fn(|i| std::array::from_fn(|j| self.matrix[j][i].clone())),
        }
    }

    /// `R v`.
    pub fn apply(&self, v: &RVector3) -> RVector3 {
        let row = |i: usize| {
            &self.matrix[i][0] * &v.x + &self.matrix[i][1] * &v.y + &self.matrix[i][2] * &v.z
        };
        RVector3::new(row(0), row(1), row(2))
    }

    /// `R v` for an exact real vector.
    pub fn apply_real(&self, v: &Vector3) -> Vector3 {
        let row = |i: usize| {
            &catallax_numbers::Real::from(self.matrix[i][0].clone()) * &v.x
                + &catallax_numbers::Real::from(self.matrix[i][1].clone()) * &v.y
                + &catallax_numbers::Real::from(self.matrix[i][2].clone()) * &v.z
        };
        Vector3::new(row(0), row(1), row(2))
    }
}

/// A rigid motion: rotate, then translate. `p' = R p + t`.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Pose {
    /// The rotation.
    pub rotation: Rotation,
    /// The translation applied after the rotation.
    pub translation: RVector3,
}

impl Pose {
    /// The identity pose.
    pub fn identity() -> Self {
        Self::new(Rotation::identity(), RVector3::zero())
    }

    /// The pose `p' = rotation p + translation`.
    pub fn new(rotation: Rotation, translation: RVector3) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    /// The pose that applies `self` first, then `next`.
    pub fn then(&self, next: &Pose) -> Pose {
        Self::new(
            self.rotation.then(&next.rotation),
            next.rotation
                .apply(&self.translation)
                .plus(&next.translation),
        )
    }

    /// The inverse pose.
    pub fn inverse(&self) -> Pose {
        let rotation = self.rotation.inverse();
        let translation = rotation
            .apply(&self.translation)
            .scaled(&Rational::from_i64(-1));
        Self::new(rotation, translation)
    }

    /// Image of a rational point.
    pub fn apply_point(&self, p: &RPoint3) -> RPoint3 {
        let v = self
            .rotation
            .apply(&RVector3::new(p.x.clone(), p.y.clone(), p.z.clone()))
            .plus(&self.translation);
        RPoint3::new(v.x, v.y, v.z)
    }

    /// Image of a rational direction (rotation only).
    pub fn apply_vector(&self, v: &RVector3) -> RVector3 {
        self.rotation.apply(v)
    }

    /// Image of an exact real point.
    pub fn apply_real_point(&self, p: &Point3) -> Point3 {
        let v = self
            .rotation
            .apply_real(&Vector3::new(p.x.clone(), p.y.clone(), p.z.clone()))
            .plus(&self.translation.to_real());
        Point3::new(v.x, v.y, v.z)
    }
}
