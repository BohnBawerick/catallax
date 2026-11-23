//! Exact grade laws shared by wall construction and field compilation.

use catallax_numbers::{Rational, Real};

/// Dot product in part coordinates.
pub fn rational_dot(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    &(&a[0] * &b[0]) + &(&(&a[1] * &b[1]) + &(&a[2] * &b[2]))
}

/// A clamped linear grade with distinct anchors.
#[derive(Clone, Debug)]
pub struct Linear {
    /// Where the ramp starts.
    pub from: [Rational; 3],
    /// The second anchor minus the first.
    pub direction: [Rational; 3],
    /// Squared length of the direction.
    pub length_sq: Rational,
    /// Value at and before the start.
    pub start: Rational,
    /// Value at and after the end.
    pub end: Rational,
}

impl Linear {
    /// Builds a ramp from distinct anchors.
    pub fn new(from: [Rational; 3], to: [Rational; 3], start: Rational, end: Rational) -> Self {
        let direction = std::array::from_fn(|i| &to[i] - &from[i]);
        Self {
            length_sq: rational_dot(&direction, &direction),
            from,
            direction,
            start,
            end,
        }
    }

    /// The unclamped fraction of the ramp at a point.
    pub fn parameter(&self, point: &[Rational; 3]) -> Rational {
        let rel = std::array::from_fn(|i| &point[i] - &self.from[i]);
        &rational_dot(&rel, &self.direction) / &self.length_sq
    }

    /// The value at a ramp parameter, including the two clamps.
    pub fn interpolate(&self, t: Rational) -> Rational {
        let t = t.max(Rational::zero()).min(Rational::one());
        &self.start + &(&t * &(&self.end - &self.start))
    }
}

/// Projects offsets perpendicular to an axis, or leaves them unchanged for a point grade.
pub fn distance_projection(direction: Option<&[Rational; 3]>) -> [[Rational; 3]; 3] {
    let mut matrix =
        std::array::from_fn(|i| std::array::from_fn(|j| Rational::from_i64(i64::from(i == j))));
    if let Some(axis) = direction {
        let squared = rational_dot(axis, axis);
        for i in 0..3 {
            for j in 0..3 {
                matrix[i][j] = &matrix[i][j] - &(&(&axis[i] * &axis[j]) / &squared);
            }
        }
    }
    matrix
}

/// A linear ramp in distance, clamped before and after two ordered radii.
#[derive(Clone, Debug)]
pub struct RadiusRamp {
    /// Distinct, nonnegative radii.
    pub radii: [Rational; 2],
    /// Squared radii used to select the clamps without square roots.
    pub squared: [Rational; 2],
    /// Start and end values.
    pub values: [Rational; 2],
    /// Change in value per unit distance.
    pub slope: Rational,
    /// Value of the extended ramp at zero distance.
    pub intercept: Rational,
}

impl RadiusRamp {
    /// Builds a ramp with validated radii.
    pub fn new(radii: [Rational; 2], values: [Rational; 2]) -> Self {
        let slope = &(&values[1] - &values[0]) / &(&radii[1] - &radii[0]);
        let intercept = &values[0] - &(&radii[0] * &slope);
        Self {
            squared: radii.each_ref().map(|r| r * r),
            radii,
            values,
            slope,
            intercept,
        }
    }

    /// The endpoint value if the squared distance lies in a clamp.
    pub fn clamp(&self, squared: &Rational) -> Option<&Rational> {
        if squared <= &self.squared[0] {
            Some(&self.values[0])
        } else if squared >= &self.squared[1] {
            Some(&self.values[1])
        } else {
            None
        }
    }

    /// The exact value at a squared distance.
    pub fn real(&self, squared: Rational) -> Real {
        if let Some(value) = self.clamp(&squared) {
            return Real::from(value.clone());
        }
        &Real::from(self.intercept.clone())
            + &(&Real::from(squared).sqrt() * &Real::from(self.slope.clone()))
    }
}
