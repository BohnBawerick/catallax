use super::{Piece, clamp01, f, q, rational_dot};
use crate::{
    geometry,
    metric::{PRECISION, Scalar},
};
use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal, Sign};
use catallax_shapes::Vector3;

/// Distance to the centre or axis is the norm of the projected offset.
#[derive(Clone, Debug)]
pub(crate) struct Radius {
    origin: [Rational; 3],
    direction: Option<[Rational; 3]>,
    projection: [[Rational; 3]; 3],
    interval_projection: [[Interval; 3]; 3],
    fast_projection: [[f64; 3]; 3],
    radii: [Rational; 2],
    squared: [Rational; 2],
    interval_radii: [Interval; 2],
    pub(super) values: [Rational; 2],
    pub(super) slope: Rational,
    intercept: Rational,
}

fn projection(direction: Option<&[Rational; 3]>) -> [[Rational; 3]; 3] {
    let mut matrix = std::array::from_fn(|i| std::array::from_fn(|j| q(i64::from(i == j))));
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

impl Radius {
    pub(super) fn new(
        origin: [Rational; 3],
        direction: Option<[Rational; 3]>,
        radii: [Rational; 2],
        values: [Rational; 2],
    ) -> Self {
        let projection = projection(direction.as_ref());
        let slope = &(&values[1] - &values[0]) / &(&radii[1] - &radii[0]);
        let intercept = &values[0] - &(&radii[0] * &slope);
        Self {
            origin,
            direction,
            interval_projection: projection
                .each_ref()
                .map(|row| row.each_ref().map(Interval::from_rational)),
            fast_projection: projection.each_ref().map(|row| row.each_ref().map(f)),
            projection,
            squared: radii.each_ref().map(|r| r * r),
            interval_radii: radii.each_ref().map(Interval::from_rational),
            radii,
            values,
            slope,
            intercept,
        }
    }

    fn offset(&self, p: &[Rational; 3]) -> [Rational; 3] {
        let delta = std::array::from_fn(|i| &p[i] - &self.origin[i]);
        self.projection
            .each_ref()
            .map(|row| rational_dot(row, &delta))
    }

    fn squared_distance(&self, p: &[Rational; 3]) -> Rational {
        let offset = self.offset(p);
        rational_dot(&offset, &offset)
    }

    fn clamp(&self, squared: &Rational) -> Option<&Rational> {
        if squared <= &self.squared[0] {
            Some(&self.values[0])
        } else if squared >= &self.squared[1] {
            Some(&self.values[1])
        } else {
            None
        }
    }

    pub(super) fn real(&self, p: &[Rational; 3]) -> Real {
        let squared = self.squared_distance(p);
        if let Some(value) = self.clamp(&squared) {
            return Real::from(value.clone());
        }
        &Real::from(self.intercept.clone())
            + &(&Real::from(squared).sqrt() * &Real::from(self.slope.clone()))
    }

    pub(super) fn ball(&self, p: &[Rational; 3]) -> Ball {
        let squared = self.squared_distance(p);
        if let Some(value) = self.clamp(&squared) {
            return Ball::from_rational(value, PRECISION);
        }
        let distance = Ball::from_rational(&squared, PRECISION)
            .sqrt()
            .expect("nonnegative squared distance");
        &Ball::from_rational(&self.intercept, PRECISION)
            + &(&distance * &Ball::from_rational(&self.slope, PRECISION))
    }

    pub(super) fn fast(&self, p: [f64; 3]) -> f64 {
        let delta = std::array::from_fn::<_, 3, _>(|i| p[i] - f(&self.origin[i]));
        let squared = self
            .fast_projection
            .iter()
            .map(|row| {
                let component = (0..3).fold(0.0, |sum, i| sum + row[i] * delta[i]);
                component * component
            })
            .sum::<f64>();
        let distance = Scalar::sqrt(&squared);
        let [r0, r1] = self.radii.each_ref().map(f);
        let [start, end] = self.values.each_ref().map(f);
        if distance <= r0 {
            start
        } else if distance >= r1 {
            end
        } else {
            f(&self.intercept) + distance * f(&self.slope)
        }
    }

    fn interval_offset(&self, p: &[Interval; 3]) -> [Interval; 3] {
        let delta =
            std::array::from_fn::<_, 3, _>(|i| p[i] - Interval::from_rational(&self.origin[i]));
        self.interval_projection.map(|row| {
            (0..3).fold(Interval::point(0.0).unwrap(), |sum, i| {
                sum + row[i] * delta[i]
            })
        })
    }

    fn norm(offset: [Interval; 3]) -> Interval {
        (offset[0].square() + offset[1].square() + offset[2].square())
            .sqrt()
            .expect("nonnegative squared distance")
    }

    pub(super) fn interval(&self, p: &[Interval; 3]) -> Interval {
        let distance = Self::norm(self.interval_offset(p));
        let t = clamp01(
            (distance - self.interval_radii[0])
                / Interval::from_rational(&(&self.radii[1] - &self.radii[0])),
        );
        Interval::from_rational(&self.values[0])
            + t * Interval::from_rational(&(&self.values[1] - &self.values[0]))
    }

    pub(super) fn piece(&self, p: &[Interval; 3]) -> Option<Piece> {
        let distance = Self::norm(self.interval_offset(p));
        let [r0, r1] = self.interval_radii;
        if distance.hi() <= r0.lo() {
            Some(Piece::Start)
        } else if distance.lo() >= r1.hi() {
            Some(Piece::End)
        } else if distance.lo() > 0.0 && distance.lo() >= r0.hi() && distance.hi() <= r1.lo() {
            Some(Piece::Curved)
        } else {
            None
        }
    }

    pub(super) fn derivatives(
        &self,
        p: &[Interval; 3],
    ) -> Option<([Interval; 3], [[Interval; 3]; 3])> {
        let offset = self.interval_offset(p);
        let distance = Self::norm(offset);
        if distance.lo() <= 0.0 || !distance.hi().is_finite() {
            return None;
        }
        let slope = Interval::from_rational(&self.slope);
        let unit = offset.map(|v| v / distance);
        // For u = Q(p - origin) and r = |u|, Q = Q^T = Q².
        // The gradient is u/r and the Hessian is (Q - (u/r)(u/r)^T)/r.
        let hessian = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                slope * (self.interval_projection[i][j] - unit[i] * unit[j]) / distance
            })
        });
        Some((unit.map(|v| v * slope), hessian))
    }

    pub(super) fn rate(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let offset = self.offset(p);
        let squared = rational_dot(&offset, &offset);
        if squared < self.squared[0] || squared > self.squared[1] {
            return Ok(None);
        }
        let direction = [&d.x, &d.y, &d.z];
        let dot = |a: &[Rational; 3]| {
            (0..3).fold(Real::zero(), |sum, i| {
                &sum + &(&Real::from(a[i].clone()) * direction[i])
            })
        };
        let rate = if squared.is_zero() {
            self.projection
                .iter()
                .fold(Real::zero(), |sum, row| &sum + &dot(row).square())
                .sqrt()
        } else {
            let along = dot(&offset);
            if (squared == self.squared[0] && geometry::sign(&along, dec)? != Sign::Positive)
                || (squared == self.squared[1] && geometry::sign(&along, dec)? != Sign::Negative)
            {
                return Ok(None);
            }
            &along / &Real::from(squared).sqrt()
        };
        Ok(Some(&rate * &Real::from(self.slope.clone())))
    }

    pub(super) fn transform(&mut self, matrix: &[[Rational; 3]; 3], translation: &[Rational; 3]) {
        let apply = |v: &[Rational; 3]| matrix.each_ref().map(|row| rational_dot(row, v));
        let moved = apply(&self.origin);
        self.origin = std::array::from_fn(|i| &moved[i] + &translation[i]);
        if let Some(direction) = &mut self.direction {
            *direction = apply(direction);
            self.projection = projection(Some(direction));
            self.interval_projection = self
                .projection
                .each_ref()
                .map(|row| row.each_ref().map(Interval::from_rational));
            self.fast_projection = self.projection.each_ref().map(|row| row.each_ref().map(f));
        }
    }
}
