use super::{Piece, clamp01, f};
use crate::{
    geometry,
    metric::{PRECISION, Scalar},
};
use catallax_features::grade::{RadiusRamp, distance_projection, rational_dot};
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
    interval_radii: [Interval; 2],
    pub(super) ramp: RadiusRamp,
}

impl Radius {
    pub(super) fn new(
        origin: [Rational; 3],
        direction: Option<[Rational; 3]>,
        radii: [Rational; 2],
        values: [Rational; 2],
    ) -> Self {
        let projection = distance_projection(direction.as_ref());
        Self {
            origin,
            direction,
            interval_projection: projection
                .each_ref()
                .map(|row| row.each_ref().map(Interval::from_rational)),
            fast_projection: projection.each_ref().map(|row| row.each_ref().map(f)),
            projection,
            interval_radii: radii.each_ref().map(Interval::from_rational),
            ramp: RadiusRamp::new(radii, values),
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

    pub(super) fn real(&self, p: &[Rational; 3]) -> Real {
        self.ramp.real(self.squared_distance(p))
    }

    pub(super) fn ball(&self, p: &[Rational; 3]) -> Ball {
        let squared = self.squared_distance(p);
        if let Some(value) = self.ramp.clamp(&squared) {
            return Ball::from_rational(value, PRECISION);
        }
        let distance = Ball::from_rational(&squared, PRECISION)
            .sqrt()
            .expect("nonnegative squared distance");
        &Ball::from_rational(&self.ramp.intercept, PRECISION)
            + &(&distance * &Ball::from_rational(&self.ramp.slope, PRECISION))
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
        let [r0, r1] = self.ramp.radii.each_ref().map(f);
        let [start, end] = self.ramp.values.each_ref().map(f);
        if distance <= r0 {
            start
        } else if distance >= r1 {
            end
        } else {
            f(&self.ramp.intercept) + distance * f(&self.ramp.slope)
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
                / Interval::from_rational(&(&self.ramp.radii[1] - &self.ramp.radii[0])),
        );
        Interval::from_rational(&self.ramp.values[0])
            + t * Interval::from_rational(&(&self.ramp.values[1] - &self.ramp.values[0]))
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
        let slope = Interval::from_rational(&self.ramp.slope);
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
        if squared < self.ramp.squared[0] || squared > self.ramp.squared[1] {
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
            if (squared == self.ramp.squared[0] && geometry::sign(&along, dec)? != Sign::Positive)
                || (squared == self.ramp.squared[1]
                    && geometry::sign(&along, dec)? != Sign::Negative)
            {
                return Ok(None);
            }
            &along / &Real::from(squared).sqrt()
        };
        Ok(Some(&rate * &Real::from(self.ramp.slope.clone())))
    }

    pub(super) fn transform(&mut self, matrix: &[[Rational; 3]; 3], translation: &[Rational; 3]) {
        let apply = |v: &[Rational; 3]| matrix.each_ref().map(|row| rational_dot(row, v));
        let moved = apply(&self.origin);
        self.origin = std::array::from_fn(|i| &moved[i] + &translation[i]);
        if let Some(direction) = &mut self.direction {
            *direction = apply(direction);
            self.projection = distance_projection(Some(direction));
            self.interval_projection = self
                .projection
                .each_ref()
                .map(|row| row.each_ref().map(Interval::from_rational));
            self.fast_projection = self.projection.each_ref().map(|row| row.each_ref().map(f));
        }
    }
}
