//! Exact finite tensor geometry for STEP extrusions of rational spline curves.

use super::*;
use crate::exact::dot;

fn include(range: &mut Option<[Rational; 2]>, value: Rational) {
    *range = Some(match range.take() {
        Some([lo, hi]) => [lo.min(value.clone()), hi.max(value)],
        None => [value.clone(), value],
    });
}

fn projection_range<'a>(
    points: impl IntoIterator<Item = &'a Vec3>,
    direction: &Vec3,
    dec: &mut Decider,
) -> Result<Option<[Rational; 2]>, ImportRefusal> {
    let mut range = None;
    for point in points {
        dec.charge_work(1).map_err(spline::source_error)?;
        include(&mut range, dot(point, direction));
    }
    Ok(range)
}

fn curve_range(
    edge: &OrientedEdge,
    direction: &Vec3,
    dec: &mut Decider,
) -> Result<[Rational; 2], ImportRefusal> {
    dec.charge_work(1).map_err(spline::source_error)?;
    Ok(match &edge.curve {
        Curve::Line(..) => projection_range([&edge.start.point, &edge.end.point], direction, dec)?
            .expect("two endpoints"),
        Curve::BSpline(curve) => projection_range(curve.geometry.controls(), direction, dec)?
            .expect("validated source controls"),
        Curve::Circle(placement, radius) | Curve::Ellipse(placement, radius, _) => {
            let radius = if let Curve::Ellipse(_, a, b) = &edge.curve {
                a.max(b)
            } else {
                radius
            };
            // The L1 norm bounds the Euclidean support radius in every frame.
            dec.charge_work(3).map_err(spline::source_error)?;
            let norm = direction
                .iter()
                .fold(Rational::zero(), |sum, component| sum + component.abs());
            let spread = radius * &norm;
            let center = dot(&placement.origin, direction);
            [&center - &spread, center + spread]
        }
    })
}

impl Reader<'_> {
    pub(super) fn linear_extrusion(
        &self,
        id: u64,
        loops: &[FaceLoop],
    ) -> Result<Surface, ImportRefusal> {
        let args = self.args(id, "SURFACE_OF_LINEAR_EXTRUSION")?;
        let profile = self.curve(reference(args, 1)?)?;
        let Curve::BSpline(profile) = profile else {
            return Err(unsupported(id, self.entity(id)?));
        };
        let axis = self.args(reference(args, 2)?, "VECTOR")?;
        let mut direction = self.direction(reference(axis, 1)?)?;
        let magnitude = self.length(axis, 2)?;
        match magnitude.sign() {
            catallax_numbers::Sign::Zero => {
                return Err(malformed("extrusion vector must not be zero"));
            }
            catallax_numbers::Sign::Negative => direction = direction.map(|v| -v),
            catallax_numbers::Sign::Positive => {}
        }
        let closed = profile.closed;
        let profile = profile.geometry;
        let mut meter = Decider::new(*self.dec.policy());
        meter.set_deadline(self.dec.deadline());
        let [profile_lo, profile_hi] =
            projection_range(profile.controls(), &direction, &mut meter)?
                .expect("validated source controls");
        let mut boundary = None;
        for bound in loops {
            meter.charge_work(1).map_err(spline::source_error)?;
            for edge in &bound.edges {
                for limit in curve_range(edge, &direction, &mut meter)? {
                    include(&mut boundary, limit);
                }
            }
            if let Some(vertex) = &bound.vertex {
                meter.charge_work(1).map_err(spline::source_error)?;
                include(&mut boundary, dot(&vertex.point, &direction));
            }
        }
        let [boundary_lo, boundary_hi] = boundary.ok_or_else(|| ImportRefusal::NotSolid {
            detail: format!("extrusion #{id} has no finite face boundary"),
        })?;
        // Direction ratios rescale v positively without changing the carrier or its normal.
        let norm_sq = dot(&direction, &direction);
        let lo = &(boundary_lo - profile_hi) / &norm_sq;
        let hi = &(boundary_hi - profile_lo) / &norm_sq;
        if lo >= hi {
            return Err(malformed(format!(
                "extrusion #{id} has no finite surface interval"
            )));
        }
        // Positive weights bound the whole profile and each whole trim curve.
        // Projection along the sweep is monotone, so its finite face lies in this strip.
        let mut controls = Vec::new();
        let mut weights = Vec::new();
        for (point, weight) in profile.controls().iter().zip(profile.weights()) {
            for t in [&lo, &hi] {
                meter.charge_work(1).map_err(spline::source_error)?;
                controls.push(std::array::from_fn(|i| &point[i] + &(t * &direction[i])));
                weights.push(weight.clone());
            }
        }
        let mut knots = Vec::new();
        for knot in profile.knots() {
            meter.charge_work(1).map_err(spline::source_error)?;
            knots.push(knot.clone());
        }
        let geometry = RationalSplineSurface::new(
            [profile.degree(), 1],
            [
                knots,
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    Rational::one(),
                    Rational::one(),
                ],
            ],
            [profile.controls().len(), 2],
            controls,
            weights,
            self.dec,
        )
        .map_err(spline::source_error)?;
        Ok(Surface::BSpline(StepSplineSurface {
            geometry,
            closed: [closed, Some(false)],
            self_intersect: None,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;
    use std::time::{Duration, Instant};

    #[test]
    fn projection_stops_on_the_shared_allowance_before_visiting_the_whole_source() {
        let point = [Rational::one(), Rational::zero(), Rational::zero()];
        let visited = Cell::new(0);
        let points = (0..100_000).map(|_| {
            visited.set(visited.get() + 1);
            &point
        });
        let (result, work) = catallax_numbers::with_operation_budget(1, || {
            projection_range(points, &point, &mut Decider::default())
        });
        assert!(matches!(result, Err(ImportRefusal::ResourceLimit)));
        assert!(work.exhausted);
        assert!(visited.get() <= 2);
    }

    #[test]
    fn projection_obeys_the_callers_deadline() {
        let point = [Rational::one(), Rational::zero(), Rational::zero()];
        let mut dec = Decider::default();
        dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
        assert!(matches!(
            projection_range([&point], &point, &mut dec),
            Err(ImportRefusal::ResourceLimit)
        ));
    }
}
