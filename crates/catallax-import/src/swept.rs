//! Exact tensor geometry for STEP swept rational spline curves.

use super::*;
use crate::exact::{cross, dot, sub};

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

// Four rational quadratic arcs sweep each profile control through one turn.
// Axis ratios keep the resulting controls and weights rational.
fn revolved(
    profile: StepSplineCurve,
    origin: &Vec3,
    axis: &Vec3,
    dec: &Decider,
) -> Result<StepSplineSurface, ImportRefusal> {
    let mut meter = Decider::new(*dec.policy());
    meter.set_deadline(dec.deadline());
    meter.charge_work(1).map_err(spline::source_error)?;
    let s = dot(axis, axis);
    if s.is_zero() {
        return Err(malformed("revolution axis must not be zero"));
    }
    let one = Rational::one();
    let two = Rational::from_i64(2);
    let next = &one + &s;
    let rotation_weights = [&one, &one, &next, &s, &s, &s, &next, &one, &one];
    let mut controls: [Vec<Vec3>; 9] = std::array::from_fn(|_| Vec::new());
    let mut weights: [Vec<Rational>; 9] = std::array::from_fn(|_| Vec::new());
    for (point, weight) in profile
        .geometry
        .controls()
        .iter()
        .zip(profile.geometry.weights())
    {
        meter.charge_work(1).map_err(spline::source_error)?;
        let k = cross(axis, &sub(point, origin));
        let l = cross(axis, &k);
        let negative_k = k.each_ref().map(|value| -value.clone());
        for row in 0..9 {
            meter.charge_work(1).map_err(spline::source_error)?;
            controls[row].push(std::array::from_fn(|i| {
                let offset = match row {
                    0 | 8 => Rational::zero(),
                    1 => k[i].clone(),
                    2 => &(&two * &(&k[i] + &l[i])) / &next,
                    3 => &(&k[i] + &(&two * &l[i])) / &s,
                    4 => &(&two * &l[i]) / &s,
                    5 => &(&negative_k[i] + &(&two * &l[i])) / &s,
                    6 => &(&two * &(&negative_k[i] + &l[i])) / &next,
                    7 => negative_k[i].clone(),
                    _ => unreachable!(),
                };
                &point[i] + &offset
            }));
            weights[row].push(weight * rotation_weights[row]);
        }
    }
    let mut u_knots = Vec::new();
    for value in [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4] {
        meter.charge_work(1).map_err(spline::source_error)?;
        u_knots.push(Rational::from_i64(value));
    }
    let mut v_knots = Vec::new();
    for knot in profile.geometry.knots() {
        meter.charge_work(1).map_err(spline::source_error)?;
        v_knots.push(knot.clone());
    }
    let geometry = RationalSplineSurface::new(
        [2, profile.geometry.degree()],
        [u_knots, v_knots],
        [9, profile.geometry.controls().len()],
        controls.into_iter().flatten().collect(),
        weights.into_iter().flatten().collect(),
        dec,
    )
    .map_err(spline::source_error)?;
    Ok(StepSplineSurface {
        geometry,
        closed: [Some(true), profile.closed],
        self_intersect: None,
    })
}

impl Reader<'_> {
    pub(super) fn revolution(&self, id: u64) -> Result<Surface, ImportRefusal> {
        let args = self.args(id, "SURFACE_OF_REVOLUTION")?;
        if args.len() != 3 {
            return Err(malformed("SURFACE_OF_REVOLUTION requires three fields"));
        }
        let Curve::BSpline(profile) = self.curve(reference(args, 1)?)? else {
            return Err(unsupported(id, self.entity(id)?));
        };
        let placement = self.args(reference(args, 2)?, "AXIS1_PLACEMENT")?;
        if placement.len() != 3 {
            return Err(malformed("AXIS1_PLACEMENT requires three fields"));
        }
        let origin = self.point(reference(placement, 1)?)?;
        let axis = match optional_ref(placement, 2)? {
            Some(id) => self.direction(id)?,
            None => [Rational::zero(), Rational::zero(), Rational::one()],
        };
        Ok(Surface::BSpline(revolved(
            profile, &origin, &axis, self.dec,
        )?))
    }

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

    fn meridian() -> StepSplineCurve {
        let zero = Rational::zero();
        let one = Rational::one();
        StepSplineCurve {
            geometry: RationalSplineCurve::new(
                1,
                vec![zero.clone(), zero.clone(), one.clone(), one.clone()],
                vec![
                    [one.clone(), zero.clone(), zero.clone()],
                    [one.clone(), zero.clone(), one.clone()],
                ],
                vec![one.clone(), one],
                &Decider::default(),
            )
            .unwrap(),
            closed: Some(false),
            self_intersect: Some(false),
        }
    }

    #[test]
    fn revolution_construction_uses_the_shared_allowance_and_deadline() {
        let origin = [Rational::zero(), Rational::zero(), Rational::zero()];
        let axis = [Rational::zero(), Rational::zero(), Rational::one()];
        let profile = meridian();
        let (result, work) = catallax_numbers::with_operation_budget(1, || {
            revolved(profile, &origin, &axis, &Decider::default())
        });
        assert!(matches!(result, Err(ImportRefusal::ResourceLimit)));
        assert!(work.exhausted);
        let mut dec = Decider::default();
        dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
        assert!(matches!(
            revolved(meridian(), &origin, &axis, &dec),
            Err(ImportRefusal::ResourceLimit)
        ));
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
