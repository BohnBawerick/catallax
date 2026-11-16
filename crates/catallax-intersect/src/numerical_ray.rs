//! Numerical finite ray meetings. No exact incidence or complete hit-count proof.

use catallax_numbers::{
    Ball, Decider, NumericalCheck, NumericalEvidence, NumericalSample, Rational, Real, Refusal,
    Sign,
};
use catallax_shapes::{Point3, RationalPatch, Vector3};
use catallax_trace::{
    curve_surface_global::global_curve_intersect_with_work,
    extent::{Extent1, Extent2},
    solution::SolutionType,
    spline::{BSplineCurve, BSplineSurface, ControlPoint},
    surface_work::WorkReport,
    vector::Vector2,
};

/// One reported meeting checked on the original finite patch and ray.
#[derive(Clone, Debug)]
pub struct NumericalRayHit {
    /// Distance along the internally normalized ray, in millimetres.
    pub distance_mm: Rational,
    /// Retained native patch parameters.
    pub uv: [Rational; 2],
    /// Natural patch normal dotted with the normalized ray direction.
    pub normal_dot: Ball,
}

/// Every reported finite meeting and its numerical checks.
#[derive(Clone, Debug)]
pub struct NumericalRayMeetings {
    /// All checked reports, including chart-boundary duplicates and tangencies.
    pub hits: Vec<NumericalRayHit>,
    /// Preset tolerance and residuals checked against the original definitions.
    pub evidence: NumericalEvidence,
    /// Floating search checkpoint work, charged to the caller.
    pub search_work: WorkReport,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("intersect.numerical-ray", detail)
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value).ok_or_else(|| bad("non-finite numerical ray value"))
}

fn approximate(value: &Real, dec: &mut Decider) -> Result<f64, Refusal> {
    dec.charge_work(1)?;
    let value = dec.enclose_budgeted_at(value, 128)?.mid().to_f64_nearest();
    exact(value)?;
    Ok(value)
}

fn pay(units: u64, deadline: Option<std::time::Instant>) -> bool {
    let mut dec = Decider::default();
    dec.set_deadline(deadline);
    dec.charge_work(units).is_ok()
}

/// Searches every finite source span with SMLib's bounded line/surface solver.
/// Every reported point is checked on the original patch and normalized ray.
/// Overlap, incomplete searches and residuals outside the preset tolerance refuse.
/// An empty numerical result supplies no exact exclusion or branch-count proof.
pub fn patch_meetings(
    patch: &RationalPatch,
    origin: &Point3,
    direction: &Vector3,
    tolerance_mm: &Rational,
    dec: &mut Decider,
) -> Result<NumericalRayMeetings, Refusal> {
    if tolerance_mm.sign() != Sign::Positive {
        return Err(bad("the ray tolerance must be positive"));
    }
    let tolerance = tolerance_mm.to_f64_nearest();
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(bad("the ray tolerance is outside finite binary64"));
    }
    let norm = direction.norm_squared().sqrt();
    if dec.enclose_budgeted_at(&norm, 128)?.lo().sign() != Sign::Positive {
        return Err(bad("the numerical ray direction is unresolved or zero"));
    }
    let direction = direction.scaled(&(Real::one() / norm));
    let source = patch.carrier();
    let mut controls = Vec::new();
    let mut end = Rational::zero();
    for (point, weight) in source.controls().iter().zip(source.weights()) {
        dec.charge_work(1)?;
        let projected = point.minus(origin).dot(&direction);
        end = end.max(dec.enclose_budgeted_at(&projected, 128)?.hi());
        let weight = approximate(weight, dec)?;
        if weight <= 0. {
            return Err(bad("source adaptation lost a positive patch weight"));
        }
        let [x, y, z] = [&point.x, &point.y, &point.z]
            .map(|value| approximate(value, dec))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .unwrap();
        let coordinates = [x * weight, y * weight, z * weight, weight];
        if coordinates.iter().any(|v| !v.is_finite()) {
            return Err(bad("a homogeneous patch control overflowed binary64"));
        }
        controls.push(ControlPoint { coordinates });
    }
    let mut knots = Vec::new();
    for row in source.knots() {
        let mut values = Vec::new();
        for knot in row {
            dec.charge_work(1)?;
            let value = knot.to_f64_nearest();
            exact(value)?;
            values.push(value);
        }
        knots.push(values);
    }
    let surface = BSplineSurface::new(
        source.degrees(),
        knots.try_into().unwrap(),
        source.control_count(),
        controls,
    )
    .map_err(|error| bad(format!("finite patch adaptation failed: {error:?}")))?;
    // Positive weights bound the complete source by its Cartesian control hull.
    // The extra millimetre keeps a closed end away from the solver's endpoint.
    let end = end + Rational::one();
    let endpoint = origin.plus(&direction.scaled(&Real::from(end.clone())));
    let mut points = Vec::new();
    for point in [origin, &endpoint] {
        let coordinates = [
            approximate(&point.x, dec)?,
            approximate(&point.y, dec)?,
            approximate(&point.z, dec)?,
            1.,
        ];
        points.push(ControlPoint { coordinates });
    }
    let hi = end.to_f64_nearest();
    let curve = BSplineCurve::new(1, vec![0., 0., hi, hi], points)
        .map_err(|error| bad(format!("finite ray adaptation failed: {error:?}")))?;
    let domain = patch.domain();
    let float_domain = Extent2::new(
        Vector2::new(domain[0][0].to_f64_nearest(), domain[1][0].to_f64_nearest()),
        Vector2::new(domain[0][1].to_f64_nearest(), domain[1][1].to_f64_nearest()),
    );
    let (solutions, search_work) = global_curve_intersect_with_work(
        &surface,
        float_domain,
        &curve,
        Extent1::new(0., hi),
        tolerance / 4.,
        dec.deadline(),
        pay,
    );
    if search_work.cancelled {
        dec.charge_work(0)?;
        return Err(bad("finite ray search cancelled without a retained result"));
    }
    let solutions =
        solutions.map_err(|error| bad(format!("finite ray search failed: {error:?}")))?;
    let mut hits = Vec::new();
    let mut samples = Vec::new();
    for solution in solutions.as_slice() {
        dec.charge_work(1)?;
        if !matches!(
            solution.kind,
            SolutionType::SingleValue | SolutionType::PolySingleValue
        ) {
            return Err(bad("a finite ray overlap needs interval ownership"));
        }
        let [t, u, v] = solution.start.parameters[..3].try_into().unwrap();
        let t = exact(t)?;
        let uv = [exact(u)?, exact(v)?];
        if t < Rational::zero()
            || t > end
            || (0..2).any(|axis| uv[axis] < domain[axis][0] || uv[axis] > domain[axis][1])
        {
            return Err(bad("a numerical ray meeting leaves its finite domains"));
        }
        let at = uv.clone().map(Real::from);
        let actual = patch.point_at(&at, dec)?;
        let on_ray = origin.plus(&direction.scaled(&Real::from(t.clone())));
        let residual = actual.minus(&on_ray).norm_squared().sqrt();
        let distance_upper_mm = dec.enclose_budgeted_at(&residual, 128)?.hi();
        let [du, dv] = patch.derivatives_at(&at, dec)?;
        let normal_dot = dec.enclose_budgeted_at(&du.cross(&dv).dot(&direction), 128)?;
        samples.push(NumericalSample {
            parameters: vec![t.clone(), uv[0].clone(), uv[1].clone()],
            distance_upper_mm,
        });
        hits.push(NumericalRayHit {
            distance_mm: t,
            uv,
            normal_dot,
        });
    }
    let mut checks = vec![NumericalCheck::Passed {
        name: "ray.retained-every-reported-finite-hit".into(),
        scope: "numerical-search-no-complete-hit-count-proof".into(),
    }];
    if !samples.is_empty() {
        checks.push(NumericalCheck::Samples {
            name: "ray.original-patch-incidence".into(),
            scope: "finite-source-chart".into(),
            samples,
        });
    }
    let evidence = NumericalEvidence::new(tolerance_mm.clone(), checks, vec![])?;
    dec.check_deadline()?;
    Ok(NumericalRayMeetings {
        hits,
        evidence,
        search_work,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;
    use catallax_shapes::BSplineSurface;

    #[test]
    fn finite_numerical_ray_keeps_hits_and_checks_original_patch_residuals() {
        let q = |s| Rational::parse(s).unwrap();
        let mut dec = Decider::default();
        let controls = [(0, 0), (1, 2), (2, 0)]
            .into_iter()
            .flat_map(|(x, y)| {
                (0..2).map(move |z| {
                    Point3::new(
                        Rational::from_ratio(x, 2).unwrap().into(),
                        Real::from_i64(y),
                        Real::from_i64(z),
                    )
                })
            })
            .collect();
        let patch = RationalPatch::new(
            BSplineSurface::new(
                [2, 1],
                [
                    vec![q("0"), q("0"), q("0"), q("1"), q("1"), q("1")],
                    vec![q("0"), q("0"), q("1"), q("1")],
                ],
                [3, 2],
                controls,
                vec![Real::one(); 6],
                &mut dec,
            )
            .unwrap(),
            &mut dec,
        )
        .unwrap();
        let direction = Vector3::new(Real::zero(), Real::one(), Real::zero());
        let ((), work) = with_operation_budget(20_000_000, || {
            for (x, y, count) in [("1/2", "1/2", 1), ("1/2", "3/2", 0), ("3/2", "1/2", 0)] {
                let origin = Point3::new(q(x).into(), q(y).into(), q("1/2").into());
                let result =
                    patch_meetings(&patch, &origin, &direction, &q("1/10000000"), &mut dec)
                        .unwrap();
                assert_eq!(result.hits.len(), count);
                for hit in &result.hits {
                    assert!((&hit.distance_mm - &q("1/2")).abs() <= q("1/10000000"));
                    assert!(hit.normal_dot.hi() < q("0"));
                }
                assert!(!result.search_work.cancelled);
                println!(
                    "NUMERICAL_RAY {x} {y} hits={count} work={:?}",
                    result.search_work
                );
            }
        });
        assert!(!work.exhausted, "{work:?}");
        let origin = Point3::new(q("1/2").into(), q("1/2").into(), q("1/2").into());
        let (result, work) = with_operation_budget(16, || {
            patch_meetings(&patch, &origin, &direction, &q("1/10000000"), &mut dec)
        });
        assert!(result.is_err() && work.exhausted);
    }
}
