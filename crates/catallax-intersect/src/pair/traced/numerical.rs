//! Approximate native windows, without exact carrier or incidence authority.

use super::*;
use catallax_numbers::{
    Ball, NumericalCheck, NumericalEvidence, NumericalSample, operation_budget_active,
    with_operation_budget,
};

const BITS: u32 = 128;

mod checked;
mod cylinder;
pub use cylinder::{NumericalNativeCylinderSource, numerical_native_cylinder_source};
mod finite;
pub use checked::{
    NumericalIntersectionPiece, NumericalNativeIntersection, NumericalNativeWork,
    numerical_native_intersection,
};
pub use finite::numerical_native_finite_intersection;

#[cfg(test)]
mod fit_tests;

/// Native cone/cylinder parameter windows with their original analytic sources.
/// Checks and rounded controls remain numerical. This is not an exact pair plan,
/// a branch count, a Boolean result or a complete body deviation bound.
#[derive(Clone, Debug)]
pub struct NumericalNativeWindows {
    originals: [Surface; 2],
    sources: [BSplineSurface; 2],
    evidence: NumericalEvidence,
    cone_heights: [Rational; 2],
    cylinder_parameters: [Rational; 2],
    cylinder_origin: Point3,
    control_rounding_mm: Rational,
}

impl NumericalNativeWindows {
    /// Original carriers, in the caller's operand order.
    pub fn originals(&self) -> &[Surface; 2] {
        &self.originals
    }
    /// Rounded native definitions over their complete finite parameter domains.
    pub fn sources(&self) -> &[BSplineSurface; 2] {
        &self.sources
    }
    /// Preset tolerance and actual carrier/control checks, not exact incidence.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
    /// Strict positive cone-height window searched by this query.
    pub fn cone_heights(&self) -> &[Rational; 2] {
        &self.cone_heights
    }
    /// Signed millimetres along the unit cylinder axis from `cylinder_origin`.
    pub fn cylinder_parameters(&self) -> &[Rational; 2] {
        &self.cylinder_parameters
    }
    /// Origin used by the retained cylinder source. Global support searches use
    /// the cone apex projection. Finite face searches use the original cylinder
    /// origin so different cone faces share one cylinder parameter frame.
    pub fn cylinder_origin(&self) -> &Point3 {
        &self.cylinder_origin
    }
    /// Maximum L1 rounding error of one Cartesian control, in millimetres.
    /// The original positive weights and native knots are unchanged.
    pub fn control_rounding_mm(&self) -> &Rational {
        &self.control_rounding_mm
    }
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("intersect.numerical-native-window", detail)
}

fn enclosure(value: &Real, dec: &mut Decider) -> Result<Ball, Refusal> {
    dec.enclose_budgeted_at(value, BITS)
}

fn absolute_upper(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}

fn absolute_lower(value: &Ball) -> Rational {
    if value.lo().sign() == Sign::Positive {
        value.lo()
    } else if value.hi().sign() == Sign::Negative {
        -value.hi()
    } else {
        r(0)
    }
}

fn unit(vector: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    let norm = vector.norm_squared().sqrt();
    if sign(&norm, dec)? != Sign::Positive {
        return Err(bad("a native axis or reference direction vanished"));
    }
    Ok(vector.scaled(&(&Real::one() / &norm)))
}

struct Frame {
    origin: Point3,
    axis: Vector3,
    x: Vector3,
    y: Vector3,
}

impl Frame {
    fn new(
        origin: Point3,
        axis: &Vector3,
        reference: &Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let axis = unit(axis, dec)?;
        let x = unit(&reference.minus(&axis.scaled(&reference.dot(&axis))), dec)?;
        let y = axis.cross(&x);
        Ok(Self { origin, axis, x, y })
    }
    fn point(&self, x: &Real, y: &Real, height: &Real) -> Point3 {
        self.origin
            .plus(&self.x.scaled(x))
            .plus(&self.y.scaled(y))
            .plus(&self.axis.scaled(height))
    }
}

struct Source {
    surface: BSplineSurface,
    rounding: Rational,
    check: NumericalCheck,
}

fn source(
    frame: &Frame,
    heights: &[Rational; 2],
    radii: [&Real; 2],
    scope: &str,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Source, Refusal> {
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    let mut samples = Vec::new();
    let mut rounding = r(0);
    for (u, (x, y, weight)) in circle().into_iter().enumerate() {
        for v in 0..2 {
            let point = frame.point(
                &(radii[v] * &x),
                &(radii[v] * &y),
                &Real::from(heights[v].clone()),
            );
            let mut values = Vec::new();
            let mut error = r(0);
            for coordinate in [&point.x, &point.y, &point.z] {
                let bound = enclosure(coordinate, dec)?;
                let value = bound.mid();
                error = error
                    + (&value - &bound.lo())
                        .abs()
                        .max((&bound.hi() - &value).abs());
                values.push(Real::from(value));
            }
            if error > tolerance.clone() / r(16) {
                return Err(bad(
                    "native control rounding exceeds its share of the preset tolerance",
                ));
            }
            rounding = rounding.max(error.clone());
            samples.push(NumericalSample {
                parameters: vec![Rational::from_i64(u as i64), Rational::from_i64(v as i64)],
                distance_upper_mm: error,
            });
            controls.push(Point3::new(
                values[0].clone(),
                values[1].clone(),
                values[2].clone(),
            ));
            weights.push(weight.clone());
        }
    }
    let check = NumericalCheck::Samples {
        name: "native.control-cartesian-rounding".into(),
        scope: scope.into(),
        samples,
    };
    let surface = BSplineSurface::new(
        [2, 1],
        [knots(), [0, 0, 1, 1].map(Rational::from_i64).to_vec()],
        [9, 2],
        controls,
        weights,
        dec,
    )?;
    Ok(Source {
        surface,
        rounding,
        check,
    })
}

fn carrier_checks(
    source: &BSplineSurface,
    frame: &Frame,
    radius: Option<&Real>,
    slope: Option<&Real>,
    scope: &str,
    checks: &mut Vec<NumericalCheck>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let scalar = |q: &Rational| Ball::from_rational(q, BITS);
    let origin = [&frame.origin.x, &frame.origin.y, &frame.origin.z]
        .map(|v| enclosure(v, dec))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let axis = [&frame.axis.x, &frame.axis.y, &frame.axis.z]
        .map(|v| enclosure(v, dec))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let radius = radius.map(|v| enclosure(v, dec)).transpose()?;
    let slope = slope.map(|v| enclosure(v, dec)).transpose()?;
    let mut samples = Vec::new();
    for i in 0..=8 {
        let u = Rational::from_ratio(i, 2).unwrap();
        for j in 0..=8 {
            dec.charge_work(1)?;
            let v = Rational::from_ratio(j, 8).unwrap();
            let point = source.evaluate(&scalar(&u), &scalar(&v), dec)?;
            let offset: [Ball; 3] = std::array::from_fn(|k| &point[k] - &origin[k]);
            let height =
                (&offset[0] * &axis[0]) + (&offset[1] * &axis[1]) + (&offset[2] * &axis[2]);
            let radial: [Ball; 3] = std::array::from_fn(|k| &offset[k] - &(&height * &axis[k]));
            let distance = (radial[0].square() + radial[1].square() + radial[2].square())
                .sqrt()
                .ok_or_else(|| bad("native radial distance was unresolved"))?;
            let residual = if let Some(radius) = &radius {
                &distance - radius
            } else {
                &distance - &(&height * slope.as_ref().expect("cone slope"))
            };
            samples.push(NumericalSample {
                parameters: vec![u.clone(), v],
                distance_upper_mm: absolute_upper(&residual),
            });
        }
    }
    checks.push(NumericalCheck::Samples {
        name: "native.sampled-analytic-carrier-residual".into(),
        scope: scope.into(),
        samples,
    });
    Ok(())
}

/// Builds numerical cone/cylinder windows in either operand order at 1e-7 mm.
/// Algebraic axes and coefficients remain supported through interval evaluation.
/// Full-support windows that can reach a parallel generator or the cone apex
/// require finite-face/pole handling. This query publishes no exact model facts.
pub fn numerical_native_windows(
    surfaces: [&Surface; 2],
    dec: &mut Decider,
) -> Result<Option<NumericalNativeWindows>, Refusal> {
    if operation_budget_active() {
        return build(surfaces, dec);
    }
    let (result, _) = with_operation_budget(20_000_000, || build(surfaces, dec));
    result
}

fn build(
    surfaces: [&Surface; 2],
    dec: &mut Decider,
) -> Result<Option<NumericalNativeWindows>, Refusal> {
    let (cone, cylinder) = match surfaces {
        [Surface::Cone(cone), Surface::Cylinder(cylinder)] => (cone, cylinder),
        [Surface::Cylinder(cylinder), Surface::Cone(cone)] => (cone, cylinder),
        _ => return Ok(None),
    };
    dec.charge_work(1)?;
    if dec.policy().max_precision_bits < BITS {
        return Err(bad("native numerical windows require 128-bit enclosures"));
    }
    let radius = cylinder.radius();
    let radius_bound = enclosure(&radius, dec)?;
    let slope_bound = enclosure(&cone.slope, dec)?;
    if radius_bound.lo().sign() != Sign::Positive || slope_bound.lo().sign() != Sign::Positive {
        return Err(bad(
            "native windows require positive cone slope and cylinder radius",
        ));
    }
    let cone_frame = Frame::new(cone.apex.clone(), &cone.axis, &cone.ref_dir, dec)?;
    let cylinder_axis = unit(&cylinder.axis, dec)?;
    let projected = cylinder
        .origin
        .plus(&cylinder_axis.scaled(&cone.apex.minus(&cylinder.origin).dot(&cylinder_axis)));
    let cylinder_frame = Frame::new(projected.clone(), &cylinder.axis, &cylinder.ref_dir, dec)?;
    let radial = cone.apex.minus(&projected);
    let distance_upper = [&radial.x, &radial.y, &radial.z]
        .into_iter()
        .map(|v| enclosure(v, dec).map(|b| absolute_upper(&b)))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .fold(r(0), |a, b| a + b);
    let cosine = cone_frame.axis.dot(&cylinder_frame.axis);
    let sine = cone_frame
        .axis
        .cross(&cylinder_frame.axis)
        .norm_squared()
        .sqrt();
    let cosine_bound = enclosure(&cosine, dec)?;
    let cosine_abs = Ball::from_rational(&absolute_lower(&cosine_bound), BITS)
        .hull(&Ball::from_rational(&absolute_upper(&cosine_bound), BITS));
    let gap = enclosure(&sine, dec)? - (&slope_bound * &cosine_abs);
    let gap_lower = absolute_lower(&gap);
    if gap_lower.sign() != Sign::Positive {
        return Err(bad(
            "near a parallel cone generator, finite face domains are required",
        ));
    }
    let apex_gap = &radial.norm_squared() - &cylinder.radius_sq;
    let apex_gap_lower = absolute_lower(&enclosure(&apex_gap, dec)?);
    if apex_gap_lower.sign() != Sign::Positive {
        return Err(bad(
            "a window reaching the cone apex needs singular-point handling",
        ));
    }
    let upper = (&distance_upper + &radius_bound.hi()) / gap_lower;
    let lower = apex_gap_lower
        / (r(2) * (&(&distance_upper + &radius_bound.hi()) + &r(1)) * (&r(1) + &slope_bound.hi()));
    let heights = [lower, &upper + &r(1)];
    if heights[0] >= heights[1] {
        return Err(bad("native cone-height window collapsed"));
    }
    let cylinder_bound = (&r(1) + &slope_bound.hi()) * heights[1].clone() + r(1);
    let cylinder_parameters = [-&cylinder_bound, cylinder_bound];
    let checks = vec![NumericalCheck::Passed {
        name: "native.retained-full-parameter-window".into(),
        scope: "cone-cylinder".into(),
    }];
    assemble(
        surfaces,
        heights,
        cylinder_parameters,
        cone_frame,
        cylinder_frame,
        checks,
        dec,
    )
}

fn assemble(
    surfaces: [&Surface; 2],
    heights: [Rational; 2],
    cylinder_parameters: [Rational; 2],
    cone_frame: Frame,
    cylinder_frame: Frame,
    mut checks: Vec<NumericalCheck>,
    dec: &mut Decider,
) -> Result<Option<NumericalNativeWindows>, Refusal> {
    let (cone, cylinder, swap) = match surfaces {
        [Surface::Cone(cone), Surface::Cylinder(cylinder)] => (cone, cylinder, false),
        [Surface::Cylinder(cylinder), Surface::Cone(cone)] => (cone, cylinder, true),
        _ => return Ok(None),
    };
    let radius = cylinder.radius();
    let tolerance = Rational::from_ratio(1, 10_000_000).unwrap();
    let radii = heights
        .each_ref()
        .map(|h| &cone.slope * &Real::from(h.clone()));
    let cone_source = source(
        &cone_frame,
        &heights,
        [&radii[0], &radii[1]],
        "cone",
        &tolerance,
        dec,
    )?;
    let cylinder_source = source(
        &cylinder_frame,
        &cylinder_parameters,
        [&radius, &radius],
        "cylinder",
        &tolerance,
        dec,
    )?;
    let rounding = cone_source.rounding.max(cylinder_source.rounding);
    checks.extend([cone_source.check, cylinder_source.check]);
    let cone_source = cone_source.surface;
    let cylinder_source = cylinder_source.surface;
    carrier_checks(
        &cone_source,
        &cone_frame,
        None,
        Some(&cone.slope),
        "cone",
        &mut checks,
        dec,
    )?;
    carrier_checks(
        &cylinder_source,
        &cylinder_frame,
        Some(&radius),
        None,
        "cylinder",
        &mut checks,
        dec,
    )?;
    let evidence = NumericalEvidence::new(tolerance, checks, vec![])?;
    let sources = if swap {
        [cylinder_source, cone_source]
    } else {
        [cone_source, cylinder_source]
    };
    Ok(Some(NumericalNativeWindows {
        originals: surfaces.map(Clone::clone),
        sources,
        evidence,
        cone_heights: heights,
        cylinder_parameters,
        cylinder_origin: cylinder_frame.origin,
        control_rounding_mm: rounding,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::current_operation_work;

    fn carriers(axis: Vector3, offset: Rational) -> [Surface; 2] {
        [
            Surface::Cone(Cone {
                apex: zero_point(),
                axis: z_axis(),
                slope: s(1),
                ref_dir: Vector3::new(s(1), s(0), s(0)),
            }),
            Surface::Cylinder(Cylinder {
                origin: Point3::new(offset.into(), s(0), s(0)),
                axis,
                radius_sq: s(1),
                ref_dir: Vector3::new(s(0), s(1), s(0)),
            }),
        ]
    }

    #[test]
    fn axial_and_algebraic_native_windows_remain_numerical_in_both_orders() {
        let ((), work) = with_operation_budget(20_000_000, || {
            for axis in [z_axis(), Vector3::new(s(1), s(-1), s(1))] {
                for offset in [r(10), Rational::from_ratio(21, 2).unwrap()] {
                    let surfaces = carriers(axis.clone(), offset);
                    for reverse in [false, true] {
                        let ordered = if reverse {
                            [&surfaces[1], &surfaces[0]]
                        } else {
                            [&surfaces[0], &surfaces[1]]
                        };
                        let windows = numerical_native_windows(ordered, &mut Decider::default())
                            .unwrap()
                            .unwrap();
                        assert_eq!(
                            windows.evidence().tolerance_mm(),
                            &Rational::from_ratio(1, 10_000_000).unwrap()
                        );
                        assert!(windows.cone_heights()[0] > r(0));
                        assert!(windows.cone_heights()[1] > windows.cone_heights()[0]);
                        assert!(
                            windows.cylinder_parameters()[0] < r(0)
                                && windows.cylinder_parameters()[1] > r(0)
                        );
                        assert!(
                            windows.control_rounding_mm()
                                <= &Rational::from_ratio(1, 160_000_000).unwrap()
                        );
                        assert_eq!(
                            windows.sources()[0].domain(),
                            [[&r(0), &r(4)], [&r(0), &r(1)]]
                        );
                        assert_eq!(
                            matches!(windows.originals()[0], Surface::Cylinder(_)),
                            reverse
                        );
                        let sample_count = windows
                            .evidence()
                            .checks()
                            .iter()
                            .filter_map(|c| match c {
                                NumericalCheck::Samples { samples, .. } => Some(samples.len()),
                                _ => None,
                            })
                            .sum::<usize>();
                        assert_eq!(sample_count, 2 * 18 + 2 * 81);
                    }
                }
            }
            eprintln!(
                "native numerical window controls total_work={}",
                current_operation_work().used
            );
        });
        assert!(!work.exhausted, "{work:?}");
    }

    #[test]
    fn a_numerical_window_cannot_bypass_the_request_allowance() {
        let surfaces = carriers(z_axis(), r(10));
        let (result, work) = with_operation_budget(0, || {
            numerical_native_windows([&surfaces[0], &surfaces[1]], &mut Decider::default())
        });
        assert!(result.is_err());
        assert!(work.exhausted);
    }
}
