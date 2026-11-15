//! Cap curves in the retained native parameter of a cylindrical source quarter.

use super::CylinderSpan;
use catallax_numbers::{
    Ball, Decider, NumericalCheck, NumericalEvidence, NumericalSample, Rational, Real, Refusal,
    Sign, operation_budget_active, with_operation_budget,
};
use catallax_shapes::{BSplineCurve, BSplineSurface, Plane, Point3, Vector3};
use catallax_topology::numerical::{DEFAULT_WORK, NumericalOptions};

const BITS: u32 = 128;
const SUBDIVISIONS: i64 = 16;

/// A checked numerical cap edge and its two retained parameter curves.
/// Samples grant no exact incidence or whole-curve deviation bound.
#[derive(Clone, Debug)]
pub struct CapCurve {
    key: String,
    curve: BSplineCurve,
    plane_uv: BSplineCurve,
    cylinder_uv: BSplineCurve,
    evidence: NumericalEvidence,
}

impl CapCurve {
    /// Stable caller-supplied identity for assembly.
    pub fn key(&self) -> &str {
        &self.key
    }
    /// Three-dimensional edge in the source quarter's native U parameter.
    pub fn curve(&self) -> &BSplineCurve {
        &self.curve
    }
    /// Plane coordinates at that same native edge parameter.
    pub fn plane_uv(&self) -> &BSplineCurve {
        &self.plane_uv
    }
    /// Cylinder-quarter coordinates at that same native edge parameter.
    pub fn cylinder_uv(&self) -> &BSplineCurve {
        &self.cylinder_uv
    }
    /// Preset tolerance and actual sample checks against both original carriers
    /// and both retained parameter compositions.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("hole.numerical-cap", detail)
}

/// Projects one canonical quadratic/linear cylinder quarter along its unit
/// axis onto a cap plane. The quarter's native knots and finite U interval stay
/// fixed; V is computed in its supplied signed-mm cylinder interval.
/// Parallel cap planes return `None`. Invalid sources and failed checks refuse.
pub fn cap_curve(
    key: &str,
    span: &CylinderSpan,
    quarter: &BSplineSurface,
    plane: &Plane,
    options: &NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<CapCurve>, Refusal> {
    let mut run = || build(key, span, quarter, plane, options, dec);
    if operation_budget_active() {
        return run();
    }
    let (result, work) = with_operation_budget(DEFAULT_WORK, run);
    if work.exhausted {
        return Err(bad("cap construction exhausted its shared work allowance"));
    }
    result
}

fn lift(values: &[Real; 3]) -> [Real; 4] {
    let two = Real::from_i64(2);
    let three = Real::from_i64(3);
    [
        values[0].clone(),
        &(&values[0] + &(&two * &values[1])) / &three,
        &(&(&two * &values[1]) + &values[2]) / &three,
        values[2].clone(),
    ]
}

fn bounds(point: &Point3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_budgeted_at(&point.x, BITS)?,
        dec.enclose_budgeted_at(&point.y, BITS)?,
        dec.enclose_budgeted_at(&point.z, BITS)?,
    ])
}

fn vector_bounds(vector: &Vector3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    bounds(
        &Point3::new(vector.x.clone(), vector.y.clone(), vector.z.clone()),
        dec,
    )
}

fn magnitude(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}

fn difference(a: &[Ball; 3], b: &[Ball; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| sum + magnitude(&(&a[i] - &b[i])))
}

fn build(
    key: &str,
    span: &CylinderSpan,
    quarter: &BSplineSurface,
    plane: &Plane,
    options: &NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<CapCurve>, Refusal> {
    dec.charge_work(1)?;
    if key.trim().is_empty() || options.tolerance_mm.sign() != Sign::Positive {
        return Err(bad(
            "cap identity and positive preset tolerance are required",
        ));
    }
    if quarter.degrees() != [2, 1] || quarter.control_count() != [3, 2] {
        return Err(bad(
            "cap construction requires one quadratic/linear cylinder quarter",
        ));
    }
    let height = &span.interval_mm[1] - &span.interval_mm[0];
    if height.sign() != Sign::Positive {
        return Err(bad("the signed-mm cylinder interval must increase"));
    }
    Plane::new(
        plane.origin.clone(),
        plane.normal.clone(),
        plane.ref_dir.clone(),
        dec,
    )?;
    let axis = span
        .cylinder
        .axis
        .scaled(&(&Real::one() / &span.cylinder.axis.norm_squared().sqrt()));
    let denominator = axis.dot(&plane.normal);
    if denominator
        .as_rational()
        .is_some_and(|value| value.is_zero())
    {
        return Ok(None);
    }
    let denominator_bound = dec.enclose_budgeted_at(&denominator, BITS)?;
    if denominator_bound.lo().is_zero() && denominator_bound.hi().is_zero() {
        return Ok(None);
    }
    if denominator_bound.lo().sign() != Sign::Positive
        && denominator_bound.hi().sign() != Sign::Negative
    {
        return Err(bad(
            "the projection direction is unresolved at the working precision",
        ));
    }
    let normal = plane
        .normal
        .scaled(&(&Real::one() / &plane.normal.norm_squared().sqrt()));
    let x = plane
        .ref_dir
        .scaled(&(&Real::one() / &plane.ref_dir.norm_squared().sqrt()));
    let y = normal.cross(&x);
    let [[u0, u1], [v0, v1]] = quarter.domain();
    let weights: [Real; 3] = std::array::from_fn(|i| quarter.weights()[2 * i].clone());
    let mut points = Vec::new();
    let mut plane_points = Vec::new();
    let mut v = Vec::new();
    for (i, weight) in weights.iter().enumerate() {
        dec.charge_work(1)?;
        if !quarter.weights()[2 * i + 1].same_definition(weight) {
            return Err(bad("axial source rows must retain the same weights"));
        }
        let point = &quarter.controls()[2 * i];
        let offset = &plane.origin.minus(point).dot(&plane.normal) / &denominator;
        let projected = point.plus(&axis.scaled(&offset));
        let relative = projected.minus(&plane.origin);
        plane_points.push(Point3::new(
            relative.dot(&x),
            relative.dot(&y),
            Real::zero(),
        ));
        v.push(
            &Real::from(v0.clone())
                + &(&(&offset / &Real::from(height.clone())) * &Real::from(v1 - v0)),
        );
        points.push(projected);
    }
    let curve = BSplineCurve::new(
        2,
        quarter.knots()[0].to_vec(),
        points,
        weights.to_vec(),
        dec,
    )?;
    let plane_uv = BSplineCurve::new(
        2,
        quarter.knots()[0].to_vec(),
        plane_points,
        weights.to_vec(),
        dec,
    )?;
    let weighted_v = std::array::from_fn(|i| &v[i] * &weights[i]);
    let elevated_weights = lift(&weights);
    let elevated_v = lift(&weighted_v);
    let a = Real::from(u0.clone());
    let b = Real::from(u1.clone());
    let two = Real::from_i64(2);
    let three = Real::from_i64(3);
    // Multiplying the quadratic denominator by native U needs degree three.
    let homogeneous_u = [
        &a * &weights[0],
        &(&(&b * &weights[0]) + &(&(&two * &a) * &weights[1])) / &three,
        &(&(&(&two * &b) * &weights[1]) + &(&a * &weights[2])) / &three,
        &b * &weights[2],
    ];
    let cylinder_uv = BSplineCurve::new(
        3,
        [vec![u0.clone(); 4], vec![u1.clone(); 4]].concat(),
        (0..4)
            .map(|i| {
                Point3::new(
                    match i {
                        0 => a.clone(),
                        3 => b.clone(),
                        _ => &homogeneous_u[i] / &elevated_weights[i],
                    },
                    &elevated_v[i] / &elevated_weights[i],
                    Real::zero(),
                )
            })
            .collect(),
        elevated_weights.to_vec(),
        dec,
    )?;
    let evidence = checks(
        key,
        span,
        quarter,
        (plane, [&x, &y]),
        [&curve, &plane_uv, &cylinder_uv],
        options,
        dec,
    )?;
    Ok(Some(CapCurve {
        key: key.into(),
        curve,
        plane_uv,
        cylinder_uv,
        evidence,
    }))
}

fn checks(
    key: &str,
    span: &CylinderSpan,
    quarter: &BSplineSurface,
    frame: (&Plane, [&Vector3; 2]),
    curves: [&BSplineCurve; 3],
    options: &NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalEvidence, Refusal> {
    let (plane, basis) = frame;
    let [curve, plane_uv, cylinder_uv] = curves;
    let origin = bounds(&plane.origin, dec)?;
    let x = vector_bounds(basis[0], dec)?;
    let y = vector_bounds(basis[1], dec)?;
    let axis = span
        .cylinder
        .axis
        .scaled(&(&Real::one() / &span.cylinder.axis.norm_squared().sqrt()));
    let normal = plane
        .normal
        .scaled(&(&Real::one() / &plane.normal.norm_squared().sqrt()));
    let original_origin = bounds(&span.cylinder.origin, dec)?;
    let axis_box = vector_bounds(&axis, dec)?;
    let normal_box = vector_bounds(&normal, dec)?;
    let radius_sq = dec.enclose_budgeted_at(&span.cylinder.radius_sq, BITS)?;
    let radius_lower = dec.enclose_budgeted_at(&span.cylinder.radius(), BITS)?.lo();
    if radius_lower.sign() != Sign::Positive {
        return Err(bad("a positive cylinder radius is required"));
    }
    let mut control_samples = Vec::new();
    for i in 0..3 {
        for j in 0..2 {
            dec.charge_work(1)?;
            let point = bounds(&quarter.controls()[2 * i + j], dec)?;
            let height = (0..3).fold(Ball::from_rational(&Rational::zero(), BITS), |sum, k| {
                sum + &(&point[k] - &original_origin[k]) * &axis_box[k]
            });
            let residual = &height - &Ball::from_rational(&span.interval_mm[j], BITS);
            control_samples.push(NumericalSample {
                parameters: vec![Rational::from_i64(i as i64), Rational::from_i64(j as i64)],
                distance_upper_mm: magnitude(&residual),
            });
        }
    }
    let mut samples: [Vec<NumericalSample>; 4] = std::array::from_fn(|_| Vec::new());
    let [[u0, u1], [v0, v1]] = quarter.domain();
    for i in 0..=SUBDIVISIONS {
        dec.charge_work(1)?;
        let t = u0 + &((u1 - u0) * Rational::from_ratio(i, SUBDIVISIONS).unwrap());
        let parameter = Ball::from_rational(&t, BITS);
        let p = curve.evaluate(&parameter, dec)?;
        let uv = plane_uv.evaluate(&parameter, dec)?;
        // Clamped spline endpoints are their retained control points. Evaluating
        // the homogeneous quotient there can widen an integer U past its domain.
        let native = if i == 0 {
            bounds(cylinder_uv.controls().next().unwrap(), dec)?
        } else if i == SUBDIVISIONS {
            bounds(cylinder_uv.controls().last().unwrap(), dec)?
        } else {
            cylinder_uv.evaluate(&parameter, dec)?
        };
        for (value, range) in [(&native[0], [u0, u1]), (&native[1], [v0, v1])] {
            if value.lo() < *range[0] || value.hi() > *range[1] {
                return Err(bad("a cap UV sample left its finite cylinder quarter"));
            }
        }
        let on_plane =
            std::array::from_fn(|j| &(&origin[j] + &(&x[j] * &uv[0])) + &(&y[j] * &uv[1]));
        let on_cylinder = quarter.evaluate(&native[0], &native[1], dec)?;
        let delta: [Ball; 3] = std::array::from_fn(|j| &p[j] - &original_origin[j]);
        let radial = [
            &(&delta[1] * &axis_box[2]) - &(&delta[2] * &axis_box[1]),
            &(&delta[2] * &axis_box[0]) - &(&delta[0] * &axis_box[2]),
            &(&delta[0] * &axis_box[1]) - &(&delta[1] * &axis_box[0]),
        ];
        let radial_sq = radial.iter().fold(
            Ball::from_rational(&Rational::zero(), BITS),
            |sum, value| sum + value * value,
        );
        // |rho-r| <= |rho^2-r^2|/r is a physical distance bound at this sample.
        let distance = &magnitude(&(&radial_sq - &radius_sq)) / &radius_lower;
        let plane_distance = (0..3).fold(Ball::from_rational(&Rational::zero(), BITS), |sum, j| {
            sum + &(&p[j] - &origin[j]) * &normal_box[j]
        });
        for (row, error) in samples.iter_mut().zip([
            magnitude(&plane_distance),
            distance,
            difference(&p, &on_plane),
            difference(&p, &on_cylinder),
        ]) {
            row.push(NumericalSample {
                parameters: vec![t.clone()],
                distance_upper_mm: error,
            });
        }
    }
    let mut checks = vec![NumericalCheck::Samples {
        name: "hole.cap.source-axial-coordinates".into(),
        scope: key.into(),
        samples: control_samples,
    }];
    checks.extend(
        [
            "hole.cap.original-plane",
            "hole.cap.original-cylinder",
            "hole.cap.plane-uv-composition",
            "hole.cap.cylinder-uv-composition",
        ]
        .into_iter()
        .zip(samples)
        .map(|(name, samples)| NumericalCheck::Samples {
            name: name.into(),
            scope: key.into(),
            samples,
        }),
    );
    NumericalEvidence::new(options.tolerance_mm.clone(), checks, Vec::new())
}
