//! Complete finite spline-line meetings and exact support nonidentity.

use crate::spline_arithmetic::{arithmetic, cut, difference, quotient, value};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{BSplineCurve, Curve, CurveId, GeomStore, Line, Point3};

fn missing(detail: &str) -> Refusal {
    Refusal::new("intersect.spline-line", detail)
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    crate::spline_arithmetic::sign(value, dec, "intersect.spline-line")
}

fn pair(geom: &GeomStore, a: CurveId, b: CurveId) -> Option<(&BSplineCurve, &Line)> {
    match (geom.curve(a), geom.curve(b)) {
        (Curve::BSpline(curve), Curve::Line(line)) | (Curve::Line(line), Curve::BSpline(curve)) => {
            Some((curve, line))
        }
        _ => None,
    }
}

// For nonzero d[k], the two equations
// d[k](X[j] - O[j]W) - d[j](X[k] - O[k]W) = 0 characterize the line.
fn constraints(line: &Line, dec: &mut Decider) -> Result<[[Real; 4]; 2], Refusal> {
    let origin = [&line.origin.x, &line.origin.y, &line.origin.z];
    let direction = [&line.dir.x, &line.dir.y, &line.dir.z];
    if origin
        .iter()
        .chain(direction.iter())
        .any(|value| value.rational_ref().is_none())
    {
        return Err(missing(
            "finite spline-line queries need rational line coefficients",
        ));
    }
    let axis = direction
        .iter()
        .position(|value| !value.rational_ref().unwrap().is_zero())
        .ok_or_else(|| missing("a line direction must be nonzero"))?;
    let mut result = Vec::new();
    for other in (0..3).filter(|&other| other != axis) {
        let mut coefficients = std::array::from_fn(|_| Real::zero());
        coefficients[other] = direction[axis].clone();
        coefficients[axis] = -direction[other];
        let first = arithmetic(direction[other], origin[axis], true, dec)?;
        let second = arithmetic(direction[axis], origin[other], true, dec)?;
        coefficients[3] = difference(&first, &second, dec)?;
        result.push(coefficients);
    }
    Ok(result.try_into().expect("two transverse line constraints"))
}

pub(super) fn different(
    geom: &GeomStore,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let Some((curve, line)) = pair(geom, a, b) else {
        return Ok(false);
    };
    let constraints = constraints(line, dec)?;
    for (span, _) in curve.spans() {
        let coordinates = curve.query_span_polynomials(span, dec)?;
        for constraint in &constraints {
            let equation = cut(&coordinates, constraint, dec)?;
            for (_, coefficient) in equation.terms() {
                if sign(&coefficient, dec)? != Sign::Zero {
                    // W is positive and the span has an interval with interior.
                    // A nonzero numerator cannot vanish throughout that interval.
                    return Ok(true);
                }
            }
        }
    }
    // Containment in an infinite line does not prove equality of finite sets.
    Ok(false)
}

fn same_position(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    for (a, b) in [(&a.x, &b.x), (&a.y, &b.y), (&a.z, &b.z)] {
        if sign(&difference(a, b, dec)?, dec)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Complete isolated meetings of a finite native B-spline and a rational line,
/// in either input order. Other curve pairs return None. An overlap interval
/// refuses; it is never replaced by a finite selection of points.
/// No incidence facts are published by this query.
pub fn native_spline_line_points(
    geom: &GeomStore,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<Vec<Point3>>, Refusal> {
    let Some((curve, line)) = pair(geom, a, b) else {
        return Ok(None);
    };
    let constraints = constraints(line, dec)?;
    let mut points = Vec::new();
    for (span, [lo, hi]) in curve.spans() {
        let coordinates = curve.query_span_polynomials(span, dec)?;
        let equations = [
            cut(&coordinates, &constraints[0], dec)?,
            cut(&coordinates, &constraints[1], dec)?,
        ];
        let mut selected = None;
        for (index, equation) in equations.iter().enumerate() {
            let rational = equation.rational().ok_or_else(|| {
                missing("finite spline-line meetings need rational homogeneous cut coefficients")
            })?;
            if !rational.is_zero() && selected.is_none() {
                selected = Some(index);
            }
        }
        let Some(selected) = selected else {
            return Err(Refusal::new(
                "intersect.spline-line-overlap",
                format!(
                    "knot span {span} [{lo}, {hi}] lies on the line; finite overlap endpoints and ownership need a proof"
                ),
            ));
        };
        for parameter in equations[selected].real_roots(dec)? {
            if sign(
                &difference(&parameter, &Real::from_rational(lo.clone()), dec)?,
                dec,
            )? == Sign::Negative
                || sign(
                    &difference(&Real::from_rational(hi.clone()), &parameter, dec)?,
                    dec,
                )? == Sign::Negative
            {
                continue;
            }
            let mut on_line = true;
            for equation in &equations {
                on_line &= sign(&value(equation, &parameter, false, dec)?, dec)? == Sign::Zero;
            }
            if !on_line {
                continue;
            }
            let mut position = Vec::new();
            for coordinate in &coordinates {
                position.push(value(coordinate, &parameter, false, dec)?);
            }
            // Admitted splines are continuous at internal knots. Positive
            // weights prove W nonzero on every retained closed span.
            let point = Point3::new(
                quotient(&position[0], &position[3], dec)?,
                quotient(&position[1], &position[3], dec)?,
                quotient(&position[2], &position[3], dec)?,
            );
            let mut duplicate = false;
            for previous in &points {
                if same_position(&point, previous, dec)? {
                    duplicate = true;
                    break;
                }
            }
            if !duplicate {
                points.push(point);
            }
        }
    }
    Ok(Some(points))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Outcome, Rational, with_operation_budget};

    fn point(x: i64, y: i64, z: i64) -> Point3 {
        Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
    }

    fn roof(weight: i64, dec: &mut Decider) -> (GeomStore, CurveId) {
        let mut geom = GeomStore::new();
        let curve = BSplineCurve::new(
            2,
            [0, 0, 0, 1, 1, 2, 2, 2].map(Rational::from_i64).to_vec(),
            vec![
                point(0, 0, 0),
                point(1, 2, 0),
                point(2, 2, 0),
                point(3, 2, 0),
                point(4, 0, 0),
            ],
            vec![
                Real::one(),
                Real::from_i64(weight),
                Real::one(),
                Real::from_i64(weight),
                Real::one(),
            ],
            dec,
        )
        .unwrap();
        let curve = geom.add_curve(Curve::BSpline(curve));
        (geom, curve)
    }

    fn line(geom: &mut GeomStore, y: i64, z: i64) -> CurveId {
        geom.add_curve(Curve::Line(Line {
            origin: point(0, y, z),
            dir: catallax_shapes::Vector3::new(Real::one(), Real::zero(), Real::zero()),
        }))
    }

    #[test]
    fn equal_endpoints_do_not_make_a_quadratic_spline_the_floor_line() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let (mut geom, curve) = roof(1, &mut dec);
            let line = line(&mut geom, 0, 0);
            for (a, b) in [(curve, line), (line, curve)] {
                assert!(matches!(
                    crate::same_curve(&geom, a, b, &mut dec),
                    Outcome::Proven(_)
                ));
                let points = native_spline_line_points(&geom, a, b, &mut dec)
                    .unwrap()
                    .unwrap();
                assert_eq!(points.len(), 2);
                for point in points {
                    assert_eq!(sign(&point.y, &mut dec).unwrap(), Sign::Zero);
                    assert!(
                        sign(&point.x, &mut dec).unwrap() == Sign::Zero
                            || sign(&(&point.x - &Real::from_i64(4)), &mut dec).unwrap()
                                == Sign::Zero
                    );
                }
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn complete_meetings_cover_roots_knots_tangency_and_the_other_line_constraint() {
        for weight in [1, 2] {
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let (mut geom, curve) = roof(weight, &mut dec);
                for (height, z, count) in [(1, 0, 2), (2, 0, 1), (1, 1, 0)] {
                    let line = line(&mut geom, height, z);
                    for (a, b) in [(curve, line), (line, curve)] {
                        let points = native_spline_line_points(&geom, a, b, &mut dec)
                            .unwrap()
                            .unwrap();
                        assert_eq!(points.len(), count);
                        for point in points {
                            assert_eq!(
                                sign(&(&point.y - &Real::from_i64(height)), &mut dec).unwrap(),
                                Sign::Zero
                            );
                            assert_eq!(
                                sign(&(&point.z - &Real::from_i64(z)), &mut dec).unwrap(),
                                Sign::Zero
                            );
                            if weight == 1 && height == 1 {
                                assert_eq!(
                                    sign(
                                        &((&point.x - &Real::from_i64(2)).square()
                                            - Real::from_i64(2)),
                                        &mut dec
                                    )
                                    .unwrap(),
                                    Sign::Zero
                                );
                            }
                        }
                    }
                }
            });
            assert!(!work.exhausted, "{} units", work.used);
        }
    }

    #[test]
    fn line_constraints_cover_each_pivot_axis_and_skew_directions() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let half = Real::from_rational(Rational::from_ratio(1, 2).unwrap());
            let quarter = Real::from_rational(Rational::from_ratio(1, 4).unwrap());
            let spline = BSplineCurve::new(
                2,
                [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec(),
                vec![
                    point(0, 0, 0),
                    Point3::new(half.clone(), Real::zero(), Real::zero()),
                    point(1, 1, 0),
                ],
                vec![Real::one(); 3],
                &mut dec,
            )
            .unwrap();
            let mut geom = GeomStore::new();
            let curve = geom.add_curve(Curve::BSpline(spline));
            for [x, y, z] in [[0, 2, 0], [0, 0, -3], [1, 2, 3]] {
                let line = geom.add_curve(Curve::Line(Line {
                    origin: Point3::new(half.clone(), quarter.clone(), Real::zero()),
                    dir: catallax_shapes::Vector3::new(
                        Real::from_i64(x),
                        Real::from_i64(y),
                        Real::from_i64(z),
                    ),
                }));
                let points = native_spline_line_points(&geom, curve, line, &mut dec)
                    .unwrap()
                    .unwrap();
                assert_eq!(points.len(), 1);
                let point = &points[0];
                assert_eq!(sign(&(&point.x - &half), &mut dec).unwrap(), Sign::Zero);
                assert_eq!(sign(&(&point.y - &quarter), &mut dec).unwrap(), Sign::Zero);
                assert_eq!(sign(&point.z, &mut dec).unwrap(), Sign::Zero);
            }
        });
        assert!(!work.exhausted, "{} units", work.used);
    }

    #[test]
    fn overlap_remains_unproved_and_zero_budget_cannot_skip_its_arithmetic() {
        let mut dec = Decider::default();
        for controls in [
            vec![point(0, 0, 0), point(1, 0, 0)],
            vec![point(0, 0, 0), point(1, 0, 0), point(2, 1, 0)],
        ] {
            let knots = if controls.len() == 2 {
                vec![0, 0, 1, 1]
            } else {
                vec![0, 0, 1, 2, 2]
            };
            let weights = vec![Real::one(); controls.len()];
            let curve = BSplineCurve::new(
                1,
                knots.into_iter().map(Rational::from_i64).collect(),
                controls,
                weights,
                &mut dec,
            )
            .unwrap();
            let mut geom = GeomStore::new();
            let curve = geom.add_curve(Curve::BSpline(curve));
            let line = line(&mut geom, 0, 0);
            let (result, work) =
                with_operation_budget(0, || different(&geom, curve, line, &mut dec));
            assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
            assert!(work.exhausted);
            let (result, work) = with_operation_budget(0, || {
                native_spline_line_points(&geom, curve, line, &mut dec)
            });
            assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
            assert!(work.exhausted);
            assert_eq!(
                native_spline_line_points(&geom, curve, line, &mut dec)
                    .unwrap_err()
                    .code,
                "intersect.spline-line-overlap"
            );
        }
    }
}
