//! Bernstein projection bounds for rational quadratic ruled charts.

use super::*;

type V3 = [Rational; 3];

fn coordinates(values: [&Real; 3]) -> Option<V3> {
    Some([
        values[0].as_rational()?,
        values[1].as_rational()?,
        values[2].as_rational()?,
    ])
}

fn dot(a: &V3, b: &V3) -> Rational {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}

fn minus(a: &V3, b: &V3) -> V3 {
    std::array::from_fn(|i| &a[i] - &b[i])
}

fn cross(a: &V3, b: &V3) -> V3 {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

struct RuledControls {
    points: Vec<V3>,
    weights: Vec<Rational>,
}

fn controls(carrier: &BSplineSurface, dec: &mut Decider) -> Result<Option<RuledControls>, Refusal> {
    if carrier.degrees() != [2, 1] || carrier.control_count() != [3, 2] {
        return Ok(None);
    }
    dec.charge_work(64)?;
    let [ku, kv] = carrier.knots();
    if ku.len() != 6
        || kv.len() != 4
        || ku[..3].iter().any(|q| q != &ku[0])
        || ku[3..].iter().any(|q| q != &ku[3])
        || kv[..2].iter().any(|q| q != &kv[0])
        || kv[2..].iter().any(|q| q != &kv[2])
    {
        return Ok(None);
    }
    let points: Option<Vec<_>> = carrier
        .controls()
        .iter()
        .map(|p| coordinates([&p.x, &p.y, &p.z]))
        .collect();
    let weights: Option<Vec<_>> = carrier.weights().iter().map(Real::as_rational).collect();
    let Some((points, weights)) = points.zip(weights) else {
        return Ok(None);
    };
    for u in 0..3 {
        dec.charge_work(1)?;
        if weights[2 * u] != weights[2 * u + 1] || weights[2 * u].sign() != Sign::Positive {
            return Ok(None);
        }
    }
    Ok(Some(RuledControls { points, weights }))
}

struct Projected {
    du: [Rational; 2],
    dv: [Rational; 2],
}

impl Projected {
    fn new(data: &RuledControls, axis: &V3, dec: &mut Decider) -> Result<Self, Refusal> {
        dec.charge_work(48)?;
        let p: [[Rational; 2]; 3] =
            std::array::from_fn(|u| std::array::from_fn(|v| dot(axis, &data.points[2 * u + v])));
        let w = [&data.weights[0], &data.weights[2], &data.weights[4]];
        let two = Rational::from_i64(2);
        let factors = [&two * &(w[0] * w[1]), w[0] * w[2], &two * &(w[1] * w[2])];
        let mut coefficients = Vec::with_capacity(6);
        for v in 0..2 {
            coefficients.extend([
                &factors[0] * &(&p[1][v] - &p[0][v]),
                &factors[1] * &(&p[2][v] - &p[0][v]),
                &factors[2] * &(&p[2][v] - &p[1][v]),
            ]);
        }
        let minimum = coefficients.iter().min().unwrap();
        let maximum = coefficients.iter().max().unwrap();
        let minimum_weight = w.iter().min().unwrap();
        let maximum_weight = w.iter().max().unwrap();
        let square_low = *minimum_weight * *minimum_weight;
        let square_high = *maximum_weight * *maximum_weight;
        let lower_denominator = if minimum.sign() == Sign::Negative {
            &square_low
        } else {
            &square_high
        };
        let upper_denominator = if maximum.sign() == Sign::Negative {
            &square_high
        } else {
            &square_low
        };
        let differences = p.each_ref().map(|row| &row[1] - &row[0]);
        Ok(Self {
            du: [minimum / lower_denominator, maximum / upper_denominator],
            dv: [
                differences.iter().min().unwrap().clone(),
                differences.iter().max().unwrap().clone(),
            ],
        })
    }
}

fn absolute_upper(interval: &[Rational; 2]) -> Rational {
    interval[0].abs().max(interval[1].abs())
}

/// Bounds the full symmetric projected Jacobian with rational Bernstein controls.
/// Positive weights repeated across v make each projected v derivative a convex
/// combination of control differences. Scaling its second row can then prove
/// a positive determinant even when the centre-dual projection is indefinite.
pub(super) fn certificate(
    carrier: &BSplineSurface,
    first: &Vector3,
    dec: &mut Decider,
) -> Result<Option<(Vector3, Rational)>, Refusal> {
    let Some(data) = controls(carrier, dec)? else {
        return Ok(None);
    };
    let Some(first) = coordinates([&first.x, &first.y, &first.z]) else {
        return Ok(None);
    };
    let normal = cross(
        &minus(&data.points[2], &data.points[0]),
        &minus(&data.points[4], &data.points[0]),
    );
    let change = dot(&normal, &minus(&data.points[1], &data.points[0]));
    if change.sign() == Sign::Zero {
        return Ok(None);
    }
    let height: V3 = std::array::from_fn(|i| &normal[i] / &change);
    let f = Projected::new(&data, &first, dec)?;
    let g = Projected::new(&data, &height, dec)?;
    if f.du[0].sign() != Sign::Positive || g.dv[0].sign() != Sign::Positive {
        return Ok(None);
    }
    let along = &f.du[0] * &g.dv[0];
    let cross = absolute_upper(&f.dv);
    let scale = Rational::one() + &(&cross * &cross) / &along;
    let off = &cross + &(&scale * &absolute_upper(&g.du));
    let determinant = &along * &scale - &(&off * &off) / &Rational::from_i64(4);
    if determinant.sign() != Sign::Positive {
        return Ok(None);
    }
    let margin = &determinant / &(&f.du[1] + &(&scale * &g.dv[1]));
    let axis = height.map(|q| Real::from_rational(&q * &scale));
    Ok(Some((
        Vector3::new(axis[0].clone(), axis[1].clone(), axis[2].clone()),
        margin,
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    fn cone(low: i64, high: i64, moved: bool, dec: &mut Decider) -> BSplineSurface {
        let r = Rational::from_i64;
        let mut controls = Vec::new();
        for (x, y) in [(1, 0), (1, 1), (0, 1)] {
            for h in [low, high] {
                let p = if moved {
                    [x * h + 2 * h + 3, y * h - 3 * h - 7, 4 * h + 2]
                } else {
                    [x * h, y * h, h]
                };
                controls.push(Point3::new(p[0].into(), p[1].into(), p[2].into()));
            }
        }
        BSplineSurface::new(
            [2, 1],
            [
                vec![r(3), r(3), r(3), r(7), r(7), r(7)],
                vec![r(-2), r(-2), r(5), r(5)],
            ],
            [3, 2],
            controls,
            [1, 1, 1, 1, 2, 2]
                .map(|q| Real::from_rational(r(q)))
                .to_vec(),
            dec,
        )
        .unwrap()
    }

    #[test]
    fn a_wide_ruled_quarter_has_a_global_chart_and_unique_finite_preimages() {
        for moved in [false, true] {
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let patch = RationalPatch::new(cone(1, 20, moved, &mut dec), &mut dec).unwrap();
                assert!(patch.embedding_lower_bound().sign() == Sign::Positive);
                assert_eq!(patch.embedding_cover().count(), 1);
                for (u, v) in [
                    (Rational::from_i64(3), Rational::from_i64(-2)),
                    (
                        Rational::from_ratio(13, 3).unwrap(),
                        Rational::from_ratio(4, 5).unwrap(),
                    ),
                    (Rational::from_i64(7), Rational::from_i64(5)),
                ] {
                    let expected = [Real::from_rational(u), Real::from_rational(v)];
                    let point = patch.point_at(&expected, &mut dec).unwrap();
                    let actual = patch.parameters_at(&point, &mut dec).unwrap().unwrap();
                    for k in 0..2 {
                        assert_eq!(
                            exact_sign(&(&actual[k] - &expected[k]), &mut dec).unwrap(),
                            Sign::Zero
                        );
                    }
                }
            });
            assert!(!work.exhausted, "moved={moved}: {work:?}");
        }
    }

    #[test]
    fn rounded_nonparallel_control_rows_keep_the_retained_chart() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let original = cone(1, 20, true, &mut dec);
            let mut points = original.controls().to_vec();
            points[3].z = &points[3].z + &Real::from(Rational::parse("1e-12").unwrap());
            let source = BSplineSurface::new(
                original.degrees(),
                original.knots().map(<[_]>::to_vec),
                original.control_count(),
                points,
                original.weights().to_vec(),
                &mut dec,
            )
            .unwrap();
            let patch = RationalPatch::new(source.clone(), &mut dec).unwrap();
            assert!(patch.carrier().same_definition(&source));
            assert_eq!(
                patch.domain(),
                &source.domain().map(|row| row.map(Clone::clone))
            );
            let at = [Real::from(5), Real::from(1)];
            let point = patch.point_at(&at, &mut dec).unwrap();
            let inverse = patch.parameters_at(&point, &mut dec).unwrap().unwrap();
            for axis in 0..2 {
                assert_eq!(
                    exact_sign(&(&at[axis] - &inverse[axis]), &mut dec).unwrap(),
                    Sign::Zero
                );
            }
        });
        assert!(!work.exhausted, "{work:?}");
    }

    #[test]
    fn a_collapsed_apex_does_not_receive_a_regular_chart() {
        let mut dec = Decider::default();
        assert!(RationalPatch::new(cone(-1, 1, false, &mut dec), &mut dec).is_err());
    }
}
