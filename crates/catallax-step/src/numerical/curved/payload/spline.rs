//! Native spline definitions with the actual written coordinates and parameters.

use super::*;

fn knots(
    source: &[Rational],
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<Vec<Rational>, Refusal> {
    let values = source
        .iter()
        .map(|q| parameter(q, errors, options, dec))
        .collect::<Result<Vec<_>, _>>()?;
    for (old, new) in source.windows(2).zip(values.windows(2)) {
        if old[0] < old[1] && new[0] >= new[1] {
            return Err(bad("rounding-topology", "distinct spline knots collapsed"));
        }
    }
    Ok(values)
}

fn weight(
    value: &Real,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let result = rounded(
        value,
        &mut errors.spline.max_weight_error,
        options,
        false,
        dec,
    )?;
    if result.sign() != Sign::Positive {
        return Err(bad("rounding-topology", "positive spline weight vanished"));
    }
    Ok(result.into())
}

pub(super) fn curve(
    source: &BSplineCurve,
    parameters: &[Rational; 2],
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let source = finite_curve(source, parameters, false, dec)?;
    let knots = knots(source.knots(), errors, options, dec)?;
    let controls = source
        .controls()
        .map(|p| point(p, errors, options, dec).map(|p| point_from(&p)))
        .collect::<Result<Vec<_>, _>>()?;
    let weights = source
        .weights()
        .map(|w| weight(w, errors, options, dec))
        .collect::<Result<Vec<_>, _>>()?;
    BSplineCurve::new(source.degree(), knots, controls, weights, dec)
}

fn finite_curve(
    source: &BSplineCurve,
    parameters: &[Rational; 2],
    uv: bool,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let mut interval = parameters.clone();
    interval.sort();
    if source.domain() == interval.each_ref() {
        return Ok(source.clone());
    }
    let degree = source.degree();
    let mut knots = Vec::new();
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    for (span, domain) in source.spans() {
        dec.charge_work(1)?;
        let lo = domain[0].clone().max(interval[0].clone());
        let hi = domain[1].clone().min(interval[1].clone());
        if lo >= hi {
            continue;
        }
        let net = source.span_bezier_controls_on(span, [&lo, &hi], dec)?;
        let first = controls.is_empty();
        knots.extend(std::iter::repeat_n(
            lo,
            if first { degree + 1 } else { degree },
        ));
        for [x, y, z, w] in net.into_iter().skip(usize::from(!first)) {
            dec.charge_work(1)?;
            controls.push(Point3::new(
                &x / &w,
                &y / &w,
                if uv { 0.into() } else { &z / &w },
            ));
            weights.push(w);
        }
    }
    if controls.is_empty() {
        return Err(bad(
            "rounding-topology",
            "finite UV interval has no retained span",
        ));
    }
    knots.extend(std::iter::repeat_n(interval[1].clone(), degree + 1));
    BSplineCurve::new(degree, knots, controls, weights, dec)
}

pub(super) fn uv_curve(
    source: &BSplineCurve,
    parameters: &[Rational; 2],
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let source = finite_curve(source, parameters, true, dec)?;
    let knots = knots(source.knots(), errors, options, dec)?;
    let controls = source
        .controls()
        .map(|p| {
            let x = rounded(
                &p.x,
                &mut errors.spline.max_parameter_error,
                options,
                false,
                dec,
            )?;
            let y = rounded(
                &p.y,
                &mut errors.spline.max_parameter_error,
                options,
                false,
                dec,
            )?;
            Ok(Point3::new(x.into(), y.into(), 0.into()))
        })
        .collect::<Result<Vec<_>, Refusal>>()?;
    let weights = source
        .weights()
        .map(|w| weight(w, errors, options, dec))
        .collect::<Result<Vec<_>, _>>()?;
    BSplineCurve::new(source.degree(), knots, controls, weights, dec)
}

pub(super) fn surface(
    source: &BSplineSurface,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let knots = source
        .knots()
        .map(|row| knots(row, errors, options, dec))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .expect("two knot axes");
    let controls = source
        .controls()
        .iter()
        .map(|p| point(p, errors, options, dec).map(|p| point_from(&p)))
        .collect::<Result<Vec<_>, _>>()?;
    let weights = source
        .weights()
        .iter()
        .map(|w| weight(w, errors, options, dec))
        .collect::<Result<Vec<_>, _>>()?;
    BSplineSurface::new(
        source.degrees(),
        knots,
        source.control_count(),
        controls,
        weights,
        dec,
    )
}

fn rational_homogeneous(source: &BSplineSurface, index: usize) -> Option<[Rational; 4]> {
    let p = source.controls()[index].as_rational()?;
    let w = source.weights()[index].as_rational()?;
    Some([&p.x * &w, &p.y * &w, &p.z * &w, w])
}

/// Encode exact C1 joins in the knot multiplicity before decimal rounding.
/// At a degree-fold knot, inverse knot insertion removes just the redundant
/// pole strip. The homogeneous equality is exact; non-removable joins stay.
fn continuity_basis<'a>(
    source: &'a BSplineSurface,
    dec: &mut Decider,
) -> Result<std::borrow::Cow<'a, BSplineSurface>, Refusal> {
    let mut carrier = std::borrow::Cow::Borrowed(source);
    for axis in 0..2 {
        let degree = carrier.degrees()[axis];
        if degree < 2 {
            continue;
        }
        let counts = carrier.control_count();
        let knots = carrier.knots()[axis];
        let mut dropped = Vec::new();
        let mut begin = degree + 1;
        while begin < counts[axis] {
            dec.charge_work(1)?;
            let mut end = begin + 1;
            while end < knots.len() && knots[end] == knots[begin] {
                end += 1;
            }
            if knots[begin] > knots[degree] && end - begin == degree {
                let strip = begin - 1;
                let alpha = &(&knots[begin] - &knots[strip]) / &(&knots[end] - &knots[strip]);
                let beta = &Rational::one() - &alpha;
                let mut redundant = true;
                for cross in 0..counts[1 - axis] {
                    let indices = [strip - 1, strip, strip + 1].map(|at| {
                        if axis == 0 {
                            at * counts[1] + cross
                        } else {
                            cross * counts[1] + at
                        }
                    });
                    let [Some(a), Some(b), Some(c)] =
                        indices.map(|i| rational_homogeneous(&carrier, i))
                    else {
                        redundant = false;
                        break;
                    };
                    for k in 0..4 {
                        dec.charge_work(1)?;
                        if b[k] != &(&beta * &a[k]) + &(&alpha * &c[k]) {
                            redundant = false;
                            break;
                        }
                    }
                    if !redundant {
                        break;
                    }
                }
                if redundant {
                    dropped.push((end - 1, strip));
                }
            }
            begin = end;
        }
        if dropped.is_empty() {
            continue;
        }
        let output_knots = std::array::from_fn(|other| {
            let row = carrier.knots()[other];
            if other != axis {
                return row.to_vec();
            }
            row.iter()
                .enumerate()
                .filter(|(i, _)| dropped.binary_search_by_key(i, |d| d.0).is_err())
                .map(|(_, q)| q.clone())
                .collect()
        });
        let retained = |index: &usize| {
            let strip = if axis == 0 {
                index / counts[1]
            } else {
                index % counts[1]
            };
            dropped.binary_search_by_key(&strip, |d| d.1).is_err()
        };
        let controls = carrier
            .controls()
            .iter()
            .enumerate()
            .filter(|(i, _)| retained(i))
            .map(|(_, p)| p.clone())
            .collect();
        let weights = carrier
            .weights()
            .iter()
            .enumerate()
            .filter(|(i, _)| retained(i))
            .map(|(_, w)| w.clone())
            .collect();
        let mut output_counts = counts;
        output_counts[axis] -= dropped.len();
        carrier = std::borrow::Cow::Owned(BSplineSurface::new(
            carrier.degrees(),
            output_knots,
            output_counts,
            controls,
            weights,
            dec,
        )?);
    }
    Ok(carrier)
}

/// The written atlas face: the rounded source, partitioned again, keeping
/// the component whose charts have the original domains after rounding.
pub(super) fn atlas(
    source: &RationalAtlas,
    errors: &mut Errors,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<RationalAtlas, Refusal> {
    let basis = continuity_basis(source.source(), dec)?;
    let carrier = surface(&basis, errors, options, dec)?;
    let mut wanted = Vec::with_capacity(source.charts().len());
    for chart in source.charts() {
        let mut domain = chart.patch().domain().clone();
        for value in domain.iter_mut().flatten() {
            dec.charge_work(1)?;
            *value = written(value, options)?;
        }
        wanted.push(domain);
    }
    let periodic = [source.periodic_u(), source.periodic_v()];
    for candidate in RationalAtlas::partition(carrier, periodic, dec)? {
        dec.charge_work(1)?;
        if candidate
            .charts()
            .iter()
            .map(|chart| chart.patch().domain())
            .eq(wanted.iter())
        {
            return Ok(candidate);
        }
    }
    Err(bad(
        "rounding-topology",
        "the written atlas source no longer has the original charts",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    #[test]
    fn decimal_export_keeps_an_exact_c1_join_as_one_face() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let controls = (0..2)
                .flat_map(|x| {
                    (0..5).map(move |j| {
                        Point3::new(
                            x.into(),
                            Rational::from_ratio(j as i64, 3).unwrap().into(),
                            Rational::from_ratio([0, 1, 0, -1, 0][j], 3).unwrap().into(),
                        )
                    })
                })
                .collect();
            let source = BSplineSurface::new(
                [1, 2],
                [
                    [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
                    [0, 0, 0, 1, 1, 2, 2, 2].map(Rational::from_i64).to_vec(),
                ],
                [2, 5],
                controls,
                vec![Real::one(); 10],
                &mut dec,
            )
            .unwrap();
            let [source]: [_; 1] = RationalAtlas::partition(source, [false; 2], &mut dec)
                .unwrap()
                .try_into()
                .unwrap();
            let written = atlas(
                &source,
                &mut Errors::default(),
                &NumericalStepOptions::default(),
                &mut dec,
            )
            .unwrap();
            let tolerance = Rational::from_ratio(1, 10_000_000).unwrap();
            for u in [Rational::zero(), Rational::one()] {
                for i in 0..=8 {
                    let v = Rational::from_ratio(i, 4).unwrap();
                    let a = source
                        .source()
                        .evaluate_box([&u, &u], [&v, &v], 128, &dec)
                        .unwrap();
                    let b = written
                        .source()
                        .evaluate_box([&u, &u], [&v, &v], 128, &dec)
                        .unwrap();
                    for k in 0..3 {
                        let error = &a[k] - &b[k];
                        assert!(error.lo().abs() <= tolerance);
                        assert!(error.hi().abs() <= tolerance);
                    }
                }
            }
        });
        assert!(!work.exhausted, "{work:?}");
    }

    #[test]
    fn a_finite_uv_basis_keeps_native_parameters_across_partial_spans() {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let source = BSplineCurve::new(
                2,
                [0, 0, 0, 1, 1, 2, 2, 2].map(Rational::from_i64).to_vec(),
                [[0, 0], [1, 1], [2, -1], [3, 2], [4, 0]]
                    .map(|[x, y]| Point3::new(x.into(), y.into(), Real::from(2).sqrt()))
                    .to_vec(),
                [1, 2, 3, 2, 1].map(Real::from).to_vec(),
                &mut dec,
            )
            .unwrap();
            let interval = [
                Rational::from_ratio(1, 4).unwrap(),
                Rational::from_ratio(7, 4).unwrap(),
            ];
            for reversed in [false, true] {
                let mut authored = interval.clone();
                if reversed {
                    authored.reverse();
                }
                let finite = finite_curve(&source, &authored, true, &mut dec).unwrap();
                assert_eq!(finite.domain(), interval.each_ref());
                assert_eq!(finite.degree(), source.degree());
                assert_eq!(
                    finite
                        .knots()
                        .iter()
                        .filter(|t| **t == Rational::one())
                        .count(),
                    2
                );
                for i in 0..=8 {
                    let share = Rational::from_ratio(i, 8).unwrap();
                    let t = Real::from(&interval[0] + &(&share * &(&interval[1] - &interval[0])));
                    let a = source.point_at(&t, &mut dec).unwrap();
                    let b = finite.point_at(&t, &mut dec).unwrap();
                    assert_eq!(a.x.as_rational(), b.x.as_rational());
                    assert_eq!(a.y.as_rational(), b.y.as_rational());
                    assert_eq!(b.z.as_rational(), Some(Rational::zero()));
                }
                assert!(finite.point_at(&Real::zero(), &mut dec).is_err());
                assert!(finite.point_at(&Real::from(2), &mut dec).is_err());
            }
        });
        assert!(!work.exhausted, "{work:?}");
    }
}
