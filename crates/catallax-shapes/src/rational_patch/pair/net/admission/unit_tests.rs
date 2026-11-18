use super::*;
use crate::{BSplineSurface, Point3, RationalAtlas, RationalPatch};
use catallax_numbers::with_operation_budget;

fn q(n: i64, d: i64) -> Rational {
    Rational::from_ratio(n, d).unwrap()
}

fn equivalent(
    source: &BSplineSurface,
    span: [usize; 2],
    dec: &mut Decider,
) -> (PolynomialBounds, PolynomialBounds) {
    let polynomials = source.span_polynomials(span, dec).unwrap();
    let coefficients = polynomials
        .each_ref()
        .map(|polynomial| polynomial.as_rational().unwrap());
    let knots = source.knots();
    let domain = std::array::from_fn(|axis| {
        [
            knots[axis][span[axis]].clone(),
            knots[axis][span[axis] + 1].clone(),
        ]
    });
    let old = PolynomialBounds::new(&coefficients, domain, dec).unwrap();
    let new = PolynomialBounds::from_span(source, span, dec).unwrap();
    assert_eq!(old.domain, new.domain);
    assert_eq!(old.widths, new.widths);
    assert_eq!(old.weight.deg, new.weight.deg);
    assert_eq!(old.weight.v, new.weight.v);
    for axis in 0..2 {
        for coordinate in 0..3 {
            assert_eq!(
                old.tangents[axis][coordinate].deg,
                new.tangents[axis][coordinate].deg
            );
            assert_eq!(
                old.tangents[axis][coordinate].v,
                new.tangents[axis][coordinate].v
            );
        }
    }
    (old, new)
}

fn compare_ranges(
    reference: &PolynomialBounds,
    candidate: &PolynomialBounds,
    cell: &Domain,
    dec: &mut Decider,
) {
    fn grids(bounds: &PolynomialBounds) -> [&Grid; 7] {
        [
            &bounds.weight,
            &bounds.tangents[0][0],
            &bounds.tangents[0][1],
            &bounds.tangents[0][2],
            &bounds.tangents[1][0],
            &bounds.tangents[1][1],
            &bounds.tangents[1][2],
        ]
    }
    let cuts = restriction(&reference.domain, &reference.widths, cell, dec).unwrap();
    let expected = grids(reference).map(|grid| restricted_range(grid, &cuts, dec).unwrap());
    let mut cache = candidate.cache.borrow_mut();
    let prepared = cache.prepare(&candidate.domain, cell, cuts, dec).unwrap();
    for (component, (grid, expected)) in grids(candidate).into_iter().zip(expected).enumerate() {
        assert_eq!(
            cache.range(grid, component, &prepared, dec).unwrap(),
            expected
        );
    }
}

#[test]
fn unit_span_bounds_preserve_derivatives_on_translated_tiny_parameter_domains() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let shift = q(1_000_000_000_000, 1);
        for domain in [
            [[q(2, 1), q(4, 1)], [q(-3, 1), q(1, 1)]],
            [
                [shift.clone(), &shift + &Rational::one().mul_pow2(-20)],
                [q(-3, 1), q(1, 1)],
            ],
        ] {
            let mut dec = Decider::default();
            let controls = [0, 1]
                .into_iter()
                .flat_map(|a| {
                    [-3, 1].map(move |v| {
                        let w = q(1 + a, 1);
                        let x = &(&q(1000, 1) * &w) + &q(a, 1);
                        let y = &w * &q(v, 1);
                        [x.into(), y.into(), Real::zero(), w.into()]
                    })
                })
                .collect();
            let source = BSplineSurface::from_homogeneous(
                [1, 1],
                domain.each_ref().map(|axis| {
                    vec![
                        axis[0].clone(),
                        axis[0].clone(),
                        axis[1].clone(),
                        axis[1].clone(),
                    ]
                }),
                [2, 2],
                controls,
                &mut dec,
            )
            .unwrap();
            let (old, new) = equivalent(&source, [1, 1], &mut dec);
            let width = &domain[0][1] - &domain[0][0];
            let low = &q(1, 4) / &width;
            let high = &Rational::one() / &width;
            for certificate in [&old, &new] {
                let [du, dv] = certificate.derivatives(&domain, 128, &mut dec).unwrap();
                assert!(du[0].contains(&low) && du[0].contains(&high));
                assert!(du[0].lo() > Rational::zero());
                assert!(du[1].is_exact_zero() && du[2].is_exact_zero());
                assert!(dv[0].is_exact_zero() && dv[2].is_exact_zero());
                assert!(dv[1].contains(&Rational::one()));
                RationalPatch::new_cell_with_polynomial_bounds(
                    source.clone(),
                    certificate,
                    &mut dec,
                )
                .unwrap();
            }
            let initial = [
                [[q(0, 1), q(1, 2)], [q(0, 1), q(1, 1)]],
                [[q(1, 2), q(1, 1)], [q(0, 1), q(1, 1)]],
                [[q(1, 4), q(3, 4)], [q(0, 1), q(1, 1)]],
                [[q(1, 3), q(2, 3)], [q(1, 5), q(4, 5)]],
                [[q(0, 1), q(1, 4)], [q(0, 1), q(1, 2)]],
                [[q(1, 4), q(1, 2)], [q(1, 2), q(1, 1)]],
                [[q(0, 1), q(1, 1)], [q(0, 1), q(1, 1)]],
            ];
            let strips = (0..64).map(|i| [[q(i, 64), q(i + 1, 64)], [q(0, 1), q(1, 1)]]);
            for parts in initial.into_iter().chain(strips).chain([
                [[q(1, 3), q(2, 3)], [q(1, 5), q(4, 5)]],
                [[q(1, 4), q(3, 4)], [q(1, 4), q(3, 4)]],
            ]) {
                let cell = std::array::from_fn(|axis| {
                    let width = &domain[axis][1] - &domain[axis][0];
                    parts[axis]
                        .each_ref()
                        .map(|share| &domain[axis][0] + &(&width * share))
                });
                compare_ranges(&old, &new, &cell, &mut dec);
                let [du, dv] = new.derivatives(&cell, 128, &mut dec).unwrap();
                for end in &parts[0] {
                    let weight = &Rational::one() + end;
                    let expected = &Rational::one() / &(&width * &(&weight * &weight));
                    assert!(du[0].contains(&expected));
                }
                assert!(du[0].lo() > Rational::zero());
                assert!(du[1].is_exact_zero() && du[2].is_exact_zero());
                assert!(dv[0].is_exact_zero() && dv[2].is_exact_zero());
                assert!(dv[1].contains(&Rational::one()));
            }
        }
    });
    assert!(!work.exhausted);
}

#[test]
fn unit_span_bounds_preserve_split_half_cylinder_admission() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let controls = [(1, 0, 1), (1, 1, 1), (0, 2, 2), (-2, 2, 2), (-4, 0, 4)]
            .into_iter()
            .flat_map(|(x, y, w)| {
                [0, 20].map(move |z| {
                    [
                        Real::from_i64(10 * x),
                        Real::from_i64(10 * y),
                        Real::from_i64(z * w),
                        Real::from_i64(w),
                    ]
                })
            })
            .collect();
        let source = BSplineSurface::from_homogeneous(
            [2, 1],
            [
                vec![
                    q(0, 1),
                    q(0, 1),
                    q(0, 1),
                    q(1, 2),
                    q(1, 2),
                    q(1, 1),
                    q(1, 1),
                    q(1, 1),
                ],
                vec![q(0, 1), q(0, 1), q(1, 1), q(1, 1)],
            ],
            [5, 2],
            controls,
            &mut dec,
        )
        .unwrap();
        for span in source.span_indices() {
            equivalent(&source, span, &mut dec);
        }
        let atlas = RationalAtlas::partition_with_polynomial_bounds(
            source.clone(),
            [false, false],
            &mut dec,
        )
        .unwrap();
        assert_eq!(atlas.len(), 1);
        for i in 0..=16 {
            let u = q(i, 16);
            let v = q(1, 3);
            let point = source.evaluate_box([&u, &u], [&v, &v], 128, &dec).unwrap();
            let radius = (point[0].square() + point[1].square()).sqrt().unwrap();
            assert!(radius.contains(&q(10, 1)));
            assert!(point[1].lo() >= Rational::zero());
            assert!(point[2].contains(&q(20, 3)));
        }
    });
    assert!(!work.exhausted);
}

#[test]
fn unit_span_bounds_keep_fold_and_nonrational_refusals() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let source = BSplineSurface::new(
            [2, 1],
            [
                [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec(),
                [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
            ],
            [3, 2],
            [1, -3, 9]
                .into_iter()
                .flat_map(|x| {
                    [0, 1]
                        .map(move |y| Point3::new(q(x, 16).into(), Real::from_i64(y), Real::zero()))
                })
                .collect(),
            vec![Real::one(); 6],
            &mut dec,
        )
        .unwrap();
        let (old, new) = equivalent(&source, [2, 1], &mut dec);
        for certificate in [&old, &new] {
            let error = RationalPatch::new_cell_with_polynomial_bounds(
                source.clone(),
                certificate,
                &mut dec,
            )
            .expect_err("an interior fold must refuse");
            assert_eq!(error.code, "shapes.rational-patch-unproven");
        }
        let source = BSplineSurface::new(
            [1, 1],
            std::array::from_fn(|_| [0, 0, 1, 1].map(Rational::from_i64).to_vec()),
            [2, 2],
            [0, 1]
                .into_iter()
                .flat_map(|x| {
                    [0, 1].map(move |y| {
                        Point3::new(Real::from_i64(x), Real::from_i64(y), Real::zero())
                    })
                })
                .collect(),
            vec![
                Real::one(),
                Real::from_i64(2).sqrt(),
                Real::one(),
                Real::one(),
            ],
            &mut dec,
        )
        .unwrap();
        let Err(error) = PolynomialBounds::from_span(&source, [1, 1], &mut dec) else {
            panic!("irrational coefficients must not be rounded into a certificate");
        };
        assert_eq!(error.code, "shapes.rational-atlas-polynomial-coefficients");
    });
    assert!(!work.exhausted);
}

#[test]
fn degree_elevated_linear_spans_keep_the_actual_polynomial_degree() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let source = BSplineSurface::new(
            [67, 1],
            [
                [vec![Rational::zero(); 68], vec![Rational::one(); 68]].concat(),
                [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
            ],
            [68, 2],
            (0..=67)
                .flat_map(|x| {
                    [0, 1]
                        .map(move |y| Point3::new(q(x, 67).into(), Real::from_i64(y), Real::zero()))
                })
                .collect(),
            vec![Real::one(); 136],
            &mut dec,
        )
        .unwrap();
        let (old, new) = equivalent(&source, [67, 1], &mut dec);
        let domain = [
            [Rational::zero(), Rational::one()],
            [Rational::zero(), Rational::one()],
        ];
        for certificate in [&old, &new] {
            let [du, dv] = certificate.derivatives(&domain, 128, &mut dec).unwrap();
            assert!(du[0].contains(&Rational::one()) && dv[1].contains(&Rational::one()));
            assert!(du[1].is_exact_zero() && du[2].is_exact_zero());
            assert!(dv[0].is_exact_zero() && dv[2].is_exact_zero());
        }
    });
    assert!(!work.exhausted);
}
