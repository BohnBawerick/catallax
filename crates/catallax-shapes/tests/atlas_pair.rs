//! Exact chart and periodic ownership through the existing cell-pair engine.

use catallax_numbers::{
    Decider, Outcome, Rational, Real, Refusal, ZeroProof, with_operation_budget,
};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Plane, Point3, RationalAtlas, SectionHits, Vector3,
};
use std::num::NonZeroUsize;

fn q(n: i64, d: i64) -> Real {
    Real::from_rational(Rational::from_ratio(n, d).unwrap())
}

fn split_vertical_plane(first: bool, dec: &mut Decider) -> RationalAtlas {
    let controls = [-1, 1]
        .into_iter()
        .flat_map(|horizontal| {
            [0, 1, 2].map(|height| {
                if first {
                    Point3::new(q(horizontal, 1), q(0, 1), q(height, 1))
                } else {
                    Point3::new(q(0, 1), q(horizontal, 1), q(height, 1))
                }
            })
        })
        .collect();
    let source = BSplineSurface::new(
        [1, 1],
        [
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
            [0, 0, 1, 2, 2].map(Rational::from_i64).to_vec(),
        ],
        [2, 3],
        controls,
        vec![Real::one(); 6],
        dec,
    )
    .unwrap();
    let [atlas] =
        <[_; 1]>::try_from(RationalAtlas::partition(source, [false, false], dec).unwrap()).unwrap();
    assert_eq!(atlas.charts().len(), 2);
    atlas
}

#[test]
fn offset_spine_crosses_both_atlas_seams_without_fragmenting() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let first = split_vertical_plane(true, &mut dec);
        let second = split_vertical_plane(false, &mut dec);
        let tolerance = Rational::parse("0.0000001").unwrap();
        let spine = catallax_shapes::NumericalFilletSpine::trace(
            &[first, second],
            [
                Rational::from_ratio(1, 10).unwrap(),
                Rational::from_ratio(1, 5).unwrap(),
            ],
            1,
            [Rational::zero(), Rational::from_i64(2)],
            [
                Rational::from_ratio(3, 5).unwrap(),
                Rational::from_ratio(9, 20).unwrap(),
                Rational::zero(),
            ],
            tolerance.clone(),
            &mut dec,
        )
        .unwrap();
        let charts: std::collections::BTreeSet<_> =
            spine.nodes().iter().map(|node| node.charts).collect();
        assert_eq!(charts, [[0, 0], [1, 1]].into_iter().collect());
        assert_eq!(spine.tolerance_mm(), &tolerance);
        assert!(spine.sampled_error_mm() <= &tolerance);
        let curve = spine.centre_curve();
        for i in 0..=32 {
            let t = Rational::from_ratio(i, 16).unwrap();
            let point = curve.evaluate_box([&t, &t], 128, &dec).unwrap();
            for (value, expected) in point.iter().zip([
                Rational::from_ratio(1, 5).unwrap(),
                Rational::from_ratio(-1, 10).unwrap(),
                t.clone(),
            ]) {
                assert!((&value.lo() - &expected).abs() <= tolerance);
                assert!((&value.hi() - &expected).abs() <= tolerance);
            }
        }
    });
    assert!(!work.exhausted, "{work:?}");
}

#[test]
fn numerical_fillet_crosses_chart_seams_with_circular_sections() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let supports = [
            split_vertical_plane(true, &mut dec),
            split_vertical_plane(false, &mut dec),
        ];
        let supports = supports.map(|support| {
            catallax_shapes::NumericalFilletSpine::extended_support(
                &support,
                1,
                &Rational::from_ratio(1, 4).unwrap(),
                &mut dec,
            )
            .unwrap()
        });
        let radius = Rational::from_ratio(1, 10).unwrap();
        let tolerance = Rational::parse("0.0000001").unwrap();
        let spine = catallax_shapes::NumericalFilletSpine::trace(
            &supports,
            [radius.clone(), radius.clone()],
            1,
            [
                Rational::from_ratio(-1, 2).unwrap(),
                Rational::from_ratio(5, 2).unwrap(),
            ],
            [
                Rational::from_ratio(11, 20).unwrap(),
                Rational::from_ratio(9, 20).unwrap(),
                Rational::from_ratio(-1, 2).unwrap(),
            ],
            tolerance.clone(),
            &mut dec,
        )
        .unwrap();
        let fillet = catallax_shapes::NumericalSplineFillet::fit(spine, &mut dec).unwrap();
        assert!(fillet.sampled_error_mm() <= &tolerance);
        let source = fillet.surface().source();
        for i in -4..=20 {
            let tau = Rational::from_ratio(i, 8).unwrap();
            for j in 0..=8 {
                let across = Rational::from_ratio(j, 8).unwrap();
                let point = source
                    .evaluate_box([&tau, &tau], [&across, &across], 128, &dec)
                    .unwrap();
                let x = &point[0] - &catallax_numbers::Ball::from_rational(&radius, 128);
                let y = &point[1] + &catallax_numbers::Ball::from_rational(&radius, 128);
                let square = &x.pow_u32(2) + &y.pow_u32(2);
                let error = (&square.lo() - &(&radius * &radius))
                    .abs()
                    .max((&square.hi() - &(&radius * &radius)).abs());
                assert!(error <= &tolerance * &radius);
                assert!((&point[2].lo() - &tau).abs() <= tolerance);
                assert!((&point[2].hi() - &tau).abs() <= tolerance);
                if j == 0 {
                    assert!(point[1].lo().abs().max(point[1].hi().abs()) <= tolerance);
                }
                if j == 8 {
                    assert!(point[0].lo().abs().max(point[0].hi().abs()) <= tolerance);
                }
            }
        }
        for normal in [[0, 0, 1], [1, 0, 1]] {
            let cut = fillet
                .plane_cut(normal.map(Rational::from_i64), Rational::zero(), &mut dec)
                .unwrap();
            for i in 0..=32 {
                let s = Rational::from_ratio(i, 32).unwrap();
                let point = cut.curve().evaluate_box([&s, &s], 128, &dec).unwrap();
                let uv = cut.parameters().evaluate_box([&s, &s], 128, &dec).unwrap();
                let tau = uv[0].mid();
                let across = uv[1].mid();
                let on_fillet = source
                    .evaluate_box([&tau, &tau], [&across, &across], 128, &dec)
                    .unwrap();
                let residual = &point[2]
                    + &(&point[0]
                        * &catallax_numbers::Ball::from_rational(
                            &Rational::from_i64(normal[0]),
                            128,
                        ));
                assert!(residual.lo().abs().max(residual.hi().abs()) <= tolerance);
                for k in 0..3 {
                    assert!((&point[k].lo() - &on_fillet[k].hi()).abs() <= tolerance);
                    assert!((&point[k].hi() - &on_fillet[k].lo()).abs() <= tolerance);
                }
            }
        }
    });
    assert!(!work.exhausted, "{work:?}");
}

#[test]
fn retained_source_clearance_keeps_native_initial_and_later_contacts() {
    for reversed in [false, true] {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let cylinder = cylinder(true, false, false, &mut dec).unwrap();
            let plane = plane(&mut dec);
            let meeting = if reversed {
                plane.patch_meeting_with_cached_spans(&cylinder, &mut dec)
            } else {
                cylinder.patch_meeting_with_cached_spans(&plane, &mut dec)
            }
            .unwrap();
            let component = meeting.curve_branch(0, &mut dec).unwrap();
            for (x, heights, allowed, expected) in [
                (1, vec![q(1, 5), q(1, 1)], true, true),
                (1, vec![q(1, 5), q(1, 1)], false, false),
                (1, vec![q(0, 1), q(1, 1)], true, false),
                (1, vec![q(0, 1), q(1, 5)], true, false),
                (1, vec![q(1, 1), q(2, 1)], false, true),
                (3, vec![q(0, 1), q(1, 1)], false, true),
                (1, vec![q(1, 5), q(6, 5), q(1, 5)], true, false),
            ] {
                let degree = heights.len() - 1;
                let knots = std::iter::repeat_n(Rational::from_i64(2), degree + 1)
                    .chain(std::iter::repeat_n(Rational::from_i64(5), degree + 1))
                    .collect();
                let path = BSplineCurve::new(
                    degree,
                    knots,
                    heights
                        .into_iter()
                        .map(|z| Point3::new(q(x, 1), Real::zero(), z))
                        .collect(),
                    vec![Real::one(); degree + 1],
                    &mut dec,
                )
                .unwrap();
                assert_eq!(
                    component
                        .clears_spline_sources(&path, allowed, &mut dec)
                        .unwrap(),
                    expected,
                    "reversed={reversed},x={x},degree={degree},allowed={allowed}"
                );
            }
        });
        assert!(!work.exhausted, "{work:?}");
    }
}

#[test]
fn retained_pair_components_keep_complete_plane_points_and_near_contacts() {
    for reversed in [false, true] {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let cylinder = cylinder(true, false, false, &mut dec).unwrap();
            let support = plane(&mut dec);
            let meeting = if reversed {
                support.patch_meeting_with_cached_spans(&cylinder, &mut dec)
            } else {
                cylinder.patch_meeting_with_cached_spans(&support, &mut dec)
            }
            .unwrap();
            assert_eq!(meeting.branches().len(), 1);
            let curve = meeting.curve_branch(0, &mut dec).unwrap();
            for (coordinate, expected) in [
                (q(0, 1), 2),
                (q(1, 2), 2),
                (q(1, 1), 1),
                (q(999_999, 1_000_000), 2),
                (q(1_000_001, 1_000_000), 0),
            ] {
                let cut = Plane::new(
                    Point3::new(coordinate.clone(), Real::zero(), Real::zero()),
                    Vector3::new(Real::one(), Real::zero(), Real::zero()),
                    Vector3::new(Real::zero(), Real::one(), Real::zero()),
                    &mut dec,
                )
                .unwrap();
                let SectionHits::Points(points) =
                    curve.plane_hits(&cut, &mut dec).unwrap_or_else(|error| {
                        panic!("reversed={reversed}, coordinate={coordinate:?}: {error:?}")
                    })
                else {
                    panic!("the circular component cannot lie in a vertical plane")
                };
                assert_eq!(points.len(), expected);
                for point in points {
                    assert!(curve.contains(&point, &mut dec).unwrap());
                    assert!(
                        matches!(dec.sign(&point.minus(&cut.origin).dot(&cut.normal)),
                        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound)
                    );
                }
            }
            let cut = Plane::new(
                Point3::new(Real::zero(), Real::zero(), q(1, 5)),
                Vector3::new(Real::zero(), Real::zero(), Real::one()),
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
                &mut dec,
            )
            .unwrap();
            assert!(matches!(
                curve.plane_hits(&cut, &mut dec).unwrap(),
                SectionHits::Contained
            ));
        });
        assert!(!work.exhausted, "{reversed}: {} work units", work.used);
    }
}

fn cylinder(
    periodic: bool,
    transpose: bool,
    gap: bool,
    dec: &mut Decider,
) -> Result<RationalAtlas, Refusal> {
    let circle = [
        (1, 0),
        (1, 1),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (-1, -1),
        (0, -1),
        (1, -1),
        (1, 0),
    ];
    let mut rows = Vec::new();
    let mut scale = Real::one();
    for (i, (x, y)) in circle.into_iter().enumerate() {
        let weight = if i % 2 == 0 {
            scale.clone()
        } else {
            &scale * &q(7, 10)
        };
        let x = if gap && i == 8 {
            &Real::from_i64(x) + &q(1, 1_000_000)
        } else {
            x.into()
        };
        rows.push((x, Real::from_i64(y), weight));
        if i % 2 == 1 {
            scale = &scale * &q(49, 50);
        }
    }
    let mut controls = Vec::new();
    let mut weights = Vec::new();
    for i in 0..18 {
        let (circle, axial) = if transpose {
            (i % 9, i / 9)
        } else {
            (i / 2, i % 2)
        };
        let (x, y, w) = &rows[circle];
        controls.push(Point3::new(
            x.clone(),
            y.clone(),
            Real::from_i64(if axial == 0 { -1 } else { 1 }),
        ));
        weights.push(w.clone());
    }
    let circle_knots = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
        .map(Rational::from_i64)
        .to_vec();
    let linear_knots = [0, 0, 1, 1].map(Rational::from_i64).to_vec();
    let (degree, count, knots, periodic) = if transpose {
        (
            [1, 2],
            [2, 9],
            [linear_knots, circle_knots],
            [false, periodic],
        )
    } else {
        (
            [2, 1],
            [9, 2],
            [circle_knots, linear_knots],
            [periodic, false],
        )
    };
    let source = BSplineSurface::new(degree, knots, count, controls, weights, dec)?;
    let mut atlases = RationalAtlas::partition(source, periodic, dec)?;
    assert_eq!(atlases.len(), 1);
    Ok(atlases.pop().unwrap())
}

fn plane(dec: &mut Decider) -> RationalAtlas {
    let source = BSplineSurface::new(
        [1, 1],
        std::array::from_fn(|_| [0, 0, 1, 1].map(Rational::from_i64).to_vec()),
        [2, 2],
        [[-2, -2], [-2, 2], [2, -2], [2, 2]]
            .map(|p| Point3::new(p[0].into(), p[1].into(), q(1, 5)))
            .to_vec(),
        vec![Real::one(); 4],
        dec,
    )
    .unwrap();
    RationalAtlas::partition(source, [false, false], dec)
        .unwrap()
        .pop()
        .unwrap()
}

fn run(periodic: bool, transpose: bool, reversed: bool, cached: bool) -> (String, u64) {
    let (result, work) = with_operation_budget(20_000_000, || -> Result<String, Refusal> {
        let mut dec = Decider::default();
        let cylinder = cylinder(periodic, transpose, false, &mut dec)?;
        let plane = plane(&mut dec);
        let (first, second) = if reversed {
            (&plane, &cylinder)
        } else {
            (&cylinder, &plane)
        };
        let found = if cached {
            first.patch_meeting_with_cached_spans(second, &mut dec)?
        } else {
            first.patch_meeting(second, &mut dec)?
        };
        assert_eq!(
            found.branches().len(),
            1,
            "chart cuts cannot fragment the circle"
        );
        let branch = &found.branches()[0];
        assert_eq!(
            branch.is_closed(),
            periodic,
            "only a proved periodic identification closes the parameter graph"
        );
        let mut summary = format!("closed={};", branch.is_closed());
        for arc in branch.arcs() {
            summary.push_str(&format!(
                "{:?};",
                (
                    arc.first_chart,
                    arc.first_span,
                    &arc.first_cell,
                    arc.second_chart,
                    arc.second_span,
                    &arc.second_cell,
                    arc.certificate
                )
            ));
            for end in &arc.ends {
                summary.push_str(&format!("{:?};", end.enclosure()));
                let point = found.point(end);
                let on_cylinder = &point.x.square() + &point.y.square() - Real::one();
                let on_plane = &point.z - &q(1, 5);
                for equation in [on_cylinder, on_plane] {
                    let bound = dec.enclose_budgeted_at(&equation, 128)?;
                    assert!(
                        bound.contains_zero(),
                        "retained endpoint leaves an original support: {bound:?}"
                    );
                    assert!(bound.hi() - bound.lo() < Rational::one().mul_pow2(-80));
                }
            }
        }
        Ok(summary)
    });
    assert!(!work.exhausted);
    (result.unwrap(), work.used)
}

#[test]
fn periodic_rows_join_exactly_in_either_parameter_and_either_pair_order() {
    let mut baseline = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let records: Vec<_> = [false, true]
            .into_iter()
            .flat_map(|transpose| {
                [false, true]
                    .into_iter()
                    .map(move |reversed| run(true, transpose, reversed, false))
            })
            .collect();
        println!(
            "atlas pair threads={threads:?} work={:?}",
            records.iter().map(|(_, work)| work).collect::<Vec<_>>()
        );
        if let Some(baseline) = &baseline {
            assert_eq!(baseline, &records);
        } else {
            baseline = Some(records);
        }
    }
    catallax_numbers::set_threads(None);
}

#[test]
fn unjoined_rows_remain_physical_even_when_their_world_points_agree() {
    for transpose in [false, true] {
        run(false, transpose, false, false);
    }
}

#[test]
fn a_nearly_closing_row_is_not_a_periodic_identification() {
    let (result, work) = with_operation_budget(20_000_000, || {
        cylinder(true, false, true, &mut Decider::default())
    });
    assert!(!work.exhausted);
    assert_eq!(result.unwrap_err().code, "shapes.rational-atlas-open");
}

#[test]
fn analytic_sections_use_the_same_periodic_row_proofs_in_both_parameters() {
    for transpose in [false, true] {
        let (result, work) = with_operation_budget(20_000_000, || -> Result<(), Refusal> {
            let mut dec = Decider::default();
            let atlas = cylinder(true, transpose, false, &mut dec)?;
            let cutting = catallax_shapes::SectionSurface::Plane(catallax_shapes::Plane {
                origin: Point3::new(Real::zero(), Real::zero(), q(1, 5)),
                normal: catallax_shapes::Vector3::new(Real::zero(), Real::zero(), Real::one()),
                ref_dir: catallax_shapes::Vector3::new(Real::one(), Real::zero(), Real::zero()),
            });
            let section = atlas.sections(&cutting, &mut dec)?;
            assert_eq!(section.branches().len(), 1);
            assert!(section.branches()[0].is_closed());
            Ok(())
        });
        assert!(!work.exhausted);
        result.unwrap();
    }
}

#[test]
fn cached_preparation_keeps_branch_bytes_and_periodic_ownership() {
    let mut baseline = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let mut records = Vec::new();
        for periodic in [false, true] {
            for transpose in [false, true] {
                for reversed in [false, true] {
                    let legacy = run(periodic, transpose, reversed, false);
                    let cached = run(periodic, transpose, reversed, true);
                    assert_eq!(legacy.0, cached.0);
                    records.push(cached);
                }
            }
        }
        println!(
            "cached pair threads={threads:?} work={:?}",
            records.iter().map(|(_, work)| *work).collect::<Vec<_>>()
        );
        if let Some(baseline) = &baseline {
            assert_eq!(baseline, &records);
        } else {
            baseline = Some(records);
        }
    }
    catallax_numbers::set_threads(None);
}

#[test]
fn extracted_periodic_charts_retain_exact_analytic_incidence() {
    use catallax_shapes::{Cylinder, SectionSurface, Vector3};

    for transpose in [false, true] {
        let mut dec = Decider::default();
        let atlas = cylinder(true, transpose, false, &mut dec).unwrap();
        let mut carrier = Cylinder {
            origin: Point3::new(0.into(), 0.into(), 0.into()),
            axis: Vector3::new(0.into(), 0.into(), 1.into()),
            radius_sq: Real::one(),
            ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
        };
        atlas
            .certify_analytic_carrier(&SectionSurface::Cylinder(carrier.clone()), &mut dec)
            .unwrap();
        carrier.radius_sq = 4.into();
        let error = atlas
            .certify_analytic_carrier(&SectionSurface::Cylinder(carrier), &mut dec)
            .unwrap_err();
        assert_eq!(error.code, "shapes.rational-atlas-carrier");
    }
}

#[test]
fn a_later_bezier_span_cannot_borrow_an_earlier_spans_carrier_identity() {
    use catallax_shapes::{Plane, SectionSurface, Vector3};

    let mut dec = Decider::default();
    let controls = (0..5)
        .flat_map(|u| {
            (0..2).map(move |v| {
                Point3::new(
                    q(u, 2),
                    Real::from_i64(v),
                    Real::from_i64(i64::from(u == 4)),
                )
            })
        })
        .collect();
    // The first span is planar; the second leaves it quadratically with a
    // genuine G1 join. Each extracted chart has the same local span indices.
    let source = BSplineSurface::new(
        [2, 1],
        [
            [0, 0, 0, 1, 1, 2, 2, 2].map(Rational::from_i64).to_vec(),
            [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        ],
        [5, 2],
        controls,
        vec![Real::one(); 10],
        &mut dec,
    )
    .unwrap();
    let atlases = RationalAtlas::partition(source, [false, false], &mut dec).unwrap();
    assert_eq!(atlases.len(), 1);
    let carrier = SectionSurface::Plane(Plane {
        origin: Point3::new(0.into(), 0.into(), 0.into()),
        normal: Vector3::new(0.into(), 0.into(), 1.into()),
        ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
    });
    let error = atlases[0]
        .certify_analytic_carrier(&carrier, &mut dec)
        .unwrap_err();
    assert_eq!(error.code, "shapes.rational-atlas-carrier");
}

#[test]
fn original_spans_keep_complete_circle_components_and_exact_periodic_ownership() {
    let mut baseline = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let mut records = Vec::new();
        for periodic in [false, true] {
            for transpose in [false, true] {
                for swapped in [false, true] {
                    let (result, work) =
                        with_operation_budget(20_000_000, || -> Result<String, Refusal> {
                            let mut dec = Decider::default();
                            let cylinder = cylinder(periodic, transpose, false, &mut dec)?;
                            let plane = plane(&mut dec);
                            let meeting = if swapped {
                                plane.patch_meeting_with_original_spans(&cylinder, &mut dec)?
                            } else {
                                cylinder.patch_meeting_with_original_spans(&plane, &mut dec)?
                            };
                            assert_eq!(meeting.branches().len(), 1);
                            let branch = &meeting.branches()[0];
                            assert_eq!(branch.is_closed(), periodic);
                            let mut record = format!("closed={periodic};");
                            for arc in branch.arcs() {
                                record.push_str(&format!(
                                    "{:?};",
                                    (
                                        arc.first_chart,
                                        arc.first_span,
                                        &arc.first_cell,
                                        arc.second_chart,
                                        arc.second_span,
                                        &arc.second_cell,
                                        arc.certificate,
                                    )
                                ));
                                for end in &arc.ends {
                                    record.push_str(&format!("{:?};", end.enclosure()));
                                }
                            }
                            let curve = meeting.curve_branch(0, &mut dec)?;
                            for (x, y) in [
                                (q(1, 1), q(0, 1)),
                                (q(0, 1), q(1, 1)),
                                (q(-1, 1), q(0, 1)),
                                (q(0, 1), q(-1, 1)),
                                (q(3, 5), q(4, 5)),
                            ] {
                                assert!(curve.contains(&Point3::new(x, y, q(1, 5)), &mut dec)?);
                            }
                            for piece in curve.pieces() {
                                let point = piece
                                    .point_at(&Rational::from_ratio(1, 2).unwrap(), &mut dec)?;
                                let circle =
                                    &(&point.x.square() + &point.y.square()) - &Real::one();
                                assert_eq!(
                                    dec.rational_identity(&circle)?,
                                    Some(Rational::from_i64(0))
                                );
                                assert_eq!(
                                    dec.rational_identity(&(&point.z - &q(1, 5)))?,
                                    Some(Rational::from_i64(0))
                                );
                            }
                            Ok(record)
                        });
                    assert!(
                        !work.exhausted,
                        "original-span circle proof exhausted its unchanged allowance"
                    );
                    records.push((result.unwrap(), work.used));
                }
            }
        }
        println!(
            "original span threads={threads:?} work={:?}",
            records.iter().map(|(_, work)| *work).collect::<Vec<_>>()
        );
        if let Some(baseline) = &baseline {
            assert_eq!(baseline, &records);
        } else {
            baseline = Some(records);
        }
    }
    catallax_numbers::set_threads(None);
}
