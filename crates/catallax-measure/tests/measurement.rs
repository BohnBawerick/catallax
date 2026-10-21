//! Blind cases from docs/pieces/06-measurement.md. Expected values are closed forms.
mod common;
use catallax_topology::Body;
use common::*;
use std::time::{Duration, Instant};

blind_test!(cases_1_to_14_primitives_hand_bodies_and_every_face, {
    let a = default_accuracy();
    for case in CASES {
        let f = build(case, &r(1), &Pose::identity());
        let result = run_mass(&f, &a);
        assert_mass(&result, &f.expected, &a, true);
        let volume = measure::volume(&f.model, f.body, &a, &mut Decider::default()).unwrap();
        let area = measure::area(&f.model, f.body, &a, &mut Decider::default()).unwrap();
        assert_scalar(&volume, &f.expected.volume, &a, true, f.expected.exact);
        assert_scalar(&area, &f.expected.area, &a, true, f.expected.exact);
        assert_eq!(
            volume.value, result.value.volume,
            "{case:?}: volume entry points diverge"
        );
        assert_eq!(
            area.value, result.value.area,
            "{case:?}: area entry points diverge"
        );
        assert_faces(&f, &a);
    }
});

blind_test!(
    case_16_accuracy_ladder_including_zero_and_signed_coordinates,
    {
        for case in CURVED {
            // Tightest accuracies only on three bodies, as the runtime clause permits.
            let digits: &[i32] =
                if matches!(case, Case::Sphere | Case::Frustum | Case::HalfPositive) {
                    &[6, 12, 18, 24, 30]
                } else {
                    &[6, 12, 18]
                };
            for &n in digits {
                let a = accuracy(n);
                let f = build(case, &r(1), &Pose::identity());
                assert_mass(&run_mass(&f, &a), &f.expected, &a, true);
                assert_faces(&f, &a);
            }
        }
    }
);

blind_test!(cases_18_19_scaled_placed_and_threaded_math_enclosures, {
    let a = default_accuracy();
    let mut count = 0;
    for case in CASES {
        for k in -30..=30 {
            let f = build(case, &ten(k), &Pose::identity());
            count += assert_mass(&run_mass(&f, &a), &f.expected, &a, true);
        }
        for pose in poses() {
            let f = build(case, &r(1), &pose);
            count += assert_mass(&run_mass(&f, &a), &f.expected, &a, true);
        }
        let pose = poses().remove(0);
        let f = build(case, &r(1), &pose);
        let baseline = run_mass(&f, &a);
        count += assert_mass(&baseline, &f.expected, &a, true);
        let repeat = run_mass(&f, &a);
        assert_eq!(baseline.value, repeat.value);
        assert_eq!(baseline.marks, repeat.marks);
        count += assert_mass(&repeat, &f.expected, &a, true);
        let results = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..8)
                .map(|_| {
                    scope.spawn(|| {
                        // Each thread gets a fresh model and Decider; no shared caches required.
                        let f = build(case, &r(1), &pose);
                        let result = run_mass(&f, &a);
                        let count = assert_mass(&result, &f.expected, &a, true);
                        let volume =
                            measure::volume(&f.model, f.body, &a, &mut Decider::default()).unwrap();
                        let area =
                            measure::area(&f.model, f.body, &a, &mut Decider::default()).unwrap();
                        assert_eq!(result.value.volume, volume.value);
                        assert_eq!(result.value.area, area.value);
                        let faces: Vec<_> = f
                            .faces
                            .iter()
                            .map(|(id, expected)| {
                                let answer =
                                    measure::face_area(&f.model, *id, &a, &mut Decider::default())
                                        .unwrap();
                                assert_scalar(&answer, expected, &a, true, f.expected.exact);
                                answer.value
                            })
                            .collect();
                        (result, count, faces)
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        let expected_faces: Vec<_> = f
            .faces
            .iter()
            .map(|(id, _)| {
                measure::face_area(&f.model, *id, &a, &mut Decider::default())
                    .unwrap()
                    .value
            })
            .collect();
        for (result, n, faces) in results {
            assert_eq!(
                baseline.value, result.value,
                "{case:?}: thread-dependent endpoints"
            );
            assert_eq!(baseline.marks, result.marks);
            assert_eq!(expected_faces, faces);
            count += n;
        }
    }
    // 16 variants: fifteen solids with 5 quantities, Empty with 2.
    // (61 scales + 20 poses + 2 repeats + 8 threads) * 77 = 7007.
    assert_eq!(count, 7007);
    assert!(count >= 3000);
});

blind_test!(case_17_budget_marks_bounds_and_zero_work, {
    let a = default_accuracy();
    for case in CASES {
        let f = build(case, &r(1), &poses().remove(1));
        for work in [10_000, 0] {
            let policy = BudgetPolicy {
                max_precision_bits: 64,
                max_work: work,
            };
            type Property = fn(
                &catallax_topology::Model,
                catallax_topology::BodyId,
                &Accuracy,
                &mut Decider,
            ) -> Answered<Enclosure>;
            let calls: [(Property, &Expr); 2] = [
                (measure::volume, &f.expected.volume),
                (measure::area, &f.expected.area),
            ];
            for (call, expected) in calls {
                match call(&f.model, f.body, &a, &mut Decider::new(policy)) {
                    Ok(answer) => {
                        assert!(
                            work != 0 || f.expected.exact,
                            "non-rational property answered without Ball work"
                        );
                        assert_scalar(&answer, expected, &a, false, f.expected.exact);
                    }
                    Err(err) => {
                        assert!(!f.expected.exact, "rational body refused: {err}");
                        if work == 0 && err.code == "measure.budget-exhausted" {
                            continue;
                        }
                        upstream_budget_refusal(err.code);
                    }
                }
            }
            match measure::mass_properties(&f.model, f.body, &a, &mut Decider::new(policy)) {
                Ok(answer) => {
                    assert!(work != 0 || f.expected.exact);
                    assert_mass(&answer, &f.expected, &a, false);
                }
                Err(err) => {
                    assert!(!f.expected.exact, "rational body refused: {err}");
                    if err.code != "measure.centroid-not-separated"
                        && !(work == 0 && err.code == "measure.budget-exhausted")
                    {
                        upstream_budget_refusal(err.code);
                    }
                }
            }
            // Rational planar faces of curved bodies may still answer with zero work.
            for (id, expected) in &f.faces {
                match measure::face_area(&f.model, *id, &a, &mut Decider::new(policy)) {
                    Ok(answer) => assert_scalar(&answer, expected, &a, false, f.expected.exact),
                    Err(err) => {
                        assert!(!f.expected.exact);
                        if !(work == 0 && err.code == "measure.budget-exhausted") {
                            upstream_budget_refusal(err.code);
                        }
                    }
                }
            }
        }
    }
});

blind_test!(
    budget_marks_use_actual_width_and_centroid_maximum_in_property_order,
    {
        // Shift and rotate the asymmetric cap so all three centroid coordinates are
        // nonzero and have distinct scales. The absolute floor is deliberately tiny.
        let a = Accuracy {
            relative: ten(-80),
            absolute: ten(-90),
        };
        let policy = BudgetPolicy {
            max_precision_bits: 64,
            max_work: 10_000,
        };
        for case in [Case::Frustum, Case::Cap, Case::HalfPositive, Case::Torus] {
            let f = build(case, &r(1), &poses().remove(2));
            match measure::mass_properties(&f.model, f.body, &a, &mut Decider::new(policy)) {
                Ok(answer) => {
                    assert_mass(&answer, &f.expected, &a, false);
                }
                Err(err) if err.code == "measure.centroid-not-separated" => {}
                Err(err) => upstream_budget_refusal(err.code),
            }
        }
    }
);

blind_test!(
    case_20_expired_deadline_preserves_default_answer_or_exact_code,
    {
        let a = default_accuracy();
        let expired = || {
            let mut d = Decider::default();
            d.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
            d
        };
        for case in CASES {
            let f = build(case, &r(1), &Pose::identity());
            let baseline = run_mass(&f, &a);
            match measure::mass_properties(&f.model, f.body, &a, &mut expired()) {
                Ok(answer) => {
                    assert_eq!(answer.value, baseline.value);
                    assert_eq!(answer.marks, baseline.marks);
                }
                Err(err) => assert_eq!(err.code, "numbers.deadline"),
            }
            for call in [measure::volume, measure::area] {
                let baseline = call(&f.model, f.body, &a, &mut Decider::default()).unwrap();
                match call(&f.model, f.body, &a, &mut expired()) {
                    Ok(answer) => {
                        assert_eq!(answer.value, baseline.value);
                        assert_eq!(answer.marks, baseline.marks);
                    }
                    Err(err) => assert_eq!(err.code, "numbers.deadline"),
                }
            }
            for (face, _) in &f.faces {
                let baseline =
                    measure::face_area(&f.model, *face, &a, &mut Decider::default()).unwrap();
                match measure::face_area(&f.model, *face, &a, &mut expired()) {
                    Ok(answer) => {
                        assert_eq!(answer.value, baseline.value);
                        assert_eq!(answer.marks, baseline.marks);
                    }
                    Err(err) => assert_eq!(err.code, "numbers.deadline"),
                }
            }
        }
    }
);

blind_test!(body_marks_stay_on_body_and_do_not_inflate_measurement, {
    let a = default_accuracy();
    for case in [Case::Block, Case::Sphere, Case::HalfNegative] {
        let mut f = build(case, &r(1), &Pose::identity());
        let baseline = run_mass(&f, &a);
        let original = f.model.body(f.body).clone();
        let mut marks = MarkLog::new();
        marks.push(Mark::new(
            "numbers.precision-limit",
            BoundKind::Distance,
            q(1, 7),
            "constructed body uncertainty",
        ));
        let marked = f.model.add_body(Body {
            marks: marks.clone(),
            ..original
        });
        let facts = f.model.geom().facts().len();
        let result =
            measure::mass_properties(&f.model, marked, &a, &mut Decider::default()).unwrap();
        assert_mass(&result, &f.expected, &a, true);
        assert_eq!(baseline.value, result.value);
        assert_eq!(baseline.marks, result.marks);
        assert_eq!(f.model.body(marked).marks, marks);
        assert_eq!(f.model.geom().facts().len(), facts);
        for call in [measure::volume, measure::area] {
            let original = call(&f.model, f.body, &a, &mut Decider::default()).unwrap();
            let result = call(&f.model, marked, &a, &mut Decider::default()).unwrap();
            assert_eq!(original.value, result.value);
            assert_eq!(original.marks, result.marks);
        }
    }
});

blind_test!(enclosure_closed_endpoints_width_and_exactness, {
    for (lo, hi, inside, outside) in [
        (q(-7, 3), q(-1, 3), q(-4, 3), r(0)),
        (r(-2), r(3), r(0), r(-3)),
        (q(2, 7), q(2, 7), q(2, 7), q(3, 7)),
    ] {
        let e = Enclosure {
            lo: lo.clone(),
            hi: hi.clone(),
        };
        assert_eq!(e.width(), &hi - &lo);
        assert_eq!(e.is_exact(), lo == hi);
        for p in [lo, hi, inside] {
            assert!(e.contains(&p));
        }
        assert!(!e.contains(&outside));
    }
});

blind_test!(default_accuracy_controls_actual_measurements, {
    let a = Accuracy::default();
    assert_eq!(a.relative, ten(-12));
    assert_eq!(a.absolute, ten(-30));
    let f = build(Case::Sphere, &r(1), &Pose::identity());
    let answer = run_mass(&f, &a);
    assert_mass(&answer, &f.expected, &default_accuracy(), true);
});

blind_test!(
    case_21_boolean_results_against_closed_forms,
    "comparison 06: Steinmetz prerequisite refuses boolean.imprint-ambiguous before measurement",
    {
        use common::boolean_cases::*;
        let a = default_accuracy();
        for case in BOOLEAN_CASES {
            let f = build_boolean(case);
            let vol = measure::volume(&f.model, f.body, &a, &mut Decider::default()).unwrap();
            assert_scalar(&vol, &f.volume, &a, true, false);
            if let Some(expected) = f.mass {
                let mass = measure::mass_properties(&f.model, f.body, &a, &mut Decider::default())
                    .unwrap();
                assert_mass(&mass, &expected, &a, true);
                assert_eq!(vol.value, mass.value.volume);
                let ar = measure::area(&f.model, f.body, &a, &mut Decider::default()).unwrap();
                assert_scalar(&ar, &expected.area, &a, true, false);
                assert_eq!(ar.value, mass.value.area);
            }
        }
    }
);

blind_test!(
    centroid_separation_ignores_absolute_floor_or_refuses_by_named_code,
    {
        let a = Accuracy {
            relative: r(0),
            absolute: r(1),
        };
        // Volume is 4*pi/3 * (1-(1-epsilon)^3), strictly positive even when
        // cancellation leaves a low-precision volume interval straddling zero.
        // A creator may avoid cancellation analytically, so refusal is not mandatory.
        for exponent in [-30, -100, -300] {
            let f = thin_spherical_shell(&ten(exponent));
            assert_mass(&run_mass(&f, &a), &f.expected, &a, true);
            for work in [0, 1, 16, 10_000] {
                let policy = BudgetPolicy {
                    max_precision_bits: 64,
                    max_work: work,
                };
                match measure::mass_properties(&f.model, f.body, &a, &mut Decider::new(policy)) {
                    Ok(answer) => {
                        assert_mass(&answer, &f.expected, &a, false);
                    }
                    Err(err) => match err.code {
                        "measure.centroid-not-separated" | "measure.budget-exhausted" => {}
                        code => upstream_budget_refusal(code),
                    },
                }
            }
        }
    }
);
