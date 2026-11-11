//! Complete patch-plane search, exact branch queries and candidate tubes.
use catallax_intersect::{
    Component, Contact, CurveForm, Cyclic, SurfacePair, curve_form, curve_is_closed, cyclic_order,
    intersect_surfaces, on_curve, order_on_curve, tangent_at,
};
use catallax_numbers::{Decider, Outcome, Rational, Real, StrictOrder, with_operation_budget};
use catallax_shapes::{
    Curve, GeomStore, Plane, Point3, PointDef, RationalPatch, SectionSurface, Surface, Vector3,
};

#[allow(dead_code)]
mod patch_fixtures;
use patch_fixtures::{candidates, equal, height_patch, loft_like, p, r};

fn cutting_plane(case: &str) -> Plane {
    if case == "loft" {
        return Plane {
            origin: p("5/2", "0", "0"),
            normal: Vector3::new(Real::one(), Real::one(), Real::one()),
            ref_dir: Vector3::new(Real::one(), Real::from_i64(-1), Real::zero()),
        };
    }
    let height = match case {
        "wavy" => r("7/20"),
        "hill" => r("3/10"),
        "tiny" => r("225/256") - r("1/1000000"),
        "tangent" => r("225/256"),
        _ => r("225/256") + r("1/1000000"),
    };
    Plane {
        origin: Point3::new(Real::zero(), Real::zero(), height.into()),
        normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
        ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
    }
}
fn patch(case: &str, dec: &mut Decider) -> RationalPatch {
    if case == "loft" {
        loft_like(dec)
    } else {
        height_patch(case != "wavy", dec)
    }
}
#[test]
fn complete_search_finds_every_open_arc_loop_and_tangent_contact() {
    for (case, branches, arcs, empty, depth) in [
        ("loft", 1, 2, 0, 1),
        ("wavy", 2, 6, 4, 1),
        ("hill", 1, 10, 15, 2),
        ("tiny", 1, 12, 76, 21),
        ("tangent", 0, 0, 0, 0),
        ("above", 0, 0, 12, 2),
    ] {
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = patch(case, &mut dec);
            let result = patch.sections(&SectionSurface::Plane(cutting_plane(case)), &mut dec);
            if let Ok(section) = &result
                && case == "loft"
            {
                let [a, b] = section.branches()[0].endpoints();
                let ends = [p("0", "0", "5/2"), p("2", "0", "1/2")];
                if a.x.rational_ref() == Some(&r("0")) {
                    equal(&a, &ends[0], &mut dec);
                    equal(&b, &ends[1], &mut dec);
                } else {
                    equal(&b, &ends[0], &mut dec);
                    equal(&a, &ends[1], &mut dec);
                }
            }
            result
        });
        assert!(!work.exhausted, "{case}: {} units", work.used);
        if case == "tangent" {
            let refusal = result.unwrap_err();
            assert_eq!(refusal.code, "intersect.patch-plane-tangent");
            assert!(refusal.detail.contains("(5/2, 5/2, 225/256)"));
            println!("WORK {case}: {} units; {}", work.used, refusal.detail);
            continue;
        }
        let section = result.unwrap();
        let search = section.search().unwrap();
        assert_eq!(section.branches().len(), branches);
        assert_eq!(
            (search.arc_cells, search.empty_cells, search.depth),
            (arcs, empty, depth)
        );
        for b in section.branches() {
            assert_eq!(b.is_closed(), case == "hill" || case == "tiny");
        }
        println!(
            "WORK {case}: {} operation units, {} decision units; {branches} branches, {arcs} arc and {empty} empty cells, depth {depth}",
            work.used, search.work
        );
    }
}

#[test]
fn intersection_pair_returns_exact_branches_and_answers_parameter_queries() {
    for case in ["loft", "wavy", "hill", "tiny", "above", "tangent"] {
        for swap in [false, true] {
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let patch = patch(case, &mut dec);
                let plane = cutting_plane(case);
                let mut geom = GeomStore::new();
                let a = geom.add_surface(Surface::Plane(plane.clone()));
                let b = geom.add_surface(Surface::BSpline(patch.into()));
                let result = if swap {
                    intersect_surfaces(&mut geom, b, a, &mut dec)
                } else {
                    intersect_surfaces(&mut geom, a, b, &mut dec)
                };
                if case == "tangent" {
                    assert_eq!(result.unwrap_err().code, "intersect.patch-plane-tangent");
                    return;
                }
                let answer = result.unwrap();
                assert!(answer.is_exact());
                if case == "above" {
                    assert_eq!(answer.value, SurfacePair::Disjoint);
                    return;
                }
                let SurfacePair::Meets(components) = answer.value else {
                    panic!("expected exact branch components")
                };
                assert_eq!(components.len(), if case == "wavy" { 2 } else { 1 });
                for component in components {
                    let Component::Curve {
                        curve,
                        contact,
                        singular,
                    } = component
                    else {
                        panic!("expected curve")
                    };
                    assert_eq!(contact, Contact::Crossing);
                    assert!(singular.is_empty());
                    let Curve::PatchSection(branch) = geom.curve(curve) else {
                        panic!("candidate must not replace the exact branch")
                    };
                    let branch = branch.clone();
                    assert_eq!(branch.follows_surface_first(&mut dec).unwrap(), !swap);
                    assert!(matches!(
                        curve_form(&geom, curve, &mut dec).unwrap(),
                        CurveForm::PatchSection(_)
                    ));
                    assert_eq!(
                        curve_is_closed(&geom, curve, &mut dec).unwrap(),
                        branch.is_closed()
                    );
                    let [start, end] = branch.endpoints();
                    let p = geom.add_point(PointDef::Explicit(start.clone()));
                    let q = geom.add_point(PointDef::Explicit(end));
                    assert!(on_curve(&geom, curve, p, &mut dec).is_coincident());
                    let raw =
                        geom.add_curve(Curve::Intersection(catallax_shapes::IntersectionCurve {
                            surfaces: [a, b],
                            witness: p,
                        }));
                    assert_eq!(
                        curve_form(&geom, raw, &mut dec).unwrap_err().code,
                        "intersect.patch-plane-selector"
                    );
                    let outside = branch
                        .patch()
                        .point_at(&[Real::zero(), Real::zero()], &mut dec)
                        .unwrap();
                    let outside = geom.add_point(PointDef::Explicit(outside));
                    assert!(matches!(
                        on_curve(&geom, curve, outside, &mut dec),
                        Outcome::Proven(_)
                    ));
                    if !branch.is_closed() {
                        assert!(matches!(
                            order_on_curve(&geom, curve, p, q, &mut dec),
                            Outcome::Proven(StrictOrder::Less)
                        ));
                        assert!(matches!(
                            order_on_curve(&geom, curve, q, p, &mut dec),
                            Outcome::Proven(StrictOrder::Greater)
                        ));
                        assert!(order_on_curve(&geom, curve, p, p, &mut dec).is_coincident());
                    } else {
                        let ids = branch
                            .witnesses()
                            .take(3)
                            .map(|p| geom.add_point(PointDef::Explicit(p)))
                            .collect::<Vec<_>>();
                        assert_eq!(ids.len(), 3);
                        assert!(matches!(
                            cyclic_order(&geom, curve, [ids[0], ids[1], ids[2]], &mut dec),
                            Outcome::Proven(Cyclic::Forward)
                        ));
                        assert!(matches!(
                            cyclic_order(&geom, curve, [ids[0], ids[2], ids[1]], &mut dec),
                            Outcome::Proven(Cyclic::Backward)
                        ));
                    }
                    let tangents = tangent_at(&geom, curve, p, &mut dec).unwrap();
                    for tangent in tangents {
                        assert!(dec.sign(&tangent.dot(&plane.normal)).is_coincident());
                        assert!(matches!(
                            dec.sign(&tangent.norm_squared()),
                            Outcome::Proven(_)
                        ));
                    }
                }
            });
            assert!(!work.exhausted, "{case} swap={swap}: {} units", work.used);
            println!("WORK pair {case} swap={swap}: {} units", work.used);
        }
    }
}

#[test]
fn oracle_tubes_respect_the_target_cap_end_caps_and_work_allowance() {
    let tight = candidates(include_str!("fixtures/patch-plane/lofttight.txt"))
        .pop()
        .unwrap();
    let default = candidates(include_str!("fixtures/patch-plane/loft.txt"))
        .pop()
        .unwrap();
    assert_eq!(tight.pieces.len(), 21);
    assert_eq!(default.pieces.len(), 9);
    assert!(candidates(include_str!("fixtures/patch-plane/tiny.txt")).is_empty());
    assert!(candidates(include_str!("fixtures/patch-plane/tangent.txt")).is_empty());
    for name in ["tight", "default", "shifted", "short"] {
        let mut candidate = if name == "default" {
            default.clone()
        } else {
            tight.clone()
        };
        if name == "shifted" {
            for piece in &mut candidate.pieces {
                for p in piece {
                    p[1] = &p[1] + &r("1/1000");
                }
            }
        }
        if name == "short" {
            candidate.pieces[0][0][0] = r("999/1000");
        }
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = loft_like(&mut dec);
            let section = patch
                .sections(&SectionSurface::Plane(cutting_plane("loft")), &mut dec)
                .unwrap();
            section.branches()[0].certify_candidate(&candidate, &mut dec)
        });
        println!(
            "WORK tube {name}: {} operation units, exhausted={}",
            work.used, work.exhausted
        );
        assert!(
            !work.exhausted,
            "the first M10 tube exceeded the unchanged 20M allowance"
        );
        if name == "tight" {
            let proof = result.unwrap();
            assert!(proof.bound_mm() <= &r("1/10000"));
            assert_eq!(proof.cap_mm(), &r("1/1000"));
            assert!(proof.transverse_slope_lower() > &Rational::zero());
            assert_eq!(proof.candidate().pieces.len(), 21);
            assert!(!proof.definition().is_closed());
            assert!(matches!(
                proof.proof(),
                catallax_numbers::GeometricProof::Approximation { .. }
            ));
            println!(
                "BOUND tube tight: U={} mm, rho={}, speed={}, slope_lower={}, decision_work={}",
                proof.bound_mm(),
                proof.half_width(),
                proof.speed_bound(),
                proof.transverse_slope_lower(),
                proof.work()
            );
        } else {
            let refusal = result.unwrap_err();
            assert_eq!(
                refusal.code,
                if name == "short" {
                    "intersect.patch-plane-candidate-endpoint"
                } else {
                    "intersect.patch-plane-candidate-cap"
                }
            );
        }
    }
    let mut dec = Decider::default();
    let patch = loft_like(&mut dec);
    let section = patch
        .sections(&SectionSurface::Plane(cutting_plane("loft")), &mut dec)
        .unwrap();
    let (result, work) = with_operation_budget(100, || {
        section.branches()[0].certify_candidate(&tight, &mut dec)
    });
    assert!(work.exhausted);
    assert!(result.is_err());
}

#[test]
fn same_exact_answers_with_one_two_and_four_workers() {
    let cases = ["loft", "wavy", "hill", "tiny", "tangent", "above"];
    let mut baseline = None;
    for workers in [1, 2, 4] {
        let mut rows = std::thread::scope(|scope| {
            let jobs = (0..workers)
                .map(|worker| {
                    scope.spawn(move || {
                        cases
                            .iter()
                            .enumerate()
                            .skip(worker)
                            .step_by(workers)
                            .map(|(index, &case)| {
                                let (answer, work) = with_operation_budget(20_000_000, || {
                                    let mut dec = Decider::default();
                                    let patch = patch(case, &mut dec);
                                    match patch.sections(
                                        &SectionSurface::Plane(cutting_plane(case)),
                                        &mut dec,
                                    ) {
                                        Err(error) => format!("{} {}", error.code, error.detail),
                                        Ok(section) => {
                                            let search = section.search().unwrap();
                                            let mut answer = format!(
                                                "{case} {} {} {}",
                                                search.arc_cells, search.empty_cells, search.depth
                                            );
                                            for branch in section.branches() {
                                                answer.push_str(&format!(
                                                    " closed={}",
                                                    branch.is_closed()
                                                ));
                                                for end in branch.endpoints() {
                                                    for coordinate in [&end.x, &end.y, &end.z] {
                                                        let b = dec
                                                            .enclose_budgeted_at(coordinate, 64)
                                                            .unwrap();
                                                        answer.push_str(&format!(
                                                            " [{},{}]",
                                                            b.lo(),
                                                            b.hi()
                                                        ));
                                                    }
                                                }
                                            }
                                            if case == "loft" {
                                                let candidate = candidates(include_str!(
                                                    "fixtures/patch-plane/lofttight.txt"
                                                ))
                                                .pop()
                                                .unwrap();
                                                let proof = section.branches()[0]
                                                    .certify_candidate(&candidate, &mut dec)
                                                    .unwrap();
                                                answer
                                                    .push_str(&format!(" U={}", proof.bound_mm()));
                                            }
                                            answer
                                        }
                                    }
                                });
                                assert!(
                                    !work.exhausted,
                                    "{case}, {workers} workers: {} units",
                                    work.used
                                );
                                (index, answer, work.used)
                            })
                            .collect::<Vec<_>>()
                    })
                })
                .collect::<Vec<_>>();
            jobs.into_iter()
                .flat_map(|job| job.join().unwrap())
                .collect::<Vec<_>>()
        });
        rows.sort_by_key(|row| row.0);
        let answers = rows
            .iter()
            .map(|row| (row.1.clone(), row.2))
            .collect::<Vec<_>>();
        let units = rows
            .iter()
            .map(|row| (cases[row.0], row.2))
            .collect::<Vec<_>>();
        println!("WORK workers={workers}: {units:?}");
        if let Some(baseline) = &baseline {
            assert_eq!(&answers, baseline);
        } else {
            baseline = Some(answers);
        }
    }
}
