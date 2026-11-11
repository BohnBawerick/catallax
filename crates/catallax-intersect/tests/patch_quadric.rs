//! Complete patch-cylinder and patch-sphere search, exact branch queries and
//! candidate tubes. Only the substituted cut polynomial differs from the plane.
use catallax_intersect::{
    Component, Contact, CurveForm, Cyclic, SurfacePair, curve_form, curve_is_closed, cyclic_order,
    intersect_surfaces, on_curve, order_on_curve, tangent_at,
};
use catallax_numbers::{Decider, Outcome, Rational, Real, StrictOrder, with_operation_budget};
use catallax_shapes::{
    Curve, Cylinder, GeomStore, PointDef, RationalPatch, SectionCandidate, SectionSurface, Sphere,
    Surface, Vector3,
};

#[allow(dead_code)]
mod patch_fixtures;
use patch_fixtures::{candidates, equal, height_patch, loft_like, p, r};

const SURFACES: [&str; 2] = ["cylinder", "sphere"];
const CUTS: [&str; 5] = ["loft", "hill", "tiny", "tangent", "above"];

fn axis(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}

/// Loft cuts have radius 5/2 about the origin, with exact ends (0,0,5/2)
/// and (2,0,3/2). Hill cuts have radius 1 and their lowest point at 3/10,
/// one millionth below the hilltop, at it, or one millionth above it.
fn cutter(surface: &str, cut: &str) -> SectionSurface {
    let peak = r("225/256");
    let low = match cut {
        "loft" => Rational::zero(),
        "hill" => r("3/10"),
        "tiny" => &peak - &r("1/1000000"),
        "tangent" => peak,
        _ => &peak + &r("1/1000000"),
    };
    let center_z = Real::from_rational(&low + &Rational::one());
    match (surface, cut) {
        ("cylinder", "loft") => SectionSurface::Cylinder(Cylinder {
            origin: p("0", "0", "0"),
            axis: axis(0, 1, 0),
            radius_sq: Real::from_rational(r("25/4")),
            ref_dir: axis(0, 0, 1),
        }),
        ("cylinder", _) => SectionSurface::Cylinder(Cylinder {
            origin: catallax_shapes::Point3::new(Real::zero(), r("5/2").into(), center_z),
            axis: axis(1, 0, 0),
            radius_sq: Real::one(),
            ref_dir: axis(0, 1, 0),
        }),
        ("sphere", "loft") => SectionSurface::Sphere(Sphere {
            center: p("0", "0", "0"),
            radius: Real::from_rational(r("5/2")),
            axis: axis(0, 0, 1),
            ref_dir: axis(1, 0, 0),
        }),
        _ => SectionSurface::Sphere(Sphere {
            center: catallax_shapes::Point3::new(r("5/2").into(), r("5/2").into(), center_z),
            radius: Real::one(),
            axis: axis(0, 0, 1),
            ref_dir: axis(1, 0, 0),
        }),
    }
}
fn patch(cut: &str, dec: &mut Decider) -> RationalPatch {
    if cut == "loft" {
        loft_like(dec)
    } else {
        height_patch(true, dec)
    }
}
fn fixture(case: &str) -> &'static str {
    match case {
        "cylinder-loft" => include_str!("fixtures/patch-quadric/cylinder-loft.txt"),
        "cylinder-lofttight" => include_str!("fixtures/patch-quadric/cylinder-lofttight.txt"),
        "cylinder-hill" => include_str!("fixtures/patch-quadric/cylinder-hill.txt"),
        "cylinder-tiny" => include_str!("fixtures/patch-quadric/cylinder-tiny.txt"),
        "cylinder-tangent" => include_str!("fixtures/patch-quadric/cylinder-tangent.txt"),
        "sphere-loft" => include_str!("fixtures/patch-quadric/sphere-loft.txt"),
        "sphere-lofttight" => include_str!("fixtures/patch-quadric/sphere-lofttight.txt"),
        "sphere-hill" => include_str!("fixtures/patch-quadric/sphere-hill.txt"),
        "sphere-tiny" => include_str!("fixtures/patch-quadric/sphere-tiny.txt"),
        "sphere-tangent" => include_str!("fixtures/patch-quadric/sphere-tangent.txt"),
        _ => unreachable!("no fixture for {case}"),
    }
}

/// A tracer may end a hair off the patch border; the checker refuses that.
/// The caller may move an end control within 1e-12 of the border in `u`
/// onto it. The result is still an untrusted candidate for the checker.
fn on_border(mut candidate: SectionCandidate) -> SectionCandidate {
    let last = candidate.pieces.len() - 1;
    for (piece, control) in [(0, 0), (last, 3)] {
        let u = &mut candidate.pieces[piece][control][0];
        for border in [Rational::zero(), Rational::one()] {
            if (&*u - &border).abs() < r("1e-12") {
                *u = border;
            }
        }
    }
    candidate
}

/// One deterministic answer line: cell counts, closure, endpoint enclosures
/// and, for loft cuts, the tight candidate bound. Refusals keep their detail.
fn answer(surface: &str, cut: &str, dec: &mut Decider) -> String {
    let patch = patch(cut, dec);
    let section = match patch.sections(&cutter(surface, cut), dec) {
        Err(error) => return format!("{} {}", error.code, error.detail),
        Ok(section) => section,
    };
    let search = section.search().unwrap();
    let mut answer = format!(
        "{surface}-{cut} {} {} {}",
        search.arc_cells, search.empty_cells, search.depth
    );
    for branch in section.branches() {
        answer.push_str(&format!(" closed={}", branch.is_closed()));
        for end in branch.endpoints() {
            for coordinate in [&end.x, &end.y, &end.z] {
                let b = dec.enclose_budgeted_at(coordinate, 64).unwrap();
                answer.push_str(&format!(" [{},{}]", b.lo(), b.hi()));
            }
        }
    }
    if cut == "loft" {
        let candidate = on_border(
            candidates(fixture(&format!("{surface}-lofttight")))
                .pop()
                .unwrap(),
        );
        let proof = section.branches()[0]
            .certify_candidate(&candidate, dec)
            .unwrap();
        answer.push_str(&format!(" U={}", proof.bound_mm()));
    }
    answer
}

#[test]
fn complete_search_finds_open_arcs_missed_loops_and_tangent_contacts() {
    for surface in SURFACES {
        for cut in CUTS {
            let (result, work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let patch = patch(cut, &mut dec);
                let result = patch.sections(&cutter(surface, cut), &mut dec);
                if let Ok(section) = &result
                    && cut == "loft"
                {
                    assert_eq!(section.branches().len(), 1);
                    let [a, b] = section.branches()[0].endpoints();
                    let ends = [p("0", "0", "5/2"), p("2", "0", "3/2")];
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
            assert!(!work.exhausted, "{surface}-{cut}: {} units", work.used);
            if cut == "tangent" {
                let refusal = result.unwrap_err();
                assert_eq!(refusal.code, format!("intersect.patch-{surface}-tangent"));
                assert!(
                    refusal
                        .detail
                        .contains("(5/2, 5/2, 225/256), parameters (3/2, 3/2)"),
                    "{}",
                    refusal.detail
                );
                println!(
                    "WORK {surface}-{cut}: {} units; {}",
                    work.used, refusal.detail
                );
                continue;
            }
            let section = result.unwrap();
            let search = section.search().unwrap();
            let branches = section.branches().len();
            println!(
                "WORK {surface}-{cut}: {} operation units, {} decision units; {branches} branches, {} arc and {} empty cells, depth {}",
                work.used, search.work, search.arc_cells, search.empty_cells, search.depth
            );
            assert_eq!(branches, usize::from(cut != "above"));
            let expected = match (surface, cut) {
                ("cylinder", "loft") | ("sphere", "loft") => (2, 1, 2),
                ("cylinder", "hill") => (10, 5, 2),
                ("cylinder", "tiny") => (12, 80, 22),
                ("sphere", "hill") => (12, 10, 3),
                ("sphere", "tiny") => (8, 84, 22),
                _ => (0, 12, 2),
            };
            assert_eq!(
                (search.arc_cells, search.empty_cells, search.depth),
                expected
            );
            for b in section.branches() {
                assert_eq!(b.is_closed(), cut == "hill" || cut == "tiny");
            }
        }
    }
}

#[test]
fn the_tracer_misses_the_tiny_loops_and_the_tangent_contacts() {
    for surface in SURFACES {
        for cut in ["tiny", "tangent"] {
            assert!(candidates(fixture(&format!("{surface}-{cut}"))).is_empty());
        }
        assert_eq!(candidates(fixture(&format!("{surface}-hill"))).len(), 1);
    }
}

#[test]
fn intersection_pair_returns_exact_branches_and_answers_parameter_queries() {
    for surface in SURFACES {
        for cut in CUTS {
            for swap in [false, true] {
                let (tangents, work) = with_operation_budget(20_000_000, || {
                    let mut tangents = Vec::new();
                    let mut dec = Decider::default();
                    let patch = patch(cut, &mut dec);
                    let analytic = cutter(surface, cut).surface();
                    let mut geom = GeomStore::new();
                    let a = geom.add_surface(analytic.clone());
                    let b = geom.add_surface(Surface::BSpline(patch.into()));
                    let result = if swap {
                        intersect_surfaces(&mut geom, b, a, &mut dec)
                    } else {
                        intersect_surfaces(&mut geom, a, b, &mut dec)
                    };
                    if cut == "tangent" {
                        assert_eq!(
                            result.unwrap_err().code,
                            format!("intersect.patch-{surface}-tangent")
                        );
                        return tangents;
                    }
                    let answer = result.unwrap();
                    assert!(answer.is_exact());
                    if cut == "above" {
                        assert_eq!(answer.value, SurfacePair::Disjoint);
                        return tangents;
                    }
                    let SurfacePair::Meets(components) = answer.value else {
                        panic!("expected exact branch components")
                    };
                    assert_eq!(components.len(), 1);
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
                        let raw = geom.add_curve(Curve::Intersection(
                            catallax_shapes::IntersectionCurve {
                                surfaces: [a, b],
                                witness: p,
                            },
                        ));
                        assert_eq!(
                            curve_form(&geom, raw, &mut dec).unwrap_err().code,
                            format!("intersect.patch-{surface}-selector")
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
                        let gradient = analytic.gradient(&start).unwrap();
                        for tangent in tangent_at(&geom, curve, p, &mut dec).unwrap() {
                            tangents.push((tangent, gradient.clone()));
                        }
                    }
                    tangents
                });
                assert!(
                    !work.exhausted,
                    "{surface}-{cut} swap={swap}: {} units",
                    work.used
                );
                println!("WORK pair {surface}-{cut} swap={swap}: {} units", work.used);
                // This check is the test's work, not the query's. At algebraic
                // loop witnesses the default decider may stop before an exact
                // zero proof; it must then bound the product below 1e-40.
                let mut dec = Decider::default();
                for (tangent, gradient) in tangents {
                    match dec.sign(&tangent.dot(&gradient)) {
                        Outcome::ProvenCoincident(_) => (),
                        Outcome::NotExact(near) => assert!(near.bound < r("1e-40")),
                        other => panic!("tangent leaves the cutting surface: {other:?}"),
                    }
                    assert!(matches!(
                        dec.sign(&tangent.norm_squared()),
                        Outcome::Proven(_)
                    ));
                }
            }
        }
    }
}

#[test]
fn oracle_tubes_respect_the_target_cap_end_caps_and_work_allowance() {
    for surface in SURFACES {
        let raw = candidates(fixture(&format!("{surface}-lofttight")))
            .pop()
            .unwrap();
        let tight = on_border(raw.clone());
        let default = on_border(
            candidates(fixture(&format!("{surface}-loft")))
                .pop()
                .unwrap(),
        );
        for name in ["raw", "tight", "default", "shifted", "short"] {
            let mut candidate = match name {
                "raw" => raw.clone(),
                "default" => default.clone(),
                _ => tight.clone(),
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
                let section = patch.sections(&cutter(surface, "loft"), &mut dec).unwrap();
                section.branches()[0].certify_candidate(&candidate, &mut dec)
            });
            println!(
                "WORK tube {surface} {name}: {} operation units, {} pieces, exhausted={}",
                work.used,
                candidate.pieces.len(),
                work.exhausted
            );
            assert!(
                !work.exhausted,
                "the tube exceeded the unchanged 20M allowance"
            );
            match (name, result) {
                ("tight", Ok(proof)) => {
                    assert!(proof.bound_mm() <= &r("1/10000"));
                    assert_eq!(proof.cap_mm(), &r("1/1000"));
                    assert!(proof.transverse_slope_lower() > &Rational::zero());
                    assert!(!proof.definition().is_closed());
                    assert!(matches!(
                        proof.proof(),
                        catallax_numbers::GeometricProof::Approximation { .. }
                    ));
                    println!(
                        "BOUND tube {surface} tight: U={} mm ({:.3e}), rho={}, speed={}, decision_work={}",
                        proof.bound_mm(),
                        proof.bound_mm().to_f64_nearest(),
                        proof.half_width(),
                        proof.speed_bound(),
                        proof.work()
                    );
                }
                ("tight", Err(error)) => panic!("{surface} tight: {} {}", error.code, error.detail),
                ("default", Ok(proof)) => {
                    // Accepted only under the hard cap, never under the target.
                    assert!(proof.bound_mm() > &r("1/10000"));
                    assert!(proof.bound_mm() <= &r("1/1000"));
                    println!(
                        "BOUND tube {surface} default: U={} mm ({:.3e}), rho={}",
                        proof.bound_mm(),
                        proof.bound_mm().to_f64_nearest(),
                        proof.half_width()
                    );
                }
                ("raw", Ok(proof)) => {
                    // The cylinder trace already ends exactly on both borders.
                    assert_eq!(surface, "cylinder");
                    assert!(proof.bound_mm() <= &r("1/10000"));
                }
                ("raw", Err(refusal)) => {
                    // The sphere trace ends 4.3e-19 off the u = 0 border.
                    assert_eq!(surface, "sphere");
                    assert_eq!(refusal.code, "intersect.patch-sphere-candidate-endpoint");
                }
                (_, Ok(proof)) => panic!("{surface} {name} accepted with U={}", proof.bound_mm()),
                (_, Err(refusal)) => assert_eq!(
                    refusal.code,
                    if name == "short" {
                        format!("intersect.patch-{surface}-candidate-endpoint")
                    } else {
                        format!("intersect.patch-{surface}-candidate-cap")
                    }
                ),
            }
        }
    }
}

#[test]
fn same_exact_answers_with_one_two_and_four_workers() {
    let cases = SURFACES
        .iter()
        .flat_map(|s| CUTS.iter().map(move |c| (*s, *c)))
        .collect::<Vec<_>>();
    let mut baseline = None;
    for workers in [1, 2, 4] {
        let mut rows = std::thread::scope(|scope| {
            let jobs = (0..workers)
                .map(|worker| {
                    let cases = &cases;
                    scope.spawn(move || {
                        cases
                            .iter()
                            .enumerate()
                            .skip(worker)
                            .step_by(workers)
                            .map(|(index, &(surface, cut))| {
                                let (answer, work) = with_operation_budget(20_000_000, || {
                                    answer(surface, cut, &mut Decider::default())
                                });
                                assert!(
                                    !work.exhausted,
                                    "{surface}-{cut}, {workers} workers: {} units",
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
            .map(|row| (format!("{}-{}", cases[row.0].0, cases[row.0].1), row.2))
            .collect::<Vec<_>>();
        println!("WORK workers={workers}: {units:?}");
        if let Some(baseline) = &baseline {
            assert_eq!(&answers, baseline);
        } else {
            baseline = Some(answers);
        }
    }
}
