//! Complete cone and ring-torus cuts, including designer contacts and copies.
use catallax_intersect::{
    Component, Contact, CurveForm, Cyclic, SurfacePair, curve_form, curve_is_closed, cyclic_order,
    intersect_surfaces, on_curve, order_on_curve, supports, tangent_at,
};
use catallax_numbers::{Decider, Outcome, Rational, Real, StrictOrder, with_operation_budget};
use catallax_shapes::{
    BSplineSurface, Cone, Curve, GeomStore, Point3, PointDef, Pose, RationalPatch, Rotation,
    SectionCandidate, SectionHits, SectionSurface, Surface, Torus, Vector3,
};

#[allow(dead_code)]
mod patch_fixtures;
use patch_fixtures::{equal, p, r};

const FAMILIES: [&str; 2] = ["cone", "torus"];
const CUTS: [&str; 6] = ["open", "axial", "loop", "tiny", "tangent", "above"];

fn expected_branches(family: &str, cut: &str) -> usize {
    match (family, cut) {
        (_, "above") => 0,
        ("torus", "axial") => 2,
        _ => 1,
    }
}

fn axis(x: i64, y: i64, z: i64) -> Vector3 {
    Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}

fn copies() -> Vec<(Real, Pose)> {
    let turned = Pose::new(
        Rotation::from_quaternion(r("1"), r("2"), r("3"), r("4")).unwrap(),
        axis(7, -3, 5),
        &mut Decider::default(),
    )
    .unwrap();
    vec![
        (Real::one(), Pose::identity()),
        (
            Real::one(),
            Pose::new(
                Rotation::identity(),
                axis(7, -3, 5),
                &mut Decider::default(),
            )
            .unwrap(),
        ),
        (Real::one(), turned.clone()),
        (r("1/10").into(), turned.clone()),
        (r("10").into(), turned),
    ]
}

fn signature(family: &str, cut: &str, copy: usize) -> (String, u64) {
    let (answer, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let (scale, pose) = &copies()[copy];
        let patch = patch(cut, scale, pose, &mut dec);
        match patch.sections(&cutter(family, cut, scale, pose), &mut dec) {
            Err(error) => {
                assert_eq!(
                    cut, "tangent",
                    "{family}-{cut} copy={copy}: {} {}",
                    error.code, error.detail
                );
                assert_eq!(error.code, format!("intersect.patch-{family}-tangent"));
                let point =
                    pose.apply_point(&Point3::new(Real::zero(), Real::zero(), scale.clone()));
                assert_eq!(
                    error.detail,
                    format!(
                        "nonregular contact at ({}, {}, {}), parameters (1/2, 1/2)",
                        point.x.rational_ref().unwrap(),
                        point.y.rational_ref().unwrap(),
                        point.z.rational_ref().unwrap()
                    )
                );
                error.code.to_owned()
            }
            Ok(section) => {
                assert_eq!(section.branches().len(), expected_branches(family, cut));
                let search = section.search().unwrap();
                let mut answer = format!(
                    "{} {} {}",
                    search.arc_cells, search.empty_cells, search.depth
                );
                for branch in section.branches() {
                    answer.push_str(&format!(" closed={}", branch.is_closed()));
                    for end in branch.parameter_endpoints() {
                        for value in end {
                            let b = dec.enclose_budgeted_at(&value, 64).unwrap();
                            answer.push_str(&format!(" [{},{}]", b.lo(), b.hi()));
                        }
                    }
                    if cut == "open" {
                        for end in branch.endpoints() {
                            let local = pose.inverse().apply_point(&end);
                            let local =
                                Point3::new(&local.x / scale, &local.y / scale, &local.z / scale);
                            let x = if local.x.rational_ref().unwrap().sign()
                                == catallax_numbers::Sign::Negative
                            {
                                "-3/5"
                            } else {
                                "3/5"
                            };
                            equal(&local, &p(x, "4/5", "2"), &mut dec);
                        }
                    }
                }
                answer
            }
        }
    });
    assert!(
        !work.exhausted,
        "{family}-{cut} copy={copy}: {} units",
        work.used
    );
    println!("WORK copy {family}-{cut} copy={copy}: {} units", work.used);
    (answer, work.used)
}

#[test]
fn copies_and_one_two_four_workers_keep_the_same_answers() {
    let cases = FAMILIES
        .into_iter()
        .flat_map(|f| {
            CUTS.into_iter()
                .flat_map(move |c| (0..copies().len()).map(move |i| (f, c, i)))
        })
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
                            .map(|(index, &(f, c, i))| (index, signature(f, c, i)))
                            .collect::<Vec<_>>()
                    })
                })
                .collect::<Vec<_>>();
            jobs.into_iter()
                .flat_map(|job| job.join().unwrap())
                .collect::<Vec<_>>()
        });
        rows.sort_by_key(|row| row.0);
        for group in rows.chunks(copies().len()) {
            for row in group {
                assert_eq!(row.1.0, group[0].1.0);
            }
        }
        println!(
            "WORK workers={workers}: {:?}",
            rows.iter()
                .map(|row| (cases[row.0], row.1.1))
                .collect::<Vec<_>>()
        );
        if let Some(baseline) = &baseline {
            assert_eq!(&rows, baseline)
        } else {
            baseline = Some(rows)
        }
    }
}

#[test]
fn pair_queries_keep_exact_branches_and_both_orientations() {
    for family in FAMILIES {
        for cut in CUTS {
            for swap in [false, true] {
                let (tangents, work) = with_operation_budget(20_000_000, || {
                    let mut dec = Decider::default();
                    let patch = patch(cut, &Real::one(), &Pose::identity(), &mut dec);
                    let analytic = cutter(family, cut, &Real::one(), &Pose::identity()).surface();
                    let patch_surface = Surface::BSpline(patch.into());
                    assert!(supports(&analytic, &patch_surface));
                    assert!(supports(&patch_surface, &analytic));
                    let mut geom = GeomStore::new();
                    let a = geom.add_surface(analytic.clone());
                    let b = geom.add_surface(patch_surface);
                    let result = if swap {
                        intersect_surfaces(&mut geom, b, a, &mut dec)
                    } else {
                        intersect_surfaces(&mut geom, a, b, &mut dec)
                    };
                    if cut == "tangent" {
                        assert_eq!(
                            result.unwrap_err().code,
                            format!("intersect.patch-{family}-tangent")
                        );
                        return Vec::new();
                    }
                    let answer = result.unwrap();
                    assert!(answer.is_exact());
                    if cut == "above" {
                        assert_eq!(answer.value, SurfacePair::Disjoint);
                        return Vec::new();
                    }
                    let SurfacePair::Meets(components) = answer.value else {
                        panic!("expected complete components")
                    };
                    assert_eq!(components.len(), expected_branches(family, cut));
                    let mut tangents = Vec::new();
                    for component in components {
                        let Component::Curve {
                            curve,
                            contact,
                            singular,
                        } = component
                        else {
                            panic!("expected regular curve")
                        };
                        assert_eq!(contact, Contact::Crossing);
                        assert!(singular.is_empty());
                        let Curve::PatchSection(branch) = geom.curve(curve) else {
                            panic!("expected retained exact definition")
                        };
                        let branch = branch.clone();
                        assert_eq!(branch.follows_surface_first(&mut dec).unwrap(), !swap);
                        assert!(matches!(
                            curve_form(&geom, curve, &mut dec).unwrap(),
                            CurveForm::PatchSection(_)
                        ));
                        assert_eq!(
                            curve_is_closed(&geom, curve, &mut dec).unwrap(),
                            cut != "open"
                        );
                        let [start, end] = branch.endpoints();
                        let p = geom.add_point(PointDef::Explicit(start.clone()));
                        let q = geom.add_point(PointDef::Explicit(end));
                        assert!(on_curve(&geom, curve, p, &mut dec).is_coincident());
                        let outside = branch.charts()[0]
                            .point_at(&[Real::zero(), Real::zero()], &mut dec)
                            .unwrap();
                        let outside = geom.add_point(PointDef::Explicit(outside));
                        assert!(matches!(
                            on_curve(&geom, curve, outside, &mut dec),
                            Outcome::Proven(_)
                        ));
                        if cut == "open" {
                            assert!(matches!(
                                order_on_curve(&geom, curve, p, q, &mut dec),
                                Outcome::Proven(StrictOrder::Less)
                            ));
                            assert!(matches!(
                                order_on_curve(&geom, curve, q, p, &mut dec),
                                Outcome::Proven(StrictOrder::Greater)
                            ));
                        } else {
                            assert_eq!(
                                branch
                                    .certify_candidate(
                                        &SectionCandidate { pieces: Vec::new() },
                                        &mut dec
                                    )
                                    .unwrap_err()
                                    .code,
                                format!("intersect.patch-{family}-candidate-domain")
                            );
                            let ids = branch
                                .witnesses()
                                .take(3)
                                .map(|p| geom.add_point(PointDef::Explicit(p)))
                                .collect::<Vec<_>>();
                            assert!(matches!(
                                cyclic_order(&geom, curve, [ids[0], ids[1], ids[2]], &mut dec),
                                Outcome::Proven(Cyclic::Forward)
                            ));
                            assert!(matches!(
                                cyclic_order(&geom, curve, [ids[0], ids[2], ids[1]], &mut dec),
                                Outcome::Proven(Cyclic::Backward)
                            ));
                        }
                        let raw = geom.add_curve(Curve::Intersection(
                            catallax_shapes::IntersectionCurve {
                                surfaces: [a, b],
                                witness: p,
                            },
                        ));
                        assert_eq!(
                            curve_form(&geom, raw, &mut dec).unwrap_err().code,
                            format!("intersect.patch-{family}-selector")
                        );
                        let gradient = analytic.gradient(&start).unwrap();
                        for tangent in tangent_at(&geom, curve, p, &mut dec).unwrap() {
                            tangents.push((tangent, gradient.clone()));
                        }
                    }
                    tangents
                });
                assert!(
                    !work.exhausted,
                    "{family}-{cut} swap={swap}: {} units",
                    work.used
                );
                println!("WORK pair {family}-{cut} swap={swap}: {} units", work.used);
                for (tangent, gradient) in tangents {
                    let mut dec = Decider::default();
                    let dot = dec
                        .enclose_budgeted_at(&tangent.dot(&gradient), 192)
                        .unwrap();
                    assert!(dot.contains_zero());
                    assert!(dot.lo().abs().max(dot.hi().abs()) < r("1e-40"));
                    let norm = dec
                        .enclose_budgeted_at(&tangent.norm_squared(), 192)
                        .unwrap();
                    assert!(norm.lo() > Rational::zero());
                }
            }
        }
    }
}

fn patch(cut: &str, scale: &Real, pose: &Pose, dec: &mut Decider) -> RationalPatch {
    let open = cut == "open";
    let mut controls = Vec::new();
    if open {
        for (x, z) in [("-3/5", "2"), ("0", "191/100"), ("3/5", "2")] {
            for y in ["7/10", "3/2"] {
                controls.push(p(x, y, z));
            }
        }
    } else {
        for (x, x2) in [("-1/4", "1/16"), ("0", "-1/16"), ("1/4", "1/16")] {
            for (y, y2) in [("-1/4", "1/16"), ("0", "-1/16"), ("1/4", "1/16")] {
                controls.push(Point3::new(
                    r(x).into(),
                    r(y).into(),
                    (r("1") - r(x2) - r(y2)).into(),
                ));
            }
        }
    }
    let controls = controls
        .into_iter()
        .map(|p| pose.apply_point(&Point3::new(&p.x * scale, &p.y * scale, &p.z * scale)))
        .collect::<Vec<_>>();
    let count = controls.len();
    let u = ["0", "0", "0", "1", "1", "1"].map(r).to_vec();
    let v = if open {
        ["0", "0", "1", "1"].map(r).to_vec()
    } else {
        u.clone()
    };
    RationalPatch::new(
        BSplineSurface::new(
            [2, if open { 1 } else { 2 }],
            [u, v],
            [3, if open { 2 } else { 3 }],
            controls,
            vec![Real::one(); count],
            dec,
        )
        .unwrap(),
        dec,
    )
    .unwrap()
}

fn cutter(family: &str, cut: &str, scale: &Real, pose: &Pose) -> SectionSurface {
    let low = match cut {
        "loop" => r("99/100"),
        "tiny" => r("999999/1000000"),
        "above" => r("1000001/1000000"),
        _ => r("1"),
    };
    let point =
        |p: Point3| pose.apply_point(&Point3::new(&p.x * scale, &p.y * scale, &p.z * scale));
    let surface = if family == "cone" {
        Surface::Cone(Cone {
            apex: point(if cut == "open" || cut == "axial" {
                p("0", "0", "0")
            } else {
                Point3::new(-Real::one(), Real::zero(), low.into())
            }),
            axis: pose.apply_vector(&if cut == "open" || cut == "axial" {
                axis(0, 0, 1)
            } else {
                axis(1, 0, 1)
            }),
            slope: if cut == "open" {
                r("1/2").into()
            } else if cut == "axial" {
                r("1/10").into()
            } else {
                Real::one()
            },
            ref_dir: pose.apply_vector(&axis(0, 1, 0)),
        })
    } else {
        Surface::Torus(Torus {
            center: point(if cut == "open" {
                p("0", "0", "2")
            } else if cut == "axial" {
                p("0", "0", "1")
            } else {
                Point3::new(Real::one(), Real::zero(), (low + r("1/4")).into())
            }),
            axis: pose.apply_vector(&axis(0, 0, 1)),
            major: scale
                * &if cut == "open" {
                    r("3/4").into()
                } else if cut == "axial" {
                    r("1/8").into()
                } else {
                    Real::one()
                },
            minor: scale * &Real::from_rational(r(if cut == "axial" { "1/32" } else { "1/4" })),
            ref_dir: pose.apply_vector(&axis(1, 0, 0)),
        })
    };
    SectionSurface::of(&surface).unwrap()
}

fn candidate(family: &str) -> SectionCandidate {
    let graph = |u: f64| {
        let x = -0.6 + 1.2 * u;
        let z = 2.0 + (x * x - 0.36) / 8.0;
        let dz = 0.3 * x;
        let (radius, dr) = if family == "cone" {
            (z / 2.0, dz / 2.0)
        } else {
            let h = z - 2.0;
            let b = (0.0625 - h * h).sqrt();
            (0.75 + b, -h * dz / b)
        };
        let y = (radius * radius - x * x).sqrt();
        ((y - 0.7) / 0.8, (radius * dr - 1.2 * x) / (y * 0.8))
    };
    let real = |v| Rational::from_f64_exact(v).unwrap();
    SectionCandidate {
        pieces: (0..16)
            .map(|i| {
                let a = Rational::from_i64(i) / r("16");
                let b = Rational::from_i64(i + 1) / r("16");
                let (v, dv) = graph(a.to_f64_nearest());
                let (w, dw) = graph(b.to_f64_nearest());
                [
                    [a.clone(), if i == 0 { r("1/8") } else { real(v) }],
                    [&a + &r("1/48"), real(v + dv / 48.0)],
                    [&b - &r("1/48"), real(w - dw / 48.0)],
                    [b, if i == 15 { r("1/8") } else { real(w) }],
                ]
            })
            .collect(),
    }
}

#[test]
fn candidate_tubes_are_bound_to_the_complete_branch() {
    for family in FAMILIES {
        for kind in ["tight", "shifted", "short"] {
            let mut candidate = candidate(family);
            if kind == "shifted" {
                for piece in &mut candidate.pieces {
                    for p in piece {
                        p[1] = &p[1] + &r("1/100");
                    }
                }
            }
            if kind == "short" {
                candidate.pieces[0][0][0] = r("1/1000");
            }
            let (result, work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let patch = patch("open", &Real::one(), &Pose::identity(), &mut dec);
                let section = patch
                    .sections(
                        &cutter(family, "open", &Real::one(), &Pose::identity()),
                        &mut dec,
                    )
                    .unwrap();
                section.branches()[0].certify_candidate(&candidate, &mut dec)
            });
            assert!(!work.exhausted, "{family}-{kind}: {} units", work.used);
            println!("WORK tube {family}-{kind}: {} units", work.used);
            if kind == "tight" {
                let proof = result.unwrap();
                assert!(proof.bound_mm() <= &r("1/10000"));
                assert_eq!(proof.cap_mm(), &r("1/1000"));
                assert!(proof.transverse_slope_lower() > &Rational::zero());
                assert!(matches!(
                    proof.proof(),
                    catallax_numbers::GeometricProof::Approximation { .. }
                ));
                println!(
                    "BOUND {family}: U={} mm; rho={}",
                    proof.bound_mm(),
                    proof.half_width()
                );
            } else {
                assert_eq!(
                    result.unwrap_err().code,
                    format!(
                        "intersect.patch-{family}-candidate-{}",
                        if kind == "short" { "endpoint" } else { "cap" }
                    )
                );
            }
        }
    }
}

#[test]
fn cone_reductions_cannot_admit_the_opposite_nappe_or_apex_boundary() {
    let mut dec = Decider::default();
    let patch = patch("open", &Real::one(), &Pose::identity(), &mut dec);
    let SectionSurface::Cone(mut cone) = cutter("cone", "open", &Real::one(), &Pose::identity())
    else {
        unreachable!()
    };
    cone.axis = axis(0, 0, -1);
    let (answer, work) = with_operation_budget(20_000_000, || {
        patch.sections(&SectionSurface::Cone(cone.clone()), &mut dec)
    });
    let answer = answer.unwrap();
    assert_eq!(answer.branches().len(), 0);
    assert_eq!(answer.curves().count(), 0);
    assert!(answer.points().is_empty());
    println!("WORK cone-opposite-nappe: {} units", work.used);
    cone.axis = axis(0, 0, 1);
    cone.apex = p("0", "0", "2");
    let (answer, work) = with_operation_budget(20_000_000, || {
        patch.sections(&SectionSurface::Cone(cone), &mut dec)
    });
    assert_eq!(answer.unwrap_err().code, "intersect.patch-cone-nappe");
    assert!(!work.exhausted);
    println!("WORK cone-apex-boundary: {} units", work.used);
}

#[test]
fn cone_third_surface_queries_keep_the_retained_nappe() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let patch = patch("open", &Real::one(), &Pose::identity(), &mut dec);
        let cutter = cutter("cone", "open", &Real::one(), &Pose::identity());
        let section = patch.sections(&cutter, &mut dec).unwrap();
        let branch = &section.branches()[0];
        assert!(matches!(
            branch.surface_hits(&cutter, &mut dec).unwrap(),
            SectionHits::Contained
        ));
        let SectionSurface::Cone(mut opposite) = cutter else {
            unreachable!()
        };
        opposite.axis = axis(0, 0, -1);
        let hits = branch
            .surface_hits(&SectionSurface::Cone(opposite), &mut dec)
            .unwrap();
        assert!(
            matches!(hits, SectionHits::Points(points) if points.is_empty()),
            "the opposite cone nappe must be excluded"
        );
    });
    assert!(!work.exhausted);
    println!("WORK cone-third-surface-nappes: {} units", work.used);
}

fn planar_patch(points: [Point3; 4], dec: &mut Decider) -> RationalPatch {
    RationalPatch::new(
        BSplineSurface::new(
            [1, 1],
            [
                ["0", "0", "1", "1"].map(r).to_vec(),
                ["0", "0", "1", "1"].map(r).to_vec(),
            ],
            [2, 2],
            points.to_vec(),
            vec![Real::one(); 4],
            dec,
        )
        .unwrap(),
        dec,
    )
    .unwrap()
}

fn contact_answer(family: &str, boundary: bool, copy: usize) -> (String, String, u64) {
    let (result, work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let (points, mut surface) = if family == "cone" {
            let y = if boundary { "0" } else { "-1" };
            (
                [
                    p("1", y, "1"),
                    p("1", "1", "1"),
                    p("2", y, "2"),
                    p("2", "1", "2"),
                ],
                SectionSurface::Cone(Cone {
                    apex: p("0", "0", "0"),
                    axis: axis(0, 0, 1),
                    slope: Real::one(),
                    ref_dir: axis(1, 0, 0),
                }),
            )
        } else {
            let points = if boundary {
                [
                    p("5/4", "0", "-1/2"),
                    p("5/4", "0", "1/2"),
                    p("5/4", "1/2", "-1/2"),
                    p("5/4", "1/2", "1/2"),
                ]
            } else {
                [
                    p("3/4", "-1/4", "1/4"),
                    p("3/4", "1/4", "1/4"),
                    p("5/4", "-1/4", "1/4"),
                    p("5/4", "1/4", "1/4"),
                ]
            };
            (
                points,
                SectionSurface::Torus(Torus {
                    center: p("0", "0", "0"),
                    axis: axis(0, 0, 1),
                    major: Real::one(),
                    minor: r("1/4").into(),
                    ref_dir: axis(1, 0, 0),
                }),
            )
        };

        let (scale, pose) = &copies()[copy];
        let points = points
            .map(|p| pose.apply_point(&Point3::new(&p.x * scale, &p.y * scale, &p.z * scale)));
        if let SectionSurface::Torus(torus) = &mut surface {
            torus.major = &torus.major * scale;
            torus.minor = &torus.minor * scale;
        }
        let surface = SectionSurface::of(&surface.surface().transformed(pose)).unwrap();
        planar_patch(points, &mut dec).sections(&surface, &mut dec)
    });
    assert!(!work.exhausted);
    let error = result.unwrap_err();
    assert_eq!(
        error.code,
        format!(
            "intersect.patch-{family}-{}",
            if family == "torus" && boundary {
                "nonregular"
            } else {
                "tangent"
            }
        )
    );
    println!(
        "WORK {family}-flush boundary={boundary} copy={copy}: {} units; {}",
        work.used, error.code
    );
    (error.code.to_owned(), error.detail, work.used)
}

#[test]
fn flush_and_boundary_contacts_refuse_on_copies_and_at_every_worker_count() {
    let cases = FAMILIES
        .into_iter()
        .flat_map(|f| {
            [false, true]
                .into_iter()
                .flat_map(move |b| (0..copies().len()).map(move |i| (f, b, i)))
        })
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
                            .map(|(index, &(f, b, i))| (index, contact_answer(f, b, i)))
                            .collect::<Vec<_>>()
                    })
                })
                .collect::<Vec<_>>();
            jobs.into_iter()
                .flat_map(|job| job.join().unwrap())
                .collect::<Vec<_>>()
        });
        rows.sort_by_key(|row| row.0);
        for group in rows.chunks(copies().len()) {
            for row in group {
                assert_eq!(row.1.0, group[0].1.0);
            }
        }
        if let Some(baseline) = &baseline {
            assert_eq!(&rows, baseline)
        } else {
            baseline = Some(rows)
        }
    }
}

#[test]
fn work_exhaustion_never_yields_a_partial_cut() {
    for family in FAMILIES {
        let mut dec = Decider::default();
        let patch = patch("open", &Real::one(), &Pose::identity(), &mut dec);
        let (result, work) = with_operation_budget(1, || {
            patch.sections(
                &cutter(family, "open", &Real::one(), &Pose::identity()),
                &mut Decider::default(),
            )
        });
        assert!(work.exhausted);
        assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    }
}

#[test]
fn complete_search_proves_open_cuts_loops_empty_cuts_and_named_contacts() {
    for family in FAMILIES {
        for cut in CUTS {
            let (result, work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let patch = patch(cut, &Real::one(), &Pose::identity(), &mut dec);
                let section = patch.sections(
                    &cutter(family, cut, &Real::one(), &Pose::identity()),
                    &mut dec,
                )?;
                assert_eq!(section.curves().count(), 0);
                assert!(section.points().is_empty());
                assert_eq!(section.branches().len(), expected_branches(family, cut));
                for branch in section.branches() {
                    assert_eq!(branch.is_closed(), cut != "open");
                    if cut == "open" {
                        let ends = branch.endpoints();
                        for end in ends {
                            let x = if end.x.rational_ref().unwrap().sign()
                                == catallax_numbers::Sign::Negative
                            {
                                "-3/5"
                            } else {
                                "3/5"
                            };
                            equal(&end, &p(x, "4/5", "2"), &mut dec);
                        }
                    }
                }
                Ok::<_, catallax_numbers::Refusal>(section)
            });
            assert!(!work.exhausted, "{family}-{cut}: {}", work.used);
            match result {
                Ok(section) => println!(
                    "WORK {family}-{cut}: {} units; {:?}",
                    work.used,
                    section.search()
                ),
                Err(error) => {
                    assert_eq!(
                        cut, "tangent",
                        "{family}-{cut}: {} {} ({} units)",
                        error.code, error.detail, work.used
                    );
                    assert_eq!(error.code, format!("intersect.patch-{family}-tangent"));
                    assert!(error.detail.contains("(0, 0, 1), parameters (1/2, 1/2)"));
                    println!(
                        "WORK {family}-{cut}: {} units; {} {}",
                        work.used, error.code, error.detail
                    );
                }
            }
        }
    }
}
