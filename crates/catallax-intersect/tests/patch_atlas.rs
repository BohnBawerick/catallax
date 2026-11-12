//! Sections of a joined atlas face: chart-seam joins, branch identity and
//! meetings with a third analytic surface.
use catallax_intersect::{
    Component, Contact, Cyclic, SurfacePair, cyclic_order, intersect_surfaces, on_curve,
    same_curve, tangent_at,
};
use catallax_numbers::{Decider, Outcome, Real, with_operation_budget};
use catallax_shapes::{
    Curve, GeomStore, Plane, Point3, PointDef, RationalAtlas, SectionHits, SectionSurface, Surface,
    Vector3,
};

#[allow(dead_code)]
mod patch_fixtures;
use patch_fixtures::{cutter, height_patch, r};

/// The bicubic hill of the plane slice, as one atlas with a chart per span.
fn hill_atlas(dec: &mut Decider) -> RationalAtlas {
    let source = height_patch(true, dec).carrier().clone();
    let mut atlases = RationalAtlas::partition(source, [false, false], dec).unwrap();
    assert_eq!(atlases.len(), 1);
    atlases.pop().unwrap()
}

fn plane(axis: usize, at: &str, sign: i64) -> SectionSurface {
    let mut origin = [Real::zero(), Real::zero(), Real::zero()];
    origin[axis] = r(at).into();
    let mut normal = [Real::zero(), Real::zero(), Real::zero()];
    normal[axis] = Real::from_i64(sign);
    let mut reference = [Real::zero(), Real::zero(), Real::zero()];
    reference[(axis + 1) % 3] = Real::one();
    let [ox, oy, oz] = origin;
    let [nx, ny, nz] = normal;
    let [rx, ry, rz] = reference;
    SectionSurface::Plane(Plane {
        origin: Point3::new(ox, oy, oz),
        normal: Vector3::new(nx, ny, nz),
        ref_dir: Vector3::new(rx, ry, rz),
    })
}

/// Cutting surfaces on the hill: the plane z = 3/10 and the plane slice's
/// cylinder and sphere cuts, at 0.3, one millionth below the top and at it.
fn cut(name: &str) -> SectionSurface {
    match name {
        "plane-hill" => plane(2, "3/10", 1),
        "plane-tiny" => plane(2, &(r("225/256") - r("1/1000000")).to_string(), 1),
        "plane-tangent" => plane(2, "225/256", 1),
        _ => {
            let (surface, cut) = name.split_once('-').unwrap();
            cutter(surface, cut)
        }
    }
}

const CUTS: [&str; 9] = [
    "plane-hill",
    "plane-tiny",
    "plane-tangent",
    "cylinder-hill",
    "cylinder-tiny",
    "cylinder-tangent",
    "sphere-hill",
    "sphere-tiny",
    "sphere-tangent",
];

fn enclose(value: &Real, dec: &mut Decider) -> f64 {
    dec.enclose_budgeted_at(value, 64)
        .unwrap()
        .lo()
        .to_f64_nearest()
}

#[test]
fn atlas_sections_join_chart_seams_and_match_the_single_patch() {
    for name in CUTS {
        let (single, _) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let patch = height_patch(true, &mut dec);
            patch.sections(&cut(name), &mut dec).map(|s| {
                s.branches()
                    .iter()
                    .map(|b| b.is_closed())
                    .collect::<Vec<_>>()
            })
        });
        // Atlas admission and search, measured on their own.
        let (joined, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let atlas = hill_atlas(&mut dec);
            atlas.sections(&cut(name), &mut dec).map(|s| {
                let charts = s
                    .branches()
                    .iter()
                    .map(|b| b.charts().len())
                    .collect::<Vec<_>>();
                let closed = s
                    .branches()
                    .iter()
                    .map(|b| b.is_closed())
                    .collect::<Vec<_>>();
                (closed, charts, s.search().clone())
            })
        });
        assert!(!work.exhausted, "{name}: {} units", work.used);
        println!("WORK atlas {name}: {} units; {:?}", work.used, joined);
        match (single, joined) {
            (Ok(single), Ok((closed, charts, _))) => {
                assert_eq!(single, closed, "{name}");
                assert!(charts.iter().all(|&c| c == 9));
            }
            (Err(a), Err(b)) => {
                assert!(name.ends_with("tangent"));
                assert_eq!(a.code, b.code);
                assert_eq!(a.detail, b.detail);
            }
            (a, b) => panic!("{name}: single {a:?} but atlas {b:?}"),
        }
    }
}

#[test]
fn atlas_pairs_answer_queries_and_identify_their_branches() {
    for name in ["plane-hill", "cylinder-hill", "sphere-hill", "sphere-tiny"] {
        for swap in [false, true] {
            let ((), work) = with_operation_budget(20_000_000, || {
                let mut dec = Decider::default();
                let atlas = hill_atlas(&mut dec);
                let mut geom = GeomStore::new();
                let a = geom.add_surface(cut(name).surface());
                let b = geom.add_surface(Surface::BSpline(atlas.into()));
                let (first, second) = if swap { (b, a) } else { (a, b) };
                let answer = intersect_surfaces(&mut geom, first, second, &mut dec).unwrap();
                assert!(answer.is_exact());
                let SurfacePair::Meets(components) = answer.value else {
                    panic!("{name}: expected one loop")
                };
                assert_eq!(components.len(), 1);
                let Component::Curve { curve, contact, .. } = &components[0] else {
                    panic!("expected a curve")
                };
                assert_eq!(*contact, Contact::Crossing);
                let Curve::PatchSection(branch) = geom.curve(*curve) else {
                    panic!("expected an exact branch")
                };
                let branch = branch.clone();
                assert!(branch.is_closed());
                assert_eq!(branch.follows_surface_first(&mut dec).unwrap(), !swap);
                let ids = branch
                    .witnesses()
                    .take(3)
                    .map(|p| geom.add_point(PointDef::Explicit(p)))
                    .collect::<Vec<_>>();
                assert!(on_curve(&geom, *curve, ids[0], &mut dec).is_coincident());
                assert!(matches!(
                    cyclic_order(&geom, *curve, [ids[0], ids[1], ids[2]], &mut dec),
                    Outcome::Proven(Cyclic::Forward)
                ));
                assert_eq!(
                    tangent_at(&geom, *curve, ids[0], &mut dec).unwrap().len(),
                    2
                );
                // The same section computed again is the same component; the
                // reversed copy too. A different cut on the same charts is not.
                let again = geom.add_curve(Curve::PatchSection(branch.reversed()));
                assert!(same_curve(&geom, *curve, again, &mut dec).is_coincident());
                let other = hill_atlas(&mut dec)
                    .sections(&plane(2, "1/2", 1), &mut dec)
                    .unwrap()
                    .branches()[0]
                    .clone();
                let other = geom.add_curve(Curve::PatchSection(other));
                assert!(matches!(
                    same_curve(&geom, *curve, other, &mut dec),
                    Outcome::Proven(_)
                ));
            });
            assert!(!work.exhausted, "{name} swap={swap}: {} units", work.used);
            println!("WORK atlas pair {name} swap={swap}: {} units", work.used);
        }
    }
}

/// The meetings of the z = 3/10 loop with a third surface, as enclosed
/// coordinates, or `None` when the loop lies on it.
fn hits(third: &SectionSurface, dec: &mut Decider) -> Option<Vec<[f64; 3]>> {
    let atlas = hill_atlas(dec);
    let section = atlas.sections(&plane(2, "3/10", 1), dec).unwrap();
    assert_eq!(section.branches().len(), 1);
    match section.branches()[0].surface_hits(third, dec).unwrap() {
        SectionHits::Contained => None,
        SectionHits::Points(points) => Some(
            points
                .iter()
                .map(|p| [enclose(&p.x, dec), enclose(&p.y, dec), enclose(&p.z, dec)])
                .collect(),
        ),
    }
}

#[test]
fn third_surface_meetings_are_exact_including_chart_seams() {
    // x(u) = 23/12 exactly on the chart seam u = 1, so those meetings are
    // cell ends. x = 5/2 and y = 5/2 cross the loop inside charts.
    for (name, third, count) in [
        ("x=5/2", plane(0, "5/2", 1), 2),
        ("y=5/2", plane(1, "5/2", 1), 2),
        ("x=23/12 on a seam", plane(0, "23/12", 1), 2),
        ("z=1 above", plane(2, "1", 1), 0),
        (
            "sphere touching z=3/10 inside the loop",
            cutter("sphere", "hill"),
            0,
        ),
    ] {
        let (result, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            hits(&third, &mut dec)
        });
        assert!(!work.exhausted, "{name}: {} units", work.used);
        println!("WORK hits {name}: {} units; {result:?}", work.used);
        let points = result.unwrap();
        assert_eq!(points.len(), count, "{name}");
        for p in &points {
            assert!((p[2] - 0.3).abs() < 1e-15, "{name}: {p:?}");
        }
    }
    // The cutting plane itself, with the opposite normal, contains the loop.
    let (contained, work) = with_operation_budget(20_000_000, || {
        hits(&plane(2, "3/10", -1), &mut Decider::default())
    });
    assert!(!work.exhausted);
    assert!(contained.is_none());
}

#[test]
fn same_atlas_answers_with_one_two_and_four_workers() {
    let jobs: Vec<(&str, SectionSurface)> = vec![
        ("plane-hill", cut("plane-hill")),
        ("sphere-tiny", cut("sphere-tiny")),
        ("cylinder-tangent", cut("cylinder-tangent")),
        ("hits x=5/2", plane(0, "5/2", 1)),
        ("hits x=23/12", plane(0, "23/12", 1)),
    ];
    let mut baseline = None;
    for workers in [1, 2, 4] {
        let mut rows = std::thread::scope(|scope| {
            let handles = (0..workers)
                .map(|worker| {
                    let jobs = &jobs;
                    scope.spawn(move || {
                        jobs.iter()
                            .enumerate()
                            .skip(worker)
                            .step_by(workers)
                            .map(|(index, (name, surface))| {
                                let (answer, work) = with_operation_budget(20_000_000, || {
                                    let mut dec = Decider::default();
                                    if name.starts_with("hits ") {
                                        format!("{:?}", hits(surface, &mut dec))
                                    } else {
                                        match hill_atlas(&mut dec).sections(surface, &mut dec) {
                                            Err(e) => format!("{} {}", e.code, e.detail),
                                            Ok(s) => format!(
                                                "{:?} {:?}",
                                                s.branches()
                                                    .iter()
                                                    .map(|b| b.is_closed())
                                                    .collect::<Vec<_>>(),
                                                s.search()
                                            ),
                                        }
                                    }
                                });
                                assert!(!work.exhausted);
                                (index, answer, work.used)
                            })
                            .collect::<Vec<_>>()
                    })
                })
                .collect::<Vec<_>>();
            handles
                .into_iter()
                .flat_map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        rows.sort_by_key(|row| row.0);
        let answers = rows
            .iter()
            .map(|row| (row.1.clone(), row.2))
            .collect::<Vec<_>>();
        println!(
            "WORK atlas workers={workers}: {:?}",
            rows.iter()
                .map(|row| (jobs[row.0].0, row.2))
                .collect::<Vec<_>>()
        );
        if let Some(baseline) = &baseline {
            assert_eq!(&answers, baseline);
        } else {
            baseline = Some(answers);
        }
    }
}

/// Two charts joined C1 at u = 0: z = u - e + 2u²/e on the left and
/// z = u - e on the right, with e = 1/1024, over y in [-1/2, 1/2]. The
/// section y = 1/10 meets z = 0 at u = -e and u = +e, one in each span.
fn close_cross_span_atlas(dec: &mut Decider) -> RationalAtlas {
    let e = r("1/1024");
    let z = [
        &(&r("-1") - &e) + &(&r("2") / &e),
        &r("-1/2") - &e,
        &r("1/2") - &e,
        &r("1") - &e,
    ];
    let x = ["-1", "-1/2", "1/2", "1"].map(r);
    let mut controls = Vec::new();
    for (x, z) in x.iter().zip(&z) {
        for y in ["-1/2", "1/2"].map(r) {
            controls.push(Point3::new(x.clone().into(), y.into(), z.clone().into()));
        }
    }
    let source = catallax_shapes::BSplineSurface::new(
        [2, 1],
        [
            ["-1", "-1", "-1", "0", "1", "1", "1"].map(r).to_vec(),
            ["0", "0", "1", "1"].map(r).to_vec(),
        ],
        [4, 2],
        controls,
        vec![Real::one(); 8],
        dec,
    )
    .unwrap();
    let mut atlases = RationalAtlas::partition(source, [false, false], dec).unwrap();
    assert_eq!(atlases.len(), 1);
    assert!(atlases[0].charts().len() >= 2);
    atlases.pop().unwrap()
}

#[test]
fn close_meetings_in_neighbouring_spans_are_both_kept_in_either_direction() {
    for reverse in [false, true] {
        let (xs, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let atlas = close_cross_span_atlas(&mut dec);
            let section = atlas.sections(&plane(1, "1/10", 1), &mut dec).unwrap();
            assert_eq!(section.branches().len(), 1);
            let branch = &section.branches()[0];
            let branch = if reverse {
                branch.reversed()
            } else {
                branch.clone()
            };
            let SectionHits::Points(points) =
                branch.surface_hits(&plane(2, "0", 1), &mut dec).unwrap()
            else {
                panic!("the branch is not on z = 0")
            };
            points
                .iter()
                .map(|p| enclose(&p.x, &mut dec))
                .collect::<Vec<_>>()
        });
        assert!(!work.exhausted, "reverse={reverse}: {} units", work.used);
        println!(
            "WORK close cross-span reverse={reverse}: {} units; {xs:?}",
            work.used
        );
        let e = 1.0 / 1024.0;
        let mut sorted = xs.clone();
        sorted.sort_by(f64::total_cmp);
        assert_eq!(sorted.len(), 2, "reverse={reverse}: {xs:?}");
        assert!(
            (sorted[0] + e).abs() < 1e-15 && (sorted[1] - e).abs() < 1e-15,
            "{xs:?}"
        );
        // Traversal order: the reversed branch meets them the other way round.
        assert_eq!(xs[0] < xs[1], xs == sorted);
    }
}
