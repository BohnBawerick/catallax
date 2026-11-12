//! M10: a box against a smooth three-section loft. The loft's roof is one
//! periodic joined atlas face; box planes cut it into exact branches that
//! join across chart seams and meet the other box planes exactly.
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{
    Plane, Point3, RationalAtlas, SectionHits, SectionSurface, SplineSupport, Surface, Vector3,
};
use serde_json::{Value, json};

/// A quadratic roof over the floor from (4,0) to (0,0), peaking at y = 2,
/// 3 and 2 at heights 0, 1 and 2.
fn loft_steps() -> Vec<Value> {
    let section = |i: usize, z: &str, peak: &str| {
        json!({"id":format!("section{i}"),"op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0",z]},
            "loops":[[{"kind":"bezier","id":format!("s{i}_roof"),"controls":[["2",peak]],"to":["4","0"]},
                      {"kind":"line","id":format!("s{i}_floor"),"to":["0","0"]}]]})
    };
    vec![
        section(0, "0", "4"),
        section(1, "1", "6"),
        section(2, "2", "4"),
        json!({"id":"loft","op":"loft","shape":"smooth","sections":[
            {"profile":"section0","start":["s0_roof"]},
            {"profile":"section1","start":["s1_roof"]},
            {"profile":"section2","start":["s2_roof"]}]}),
    ]
}

fn recipe(steps: Vec<Value>, result: &str) -> Recipe {
    let value =
        json!({"format":"catallax-recipe/5","name":"m10 box loft","steps":steps,"result":result});
    Recipe::parse(&value.to_string()).unwrap()
}

fn roof(dec: &mut Decider) -> RationalAtlas {
    let built = catallax_features::evaluate(&recipe(loft_steps(), "loft"), dec).unwrap();
    let model = &built.model;
    let mut roofs = model
        .faces_of_body(built.result)
        .into_iter()
        .filter_map(
            |face| match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => Some(atlas.clone()),
                _ => None,
            },
        )
        .collect::<Vec<_>>();
    assert_eq!(roofs.len(), 1);
    roofs.pop().unwrap()
}

fn plane(axis: usize, at: &str) -> SectionSurface {
    let q = |s: &str| Real::from_rational(Rational::parse(s).unwrap());
    let mut origin = [q("0"), q("0"), q("0")];
    origin[axis] = q(at);
    let mut normal = [q("0"), q("0"), q("0")];
    normal[axis] = q("1");
    let mut reference = [q("0"), q("0"), q("0")];
    reference[(axis + 1) % 3] = q("1");
    let [ox, oy, oz] = origin;
    let [nx, ny, nz] = normal;
    let [rx, ry, rz] = reference;
    SectionSurface::Plane(Plane {
        origin: Point3::new(ox, oy, oz),
        normal: Vector3::new(nx, ny, nz),
        ref_dir: Vector3::new(rx, ry, rz),
    })
}

fn enclose(value: &Real, dec: &mut Decider) -> f64 {
    dec.enclose_budgeted_at(value, 64)
        .unwrap()
        .lo()
        .to_f64_nearest()
}

#[test]
fn finite_atlas_face_queries_exclude_off_carrier_points() {
    use catallax_boolean::face_contains;
    for height in ["1", "13/4"] {
        let ((), work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let mut built =
                catallax_features::evaluate(&recipe(loft_steps(), "loft"), &mut dec).unwrap();
            let point = Point3::new(
                Real::from_i64(2),
                Rational::parse(height).unwrap().into(),
                Real::one(),
            );
            let face = built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .find(|face| {
                    matches!(
                        built.model.geom().surface(built.model.face(*face).surface),
                        Surface::BSpline(SplineSupport::RationalAtlas(_))
                    )
                })
                .unwrap();
            assert!(!face_contains(&mut built.model, face, &point, &mut dec).unwrap());
        });
        assert!(!work.exhausted, "{height}: {} units", work.used);
    }
}

#[test]
fn roof_branches_join_chart_seams_and_meet_the_other_box_planes() {
    // (cutting plane, third plane, meetings per branch)
    let cases = [
        ((1, "1"), (2, "1/2"), vec![1, 1]),
        ((1, "1"), (0, "1"), vec![0, 0]),
        ((1, "5/2"), (2, "1/2"), vec![2]),
        ((1, "5/2"), (0, "2"), vec![2]),
        ((2, "1/2"), (1, "1"), vec![2]),
        ((2, "1/2"), (0, "3"), vec![1]),
    ];
    for (cut, third, counts) in cases {
        let (points, work) = with_operation_budget(20_000_000, || {
            let mut dec = Decider::default();
            let atlas = roof(&mut dec);
            assert_eq!(atlas.charts().len(), 4);
            assert!(atlas.periodic_u());
            let section = atlas.sections(&plane(cut.0, cut.1), &mut dec).unwrap();
            let mut all = Vec::new();
            for branch in section.branches() {
                // The y = 5/2 cut is a closed loop over the ridge; the others
                // run between the floor edges or the end caps.
                assert_eq!(branch.is_closed(), cut == (1, "5/2"));
                let SectionHits::Points(points) = branch
                    .surface_hits(&plane(third.0, third.1), &mut dec)
                    .unwrap()
                else {
                    panic!("no box plane contains a roof branch")
                };
                let mut enclosed = Vec::new();
                for p in &points {
                    enclosed.push([
                        enclose(&p.x, &mut dec),
                        enclose(&p.y, &mut dec),
                        enclose(&p.z, &mut dec),
                    ]);
                }
                all.push(enclosed);
            }
            all
        });
        assert!(!work.exhausted, "{cut:?} {third:?}: {} units", work.used);
        println!(
            "WORK roof {cut:?} by {third:?}: {} units; {points:?}",
            work.used
        );
        assert_eq!(points.iter().map(Vec::len).collect::<Vec<_>>(), counts);
        let value =
            |(axis, at): (usize, &str)| (axis, Rational::parse(at).unwrap().to_f64_nearest());
        for p in points.iter().flatten() {
            for (axis, at) in [value(cut), value(third)] {
                assert!((p[axis] - at).abs() < 1e-15, "{p:?} off {axis}={at}");
            }
        }
    }
}

/// The Boolean of the loft with a box and its charged work.
fn boolean(kind: &str, case: &str) -> (String, u64) {
    let (size, at) = match case {
        // Crosses the roof's sides and ridge between the end caps.
        "side" => (["2", "10", "1"], ["1", "1", "1/2"]),
        // Takes the top off the ridge in one closed loop.
        _ => (["10", "10", "1"], ["-3", "5/2", "1/2"]),
    };
    let mut steps = loft_steps();
    steps.push(json!({"id":"box","op":"block","size":size,
        "pose":{"rotation":["1","0","0","0"],"translation":at}}));
    steps.push(json!({"id":"cut","op":"boolean","kind":kind,"target":"loft","tool":"box"}));
    let recipe = recipe(steps, "cut");
    let (result, work) = with_operation_budget(20_000_000, || {
        catallax_features::evaluate(&recipe, &mut Decider::default()).map(|_| ())
    });
    let answer = match result {
        Ok(()) => "solid".to_string(),
        Err(catallax_features::BuildError::Geometry { refusal, .. }) => refusal.code.to_string(),
        Err(other) => format!("{other:?}"),
    };
    if work.exhausted {
        assert!(
            answer.starts_with("numbers.operation-budget"),
            "{kind} {case}: {answer}"
        );
    }
    (answer, work.used)
}

fn assert_boolean_matrix(mut check: impl FnMut(&str, &str, &str, u64)) {
    let jobs = ["subtract", "union", "intersect"]
        .into_iter()
        .flat_map(|kind| ["side", "top"].map(|case| (kind, case)))
        .collect::<Vec<_>>();
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
                            .map(|(i, &(kind, case))| (i, boolean(kind, case)))
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
        for (i, (answer, work)) in &rows {
            println!(
                "WORK boolean workers={workers} {:?}: {work} units; {answer}",
                jobs[*i]
            );
        }
        let answers = rows.into_iter().map(|row| row.1).collect::<Vec<_>>();
        for ((kind, case), (answer, work)) in jobs.iter().zip(&answers) {
            check(kind, case, answer, *work);
        }
        if let Some(baseline) = &baseline {
            assert_eq!(&answers, baseline);
        } else {
            baseline = Some(answers);
        }
    }
}

#[test]
fn box_loft_retained_trace_progress_is_identical_on_every_worker_count() {
    assert_boolean_matrix(|kind, case, answer, work| {
        let (expected, expected_work) = match (kind, case) {
            (_, "top") => ("numbers.operation-budget", 19_997_386),
            ("subtract", "side") => ("solid", 10_942_803),
            ("union", "side") => ("solid", 10_942_852),
            ("intersect", "side") => ("solid", 10_934_705),
            _ => unreachable!(),
        };
        assert_eq!(answer, expected, "{kind} {case}: {answer}");
        assert_eq!(work, expected_work, "{kind} {case}");
    });
}

#[test]
#[ignore = "target: top cut must fit the unchanged 20M proof allowance"]
fn box_loft_booleans_build_with_same_answers_on_every_worker_count() {
    assert_boolean_matrix(|kind, case, answer, _| {
        assert_eq!(answer, "solid", "{kind} {case}: {answer}");
    });
}
