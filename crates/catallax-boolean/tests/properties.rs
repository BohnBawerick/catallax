//! Exact generated families, placement invariance, determinism and budgets.
mod common;
use common::*;
use std::time::{Duration, Instant};

fn refusal_placement(
    id: &'static str,
    op: BoolOp,
    scale: &str,
    quaternion: [i64; 4],
    translation: [&str; 3],
) {
    let [w, x, y, z] = quaternion.map(r);
    let pose = Pose::new(
        Rotation::from_quaternion(w, x, y, z).unwrap(),
        (RVector3::new(q(translation[0]), q(translation[1]), q(translation[2]))).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    let scale = q(scale);
    let c = case(id);
    let mut model = Model::new();
    let pa = c.a.transformed(&scale, &pose);
    let pb = c.b.transformed(&scale, &pose);
    let a = pa.build(&mut model);
    let b = pb.build(&mut model);
    let mut dec = Decider::default();
    let answer = boolean(&mut model, a, b, op, &mut dec).unwrap();
    assert!(answer.marks.is_empty());
    topology_rules(&model, answer.body);
    let (_, _, expected) = expectations(id)
        .into_iter()
        .find(|(operation, reverse, _)| *operation == op && !reverse)
        .unwrap();
    let Expected::Counts(want) = expected else {
        panic!("placement must have a manifold result");
    };
    check_counts(&model, answer.body, &want);
    for (point, want) in witnesses(id, op) {
        query(
            &mut model,
            answer.body,
            &moved_point(&point, &scale, &pose),
            want,
            &mut dec,
        );
    }
    let point = moved_point(&p("1/7", "2/7", "3/7"), &scale, &pose);
    if let Some(want) = set_rule(op, pa.locate(&point), pb.locate(&point)) {
        query(&mut model, answer.body, &point, want, &mut dec);
    }
}

#[test]
fn crossing_cylinders_tiny_rotated_union() {
    refusal_placement(
        "23",
        BoolOp::Unite,
        "1e-20",
        [40, -42, -1, 45],
        ["434139/7", "-354933/5", "24790"],
    );
}

#[test]
fn crossing_cylinders_rotated_intersection() {
    refusal_placement(
        "23",
        BoolOp::Intersect,
        "1/10",
        [42, 66, -43, 64],
        ["89983/4", "-34922", "-877620/13"],
    );
}

#[test]
fn crossing_cylinders_tiny_intersection() {
    refusal_placement("23", BoolOp::Intersect, "1e-30", [1, 0, 0, 0], ["0"; 3]);
}

#[test]
fn near_tangent_cylinders_tiny_union() {
    refusal_placement("29-bite", BoolOp::Unite, "1e-30", [1, 0, 0, 0], ["0"; 3]);
}

#[test]
fn tiny_rotated_block_cylinder_cap() {
    refusal_placement(
        "13-bite",
        BoolOp::Intersect,
        "1e-28",
        [78, 44, -89, -60],
        ["794523/17", "-492109/15", "-667617/17"],
    );
}

#[test]
fn coincident_sections_preserve_a_cavity_and_curve_contact() {
    let global = Pose::new(
        Rotation::from_quaternion(r(1), r(1), r(1), r(1)).unwrap(),
        (v(0, 0, 0)).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    let pose = |[x, y, z]: [&str; 3]| {
        Pose::new(
            global.rotation.clone(),
            (RVector3::new(q(x), q(y), q(z))).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap()
    };
    let primitives = [
        Primitive {
            shape: Shape::Block([r(6), r(6), r(6)]),
            pose: pose(["-3", "-3", "-3"]),
        },
        Primitive {
            shape: Shape::Cylinder(r(1), r(2)),
            pose: pose(["-1", "4", "0"]),
        },
        Primitive {
            shape: Shape::Ball(r(2)),
            pose: pose(["1/5", "1/2", "1/3"]),
        },
        Primitive {
            shape: Shape::Block([r(6), r(6), r(6)]),
            pose: pose(["-14/5", "-2", "-3"]),
        },
    ];
    let mut model = Model::new();
    let [a, b, cavity, clip] = primitives.each_ref().map(|p| p.build(&mut model));
    let mut dec = Decider::default();
    let joined = boolean(&mut model, a, b, BoolOp::Unite, &mut dec)
        .unwrap()
        .body;
    let hollow = boolean(&mut model, joined, cavity, BoolOp::Subtract, &mut dec)
        .unwrap()
        .body;
    let answer = boolean(&mut model, hollow, clip, BoolOp::Intersect, &mut dec).unwrap();
    assert!(answer.marks.is_empty());
    topology_rules(&model, answer.body);
    let counts = model.counts(answer.body);
    // In the unrotated frame: a clipped box with one spherical cavity, plus
    // half a cylinder touching the box along its innermost generator.
    assert_eq!((counts.lumps, counts.shells), (2, 3));
    for (point, want) in [
        (["-19/10", "-5/2", "-5/2"], Location::Inside),
        (["1/2", "1/3", "1/5"], Location::Outside),
        (["5/2", "1/3", "1/5"], Location::On),
        (["7/2", "0", "0"], Location::Inside),
        (["3", "0", "0"], Location::On),
        (["4", "0", "0"], Location::On),
        (["9/2", "0", "0"], Location::Outside),
    ] {
        let [x, y, z] = point;
        query(
            &mut model,
            answer.body,
            &global
                .apply_point(&p(x, y, z).to_real())
                .as_rational()
                .expect("rational Boolean fixture"),
            want,
            &mut dec,
        );
    }
}

blind_test!(cases_41_to_44_chains, {
    for number in 41..=44 {
        run_chain(number, &r(1), &Pose::identity());
    }
});

fn near_points(shape: &Primitive) -> Vec<RPoint3> {
    let mut out = Vec::new();
    for sign in [-1, 1] {
        let eps = r(sign) * q("1e-30");
        let local = match &shape.shape {
            Shape::Block(s) => RPoint3::new(&s[0] + &eps, &s[1] / &r(2), &s[2] / &r(2)),
            Shape::Cylinder(radius, height) => RPoint3::new(
                (radius + &eps) * q("3/5"),
                (radius + &eps) * q("4/5"),
                height / &r(2),
            ),
            Shape::Ball(radius) => {
                RPoint3::new((radius + &eps) * q("3/5"), (radius + &eps) * q("4/5"), r(0))
            }
        };
        out.push(
            shape
                .pose
                .apply_point(&local.to_real())
                .as_rational()
                .expect("rational Boolean fixture"),
        );
    }
    out
}

blind_test!(case_54_primitive_exact_points_edges_and_corners, {
    let primitives = [
        block_at([2, 3, 5], [0; 3]),
        cyl_at(2, 5, [0; 3], false),
        ball_at(5, [0; 3]),
    ];
    let queries = [
        vec![
            (p("1", "1", "1"), Location::Inside),
            (p("2", "1", "1"), Location::On),
            (p("3", "1", "1"), Location::Outside),
            (p("2", "3", "5"), Location::On),
        ],
        vec![
            (p("0", "0", "1"), Location::Inside),
            (p("2", "0", "1"), Location::On),
            (p("6/5", "8/5", "5"), Location::On),
            (p("2", "2", "1"), Location::Outside),
        ],
        vec![
            (p("3", "4", "0"), Location::On),
            (p("3", "4", "1e-30"), Location::Outside),
            (
                RPoint3::new(r(3), r(4) - q("1e-30"), r(0)),
                Location::Inside,
            ),
        ],
    ];
    for (primitive, queries) in primitives.into_iter().zip(queries) {
        let mut m = Model::new();
        let b = primitive.build(&mut m);
        let before = topology_snapshot(&m, &[b]);
        for (point, want) in queries {
            query(&mut m, b, &point, want, &mut Decider::default());
        }
        assert_eq!(before, topology_snapshot(&m, &[b]));
    }
});

blind_test!(
    cases_55_56_generated_material_identities_and_commutativity,
    {
        let mut questions = 0;
        for &id in IDS {
            let c = case(id);
            for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
                let expected = expectations(id)
                    .into_iter()
                    .find(|(o, reverse, _)| *o == op && !*reverse)
                    .map(|(_, _, e)| e);
                if matches!(expected, Some(Expected::Refused(..))) {
                    run_case(&c, op, false, expected.as_ref(), 0);
                    continue;
                }
                let mut m = Model::new();
                let a = c.a.build(&mut m);
                let b = c.b.build(&mut m);
                let mut dec = Decider::default();
                let answer = boolean(&mut m, a, b, op, &mut dec).unwrap();
                assert!(answer.marks.is_empty());
                topology_rules(&m, answer.body);
                if let Some(Expected::Counts(want)) = &expected {
                    check_counts(&m, answer.body, want);
                }
                let mut rng = Rng(0x550056);
                let mut n = 0;
                let near: Vec<_> = near_points(&c.a)
                    .into_iter()
                    .chain(near_points(&c.b))
                    .collect();
                for point in near {
                    if let Some(want) = set_rule(op, c.a.locate(&point), c.b.locate(&point)) {
                        query(&mut m, answer.body, &point, want, &mut dec);
                        questions += 1;
                    }
                }
                while n < 256 {
                    let point = rng.point();
                    if let Some(want) = set_rule(op, c.a.locate(&point), c.b.locate(&point)) {
                        query(&mut m, answer.body, &point, want, &mut dec);
                        questions += 1;
                        n += 1;
                    }
                }
                if op != BoolOp::Subtract {
                    let swapped = boolean(&mut m, b, a, op, &mut Decider::default()).unwrap();
                    assert!(swapped.marks.is_empty());
                    assert_eq!(counts(&m, swapped.body), counts(&m, answer.body));
                    topology_rules(&m, swapped.body);
                }
            }
        }
        assert!(
            questions >= 30_000,
            "only {questions} maths-answered queries executed"
        );
        eprintln!("cases 55/56: {questions} exact result-containment questions");
    }
);

blind_test!(case_55_random_rational_primitive_poses_near_and_far, {
    let mut rng = Rng(0x55a11);
    for _ in 0..100 {
        let pose = rng.pose();
        for primitive in [
            block_at([2, 3, 5], [0; 3]),
            cyl_at(2, 5, [0; 3], false),
            ball_at(5, [0; 3]),
        ] {
            let posed = primitive.transformed(&r(1), &pose);
            let mut m = Model::new();
            let b = posed.build(&mut m);
            let mut dec = Decider::default();
            for point in near_points(&posed) {
                query(&mut m, b, &point, posed.locate(&point), &mut dec);
            }
            for _ in 0..64 {
                let local = rng.point();
                let point = pose
                    .apply_point(&local.to_real())
                    .as_rational()
                    .expect("rational Boolean fixture");
                let want = primitive.locate(&local);
                query(&mut m, b, &point, want, &mut dec);
            }
        }
    }
});

fn placement_poses() -> [Pose; 20] {
    let mut rng = Rng(0x570020);
    std::array::from_fn(|_| rng.pose())
}

fn placed_pair(
    c: &Case,
    op: BoolOp,
    reverse: bool,
    expected: &Expected,
    scale: &Rational,
    pose: &Pose,
) {
    let pa = c.a.transformed(scale, pose);
    let pb = c.b.transformed(scale, pose);
    let mut m = Model::new();
    let a = pa.build(&mut m);
    let b = pb.build(&mut m);
    let (a, b) = if reverse { (b, a) } else { (a, b) };
    let mut dec = Decider::default();
    let answer = boolean(&mut m, a, b, op, &mut dec);
    match expected {
        Expected::Refused(code, stage) => {
            let e = answer.unwrap_err();
            assert_eq!(e.refusal.code, *code);
            assert_eq!(e.stage, *stage);
        }
        Expected::Counts(want) => {
            let answer = answer.unwrap_or_else(|error| {
                panic!(
                    "case {} {op:?}, reverse={reverse}, scale={scale}, pose={pose:?}: {error:?}; statistics: {:?}",
                    c.id,
                    dec.stats()
                )
            });
            assert!(answer.marks.is_empty());
            check_counts(&m, answer.body, want);
            if !reverse {
                for (point, want) in witnesses(c.id, op) {
                    query(
                        &mut m,
                        answer.body,
                        &moved_point(&point, scale, pose),
                        want,
                        &mut dec,
                    );
                }
            }
            let point = moved_point(&p("1/7", "2/7", "3/7"), scale, pose);
            let (a, b) = if reverse { (&pb, &pa) } else { (&pa, &pb) };
            if let Some(want) = set_rule(op, a.locate(&point), b.locate(&point)) {
                query(&mut m, answer.body, &point, want, &mut dec);
            }
        }
    }
}

fn placed_family(scale: &Rational, pose: &Pose) -> usize {
    let mut operations = 0;
    for &id in IDS {
        let c = case(id);
        for (op, reverse, expected) in expectations(id) {
            placed_pair(&c, op, reverse, &expected, scale, pose);
            operations += 1;
        }
    }
    for number in 41..=44 {
        let (_, _, steps, _) = run_chain(number, scale, pose);
        operations += steps.len();
    }
    operations
}

blind_test!(case_57_default_scale_and_placement_variants, {
    let poses = placement_poses();
    let identity = Pose::identity();
    let small = q("1e-30");
    let large = q("1e30");
    let mut rng = Rng(0x570057);
    let mut operations = 0;
    for &id in IDS {
        let c = case(id);
        for (op, reverse, expected) in expectations(id) {
            for scale in [&small, &large] {
                placed_pair(&c, op, reverse, &expected, scale, &identity);
            }
            let k = (rng.next() % 59) as i32 - 29;
            let scale = q(&format!("1e{k}"));
            let pose = &poses[operations % poses.len()];
            placed_pair(&c, op, reverse, &expected, &scale, pose);
            operations += 1;
        }
    }
    for number in 41..=44 {
        for scale in [&small, &large] {
            run_chain(number, scale, &identity);
        }
        let k = (rng.next() % 59) as i32 - 29;
        let scale = q(&format!("1e{k}"));
        let (_, _, steps, _) = run_chain(number, &scale, &poses[operations % poses.len()]);
        operations += steps.len();
    }
    eprintln!("case 57 default: {} booleans", operations * 3);
});

blind_test!(
    #[ignore = "long-run qualification: 18,392 booleans; run on demand"]
    case_57_long_run_scale_and_placement_qualification,
    {
        let start = Instant::now();
        let poses = placement_poses();
        let identity = Pose::identity();
        let mut operations = 0;
        for k in -30..=30 {
            operations += placed_family(&q(&format!("1e{k}")), &identity);
        }
        for k in [-30, 0, 30] {
            let scale = q(&format!("1e{k}"));
            for pose in &poses {
                operations += placed_family(&scale, pose);
            }
        }
        eprintln!(
            "case 57 long-run: {operations} booleans in {:?}",
            start.elapsed()
        );
    }
);

blind_test!(case_58_two_fresh_runs_and_eight_threads, {
    for &id in IDS {
        for (op, reverse, expected) in expectations(id) {
            let first = run_case(&case(id), op, reverse, Some(&expected), 0);
            assert_eq!(
                first,
                run_case(&case(id), op, reverse, Some(&expected), 0),
                "{id} {op:?} reverse={reverse}"
            );
            let results = std::thread::scope(|scope| {
                let workers: Vec<_> = (0..8)
                    .map(|_| {
                        let expected = &expected;
                        scope.spawn(move || run_case(&case(id), op, reverse, Some(expected), 0))
                    })
                    .collect();
                workers
                    .into_iter()
                    .map(|w| w.join().unwrap())
                    .collect::<Vec<_>>()
            });
            for result in results {
                assert_eq!(first, result, "{id} {op:?}");
            }
        }
    }
    for number in 41..=44 {
        let first = run_chain(number, &r(1), &Pose::identity());
        assert_eq!(first, run_chain(number, &r(1), &Pose::identity()));
        std::thread::scope(|s| {
            let workers: Vec<_> = (0..8)
                .map(|_| s.spawn(move || run_chain(number, &r(1), &Pose::identity())))
                .collect();
            for w in workers {
                assert_eq!(first, w.join().unwrap());
            }
        });
    }
});

fn check_budget_refusal(e: &BooleanRefusal) {
    match e.refusal.code {
        "boolean.imprint-ambiguous" => {
            assert!(matches!(e.stage, Stage::Intersect | Stage::Imprint))
        }
        "boolean.merge-ambiguous" | "boolean.self-touching-shell" => {
            assert_eq!(e.stage, Stage::Sew)
        }
        "boolean.classify-ambiguous" | "boolean.no-rational-point" => {}
        "numbers.budget-exhausted"
        | "numbers.precision-limit"
        | "numbers.root-not-certified"
        | "numbers.sqrt-not-separated" => {}
        code => panic!("unexpected refusal for a supported, valid pair: {code}"),
    }
}

blind_test!(
    case_59_low_budget_answers_are_sound_and_deadlines_are_honest,
    {
        for &id in IDS {
            let c = case(id);
            for (op, reverse, expected) in expectations(id) {
                for expired in [false, true] {
                    let mut m = Model::new();
                    let a = c.a.build(&mut m);
                    let b = c.b.build(&mut m);
                    let (a, b) = if reverse { (b, a) } else { (a, b) };
                    let before = topology_snapshot(&m, &[a, b]);
                    let mut dec = if expired {
                        Decider::default()
                    } else {
                        Decider::new(BudgetPolicy {
                            max_precision_bits: 64,
                            max_work: 10_000,
                        })
                    };
                    if expired {
                        dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
                    }
                    match boolean(&mut m, a, b, op, &mut dec) {
                        Err(e) => {
                            if expired {
                                if e.refusal.code != "numbers.deadline" {
                                    let Expected::Refused(code, stage) = &expected else {
                                        panic!("unexpected deadline refusal: {e:?}")
                                    };
                                    assert_eq!(e.refusal.code, *code);
                                    assert_eq!(e.stage, *stage);
                                }
                            } else {
                                check_budget_refusal(&e);
                            }
                        }
                        Ok(answer) => {
                            topology_rules(&m, answer.body);
                            assert_eq!(answer.marks, m.body(answer.body).marks);
                            for mark in answer.marks.iter() {
                                assert_eq!(mark.code, "numbers.budget-exhausted");
                                assert!(mark.bound >= r(0));
                            }
                            if expired || answer.marks.is_empty() {
                                let Expected::Counts(want) = &expected else {
                                    panic!("answered a mandatory refusal")
                                };
                                check_counts(&m, answer.body, want);
                                let (pa, pb) = if reverse { (&c.b, &c.a) } else { (&c.a, &c.b) };
                                let mut rng = Rng(0x590001);
                                for _ in 0..32 {
                                    let p = rng.point();
                                    if let Some(want) = set_rule(op, pa.locate(&p), pb.locate(&p)) {
                                        query(
                                            &mut m,
                                            answer.body,
                                            &p,
                                            want,
                                            &mut Decider::default(),
                                        );
                                    }
                                }
                            }
                        }
                    }
                    assert_eq!(before, topology_snapshot(&m, &[a, b]));
                }
            }
        }
    }
);

fn within_radial_distance(s: &Rational, radius: &Rational, height_sq: &Rational) -> bool {
    if height_sq < &r(0) {
        return false;
    }
    let left = &(s + &sq(radius)) - height_sq;
    left <= r(0) || sq(&left) <= &(r(4) * sq(radius)) * s
}
fn distance_bound_contains_boundary(shape: &Primitive, point: &RPoint3, bound: &Rational) -> bool {
    let p = shape
        .pose
        .inverse()
        .apply_point(&point.to_real())
        .as_rational()
        .expect("rational Boolean fixture");
    let bound_sq = sq(bound);
    match &shape.shape {
        Shape::Ball(radius) => {
            within_radial_distance(&(sq(&p.x) + sq(&p.y) + sq(&p.z)), radius, &bound_sq)
        }
        Shape::Cylinder(radius, h) => {
            let rho_sq = sq(&p.x) + sq(&p.y);
            let axial_gap = if p.z < r(0) {
                -&p.z
            } else if &p.z > h {
                &p.z - h
            } else {
                r(0)
            };
            let side = within_radial_distance(&rho_sq, radius, &(&bound_sq - &sq(&axial_gap)));
            let cap = |z: Rational| {
                let remaining = &bound_sq - &sq(&(&p.z - &z));
                if rho_sq <= sq(radius) {
                    remaining >= r(0)
                } else {
                    within_radial_distance(&rho_sq, radius, &remaining)
                }
            };
            side || cap(r(0)) || cap(h.clone())
        }
        Shape::Block(size) => {
            let coordinates = [&p.x, &p.y, &p.z];
            (0..3).any(|axis| {
                [r(0), size[axis].clone()].iter().any(|plane| {
                    let mut distance = sq(&(coordinates[axis] - plane));
                    for i in 0..3 {
                        if i == axis {
                            continue;
                        }
                        let gap = if coordinates[i] < &r(0) {
                            -coordinates[i]
                        } else if coordinates[i] > &size[i] {
                            coordinates[i] - &size[i]
                        } else {
                            r(0)
                        };
                        distance = distance + sq(&gap);
                    }
                    distance <= bound_sq
                })
            })
        }
    }
}
blind_test!(classify_point_not_exact_has_a_true_distance_bound, {
    let mut rng = Rng(0xb0a0d);
    for shape in [
        block_at([2, 3, 5], [0; 3]),
        cyl_at(2, 5, [0; 3], false),
        ball_at(5, [0; 3]),
    ] {
        let mut m = Model::new();
        let b = shape.build(&mut m);
        for work in [0, 1, 10_000] {
            let mut points = near_points(&shape);
            points.extend((0..100).map(|_| rng.point()));
            for point in points {
                let mut dec = Decider::new(BudgetPolicy {
                    max_precision_bits: 64,
                    max_work: work,
                });
                let want = shape.locate(&point);
                match classify_point(&mut m, b, &point, &mut dec) {
                    Outcome::Proven(Containment::Inside) => assert_eq!(want, Location::Inside),
                    Outcome::Proven(Containment::Outside) => assert_eq!(want, Location::Outside),
                    Outcome::ProvenCoincident(_) => assert_eq!(want, Location::On),
                    Outcome::NotExact(n) => {
                        assert_eq!(n.code, "numbers.budget-exhausted");
                        assert!(n.bound >= r(0));
                        assert!(
                            distance_bound_contains_boundary(&shape, &point, &n.bound),
                            "false Distance bound {n:?} at {point:?}"
                        );
                    }
                    Outcome::Refused(e) => assert!(
                        matches!(
                            e.code,
                            "boolean.classify-ambiguous"
                                | "numbers.budget-exhausted"
                                | "numbers.root-not-certified"
                        ),
                        "{e:?}"
                    ),
                }
            }
        }
    }
});
