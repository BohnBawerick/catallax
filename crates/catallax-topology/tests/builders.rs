//! Blind primitive, refusal, scale, placement and determinism contracts.
#[path = "common/assertions.rs"]
mod assertions;
mod common;
use assertions::*;
use common::*;

blind_test!(cases_01_to_05_exact_builders, {
    let mut m = Model::new();
    for shape in primitives() {
        let b = shape.build(&mut m, &Pose::identity()).unwrap();
        assert_builder(&m, b, &shape, &Pose::identity());
    }
});

blind_test!(case_06_empty_adds_only_a_body, {
    let mut m = Model::new();
    primitives()[0].build(&mut m, &Pose::identity()).unwrap();
    let before = lengths(&m);
    let b = empty(&mut m);
    let mut want = before;
    want[10] += 1;
    assert_eq!(lengths(&m), want);
    assert_eq!(b.index() as usize, before[10]);
    assert_eq!(m.body(b).kind, BodyKind::Empty);
    assert!(m.body(b).lumps.is_empty());
    assert!(m.body(b).marks.is_empty());
    assert_eq!(m.counts(b), EntityCounts::default());
    assert_eq!(walk(&m, b), (vec![], vec![], vec![]));
    assert!(m.faces_of_body(b).is_empty());
    assert!(m.edges_of_body(b).is_empty());
    assert!(m.vertices_of_body(b).is_empty());
    sound(&m, b);
});

blind_test!(case_07_exact_arena_growth_and_no_reuse, {
    let mut m = Model::new();
    // Surface, curve, point, vertex, edge, coedge, loop, face, shell, lump, body.
    let rows = [
        [6, 12, 8, 8, 12, 24, 6, 6, 1, 1, 1],
        [3, 2, 0, 0, 2, 4, 4, 3, 1, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        [3, 2, 0, 0, 2, 4, 4, 3, 1, 1, 1],
        [3, 2, 0, 0, 2, 4, 4, 3, 1, 1, 1],
        [2, 1, 1, 1, 1, 2, 3, 2, 1, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
    ];
    assert_eq!(lengths(&m), [0; 12]);
    for _ in 0..2 {
        for (shape, row) in primitives().iter().zip(rows) {
            let before = lengths(&m);
            let b = shape.build(&mut m, &Pose::identity()).unwrap();
            let after = lengths(&m);
            for i in 0..11 {
                assert_eq!(after[i] - before[i], row[i], "{} arena {i}", shape.name());
            }
            let (fs, es, vs) = walk(&m, b);
            for f in fs {
                let f = m.face(f);
                assert!(f.surface.index() as usize >= before[0]);
            }
            for e in es {
                assert!(e.index() as usize >= before[4]);
                assert!(m.edge(e).curve.index() as usize >= before[1]);
            }
            for v in vs {
                assert!(v.index() as usize >= before[3]);
                assert!(m.vertex(v).point.index() as usize >= before[2]);
            }
            assert_builder(&m, b, shape, &Pose::identity());
        }
    }
});

blind_test!(case_08_quarter_turn_has_named_corner_and_plane, {
    let pose = Pose::new(
        Rotation::from_quaternion(r(1), r(0), r(0), r(1)).unwrap(),
        v(10, 0, 0).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    let shape = Primitive::Block([r(1), r(1), r(1)]);
    let mut m = Model::new();
    let b = shape.build(&mut m, &pose).unwrap();
    assert_builder(&m, b, &shape, &pose);
    assert!(
        m.vertices_of_body(b)
            .iter()
            .any(|v| explicit(&m, *v) == p(10, 1, 0))
    );
    let f = m.faces_of_body(b)[0];
    assert!(m.geom().surface(m.face(f).surface).same_definition(&plane(
        p(10, 0, 0),
        v(0, -1, 0),
        v(-1, 0, 0)
    )));
});

blind_test!(cases_08_09_generated_rational_poses_and_orientations, {
    let mut rng = Rng(0x3b7e_2026_1019);
    // 1,100 * 7 = 7,700 builder-and-check questions, each with exact geometry,
    // outward implicit signs, block left-hand rule, topology and incidence checks.
    for _ in 0..1100 {
        let pose = rng.pose();
        let scale = q((rng.next() % 17 + 1) as i64, (rng.next() % 13 + 1) as i64);
        for shape in primitives() {
            let shape = shape.scaled(&scale);
            let mut m = Model::new();
            let b = shape.build(&mut m, &pose).unwrap();
            assert_builder(&m, b, &shape, &pose);
        }
    }
});

blind_test!(case_11_refusal_codes_precedence_and_atomicity, {
    let mut inputs = vec![];
    for axis in 0..3 {
        for bad in [0, -1] {
            let mut sizes = [r(2), r(3), r(5)];
            sizes[axis] = r(bad);
            inputs.push((Primitive::Block(sizes), "topology.non-positive-size"));
        }
    }
    for bad in [0, -1] {
        inputs.extend([
            (
                Primitive::Cylinder(r(bad), r(5)),
                "topology.non-positive-size",
            ),
            (
                Primitive::Cylinder(r(2), r(bad)),
                "topology.non-positive-size",
            ),
            (Primitive::Sphere(r(bad)), "topology.non-positive-size"),
            (
                Primitive::Cone(r(bad), r(1), r(2)),
                "topology.non-positive-size",
            ),
            (
                Primitive::Cone(r(2), r(1), r(bad)),
                "topology.non-positive-size",
            ),
            (Primitive::Torus(r(bad), r(1)), "topology.non-positive-size"),
            (Primitive::Torus(r(3), r(bad)), "topology.non-positive-size"),
        ]);
    }
    inputs.extend([
        (
            Primitive::Cone(r(2), r(-1), r(1)),
            "topology.non-positive-size",
        ),
        (
            Primitive::Cone(r(2), r(2), r(1)),
            "topology.cone-is-cylinder",
        ),
        (
            Primitive::Cone(r(-1), r(-1), r(1)),
            "topology.non-positive-size",
        ),
        (
            Primitive::Cone(r(2), r(2), r(0)),
            "topology.non-positive-size",
        ),
        (
            Primitive::Torus(r(1), r(1)),
            "topology.torus-self-intersecting",
        ),
        (
            Primitive::Torus(r(1), r(2)),
            "topology.torus-self-intersecting",
        ),
        (Primitive::Torus(r(-2), r(-1)), "topology.non-positive-size"),
    ]);
    for populated in [false, true] {
        let mut m = Model::new();
        let mut bodies = vec![];
        if populated {
            for shape in primitives() {
                bodies.push(shape.build(&mut m, &Pose::identity()).unwrap());
            }
        }
        for (shape, code) in &inputs {
            let before = snapshot(&m, &bodies);
            let sizes = lengths(&m);
            // A panic fails this test directly. No catch_unwind converts it to a refusal.
            assert_eq!(
                shape.build(&mut m, &Pose::identity()).unwrap_err().code,
                *code
            );
            assert_eq!(lengths(&m), sizes);
            assert_eq!(snapshot(&m, &bodies), before);
        }
    }
});

blind_test!(case_12_scale_family_preserves_exact_fields, {
    for k in -30..=30 {
        let factor = Rational::parse(&format!("1e{k}")).unwrap();
        for shape in primitives() {
            let shape = shape.scaled(&factor);
            let mut m = Model::new();
            let b = shape.build(&mut m, &Pose::identity()).unwrap();
            // The oracle reconstructs scaled positions/radii, squares each circle radius,
            // and divides scaled radius differences by scaled height for unchanged slopes.
            assert_builder(&m, b, &shape, &Pose::identity());
        }
    } // 61 * 7 = 427 questions.
});

blind_test!(case_33_determinism_across_models_and_eight_threads, {
    fn run() -> Vec<Snapshot> {
        let mut rng = Rng(0x00de_cafb_ad33);
        let mut results = vec![];
        for _ in 0..48 {
            let mut m = Model::new();
            let mut bodies = vec![];
            let pose = rng.pose();
            for shape in primitives() {
                let b = shape.build(&mut m, &pose).unwrap();
                sound(&m, b);
                bodies.push(b);
            }
            results.push(snapshot(&m, &bodies));
        }
        results
    }
    let baseline = run();
    assert_eq!(baseline, run());
    let workers: Vec<_> = (0..8).map(|_| std::thread::spawn(run)).collect();
    for worker in workers {
        assert_eq!(baseline, worker.join().unwrap());
    }
    // Ten runs * 48 sequences * seven shapes = 3,360 questions.
    // Generated total: 7,700 + 427 + 3,360 = 11,487, independent of test scheduling.
});
