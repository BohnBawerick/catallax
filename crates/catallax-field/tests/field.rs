//! Independent geometric answers, exact supports, Boolean trimming and reproducibility.

use catallax_field::{Batch, Direction, FaceState, Location, Plan, Query, Value, Workspace};
use catallax_numbers::{Interval, Rational, Sign};
use catallax_recipe::{
    BooleanKind, ElementId, Number, Operation, Origin, Pose, Recipe, Role, Step, StepId,
};
use catallax_shapes::{RPoint3, RVector3, Surface};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn point(p: [i64; 3]) -> RPoint3 {
    let [x, y, z] = p.map(r);
    RPoint3::new(x, y, z)
}
fn fraction(p: [&str; 3]) -> RPoint3 {
    let [x, y, z] = p.map(|s| Rational::parse(s).unwrap());
    RPoint3::new(x, y, z)
}
fn id(s: &str) -> StepId {
    StepId::new(s).unwrap()
}
fn block(size: i64, shift: [i64; 3]) -> Operation {
    let mut pose = Pose::identity();
    pose.translation.0 = shift.map(Number::from);
    Operation::Block {
        size: [size.into(), size.into(), size.into()],
        pose,
    }
}
fn cylinder(radius: i64, height: i64) -> Operation {
    Operation::Cylinder {
        radius: radius.into(),
        height: height.into(),
        pose: Pose::identity(),
    }
}
fn sphere(radius: i64) -> Operation {
    Operation::Sphere {
        radius: radius.into(),
        pose: Pose::identity(),
    }
}
fn cone() -> Operation {
    Operation::Cone {
        base_radius: 3.into(),
        top_radius: 0.into(),
        height: 4.into(),
        pose: Pose::identity(),
    }
}
fn torus() -> Operation {
    Operation::Torus {
        major: 4.into(),
        minor: 1.into(),
        pose: Pose::identity(),
    }
}
fn primitive(op: Operation) -> Plan {
    Plan::compile(&Recipe::new(
        "solid",
        vec![Step {
            id: id("solid"),
            op,
        }],
        id("solid"),
    ))
    .unwrap()
}
fn boolean(a: Operation, b: Operation, kind: BooleanKind) -> Plan {
    Plan::compile(&Recipe::new(
        "pair",
        vec![
            Step { id: id("a"), op: a },
            Step { id: id("b"), op: b },
            Step {
                id: id("part"),
                op: Operation::Boolean {
                    kind,
                    target: id("a"),
                    tool: id("b"),
                },
            },
        ],
        id("part"),
    ))
    .unwrap()
}
fn query(plan: &Plan, p: &RPoint3) -> Query {
    plan.proven(p, &mut Workspace::default())
}
fn contains(interval: Interval, exact: &Rational) {
    if interval.lo().is_finite() {
        assert!(
            Rational::from_f64_exact(interval.lo()).unwrap() <= *exact,
            "{interval:?} excludes {exact}"
        );
    }
    if interval.hi().is_finite() {
        assert!(
            Rational::from_f64_exact(interval.hi()).unwrap() >= *exact,
            "{interval:?} excludes {exact}"
        );
    }
}
fn distance(plan: &Plan, p: &RPoint3, location: Location, exact: &str) {
    let result = query(plan, p);
    assert_eq!(result.location, location, "at {p:?}");
    let exact = Rational::parse(exact).unwrap();
    let Value::TrueDistance(value) = result.value else {
        panic!("primitive distance mislabelled");
    };
    contains(value, &exact);
    contains(result.distance.millimetres, &exact.abs());
    let fast = plan
        .fast(
            [
                p.x.to_f64_nearest(),
                p.y.to_f64_nearest(),
                p.z.to_f64_nearest(),
            ],
            &mut Workspace::default(),
        )
        .unwrap();
    assert!((fast.value - exact.to_f64_nearest()).abs() < 1e-12);
    assert!(fast.true_distance_formula);
}
fn unit(query: Direction, expected: [i64; 3]) {
    let Direction::Unique(components) = query else {
        panic!("expected unique normal, got {query:?}");
    };
    for (component, value) in components.into_iter().zip(expected) {
        contains(component, &r(value));
    }
}

#[test]
fn every_primitive_inside_outside_boundary_and_true_distance() {
    for (op, cases) in [
        (
            block(2, [0, 0, 0]),
            vec![
                (point([1, 1, 1]), Location::Inside, "-1"),
                (point([5, 6, 2]), Location::Outside, "5"),
                (point([0, 1, 1]), Location::Boundary, "0"),
            ],
        ),
        (
            cylinder(2, 4),
            vec![
                (point([0, 0, 2]), Location::Inside, "-2"),
                (point([5, 0, 8]), Location::Outside, "5"),
                (point([2, 0, 2]), Location::Boundary, "0"),
            ],
        ),
        (
            sphere(2),
            vec![
                (point([0, 0, 0]), Location::Inside, "-2"),
                (point([0, 0, 5]), Location::Outside, "3"),
                (point([0, 2, 0]), Location::Boundary, "0"),
            ],
        ),
        (
            cone(),
            vec![
                (point([0, 0, 1]), Location::Inside, "-1"),
                (fraction(["11/2", "0", "5"]), Location::Outside, "5"),
                (fraction(["3/2", "0", "2"]), Location::Boundary, "0"),
            ],
        ),
        (
            torus(),
            vec![
                (point([4, 0, 0]), Location::Inside, "-1"),
                (point([4, 0, 3]), Location::Outside, "2"),
                (point([0, 0, 0]), Location::Outside, "3"),
                (point([5, 0, 0]), Location::Boundary, "0"),
            ],
        ),
    ] {
        let plan = primitive(op);
        for (p, loc, exact) in cases {
            distance(&plan, &p, loc, exact);
        }
    }
}

#[test]
fn irrational_answers_are_enclosed_not_rounded_into_exact_distances() {
    let q = query(&primitive(block(2, [0, 0, 0])), &point([3, 3, 1]));
    let Value::TrueDistance(d) = q.value else {
        unreachable!()
    };
    let lo = Rational::from_f64_exact(d.lo()).unwrap();
    let hi = Rational::from_f64_exact(d.hi()).unwrap();
    assert!(&lo * &lo <= r(2));
    assert!(&hi * &hi >= r(2));
    assert!(d.width() > 0.0 && d.width() < 1e-6);
    let q = query(&primitive(sphere(1)), &point([1, 1, 1]));
    let d = q.value.enclosure();
    let lo = Rational::from_f64_exact(d.lo()).unwrap() + r(1);
    let hi = Rational::from_f64_exact(d.hi()).unwrap() + r(1);
    assert!(&lo * &lo <= r(3));
    assert!(&hi * &hi >= r(3));
}

#[test]
fn primitive_boundaries_agree_with_independent_exact_shape_definitions() {
    let axis = RVector3::new(r(0), r(0), r(1));
    let reference = RVector3::new(r(1), r(0), r(0));
    let cases = [
        (
            block(2, [0, 0, 0]),
            Surface::Plane(
                catallax_shapes::Plane::new(
                    point([2, 0, 0]).into(),
                    reference.clone().into(),
                    axis.clone().into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ),
            point([2, 1, 1]),
            point([1, 1, 1]),
            point([3, 1, 1]),
        ),
        (
            cylinder(2, 4),
            Surface::Cylinder(
                catallax_shapes::Cylinder::new(
                    point([0, 0, 0]).into(),
                    axis.clone().into(),
                    r(2).into(),
                    reference.clone().into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ),
            point([2, 0, 2]),
            point([1, 0, 2]),
            point([3, 0, 2]),
        ),
        (
            sphere(2),
            Surface::Sphere(
                catallax_shapes::Sphere::new(
                    point([0, 0, 0]).into(),
                    r(2).into(),
                    axis.clone().into(),
                    reference.clone().into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ),
            point([2, 0, 0]),
            point([1, 0, 0]),
            point([3, 0, 0]),
        ),
        (
            cone(),
            Surface::Cone(
                catallax_shapes::Cone::new(
                    point([0, 0, 4]).into(),
                    axis.scaled(&r(-1)).into(),
                    Rational::parse("3/4").unwrap().into(),
                    reference.clone().into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ),
            fraction(["3/2", "0", "2"]),
            point([1, 0, 2]),
            point([2, 0, 2]),
        ),
        (
            torus(),
            Surface::Torus(
                catallax_shapes::Torus::new(
                    point([0, 0, 0]).into(),
                    axis.into(),
                    r(4).into(),
                    r(1).into(),
                    reference.into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ),
            point([5, 0, 0]),
            fraction(["9/2", "0", "0"]),
            point([6, 0, 0]),
        ),
    ];
    for (op, surface, boundary, inside, outside) in cases {
        let plan = primitive(op);
        for (p, sign, loc) in [
            (boundary, Sign::Zero, Location::Boundary),
            (inside, Sign::Negative, Location::Inside),
            (outside, Sign::Positive, Location::Outside),
        ] {
            assert_eq!(
                surface
                    .implicit()
                    .rational()
                    .unwrap()
                    .eval(&[p.x.clone(), p.y.clone(), p.z.clone()])
                    .sign(),
                sign
            );
            assert_eq!(query(&plan, &p).location, loc);
        }
    }
}

#[test]
fn exact_rigid_placement_preserves_metric_and_rotates_surface_direction() {
    let mut pose = Pose::identity();
    pose.rotation = [1.into(), 0.into(), 1.into(), 0.into()];
    pose.translation.0 = [10.into(), 20.into(), 30.into()];
    let plan = primitive(Operation::Cylinder {
        radius: 2.into(),
        height: 4.into(),
        pose,
    });
    distance(&plan, &point([12, 20, 30]), Location::Inside, "-2");
    distance(&plan, &point([12, 20, 32]), Location::Boundary, "0");
    unit(query(&plan, &point([12, 20, 32])).direction, [0, 0, 1]);
    distance(&plan, &point([18, 20, 35]), Location::Outside, "5");
}

#[test]
fn sharp_edges_corners_and_apex_never_choose_an_arbitrary_normal() {
    let plan = primitive(block(2, [0, 0, 0]));
    for p in [[0, 0, 1], [0, 0, 0]] {
        assert_eq!(query(&plan, &point(p)).direction, Direction::Ambiguous);
    }
    unit(query(&plan, &point([0, 1, 1])).direction, [-1, 0, 0]);
    assert_eq!(
        query(&primitive(cylinder(2, 4)), &point([2, 0, 0])).direction,
        Direction::Ambiguous
    );
    assert_eq!(
        query(&primitive(cone()), &point([0, 0, 4])).direction,
        Direction::Ambiguous
    );
    assert_eq!(
        query(&primitive(sphere(2)), &point([0, 0, 0])).direction,
        Direction::Undefined
    );
}

#[test]
fn all_booleans_have_exact_signs_and_bounds_containing_the_known_box_answer() {
    for (kind, lo, hi) in [
        (BooleanKind::Union, 0, 6),
        (BooleanKind::Subtract, 0, 2),
        (BooleanKind::Intersect, 2, 4),
    ] {
        let plan = boolean(block(4, [0, 0, 0]), block(4, [2, 0, 0]), kind);
        for x in -1..=7 {
            let result = query(&plan, &point([x, 1, 1]));
            let loc = if x < lo || x > hi {
                Location::Outside
            } else if x == lo || x == hi {
                Location::Boundary
            } else {
                Location::Inside
            };
            assert_eq!(result.location, loc, "{kind:?}, x={x}");
            let exact = if x < lo {
                lo - x
            } else if x > hi {
                x - hi
            } else {
                (x - lo).min(hi - x).min(1)
            };
            contains(result.distance.millimetres, &r(exact));
            assert!(matches!(result.value, Value::Residual(_)));
            assert!(
                !plan
                    .fast([x as f64, 1.0, 1.0], &mut Workspace::default())
                    .unwrap()
                    .true_distance_formula
            );
        }
    }
}

#[test]
fn boolean_residual_is_not_true_distance_at_a_reentrant_corner() {
    let plan = boolean(
        block(4, [-4, -2, -2]),
        block(4, [-2, -4, -2]),
        BooleanKind::Union,
    );
    let result = query(&plan, &point([-1, -1, 0]));
    assert_eq!(result.location, Location::Inside);
    let Value::Residual(residual) = result.value else {
        panic!("union cannot claim true distance");
    };
    contains(residual, &r(-1));
    let bound = result.distance.millimetres;
    assert!(bound.lo() <= 1.0 && bound.hi() == f64::INFINITY);
    // The nearest exterior point is (0,0,0), at sqrt(2), not |residual|=1.
    assert!(
        Rational::from_f64_exact(bound.lo())
            .unwrap()
            .pow(2)
            .unwrap()
            <= r(2)
    );
}

#[test]
fn recipe_names_survive_subtraction_and_the_tool_normal_reverses() {
    let mut tool = cylinder(1, 4);
    let Operation::Cylinder { pose, .. } = &mut tool else {
        unreachable!()
    };
    pose.translation.0 = [2.into(), 2.into(), 0.into()];
    let plan = boolean(block(4, [0, 0, 0]), tool, BooleanKind::Subtract);
    let wall = Origin {
        step: id("b"),
        role: Role::Side {
            segment: ElementId::new("wall").unwrap(),
            path: None,
        },
    };
    let wall = plan
        .select(&wall)
        .expect("tool wall name must remain selectable");
    let hit = plan
        .face_query(wall, &point([3, 2, 2]), &mut Workspace::default())
        .unwrap();
    assert_eq!(hit.state, FaceState::OnFace);
    unit(hit.direction, [-1, 0, 0]);
    let hit = plan
        .face_query(wall, &point([3, 2, 5]), &mut Workspace::default())
        .unwrap();
    assert_eq!(hit.state, FaceState::NotOnFace);
    let face = Origin {
        step: id("a"),
        role: Role::Side {
            segment: ElementId::new("x_max").unwrap(),
            path: None,
        },
    };
    let face = plan.select(&face).unwrap();
    assert_eq!(
        plan.face_query(face, &point([4, 1, 1]), &mut Workspace::default())
            .unwrap()
            .state,
        FaceState::OnFace
    );
}

#[test]
fn optimisation_keeps_all_names_but_never_resurrects_removed_faces() {
    let plan = boolean(
        block(2, [0, 0, 0]),
        block(2, [0, 0, 0]),
        BooleanKind::Subtract,
    );
    let origin = Origin {
        step: id("a"),
        role: Role::Start,
    };
    let face = plan.select(&origin).unwrap();
    let result = query(&plan, &point([0, 0, 0]));
    assert_eq!(result.location, Location::Outside);
    assert_eq!(
        result.distance.millimetres.lo().to_bits(),
        f64::INFINITY.to_bits()
    );
    assert_eq!(
        plan.face_query(face, &point([1, 1, 0]), &mut Workspace::default())
            .unwrap()
            .state,
        FaceState::NotOnFace
    );
    let union = boolean(block(2, [0, 0, 0]), block(2, [0, 0, 0]), BooleanKind::Union);
    for step in ["a", "b"] {
        let name = Origin {
            step: id(step),
            role: Role::Start,
        };
        let face = union.select(&name).unwrap();
        assert_eq!(
            union
                .face_query(face, &point([1, 1, 0]), &mut Workspace::default())
                .unwrap()
                .state,
            FaceState::OnFace
        );
    }
}

#[test]
fn tangent_regularisation_is_unresolved_rather_than_a_false_boundary() {
    let plan = boolean(
        block(2, [0, 0, 0]),
        block(2, [2, 0, 0]),
        BooleanKind::Intersect,
    );
    let result = query(&plan, &point([2, 1, 1]));
    assert_eq!(result.location, Location::Indeterminate);
    assert_eq!(result.direction, Direction::Undefined);
    assert!(matches!(result.value, Value::Residual(_)));
    contains(result.distance.millimetres, &r(0));
}

#[test]
fn nonfinite_fast_inputs_and_impossible_geometry_are_explicit_errors() {
    let plan = primitive(sphere(2));
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            plan.fast([x, 0.0, 0.0], &mut Workspace::default())
                .unwrap_err()
                .code,
            "field.non-finite-point"
        );
    }
    for op in [
        block(0, [0, 0, 0]),
        cylinder(-1, 2),
        sphere(0),
        Operation::Torus {
            major: 1.into(),
            minor: 1.into(),
            pose: Pose::identity(),
        },
    ] {
        let error = Plan::compile(&Recipe::new(
            "invalid",
            vec![Step { id: id("a"), op }],
            id("a"),
        ))
        .unwrap_err();
        assert_eq!(error.code, "field.impossible-primitive");
    }
    assert_eq!(
        plan.fast([f64::MAX, 0.0, 0.0], &mut Workspace::default())
            .unwrap_err()
            .code,
        "field.fast-range"
    );
    let zero = plan
        .fast([2.0, -0.0, 0.0], &mut Workspace::default())
        .unwrap();
    assert_eq!(zero.value.to_bits(), 0);
}

#[test]
fn exact_predicates_retain_sub_nanometre_features_at_large_coordinates() {
    let mut pose = Pose::identity();
    pose.translation.0[0] = Number::parse("10000000000000000").unwrap();
    let plan = primitive(Operation::Sphere {
        radius: Number::parse("0.000000000001").unwrap(),
        pose,
    });
    distance_proven_only(
        &plan,
        &fraction(["10000000000000000.000000000001", "0", "0"]),
        Location::Boundary,
        "0",
    );
    distance_proven_only(
        &plan,
        &fraction(["10000000000000000.000000000002", "0", "0"]),
        Location::Outside,
        "0.000000000001",
    );
}
fn distance_proven_only(plan: &Plan, p: &RPoint3, location: Location, exact: &str) {
    let q = query(plan, p);
    assert_eq!(q.location, location);
    contains(q.value.enclosure(), &Rational::parse(exact).unwrap());
}

fn bits(q: Query) -> Vec<u64> {
    let loc = match q.location {
        Location::Inside => 0,
        Location::Outside => 1,
        Location::Boundary => 2,
        Location::Indeterminate => 3,
        _ => panic!("new classification needs a fingerprint"),
    };
    let value = q.value.enclosure();
    let mut bits = vec![
        loc,
        u64::from(matches!(q.value, Value::TrueDistance(_))),
        value.lo().to_bits(),
        value.hi().to_bits(),
        q.distance.millimetres.lo().to_bits(),
        q.distance.millimetres.hi().to_bits(),
    ];
    match q.direction {
        Direction::Undefined => bits.push(0),
        Direction::Ambiguous => bits.push(1),
        Direction::Unique(n) => {
            bits.push(2);
            for i in n {
                bits.extend([i.lo().to_bits(), i.hi().to_bits()]);
            }
        }
        _ => panic!("new direction needs a fingerprint"),
    }
    bits
}

#[test]
fn batches_preserve_order_and_identical_bits_with_one_two_and_four_threads() {
    let plans = [
        primitive(block(2, [0, 0, 0])),
        primitive(cylinder(2, 4)),
        primitive(sphere(2)),
        primitive(cone()),
        primitive(torus()),
        boolean(block(4, [0, 0, 0]), sphere(2), BooleanKind::Union),
        boolean(block(4, [0, 0, 0]), sphere(2), BooleanKind::Subtract),
        boolean(block(4, [0, 0, 0]), sphere(2), BooleanKind::Intersect),
    ];
    let points: Vec<_> = (-2..=5)
        .flat_map(|x| [point([x, 0, 0]), point([x, 1, 1]), point([x, 2, 3])])
        .collect();
    let floats: Vec<_> = points
        .iter()
        .map(|p| {
            [
                p.x.to_f64_nearest(),
                p.y.to_f64_nearest(),
                p.z.to_f64_nearest(),
            ]
        })
        .collect();
    for plan in plans {
        let expected: Vec<_> = plan
            .proven_batch(&points, &mut Batch::default())
            .iter()
            .copied()
            .map(bits)
            .collect();
        let fast: Vec<_> = plan
            .fast_batch(&floats, &mut Batch::default())
            .unwrap()
            .iter()
            .map(|q| q.value.to_bits())
            .collect();
        for count in [1, 2, 4] {
            std::thread::scope(|scope| {
                let handles: Vec<_> = points
                    .chunks(points.len().div_ceil(count))
                    .zip(floats.chunks(floats.len().div_ceil(count)))
                    .map(|(points, floats)| {
                        let plan = &plan;
                        scope.spawn(move || {
                            let mut batch = Batch::default();
                            let proven = plan
                                .proven_batch(points, &mut batch)
                                .iter()
                                .copied()
                                .map(bits)
                                .collect::<Vec<_>>();
                            let fast = plan
                                .fast_batch(floats, &mut batch)
                                .unwrap()
                                .iter()
                                .map(|q| q.value.to_bits())
                                .collect::<Vec<_>>();
                            // Reuse both result and register storage with a shorter batch.
                            assert_eq!(
                                bits(plan.proven_batch(&points[..1], &mut batch)[0]),
                                proven[0]
                            );
                            (proven, fast)
                        })
                    })
                    .collect();
                let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
                assert_eq!(
                    results.iter().flat_map(|r| r.0.clone()).collect::<Vec<_>>(),
                    expected
                );
                assert_eq!(
                    results.iter().flat_map(|r| r.1.clone()).collect::<Vec<_>>(),
                    fast
                );
            });
        }
    }
}

#[test]
fn expanding_truncated_and_equal_radius_cones_retain_their_exact_distance() {
    for (base, top) in [(3, 1), (1, 3)] {
        let plan = primitive(Operation::Cone {
            base_radius: base.into(),
            top_radius: top.into(),
            height: 4.into(),
            pose: Pose::identity(),
        });
        let inside = query(&plan, &point([0, 0, 2]));
        assert_eq!(inside.location, Location::Inside);
        let bound = inside.distance.millimetres;
        let lo = Rational::from_f64_exact(bound.lo()).unwrap();
        let hi = Rational::from_f64_exact(bound.hi()).unwrap();
        let squared = Rational::parse("16/5").unwrap();
        assert!(&lo * &lo <= squared && &hi * &hi >= squared);
        distance(&plan, &point([2, 0, 2]), Location::Boundary, "0");
        unit(query(&plan, &point([0, 0, 4])).direction, [0, 0, 1]);
    }
    let plan = primitive(Operation::Cone {
        base_radius: 2.into(),
        top_radius: 2.into(),
        height: 4.into(),
        pose: Pose::identity(),
    });
    distance(&plan, &point([0, 0, 2]), Location::Inside, "-2");
}

#[test]
fn discarded_steps_and_tangent_hidden_supports_cannot_be_selected_as_faces() {
    let recipe = Recipe::new(
        "unused",
        vec![
            Step {
                id: id("a"),
                op: block(2, [0, 0, 0]),
            },
            Step {
                id: id("b"),
                op: sphere(2),
            },
        ],
        id("a"),
    );
    let plan = Plan::compile(&recipe).unwrap();
    let name = Origin {
        step: id("b"),
        role: Role::Side {
            segment: ElementId::new("wall").unwrap(),
            path: None,
        },
    };
    let face = plan.select(&name).unwrap();
    assert_eq!(
        plan.face_query(face, &point([2, 0, 0]), &mut Workspace::default())
            .unwrap()
            .state,
        FaceState::NotOnFace
    );
    let mut outer = sphere(2);
    let mut inner = sphere(1);
    let Operation::Sphere { pose, .. } = &mut inner else {
        unreachable!()
    };
    pose.translation.0[0] = 1.into();
    let plan = boolean(outer.clone(), inner, BooleanKind::Union);
    let face = plan.select(&name).unwrap();
    assert_eq!(
        plan.face_query(face, &point([2, 0, 0]), &mut Workspace::default())
            .unwrap()
            .state,
        FaceState::Indeterminate
    );
    let Operation::Sphere { pose, .. } = &mut outer else {
        unreachable!()
    };
    pose.translation.0[0] = 2.into();
    let plan = boolean(block(4, [0, 0, 0]), outer, BooleanKind::Intersect);
    assert_eq!(
        query(&plan, &point([4, 0, 0])).direction,
        Direction::Ambiguous
    );
}

#[test]
fn algebraic_rigid_motion_agrees_with_brep_and_preserves_bits_and_work() {
    fn run() -> (Vec<Vec<u64>>, u64) {
        use catallax_numbers::{Decider, Outcome, Real, sin_cos_degrees};
        use catallax_shapes::{Pose as GeometryPose, Rotation, Vector3};

        let mut dec = Decider::default();
        let [sine, cosine] =
            sin_cos_degrees(&Rational::from_ratio(45, 2).unwrap(), &mut dec).unwrap();
        let rotation =
            Rotation::from_real_quaternion(cosine, Real::zero(), sine, Real::zero(), &mut dec)
                .unwrap();
        let pose = GeometryPose::new(rotation, Vector3::zero(), &mut dec).unwrap();
        let axis = pose.apply_vector(&Vector3::new(Real::zero(), Real::zero(), Real::one()));
        assert!(
            dec.sign(&(&axis.x - &(&Real::from_i64(2).sqrt() / &Real::from_i64(2))))
                .is_coincident()
        );
        assert!(dec.compare(&axis.x, &axis.z).is_coincident());

        let recipe = Recipe::new(
            "solid",
            vec![Step {
                id: id("solid"),
                op: cylinder(1, 4),
            }],
            id("solid"),
        );
        let field = Plan::compile_checked(&recipe, &mut dec)
            .unwrap()
            .transformed(&pose, &mut dec)
            .unwrap();
        let joined = boolean(
            block(1, [0, 0, 0]),
            block(1, [0, -1, 0]),
            BooleanKind::Union,
        )
        .transformed(&pose, &mut dec)
        .unwrap();
        let mut first = None;
        for _ in 0..2 {
            let mut model = catallax_topology::Model::new();
            let body = catallax_topology::cylinder(&mut model, r(1), r(4), &pose).unwrap();
            let mut workspace = Workspace::default();
            dec.take_stats();
            let mut answers = Vec::new();
            for (p, expected, distance) in [
                ([1, 0, 1], Location::Inside, -1.0),
                ([0, 2, 0], Location::Outside, 1.0),
                ([1, 1, 1], Location::Boundary, 0.0),
            ] {
                let p = point(p);
                let query = field.proven_checked(&p, &mut workspace, &mut dec).unwrap();
                assert_eq!(query.location, expected);
                let Value::TrueDistance(value) = query.value else {
                    panic!("rigid motion lost the metric certificate")
                };
                assert!(value.lo() <= distance && value.hi() >= distance);
                let brep = catallax_boolean::classify_point(&mut model, body, &p, &mut dec);
                match expected {
                    Location::Inside => assert!(matches!(
                        brep,
                        Outcome::Proven(catallax_boolean::Containment::Inside)
                    )),
                    Location::Outside => assert!(matches!(
                        brep,
                        Outcome::Proven(catallax_boolean::Containment::Outside)
                    )),
                    Location::Boundary => {
                        assert!(matches!(brep, Outcome::ProvenCoincident(_)));
                        let Direction::Unique(normal) = query.direction else {
                            panic!("smooth side has no unique normal")
                        };
                        for (value, expected) in normal.into_iter().zip([0.0, 1.0, 0.0]) {
                            assert!(value.lo() <= expected && value.hi() >= expected);
                        }
                    }
                    _ => unreachable!(),
                }
                answers.push(bits(query));
            }
            // The joined blocks have a genuine sharp edge through the origin.
            // Its algebraic supports require ray germs, not a tolerance test.
            let query = joined
                .proven_checked(&point([0, 0, 0]), &mut workspace, &mut dec)
                .unwrap();
            assert_eq!(query.location, Location::Boundary);
            assert!(matches!(query.value, Value::Residual(_)));
            assert_eq!(query.direction, Direction::Ambiguous);
            answers.push(bits(query));
            let observed = (answers, dec.take_stats().work);
            if let Some(expected) = &first {
                assert_eq!(&observed, expected);
            } else {
                first = Some(observed);
            }
        }
        first.unwrap()
    }
    let expected = run();
    for count in [1, 2, 4] {
        std::thread::scope(|scope| {
            let workers: Vec<_> = (0..count).map(|_| scope.spawn(run)).collect();
            for worker in workers {
                assert_eq!(worker.join().unwrap(), expected);
            }
        });
    }
}

#[test]
fn checked_metric_enclosures_respect_the_callers_precision_policy() {
    use catallax_numbers::{BudgetPolicy, Decider, Real, sin_cos_degrees};
    use catallax_shapes::{Pose as GeometryPose, Rotation, Vector3};
    let id = StepId::new("sphere").unwrap();
    let recipe = Recipe::new(
        "sphere",
        vec![Step {
            id: id.clone(),
            op: Operation::Sphere {
                radius: 1.into(),
                pose: Pose::identity(),
            },
        }],
        id,
    );
    let mut dec = Decider::default();
    let [sine, cosine] = sin_cos_degrees(&Rational::from_ratio(45, 2).unwrap(), &mut dec).unwrap();
    let pose = GeometryPose::new(
        Rotation::from_real_quaternion(cosine, Real::zero(), sine, Real::zero(), &mut dec).unwrap(),
        Vector3::zero(),
        &mut dec,
    )
    .unwrap();
    let native = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let rotated = native.transformed(&pose, &mut dec).unwrap();
    let p = point([2, 0, 0]);
    let policy = BudgetPolicy {
        max_precision_bits: 64,
        ..BudgetPolicy::default()
    };
    let query = native
        .proven_checked(&p, &mut Workspace::default(), &mut Decider::new(policy))
        .unwrap();
    assert_eq!(query.location, Location::Outside);
    assert!(matches!(query.value, Value::TrueDistance(value) if value.contains(1.0)));
    let query = rotated
        .proven_checked(&p, &mut Workspace::default(), &mut dec)
        .unwrap();
    assert_eq!(query.location, Location::Outside);
    assert!(matches!(query.value, Value::TrueDistance(value) if value.contains(1.0)));
    let refusal = rotated
        .proven_checked(&p, &mut Workspace::default(), &mut Decider::new(policy))
        .unwrap_err();
    assert_eq!(refusal.code, "numbers.precision-limit");
}
