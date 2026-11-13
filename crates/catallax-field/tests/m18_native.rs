//! Native field motions, circular blends and their explicit remaining limits.
use catallax_field::{FaceState, Location, Plan, Value, Workspace};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::RPoint3;
use serde_json::{Value as Json, json};

fn pose(x: &str, y: &str, z: &str) -> Json {
    json!({"rotation":["1","0","0","0"],"translation":[x,y,z]})
}
fn rectangle() -> Json {
    json!([[{"kind":"line","id":"bottom","to":["4","0"]},{"kind":"line","id":"right","to":["4","3"]},{"kind":"line","id":"top","to":["0","3"]},{"kind":"line","id":"left","to":["0","0"]}]])
}
fn recipe(steps: Vec<Json>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"m18-native","steps":steps,"result":result})
            .to_string(),
    )
    .unwrap()
}
fn p(x: &str, y: &str, z: &str) -> RPoint3 {
    RPoint3::new(
        Rational::parse(x).unwrap(),
        Rational::parse(y).unwrap(),
        Rational::parse(z).unwrap(),
    )
}
fn location(plan: &Plan, point: &RPoint3) -> Location {
    plan.proven_checked(point, &mut Workspace::default(), &mut Decider::default())
        .unwrap()
        .location
}

#[test]
fn straight_sweep_keeps_world_displacement_and_path_names() {
    let r = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0","0","0"),"loops":rectangle()}),
        json!({"id":"path","op":"path","start":["100","200","300"],"segments":[{"kind":"line","id":"run","to":["101","202","303"]}]}),
        json!({"id":"sweep","op":"sweep","profile":"profile","path":"path","orientation":"fixed"}),
    ]);
    let plan = Plan::compile(&r).unwrap();
    assert_eq!(location(&plan, &p("5/2", "5/2", "3/2")), Location::Inside);
    assert_eq!(location(&plan, &p("0", "3/2", "0")), Location::Boundary);
    assert_eq!(location(&plan, &p("100", "200", "300")), Location::Outside);
    let origin = Origin {
        step: StepId::new("sweep").unwrap(),
        role: Role::Side {
            segment: catallax_recipe::ElementId::new("right").unwrap(),
            path: Some(catallax_recipe::ElementId::new("run").unwrap()),
        },
    };
    assert_eq!(
        plan.face_query(
            plan.select(&origin).unwrap(),
            &p("9/2", "5/2", "3/2"),
            &mut Workspace::default()
        )
        .unwrap()
        .state,
        FaceState::OnFace
    );
    let mesh =
        catallax_field::display::mesh_plan(&plan, &Default::default(), &mut Decider::default())
            .unwrap();
    assert!(!mesh.triangles.is_empty());
}

#[test]
fn fixed_line_chain_unites_spans_and_rm_corner_refuses_with_join() {
    let mut r = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0","0","0"),"loops":rectangle()}),
        json!({"id":"path","op":"path","start":["0","0","0"],"segments":[{"kind":"line","id":"a","to":["0","0","2"]},{"kind":"line","id":"b","to":["1","0","4"]}]}),
        json!({"id":"sweep","op":"sweep","profile":"profile","path":"path","orientation":"fixed"}),
    ]);
    let plan = Plan::compile(&r).unwrap();
    assert_eq!(location(&plan, &p("2", "1", "2")), Location::Inside);
    assert_eq!(location(&plan, &p("3", "1", "3")), Location::Inside);
    if let catallax_recipe::Operation::Sweep { orientation, .. } = &mut r.steps[2].op {
        *orientation = catallax_recipe::SweepOrientation::RotationMinimising;
    }
    let e = Plan::compile(&r).unwrap_err();
    assert_eq!(e.code, "field.sweep-corner-unsupported");
    assert!(e.detail.contains("join 1"));
}

#[test]
fn cardinal_arc_sweep_uses_the_field_revolution_and_refuses_fixed_motion() {
    let mut r = recipe(vec![
        json!({"id":"profile","op":"sketch","plane":pose("0","0","0"),"loops":[[{"kind":"circle","id":"tube","centre":["4","0"],"radius":"1"}]]}),
        json!({"id":"path","op":"path","start":["4","0","0"],"segments":[{"kind":"arc","id":"turn","centre":["0","0","0"],"to":["0","0","-4"],"axis":["0","1","0"]}]}),
        json!({"id":"sweep","op":"sweep","profile":"profile","path":"path","orientation":"rotation_minimising"}),
    ]);
    let plan = Plan::compile(&r).unwrap();
    assert_eq!(location(&plan, &p("3", "0", "-4")), Location::Boundary);
    assert_eq!(location(&plan, &p("12/5", "0", "-16/5")), Location::Inside);
    assert_eq!(location(&plan, &p("0", "0", "0")), Location::Outside);
    if let catallax_recipe::Operation::Sweep { orientation, .. } = &mut r.steps[2].op {
        *orientation = catallax_recipe::SweepOrientation::Fixed;
    }
    assert_eq!(
        Plan::compile(&r).unwrap_err().code,
        "field.sweep-motion-unsupported"
    );
}

#[test]
fn translated_lofts_build_and_different_circles_reduce_to_a_cone() {
    for shape in ["ruled", "smooth"] {
        let r = recipe(vec![
            json!({"id":"a","op":"sketch","plane":pose("0","0","0"),"loops":rectangle()}),
            json!({"id":"b","op":"sketch","plane":pose("1","0","3"),"loops":rectangle()}),
            json!({"id":"c","op":"sketch","plane":pose("2","0","6"),"loops":rectangle()}),
            json!({"id":"loft","op":"loft","sections":[{"profile":"a","start":["bottom"]},{"profile":"b","start":["bottom"]},{"profile":"c","start":["bottom"]}],"shape":shape}),
        ]);
        let plan = Plan::compile(&r).unwrap();
        assert_eq!(location(&plan, &p("3", "1", "3")), Location::Inside);
        assert_eq!(location(&plan, &p("0", "1", "3")), Location::Outside);
    }
    let r = recipe(vec![
        json!({"id":"a","op":"sketch","plane":pose("2","3","4"),"loops":[[{"kind":"circle","id":"lower","centre":["1","0"],"radius":"2"}]]}),
        json!({"id":"b","op":"sketch","plane":pose("2","3","10"),"loops":[[{"kind":"circle","id":"upper","centre":["1","0"],"radius":"4"}]]}),
        json!({"id":"loft","op":"loft","sections":[{"profile":"a","start":["lower"]},{"profile":"b","start":["upper"]}],"shape":"smooth"}),
    ]);
    let plan = Plan::compile(&r).unwrap();
    assert_eq!(location(&plan, &p("6", "3", "7")), Location::Boundary);
    assert_eq!(location(&plan, &p("3", "3", "7")), Location::Inside);
    assert_eq!(location(&plan, &p("7", "3", "7")), Location::Outside);
    assert_eq!(
        plan.faces()
            .iter()
            .filter(
                |f| matches!(&f.origin.role,Role::Side { segment,.. } if segment.as_str()=="lower")
            )
            .count(),
        1
    );
}

fn blend(kind: &str) -> Recipe {
    recipe(vec![
        json!({"id":"a","op":"block","size":["20","40","40"],"pose":pose("-20","-20","-20")}),
        json!({"id":"b","op":"block","size":["40","20","40"],"pose":pose("-20","-20","-20")}),
        json!({"id":"blend","op":"smooth_blend","kind":kind,"target":"a","tool":"b","radius":"5"}),
    ])
}

#[test]
fn circular_blends_use_true_distance_and_name_the_changed_patch() {
    for (kind, boundary, inside, outside) in [
        (
            "union",
            p("2", "1", "0"),
            p("1", "1", "0"),
            p("4", "4", "0"),
        ),
        (
            "intersect",
            p("-2", "-1", "0"),
            p("-4", "-4", "0"),
            p("-1", "-1", "0"),
        ),
        (
            "subtract",
            p("-2", "1", "0"),
            p("-4", "4", "0"),
            p("-1", "1", "0"),
        ),
    ] {
        let plan = Plan::compile(&blend(kind)).unwrap();
        assert_eq!(location(&plan, &boundary), Location::Boundary, "{kind}");
        assert_eq!(location(&plan, &inside), Location::Inside, "{kind}");
        assert_eq!(location(&plan, &outside), Location::Outside, "{kind}");
        assert_eq!(
            plan.smooth_origins_checked(&boundary, &mut Decider::default())
                .unwrap(),
            vec![Origin {
                step: StepId::new("blend").unwrap(),
                role: Role::Smooth
            }]
        );
        let q = plan
            .proven_checked(&outside, &mut Workspace::default(), &mut Decider::default())
            .unwrap();
        assert!(matches!(q.value, Value::Residual(_)));
        assert!(
            q.distance.millimetres.lo()
                <= (outside.minus(&boundary).norm_squared().to_f64_nearest()).sqrt()
        );
        assert!(
            !plan
                .fast([0., 0., 0.], &mut Workspace::default())
                .unwrap()
                .true_distance_formula
        );
    }
}

#[test]
fn unchanged_blend_boundary_keeps_only_the_exposed_stock_name() {
    let plan = Plan::compile(&blend("union")).unwrap();
    let names = plan
        .smooth_origins_checked(&p("0", "10", "0"), &mut Decider::default())
        .unwrap();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0].step, StepId::new("a").unwrap());
    assert!(
        plan.smooth_origins_checked(&p("0", "0", "0"), &mut Decider::default())
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        catallax_field::display::mesh_plan(&plan, &Default::default(), &mut Decider::default())
            .unwrap_err()
            .code,
        "smooth.boundary-unsupported"
    );
    assert_eq!(
        catallax_field::printing::classify_box(
            &plan,
            &std::array::from_fn(|_| Rational::from_i64(-1)),
            &std::array::from_fn(|_| Rational::from_i64(1)),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "smooth.boundary-unsupported"
    );
}

#[test]
fn composite_distance_and_nonpositive_radius_refuse_without_residual_substitution() {
    let mut r = blend("union");
    if let catallax_recipe::Operation::SmoothBlend { radius, .. } = &mut r.steps[2].op {
        *radius = 0.into();
    }
    assert_eq!(
        Plan::compile(&r).unwrap_err().code,
        "smooth.non-positive-radius"
    );
    let mut r = blend("union");
    r.steps.insert(
        2,
        catallax_recipe::Step {
            id: StepId::new("composite").unwrap(),
            op: catallax_recipe::Operation::Boolean {
                kind: catallax_recipe::BooleanKind::Union,
                target: StepId::new("a").unwrap(),
                tool: StepId::new("b").unwrap(),
            },
        },
    );
    if let catallax_recipe::Operation::SmoothBlend { target, .. } = &mut r.steps[3].op {
        *target = StepId::new("composite").unwrap();
    }
    assert_eq!(
        Plan::compile(&r).unwrap_err().code,
        "smooth.distance-unsupported"
    );
}

#[test]
fn smooth_answers_and_work_are_identical_at_one_two_and_four_workers() {
    let run = |threads| {
        catallax_field::workers::with_workers(threads, || {
            let (q, w) = with_operation_budget(20_000_000, || {
                let plan = Plan::compile(&blend("union")).unwrap();
                [p("2", "1", "0"), p("1", "1", "0"), p("4", "4", "0")].map(|p| {
                    plan.proven_checked(&p, &mut Workspace::default(), &mut Decider::default())
                        .unwrap()
                })
            });
            assert!(!w.exhausted);
            (q, w.used)
        })
        .unwrap()
    };
    let baseline = run(1);
    assert_eq!(run(2), baseline);
    assert_eq!(run(4), baseline);
}
