//! Body copies retain exact definitions, regularisation and per-copy origins.
use catallax_features::{Built, evaluate};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{Recipe, Role, StepId};
use serde_json::{Value, json};

fn recipe(mut steps: Vec<Value>, operation: Value) -> Recipe {
    steps.push(operation);
    Recipe::parse(&json!({"format":catallax_recipe::FORMAT,"name":"copies","steps":steps,"result":steps.last().unwrap()["id"]}).to_string()).unwrap()
}
fn seed() -> Vec<Value> {
    vec![
        json!({"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
            {"kind":"line","id":"a","to":["1","0"]},{"kind":"line","id":"b","to":["1","2"]},{"kind":"line","id":"c","to":["0","2"]},{"kind":"line","id":"d","to":["0","0"]}
        ]]}),
        json!({"id":"seed","op":"extrude","profile":"profile","distance":"3"}),
    ]
}
fn build(recipe: &Recipe) -> Built {
    let built = evaluate(recipe, &mut Decider::default()).unwrap();
    assert_eq!(built.model.check(built.result), vec![]);
    built
}
fn volume(b: &Built) -> catallax_measure::Enclosure {
    catallax_measure::mass_properties(
        &b.model,
        b.result,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    .volume
}
fn step(b: &Built) -> Vec<u8> {
    let mut out = vec![];
    let report = catallax_step::write_step(
        &b.model,
        &[b.result],
        &Default::default(),
        &mut Decider::default(),
        &mut out,
    )
    .unwrap();
    assert!(report.marks.is_empty());
    assert!(report.value.max_deviation <= Rational::parse("1e-7").unwrap());
    out
}

#[test]
fn linear_disjoint_touching_overlapping_and_coincident_copies() {
    for (spacing, expected, lumps) in [
        ("3", "18", 3),
        ("1", "18", 1),
        ("1/2", "12", 1),
        ("0", "6", 1),
    ] {
        let r = recipe(
            seed(),
            json!({"id":"row","op":"pattern","body":"seed","layout":{"linear":{"vector":[spacing,"0","0"],"count":3}}}),
        );
        let b = build(&r);
        assert_eq!(volume(&b).lo, Rational::parse(expected).unwrap());
        assert_eq!(volume(&b).hi, Rational::parse(expected).unwrap());
        assert_eq!(b.model.counts(b.result).lumps, lumps);
        let origins = b
            .model
            .faces_of_body(b.result)
            .into_iter()
            .flat_map(|f| b.model.face_origins(f))
            .collect::<Vec<_>>();
        for k in 1..3 {
            assert!(
                origins
                    .iter()
                    .any(|o| matches!(o.role,Role::Copy{instance,..} if instance == k))
            );
        }
        assert_eq!(step(&b), step(&build(&r)));
    }
}

#[test]
fn circular_uses_incremental_angle_and_a_non_origin_nonunit_axis() {
    let r = recipe(
        seed(),
        json!({"id":"ring","op":"pattern","body":"seed","layout":{"circular":{"axis":{"point":["10","0","0"],"direction":["0","0","7"]},"angle":"90","count":4}}}),
    );
    let b = build(&r);
    assert_eq!(volume(&b).lo, Rational::from_i64(24));
    assert_eq!(b.model.counts(b.result).lumps, 4);
    let seed = b.bodies[&StepId::new("seed").unwrap()];
    let source = b.model.faces_of_body(seed);
    for k in 1..4 {
        let map = catallax_features::pattern_placement(
            match &r.steps.last().unwrap().op {
                catallax_recipe::Operation::Pattern { layout, .. } => layout,
                _ => unreachable!(),
            },
            k,
            &mut Decider::default(),
        )
        .unwrap();
        for f in &source {
            let expected = map.surface(b.model.geom().surface(b.model.face(*f).surface));
            let origin = catallax_features::copy_origin(
                &StepId::new("ring").unwrap(),
                k,
                &b.model.face_origins(*f)[0],
            );
            let copy = b
                .model
                .faces_of_body(b.result)
                .into_iter()
                .find(|f| b.model.face_origins(*f).contains(&origin))
                .unwrap();
            assert!(
                b.model
                    .geom()
                    .surface(b.model.face(copy).surface)
                    .same_definition(&expected)
            );
        }
    }
    step(&b);
}

#[test]
fn mirror_oblique_plane_preserves_outward_normals_and_reverses_loops() {
    for keep in [false, true] {
        let r = recipe(
            seed(),
            json!({"id":"mirror","op":"mirror","body":"seed","plane":{"point":["5","0","0"],"normal":["2","3","0"]},"keep_original":keep}),
        );
        let b = build(&r);
        assert_eq!(volume(&b).lo, Rational::from_i64(if keep { 12 } else { 6 }));
        assert_eq!(b.model.counts(b.result).lumps, if keep { 2 } else { 1 });
        let origins = b
            .model
            .faces_of_body(b.result)
            .into_iter()
            .flat_map(|f| b.model.face_origins(f))
            .collect::<Vec<_>>();
        assert_eq!(
            origins
                .iter()
                .filter(|o| matches!(o.role, Role::Copy { instance: 1, .. }))
                .count(),
            6
        );
        assert_eq!(step(&b), step(&build(&r)));
    }
}

#[test]
fn revolved_tools_pattern_and_cut_preserve_each_cylinder_name() {
    let mut steps = vec![
        json!({"id":"base","op":"block","size":["12","4","6"],"pose":{"rotation":["1","0","0","0"],"translation":["-2","0","-3"]}}),
    ];
    steps.push(json!({"id":"section","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
        {"kind":"line","id":"radius","to":["1","0"]},{"kind":"line","id":"wall","to":["1","4"]},{"kind":"line","id":"top","to":["0","4"]},{"kind":"line","id":"axis","to":["0","0"]}
    ]]}));
    steps.push(json!({"id":"seed","op":"revolve","profile":"section","axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":"360"}));
    steps.push(json!({"id":"row","op":"pattern","body":"seed","layout":{"linear":{"vector":["4","0","0"],"count":3}}}));
    let b = build(&recipe(
        steps,
        json!({"id":"cut","op":"boolean","kind":"subtract","target":"base","tool":"row"}),
    ));
    let expected = Ball::from_rational(&Rational::from_i64(288), 192)
        - Ball::pi(192) * Ball::from_rational(&Rational::from_i64(12), 192);
    let v = volume(&b);
    assert!(v.lo <= expected.lo() && v.hi >= expected.hi());
    let origins = b
        .model
        .faces_of_body(b.result)
        .into_iter()
        .flat_map(|f| b.model.face_origins(f))
        .collect::<Vec<_>>();
    for k in 1..3 {
        assert!(
            origins
                .iter()
                .any(|o| matches!(o.role,Role::Copy{instance,..} if instance == k))
        );
    }
    step(&b);
}

#[test]
fn impossible_placements_fail_even_for_one_or_empty_copies() {
    for (op, code) in [
        (
            json!({"id":"bad","op":"pattern","body":"seed","layout":{"linear":{"vector":["0","0","0"],"count":0}}}),
            "pattern.zero-count",
        ),
        (
            json!({"id":"bad","op":"pattern","body":"seed","layout":{"circular":{"axis":{"point":["0","0","0"],"direction":["0","0","0"]},"angle":"0","count":1}}}),
            "pattern.zero-axis",
        ),
        (
            json!({"id":"bad","op":"mirror","body":"seed","plane":{"point":["0","0","0"],"normal":["0","0","0"]},"keep_original":false}),
            "mirror.zero-normal",
        ),
    ] {
        for body in ["seed", "empty"] {
            let mut source = seed();
            source.push(json!({"id":"empty","op":"boolean","kind":"subtract","target":"seed","tool":"seed"}));
            let mut operation = op.clone();
            operation["body"] = json!(body);
            let error = evaluate(&recipe(source, operation), &mut Decider::default()).unwrap_err();
            assert!(
                matches!(error,catallax_features::BuildError::Geometry{ref refusal,..} if refusal.code == code)
            );
        }
    }
}

#[test]
fn curved_copies_retain_source_supports_mass_names_and_step_bytes() {
    let frame = json!({"rotation":["1","0","0","0"],"translation":["3","0","0"]});
    let axis = json!({"point":["3","0","0"],"direction":["0","1","0"]});
    let circle = json!({"id":"profile","op":"sketch","plane":frame,"loops":[[
        {"kind":"circle","id":"wall","centre":["0","0"],"radius":"1/2"}
    ]]});
    let section = json!({"id":"profile","op":"sketch","plane":frame,"loops":[[
        {"kind":"line","id":"a","to":["2","0"]},{"kind":"line","id":"b","to":["2","1"]},{"kind":"line","id":"c","to":["1","1"]},{"kind":"line","id":"d","to":["1","0"]}
    ]]});
    let sources = [
        (
            vec![
                circle,
                json!({"id":"seed","op":"extrude","profile":"profile","distance":"2"}),
            ],
            Rational::parse("1/2").unwrap(),
        ),
        (
            vec![
                section.clone(),
                json!({"id":"seed","op":"revolve","profile":"profile","axis":axis,"angle":"360"}),
            ],
            Rational::from_i64(3),
        ),
        (
            vec![
                section,
                json!({"id":"seed","op":"revolve","profile":"profile","axis":axis,"angle":"90"}),
            ],
            Rational::parse("3/4").unwrap(),
        ),
    ];
    for (source, coefficient) in sources {
        for (operation, count) in [
            (
                json!({"id":"copies","op":"pattern","body":"seed","layout":{"linear":{"vector":["5","0","0"],"count":3}}}),
                3,
            ),
            (
                json!({"id":"copies","op":"pattern","body":"seed","layout":{"circular":{"axis":{"point":["0","0","0"],"direction":["0","0","7"]},"angle":"90","count":4}}}),
                4,
            ),
            (
                json!({"id":"copies","op":"mirror","body":"seed","plane":{"point":["0","0","0"],"normal":["2","3","0"]},"keep_original":false}),
                1,
            ),
        ] {
            let recipe = recipe(source.clone(), operation);
            let built = build(&recipe);
            let source_body = built.bodies[&StepId::new("seed").unwrap()];
            let source_faces = built.model.faces_of_body(source_body);
            let expected = Ball::pi(192)
                * Ball::from_rational(&(&coefficient * &Rational::from_i64(count)), 192);
            let volume = volume(&built);
            assert!(volume.lo <= expected.lo() && volume.hi >= expected.hi());
            assert!(volume.width() <= Rational::parse("1e-9").unwrap());
            assert_eq!(built.model.counts(built.result).lumps, count as usize);
            assert_eq!(
                built.model.counts(built.result).faces,
                source_faces.len() * count as usize
            );
            let operation = &recipe.steps.last().unwrap().op;
            let instances = if matches!(operation, catallax_recipe::Operation::Mirror { .. }) {
                1..2
            } else {
                0..count as u32
            };
            for instance in instances {
                let map = match operation {
                    catallax_recipe::Operation::Pattern { layout, .. } => {
                        catallax_features::pattern_placement(
                            layout,
                            instance,
                            &mut Decider::default(),
                        )
                        .unwrap()
                    }
                    catallax_recipe::Operation::Mirror { plane, .. } => {
                        catallax_features::mirror_placement(plane).unwrap()
                    }
                    _ => unreachable!(),
                };
                for source in &source_faces {
                    let origin = catallax_features::copy_origin(
                        &StepId::new("copies").unwrap(),
                        instance,
                        &built.model.face_origins(*source)[0],
                    );
                    let copy = built
                        .model
                        .faces_of_body(built.result)
                        .into_iter()
                        .find(|f| built.model.face_origins(*f).contains(&origin))
                        .unwrap();
                    let expected = map.surface(
                        built
                            .model
                            .geom()
                            .surface(built.model.face(*source).surface),
                    );
                    assert!(
                        built
                            .model
                            .geom()
                            .surface(built.model.face(copy).surface)
                            .same_definition(&expected)
                    );
                }
            }
            assert_eq!(step(&built), step(&build(&recipe)));
        }
    }
}
