//! A face-mounted cut must share the supporting plane, not a rounded copy of it.
use catallax_app::soliddraft::convert;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Number, Operation, Role};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation};
use serde_json::{Value, json};

fn recipe(steps: Vec<Value>, result: &str) -> String {
    json!({"format":"soliddraft-recipe/1","name":"angled face cut","steps":steps,"result":result})
        .to_string()
}

fn polygon(points: &[[f64; 2]]) -> Value {
    json!({"outer": (0..points.len()).map(|i| json!({
        "type":"line", "from":points[i], "to":points[(i+1)%points.len()]
    })).collect::<Vec<_>>()})
}

fn cut_profile() -> Value {
    json!({"id":"CutSketch","op":"profile",
    "plane":{"origin":[58.52917258875233,0.0,73.20363419956955],
        "normal":[0.6244755519662776,0.0,0.7810443553322775],"x_dir":[0.0,1.0,0.0]},
    "loops":polygon(&[
        [-52.195125579833984,1.726492166519165],
        [-69.71279914190637,-8.387341380119322],
        [-69.71279914190639,-28.615008473396294],
        [-52.195125579833984,-38.72884202003479],
        [-34.6774520177616,-28.6150084733963],
        [-34.6774520177616,-8.387341380119322],
    ])})
}

#[test]
fn identity_copy_placements_keep_rational_face_axes() {
    let normal = RVector3::new(2.into(), 0.into(), 3.into());
    let x = RVector3::new(0.into(), 1.into(), 0.into());
    let plane = catallax_shapes::Plane {
        origin: RPoint3::new(3.into(), 0.into(), 2.into()).into(),
        normal: normal.clone().into(),
        ref_dir: x.clone().into(),
    };
    for scale in [1, -2] {
        let mut placement = catallax_recipe::Pose::identity();
        placement.rotation[0] = Number::from(scale);
        let pose = catallax_features::face_placement(
            std::slice::from_ref(&plane),
            &normal,
            &x,
            &Rational::zero(),
            &placement,
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(
            pose.rotation.matrix()[1]
                .each_ref()
                .map(|v| v.as_rational()),
            [Some(1.into()), Some(0.into()), Some(0.into())]
        );
        assert_eq!(pose.translation.y.as_rational(), Some(0.into()));
    }
}

#[test]
fn rounded_frame_does_not_preserve_the_supporting_plane() {
    let converted = convert(&recipe(
        vec![
            cut_profile(),
            json!({"id":"Tool","op":"extrude","profile":"CutSketch","distance":150}),
        ],
        "Tool",
    ))
    .unwrap();
    let Operation::Sketch { plane, .. } = &converted.recipe.steps[0].op else {
        panic!()
    };
    let [w, x, y, z] = plane.rotation.each_ref().map(|v| v.value().clone());
    let [tx, ty, tz] = plane.translation.0.each_ref().map(|v| v.value().clone());
    let pose = Pose::from_rational(
        Rotation::from_quaternion(w, x, y, z).unwrap(),
        RVector3::new(tx, ty, tz),
    );
    let source_normal = RVector3::new(
        Rational::from_i64(120),
        Rational::zero(),
        Rational::parse("150.08645629882812").unwrap(),
    );
    let source_point = RPoint3::new(source_normal.z.clone(), Rational::zero(), Rational::zero());
    let p = pose
        .apply_rational_point(&RPoint3::new(
            Rational::from_i64(-52),
            Rational::from_i64(-18),
            Rational::zero(),
        ))
        .unwrap();
    let residual = source_normal.dot(&p.minus(&source_point));
    eprintln!(
        "converted frame: {plane:?}\nexact residual at (-52,-18): {residual}\nresidual estimate: {}",
        residual.to_f64_nearest()
    );
    assert!(!residual.is_zero());
}

#[test]
fn face_mounted_cut_opens_supporting_plane() {
    let n = [0.6244755519662776, 0.0, 0.7810443553322775];
    let steps = vec![
        json!({"id":"BaseSketch","op":"profile","plane":{"origin":[0,0,0],"normal":[0,-1,0],"x_dir":[1,0,0]},
            "loops":polygon(&[[0.,0.],[150.08645629882812,0.],[0.,120.]])}),
        json!({"id":"Boss","op":"extrude","profile":"BaseSketch","distance":105.5}),
        cut_profile(),
        json!({"id":"Tool","op":"extrude","profile":"CutSketch","distance":150,"direction":n.map(|v| -v)}),
        json!({"id":"Cut","op":"boolean","kind":"subtract","target":"Boss","tool":"Tool"}),
    ];
    let mut converted = convert(&recipe(steps, "Cut")).unwrap();
    if let Operation::Sketch { plane, .. } = &mut converted.recipe.steps[0].op {
        plane.rotation = [1, 1, 0, 0].map(Number::from);
    }
    if let Operation::Sketch { plane, .. } = &mut converted.recipe.steps[2].op {
        *plane = catallax_recipe::Pose::identity();
    }
    converted.recipe.steps[3].op = Operation::ExtrudeOnFace {
        profile: converted.steps["CutSketch"].clone(),
        body: converted.steps["Boss"].clone(),
        face: catallax_recipe::OriginPattern {
            step: converted.steps["Boss"].clone(),
            role: Some(Role::Side {
                segment: catallax_recipe::ElementId::new("l0s1").unwrap(),
                path: None,
            }),
        },
        normal_direction: catallax_recipe::Vector3(
            n.map(|v| Number::parse(&v.to_string()).unwrap()),
        ),
        x_direction: catallax_recipe::Vector3([0, 1, 0].map(Number::from)),
        offset: Number::from(0),
        distance: Number::from(-150),
        placement: catallax_recipe::Pose::identity(),
    };
    let (built, work) =
        catallax_numbers::with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            catallax_features::evaluate(&converted.recipe, &mut Decider::default())
        });
    eprintln!("face-mounted cut work: {work:?}");
    let built = built.unwrap();
    let tool = &converted.steps["Tool"];
    let start_caps: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| {
            built
                .model
                .face_origins(*f)
                .iter()
                .any(|o| &o.step == tool && o.role == Role::Start)
        })
        .collect();
    assert!(
        start_caps.is_empty(),
        "the cut left its Start cap behind the supporting face: {start_caps:?}"
    );
}
