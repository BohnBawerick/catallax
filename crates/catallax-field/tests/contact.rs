//! Closed-patch incidence must distinguish transverse rails from tangent donors.
use catallax_field::{Direction, FaceState, Plan, Workspace};
use catallax_numbers::Decider;
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};

fn pose(x: &str) -> Value {
    json!({"rotation":["1","0","0","0"],"translation":[x,"0","0"]})
}
fn recipe(steps: Vec<Value>) -> Recipe {
    let result = steps.last().unwrap()["id"].clone();
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"closed contact","steps":steps,"result":result}).to_string()).unwrap()
}

#[test]
fn transverse_curved_boolean_faces_keep_their_closed_intersection_circle() {
    for kind in ["union", "intersect", "subtract"] {
        let plan = Plan::compile_checked(
            &recipe(vec![
                json!({"id":"a","op":"sphere","radius":"5","pose":pose("0")}),
                json!({"id":"b","op":"sphere","radius":"5","pose":pose("6")}),
                json!({"id":"result","op":"boolean","kind":kind,"target":"a","tool":"b"}),
            ]),
            &mut Decider::default(),
        )
        .unwrap();
        let p = RPoint3::new(3.into(), 4.into(), 0.into());
        for face in plan.faces() {
            let id = plan.select(&face.origin).unwrap();
            let query = plan
                .closed_face_query_checked(
                    id,
                    &p,
                    &mut Workspace::default(),
                    &mut Decider::default(),
                )
                .unwrap();
            assert_eq!(query.state, FaceState::OnFace, "{kind}: {:?}", face.origin);
            let Direction::Unique(normal) = query.direction else {
                panic!("{kind}: {query:?}");
            };
            let b = face.origin.step.as_str() == "b";
            let sense = if b && kind == "subtract" { -1 } else { 1 };
            let expected = [
                catallax_numbers::Rational::from_ratio(if b { -3 * sense } else { 3 * sense }, 5)
                    .unwrap(),
                catallax_numbers::Rational::from_ratio(4 * sense, 5).unwrap(),
                0.into(),
            ];
            for (component, value) in normal.iter().zip(expected) {
                assert!(
                    catallax_numbers::Rational::from_f64_exact(component.lo()).unwrap() <= value
                );
                assert!(
                    catallax_numbers::Rational::from_f64_exact(component.hi()).unwrap() >= value
                );
            }
        }
    }
}

#[test]
fn buried_internal_tangency_cannot_donate_a_named_patch() {
    let plan = Plan::compile_checked(
        &recipe(vec![
            json!({"id":"outer","op":"sphere","radius":"2","pose":pose("0")}),
            json!({"id":"inner","op":"sphere","radius":"1","pose":pose("1")}),
            json!({"id":"result","op":"boolean","kind":"union","target":"outer","tool":"inner"}),
        ]),
        &mut Decider::default(),
    )
    .unwrap();
    let face = plan
        .faces()
        .iter()
        .find(|f| f.origin.step.as_str() == "inner")
        .unwrap();
    let query = plan
        .closed_face_query_checked(
            plan.select(&face.origin).unwrap(),
            &RPoint3::new(2.into(), 0.into(), 0.into()),
            &mut Workspace::default(),
            &mut Decider::default(),
        )
        .unwrap();
    assert_ne!(query.state, FaceState::OnFace);
}

#[test]
fn planar_fans_admit_closed_corners_with_each_faces_own_normal() {
    let plan = Plan::compile_checked(
        &recipe(vec![
            json!({"id":"a","op":"block","size":["2","2","2"],"pose":pose("0")}),
            json!({"id":"b","op":"block","size":["1","1","1"],"pose":pose("0")}),
            json!({"id":"result","op":"boolean","kind":"intersect","target":"a","tool":"b"}),
        ]),
        &mut Decider::default(),
    )
    .unwrap();
    let mut count = 0;
    for face in plan
        .faces()
        .iter()
        .filter(|f| f.origin.step.as_str() == "a")
    {
        let query = plan
            .closed_face_query_checked(
                plan.select(&face.origin).unwrap(),
                &RPoint3::origin(),
                &mut Workspace::default(),
                &mut Decider::default(),
            )
            .unwrap();
        if query.state == FaceState::OnFace {
            count += 1;
            assert!(matches!(query.direction, Direction::Unique(_)));
        }
    }
    assert_eq!(count, 3);
}

#[test]
fn fractional_bevel_and_cap_both_own_the_contact_circle() {
    let mut v: Value = serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
    v["steps"][2]["distance"] = json!("3/4");
    let plan = Plan::compile_checked(
        &Recipe::parse(&v.to_string()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    let p = RPoint3::new(
        catallax_numbers::Rational::parse("45/4").unwrap(),
        0.into(),
        0.into(),
    );
    let mut names = Vec::new();
    for face in plan.faces() {
        let q = plan
            .closed_face_query_checked(
                plan.select(&face.origin).unwrap(),
                &p,
                &mut Workspace::default(),
                &mut Decider::default(),
            )
            .unwrap();
        if q.state == FaceState::OnFace {
            names.push(face.origin.clone());
        }
    }
    assert_eq!(names.len(), 2, "{names:?}");
    assert!(names.iter().any(|n| matches!(n.role, Role::Bevel { .. })));
    assert!(
        names
            .iter()
            .any(|n| matches!(&n.role,Role::Side{segment,..} if segment.as_str()=="bottom"))
    );
}

#[test]
fn circular_offset_tangent_rails_keep_both_incident_names() {
    let plan = Plan::compile(&recipe(vec![
        json!({"id":"stock","op":"cylinder","radius":"2","height":"4","pose":pose("0")}),
        json!({"id":"parallel","op":"offset","body":"stock","distance":"1/5"}),
    ]))
    .unwrap();
    for (p, component) in [(["11/5", "0", "0"], 0), (["2", "0", "-1/5"], 2)] {
        let [x, y, z] = p.map(|v| catallax_numbers::Rational::parse(v).unwrap());
        let point = RPoint3::new(x, y, z);
        let mut names = Vec::new();
        for face in plan.faces() {
            let query = plan
                .closed_face_query_checked(
                    plan.select(&face.origin).unwrap(),
                    &point,
                    &mut Workspace::default(),
                    &mut Decider::default(),
                )
                .unwrap();
            if query.state != FaceState::OnFace {
                continue;
            }
            let Direction::Unique(normal) = query.direction else {
                panic!("{query:?}")
            };
            assert!(normal[component].contains(if component == 0 { 1.0 } else { -1.0 }));
            if !names.contains(&face.origin) {
                names.push(face.origin.clone());
            }
        }
        assert_eq!(names.len(), 2, "{p:?}: {names:?}");
        assert!(names.iter().any(|o| matches!(o.role, Role::Round { .. })));
        assert!(names.iter().any(|o| o.step.as_str() == "stock"));
    }
}

#[test]
fn tangent_offset_cells_cannot_donate_a_buried_named_patch() {
    let plan = Plan::compile(&recipe(vec![
        json!({"id":"outer","op":"sphere","radius":"3/2","pose":pose("0")}),
        json!({"id":"outer_parallel","op":"offset","body":"outer","distance":"1/2"}),
        json!({"id":"inner","op":"sphere","radius":"1/2","pose":pose("1")}),
        json!({"id":"inner_parallel","op":"offset","body":"inner","distance":"1/2"}),
        json!({"id":"result","op":"boolean","kind":"union","target":"outer_parallel","tool":"inner_parallel"}),
    ])).unwrap();
    let origin = &plan
        .faces()
        .iter()
        .find(|f| f.origin.step.as_str() == "inner")
        .unwrap()
        .origin;
    let query = plan
        .closed_face_query_checked(
            plan.select(origin).unwrap(),
            &RPoint3::new(2.into(), 0.into(), 0.into()),
            &mut Workspace::default(),
            &mut Decider::default(),
        )
        .unwrap();
    assert_ne!(query.state, FaceState::OnFace);
}
