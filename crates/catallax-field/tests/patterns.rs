//! Coincident copies must not inherit aliases from unrelated body snapshots.
use catallax_field::{FaceState, Plan, Workspace};
use catallax_numbers::Rational;
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::json;
fn plan(op: serde_json::Value) -> Plan {
    let mut steps = vec![
        json!({"id":"seed","op":"block","size":["1","1","1"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}),
    ];
    if let serde_json::Value::Array(operations) = op {
        steps.extend(operations);
    } else {
        steps.push(op);
    }
    let recipe = Recipe::parse(&json!({"format":catallax_recipe::FORMAT,"name":"coincident-copy","steps":steps,"result":"copies"}).to_string()).unwrap();
    Plan::compile(&recipe).unwrap()
}
fn top(plan: &Plan) -> Vec<catallax_recipe::Origin> {
    let p = RPoint3::new(
        Rational::parse("1/2").unwrap(),
        Rational::parse("1/2").unwrap(),
        1.into(),
    );
    let mut workspace = Workspace::default();
    plan.faces()
        .iter()
        .filter_map(|face| {
            let id = plan.select(&face.origin).unwrap();
            (plan.face_query(id, &p, &mut workspace).unwrap().state == FaceState::OnFace)
                .then(|| face.origin.clone())
        })
        .collect()
}
#[test]
fn coincident_pattern_uses_the_seed_snapshot_for_every_copy() {
    let plan = plan(
        json!({"id":"copies","op":"pattern","body":"seed","layout":{"linear":{"vector":["0","0","0"],"count":3}}}),
    );
    let names = top(&plan);
    assert_eq!(names.len(), 3);
    assert_eq!(names.iter().filter(|o| o.role == Role::End).count(), 1);
    for k in [1, 2] {
        assert!(names.iter().any(
            |o| matches!(&o.role,Role::Copy{instance,of} if *instance == k && of.role == Role::End)
        ));
    }
}
#[test]
fn symmetric_mirror_without_original_selects_only_wrapped_origins() {
    let plan = plan(
        json!({"id":"copies","op":"mirror","body":"seed","plane":{"point":["1/2","0","0"],"normal":["1","0","0"]},"keep_original":false}),
    );
    let names = top(&plan);
    assert_eq!(names.len(), 1);
    assert!(matches!(&names[0].role,Role::Copy{instance:1,of} if of.role == Role::End));
}

#[test]
fn nested_mirrors_do_not_reactivate_the_unused_seed_alias() {
    let plan = plan(json!([
        {"id":"image","op":"mirror","body":"seed","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false},
        {"id":"copies","op":"mirror","body":"image","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}
    ]));
    let names = top(&plan);
    assert_eq!(names.len(), 1);
    assert!(matches!(&names[0].role,Role::Copy{instance:1,of}
        if of.step.as_str() == "image" && matches!(&of.role,Role::Copy{instance:1,of} if of.role == Role::End)));
}
