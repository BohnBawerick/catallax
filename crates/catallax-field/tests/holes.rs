//! Independent signs and trimmed ledge names for named drill tools.
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::Rational;
use catallax_recipe::{HolePart, Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::json;

fn point(x: &str, y: &str, z: &str) -> RPoint3 {
    RPoint3::new(
        Rational::parse(x).unwrap(),
        Rational::parse(y).unwrap(),
        Rational::parse(z).unwrap(),
    )
}
#[test]
fn short_countersink_ledge_is_boundary_not_an_untrimmed_cone_extension() {
    let recipe=Recipe::parse(&json!({"format":"catallax-recipe/3","name":"ledge","steps":[
        {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-5","-5","0"]}},
        {"id":"hole","op":"hole","target":"stock","at":[{"id":"centre","point":["0","0","0"],"direction":["0","0","1"]}],"diameter":"2","depth":{"blind":"1"},"tip":{"cone":"90"},"head":{"countersink":{"diameter":"8","angle":"90"}}}],"result":"hole"}).to_string()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (p, want) in [
        (point("2", "0", "1"), Location::Boundary),
        (point("2", "0", "1/2"), Location::Outside),
        (point("2", "0", "3/2"), Location::Inside),
        (point("0", "0", "3/2"), Location::Outside),
        (point("0", "0", "5/2"), Location::Inside),
    ] {
        assert_eq!(plan.proven(&p, &mut workspace).location, want, "{p:?}");
    }
    let floor = plan
        .faces()
        .iter()
        .find(|f| {
            matches!(
                f.origin.role,
                Role::Hole {
                    part: HolePart::Floor,
                    ..
                }
            )
        })
        .unwrap();
    let id = plan.select(&floor.origin).unwrap();
    assert_eq!(
        plan.face_query(id, &point("2", "0", "1"), &mut workspace)
            .unwrap()
            .state,
        FaceState::OnFace
    );
    assert_eq!(
        plan.face_query(id, &point("1/2", "0", "1"), &mut workspace)
            .unwrap()
            .state,
        FaceState::NotOnFace
    );
}

#[test]
fn mirrored_drill_cells_keep_their_cut_and_wrapped_floor() {
    let recipe=Recipe::parse(&json!({"format":"catallax-recipe/3","name":"mirrored drill","steps":[
        {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-5","-5","0"]}},
        {"id":"hole","op":"hole","target":"stock","at":[{"id":"centre","point":["2","0","0"],"direction":["0","0","1"]}],"diameter":"2","depth":{"blind":"1"},"tip":"flat"},
        {"id":"image","op":"mirror","body":"hole","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}
    ],"result":"image"}).to_string()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (p, want) in [
        (point("-2", "0", "1/2"), Location::Outside),
        (point("2", "0", "1/2"), Location::Inside),
        (point("-2", "0", "3/2"), Location::Inside),
    ] {
        assert_eq!(plan.proven(&p, &mut workspace).location, want, "{p:?}");
    }
    let floor=plan.faces().iter().find(|face|matches!(&face.origin.role,Role::Copy{instance:1,of} if matches!(of.role,Role::Hole{part:HolePart::Floor,..}))).unwrap();
    assert_eq!(
        plan.face_query(
            plan.select(&floor.origin).unwrap(),
            &point("-2", "0", "1"),
            &mut workspace
        )
        .unwrap()
        .state,
        FaceState::OnFace
    );
}
