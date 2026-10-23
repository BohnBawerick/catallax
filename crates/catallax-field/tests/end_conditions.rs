//! Independent first-half versus second-half cutting, including a reversed axis.
use catallax_field::{Location, Plan, Workspace};
use catallax_numbers::Rational;
use catallax_recipe::Recipe;
use catallax_shapes::RPoint3;
use serde_json::json;
#[test]
fn one_sided_tools_classify_both_sides_of_the_authoring_plane() {
    for (side, direction, removed_z) in [
        ("forward", "1", 1),
        ("backward", "1", -1),
        ("forward", "-1", -1),
        ("backward", "-1", 1),
    ] {
        let recipe=Recipe::parse(&json!({"format":"catallax-recipe/4","name":"half cut","steps":[
            {"id":"stock","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","-5"]}},
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"circle","id":"circle","centre":["5","5"],"radius":"1"}]]},
            {"id":"cut","op":"extrude_through","target":"stock","profile":"profile","direction":["0","0",direction],"side":side}],"result":"cut"}).to_string()).unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        let mut workspace = Workspace::default();
        for (z, expected) in [
            (removed_z, Location::Outside),
            (-removed_z, Location::Inside),
            (0, Location::Boundary),
        ] {
            let point = RPoint3::new(
                Rational::from_i64(5),
                Rational::from_i64(5),
                Rational::from_i64(z),
            );
            assert_eq!(
                plan.proven(&point, &mut workspace).location,
                expected,
                "{side} along {direction} at {z}"
            );
        }
    }
}

fn point(x: &str, y: &str, z: &str) -> RPoint3 {
    let [x, y, z] = [x, y, z].map(|s| Rational::parse(s).unwrap());
    RPoint3::new(x, y, z)
}

#[test]
fn stepped_first_hits_and_curtain_share_one_end_origin_across_distinct_supports() {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/stepped-limit.json"
    ))
    .unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let mut workspace = Workspace::default();
    for (p, expected) in [
        (point("3/4", "3", "4"), Location::Inside),
        (point("3/4", "3", "6"), Location::Outside),
        (point("5/4", "3", "9"), Location::Inside),
        (point("5/4", "3", "11"), Location::Outside),
        (point("1", "3", "7"), Location::Boundary),
        (point("1", "3", "2"), Location::Inside),
        (point("1/2", "3", "7"), Location::Outside),
    ] {
        assert_eq!(plan.proven(&p, &mut workspace).location, expected, "{p:?}");
    }
    let end = catallax_recipe::Origin {
        step: catallax_recipe::StepId::new("cap").unwrap(),
        role: catallax_recipe::Role::End,
    };
    let selected = plan.select(&end).unwrap();
    for p in [
        point("3/4", "3", "5"),
        point("5/4", "3", "10"),
        point("1", "3", "7"),
    ] {
        assert_eq!(
            plan.face_query(selected, &p, &mut workspace).unwrap().state,
            catallax_field::FaceState::OnFace,
            "{p:?}"
        );
    }
    assert_eq!(
        plan.face_query(selected, &point("1", "3", "2"), &mut workspace)
            .unwrap()
            .state,
        catallax_field::FaceState::NotOnFace
    );
}
