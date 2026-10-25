//! Native field compilation keeps the first toroidal branch and face ownership.
use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::Rational;
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};
fn point(x: &str, z: &str) -> RPoint3 {
    RPoint3::new(
        Rational::parse(x).unwrap(),
        Rational::zero(),
        Rational::parse(z).unwrap(),
    )
}
#[test]
fn torus_field_keeps_curved_caps_holes_and_closed_boundary_names() {
    for (start, direction, inside, outside, end) in [
        ("0", "2", "3", "9/2", "4"),
        ("5", "1", "11/2", "7", "6"),
        ("10", "-2", "7", "11/2", "6"),
    ] {
        let mut value: Value = serde_json::from_str(include_str!(
            "../../catallax-features/tests/fixtures/curved-limits/torus-entry-next.json"
        ))
        .unwrap();
        value["steps"][0]["plane"]["translation"][2] = json!(start);
        value["steps"][2]["direction"] = json!(["0", "0", direction]);
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let plan = Plan::compile(&recipe).unwrap();
        let mut workspace = Workspace::default();
        for (x, z, expected) in [
            ("3", inside, Location::Inside),
            ("3", outside, Location::Outside),
            ("0", inside, Location::Outside),
            ("4", inside, Location::Outside),
            ("3", end, Location::Boundary),
            ("3", start, Location::Boundary),
        ] {
            assert_eq!(plan.proven(&point(x, z), &mut workspace).location, expected);
        }
        for (role, p) in [
            (Role::Start, point("3", start)),
            (Role::End, point("3", end)),
        ] {
            let candidates: Vec<_> = plan
                .faces()
                .iter()
                .enumerate()
                .filter(|(_, f)| f.origin.step.as_str() == "cap" && f.origin.role == role)
                .filter_map(|(i, _)| plan.result_face(i))
                .collect();
            assert!(candidates.iter().any(
                |f| plan.face_query(*f, &p, &mut workspace).unwrap().state == FaceState::OnFace
            ));
        }
    }
}
