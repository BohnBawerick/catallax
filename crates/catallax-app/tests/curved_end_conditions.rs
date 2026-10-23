//! Real application meshes retain the curved end and all authored side names.
use catallax_numbers::Decider;
use catallax_recipe::{Recipe, Role};

#[test]
fn application_builds_and_picks_a_curved_first_hit_cap() {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/curved-limits/sphere-entry-to.json"
    ))
    .unwrap();
    let body = catallax_app::build(&recipe, &Default::default(), &mut Decider::default()).unwrap();
    assert!(body.built.model.check(body.built.result).is_empty());
    assert!(body.display.marks.is_empty());
    assert!(!body.display.value.triangles.is_empty());
    for role in [Role::Start, Role::End] {
        assert!(body.display.value.triangles.iter().any(|t| {
            body.names.faces[&t.face]
                .iter()
                .any(|o| o.step.as_str() == "cap" && o.role == role)
        }));
    }
    assert!(body.names.edges.values().all(Option::is_some));
    assert!(
        body.names
            .faces
            .values()
            .flatten()
            .all(|o| o.step.as_str() == "cap")
    );
}
