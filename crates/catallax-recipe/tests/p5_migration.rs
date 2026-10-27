//! The reviewed housing migration changes overflow meaning, not its authored graph.

use catallax_recipe::{FORMAT, FORMAT_V4, Number, Operation, Recipe};
use serde_json::Value;

#[test]
fn p5_v4_preserves_the_v3_graph_names_and_authored_radius_edit() {
    let mut legacy = Recipe::parse(include_str!("fixtures/P5-flanged-housing-v3.json")).unwrap();
    let mut migrated =
        Recipe::parse(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    assert_eq!(legacy.format, FORMAT);
    assert_eq!(migrated.format, FORMAT_V4);
    assert_eq!(migrated.steps, legacy.steps);
    assert_eq!(migrated.result, legacy.result);
    assert_eq!(migrated.name, legacy.name);

    let case: Value = serde_json::from_str(include_str!(
        "../../../.claude/skills/verify-catallax/scripts/head-to-head/cases/P5-flanged-housing.json"
    ))
    .unwrap();
    let edit = &case["edit"];
    assert_eq!(edit["step"], "filleted");
    assert_eq!(edit["path"], serde_json::json!(["radius"]));
    let edited_radius = Number::parse(edit["value"].as_str().unwrap()).unwrap();
    assert_eq!(edited_radius, Number::parse("2.5").unwrap());

    for recipe in [&mut legacy, &mut migrated] {
        let step = recipe
            .steps
            .iter_mut()
            .find(|step| step.id.as_str() == edit["step"].as_str().unwrap())
            .unwrap();
        let Operation::Fillet { radius, .. } = &mut step.op else {
            panic!("the authored edit must still address the junction fillet");
        };
        assert_eq!(*radius, Number::from(2));
        *radius = edited_radius.clone();
        recipe.check().unwrap();
        assert_eq!(Recipe::parse(&recipe.to_text()).unwrap(), *recipe);
    }
    assert_eq!(migrated.steps, legacy.steps);
    assert_eq!(migrated.result, legacy.result);
}
