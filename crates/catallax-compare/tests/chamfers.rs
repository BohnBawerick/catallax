//! Certified trimmed distances must resolve the real default comparison.

use catallax_compare::{Completion, Controller, Verdict};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;
use serde_json::json;

fn agrees(recipe: &Recipe) {
    let controller = Controller::default();
    let revision = controller.begin();
    let (completion, work) =
        with_operation_budget(100_000_000, || revision.check(recipe, &Default::default()));
    assert_eq!(completion, Completion::Published);
    assert!(!work.exhausted);
    let report = controller.latest().unwrap();
    assert_eq!(report.verdict, Verdict::Agree, "{:?}", report.findings);
    assert!(report.findings.is_empty());
    assert!(report.not_comparable.is_empty());
    assert!(report.forward_vertices > 0 && report.reverse_samples > 0);
    assert!(
        report
            .sampled_distance_upper
            .as_ref()
            .is_some_and(|upper| { upper <= &catallax_compare::Options::default().accuracy })
    );
}

#[test]
fn p2_compares_in_both_directions_with_the_default_allowance() {
    agrees(&Recipe::parse(include_str!("../../../parts/P2-collar.json")).unwrap());
}

#[test]
fn fractional_setbacks_include_the_closed_circular_contacts() {
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("../../../parts/P2-collar.json")).unwrap();
    value["steps"][2]["distance"] = json!("3/4");
    agrees(&Recipe::parse(&value.to_string()).unwrap());
}

#[test]
fn unequal_inner_rim_excludes_the_consumed_sharp_circle() {
    let pose = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
    let value = json!({"format":"catallax-recipe/3","name":"unequal inner rim","steps":[
        {"id":"outer","op":"cylinder","radius":"8","height":"9","pose":pose},
        {"id":"inner","op":"cylinder","radius":"3","height":"9","pose":pose},
        {"id":"ring","op":"boolean","kind":"subtract","target":"outer","tool":"inner"},
        {"id":"bevel","op":"chamfer","body":"ring","distance":"2","other_distance":"1","edges":[{"between":[
            {"step":"outer","role":"start"}, {"step":"inner","role":{"side":{"segment":"wall"}}}
        ]}]}
    ],"result":"bevel"});
    agrees(&Recipe::parse(&value.to_string()).unwrap());
}
