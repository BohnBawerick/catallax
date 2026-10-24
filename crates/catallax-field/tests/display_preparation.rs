//! Unsupported charts refuse before preparing restrictions for earlier supports.
use catallax_field::{
    Plan,
    display::{Diagnostics, Options, mesh_plan_reported},
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::json;

#[test]
fn unsupported_slanted_cap_does_not_prepare_cylinder_restrictions() {
    let recipe = Recipe::parse(
        &json!({
            "format":"catallax-recipe/3", "name":"slanted cap display preparation",
            "steps":[
                {"id":"base", "op":"cylinder", "radius":"1", "height":"4",
                 "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
                {"id":"drafted", "op":"draft", "body":"base",
                 "faces":[{"face":{"step":"base","role":"end"}}],
                 "pull":["-1","0","1"],
                 "neutral":{"point":["0","0","0"],"normal":["1","0","0"]},
                 "angle":"0"}
            ], "result":"drafted"
        })
        .to_string(),
    )
    .unwrap();
    let (result, work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
        let mut report = Diagnostics::default();
        let result = mesh_plan_reported(&plan, &Options::default(), &mut dec, &mut report);
        (result, report)
    });
    assert!(!work.exhausted);
    let (result, report) = result;
    assert_eq!(result.unwrap_err().code, "field.display-unsupported");
    assert_eq!(report.total.cells, 0);
    assert_eq!(report.total.work, 0);
    assert!(report.supports.iter().all(|s| s.total.cells == 0));
}
