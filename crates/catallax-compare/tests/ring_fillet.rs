//! The field does not build the ring between nonparallel cylinders yet, so the
//! comparison has one side and must say so instead of agreeing.
use catallax_compare::{Completion, Controller, Options, Verdict};
use catallax_recipe::Recipe;
use serde_json::json;

#[test]
fn a_ring_with_one_side_stays_unresolved() {
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"tee","steps":[
            {"id":"wide","op":"cylinder","radius":"10","height":"60",
                "pose":{"rotation":["1","1","0","0"],"translation":["0","30","0"]}},
            {"id":"wall","op":"cylinder","radius":"4","height":"30",
                "pose":{"rotation":["1","0","0","0"],"translation":["1","0","0"]}},
            {"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"},
            {"id":"round","op":"fillet","body":"joined","radius":"1",
                "edges":[{"between":[{"step":"wide","role":{"side":{"segment":"wall"}}},
                    {"step":"wall","role":{"side":{"segment":"wall"}}}]}]}
        ],"result":"round"})
        .to_string(),
    )
    .unwrap();
    let controller = Controller::default();
    assert_eq!(
        controller.begin().check(&recipe, &Options::default()),
        Completion::Published
    );
    let report = controller.latest().unwrap();
    assert_eq!(report.verdict, Verdict::Unresolved, "{report:?}");
    assert!(!report.findings.is_empty());
    assert_eq!(report.forward_vertices, 0);
    assert_eq!(report.reverse_samples, 0);
}
