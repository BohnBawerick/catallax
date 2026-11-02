use catallax_features::evaluate;
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn document() -> Value {
    serde_json::from_str(include_str!("fixtures/nonconvex-shell.json")).unwrap()
}
fn measure(source: Value) {
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let (answer, work) = with_operation_budget(100_000_000, || {
        evaluate(&recipe, &mut Decider::default())
    });
    eprintln!("MULTIWAY result_work={work:?} result={:?}", answer.as_ref().map(|_| "ok"));
    assert!(!work.exhausted, "{work:?}");
    let built = answer.unwrap();
    assert!(built.model.check(built.result).is_empty());
}
fn chain(points: &[[i64; 2]], prefix: &str) -> Value {
    json!(points.iter().enumerate().map(|(i, _)| json!({
        "kind":"line", "id":format!("{prefix}-{i}"),
        "to":points[(i+1)%points.len()].map(|n| n.to_string())
    })).collect::<Vec<_>>())
}
#[test]
fn measure_reentrant_tray() {
    measure(document());
}
#[test]
fn measure_holed_shell() {
    let mut source = document();
    source["steps"][0]["loops"] = json!([
        chain(&[[0,0],[10,0],[10,10],[0,10]], "outer"),
        chain(&[[3,3],[3,7],[7,7],[7,3]], "hole")
    ]);
    source["steps"][2]["open"] = json!([
        {"face":{"step":"stock","role":"start"}},
        {"face":{"step":"stock","role":"end"}}
    ]);
    source["steps"][2]["thickness"] = json!("3/2");
    measure(source);
}
