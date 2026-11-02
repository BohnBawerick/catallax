use catallax_features::evaluate;
use catallax_numbers::{Ball, Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn document() -> Value {
    serde_json::from_str(include_str!("fixtures/nonconvex-shell.json")).unwrap()
}
fn measure(source: Value, constant: &str, pi: &str) {
    let recipe = Recipe::parse(&source.to_string()).unwrap();
    let (answer, work) = with_operation_budget(100_000_000, || {
        evaluate(&recipe, &mut Decider::default())
    });
    eprintln!("MULTIWAY result_work={work:?} result={:?}", answer.as_ref().map(|_| "ok"));
    assert!(!work.exhausted, "{work:?}");
    let built = answer.unwrap();
    assert!(built.model.check(built.result).is_empty());
    let (mass, mass_work) = with_operation_budget(100_000_000, || {
        catallax_measure::mass_properties(
            &built.model,
            built.result,
            &Default::default(),
            &mut Decider::default(),
        )
    });
    eprintln!("MULTIWAY phase=Measure work={mass_work:?} result={:?}", mass.as_ref().map(|_| "ok"));
    assert!(!mass_work.exhausted, "{mass_work:?}");
    let volume = mass.unwrap().value.volume;
    let ball = |s| Ball::from_rational(&Rational::parse(s).unwrap(), 512);
    let expected = ball(constant) + Ball::pi(512) * ball(pi);
    assert!(volume.lo <= expected.lo());
    assert!(volume.hi >= expected.hi());
    assert!(volume.width() < Rational::parse("1e-6").unwrap());
}
fn chain(points: &[[i64; 2]], prefix: &str) -> Value {
    json!(points.iter().enumerate().map(|(i, _)| json!({
        "kind":"line", "id":format!("{prefix}-{i}"),
        "to":points[(i+1)%points.len()].map(|n| n.to_string())
    })).collect::<Vec<_>>())
}
#[test]
fn measure_reentrant_tray() {
    measure(document(), "120", "3/4");
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
    measure(source, "300", "9");
}
