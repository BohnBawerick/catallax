//! Spline-face offsets wait for their own metric offset and trim proof.
use catallax_features::{BuildError, evaluate};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Surface;
use serde_json::{Value, json};

fn source() -> Value {
    serde_json::from_str(include_str!("fixtures/spline-face-offset.json")).unwrap()
}
fn recipe(source: &Value) -> Recipe {
    Recipe::parse(&source.to_string()).unwrap()
}

#[test]
fn nonzero_spline_face_offsets_name_the_m10_dependency() {
    for distance in ["1/4", "-1/4"] {
        let mut source = source();
        source["steps"][2]["distance"] = json!(distance);
        let (built, work) = with_operation_budget(1_000_000, || {
            evaluate(&recipe(&source), &mut Decider::default())
        });
        let BuildError::Geometry { step, refusal } = built.unwrap_err() else {
            panic!("expected geometry refusal");
        };
        assert_eq!(step.as_str(), "result");
        assert_eq!(refusal.code, "distance.spline-offset-unsupported");
        assert!(refusal.detail.contains("M10"));
        assert!(!work.exhausted);
    }
}

#[test]
fn zero_offset_keeps_the_represented_spline_body_identity() {
    let mut source = source();
    source["steps"][2]["distance"] = json!("0");
    let offset = evaluate(&recipe(&source), &mut Decider::default()).unwrap();
    source["steps"].as_array_mut().unwrap().pop();
    source["result"] = json!("stock");
    let stock = evaluate(&recipe(&source), &mut Decider::default()).unwrap();
    assert!(
        stock
            .model
            .faces_of_body(stock.result)
            .iter()
            .any(|face| matches!(
                stock.model.geom().surface(stock.model.face(*face).surface),
                Surface::BSpline(_)
            ))
    );
    assert_eq!(offset.result, stock.result);
    assert_eq!(
        offset.model.faces_of_body(offset.result),
        stock.model.faces_of_body(stock.result)
    );
    assert!(offset.model.check(offset.result).is_empty());
}
