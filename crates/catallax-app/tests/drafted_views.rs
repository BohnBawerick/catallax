//! Drawing requests for the drafted SolidDraft scenarios.

use catallax_app::DEFAULT_BUILD_WORK;
use catallax_app::views::{VIEW_RAY_MARK, VIEW_SPLIT_MARK, hidden_line_view, stroke_drawing};
use catallax_features::evaluate;
use catallax_numbers::{Decider, Refusal, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::RVector3;
use catallax_topology::{BodyId, Model};

const DRAWING_WORK_TARGET: u64 = 20_000_000;

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

fn part(source: &str) -> (Model, BodyId) {
    let recipe = Recipe::parse(source).unwrap();
    let (built, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        evaluate(&recipe, &mut Decider::default())
    });
    assert!(!work.exhausted);
    let built = built.unwrap();
    (built.model, built.result)
}

#[test]
fn drafted_parts_fit_the_drawing_request() {
    use serde_json::{Value, json};

    // The preserved-bore fixture matches SolidDraft's draft-countersink scenario.
    let countersink: Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/draft-preserved-bore/P-N.json"
    ))
    .unwrap();
    let mut through = countersink.clone();
    through["steps"][2]["at"][0]["point"] = json!(["0", "0", "20"]);
    through["steps"][2]["diameter"] = json!("10");
    through["steps"][2].as_object_mut().unwrap().remove("head");
    let mut plain = countersink.clone();
    plain["steps"].as_array_mut().unwrap().remove(2);
    plain["steps"][2]["body"] = json!("base");
    let mut boss = plain.clone();
    boss["steps"].as_array_mut().unwrap().extend([
        json!({
            "id": "boss-profile", "op": "sketch",
            "plane": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "20"]},
            "loops": [[{"id": "wall", "kind": "circle", "centre": ["0", "0"], "radius": "5"}]]
        }),
        json!({"id": "boss", "op": "extrude", "profile": "boss-profile", "distance": "5"}),
        json!({"id": "joined", "op": "boolean", "kind": "union", "target": "drafted", "tool": "boss"}),
    ]);
    boss["result"] = json!("joined");

    let mut failures = Vec::new();
    for (name, mut recipe, draft_index, angles) in [
        ("draft", plain, 2, &["3"][..]),
        ("draft-with-boss", boss, 2, &["5"][..]),
        ("draft-countersink", countersink, 3, &["3", "5"][..]),
        ("draft-through-hole", through, 3, &["3", "5"][..]),
    ] {
        for angle in angles {
            recipe["steps"][draft_index]["angle"] = json!(angle);
            let (model, body) = part(&recipe.to_string());
            for (orientation, direction) in [
                ("front", vector(0, -1, 0)),
                ("top", vector(0, 0, 1)),
                ("right", vector(1, 0, 0)),
                ("iso", vector(1, -1, 1)),
            ] {
                let (answer, work) =
                    with_operation_budget(DEFAULT_BUILD_WORK, || -> Result<_, Refusal> {
                        let mut dec = Decider::default();
                        let view = hidden_line_view(&model, body, &direction, &mut dec)
                            .map_err(|q| q.refusal)?;
                        assert!(
                            !view
                                .marks
                                .iter()
                                .any(|m| { matches!(m.code, VIEW_RAY_MARK | VIEW_SPLIT_MARK) })
                        );
                        let mut drawn = 0;
                        for stroke in &view.value.strokes {
                            if !stroke.seam && !stroke.smooth {
                                let curve = stroke_drawing(&view.value, stroke, &mut dec)?;
                                assert!(curve.value.fit.is_none());
                                drawn += 1;
                            }
                        }
                        assert!(drawn > 0);
                        Ok(drawn)
                    });
                eprintln!(
                    "{name} {angle} {orientation}: {answer:?} at {} work",
                    work.used
                );
                if answer.is_err() || work.exhausted || work.used > DRAWING_WORK_TARGET {
                    failures.push((name, angle, orientation, answer, work));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
