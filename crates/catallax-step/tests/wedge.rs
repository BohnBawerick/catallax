//! A planar wedge closing to an edge too short to read leaves the file, not the model.
mod common;
use catallax_recipe::Recipe;
use common::*;

fn sketch(id: &str, points: &[[&str; 2]]) -> serde_json::Value {
    let segments: Vec<_> = points
        .iter()
        .enumerate()
        .map(|(i, to)| serde_json::json!({"id": format!("s{i}"), "kind": "line", "to": to}))
        .collect();
    serde_json::json!({
        "id": id, "op": "sketch", "loops": [segments],
        "plane": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]},
    })
}

// A 12 by 9 by 1 plate whose front edge leans `tilt` forward at x = 12, united
// with a flange 6 tall and 1 deep along that front, starting `inset` from the
// plate's left wall. SolidDraft's GUI frames gave such a bracket a 2e-17 lean
// and a 2e-28 inset.
fn bracket(tilt: &str, inset: &str) -> (Model, BodyId) {
    let lean = format!("-{tilt}");
    let recipe = serde_json::json!({
        "format": "catallax-recipe/4", "name": "leaning bracket", "result": "body",
        "steps": [
            sketch("plate-profile", &[["12", &lean], ["12", "9"], ["0", "9"], ["0", "0"]]),
            {"id": "plate", "op": "extrude", "profile": "plate-profile",
             "direction": ["0", "0", "1"], "distance": "1"},
            sketch("flange-profile", &[["12", "0"], ["12", "1"], [inset, "1"], [inset, "0"]]),
            {"id": "flange", "op": "extrude", "profile": "flange-profile",
             "direction": ["0", "0", "1"], "distance": "6"},
            {"id": "body", "op": "boolean", "kind": "union", "target": "plate", "tool": "flange"},
        ],
    });
    built(recipe)
}

// A plate whose left wall leans `lean` out at its back, under a flange along
// that wall, all multiplied through by `scale`: the plate's top keeps a
// three-sided wedge beside the flange.
fn wedge_face(lean: &str, scale: &str) -> (Model, BodyId) {
    let k = r(scale);
    let at = |v: &str| (&r(v) * &k).to_string();
    let back = at(&format!("-{lean}"));
    let [zero, one, six, nine, twelve] = ["0", "1", "6", "9", "12"].map(at);
    let recipe = serde_json::json!({
        "format": "catallax-recipe/4", "name": "leaning plate", "result": "body",
        "steps": [
            sketch("plate-profile", &[[&twelve, &zero], [&twelve, &nine], [&back, &nine], [&zero, &zero]]),
            {"id": "plate", "op": "extrude", "profile": "plate-profile",
             "direction": ["0", "0", "1"], "distance": one},
            sketch("flange-profile", &[[&one, &zero], [&one, &nine], [&zero, &nine], [&zero, &zero]]),
            {"id": "flange", "op": "extrude", "profile": "flange-profile",
             "direction": ["0", "0", "1"], "distance": six},
            {"id": "body", "op": "boolean", "kind": "union", "target": "plate", "tool": "flange"},
        ],
    });
    built(recipe)
}

fn built(recipe: serde_json::Value) -> (Model, BodyId) {
    let recipe = Recipe::parse(&recipe.to_string()).unwrap();
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    (built.model, built.result)
}

fn written(model: &Model, body: BodyId) -> (Step, StepReport) {
    let mut bytes = Vec::new();
    let answer = write_step(
        model,
        &[body],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    (Step::parse(&bytes, "catallax"), answer.value)
}

fn assert_reconciled(model: &Model, body: BodyId, report: &StepReport) {
    let exact = model.counts(body);
    let faces: usize = report.regularizations.iter().map(|x| x.removed_faces).sum();
    let edges: usize = report.regularizations.iter().map(|x| x.removed_edges).sum();
    assert_eq!(exact.faces, report.faces + faces);
    assert_eq!(exact.edges + report.added_edges, report.edges + edges);
}

#[test]
fn a_leaning_front_and_an_inset_flange_are_written_as_flush_walls() {
    for (tilt, inset) in [("2e-17", "2e-28"), ("1e-7", "1e-30"), ("9e-8", "9e-15")] {
        let (model, body) = bracket(tilt, inset);
        let exact = model.counts(body);
        assert_eq!(exact.faces, 10, "{tilt} {inset}");
        let (file, report) = written(&model, body);
        // The L-shaped prism the flush bracket would be.
        assert_eq!(report.faces, 8, "{tilt} {inset}");
        assert_eq!(file.count("ADVANCED_FACE"), 8);
        assert_eq!(file.count("EDGE_CURVE"), report.edges);
        let [front, side] = report.regularizations.as_slice() else {
            panic!("{tilt} {inset}: {:?}", report.regularizations);
        };
        // Each closes one end and joins the far ends with one new edge.
        for record in [front, side] {
            assert_eq!((record.removed_faces, record.removed_edges), (1, 2));
            assert!(report.max_deviation >= record.max_deviation);
        }
        assert_eq!(front.max_deviation, r(tilt));
        assert_eq!(side.max_deviation, r(inset));
        assert_reconciled(&model, body, &report);
        // The exact model is untouched.
        assert_eq!(model.counts(body), exact);
        assert!(model.check(body).is_empty());
    }
}

#[test]
fn a_three_sided_wedge_leaves_with_its_wall() {
    for lean in ["1e-30", "2e-17", "1e-7"] {
        let (model, body) = wedge_face(lean, "1");
        let exact = model.counts(body);
        let (file, report) = written(&model, body);
        let [record] = report.regularizations.as_slice() else {
            panic!("{lean}: {:?}", report.regularizations);
        };
        // The wedge, its moved wall and its three edges.
        assert_eq!((record.removed_faces, record.removed_edges), (2, 3));
        assert_eq!(record.max_deviation, r(lean));
        assert_eq!(report.faces + 2, exact.faces);
        assert_eq!(file.count("ADVANCED_FACE"), report.faces);
        assert_reconciled(&model, body, &report);
        assert_eq!(model.counts(body), exact);
    }
}

#[test]
fn a_wedge_wider_than_the_floor_keeps_its_bytes() {
    for (model, body) in [
        bracket("1.00000000001e-7", "0"),
        wedge_face("1.00000000001e-7", "1"),
    ] {
        let (_, report) = written(&model, body);
        assert!(report.regularizations.is_empty());
        assert_eq!(report.faces, model.counts(body).faces);
    }
}

#[test]
fn a_uniformly_tiny_wedge_is_left_alone() {
    // Scaled with its part, a wedge stays wide against its length or every
    // edge beside it falls below the floor.
    for scale in ["1e-30", "1e-20", "1e-8"] {
        for lean in ["1", "1/2", "1e-3"] {
            let (model, body) = wedge_face(lean, scale);
            let (_, report) = written(&model, body);
            assert!(report.regularizations.is_empty(), "{scale} {lean}");
        }
    }
}
