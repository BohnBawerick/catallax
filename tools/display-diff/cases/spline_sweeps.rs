#[test]
fn display_diff_multichart_preview() {
    for reverse in [false, true] {
        let mut value = definition(true);
        let path = &mut value["steps"][1]["segments"][0];
        path["degree"] = json!(2);
        let heights = if reverse { ["-1", "-2", "-3"] } else { ["1", "2", "3"] };
        path["control_points"] = json!([
            ["0", "0", "0"],
            ["1", "0", heights[0]],
            ["1", "1", heights[1]],
            ["0", "1", heights[2]]
        ]);
        path["weights"] = json!(["1", "2", "2", "1"]);
        path["knots"] = json!(["0", "0", "0", "1/2", "1", "1", "1"]);
        let built = catallax_features::evaluate(&recipe(&value), &mut Decider::default()).unwrap();
        let _ = catallax_numbers::with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            catallax_mesh::preview_mesh_body(
                &built.model,
                built.result,
                &catallax_mesh::PreviewOptions::default(),
                &mut Decider::default(),
            )
        });
    }
}
