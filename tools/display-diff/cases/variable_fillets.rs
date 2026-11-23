#[test]
fn display_diff_angular_zero_end_fillets() {
    for (start, end, backwards) in [(0, 2, false), (2, 0, false), (0, 2, true), (2, 0, true)] {
        let built = build(&recipe(start, end, backwards));
        let options = MeshOptions {
            chord_tolerance: q("0.01"),
            angular_deviation: Some(q("1/12")),
            ..Default::default()
        };
        let _ = catallax_numbers::with_operation_budget(20_000_000, || {
            mesh_body(&built.model, built.result, &options, &mut decider())
        });
    }
}
