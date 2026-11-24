#[test]
fn display_diff_equal_crossing_within_target() {
    let built = evaluate(&parsed(&recipe("union", "1", "8", "0")), &mut decider()).unwrap();
    let options = catallax_mesh::MeshOptions {
        chord_tolerance: Rational::from_ratio(6, 125).unwrap(),
        angular_deviation: Some(Rational::from_ratio(314159265, 1200000000).unwrap()),
        ..Default::default()
    };
    let _ = catallax_numbers::with_operation_budget(20_000_000, || {
        catallax_mesh::mesh_body(&built.model, built.result, &options, &mut decider())
    });
}
