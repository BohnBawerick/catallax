#[test]
fn display_diff_short_offset_prism() {
    let certified = build(&recipe("-1/4", "5001/10000"));
    let built = certified.built();
    let _ = catallax_mesh::mesh_body(
        &built.model,
        built.result,
        &Default::default(),
        &mut Decider::default(),
    );
}
