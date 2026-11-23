#[test]
fn display_diff_m11_loft() {
    let plan = catallax_field::Plan::compile(&recipe("b_0", false)).unwrap();
    for chord in ["1/100", "1/200"] {
        let options = catallax_field::display::Options {
            tolerance: r(chord),
            ..Default::default()
        };
        let _ = catallax_numbers::with_operation_budget(20_000_000, || {
            catallax_field::display::mesh_plan(&plan, &options, &mut decider())
        });
    }
}
