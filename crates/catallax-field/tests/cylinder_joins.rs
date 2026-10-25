//! Native field display remains an explicit gate at unchanged public controls.
use catallax_field::{
    Plan,
    display::{Options, mesh_plan},
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
#[test]
fn bounded_cylinder_join_native_display_has_complete_certified_outputs() {
    let mut failures = Vec::new();
    for text in [
        include_str!("../../catallax-features/tests/fixtures/cylinder-joins/boss-plate.json"),
        include_str!("../../catallax-features/tests/fixtures/cylinder-joins/stepped-shaft.json"),
        include_str!(
            "../../catallax-features/tests/fixtures/cylinder-joins/boss-plate-thicken.json"
        ),
        include_str!(
            "../../catallax-features/tests/fixtures/cylinder-joins/stepped-shaft-thicken.json"
        ),
    ] {
        let r = Recipe::parse(text).unwrap();
        let options = Options::default();
        let started = std::time::Instant::now();
        let (answer, work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let plan = Plan::compile_checked(&r, &mut dec)?;
            mesh_plan(&plan, &options, &mut dec)
        });
        eprintln!(
            "{}: {} work, {:?} wall; {:?}",
            r.name,
            work.used,
            started.elapsed(),
            answer
                .as_ref()
                .map(|m| (m.positions.len(), m.triangles.len(), m.cells))
        );
        match answer {
            Ok(mesh) => {
                assert!(!mesh.positions.is_empty() && !mesh.triangles.is_empty());
                assert!(mesh.max_deviation <= options.tolerance);
                assert!(!work.exhausted);
            }
            Err(e) => failures.push((r.name, e)),
        }
    }
    assert!(failures.is_empty(), "{failures:?}");
}
