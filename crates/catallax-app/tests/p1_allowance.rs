//! Canonical P1 fits the eager application's unchanged request allowance.
use catallax_app::{DEFAULT_BUILD_WORK, build};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Rational, set_threads, with_operation_budget};
use catallax_recipe::Recipe;
use std::num::NonZeroUsize;

#[test]
fn canonical_p1_builds_with_display_inside_the_default_allowance_at_every_thread_count() {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let options = MeshOptions {
        chord_tolerance: Rational::parse("0.01").unwrap(),
        max_vertices: 200_000,
        ..Default::default()
    };
    assert_eq!(DEFAULT_BUILD_WORK, 20_000_000);
    let mut reference = None;
    for threads in [1, 2, 4] {
        set_threads(NonZeroUsize::new(threads));
        let (body, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
            build(&recipe, &options, &mut Decider::default()).unwrap()
        });
        assert!(!work.exhausted);
        assert!(work.used <= DEFAULT_BUILD_WORK);
        assert!(body.built.model.check(body.built.result).is_empty());
        assert!(body.display.marks.is_empty());
        assert_eq!(body.names.faces.len(), 19);
        assert_eq!(body.names.edges.len(), 42);
        assert_eq!(body.names.vertices.len(), 20);
        assert_eq!(body.display.value.positions.len(), 2800);
        assert_eq!(body.display.value.triangles.len(), 5620);
        let snapshot = (body.names, body.display.value, work.used);
        if let Some(reference) = &reference {
            assert_eq!(&snapshot, reference, "{threads} threads");
        } else {
            reference = Some(snapshot);
        }
        println!("{threads} threads: {} work units", work.used);
    }
    set_threads(None);
}
