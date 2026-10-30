//! SolidDraft's on-face slot cut through the bracket flange, with no start offset.
//! Each call rebuilds only the slot from its target's checkpoint, as the application
//! does, and must stay within the application allowance.
use catallax_app::{DEFAULT_BUILD_WORK, build_incremental_without_display};
use catallax_measure::{Accuracy, volume};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;

fn slot_cut(source: &str, faces: usize, expected_volume: &str) {
    let recipe = Recipe::parse(source).unwrap();
    let slot = recipe.steps.last().unwrap();
    let catallax_recipe::Operation::Boolean { target, .. } = &slot.op else {
        panic!("the last step is the slot subtract");
    };
    // The target body, built as an earlier application call.
    let end = recipe.steps.iter().position(|s| s.id == *target).unwrap();
    let mut prefix = recipe.clone();
    prefix.steps.truncate(end + 1);
    prefix.result = target.clone();
    let (target, work) = with_operation_budget(100_000_000, || {
        build_incremental_without_display(&prefix, None, &mut Decider::default())
    });
    assert!(!work.exhausted);
    let target = target.unwrap();

    let (body, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        build_incremental_without_display(&recipe, Some(target.checkpoint), &mut Decider::default())
    });
    eprintln!("{}: slot subtract {work:?}", recipe.name);
    assert!(
        !work.exhausted,
        "slot subtract over the allowance: {work:?}"
    );
    let body = body.unwrap();
    let (display, work) = with_operation_budget(DEFAULT_BUILD_WORK, || {
        body.display(&MeshOptions::default(), &mut Decider::default())
    });
    eprintln!("{}: display {work:?}", recipe.name);
    assert!(!work.exhausted, "display over the allowance: {work:?}");
    display.unwrap();

    let built = body.checkpoint.built();
    assert_eq!(built.model.faces_of_body(built.result).len(), faces);
    let v = volume(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    let expected = Rational::parse(expected_volume).unwrap();
    let tolerance = Rational::parse("1e-6").unwrap();
    assert!(
        v.lo > &expected - &tolerance && v.hi < &expected + &tolerance,
        "volume [{}, {}]",
        v.lo.to_f64_nearest(),
        v.hi.to_f64_nearest()
    );
}

#[test]
fn on_face_slot_through_the_flange_fits_the_application_allowance() {
    slot_cut(
        include_str!("fixtures/on-face-slot-minimal.json"),
        17,
        "174366.367281174",
    );
}

#[test]
fn on_face_slot_through_the_bored_bracket_fits_the_application_allowance() {
    slot_cut(
        include_str!("fixtures/on-face-slot-bracket.json"),
        18,
        "170521.057873181",
    );
}
