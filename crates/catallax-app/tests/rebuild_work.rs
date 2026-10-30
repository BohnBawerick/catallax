//! Certified display is requested separately from exact suffix rebuilding.
use catallax_app::{build_incremental, build_incremental_without_display};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

fn bracket_calls() -> Vec<Recipe> {
    let question: Value =
        serde_json::from_str(include_str!("fixtures/bracket-width-140.question.json")).unwrap();
    let recipe = &question["problem"];
    let steps = recipe["steps"].as_array().unwrap();
    steps
        .iter()
        .enumerate()
        .filter_map(|(i, step)| {
            let selected = match step["op"].as_str().unwrap() {
                "extrude" => vec![steps[i - 1].clone(), step.clone()],
                "boolean" => steps[..=i].to_vec(),
                _ => return None,
            };
            let mut call = recipe.clone();
            call["result"] = step["id"].clone();
            call["steps"] = selected.into();
            Some(Recipe::parse(&call.to_string()).unwrap())
        })
        .collect()
}

#[test]
fn bracket_suffix_builds_do_not_spend_unrequested_display_work() {
    let mut done = Vec::<(Recipe, catallax_features::EvaluationCheckpoint)>::new();
    let mut final_body = None;
    let ((), work) = with_operation_budget(100_000_000, || {
        for recipe in bracket_calls() {
            let prefix = done
                .iter()
                .filter(|(prior, _)| recipe.steps.starts_with(&prior.steps))
                .max_by_key(|(prior, _)| prior.steps.len())
                .map(|(_, checkpoint)| checkpoint.clone());
            let body = build_incremental_without_display(&recipe, prefix, &mut Decider::default())
                .unwrap();
            done.push((recipe, body.checkpoint.clone()));
            final_body = Some(body);
        }
    });
    eprintln!("17 build calls without display: {work:?}");
    assert!(!work.exhausted);
    assert!(
        work.used <= 796_422,
        "unrequested certified display work: {work:?}"
    );
    let body = final_body.unwrap();
    let (display, work) = with_operation_budget(100_000_000, || {
        body.display(&MeshOptions::default(), &mut Decider::default())
            .unwrap()
    });
    eprintln!("requested final display: {work:?}");
    assert!(!work.exhausted);
    assert!(
        work.used <= 2_775_075,
        "requested certified display work: {work:?}"
    );
    assert_display_identity(&body, &display, &MeshOptions::default());
    assert_eq!(
        body.names,
        build_incremental(
            &done.last().unwrap().0,
            None,
            &MeshOptions::default(),
            &mut Decider::default()
        )
        .unwrap()
        .names
    );
}

fn assert_display_identity(
    body: &catallax_app::EvaluatedApplicationBody,
    display: &catallax_numbers::Answer<catallax_mesh::Mesh>,
    options: &MeshOptions,
) {
    let built = body.checkpoint.built();
    let expected = mesh_body(&built.model, built.result, options, &mut Decider::default()).unwrap();
    assert_eq!(display.value, expected.value);
    assert_eq!(display.marks, expected.marks);
    let mut actual_obj = Vec::new();
    let mut expected_obj = Vec::new();
    display.value.write_obj(&mut actual_obj).unwrap();
    expected.value.write_obj(&mut expected_obj).unwrap();
    assert_eq!(actual_obj, expected_obj);
}

#[test]
fn deferred_display_and_names_match_eager_fixture_bodies() {
    for source in [
        include_str!("fixtures/bracket.json"),
        include_str!("fixtures/rounded-base-flange.json"),
        include_str!("../../catallax-features/tests/fixtures/algebraic-extrusions/rectangle.json"),
        include_str!(
            "../../catallax-features/tests/fixtures/circular-sphere-limits/sphere-entry-next.json"
        ),
        include_str!("../../catallax-features/tests/fixtures/curved-limits/sphere-entry-to.json"),
        include_str!("../../catallax-features/tests/fixtures/open-spherical-cap/base-open.json"),
        include_str!("../../../parts/P3-filleted-boss.json"),
    ] {
        let value: Value = serde_json::from_str(source).unwrap();
        let recipe = if value["format"] == "soliddraft-recipe/1" {
            catallax_app::soliddraft::convert(source).unwrap().recipe
        } else {
            Recipe::parse(source).unwrap()
        };
        let options = MeshOptions::default();
        let eager = build_incremental(&recipe, None, &options, &mut Decider::default()).unwrap();
        let deferred =
            build_incremental_without_display(&recipe, None, &mut Decider::default()).unwrap();
        assert_eq!(deferred.names, eager.names);
        let display = deferred.display(&options, &mut Decider::default()).unwrap();
        assert_eq!(display.value, eager.display.value);
        assert_eq!(display.marks, eager.display.marks);
        assert_display_identity(&deferred, &display, &options);
    }
}

fn cylinder(radius: i64) -> Recipe {
    catallax_app::soliddraft::convert(&json!({
        "format":"soliddraft-recipe/1", "name":"display-refusal", "result":"body", "steps":[
            {"id":"profile","op":"profile","plane":"XY","circle":{"center":[0,0],"radius":radius}},
            {"id":"body","op":"extrude","profile":"profile","distance":10}
        ]
    }).to_string()).unwrap().recipe
}

#[test]
fn deferred_display_preserves_refusal_codes_and_caller_allowances() {
    let recipe = cylinder(5);
    let options = MeshOptions {
        max_vertices: 8,
        ..Default::default()
    };
    let eager = build_incremental(&recipe, None, &options, &mut Decider::default()).unwrap_err();
    assert!(eager.code().starts_with("mesh."), "{eager:?}");
    let deferred =
        build_incremental_without_display(&recipe, None, &mut Decider::default()).unwrap();
    assert_eq!(
        deferred
            .display(&options, &mut Decider::default())
            .unwrap_err()
            .code(),
        eager.code()
    );
    let (answer, work) = with_operation_budget(0, || {
        deferred.display(&MeshOptions::default(), &mut Decider::default())
    });
    assert_eq!(answer.unwrap_err().code(), "numbers.operation-budget");
    assert!(work.exhausted);
    assert_eq!(work.used, 0);
    let (answer, work) = with_operation_budget(0, || {
        build_incremental_without_display(&recipe, None, &mut Decider::default())
    });
    assert_eq!(answer.unwrap_err().code(), "numbers.operation-budget");
    assert!(work.exhausted);
    assert_eq!(work.used, 0);
    let display = deferred
        .display(&MeshOptions::default(), &mut Decider::default())
        .unwrap();
    assert_display_identity(&deferred, &display, &MeshOptions::default());
}

#[test]
fn deferred_display_keeps_the_default_application_guard() {
    let recipe = cylinder(500_000);
    let options = MeshOptions::default();
    let eager = build_incremental(&recipe, None, &options, &mut Decider::default()).unwrap_err();
    assert_eq!(eager.code(), "application.work-limit");
    let deferred =
        build_incremental_without_display(&recipe, None, &mut Decider::default()).unwrap();
    assert_eq!(
        deferred
            .display(&options, &mut Decider::default())
            .unwrap_err()
            .code(),
        eager.code()
    );
}
