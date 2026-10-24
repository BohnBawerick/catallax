//! Prefix reuse must keep full-rebuild topology, origins and STEP bytes.
use catallax_features::{Built, EvaluationCheckpoint, evaluate, evaluate_incremental};
use catallax_numbers::{BudgetPolicy, Decider, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use serde_json::json;

fn recipe() -> Recipe {
    Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"prefix",
        "steps":[
            {"id":"stock","op":"block","size":["8","6","4"],
                "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"tool","op":"cylinder","radius":"1","height":"6",
                "pose":{"rotation":["1","0","0","0"],"translation":["3","3","-1"]}},
            {"id":"cut","op":"boolean","kind":"subtract","target":"stock","tool":"tool"}
        ],"result":"cut"})
        .to_string(),
    )
    .unwrap()
}

fn export(built: &Built) -> Vec<u8> {
    assert!(built.model.check(built.result).is_empty());
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    bytes
}

fn prefix(recipe: &Recipe, count: usize) -> EvaluationCheckpoint {
    let mut prefix = recipe.clone();
    prefix.steps.truncate(count);
    prefix.result = prefix.steps[count - 1].id.clone();
    evaluate_incremental(&prefix, None, &mut Decider::default()).unwrap()
}

fn assert_same(recipe: &Recipe, checkpoint: EvaluationCheckpoint) {
    let full = evaluate(recipe, &mut Decider::default()).unwrap();
    let resumed = evaluate_incremental(recipe, Some(checkpoint), &mut Decider::default()).unwrap();
    let resumed = resumed.built();
    assert_eq!(export(resumed), export(&full));
    assert_eq!(resumed.bodies, full.bodies);
    assert_eq!(
        resumed.model.counts(resumed.result),
        full.model.counts(full.result)
    );
    for (a, b) in resumed
        .model
        .faces_of_body(resumed.result)
        .into_iter()
        .zip(full.model.faces_of_body(full.result))
    {
        assert_eq!(resumed.model.face_origins(a), full.model.face_origins(b));
    }
    assert_eq!(resumed.recipe(), recipe);
}

#[test]
fn append_result_change_and_earlier_edits_match_fresh_builds() {
    let recipe = recipe();
    assert_same(&recipe, prefix(&recipe, 2));
    let checkpoint = prefix(&recipe, 2);
    let mut edited = recipe.clone();
    if let catallax_recipe::Operation::Block { size, .. } = &mut edited.steps[0].op {
        size[0] = catallax_recipe::Number::new(10.into());
    } else {
        panic!("stock block");
    }
    assert_same(&edited, checkpoint.clone());
    edited = recipe.clone();
    edited.steps.swap(0, 1);
    assert_same(&edited, checkpoint.clone());
    edited = recipe.clone();
    edited.steps.truncate(1);
    edited.result = StepId::new("stock").unwrap();
    assert_same(&edited, checkpoint.clone());
    edited = recipe.clone();
    edited.name = "renamed".into();
    edited.description = "A new description does not change geometry".into();
    edited.result = StepId::new("stock").unwrap();
    assert_same(&edited, checkpoint);
}

#[test]
fn a_different_numeric_policy_cannot_reuse_an_earlier_proof() {
    let recipe = recipe();
    let checkpoint = evaluate_incremental(&recipe, None, &mut Decider::default()).unwrap();
    let policy = BudgetPolicy {
        max_work: 0,
        ..Default::default()
    };
    let full = evaluate(&recipe, &mut Decider::new(policy)).unwrap_err();
    let resumed =
        evaluate_incremental(&recipe, Some(checkpoint), &mut Decider::new(policy)).unwrap_err();
    assert_eq!(full.to_string(), resumed.to_string());
}

#[test]
fn failed_suffix_does_not_change_a_retained_checkpoint() {
    let recipe = recipe();
    let checkpoint = prefix(&recipe, 2);
    let before = export(checkpoint.built());
    let (failed, work) = with_operation_budget(0, || {
        evaluate_incremental(&recipe, Some(checkpoint.clone()), &mut Decider::default())
    });
    assert!(failed.is_err());
    assert!(work.exhausted);
    assert_eq!(export(checkpoint.built()), before);
    assert_same(&recipe, checkpoint);
}

#[test]
fn cached_prefix_is_paid_once_and_invalid_suffix_still_refuses() {
    let recipe = recipe();
    let checkpoint = evaluate_incremental(&recipe, None, &mut Decider::default()).unwrap();
    let (cached, work) = with_operation_budget(0, || {
        evaluate_incremental(&recipe, Some(checkpoint.clone()), &mut Decider::default())
    });
    assert!(!work.exhausted);
    assert_eq!(export(cached.unwrap().built()), export(checkpoint.built()));
    let mut invalid = recipe.clone();
    invalid.result = StepId::new("missing").unwrap();
    let full = evaluate(&invalid, &mut Decider::default()).unwrap_err();
    let cached =
        evaluate_incremental(&invalid, Some(checkpoint), &mut Decider::default()).unwrap_err();
    assert_eq!(full.to_string(), cached.to_string());
}
