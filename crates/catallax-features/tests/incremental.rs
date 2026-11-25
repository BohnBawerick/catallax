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
    assert_eq!(records(resumed), records(&full));
}

// Feature records by origin, in build order: fillets, chamfers, removed faces.
fn records(built: &Built) -> Vec<Vec<String>> {
    vec![
        built
            .fillets
            .iter()
            .map(|d| format!("{:?}", d.origin))
            .collect(),
        built
            .chamfers
            .iter()
            .map(|d| format!("{:?}", d.origin))
            .collect(),
        built
            .consumed
            .iter()
            .map(|c| format!("{:?} by {:?}", c.source, c.by))
            .collect(),
    ]
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

#[test]
fn undo_and_independent_branches_keep_names_facts_definitions_and_step_bytes() {
    let recipe = recipe();
    let saved = prefix(&recipe, 2);
    let before = export(saved.built());
    let definition = saved.built().recipe().clone();
    let bodies = saved.built().bodies.clone();
    let facts: Vec<_> = saved.built().model.geom().facts().iter().copied().collect();
    let names: Vec<_> = saved
        .built()
        .model
        .faces_of_body(saved.built().result)
        .into_iter()
        .map(|face| (face, saved.built().model.face_origins(face).to_vec()))
        .collect();
    let edges: Vec<_> = saved
        .built()
        .model
        .edges_of_body(saved.built().result)
        .into_iter()
        .map(|edge| {
            (
                edge,
                saved.built().model.edge_name(saved.built().result, edge),
            )
        })
        .collect();
    let mut branch = recipe.clone();
    if let catallax_recipe::Operation::Boolean { kind, .. } = &mut branch.steps[2].op {
        *kind = catallax_recipe::BooleanKind::Union;
    } else {
        panic!("boolean cut");
    }
    std::thread::scope(|scope| {
        let reader = scope.spawn(|| export(saved.built()));
        assert_same(&recipe, saved.clone());
        assert_same(&branch, saved.clone());
        assert_eq!(reader.join().unwrap(), before);
    });
    assert_eq!(export(saved.built()), before);
    assert_eq!(saved.built().recipe(), &definition);
    assert_eq!(saved.built().bodies, bodies);
    assert_eq!(
        saved
            .built()
            .model
            .geom()
            .facts()
            .iter()
            .copied()
            .collect::<Vec<_>>(),
        facts
    );
    // Extracting and editing a public model cannot mutate the checkpoint either.
    let mut extracted = saved.clone().into_built();
    extracted
        .model
        .map_face_origins(extracted.result, |origin| {
            let mut changed = origin.clone();
            changed.step = StepId::new("other").unwrap();
            changed
        });
    for (face, origins) in names {
        assert_eq!(saved.built().model.face_origins(face), origins);
    }
    for (edge, name) in edges {
        assert_eq!(
            saved.built().model.edge_name(saved.built().result, edge),
            name
        );
    }
    assert_eq!(export(saved.built()), before);
}

// Pattern and mirror copy the fillet and chamfer records of their source. A
// suffix that resumes after those records were made copies them as a fresh
// build does; an earlier edit discards the checkpoint, refusal included.
#[test]
fn copies_of_checkpointed_fillets_and_chamfers_match_fresh_builds() {
    let wall = json!({"step":"base","role":{"side":{"segment":"wall"}}});
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"placed features","steps":[
            {"id":"base","op":"cylinder","radius":"8","height":"12",
                "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"round","op":"fillet","body":"base","radius":"2",
                "edges":[{"between":[{"step":"base","role":"end"},wall]}]},
            {"id":"bevel","op":"chamfer","body":"round","distance":"1",
                "edges":[{"between":[{"step":"base","role":"start"},wall]}]},
            {"id":"copies","op":"pattern","body":"bevel",
                "layout":{"linear":{"vector":["30","0","0"],"count":2}}},
            {"id":"reflected","op":"mirror","body":"copies",
                "plane":{"point":["-20","0","0"],"normal":["1","0","0"]},"keep_original":true}
        ],"result":"reflected"})
        .to_string(),
    )
    .unwrap();
    let checkpoint = prefix(&recipe, 3);
    // The checkpoint already holds the source fillet and chamfer records.
    assert_eq!(checkpoint.built().fillets.len(), 1);
    assert_eq!(checkpoint.built().chamfers.len(), 1);
    assert_same(&recipe, checkpoint.clone());
    let full = evaluate(&recipe, &mut Decider::default()).unwrap();
    // One source, its pattern copy, and the mirror of both.
    assert_eq!(full.fillets.len(), 4);
    assert_eq!(full.chamfers.len(), 4);
    let mut edited = recipe.clone();
    if let catallax_recipe::Operation::Fillet { radius, .. } = &mut edited.steps[1].op {
        *radius = catallax_recipe::Number::new(3.into());
    } else {
        panic!("fillet");
    }
    assert_same(&edited, checkpoint.clone());
    if let catallax_recipe::Operation::Chamfer { distance, .. } = &mut edited.steps[2].op {
        *distance = catallax_recipe::Number::new(0.into());
    } else {
        panic!("chamfer");
    }
    let fresh = evaluate(&edited, &mut Decider::default()).unwrap_err();
    let resumed =
        evaluate_incremental(&edited, Some(checkpoint), &mut Decider::default()).unwrap_err();
    assert!(
        matches!(&fresh, catallax_features::BuildError::Geometry { step, .. } if step.as_str() == "bevel")
    );
    assert_eq!(fresh.to_string(), resumed.to_string());
}
