//! A hole after a numerical fillet, through the build the application calls
//! (Checkpoint 2 decisions, section 4.5: "The `features.numerical-continuation`
//! refusal goes away: a hole after a numerical fillet works").
//!
//! The part is the owner's pocket in a rounded profile, with the six edges
//! of the pocket's rim rounded at 1 mm. The exact builders refuse that
//! fillet and the rolling ball builds it, so the body has fitted faces and
//! tolerant edges beside its exact ones.
use catallax_features::{BuildError, EvaluationOutcome, evaluate_incremental_outcome};
use catallax_numbers::{Decider, OperationWork, with_operation_budget};
use catallax_recipe::{HolePart, Recipe, Role};
use catallax_topology::numerical::tolerance::EntityKind;
use serde_json::{Value, json};

fn side(step: &str, segment: &str) -> Value {
    json!({"step": step, "role": {"side": {"segment": segment}}})
}

/// The pocket with its rim rounded, then `steps`.
fn rounded_pocket(steps: Vec<Value>) -> Recipe {
    let mut recipe: Value =
        serde_json::from_str(include_str!("fixtures/owner-curved-pocket.json")).unwrap();
    let rim: Vec<Value> = [
        ("top", "near"),
        ("top", "left"),
        ("top", "far"),
        ("right-round", "near"),
        ("right-round", "right"),
        ("right-round", "far"),
    ]
    .into_iter()
    .map(|(blank, tool)| json!({"between": [side("blank", blank), side("pocket-tool", tool)]}))
    .collect();
    let all = recipe["steps"].as_array_mut().unwrap();
    all.push(json!({"id": "round", "op": "fillet", "body": "pocket", "radius": "1", "edges": rim}));
    let result = steps
        .last()
        .map_or(json!("round"), |step| step["id"].clone());
    all.extend(steps);
    recipe["result"] = result;
    Recipe::parse(&recipe.to_string()).unwrap()
}

/// The application's 100M allowance, as the ordinary no-display build gives.
fn build(recipe: &Recipe) -> (Result<EvaluationOutcome, BuildError>, OperationWork) {
    with_operation_budget(100_000_000, || {
        evaluate_incremental_outcome(recipe, None, &Default::default(), &mut Decider::default())
    })
}

/// A drill along the part's length, from the end cap y = 0, 15 mm up and
/// 25 mm left of the middle: through both caps, away from the pocket.
fn lengthwise(depth: Value) -> Value {
    json!({"id": "drilled", "op": "hole", "target": "round",
        "at": [{"id": "a", "point": ["-25", "0", "15"], "direction": ["0", "-1", "0"]}],
        "diameter": "6", "depth": depth, "tip": "flat"})
}

#[test]
fn a_hole_after_a_numerical_fillet_builds() {
    let (fillet, fillet_work) = build(&rounded_pocket(Vec::new()));
    let EvaluationOutcome::Numerical(rounded) = fillet.unwrap() else {
        panic!("the rim fillet is built by the rolling ball");
    };
    let recipe = rounded_pocket(vec![lengthwise(json!("through"))]);
    let (outcome, work) = build(&recipe);
    assert!(!work.exhausted, "{work:?}");
    let EvaluationOutcome::Numerical(drilled) = outcome.unwrap() else {
        panic!("a hole in a body with tolerant entities stays one");
    };
    let body = drilled.body();
    assert_eq!(body.topology_faults(), Vec::new());
    let (before, after) = (rounded.body().counts(), body.counts());
    assert_eq!(
        (
            after.faces - before.faces,
            after.edges - before.edges,
            after.vertices - before.vertices
        ),
        (4, 12, 8)
    );
    // The walls are named as the hole's, and every other face keeps the
    // origins it had after the fillet.
    let walls: Vec<_> = drilled
        .face_origins()
        .iter()
        .filter(|(_, origins)| {
            origins
                .iter()
                .any(|origin| origin.step.as_str() == "drilled")
        })
        .collect();
    assert_eq!(walls.len(), 4);
    for (_, origins) in &walls {
        assert!(
            matches!(&origins[..], [origin] if matches!(&origin.role,
                Role::Hole { at, part: HolePart::Wall } if at.to_string() == "a")),
            "{origins:?}"
        );
    }
    for (key, origins) in rounded.face_origins() {
        assert_eq!(drilled.face_origins().get(key), Some(origins), "{key}");
    }
    assert_eq!(
        drilled.face_origins().len(),
        rounded.face_origins().len() + 4
    );
    // Exact and tolerant entities in one body.
    let edges = &body.definition().edges;
    let tolerant = edges
        .iter()
        .filter(|edge| body.stored_tolerance(EntityKind::Edge, &edge.key).is_some())
        .count();
    assert!(
        tolerant > 0 && tolerant < edges.len(),
        "{tolerant} of {}",
        edges.len()
    );
    for edge in edges
        .iter()
        .filter(|edge| edge.key.starts_with("drilled/a/"))
    {
        assert!(
            body.stored_tolerance(EntityKind::Edge, &edge.key).is_none(),
            "{}",
            edge.key
        );
    }
    println!(
        "the fillet: {} used; the fillet and the hole: {} used",
        fillet_work.used, work.used
    );
}

#[test]
fn what_has_no_construction_on_a_tolerant_body_refuses_by_its_own_name() {
    let refusal = |steps: Vec<Value>| {
        let (outcome, work) = build(&rounded_pocket(steps));
        assert!(!work.exhausted, "{work:?}");
        match outcome {
            Err(BuildError::Geometry { step, refusal }) => (step.as_str().to_owned(), refusal),
            other => panic!("expected a refusal, got {other:?}"),
        }
    };
    // A blind hole has a floor, which this construction does not make.
    let (step, blind) = refusal(vec![lengthwise(json!({"blind": "10"}))]);
    assert_eq!(
        (step.as_str(), blind.code),
        ("drilled", "hole.numerical-tool")
    );
    // An operation with no construction on such a body says so. The old
    // name is gone.
    let (step, again) = refusal(vec![json!({"id": "again", "op": "fillet", "body": "round",
        "radius": "1", "edges": [{"between": [side("blank", "top"), side("blank", "left-round")]}]})]);
    assert_eq!(
        (step.as_str(), again.code),
        ("again", "features.tolerant-input")
    );
    assert!(!again.detail.contains("numerical-continuation"));
}
