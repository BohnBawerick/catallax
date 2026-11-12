//! Canonical P1 with the Boolean entry's construction shortcuts bypassed.
use catallax_boolean::{BoolOp, classify, imprint, intersect, keep, sew};
use catallax_features::{EvaluationCheckpoint, evaluate_incremental};
use catallax_numbers::{Decider, reuse_operation_enclosures, with_operation_budget};
use catallax_recipe::{BooleanKind, Operation, Recipe, ValueKind};
use catallax_topology::{BodyId, Model};

fn export(model: &Model, body: BodyId) -> Vec<u8> {
    assert!(model.check(body).is_empty());
    let (bytes, work) = with_operation_budget(100_000_000, || {
        let mut bytes = Vec::new();
        let answer = catallax_step::write_step(
            model,
            &[body],
            &Default::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        assert!(answer.marks.is_empty());
        bytes
    });
    assert!(!work.exhausted);
    bytes
}

fn names(model: &Model, body: BodyId) -> serde_json::Value {
    serde_json::json!({
        "faces": model.faces_of_body(body).into_iter()
            .map(|face| model.face_origins(face)).collect::<Vec<_>>(),
        "edges": model.edges_of_body(body).into_iter()
            .map(|edge| model.edge_name(body, edge)).collect::<Vec<_>>(),
        "vertices": model.vertices_of_body(body).into_iter()
            .map(|vertex| model.vertex_name(body, vertex)).collect::<Vec<_>>(),
    })
}

#[test]
#[ignore = "complete P1 construction and every intermediate STEP export"]
fn every_p1_boolean_matches_the_five_stages_without_construction_shortcuts() {
    let recipe = Recipe::parse(include_str!("../../../parts/P1-bracket.json")).unwrap();
    let mut checkpoint: Option<EvaluationCheckpoint> = None;
    let mut checked = 0;
    for (index, step) in recipe.steps.iter().enumerate() {
        if step.op.kind() != ValueKind::Body {
            continue;
        }
        let before = checkpoint.clone();
        let mut prefix = recipe.clone();
        prefix.steps.truncate(index + 1);
        prefix.result = step.id.clone();
        let (normal, normal_work) = with_operation_budget(100_000_000, || {
            evaluate_incremental(&prefix, checkpoint.take(), &mut Decider::default()).unwrap()
        });
        assert!(!normal_work.exhausted, "{}", step.id);
        if let Operation::Boolean { kind, target, tool } = &step.op {
            let before = before.as_ref().unwrap().built();
            let mut model = before.model.clone();
            let op = match kind {
                BooleanKind::Union => BoolOp::Unite,
                BooleanKind::Subtract => BoolOp::Subtract,
                BooleanKind::Intersect => BoolOp::Intersect,
            };
            let (body, work) = with_operation_budget(100_000_000, || {
                reuse_operation_enclosures();
                let mut dec = Decider::default();
                let graph = intersect(
                    &mut model,
                    before.bodies[target],
                    before.bodies[tool],
                    &mut dec,
                )
                .unwrap();
                assert!(graph.marks.is_empty());
                let pieces = imprint(&mut model, &graph.value, &mut dec).unwrap();
                assert!(pieces.marks.is_empty());
                let classes = classify(&mut model, &pieces.value, &mut dec).unwrap();
                assert!(classes.marks.is_empty());
                let selection = keep(&classes.value, &pieces.value, op);
                let answer = sew(&mut model, &pieces.value, &selection, &mut dec).unwrap();
                assert!(answer.marks.is_empty());
                answer.value
            });
            assert!(!work.exhausted, "{}", step.id);
            let expected = normal.built();
            assert_eq!(
                model.counts(body),
                expected.model.counts(expected.result),
                "{}",
                step.id
            );
            assert_eq!(
                names(&model, body),
                names(&expected.model, expected.result),
                "{}",
                step.id
            );
            assert_eq!(
                export(&model, body),
                export(&expected.model, expected.result),
                "{}",
                step.id
            );
            println!(
                "{}: normal={}, five_stages={}",
                step.id, normal_work.used, work.used
            );
            checked += 1;
        }
        checkpoint = Some(normal);
    }
    assert_eq!(checked, 16);
}
