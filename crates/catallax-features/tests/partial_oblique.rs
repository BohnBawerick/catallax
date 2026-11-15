//! Open oblique rims share their end section with the vertical wall.

use catallax_features::{Built, evaluate_incremental};
use catallax_numbers::{BudgetPolicy, Decider, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Curve;
use serde_json::{Value, json};

fn face_names(built: &Built) -> Vec<Vec<catallax_recipe::Origin>> {
    let mut names: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .map(|face| built.model.face_origins(face).to_vec())
        .collect();
    names.sort();
    names
}

#[test]
fn the_public_partial_oblique_rim_builds_and_exports_after_a_radius_edit() {
    let mut value: Value = serde_json::from_str(include_str!(
        "../../../evidence/phase1/P03/inputs/partial-oblique.json"
    ))
    .unwrap();
    let mut checkpoint = None;
    let mut original_names = None;
    let mut original_edges = None;
    for radius in ["0.5", "0.75"] {
        value["steps"][6]["radius"] = json!(radius);
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 1_000_000_000,
            ..Default::default()
        });
        let ((next, bytes), work) = with_operation_budget(1_000_000_000, || {
            let next = evaluate_incremental(&recipe, checkpoint.take(), &mut dec).unwrap();
            let built = next.built();
            let mut bytes = Vec::new();
            let report = catallax_step::write_step(
                &built.model,
                &[built.result],
                &Default::default(),
                &mut dec,
                &mut bytes,
            )
            .unwrap();
            assert_eq!(report.marks.len(), 1);
            assert_eq!(report.marks[0].code, "step.approximate-pipe-end");
            assert!(report.value.definition_deviation.is_none());
            assert!(String::from_utf8_lossy(&bytes).contains("complete_body_bound=unestablished"));
            assert!(!bytes.is_empty());
            (next, bytes)
        });
        assert!(!work.exhausted, "{radius}: {work:?}");
        let built = next.built();
        let mut query_dec = Decider::default();
        assert!(built.model.check(built.result).is_empty());
        assert_eq!(built.model.faces_of_body(built.result).len(), 5);
        let mut caps = 0;
        for edge in built.model.edges_of_body(built.result) {
            if let Curve::PatchSection(section) =
                built.model.geom().curve(built.model.edge(edge).curve)
            {
                caps += 1;
                assert_eq!(built.model.coedges_of_edge(edge).len(), 2);
                for point in section.endpoints() {
                    assert!(query_dec.sign(&point.x).is_coincident());
                }
            }
        }
        assert_eq!(caps, 2);
        let names = face_names(built);
        assert!(
            names
                .iter()
                .flatten()
                .any(|origin| matches!(origin.role, catallax_recipe::Role::Blend { .. }))
        );
        let mut edges: Vec<_> = built
            .model
            .edges_of_body(built.result)
            .into_iter()
            .filter_map(|edge| built.model.edge_name(built.result, edge))
            .collect();
        edges.sort_by(|a, b| a.between().cmp(b.between()));
        assert_eq!(edges.len(), 9);
        if let Some(original) = &original_names {
            assert_eq!(original, &names);
        }
        if let Some(original) = &original_edges {
            assert_eq!(original, &edges);
        }
        original_names = Some(names);
        original_edges = Some(edges);
        assert!(
            String::from_utf8(bytes)
                .unwrap()
                .contains("RATIONAL_B_SPLINE_SURFACE")
        );
        checkpoint = Some(next);
    }
}
