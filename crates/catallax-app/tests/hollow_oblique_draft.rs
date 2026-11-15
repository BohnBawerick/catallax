//! Application build, edit, recipe persistence, picking and STEP for a hollow draft.
use catallax_app::build_incremental;
use catallax_numbers::Decider;
use catallax_recipe::Recipe;
use std::collections::BTreeSet;

#[test]
fn hollow_oblique_draft_edits_reopens_and_exports_with_physical_names() {
    let mut source: serde_json::Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/draft-annular/oblique-bound.json"
    ))
    .unwrap();
    let mut prefix = None;
    let mut original_names = None;
    let mut exports = Vec::new();
    for radius in ["8", "9", "9"] {
        source["steps"][0]["loops"][0][0]["radius"] = serde_json::json!(radius);
        let recipe = Recipe::parse(&source.to_string()).unwrap();
        let saved = recipe.to_text();
        let reopened = Recipe::parse(&saved).unwrap();
        assert_eq!(reopened, recipe);
        let app = build_incremental(
            &reopened,
            prefix,
            &Default::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(app.names.faces.len(), 4);
        assert_eq!(app.names.edges.len(), 4);
        let names: BTreeSet<_> = app.names.faces.values().flatten().cloned().collect();
        if let Some(before) = &original_names {
            assert_eq!(&names, before);
        } else {
            original_names = Some(names);
        }
        assert!(app.names.edges.values().all(Option::is_some));
        let display = &app.display.value;
        assert_eq!(
            display
                .faces
                .iter()
                .map(|f| f.face)
                .collect::<BTreeSet<_>>(),
            app.names.faces.keys().copied().collect()
        );
        assert_eq!(
            display
                .edges
                .iter()
                .map(|e| e.edge)
                .collect::<BTreeSet<_>>(),
            app.names.edges.keys().copied().collect()
        );
        assert!(
            display
                .positions
                .iter()
                .all(|p| p[0] * p[0] + p[1] * p[1] > 1.0)
        );
        let built = app.checkpoint.built();
        let mut bytes = Vec::new();
        let exported = catallax_step::write_step(
            &built.model,
            &[built.result],
            &Default::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        assert_eq!(exported.value.face_map.len(), 6);
        assert!(
            exported
                .value
                .face_map
                .iter()
                .all(|f| app.names.faces.contains_key(&f.face))
        );
        exports.push(bytes);
        prefix = Some(app.checkpoint);
    }
    assert_ne!(exports[0], exports[1]);
    assert_eq!(exports[1], exports[2]);
}
