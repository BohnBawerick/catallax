//! Independent prefix construction coverage and work/refusal reporting.
use catallax_compare::{Construction, independent_steps};
use catallax_features::ImportFiles;
use catallax_recipe::Recipe;
use serde_json::json;

fn recipe() -> Recipe {
    let pose = json!({"rotation":["1","0","0","0"],"translation":["0","0","0"]});
    Recipe::parse(
        &json!({"format":"catallax-recipe/4","name":"steps","steps":[
        {"id":"a","op":"sphere","radius":"5","pose":pose},
        {"id":"b","op":"sphere","radius":"2","pose":pose},
        {"id":"smooth","op":"smooth_blend","kind":"union","target":"a","tool":"b","radius":"1"},
        {"id":"later","op":"offset","body":"smooth","distance":"1"}
    ],"result":"later"})
        .to_string(),
    )
    .unwrap()
}

#[test]
fn report_keeps_every_step_and_the_original_first_failure() {
    let rows = independent_steps(&recipe(), &ImportFiles::new(), 20_000_000).unwrap();
    assert_eq!(rows.len(), 4);
    for row in &rows[..2] {
        assert!(matches!(row.brep.construction, Construction::Built { .. }));
        assert!(matches!(row.field.construction, Construction::Built { .. }));
    }
    assert!(matches!(
        rows[2].field.construction,
        Construction::Built {
            faces: 0,
            empty: false
        }
    ));
    let Construction::Refused {
        blocked_at,
        refusal,
    } = &rows[2].brep.construction
    else {
        panic!("B-rep has no smooth builder");
    };
    assert_eq!(blocked_at.as_deref(), Some("smooth"));
    assert!(!refusal.code.is_empty());
    assert_eq!(rows[3].brep.construction, rows[2].brep.construction);
    assert!(
        rows.iter()
            .all(|r| r.brep.work <= r.brep.limit && r.field.work <= r.field.limit)
    );
}

#[test]
fn work_exhaustion_never_counts_as_built() {
    let rows = independent_steps(&recipe(), &ImportFiles::new(), 1).unwrap();
    assert!(rows.iter().all(
        |r| matches!(r.brep.construction, Construction::Refused { .. })
            && matches!(r.field.construction, Construction::Refused { .. })
    ));
}

#[test]
fn target_parts_report_all_authored_steps_without_cross_kernel_fallback() {
    // P5's authored radius-2 refusal is deliberately retained.
    let r = Recipe::parse(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    let rows = independent_steps(&r, &ImportFiles::new(), 20_000_000).unwrap();
    assert_eq!(rows.len(), r.steps.len());
    let last = rows.last().unwrap();
    let Construction::Refused { refusal, .. } = &last.field.construction else {
        panic!("radius-2 field verdict changed");
    };
    assert_eq!(refusal.code, "fillet.unsupported-overflow");
}
