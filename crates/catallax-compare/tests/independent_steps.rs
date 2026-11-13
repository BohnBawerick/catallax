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
fn zero_work_construction_may_build_with_one_unit() {
    let mut part = recipe();
    part.steps.truncate(1);
    part.result = part.steps[0].id.clone();
    let rows = independent_steps(&part, &ImportFiles::new(), 1).unwrap();
    for build in [&rows[0].brep, &rows[0].field] {
        assert!(matches!(build.construction, Construction::Built { .. }));
        assert_eq!(build.work, 0);
        assert_eq!(build.limit, 1);
    }
}

#[test]
fn work_exhaustion_never_counts_as_built() {
    let part = Recipe::parse(&json!({
        "format":"catallax-recipe/4","name":"metered-prism","steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
                {"kind":"line","id":"bottom","to":["4","0"]},
                {"kind":"line","id":"right","to":["4","3"]},
                {"kind":"line","id":"top","to":["0","3"]},
                {"kind":"line","id":"left","to":["0","0"]}
            ]]},
            {"id":"solid","op":"extrude","profile":"profile","distance":"2"}
        ],"result":"solid"
    }).to_string()).unwrap();
    let complete = independent_steps(&part, &ImportFiles::new(), 20_000_000).unwrap();
    assert!(matches!(
        complete[1].brep.construction,
        Construction::Built { .. }
    ));
    assert!(matches!(
        complete[1].field.construction,
        Construction::Built { .. }
    ));
    for field in [false, true] {
        let part = part.clone();
        let (_, work) = std::thread::spawn(move || {
            catallax_numbers::with_operation_budget(1, || {
                let mut dec = catallax_numbers::Decider::default();
                if field {
                    catallax_field::Plan::compile_checked(&part, &mut dec).is_ok()
                } else {
                    catallax_features::evaluate(&part, &mut dec).is_ok()
                }
            })
        })
        .join()
        .unwrap();
        assert!(work.exhausted, "field={field}: {work:?}");
    }
    let rows = independent_steps(&part, &ImportFiles::new(), 1).unwrap();
    for build in [&rows[1].brep, &rows[1].field] {
        assert!(
            matches!(build.construction, Construction::Refused { .. }),
            "{build:?}"
        );
        assert_eq!(build.limit, 1);
        assert!(build.work <= build.limit);
    }
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
