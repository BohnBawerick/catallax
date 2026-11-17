//! Tier selection, retained siblings and replay from source plus controls.
mod common;

use catallax_import::{
    numerical::StepImport,
    sha256_hex,
    tolerant::{
        self,
        opened::{self, BodyView, Tier},
    },
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::StepId;
use common::{BoxFixture, box_step, boxes_step, r};

fn open(source: &str) -> opened::OpenedStep {
    opened::open_step(
        source.as_bytes(),
        &sha256_hex(source.as_bytes()),
        &StepId::new("source").unwrap(),
        Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

fn history(file: &opened::OpenedStep) -> serde_json::Value {
    if let Some(native) = file.tolerant_admission() {
        return serde_json::to_value(native.evidence()).unwrap();
    }
    if let Some(healed) = file.bodies()[0].healing() {
        return serde_json::json!({"healing":healed.report(),"tolerances":healed.tolerances()});
    }
    match file.admission(0).unwrap() {
        StepImport::Numerical(admitted) => serde_json::to_value(admitted.record()).unwrap(),
        _ => panic!("normal-gap repair must retain approximate history"),
    }
}

fn checked_artifact(file: &opened::OpenedStep) -> ((usize, usize, usize), Vec<u8>) {
    let mut dec = Decider::default();
    match file.body(0).unwrap() {
        BodyView::Certified { model, body } => {
            let counts = model.counts(body);
            let mut bytes = Vec::new();
            catallax_step::write_step(model, &[body], &Default::default(), &mut dec, &mut bytes)
                .unwrap();
            ((counts.faces, counts.edges, counts.vertices), bytes)
        }
        BodyView::Numerical(body) => {
            let counts = body.counts();
            let written =
                catallax_step::numerical::build(body, &Default::default(), &mut dec).unwrap();
            written.check_written(written.bytes(), &mut dec).unwrap();
            (
                (counts.faces, counts.edges, counts.vertices),
                written.bytes().to_vec(),
            )
        }
        _ => panic!("stricter checked upgrade required"),
    }
}

#[test]
fn a_good_sibling_remains_solid_while_an_open_body_retains_all_five_faces() {
    let good = BoxFixture::new(["2", "3", "4"]);
    let mut sheet = BoxFixture::new(["2", "3", "4"]);
    sheet.origin = [r("5"), r("0"), r("0")];
    sheet.omit = Some(0);
    let source = boxes_step(&[good, sheet]);
    let opened = open(&source);
    assert_eq!(opened.bodies().len(), 2);
    assert_eq!(opened.bodies()[0].tier(), Tier::ValidSolid);
    assert!(!opened.bodies()[0].approximate());
    assert!(matches!(opened.body(0), Some(BodyView::Certified { .. })));
    assert_eq!(opened.bodies()[1].tier(), Tier::SurfaceBody);
    let Some(BodyView::Surfaces { source, .. }) = opened.body(1) else {
        panic!("source surface representation required")
    };
    let source = source.materialize(&mut Decider::default()).unwrap();
    assert_eq!(source.shells()[0].faces.len(), 5);
    assert_ne!(opened.source().bodies()[0].name(), source.name());
}

#[test]
fn admission_exhaustion_keeps_the_complete_prepaid_surface_representation() {
    let source = box_step(BoxFixture::new(["2", "3", "4"]));
    let hash = sha256_hex(source.as_bytes());
    let (baseline, _) = with_operation_budget(20_000_000, || open(&source));
    let (result, work) = with_operation_budget(baseline.read_work_used() + 128, || {
        opened::open_step(
            source.as_bytes(),
            &hash,
            &StepId::new("source").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    let opened = result.unwrap();
    assert!(work.exhausted);
    assert_eq!(opened.source(), baseline.source());
    assert_eq!(opened.bodies()[0].tier(), Tier::SurfaceBody);
    assert!(opened.first_refusal().is_some());
    assert!(matches!(
        opened.body(0),
        Some(BodyView::Surfaces {
            candidate: None,
            ..
        })
    ));
}

#[test]
fn normal_gap_replays_from_saved_source_controls_with_approximate_history() {
    let mut fixture = BoxFixture::new(["2", "3", "4"]);
    fixture.shift = Some((0, [r("0.001"), r("0"), r("0")]));
    let source = box_step(fixture);
    let original = open(&source);
    let controls = serde_json::to_vec(&original.saved_controls()).unwrap();
    let controls = serde_json::from_slice(&controls).unwrap();
    let replayed = opened::reopen(
        source.as_bytes(),
        &controls,
        &StepId::new("source").unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    for file in [&original, &replayed] {
        assert_eq!(file.bodies()[0].tier(), Tier::RepairedSolid);
        assert!(file.bodies()[0].approximate());
        assert_eq!(checked_artifact(file).0, (6, 12, 8));
    }
    let retained = history(&original);
    assert_eq!(retained, history(&replayed));
    assert!(
        retained.get("healing").is_some()
            || ["repairs", "changes", "face_tolerances", "source_incidence"]
                .iter()
                .any(|key| retained[*key]
                    .as_array()
                    .is_some_and(|values| !values.is_empty()))
    );
}

#[test]
fn retained_surface_sewing_requires_closure_within_the_callers_cap() {
    let mut fixture = BoxFixture::new(["2", "3", "4"]);
    fixture.shift = Some((0, [r("0.001"), r("0"), r("0")]));
    let source = box_step(fixture)
        .lines()
        .map(|line| {
            if line.contains("MANIFOLD_SOLID_BREP") {
                let (id, _) = line.split_once('=').unwrap();
                let shell = line.split(",#").nth(1).unwrap().trim_end_matches(");");
                format!("{id}=SHELL_BASED_SURFACE_MODEL('',(#{shell}));")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let hash = sha256_hex(source.as_bytes());
    let step = StepId::new("source").unwrap();
    for (cap, expected) in [
        ("0.002", Tier::RepairedSolid),
        ("0.0005", Tier::SurfaceBody),
    ] {
        let options = tolerant::ToleranceOptions {
            maximum_mm: r(cap),
            ..Default::default()
        };
        let opened = opened::open_step(
            source.as_bytes(),
            &hash,
            &step,
            options,
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(opened.bodies()[0].tier(), expected);
        assert!(!opened.source().bodies()[0].declared_solid());
        assert_eq!(
            opened.source().bodies()[0]
                .materialize(&mut Decider::default())
                .unwrap()
                .shells()[0]
                .faces
                .len(),
            6
        );
        assert!(opened.source().bodies()[0].tolerance_changes().is_empty());
        if expected == Tier::SurfaceBody {
            assert!(matches!(opened.body(0), Some(BodyView::Surfaces { .. })));
            assert!(opened.bodies()[0].refusal().is_some());
            continue;
        }
        let (counts, exported) = checked_artifact(&opened);
        assert_eq!(counts, (6, 12, 8));
        let replayed = opened::reopen(
            source.as_bytes(),
            &opened.saved_controls(),
            &step,
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(replayed.bodies()[0].tier(), Tier::RepairedSolid);
        assert_eq!(history(&opened), history(&replayed));
        if let Ok(folder) = std::env::var("CATALLAX_HEALING_EVIDENCE") {
            std::fs::write(
                std::path::Path::new(&folder).join("stitched-surface.step"),
                &exported,
            )
            .unwrap();
        }
    }
}

#[test]
fn default_opening_preserves_checked_material_in_a_multi_body_document() {
    let fixtures = (0..24)
        .map(|index| {
            let mut fixture = BoxFixture::new(["2", "3", "4"]);
            fixture.origin[0] = r(&(index * 5).to_string());
            fixture
        })
        .collect::<Vec<_>>();
    let source = boxes_step(&fixtures);
    let hash = sha256_hex(source.as_bytes());
    let (checked, checked_work) = with_operation_budget(20_000_000, || {
        opened::open_checked_step(
            source.as_bytes(),
            &hash,
            &StepId::new("source").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
        .unwrap()
    });
    assert!(!checked_work.exhausted);
    let (default, default_work) = with_operation_budget(20_000_000, || open(&source));
    assert!(!default_work.exhausted);
    for index in 0..fixtures.len() {
        assert_eq!(checked.bodies()[index].tier(), Tier::ValidSolid);
        assert_eq!(
            default.bodies()[index].tier(),
            checked.bodies()[index].tier()
        );
        let Some(BodyView::Certified { model, body }) = default.body(index) else {
            panic!("default opening downgraded previously checked material");
        };
        assert_eq!(model.body(body).kind, catallax_topology::BodyKind::Solid);
    }
}
