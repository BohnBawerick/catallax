//! Geometric measurements and whole-candidate publication, not ledger-only repairs.
mod common;

use catallax_import::{
    build_file, sha256_hex,
    tolerant::{
        self, ToleranceOptions,
        healing::{self, HealingChange},
    },
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_topology::Model;
use common::{BoxFixture, box_step, r};

fn source(shift: &str) -> String {
    let mut fixture = BoxFixture::new(["2", "3", "4"]);
    fixture.shift = Some((0, [r("0"), r("0"), r(shift)]));
    box_step(fixture)
}
fn read(text: &str) -> tolerant::TolerantFile {
    tolerant::read_step(
        text.as_bytes(),
        &sha256_hex(text.as_bytes()),
        ToleranceOptions::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

#[test]
fn a_shared_edge_source_is_sewn_with_measured_geometry_and_an_approximate_result() {
    let text = source("0.001");
    let input = read(&text);
    let before = input.clone();
    assert_eq!(input.bodies()[0].edges().len(), 16);
    assert_eq!(input.bodies()[0].vertices().len(), 12);
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(input, before);
    assert_eq!(healed.source(), &before);
    assert!(healed.report().approximate);
    assert_eq!(
        healed
            .report()
            .changes
            .iter()
            .filter(|c| matches!(c, HealingChange::EdgesSewn { .. }))
            .count(),
        4
    );
    assert_eq!(
        healed
            .report()
            .changes
            .iter()
            .filter(|c| matches!(c, HealingChange::VerticesJoined { .. }))
            .count(),
        4
    );
    assert_eq!(
        healed
            .report()
            .stages
            .iter()
            .map(|s| s.stage.as_str())
            .collect::<Vec<_>>(),
        ["sewing", "wire", "face", "shell", "solid"]
    );
    assert!(!healed.report().measurements.is_empty());
    let result = healed.definition(&mut Decider::default()).unwrap();
    let mut model = Model::new();
    let bodies = build_file(&mut model, &result, None).unwrap();
    assert_eq!(bodies.len(), 1);
    let counts = model.counts(bodies[0]);
    assert_eq!((counts.faces, counts.edges, counts.vertices), (6, 12, 8));
    let record = healed.tolerances();
    assert!(!record.bodies[0].changes.is_empty());
    assert!(
        record.bodies[0]
            .changes
            .iter()
            .all(|c| r(&c.after_mm) <= r("0.01"))
    );
    assert!(before.bodies()[0].tolerance_changes().is_empty());
}

#[test]
fn a_late_geometry_cap_failure_publishes_neither_earlier_sewing_nor_growth() {
    let mut text = source("0.001");
    let line = text
        .lines()
        .find(|l| l.contains("=LINE("))
        .unwrap()
        .to_owned();
    let origin = line.split(",#").nth(1).unwrap().split(',').next().unwrap();
    let changed = line.replacen(&format!(",#{origin}"), ",#900001", 1);
    text = text.replacen(&line, &changed, 1).replacen(
        "ENDSEC;\nEND-ISO-10303-21;",
        "#900001=CARTESIAN_POINT('',(100.,100.,100.));\nENDSEC;\nEND-ISO-10303-21;",
        1,
    );
    let input = read(&text);
    let before = input.clone();
    let error = healing::heal(&input, &mut Decider::default()).unwrap_err();
    assert!(matches!(error,tolerant::Error::Numerical(ref r) if r.code=="import.healing-cap"));
    assert_eq!(input, before);
    assert!(input.bodies()[0].tolerance_changes().is_empty());
}

#[test]
fn over_cap_open_geometry_is_retained_with_explicit_shell_issues() {
    let input = read(&source("0.02"));
    let before = input.clone();
    let healed = healing::heal(&input, &mut Decider::default()).unwrap();
    assert_eq!(input, before);
    assert!(
        !healed
            .report()
            .changes
            .iter()
            .any(|c| matches!(c, HealingChange::EdgesSewn { .. }))
    );
    assert!(
        healed
            .report()
            .issues
            .iter()
            .any(|i| i.code == "shell.open-or-nonmanifold")
    );
}

#[test]
fn an_exhausted_candidate_preserves_the_original_source_and_ledger() {
    let input = read(&source("0.001"));
    let before = input.clone();
    let (complete, measured) = with_operation_budget(20_000_000, || {
        healing::heal(&input, &mut Decider::default())
    });
    complete.unwrap();
    for allowance in [0, 1, measured.used / 2, measured.used - 1] {
        let (result, work) =
            with_operation_budget(allowance, || healing::heal(&input, &mut Decider::default()));
        assert!(result.is_err() && work.exhausted);
        assert_eq!(input, before);
    }
}
