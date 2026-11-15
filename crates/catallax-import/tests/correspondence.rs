//! Replacement controls use the completed production import, not parser-only data.

mod common;

use catallax_import::{
    Decision, ImportOptions, Imported, apply,
    correspondence::{EntityKind, UnresolvedReason, compare},
    propose, read_step, rebuild, sha256_hex,
};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_topology::Model;
use common::{BoxFixture, box_step, box_with_placement, boxes_step};

fn imported(bytes: &[u8]) -> Imported {
    let (result, work) = with_operation_budget(20_000_000, || {
        apply(
            &mut Model::new(),
            bytes,
            "control.step",
            &ImportOptions::default(),
            &[],
            None,
        )
    });
    assert!(!work.exhausted);
    result.unwrap()
}

// The same public renumbering control as catallax-features/tests/import_edit.rs.
fn renumber(text: &str) -> String {
    let mut result = String::new();
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        result.push(c);
        if c == '#' {
            let mut digits = String::new();
            while chars.peek().is_some_and(char::is_ascii_digit) {
                digits.push(chars.next().unwrap());
            }
            result.push_str(&(digits.parse::<u64>().unwrap() + 10000).to_string());
        }
    }
    result
}

#[test]
fn renumbering_verifies_every_face_and_physical_edge() {
    let old = box_step(BoxFixture::new(["2", "3", "4"]));
    let new = renumber(&old);
    let a = imported(old.as_bytes());
    let b = imported(new.as_bytes());
    let report = compare(&a, &b, &mut Decider::default()).unwrap();
    assert_eq!(report.old_source_sha256, sha256_hex(old.as_bytes()));
    assert_eq!(report.new_source_sha256, sha256_hex(new.as_bytes()));
    assert_eq!(report.verified.len(), 18);
    assert_eq!(
        report
            .verified
            .iter()
            .filter(|m| m.old.kind == EntityKind::Face)
            .count(),
        6
    );
    assert_eq!(
        report
            .verified
            .iter()
            .filter(|m| m.old.kind == EntityKind::Edge)
            .count(),
        12
    );
    assert!(report.unresolved.is_empty());
    assert!(report.unmatched_new.is_empty());
    for entry in &report.verified {
        assert_eq!(entry.old.kind, entry.new.kind);
        assert_ne!(entry.old.key, entry.new.key);
        assert_ne!(entry.old.solid, entry.new.solid);
        let before: u64 = entry.old.key[1..].parse().unwrap();
        let after: u64 = entry.new.key[1..].parse().unwrap();
        assert_eq!(after, before + 10000);
    }
}

#[test]
fn reused_entity_numbers_do_not_bind_changed_geometry() {
    let old = box_step(BoxFixture::new(["2", "3", "4"]));
    let new = box_step(BoxFixture::new(["5", "3", "4"]));
    let report = compare(
        &imported(old.as_bytes()),
        &imported(new.as_bytes()),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(!report.verified.is_empty()); // unchanged x=0 face and its edges
    assert!(
        report
            .unresolved
            .iter()
            .any(|m| m.old.kind == EntityKind::Face)
    );
    assert!(
        report
            .unresolved
            .iter()
            .any(|m| m.old.kind == EntityKind::Edge)
    );
    assert!(
        report
            .unresolved
            .iter()
            .all(|m| m.reason == UnresolvedReason::ChangedOrRemoved)
    );
    let face = read_step(old.as_bytes()).unwrap().solids[0].faces[1]
        .name
        .clone();
    assert!(report.unresolved.iter().any(|m| m.old.key == face));
}

#[test]
fn repeated_geometry_in_either_file_stays_ambiguous() {
    let one = box_step(BoxFixture::new(["2", "3", "4"]));
    let two = boxes_step(&[
        BoxFixture::new(["2", "3", "4"]),
        BoxFixture::new(["2", "3", "4"]),
    ]);
    let one = imported(one.as_bytes());
    let two = imported(two.as_bytes());
    for (old, new, count, candidates) in [(&one, &two, 18, 2), (&two, &one, 36, 1)] {
        let report = compare(old, new, &mut Decider::default()).unwrap();
        assert!(report.verified.is_empty());
        assert_eq!(report.unresolved.len(), count);
        assert!(
            report
                .unresolved
                .iter()
                .all(|m| m.reason == UnresolvedReason::AmbiguousGeometry
                    && m.candidates.len() == candidates)
        );
    }
}

#[test]
fn multi_solid_reordering_uses_definitions_and_retains_membership() {
    let first = || BoxFixture::new(["2", "3", "4"]);
    let second = || {
        let mut fixture = BoxFixture::new(["1", "1", "1"]);
        fixture.origin[0] = common::r("10");
        fixture
    };
    let old = boxes_step(&[first(), second()]);
    let new = renumber(&boxes_step(&[second(), first()]));
    let old_file = read_step(old.as_bytes()).unwrap();
    let new_file = read_step(new.as_bytes()).unwrap();
    let report = compare(
        &imported(old.as_bytes()),
        &imported(new.as_bytes()),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(report.verified.len(), 36);
    assert!(report.unresolved.is_empty());
    for pair in report.verified {
        let old_index = old_file
            .solids
            .iter()
            .position(|s| s.name == pair.old.solid)
            .unwrap();
        let new_index = new_file
            .solids
            .iter()
            .position(|s| s.name == pair.new.solid)
            .unwrap();
        assert_eq!(old_index + new_index, 1);
    }
}

#[test]
fn recorded_planar_derivatives_compare_and_survive_reopen() {
    let old = box_with_placement(5, Some(["0.000001", "0", "1"]), Some(["1", "0", "0"]));
    let new = renumber(&old);
    let a = imported(old.as_bytes());
    let b = imported(new.as_bytes());
    assert!(!a.record.source_incidence.is_empty());
    let report = compare(&a, &b, &mut Decider::default()).unwrap();
    assert_eq!(report.verified.len(), 18);
    assert!(report.unresolved.is_empty());
    let reopened = rebuild(
        &mut Model::new(),
        old.as_bytes(),
        &ImportOptions::default(),
        &a.record,
        None,
    )
    .unwrap();
    assert_eq!(
        compare(&reopened, &b, &mut Decider::default()).unwrap(),
        report
    );
}

#[test]
fn a_declined_solid_has_only_unresolved_entries() {
    let source = box_with_placement(5, Some(["0.000001", "0", "1"]), Some(["1", "0", "0"]));
    let options = ImportOptions::default();
    let proposals = propose(source.as_bytes(), &options).unwrap();
    let partial = apply(
        &mut Model::new(),
        source.as_bytes(),
        "control.step",
        &options,
        &[Decision::Decline {
            id: proposals.proposals[0].id.clone(),
        }],
        None,
    )
    .unwrap();
    assert_eq!(partial.left_open.len(), 1);
    let report = compare(
        &partial,
        &imported(source.as_bytes()),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(report.verified.is_empty());
    assert_eq!(report.unresolved.len(), 18);
    assert!(
        report
            .unresolved
            .iter()
            .all(|m| m.reason == UnresolvedReason::UnverifiedGeometry)
    );
}

#[test]
fn unavailable_peers_prevent_uniqueness_in_either_import() {
    let source = boxes_step(&[
        BoxFixture::new(["10", "20", "30"]),
        BoxFixture::new(["10", "20", "30"]),
    ]);
    let placement = source
        .lines()
        .rfind(|line| line.contains("AXIS2_PLACEMENT_3D("))
        .unwrap();
    let axis = placement.split(',').nth(2).unwrap();
    let prefix = format!("{axis}=");
    let direction = source
        .lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap();
    let source = source.replace(direction, &format!("{prefix}DIRECTION('',(0.000001,0,1));"));
    let options = ImportOptions::default();
    let proposals = propose(source.as_bytes(), &options).unwrap();
    assert!(!proposals.proposals.is_empty());
    let decisions = proposals
        .proposals
        .iter()
        .map(|proposal| Decision::Decline {
            id: proposal.id.clone(),
        })
        .collect::<Vec<_>>();
    let partial = apply(
        &mut Model::new(),
        source.as_bytes(),
        "mixed.step",
        &options,
        &decisions,
        None,
    )
    .unwrap();
    assert_eq!(partial.bodies.len(), 1);
    assert_eq!(partial.left_open.len(), 1);
    let complete = imported(box_step(BoxFixture::new(["10", "20", "30"])).as_bytes());
    for (old, new, unavailable) in [(&complete, &partial, 0), (&partial, &complete, 18)] {
        let report = compare(old, new, &mut Decider::default()).unwrap();
        assert!(report.verified.is_empty());
        assert_eq!(report.unresolved.len(), 18 + unavailable);
        assert_eq!(
            report
                .unresolved
                .iter()
                .filter(|entry| entry.reason == UnresolvedReason::UnverifiedGeometry)
                .count(),
            unavailable
        );
        let mut known = report
            .unresolved
            .iter()
            .filter(|entry| entry.reason != UnresolvedReason::UnverifiedGeometry);
        assert_eq!(known.clone().count(), 18);
        assert!(known.all(|entry| {
            entry.reason == UnresolvedReason::UnverifiedMultiplicity && entry.candidates.len() == 1
        }));
    }
}

#[test]
fn quadratic_ambiguity_output_obeys_the_work_allowance() {
    let single = imported(box_step(BoxFixture::new(["2", "3", "4"])).as_bytes());
    let (result, baseline) = with_operation_budget(20_000_000, || {
        compare(&single, &single, &mut Decider::default())
    });
    assert_eq!(result.unwrap().verified.len(), 18);
    let copies = 40;
    let source = boxes_step(
        &std::iter::repeat_with(|| BoxFixture::new(["2", "3", "4"]))
            .take(copies)
            .collect::<Vec<_>>(),
    );
    let repeated = imported(source.as_bytes());
    // A linear multiple covers the linear inventories, but cannot authorize
    // the quadratic ambiguity candidate output of identical repeated solids.
    let (result, work) = with_operation_budget(baseline.used * copies as u64, || {
        compare(&repeated, &repeated, &mut Decider::default())
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}

#[test]
fn public_record_edits_cannot_change_the_compared_source_identity() {
    let source = box_step(BoxFixture::new(["2", "3", "4"]));
    let mut old = imported(source.as_bytes());
    let new = imported(renumber(&source).as_bytes());
    old.record.source.sha256 = "0".repeat(64);
    let report = compare(&old, &new, &mut Decider::default()).unwrap();
    assert_eq!(report.old_source_sha256, sha256_hex(source.as_bytes()));
}

#[test]
fn exhausted_comparison_returns_no_partial_report() {
    let source = box_step(BoxFixture::new(["2", "3", "4"]));
    let old = imported(source.as_bytes());
    let new = imported(renumber(&source).as_bytes());
    let (result, work) = with_operation_budget(1, || compare(&old, &new, &mut Decider::default()));
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert_eq!(
        compare(&old, &new, &mut Decider::default())
            .unwrap()
            .verified
            .len(),
        18
    );
}
