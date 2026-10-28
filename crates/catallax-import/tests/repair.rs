//! Repair, decisions and the change record, through `apply`.

mod common;

use catallax_import::{
    Approval, Decision, GapKind, ImportOptions, ImportRefusal, Imported, apply, propose, sha256_hex,
};
use catallax_measure::{Accuracy, volume};
use catallax_numbers::{Decider, Rational};
use catallax_topology::{BodyId, Model};
use common::*;

fn shifted_by(size: [&str; 3], face: usize, by: [&str; 3]) -> String {
    let mut b = BoxFixture::new(size);
    b.shift = Some((face, by.map(r)));
    box_step(b)
}

fn run(text: &str, decisions: &[Decision]) -> Result<(Model, Imported), ImportRefusal> {
    let mut model = Model::new();
    let imported = apply(
        &mut model,
        text.as_bytes(),
        "part.step",
        &ImportOptions::default(),
        decisions,
        None,
    )?;
    Ok((model, imported))
}

fn vol(model: &Model, body: BodyId) -> Rational {
    let v = volume(model, body, &Accuracy::default(), &mut Decider::default())
        .unwrap()
        .value;
    assert_eq!(v.lo, v.hi);
    v.lo
}

fn accept(id: &str, tol: &str) -> Decision {
    Decision::Accept {
        id: id.into(),
        tolerance_mm: tol.into(),
        approver: "jev".into(),
    }
}

fn decline(id: &str) -> Decision {
    Decision::Decline { id: id.into() }
}

#[test]
fn a_small_seam_is_repaired_automatically_and_recorded() {
    let text = shifted_by(["10", "10", "10"], 1, ["0.004", "0", "0"]);
    let (model, done) = run(&text, &[]).unwrap();
    let [(_, body)] = done.bodies.as_slice() else {
        panic!("one body")
    };
    assert!(model.check(*body).is_empty());
    // The seam closes by extending the four neighbours to the shifted face.
    assert_eq!(vol(&model, *body), r("1000.4"));
    assert!(done.left_open.is_empty());
    let [change] = done.record.changes.as_slice() else {
        panic!("one change")
    };
    assert_eq!(change.approval, Approval::Automatic);
    assert_eq!(change.kind, GapKind::SeamClose);
    assert_eq!(change.max_movement_mm, "0.004");
    assert_eq!(change.tolerance_mm, "0.01");
    assert_eq!(change.bounds_before[1][0], "10.004");
    assert_eq!(change.bounds_after[1][0], "10.004");
    assert_eq!(done.record.source.sha256, sha256_hex(text.as_bytes()));
    let derivative = done.derivative.expect("derivative");
    assert_eq!(
        done.record.derivative.unwrap().sha256,
        sha256_hex(&derivative)
    );
}

#[test]
fn the_measured_movement_is_the_true_three_dimensional_distance() {
    let text = shifted_by(["10", "10", "10"], 1, ["0.004", "0.003", "0"]);
    let (_, done) = run(&text, &[]).unwrap();
    assert_eq!(done.record.changes[0].max_movement_mm, "0.005");
}

#[test]
fn separate_bodies_stay_separate_and_unrecorded() {
    let mut b = BoxFixture::new(["1", "1", "1"]);
    b.origin = [r("1.005"), r("0"), r("0")];
    let (_, done) = run(&boxes_step(&[BoxFixture::new(["1", "1", "1"]), b]), &[]).unwrap();
    assert_eq!(done.bodies.len(), 2);
    assert!(done.record.changes.is_empty() && done.derivative.is_none());
}

#[test]
fn a_thin_opening_stays_open_and_only_a_decline_is_accepted() {
    let mut b = BoxFixture::new(["10", "10", "0.005"]);
    b.omit = Some(1);
    let text = box_step(b);
    assert!(
        matches!(run(&text, &[]).unwrap_err(), ImportRefusal::InvalidDecision { id, .. } if id == "g0")
    );
    assert!(matches!(
        run(&text, &[accept("g0", "0.01")]).unwrap_err(),
        ImportRefusal::InvalidDecision { detail, .. } if detail.contains("opening")
    ));
    let (_, done) = run(&text, &[decline("g0")]).unwrap();
    assert!(done.bodies.is_empty());
    assert_eq!(done.left_open.len(), 1);
    assert!(done.derivative.is_none() && done.record.derivative.is_none());
    let change = &done.record.changes[0];
    assert_eq!(change.approval, Approval::Declined);
    assert_eq!(change.max_movement_mm, "0");
    assert_eq!(change.bounds_before, change.bounds_after);
}

#[test]
fn an_over_cap_gap_asks_and_an_approval_is_recorded() {
    let text = shifted_by(["10", "10", "10"], 1, ["0.02", "0", "0"]);
    assert!(matches!(
        run(&text, &[]).unwrap_err(),
        ImportRefusal::InvalidDecision { .. }
    ));
    assert!(matches!(
        run(&text, &[accept("g0", "0.01")]).unwrap_err(),
        ImportRefusal::InvalidDecision { detail, .. } if detail.contains("below")
    ));
    let (model, done) = run(&text, &[accept("g0", "0.02")]).unwrap();
    assert_eq!(vol(&model, done.bodies[0].1), r("1002"));
    let change = &done.record.changes[0];
    assert_eq!(change.approval, Approval::User);
    assert_eq!(change.approver.as_deref(), Some("jev"));
    assert_eq!(change.tolerance_mm, "0.02");
    assert_eq!(change.max_movement_mm, "0.02");
}

#[test]
fn a_decline_leaves_the_original_unchanged_and_builds_nothing() {
    let text = shifted_by(["10", "10", "10"], 1, ["0.02", "0", "0"]);
    let (_, done) = run(&text, &[decline("g0")]).unwrap();
    assert!(done.bodies.is_empty() && done.derivative.is_none());
    assert_eq!(done.left_open.len(), 1);
    assert_eq!(done.record.changes[0].approval, Approval::Declined);
}

#[test]
fn an_automatic_repair_can_be_declined_but_not_accepted() {
    let text = shifted_by(["10", "10", "10"], 1, ["0.004", "0", "0"]);
    assert!(matches!(
        run(&text, &[accept("g0", "0.01")]).unwrap_err(),
        ImportRefusal::InvalidDecision { .. }
    ));
    let (_, done) = run(&text, &[decline("g0")]).unwrap();
    assert!(done.bodies.is_empty() && done.record.changes[0].approval == Approval::Declined);
}

#[test]
fn decisions_for_unknown_or_repeated_proposals_are_invalid() {
    let text = shifted_by(["10", "10", "10"], 1, ["0.02", "0", "0"]);
    assert!(matches!(
        run(&text, &[decline("g9")]).unwrap_err(),
        ImportRefusal::InvalidDecision { .. }
    ));
    assert!(matches!(
        run(&text, &[decline("g0"), decline("g0")]).unwrap_err(),
        ImportRefusal::InvalidDecision { .. }
    ));
}

#[test]
fn proposals_and_decisions_agree_on_the_question_id() {
    let text = shifted_by(["10", "10", "10"], 1, ["0.02", "0", "0"]);
    let set = propose(text.as_bytes(), &ImportOptions::default()).unwrap();
    assert_eq!(
        set.proposals[0].question,
        format!("repair:{}:g0", &set.source_sha256[..12])
    );
}

#[test]
fn a_gap_no_side_can_close_without_changing_a_surface_fails_validation() {
    let mut b = BoxFixture::new(["10", "10", "10"]);
    b.shift = Some((1, [r("0.004"), r("0"), r("0")]));
    b.also = Some((3, [r("0"), r("0.003"), r("0")]));
    let text = box_step(b);
    let set = propose(text.as_bytes(), &ImportOptions::default()).unwrap();
    assert_eq!(set.proposals.len(), 1);
    assert!(matches!(
        run(&text, &[]).unwrap_err(),
        ImportRefusal::ValidationFailed { id, .. } if id == "g0"
    ));
}
