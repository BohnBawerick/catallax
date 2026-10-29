//! Rebuilding an import from its source and persisted record.

mod common;

use catallax_import::{
    Decision, ImportOptions, ImportRefusal, Imported, RepairRecord, apply, rebuild, sha256_hex,
};
use catallax_numbers::Rational;
use catallax_topology::Model;
use common::*;

fn text(origin: [&str; 3], unit: Unit, by: &str) -> String {
    let mut b = BoxFixture::new(["10", "12", "14"]);
    b.origin = origin.map(r);
    b.unit = unit;
    b.shift = Some((1, [r(by), Rational::zero(), Rational::zero()]));
    box_step(b)
}

fn first(text: &str, decisions: &[Decision]) -> Imported {
    apply(
        &mut Model::new(),
        text.as_bytes(),
        "p.step",
        &ImportOptions::default(),
        decisions,
        None,
    )
    .unwrap()
}

fn reopen(record: &RepairRecord) -> RepairRecord {
    serde_json::from_str(&serde_json::to_string(record).unwrap()).unwrap()
}

fn again(text: &str, record: &RepairRecord) -> Result<Imported, ImportRefusal> {
    rebuild(
        &mut Model::new(),
        text.as_bytes(),
        &ImportOptions::default(),
        record,
        None,
    )
}

#[test]
fn an_approved_record_survives_reopen_and_reproduces_the_derivative() {
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.02");
    let done = first(
        &t,
        &[Decision::Accept {
            id: "g0".into(),
            tolerance_mm: "0.02".into(),
            approver: "jev".into(),
        }],
    );
    let record = reopen(&done.record);
    assert_eq!(record, done.record);
    let rebuilt = again(&t, &record).unwrap();
    assert_eq!(rebuilt.derivative, done.derivative);
    assert_eq!(rebuilt.record, done.record);
    assert_eq!(rebuilt.bodies.len(), 1);
}

#[test]
fn automatic_and_declined_records_rebuild_too() {
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.004");
    let done = first(&t, &[]);
    assert_eq!(
        again(&t, &reopen(&done.record)).unwrap().derivative,
        done.derivative
    );
    let declined = first(&t, &[Decision::Decline { id: "g0".into() }]);
    let back = again(&t, &reopen(&declined.record)).unwrap();
    assert!(back.bodies.is_empty() && back.derivative.is_none());
}

#[test]
fn other_bytes_or_a_forged_record_are_refused() {
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.004");
    let record = first(&t, &[]).record;
    let other = t.replace("'box'", "'bax'");
    assert!(matches!(
        again(&other, &record).unwrap_err(),
        ImportRefusal::HashMismatch { .. }
    ));
    let mut forged = record.clone();
    forged.derivative.as_mut().unwrap().sha256 = sha256_hex(b"x");
    assert!(matches!(
        again(&t, &forged).unwrap_err(),
        ImportRefusal::HashMismatch { expected, .. } if expected == sha256_hex(b"x")
    ));
    let mut moved = record.clone();
    moved.changes[0].max_movement_mm = "0.001".into();
    assert!(matches!(
        again(&t, &moved).unwrap_err(),
        ImportRefusal::ValidationFailed { .. }
    ));
    let mut format = record;
    format.format = "catallax-import-repair/9".into();
    assert!(matches!(
        again(&t, &format).unwrap_err(),
        ImportRefusal::Malformed { .. }
    ));
}

#[test]
fn moved_scaled_and_threaded_copies_repair_the_same_way() {
    // A translated copy and an inch copy of the same part repair by the same movement.
    let far = text(["100.5", "-30", "7"], Unit::Millimetre, "0.004");
    assert_eq!(first(&far, &[]).record.changes[0].max_movement_mm, "0.004");
    // 0.004 mm is 0.004/25.4 inch, a terminating decimal 0.00015748...: use an exact inch shift.
    let inch = text(["0", "0", "0"], Unit::Inch, "0.0001");
    let done = first(&inch, &[]);
    assert_eq!(done.record.changes[0].max_movement_mm, "0.00254");
    // Four threads give byte-identical derivatives and hashes.
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.004");
    let runs: Vec<_> = (0..4)
        .map(|_| {
            let t = t.clone();
            std::thread::spawn(move || first(&t, &[]).derivative.unwrap())
        })
        .collect();
    let all: Vec<_> = runs.into_iter().map(|h| h.join().unwrap()).collect();
    assert!(all.windows(2).all(|w| w[0] == w[1]));
}

#[test]
fn ties_at_the_cap_repair_and_just_above_ask() {
    let at = first(&text(["0", "0", "0"], Unit::Millimetre, "0.01"), &[]);
    assert_eq!(at.record.changes[0].max_movement_mm, "0.01");
    let t = text(["0", "0", "0"], Unit::Millimetre, "0.010001");
    assert!(matches!(
        apply(
            &mut Model::new(),
            t.as_bytes(),
            "p.step",
            &ImportOptions::default(),
            &[],
            None
        )
        .unwrap_err(),
        ImportRefusal::InvalidDecision { .. }
    ));
}
