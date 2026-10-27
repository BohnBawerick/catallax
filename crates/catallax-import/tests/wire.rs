//! Round trips of the import wire types.

use catallax_import::*;

fn b(v: &str) -> Bounds {
    [
        [v.into(), v.into(), v.into()],
        [v.into(), v.into(), v.into()],
    ]
}

#[test]
fn record_round_trips_and_omits_empty_optionals() {
    let record = RepairRecord {
        format: REPAIR_FORMAT.into(),
        source: SourceRef {
            name: "a.step".into(),
            sha256: "ab".repeat(32),
        },
        derivative: None,
        changes: vec![Change {
            id: "g0".into(),
            entities: vec!["#12".into(), "#40".into()],
            kind: GapKind::SeamClose,
            bounds_before: b("0"),
            bounds_after: b("0"),
            max_movement_mm: "0.004".into(),
            tolerance_mm: AUTO_CAP_MM.into(),
            reason: "boundary gap below automatic cap".into(),
            approval: Approval::Declined,
            approver: None,
        }],
    };
    let text = serde_json::to_string(&record).unwrap();
    assert!(!text.contains("derivative") && !text.contains("approver"));
    assert!(text.contains("\"approval\":\"declined\"") && text.contains("seam_close"));
    assert_eq!(serde_json::from_str::<RepairRecord>(&text).unwrap(), record);
}

#[test]
fn decisions_and_refusals_use_tagged_json() {
    let accept: Decision = serde_json::from_str(
        r#"{"decision":"accept","id":"g1","tolerance_mm":"0.02","approver":"jev"}"#,
    )
    .unwrap();
    assert!(matches!(accept, Decision::Accept { ref id, .. } if id == "g1"));
    let decline = serde_json::to_string(&Decision::Decline { id: "g1".into() }).unwrap();
    assert_eq!(decline, r#"{"decision":"decline","id":"g1"}"#);
    let refusal = serde_json::to_string(&ImportRefusal::Unsupported {
        entity: "#40".into(),
        type_name: "B_SPLINE_SURFACE_WITH_KNOTS".into(),
    })
    .unwrap();
    assert!(refusal.starts_with(r#"{"code":"unsupported""#));
}
