//! Round trips of the import wire types.

use catallax_import::*;

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
