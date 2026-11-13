use catallax_cli::{BoundedSummary, Status, StepSummary, operation_result_summary};
use catallax_numbers::{
    Answer, BoundedAnswer, BoundedEvidence, BoundedRule, CertificateIdentity, CertificateScope,
    EvidenceBound, EvidenceUnit, OperationResult, Rational, Refusal,
};
use serde_json::json;

static RULE: BoundedRule = BoundedRule {
    case_kind: "test.whole-result",
    name: "test.whole-result/1",
    level: "0.01",
    unit: EvidenceUnit::Millimetres,
};

fn evidence() -> BoundedEvidence {
    BoundedEvidence::new(
        &RULE,
        CertificateIdentity {
            revision: "recipe/1".into(),
            representation: "body/1".into(),
            target: "kept-law/1".into(),
            scope: CertificateScope::WholeSolid,
        },
        EvidenceBound {
            upper: Rational::parse("0.005").unwrap(),
            unit: EvidenceUnit::Millimetres,
        },
    )
    .unwrap()
}

#[test]
fn bounded_wire_record_keeps_exact_rationals_rule_and_definition() {
    let bounded = serde_json::to_value(BoundedSummary(evidence())).unwrap();
    assert_eq!(
        bounded,
        json!({"case_kind":"test.whole-result", "rule":"test.whole-result/1",
        "level":"1/100", "measured_bound":"1/200", "unit":"mm",
        "definition":{"revision":"recipe/1", "representation":"body/1", "target":"kept-law/1",
            "scope":{"kind":"whole-solid"}}})
    );
    let result = OperationResult::Bounded(BoundedAnswer::new(7, evidence()));
    let summary = operation_result_summary(&result, |value| json!(value));
    assert_eq!(summary["status"], "bounded");
    assert_eq!(summary["proof"], "bounded");
    assert_eq!(summary["exact"], false);
    assert_eq!(summary["bounded"], bounded);
    assert!(summary.get("refusal").is_none());
    assert_eq!(
        serde_json::to_string(&Status::Bounded).unwrap(),
        "\"bounded\""
    );
    let proven = operation_result_summary(&OperationResult::Proven(Answer::exact(7)), |v| json!(v));
    assert_eq!(proven["proof"], "proven");
    assert!(proven.get("bounded").is_none());
    let refused: OperationResult<u32> =
        OperationResult::Refused(Refusal::new("test.no-bound", "no bound"));
    let refused = operation_result_summary(&refused, |v| json!(v));
    assert_eq!(refused["proof"], "refused");
    assert!(refused.get("bounded").is_none());
}

#[test]
fn step_summary_retains_bounded_status_and_legacy_bytes() {
    let report = catallax_step::StepReport::default();
    let proven = StepSummary::from_result(
        &OperationResult::Proven(Answer::exact(report.clone())),
        "answer.step",
        None,
    )
    .unwrap();
    assert_eq!(
        serde_json::to_string(&proven).unwrap(),
        "{\"file\":\"answer.step\",\"max_deviation\":\"0\",\"approximated_curves\":0,\"added_vertices\":0,\"added_edges\":0}"
    );
    let bounded = OperationResult::Bounded(BoundedAnswer::new(report, evidence()));
    let bounded = StepSummary::from_result(&bounded, "answer.step", None).unwrap();
    assert_eq!(
        serde_json::to_value(bounded).unwrap()["bounded"]["case_kind"],
        "test.whole-result"
    );
    let refused = OperationResult::Refused(Refusal::new("step.no-bound", "no bound"));
    assert_eq!(
        StepSummary::from_result(&refused, "answer.step", None)
            .unwrap_err()
            .code,
        "step.no-bound"
    );
}

#[test]
fn old_refusal_manifest_has_identical_bytes_and_bounded_manifest_is_separate() {
    let manifest = catallax_cli::failure(None, Status::Refused, "test.no-bound", "no bound");
    let before = catallax_cli::render_manifest(&manifest);
    assert_eq!(
        serde_json::to_string(&manifest).unwrap(),
        format!(
            "{{\"question_id\":null,\"kernel\":{},\"status\":\"refused\",\"exact\":false,\"result\":null,\"properties\":null,\"step\":null,\"marks\":[],\"refusal\":{{\"code\":\"test.no-bound\",\"step\":null,\"stage\":null,\"detail\":\"no bound\"}},\"stats\":{}}}",
            serde_json::to_string(&manifest.kernel).unwrap(),
            serde_json::to_string(&manifest.stats).unwrap(),
        )
    );
    assert!(!before.contains("bounded"));
    assert!(before.contains("\"status\": \"refused\""));
    assert_eq!(
        serde_json::to_string(manifest.refusal.as_ref().unwrap()).unwrap(),
        "{\"code\":\"test.no-bound\",\"step\":null,\"stage\":null,\"detail\":\"no bound\"}"
    );
    let bounded = manifest.with_bounded_evidence(evidence());
    assert_eq!(bounded.status, Status::Bounded);
    assert!(!bounded.exact);
    assert!(bounded.refusal.is_none());
    assert!(catallax_cli::render_manifest(&bounded).contains("\"status\": \"bounded\""));
}

#[test]
fn old_proven_manifest_has_identical_bytes() {
    let mut manifest = catallax_cli::failure(None, Status::Answered, "unused", "unused");
    manifest.exact = true;
    manifest.refusal = None;
    assert_eq!(
        serde_json::to_string(&manifest).unwrap(),
        format!(
            "{{\"question_id\":null,\"kernel\":{},\"status\":\"answered\",\"exact\":true,\"result\":null,\"properties\":null,\"step\":null,\"marks\":[],\"refusal\":null,\"stats\":{}}}",
            serde_json::to_string(&manifest.kernel).unwrap(),
            serde_json::to_string(&manifest.stats).unwrap()
        )
    );
}
