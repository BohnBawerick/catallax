//! Wire records for rule-bound answers. Existing manifests omit these fields.

use catallax_numbers::{BoundedEvidence, CertificateScope, EvidenceUnit, OperationResult};
use serde_json::{Value, json};

/// A bounded answer's rule, preset level and measured kept-definition bound.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundedSummary(pub BoundedEvidence);

impl serde::Serialize for BoundedSummary {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let evidence = &self.0;
        let rule = evidence.rule();
        let definition = evidence.definition();
        let scope = match &definition.scope {
            CertificateScope::WholeSolid => json!({"kind":"whole-solid"}),
            CertificateScope::Entities(keys) => json!({"kind":"entities", "keys":keys}),
            CertificateScope::ParameterBox { chart, bounds } => json!({
                "kind":"parameter-box", "chart":chart,
                "bounds":bounds.iter().map(|(lo, hi)| [lo.to_string(), hi.to_string()]).collect::<Vec<_>>()
            }),
            CertificateScope::Query(key) => json!({"kind":"query", "key":key}),
        };
        let unit = match rule.unit {
            EvidenceUnit::Millimetres => "mm",
            EvidenceUnit::Radians => "rad",
            EvidenceUnit::Scalar => "scalar",
            EvidenceUnit::SquareMillimetres => "mm2",
            EvidenceUnit::CubicMillimetres => "mm3",
            EvidenceUnit::Raw => unreachable!("bounded evidence rejects unknown units"),
        };
        json!({
            "case_kind":rule.case_kind, "rule":rule.name,
            "level":evidence.level().to_string(), "measured_bound":evidence.measured().to_string(), "unit":unit,
            "definition":{"revision":definition.revision, "representation":definition.representation,
                "target":definition.target, "scope":scope}
        }).serialize(serializer)
    }
}

/// Summarizes operation delivery without promoting a bounded value to proven.
/// Existing adapters keep their existing wire formats until they adopt a rule.
pub fn operation_result_summary<T>(
    result: &OperationResult<T>,
    summarize: impl FnOnce(&T) -> Value,
) -> Value {
    match result {
        OperationResult::Proven(answer) => json!({
            "format":"catallax-operation-result/1", "status":"answered", "proof":"proven",
            "exact":answer.is_exact(), "result":summarize(&answer.value),
            "marks":answer.marks.iter().map(|mark| crate::mark_entry(mark, None)).collect::<Vec<_>>()
        }),
        OperationResult::Bounded(answer) => json!({
            "format":"catallax-operation-result/1", "status":"bounded", "proof":"bounded",
            "exact":false, "result":summarize(answer.value()),
            "bounded":BoundedSummary(answer.evidence().clone())
        }),
        OperationResult::Refused(refusal) => json!({
            "format":"catallax-operation-result/1", "status":"refused", "proof":"refused",
            "exact":false, "result":null,
            "refusal":{"code":refusal.code, "detail":refusal.detail}
        }),
    }
}
