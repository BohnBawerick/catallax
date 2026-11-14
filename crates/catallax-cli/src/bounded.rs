//! Separate wire records for bounded answers and checked numerical builds.

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

/// Serializes actual numerical checks without creating a complete-bound field.
/// Parameter locations and pointwise distance enclosures remain explicit.
pub fn numerical_evidence_summary(evidence: &catallax_numbers::NumericalEvidence) -> Value {
    use catallax_numbers::{NumericalCheck, NumericalRepair};
    let checks: Vec<_> = evidence.checks().iter().map(|check| match check {
        NumericalCheck::Passed { name, scope } => json!({"kind":"passed", "name":name, "scope":scope}),
        NumericalCheck::Samples { name, scope, samples } => json!({
            "kind":"samples", "name":name, "scope":scope,
            "sampled_max_mm":samples.iter().map(|s| &s.distance_upper_mm).max().map(ToString::to_string),
            "samples":samples.iter().map(|sample| json!({
                "parameters":sample.parameters.iter().map(ToString::to_string).collect::<Vec<_>>(),
                "distance_upper_mm":sample.distance_upper_mm.to_string()
            })).collect::<Vec<_>>()
        }),
    }).collect();
    let repairs: Vec<_> = evidence
        .repairs()
        .iter()
        .map(|repair| match repair {
            NumericalRepair::VertexMoved {
                vertex,
                before,
                after,
                distance_upper_mm,
            } => json!({
                "kind":"vertex-moved", "vertex":vertex,
                "before_mm":before.iter().map(ToString::to_string).collect::<Vec<_>>(),
                "after_mm":after.iter().map(ToString::to_string).collect::<Vec<_>>(),
                "distance_upper_mm":distance_upper_mm.to_string()
            }),
        })
        .collect();
    json!({"tolerance_mm":evidence.tolerance_mm().to_string(), "checks":checks, "repairs":repairs})
}

/// Public numerical build delivery. This format is separate from exact manifests
/// and from complete bounded answers. It grants no export authority; `step` is null.
/// Refusals publish neither a body nor partial acceptance evidence. A retained
/// undecided decision keeps its code, detail and exact rational magnitude bound;
/// that residual does not grant geometric accuracy or body acceptance.
pub fn numerical_build_summary(
    result: &Result<catallax_app::NumericalBody, catallax_numbers::Refusal>,
) -> Value {
    match result {
        Ok(body) => {
            let counts = body.counts();
            json!({"format":"catallax-numerical-result/1", "status":"numerical", "proof":"numerical", "exact":false,
                "definition":body.definition().name,
                "result":{"lumps":counts.lumps,"shells":counts.shells,"faces":counts.faces,"edges":counts.edges,"vertices":counts.vertices},
                "numerical":numerical_evidence_summary(body.evidence()), "step":null})
        }
        Err(refusal) => {
            let mut reason = json!({"code":refusal.code,"detail":refusal.detail});
            if let Some(decision) = refusal.undecided() {
                reason["undecided"] = json!({
                    "code":decision.code, "detail":decision.detail,
                    "bound":decision.bound.to_string()
                });
            }
            json!({"format":"catallax-numerical-result/1", "status":"refused", "proof":"refused", "exact":false,
            "result":null, "numerical":null, "step":null, "refusal":reason})
        }
    }
}
