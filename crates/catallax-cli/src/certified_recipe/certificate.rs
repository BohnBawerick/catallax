//! Wire records report checked claims; replay and geometry checking confer authority.

use std::fmt::Write as _;

use catallax_features::{CertificateKey, CertifiedBuilt};
use catallax_numbers::{
    CertificateIdentity, CertificateScope, DecisionEvidence, EvidenceBound, EvidenceUnit,
    ExceptionSource, FeatureSizeLowerBound, GeometricProof, LocalException,
};
use serde_json::{Value, json};

pub(super) fn key(value: &CertificateKey) -> String {
    let mut text = String::with_capacity(64);
    for byte in value {
        write!(text, "{byte:02x}").unwrap();
    }
    text
}
fn unit(unit: EvidenceUnit) -> &'static str {
    match unit {
        EvidenceUnit::Millimetres => "mm",
        EvidenceUnit::Radians => "rad",
        EvidenceUnit::Scalar => "scalar",
        EvidenceUnit::SquareMillimetres => "mm2",
        EvidenceUnit::CubicMillimetres => "mm3",
        EvidenceUnit::Raw => "raw",
    }
}
pub(super) fn bound(bound: &EvidenceBound) -> Value {
    json!({"lo":"0","hi":bound.upper.to_string(),"unit":unit(bound.unit)})
}
fn scope(scope: &CertificateScope<CertificateKey>) -> Value {
    match scope {
        CertificateScope::WholeSolid => json!({"kind":"whole-solid"}),
        CertificateScope::Entities(entities) => {
            json!({"kind":"complete-entities-and-trims","entities":entities.iter().map(key).collect::<Vec<_>>()})
        }
        CertificateScope::ParameterBox { chart, bounds } => {
            json!({"kind":"closed-parameter-box","chart":key(chart),
            "bounds":bounds.iter().map(|(lo, hi)| [lo.to_string(),hi.to_string()]).collect::<Vec<_>>()})
        }
        CertificateScope::Query(query) => json!({"kind":"query","query":key(query)}),
    }
}
fn identity(identity: &CertificateIdentity<CertificateKey>) -> Value {
    json!({"revision":key(&identity.revision),"representation":key(&identity.representation),"target":key(&identity.target),"scope":scope(&identity.scope)})
}
fn decision(decision: &DecisionEvidence<CertificateKey>) -> Value {
    json!({"original_code":decision.code,"enclosure":[decision.enclosure[0].to_string(),decision.enclosure[1].to_string()],
        "unit":unit(decision.unit),"entities":decision.entities.iter().map(key).collect::<Vec<_>>(),"branch":key(&decision.branch)})
}
fn proof(proof: &GeometricProof<CertificateKey>) -> Value {
    match proof {
        GeometricProof::ExactDefinition(definition) => {
            json!({"kind":"exact-definition","definition":key(definition)})
        }
        GeometricProof::Approximation {
            definition,
            theorem,
        } => json!({"kind":"approximation","definition":key(definition),"theorem":key(theorem)}),
        GeometricProof::Conversion {
            decision: original,
            theorem,
        } => {
            json!({"kind":"scoped-residual-conversion","decision":decision(original),"theorem":key(theorem)})
        }
        GeometricProof::Transfer(theorem) => {
            json!({"kind":"operation-transfer","theorem":key(theorem)})
        }
        GeometricProof::WholeResult(theorem) => {
            json!({"kind":"whole-oriented-result","theorem":key(theorem)})
        }
    }
}
fn size(size: &FeatureSizeLowerBound) -> Value {
    match size {
        FeatureSizeLowerBound::NoNeighbour => json!({"kind":"no-neighbour"}),
        FeatureSizeLowerBound::Millimetres(value) => {
            json!({"kind":"proved-lower-bound","mm":value.to_string()})
        }
    }
}
fn exception(exception: &LocalException<CertificateKey>) -> Value {
    let source = match &exception.source {
        ExceptionSource::KernelApproximation => json!({"kind":"kernel-approximation"}),
        ExceptionSource::BoundedDecision { decision } => {
            json!({"kind":"bounded-decision","decision":key(decision)})
        }
        ExceptionSource::Import {
            file,
            source_uncertainty,
            acceptance,
        } => json!({"kind":"import","file":key(file),
            "source_uncertainty":source_uncertainty.as_ref().map(bound),
            "acceptance":acceptance.as_ref().map(|a| json!({"file":key(&a.file),"limit":bound(&a.limit),"record":key(&a.record)}))}),
    };
    json!({"enclosure":{"lo":exception.lower_mm.to_string(),"hi":exception.bound.upper.to_string(),"unit":unit(exception.bound.unit)},
        "hard_limit":bound(&exception.limit),"source":source,"shortest_neighbouring_edge":size(&exception.shortest_edge),
        "narrowest_neighbouring_face":size(&exception.narrowest_face)})
}

pub(super) fn record(built: &CertifiedBuilt) -> Value {
    let certificate = built.result_certificate();
    let book = built.evidence();
    let model = &built.built().model;
    let regions: Vec<_> = certificate.faces().into_iter().flatten().enumerate().map(|(index, face)| {
        let region = certificate.region_for_face(face).expect("checked complete face-to-region binding");
        let role = ["lower-cap", "upper-cap", "side"][index];
        let domain = if let Some(spline) = certificate.spline() {
            let domain = spline.carrier().domain();
            if index == 2 {
                json!({"kind":"complete-finite-spline-strip","bounds":domain.map(|bounds| bounds.map(ToString::to_string)),
                    "periodic_identification":{"axis":"u","from":domain[0][0].to_string(),"to":domain[0][1].to_string(),"continuity":"C1"},
                    "physical_boundaries":certificate.rings().unwrap().map(|edge| edge.index())})
            } else {
                json!({"kind":"complete-convex-cap-disk","boundary":certificate.rings().unwrap()[index].index(),
                    "radial_parameter":["0","1"],"boundary_parameter":domain[0].map(ToString::to_string)})
            }
        } else {
            json!({"kind":if index == 2 {"complete-analytic-periodic-strip"} else {"complete-convex-cap-disk"},
                "physical_boundaries":certificate.rings().unwrap().map(|edge| edge.index())})
        };
        json!({"id":key(&region.identity.representation),"face":face.index(),"role":role,"origins":model.face_origins(face),
            "identity":identity(&region.identity),"evidence":key(&region.evidence),"domain":domain,"exception":region.exception.as_ref().map(exception)})
    }).collect();
    let evidence: Vec<_> = book.iter().map(|(id, record)| json!({"id":key(id),"identity":identity(&record.identity),
        "bound":record.bound.as_ref().map(bound),"proof":proof(&record.proof),
        "dependencies":record.dependencies.iter().map(key).collect::<Vec<_>>(),
        "context":{"producer":record.context.producer,"step":record.context.step.as_ref().map(key),"stage":record.context.stage}})).collect();
    let decisions: Vec<_> = built.decisions().iter().map(|(id, choice)| json!({"decision":key(id),"branch":key(&choice.branch),"evidence":key(&choice.evidence)})).collect();
    let spline = certificate.spline().map(|spline| {
        let carrier = spline.carrier();
        json!({"kind":"bspline","degrees":carrier.degrees(),"control_count":carrier.control_count(),
            "knots":carrier.knots().map(|knots| knots.iter().map(ToString::to_string).collect::<Vec<_>>()),
            "domain":carrier.domain().map(|bounds| bounds.map(ToString::to_string)),
            "u_identification":{"from":"0","to":"4","continuity":"C1","physical_edge":false},
            "boundary_curves":certificate.rings().unwrap().into_iter().enumerate().map(|(index, edge)| json!({
                "edge":edge.index(),"kind":"bspline","shared_surface_boundary":if index == 0 {"v-min"} else {"v-max"},
                "closed":true,"physical":true,"added_seam_vertices":0})).collect::<Vec<_>>()})
    });
    let record = certificate.record();
    json!({"identity":identity(&record.identity),"evidence":key(&record.evidence),"regions":regions,
        "uniform_signed_distance_bound":record.uniform_bound.as_ref().map(bound),
        "bounded_decisions":record.decisions.iter().map(key).collect::<Vec<_>>(),"decide_once":decisions,"evidence_records":evidence,
        "settings":{"modelling_tolerance_mm":book.tolerance().millimetres().to_string(),
            "kernel_fit_target_mm":book.local_settings().kernel_target_mm().to_string(),
            "producer_fit_target_mm":certificate.fit_target_mm().map(ToString::to_string),
            "kernel_hard_cap_mm":book.local_settings().kernel_cap_mm().to_string(),"import_hard_cap_mm":book.local_settings().import_cap_mm().to_string()},
        "retained_target":{"kind":"elliptic-normal-offset-prism","distance_mm":certificate.target().distance_mm().to_string(),
            "empty":certificate.target().is_empty(),"side":match certificate.target().definition().side() {
                catallax_shapes::OffsetSide::Inward=>"inward",catallax_shapes::OffsetSide::Outward=>"outward"}},
        "spline":spline})
}
