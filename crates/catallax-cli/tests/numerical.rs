//! Numerical wire records cannot promote sampled distances or authorize export.

use catallax_cli::{numerical_build_summary, numerical_evidence_summary};
use catallax_numbers::{NumericalCheck, NumericalEvidence, NumericalSample, Rational, Refusal};
use serde_json::json;

#[test]
fn zero_samples_retain_their_locations_and_do_not_become_a_complete_bound() {
    let evidence = NumericalEvidence::new(
        Rational::parse("0.01").unwrap(),
        vec![NumericalCheck::Samples {
            name: "trim.edge-surface-distance".into(),
            scope: "roof/e0:0".into(),
            samples: vec![NumericalSample {
                parameters: vec![Rational::zero(), Rational::zero(), Rational::one()],
                distance_upper_mm: Rational::zero(),
            }],
        }],
        vec![],
    )
    .unwrap();
    let record = numerical_evidence_summary(&evidence);
    assert_eq!(record["tolerance_mm"], "1/100");
    assert_eq!(record["checks"][0]["sampled_max_mm"], "0");
    assert_eq!(
        record["checks"][0]["samples"][0]["parameters"],
        json!(["0", "0", "1"])
    );
    assert_eq!(record["repairs"], json!([]));
    for key in [
        "bounded",
        "measured_bound",
        "whole_body_bound",
        "proof",
        "exact",
    ] {
        assert!(record.get(key).is_none());
    }
}

#[test]
fn a_failed_check_never_serializes_partial_geometry_or_step_authority() {
    let result = Err(Refusal::new(
        "numerical.invalid-trim",
        "crossing parameter loop",
    ));
    let record = numerical_build_summary(&result);
    assert_eq!(record["format"], "catallax-numerical-result/1");
    assert_eq!(record["status"], "refused");
    assert_eq!(record["exact"], false);
    assert!(record["result"].is_null());
    assert!(record["numerical"].is_null());
    assert!(record["step"].is_null());
}

#[path = "../../catallax-app/tests/numerical_cases/mod.rs"]
mod numerical_cases;

#[test]
fn a_real_build_serializes_checks_and_repairs_but_never_a_step_result() {
    use catallax_app::{NumericalOptions, build_numerical};
    use catallax_numbers::Decider;
    use numerical_cases::{fixture, q, roof};
    let mut definition = fixture("0.001");
    let position = definition.vertices[0].position.clone();
    let key = definition.vertices[0].key.clone();
    definition.vertices[0].position[0] = &position[0] + &q("0.000001");
    let body = build_numerical(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let body = body
        .repair_vertex(&key, position, &mut Decider::default())
        .unwrap();
    let record = numerical_build_summary(&Ok(body.clone()));
    assert_eq!(record["status"], "numerical");
    assert_eq!(record["proof"], "numerical");
    assert_eq!(record["exact"], false);
    assert_eq!(record["result"]["faces"], 6);
    assert_eq!(record["numerical"]["tolerance_mm"], "1/100");
    assert_eq!(record["numerical"]["repairs"][0]["kind"], "vertex-moved");
    assert_eq!(record["numerical"]["repairs"][0]["vertex"], key);
    let samples = record["numerical"]["checks"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["kind"] == "samples")
        .unwrap();
    assert!(samples["sampled_max_mm"].is_string());
    assert!(samples["samples"].is_array());
    assert!(record["step"].is_null());
    assert!(record.get("bounded").is_none());
    let refused = body.replace_surface("roof", roof("0.1"), &mut Decider::default());
    let record = numerical_build_summary(&refused);
    assert_eq!(record["status"], "refused");
    assert_eq!(
        record["refusal"]["code"],
        "numbers.numerical-tolerance-exceeded"
    );
    assert!(record["result"].is_null());
    assert!(record["numerical"].is_null());
    assert!(record["step"].is_null());
}

#[test]
fn malformed_records_and_underreported_repairs_do_not_form_evidence() {
    use catallax_numbers::NumericalRepair;
    use numerical_cases::q;
    let passed = || NumericalCheck::Passed {
        name: "topology.structure".into(),
        scope: "body".into(),
    };
    for (tolerance, checks) in [
        (q("0"), vec![passed()]),
        (q("0.01"), vec![]),
        (
            q("0.01"),
            vec![NumericalCheck::Passed {
                name: " ".into(),
                scope: "body".into(),
            }],
        ),
    ] {
        assert_eq!(
            NumericalEvidence::new(tolerance, checks, vec![])
                .unwrap_err()
                .code,
            "numbers.invalid-numerical-evidence"
        );
    }
    for (distance, code) in [
        ("-1", "numbers.invalid-numerical-evidence"),
        ("0.02", "numbers.numerical-tolerance-exceeded"),
    ] {
        let check = NumericalCheck::Samples {
            name: "trim.edge-surface-distance".into(),
            scope: "roof".into(),
            samples: vec![NumericalSample {
                parameters: vec![q("0")],
                distance_upper_mm: q(distance),
            }],
        };
        assert_eq!(
            NumericalEvidence::new(q("0.01"), vec![check], vec![])
                .unwrap_err()
                .code,
            code
        );
    }
    let repair = NumericalRepair::VertexMoved {
        vertex: "v0".into(),
        before: [q("0"), q("0"), q("0")],
        after: [q("0.001"), q("0"), q("0")],
        distance_upper_mm: q("0"),
    };
    assert_eq!(
        NumericalEvidence::new(q("0.01"), vec![passed()], vec![repair])
            .unwrap_err()
            .code,
        "numbers.invalid-numerical-evidence"
    );
}
