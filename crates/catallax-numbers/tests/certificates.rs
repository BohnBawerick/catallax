#![allow(missing_docs)]

use catallax_numbers::*;
use std::collections::BTreeSet;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn bound(value: &str) -> EvidenceBound {
    EvidenceBound {
        upper: q(value),
        unit: EvidenceUnit::Millimetres,
    }
}
fn identity() -> CertificateIdentity<u32> {
    CertificateIdentity {
        revision: 10,
        representation: 20,
        target: 30,
        scope: CertificateScope::WholeSolid,
    }
}
fn evidence(dependencies: &[u32]) -> GeometricEvidence<u32> {
    GeometricEvidence {
        identity: identity(),
        bound: Some(bound("0.0000008")),
        proof: GeometricProof::WholeResult(100),
        dependencies: dependencies.iter().copied().collect(),
        context: ProofContext {
            producer: "boolean",
            step: Some(40),
            stage: Some("sew"),
        },
    }
}
fn certificate() -> GapCertificate<u32> {
    GapCertificate {
        identity: identity(),
        bound: bound("0.0000009"),
        evidence: 3,
    }
}

#[test]
fn strict_cap_rejects_equality_and_never_clamps_or_converts_units() {
    let setting = ModellingTolerance::default();
    let tau = q("0.000001");
    let tiny = q("1e-100");
    assert_eq!(setting.millimetres(), &tau);
    setting
        .check(&EvidenceBound {
            upper: &tau - &tiny,
            unit: EvidenceUnit::Millimetres,
        })
        .unwrap();
    for upper in [tau.clone(), &tau + &tiny] {
        let input = EvidenceBound {
            upper: upper.clone(),
            unit: EvidenceUnit::Millimetres,
        };
        assert_eq!(
            setting.check(&input).unwrap_err().code,
            GEOMETRIC_GAP_EXCEEDED
        );
        assert_eq!(input.upper, upper);
    }
    for unit in [
        EvidenceUnit::Radians,
        EvidenceUnit::Scalar,
        EvidenceUnit::SquareMillimetres,
        EvidenceUnit::CubicMillimetres,
        EvidenceUnit::Raw,
    ] {
        assert_eq!(
            setting
                .check(&EvidenceBound {
                    upper: q("0"),
                    unit
                })
                .unwrap_err()
                .code,
            GEOMETRIC_GAP_UNPROVEN
        );
    }
    assert_eq!(
        setting.check(&bound("-1e-20")).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
}

#[test]
fn consecutive_bounds_add_instead_of_taking_a_maximum() {
    let first = bound("0.0000004");
    let second = first.add(&first).unwrap();
    assert_eq!(second.upper, q("0.0000008"));
    ModellingTolerance::default().check(&second).unwrap();
    let at_limit = second.add(&bound("0.0000002")).unwrap();
    assert_eq!(
        ModellingTolerance::default()
            .check(&at_limit)
            .unwrap_err()
            .code,
        GEOMETRIC_GAP_EXCEEDED
    );
    assert!(
        first
            .add(&EvidenceBound {
                upper: q("1"),
                unit: EvidenceUnit::Radians
            })
            .is_err()
    );
    assert!(first.add(&bound("-1")).is_err());
}

#[test]
fn shared_real_enclosure_is_outward_and_does_not_promote_raw_units() {
    let value = -Real::from_i64(2).sqrt();
    let enclosed = EvidenceBound::enclose_magnitude(
        &value,
        EvidenceUnit::Raw,
        &q("1e-12"),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(&enclosed.upper * &enclosed.upper >= q("2"));
    assert!(enclosed.upper < q("1.415"));
    assert_eq!(
        ModellingTolerance::default()
            .check(&enclosed)
            .unwrap_err()
            .code,
        GEOMETRIC_GAP_UNPROVEN
    );
    assert!(
        EvidenceBound::enclose_magnitude(
            &value,
            EvidenceUnit::Millimetres,
            &q("0"),
            &mut Decider::default(),
        )
        .is_err()
    );
}

#[test]
fn evidence_shares_identity_and_orders_records_independently_of_insertion() {
    let dec = Decider::default();
    let mut forward = EvidenceBook::new(ModellingTolerance::default());
    let mut reverse = EvidenceBook::new(ModellingTolerance::default());
    // Both parents refer to one shared leaf. Traversal must not report a cycle.
    let records = [
        (0, evidence(&[])),
        (1, evidence(&[0])),
        (2, evidence(&[0])),
        (3, evidence(&[1, 2])),
    ];
    for (id, record) in &records {
        forward.record(*id, record.clone(), &dec).unwrap();
    }
    for (id, record) in records.iter().rev() {
        reverse.record(*id, record.clone(), &dec).unwrap();
    }
    assert_eq!(
        forward.iter().collect::<Vec<_>>(),
        reverse.iter().collect::<Vec<_>>()
    );
    forward.check_gap_structure(&certificate(), &dec).unwrap();
    reverse.check_gap_structure(&certificate(), &dec).unwrap();
    assert!(!forward.record(0, evidence(&[]), &dec).unwrap());
    assert_eq!(
        forward.record(0, evidence(&[1]), &dec).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    assert_eq!(forward.get(&0).unwrap().dependencies, BTreeSet::new());
}

#[test]
fn missing_cyclic_or_incomplete_evidence_cannot_pass_structure_checking() {
    let dec = Decider::default();
    for records in [
        vec![],
        vec![(3, evidence(&[1]))],
        vec![(3, evidence(&[1])), (1, evidence(&[3]))],
    ] {
        let mut book = EvidenceBook::new(ModellingTolerance::default());
        for (id, record) in records {
            book.record(id, record, &dec).unwrap();
        }
        assert_eq!(
            book.check_gap_structure(&certificate(), &dec)
                .unwrap_err()
                .code,
            GEOMETRIC_GAP_UNPROVEN
        );
    }
    let mut record = evidence(&[]);
    record.proof = GeometricProof::Conversion {
        decision: DecisionEvidence {
            code: "numbers.precision-limit",
            enclosure: [q("1"), q("-1")],
            unit: EvidenceUnit::Raw,
            entities: BTreeSet::from([50]),
            branch: 60,
        },
        theorem: 70,
    };
    let mut book = EvidenceBook::new(ModellingTolerance::default());
    book.record(3, record, &dec).unwrap();
    assert!(book.check_gap_structure(&certificate(), &dec).is_err());
}

#[test]
fn certificate_binding_domains_and_understated_bounds_are_rejected() {
    let dec = Decider::default();
    let mut book = EvidenceBook::new(ModellingTolerance::default());
    book.record(3, evidence(&[]), &dec).unwrap();
    let original = certificate();
    let mut changed = original.clone();
    changed.identity.revision += 1;
    assert!(book.check_gap_structure(&changed, &dec).is_err());
    let mut changed = original.clone();
    changed.identity.representation += 1;
    assert!(book.check_gap_structure(&changed, &dec).is_err());
    let mut changed = original.clone();
    changed.identity.target += 1;
    assert!(book.check_gap_structure(&changed, &dec).is_err());
    let mut changed = original;
    changed.bound = bound("0.0000007");
    assert!(book.check_gap_structure(&changed, &dec).is_err());
    for scope in [
        CertificateScope::Entities(BTreeSet::new()),
        CertificateScope::ParameterBox {
            chart: 80,
            bounds: vec![],
        },
        CertificateScope::ParameterBox {
            chart: 80,
            bounds: vec![(q("1"), q("0"))],
        },
    ] {
        let mut record = evidence(&[]);
        record.identity.scope = scope.clone();
        let mut claim = certificate();
        claim.identity.scope = scope;
        let mut book = EvidenceBook::new(ModellingTolerance::default());
        book.record(3, record, &dec).unwrap();
        assert!(book.check_gap_structure(&claim, &dec).is_err());
    }
}

#[test]
fn decide_once_reuses_the_settled_branch_and_rejects_conflicting_imports() {
    let mut book = DecisionBook::default();
    let mut dec = Decider::default();
    let mut calls = 0;
    for branch in [7, 8, 9] {
        let choice = book
            .decide_once(1, &mut dec, |_| {
                calls += 1;
                Ok(RecordedChoice {
                    branch,
                    evidence: 2,
                })
            })
            .unwrap();
        assert_eq!(choice.branch, 7);
    }
    assert_eq!(calls, 1);
    assert!(
        !book
            .record(
                1,
                RecordedChoice {
                    branch: 7,
                    evidence: 2
                },
                &dec
            )
            .unwrap()
    );
    assert!(
        book.record(
            1,
            RecordedChoice {
                branch: 8,
                evidence: 2
            },
            &dec
        )
        .is_err()
    );
    assert!(
        book.record(
            1,
            RecordedChoice {
                branch: 7,
                evidence: 3
            },
            &dec
        )
        .is_err()
    );
    assert_eq!(book.get(&1).unwrap().branch, 7);
    book.record(
        0,
        RecordedChoice {
            branch: 4,
            evidence: 5,
        },
        &dec,
    )
    .unwrap();
    assert_eq!(
        book.iter().map(|(key, _)| *key).collect::<Vec<_>>(),
        vec![0, 1]
    );
}

#[test]
fn refusal_or_cancellation_never_commits_a_partial_decision() {
    let mut book = DecisionBook::default();
    let mut dec = Decider::default();
    let refused = book.decide_once(1, &mut dec, |_| {
        Err(Refusal::new("test.refused", "no choice"))
    });
    assert_eq!(refused.unwrap_err().code, "test.refused");
    assert!(book.get(&1).is_none());
    // The closure deliberately catches exhaustion. The record must still fail closed.
    let (result, work) = with_operation_budget(1, || {
        book.decide_once(1, &mut dec, |dec| {
            let _ = dec.check_deadline();
            Ok(RecordedChoice {
                branch: 7,
                evidence: 2,
            })
        })
        .map(|_| ())
    });
    assert!(result.is_err());
    assert!(work.exhausted);
    assert!(book.get(&1).is_none());
    book.record(
        1,
        RecordedChoice {
            branch: 7,
            evidence: 2,
        },
        &dec,
    )
    .unwrap();
    let (result, work) = with_operation_budget(0, || {
        book.decide_once(1, &mut dec, |_| {
            panic!("a recorded branch must not be recomputed")
        })
        .map(|_| ())
    });
    assert!(result.is_err());
    assert!(work.exhausted);
    assert_eq!(book.get(&1).unwrap().branch, 7);
}
