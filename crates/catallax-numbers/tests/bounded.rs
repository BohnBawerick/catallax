use catallax_numbers::{
    Answer, BoundedAnswer, BoundedEvidence, BoundedRule, CertificateIdentity, CertificateScope,
    EvidenceBound, EvidenceUnit, OperationResult, Rational, Refusal,
};

// This rule is a test fixture, not a supported geometric operation.
static RULE: BoundedRule = BoundedRule {
    case_kind: "test.whole-result",
    name: "test.whole-result/1",
    level: "0.01",
    unit: EvidenceUnit::Millimetres,
};

fn definition() -> CertificateIdentity<String> {
    CertificateIdentity {
        revision: "recipe/1".into(),
        representation: "body/1".into(),
        target: "kept-law/1".into(),
        scope: CertificateScope::WholeSolid,
    }
}

fn measured(text: &str) -> EvidenceBound {
    EvidenceBound {
        upper: Rational::parse(text).unwrap(),
        unit: EvidenceUnit::Millimetres,
    }
}

#[test]
fn level_is_the_static_rule_and_bounded_is_never_proven_even_at_zero() {
    for bound in ["0", "0.005", "0.01"] {
        let evidence = BoundedEvidence::new(&RULE, definition(), measured(bound)).unwrap();
        assert_eq!(evidence.level(), &Rational::parse("0.01").unwrap());
        assert_eq!(evidence.measured(), &Rational::parse(bound).unwrap());
        assert_eq!(evidence.definition(), &definition());
        assert!(
            evidence
                .to_string()
                .starts_with("bounded: test.whole-result under test.whole-result/1:")
        );
        let result = OperationResult::Bounded(BoundedAnswer::new(7, evidence));
        assert!(result.proven().is_none());
        assert!(result.to_string().starts_with("bounded:"));
    }
    let proved = OperationResult::from(Ok(Answer::exact(7)));
    assert_eq!(proved.proven().unwrap().value, 7);
    assert_eq!(proved.to_string(), "proven: 7");
    let refused: OperationResult<u32> = Err(Refusal::new("test.refused", "no bound")).into();
    assert!(refused.proven().is_none());
    assert_eq!(refused.to_string(), "refused: test.refused: no bound");
}

#[test]
fn rejects_exceeded_negative_unknown_or_mismatched_bounds_and_missing_definitions() {
    assert_eq!(
        BoundedEvidence::new(&RULE, definition(), measured("0.010001"))
            .unwrap_err()
            .code,
        "numbers.bounded-level-exceeded"
    );
    assert!(BoundedEvidence::new(&RULE, definition(), measured("-0.001")).is_err());
    for unit in [EvidenceUnit::Raw, EvidenceUnit::Radians] {
        let mut bound = measured("0.001");
        bound.unit = unit;
        assert!(BoundedEvidence::new(&RULE, definition(), bound).is_err());
    }
    let mut missing = definition();
    missing.target.clear();
    assert!(BoundedEvidence::new(&RULE, missing, measured("0.001")).is_err());
    for scope in [
        CertificateScope::Entities(Default::default()),
        CertificateScope::Query(String::new()),
        CertificateScope::ParameterBox {
            chart: "chart/1".into(),
            bounds: vec![],
        },
    ] {
        let mut missing = definition();
        missing.scope = scope;
        assert!(BoundedEvidence::new(&RULE, missing, measured("0.001")).is_err());
    }
}
