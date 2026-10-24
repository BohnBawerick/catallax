#![allow(missing_docs)]

use catallax_numbers::*;
use std::collections::{BTreeMap, BTreeSet};

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn mm(value: &str) -> EvidenceBound {
    EvidenceBound {
        upper: q(value),
        unit: EvidenceUnit::Millimetres,
    }
}
fn local(value: &str) -> LocalException<u32> {
    LocalException {
        lower_mm: q("0"),
        bound: mm(value),
        limit: mm("0.001"),
        source: ExceptionSource::KernelApproximation,
        shortest_edge: FeatureSizeLowerBound::Millimetres(q("1")),
        narrowest_face: FeatureSizeLowerBound::Millimetres(q("1")),
    }
}
fn identity(
    representation: u32,
    target: u32,
    scope: CertificateScope<u32>,
) -> CertificateIdentity<u32> {
    CertificateIdentity {
        revision: 10,
        representation,
        target,
        scope,
    }
}
fn context() -> ProofContext<u32> {
    ProofContext {
        producer: "features.offset",
        step: Some(40),
        stage: None,
    }
}
fn covered(
    value: &str,
    proof: GeometricProof<u32>,
    uniform: Option<EvidenceBound>,
) -> (EvidenceBook<u32>, GeometricCertificate<u32>) {
    let exact = identity(1, 101, CertificateScope::Entities(BTreeSet::from([1])));
    let fit = identity(2, 102, CertificateScope::Entities(BTreeSet::from([2])));
    let whole = identity(20, 30, CertificateScope::WholeSolid);
    let mut book = EvidenceBook::new(ModellingTolerance::default());
    let dec = Decider::default();
    for (key, evidence) in [
        (
            11,
            GeometricEvidence {
                identity: exact.clone(),
                bound: Some(mm("0")),
                proof: GeometricProof::ExactDefinition(101),
                dependencies: BTreeSet::new(),
                context: context(),
            },
        ),
        (
            12,
            GeometricEvidence {
                identity: fit.clone(),
                bound: Some(mm(value)),
                proof,
                dependencies: BTreeSet::new(),
                context: context(),
            },
        ),
        (
            99,
            GeometricEvidence {
                identity: whole.clone(),
                bound: uniform.clone(),
                proof: GeometricProof::WholeResult(900),
                dependencies: BTreeSet::from([11, 12]),
                context: context(),
            },
        ),
    ] {
        book.record(key, evidence, &dec).unwrap();
    }
    let result = GeometricCertificate {
        identity: whole,
        evidence: 99,
        regions: BTreeMap::from([
            (
                1,
                RegionCertificate {
                    identity: exact,
                    evidence: 11,
                    exception: None,
                },
            ),
            (
                2,
                RegionCertificate {
                    identity: fit,
                    evidence: 12,
                    exception: Some(local(value)),
                },
            ),
        ]),
        decisions: BTreeSet::new(),
        uniform_bound: uniform,
    };
    (book, result)
}
fn fit_proof() -> GeometricProof<u32> {
    GeometricProof::Approximation {
        definition: 102,
        theorem: 800,
    }
}
fn check(book: &EvidenceBook<u32>, result: &GeometricCertificate<u32>) -> Result<(), Refusal> {
    book.check_structure(result, &BTreeSet::from([1, 2]), &Decider::default())
}

#[test]
fn local_fit_target_is_not_a_strict_whole_result_cap() {
    let settings = LocalToleranceSettings::default();
    let tau = ModellingTolerance::default();
    for (value, target_met) in [
        ("0.00005", true),
        ("0.0001", true),
        ("0.00015", false),
        ("0.001", false),
    ] {
        let claim = local(value);
        claim.check(&tau, &settings).unwrap();
        assert_eq!(claim.kernel_target_met(&settings), Some(target_met));
    }
    assert_eq!(
        local("0.0010000000001")
            .check(&tau, &settings)
            .unwrap_err()
            .code,
        GEOMETRIC_GAP_EXCEEDED
    );
    let (book, result) = covered("0.00005", fit_proof(), None);
    check(&book, &result).unwrap();
    // A separately claimed full metric can also exceed 1 nm. Structure checking
    // does not prove that theorem, and never synthesises it from max(face bounds).
    let (book, result) = covered("0.00005", fit_proof(), Some(mm("0.00007")));
    check(&book, &result).unwrap();
}

#[test]
fn decision_only_exceptions_cannot_borrow_a_fit_allowance() {
    let settings = LocalToleranceSettings::default();
    let tau = ModellingTolerance::default();
    let mut claim = local("0.0000009");
    claim.source = ExceptionSource::BoundedDecision { decision: 77 };
    claim.limit = mm("0.000001");
    claim.check(&tau, &settings).unwrap();
    for value in ["0.000001", "0.000002"] {
        claim.bound = mm(value);
        assert_eq!(
            claim.check(&tau, &settings).unwrap_err().code,
            GEOMETRIC_GAP_EXCEEDED
        );
    }
    let decision = GeometricProof::Conversion {
        decision: DecisionEvidence {
            code: "numbers.precision-limit",
            enclosure: [q("-1"), q("1")],
            unit: EvidenceUnit::Raw,
            entities: BTreeSet::from([2]),
            branch: 7,
        },
        theorem: 801,
    };
    let (book, result) = covered("0.000002", decision, None);
    assert_eq!(
        check(&book, &result).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
}

#[test]
fn shorter_neighbours_tighten_caps_and_invalid_units_never_convert() {
    let settings = LocalToleranceSettings::default();
    let tau = ModellingTolerance::default();
    let mut claim = local("0.00005");
    claim.narrowest_face = FeatureSizeLowerBound::Millimetres(q("0.0005"));
    claim.check(&tau, &settings).unwrap();
    claim.narrowest_face = FeatureSizeLowerBound::Millimetres(q("0.0004"));
    assert_eq!(
        claim.check(&tau, &settings).unwrap_err().code,
        GEOMETRIC_GAP_EXCEEDED
    );
    claim = local("0.00005");
    claim.shortest_edge = FeatureSizeLowerBound::Millimetres(q("0"));
    assert_eq!(
        claim.check(&tau, &settings).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    claim = local("0.00005");
    claim.bound.unit = EvidenceUnit::SquareMillimetres;
    assert_eq!(
        claim.check(&tau, &settings).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    claim = local("0.00005");
    claim.lower_mm = q("0.00006");
    assert_eq!(
        claim.check(&tau, &settings).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
}

#[test]
fn import_override_is_explicit_file_specific_and_separate_from_uncertainty() {
    let settings = LocalToleranceSettings::default();
    let tau = ModellingTolerance::default();
    let mut claim = local("0.01");
    claim.limit = mm("0.01");
    claim.source = ExceptionSource::Import {
        file: 50,
        source_uncertainty: Some(mm("0.5")),
        acceptance: None,
    };
    claim.check(&tau, &settings).unwrap();
    assert_eq!(claim.kernel_target_met(&settings), None);
    claim.bound = mm("0.010000001");
    assert_eq!(
        claim.check(&tau, &settings).unwrap_err().code,
        GEOMETRIC_GAP_EXCEEDED
    );
    claim.limit = mm("0.02");
    claim.source = ExceptionSource::Import {
        file: 50,
        source_uncertainty: Some(mm("0.5")),
        acceptance: Some(ImportToleranceAcceptance {
            file: 50,
            limit: mm("0.02"),
            record: 51,
        }),
    };
    claim.check(&tau, &settings).unwrap();
    let ExceptionSource::Import {
        acceptance: Some(accepted),
        ..
    } = &mut claim.source
    else {
        unreachable!()
    };
    accepted.file = 60;
    assert_eq!(
        claim.check(&tau, &settings).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
}

#[test]
fn whole_result_requires_every_region_and_every_changed_trim_exception() {
    let (book, result) = covered("0.00005", fit_proof(), None);
    let mut missing = result.clone();
    missing.regions.remove(&1);
    assert_eq!(
        check(&book, &missing).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    let mut hidden = result.clone();
    hidden.regions.get_mut(&2).unwrap().exception = None;
    assert_eq!(
        check(&book, &hidden).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    let mut invented = result.clone();
    invented.regions.get_mut(&1).unwrap().exception = Some(local("0.00005"));
    assert_eq!(
        check(&book, &invented).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    let mut understated = result;
    understated
        .regions
        .get_mut(&2)
        .unwrap()
        .exception
        .as_mut()
        .unwrap()
        .bound = mm("0.00004");
    assert_eq!(
        check(&book, &understated).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
}

#[test]
fn metric_summary_and_stale_domain_require_their_own_bound_proofs() {
    let (book, result) = covered("0.00005", fit_proof(), None);
    let mut invented = result.clone();
    invented.uniform_bound = Some(mm("0.00005"));
    assert_eq!(
        check(&book, &invented).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    let mut stale = result.clone();
    stale.regions.get_mut(&2).unwrap().identity.revision += 1;
    assert_eq!(
        check(&book, &stale).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    let mut wrong_domain = result;
    wrong_domain.regions.get_mut(&2).unwrap().identity.scope = CertificateScope::Query(2);
    assert_eq!(
        check(&book, &wrong_domain).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
    let (book, result) = covered(
        "0.00005",
        GeometricProof::Approximation {
            definition: 103,
            theorem: 800,
        },
        None,
    );
    assert_eq!(
        check(&book, &result).unwrap_err().code,
        GEOMETRIC_GAP_UNPROVEN
    );
}
