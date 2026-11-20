//! Shared by the numerical STEP tests of each writer path.
use catallax_numbers::{Decider, NumericalCheck, NumericalEvidence, Rational};
use catallax_topology::numerical::{
    CopyRule, CopySource, CopyTag, NumericalBody, build_with_sampled_carriers,
};

/// The body again with every edge use, and every retained u/v curve, tagged
/// as a derived copy with the stated bound.
fn tagged(body: &NumericalBody, bound: &str) -> NumericalBody {
    let tag = CopyTag::Derived {
        from: CopySource::Edge3d,
        rule: CopyRule::ExactComposition,
        bound: Rational::parse(bound).unwrap(),
    };
    let mut definition = body.definition().clone();
    for trim in definition
        .faces
        .iter_mut()
        .flat_map(|face| face.loops.iter_mut().flatten())
    {
        trim.tag = tag.clone();
    }
    let mut native = body.native_trims().to_vec();
    for binding in &mut native {
        binding.tag = tag.clone();
    }
    build_with_sampled_carriers(
        definition,
        native,
        body.sampled_carriers().clone(),
        body.options().clone(),
        &mut Decider::default(),
    )
    .unwrap()
}

/// How many edge uses a set of checks records by a stated bound, and how many
/// by a measurement.
fn bound_and_measured(evidence: &NumericalEvidence) -> [usize; 2] {
    let count = |suffix: &str| {
        evidence
            .checks()
            .iter()
            .filter(|check| {
                matches!(check, NumericalCheck::Samples { name, .. }
                    if name.starts_with("trim.") && name.ends_with(suffix))
            })
            .count()
    };
    [
        count("edge-surface-derived-bound"),
        count("edge-surface-distance"),
    ]
}

/// STEP carries no copy tag (decision Q8, section 4.4 of the Checkpoint 2
/// decisions). The writer rounds both copies of every use, so the check of
/// the written values measures every use, also one the source body recorded
/// by its bound.
pub fn assert_the_written_body_is_measured(body: &NumericalBody, bound: &str) {
    let source = tagged(body, bound);
    let [stated, measured] = bound_and_measured(source.evidence());
    assert!(stated > 0, "no use of the source is recorded by its bound");
    let file = catallax_step::numerical::build(
        &source,
        &catallax_step::numerical::NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(
        bound_and_measured(&file.report().emitted_evidence),
        [0, stated + measured]
    );
}
