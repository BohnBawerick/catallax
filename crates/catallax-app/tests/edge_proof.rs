use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use catallax_app::{
    EdgeProofFailure, EdgeProofKey, EdgeProofOutcome, EvaluatedApplicationOutcome,
    build_incremental_with_files_without_display, prove_edge_gaps,
};
use catallax_numbers::{Decider, Rational, set_threads};
use catallax_recipe::Recipe;
use catallax_topology::Tolerance;
use serde_json::Value;

#[test]
fn whole_edge_reports_cover_exact_parts_and_a_fitted_loft() {
    let corpus: Value =
        serde_json::from_str(include_str!("../../../tools/predictor/corpus.json")).unwrap();
    let loft = corpus
        .as_object()
        .unwrap()
        .values()
        .find_map(|v| {
            v.as_array()?
                .iter()
                .find(|item| item["id"] == "loft-smooth-fillet-r1_10")
        })
        .unwrap()["inline"]
        .to_string();
    let limit = Tolerance::new(Rational::parse("0.01").unwrap()).unwrap();
    let cancel = AtomicBool::new(false);
    for (name, source) in [
        ("P1", include_str!("../../../parts/P1-bracket.json")),
        ("P5", include_str!("../../../parts/P5-flanged-housing.json")),
        ("loft-smooth-fillet-r1_10", loft.as_str()),
    ] {
        set_threads(NonZeroUsize::new(1));
        let part = build_incremental_with_files_without_display(
            &Recipe::parse(source).unwrap(),
            None,
            &BTreeMap::new(),
            &mut Decider::default(),
        )
        .unwrap_or_else(|error| panic!("{name}: {error}"));
        let expected = match &part {
            EvaluatedApplicationOutcome::Exact(body) => {
                assert_ne!(name, "loft-smooth-fillet-r1_10");
                let built = body.checkpoint.built();
                built
                    .model
                    .edges_of_body(built.result)
                    .into_iter()
                    .map(EdgeProofKey::Exact)
                    .collect::<Vec<_>>()
            }
            EvaluatedApplicationOutcome::Numerical(body) => {
                assert_ne!(name, "P1");
                body.body()
                    .definition()
                    .edges
                    .iter()
                    .map(|e| EdgeProofKey::Numerical(e.key.clone()))
                    .collect()
            }
        };
        assert!(!expected.is_empty());

        let mut reference = None;
        for threads in [1, 2, 4, 24] {
            set_threads(NonZeroUsize::new(threads));
            let start = Instant::now();
            let report = prove_edge_gaps(&part, &limit, 4096, &cancel);
            let elapsed = start.elapsed();
            assert!(report.edges.iter().map(|edge| &edge.edge).eq(&expected));
            let mut proved = 0;
            for edge in &report.edges {
                if let EdgeProofOutcome::Proved { upper_mm } = &edge.outcome {
                    assert!(upper_mm <= limit.mm());
                    if matches!(edge.edge, EdgeProofKey::Exact(_)) {
                        assert_eq!(upper_mm, &Rational::zero());
                    }
                    proved += 1;
                } else {
                    assert!(
                        matches!(edge.edge, EdgeProofKey::Numerical(_)),
                        "{name}: {edge:?}"
                    );
                }
            }
            println!(
                "{name}, {threads} threads: {proved}/{} proved, {} work, {:.6} s",
                report.edges.len(),
                report.work,
                elapsed.as_secs_f64()
            );
            if let Some(reference) = &reference {
                assert_eq!(&report, reference, "{name}, {threads} threads");
            } else {
                reference = Some(report);
            }
        }
        if name == "P1" {
            let no_work = prove_edge_gaps(&part, &limit, 0, &cancel);
            assert_eq!(no_work.work, 0);
            assert!(
                no_work.edges.iter().all(|edge| edge.outcome
                    == EdgeProofOutcome::CouldNotProve(EdgeProofFailure::WorkCap))
            );
            cancel.store(true, Ordering::Relaxed);
            let cancelled = prove_edge_gaps(&part, &limit, 4096, &cancel);
            assert_eq!(cancelled.work, 0);
            assert!(
                cancelled.edges.iter().all(|edge| edge.outcome
                    == EdgeProofOutcome::CouldNotProve(EdgeProofFailure::Cancelled))
            );
            cancel.store(false, Ordering::Relaxed);
        }
    }
    set_threads(None);
}
