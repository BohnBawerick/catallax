//! Measure the complete two-way check, including both builds and display meshing.
use catallax_compare::{Completion, Controller, Options, Reason};
use catallax_recipe::Recipe;
use std::{collections::BTreeMap, time::Instant};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() {
        return Err("usage: compare recipe.json [runs]".into());
    }
    let recipe = Recipe::parse(&std::fs::read_to_string(&args[0])?)?;
    let runs: usize = args.get(1).map(|n| n.parse()).transpose()?.unwrap_or(1);
    let controller = Controller::default();
    for run in 0..runs {
        let revision = controller.begin();
        let start = Instant::now();
        assert_eq!(
            revision.check(&recipe, &Options::default()),
            Completion::Published
        );
        let elapsed = start.elapsed();
        let report = controller.latest().unwrap();
        println!(
            "run={} elapsed_ms={:.3} verdict={:?} forward_vertices={} reverse_samples={} findings={}",
            run + 1,
            elapsed.as_secs_f64() * 1000.,
            report.verdict,
            report.forward_vertices,
            report.reverse_samples,
            report.findings.len()
        );
        let mut reasons = BTreeMap::new();
        for finding in &report.findings {
            let key = match &finding.reason {
                Reason::Unresolved(s) => format!("unresolved: {s}"),
                reason => format!("{reason:?}"),
            };
            *reasons.entry(key).or_insert(0usize) += 1;
        }
        for (reason, count) in reasons {
            println!("  {count} {reason}");
        }
        for op in &report.not_comparable {
            println!(
                "  step={} op={} reason={}",
                op.step, op.operation, op.reason
            );
        }
        if let Some(f) = report
            .findings
            .iter()
            .find(|f| !matches!(f.reason, Reason::Unresolved(_)))
        {
            println!("  first_disagreement={f:?}");
        }
    }
    Ok(())
}
