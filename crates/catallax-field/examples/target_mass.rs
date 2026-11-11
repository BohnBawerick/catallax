//! Complete public target evidence; no prefix recipes or saved output inputs.
use catallax_field::{solver::GridOptions, target_mass};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
fn main() {
    let divisions = std::env::args().nth(1).map_or(16, |s| s.parse().unwrap());
    let density = Rational::parse("0.00785").unwrap();
    let agreed = std::env::args().nth(2).is_some_and(|s| s == "agreed");
    for name in [
        "P1-bracket",
        "P2-collar",
        "P3-filleted-boss",
        "P4-shelled-tray",
        "P5-flanged-housing",
        "P6-lattice-bracket",
    ] {
        let text = std::fs::read_to_string(format!("parts/{name}.json")).unwrap();
        let recipe = Recipe::parse(&text).unwrap();
        let mut prior = None;
        for workers in [Some(1), Some(2), Some(4), None] {
            let (run, work) = with_operation_budget(100_000_000, || {
                let options = GridOptions {
                    divisions,
                    workers,
                    ..Default::default()
                };
                if agreed {
                    return target_mass::measure_agreed(
                        &recipe,
                        &density,
                        &options,
                        &Default::default(),
                        &mut Decider::default(),
                    );
                }
                target_mass::measure(
                    &recipe,
                    &density,
                    &GridOptions {
                        divisions,
                        workers,
                        ..Default::default()
                    },
                    &mut Decider::default(),
                )
            });
            println!(
                "{name} divisions={divisions} workers={workers:?} spawned={} charged={} exhausted={} classifier_work={}",
                run.spawned_workers, work.used, work.exhausted, run.classifier_work
            );
            let record = (run.answer, work.used, work.exhausted, run.classifier_work);
            if let Some(prior) = &prior {
                assert_eq!(&record, prior);
            } else {
                println!("{name} result={:?}", record.0);
            }
            prior = Some(record);
        }
    }
}
