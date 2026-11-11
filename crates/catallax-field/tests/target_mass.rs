//! Exact signed moments, worker records, refusals and independent B-rep agreement.
use catallax_field::{solver::GridOptions, target_mass};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn box_recipe() -> Recipe {
    Recipe::parse(r#"{"format":"catallax-recipe/3","name":"signed box","steps":[{"id":"b","op":"block","size":["2","4","6"],"pose":{"rotation":["1","0","0","0"],"translation":["-3","-1","2"]}}],"result":"b"}"#).unwrap()
}
#[test]
fn full_records_and_refusals_are_worker_independent() {
    let recipe = box_recipe();
    for cap in [100_000_000, 100, 1] {
        for work_per_cell in [256, 1] {
            let mut prior = None;
            for workers in [Some(1), Some(2), Some(3), Some(4), Some(99), None] {
                let (run, work) = with_operation_budget(cap, || {
                    target_mass::measure(
                        &recipe,
                        &q("2"),
                        &GridOptions {
                            divisions: 9,
                            workers,
                            work_per_cell,
                        },
                        &mut Decider::default(),
                    )
                });
                let record = (run.answer, work.used, work.exhausted, run.classifier_work);
                if cap == 100_000_000 && work_per_cell == 256 {
                    let result = record.0.as_ref().unwrap();
                    assert_eq!(result.volume.lo, q("48"));
                    assert_eq!(result.volume.hi, q("48"));
                    assert_eq!(result.grams.lo, q("96"));
                    let c = result.centre.as_ref().unwrap();
                    for (v, expected) in c.iter().zip(["-2", "1", "5"]) {
                        assert_eq!(v.lo, q(expected));
                        assert_eq!(v.hi, q(expected));
                    }
                    assert!(!record.2);
                    let expected_workers = workers
                        .unwrap_or_else(|| {
                            std::thread::available_parallelism().map_or(1, usize::from)
                        })
                        .min(6);
                    assert_eq!(run.spawned_workers, expected_workers);
                } else {
                    assert!(record.0.is_err());
                }
                if let Some(prior) = &prior {
                    assert_eq!(&record, prior);
                }
                prior = Some(record);
            }
        }
    }
}
#[test]
fn density_is_explicit_and_positive() {
    for density in ["0", "-1"] {
        let (run, _) = with_operation_budget(1_000_000, || {
            target_mass::measure(
                &box_recipe(),
                &q(density),
                &GridOptions::default(),
                &mut Decider::default(),
            )
        });
        assert_eq!(run.answer.unwrap_err().code, "target-mass.density");
        assert_eq!(run.spawned_workers, 0);
    }
}

#[test]
fn independent_box_agreement_returns_one_exact_result() {
    let (run, work) = with_operation_budget(1_000_000, || {
        target_mass::measure_agreed(
            &box_recipe(),
            &q("2"),
            &GridOptions {
                divisions: 4,
                ..Default::default()
            },
            &Default::default(),
            &mut Decider::default(),
        )
    });
    let result = run.answer.unwrap();
    assert!(result.brep_agreed);
    assert_eq!(result.grams.lo, q("96"));
    assert_eq!(result.grams.hi, q("96"));
    assert!(!work.exhausted);
}

#[test]
fn translated_sphere_encloses_analytic_mass_and_centre() {
    let recipe = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"sphere","steps":[{"id":"s","op":"sphere","radius":"2","pose":{"rotation":["1","0","0","0"],"translation":["-5","2","-1"]}}],"result":"s"}"#).unwrap();
    let (run, work) = with_operation_budget(20_000_000, || {
        target_mass::measure(
            &recipe,
            &q("3"),
            &GridOptions {
                divisions: 16,
                ..Default::default()
            },
            &mut Decider::default(),
        )
    });
    let result = run.answer.unwrap();
    // 3.14 < pi < 22/7: both bounds of the analytic volume must be covered.
    assert!(result.volume.lo < q("32/3") * q("3.14"));
    assert!(result.volume.hi > q("32/3") * q("22/7"));
    for (c, expected) in result.centre.unwrap().iter().zip(["-5", "2", "-1"]) {
        assert!(c.contains(&q(expected)));
    }
    assert!(!work.exhausted);
}

#[test]
fn expired_deadline_never_returns_a_mass() {
    let mut dec = Decider::default();
    dec.set_deadline(Some(std::time::Instant::now()));
    let (run, _) = with_operation_budget(1_000_000, || {
        target_mass::measure(&box_recipe(), &q("1"), &GridOptions::default(), &mut dec)
    });
    assert!(run.answer.is_err());
}
