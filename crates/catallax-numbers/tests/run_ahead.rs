#![allow(missing_docs)]
//! Workers must change nothing: the outcomes, the statistics and the request's
//! work are those of one thread, with every record adopted or none.
mod common;
use catallax_numbers::*;
use common::*;
use std::{num::NonZeroUsize, sync::Mutex};

// The thread count is one setting for the process.
static SETTING: Mutex<()> = Mutex::new(());

const ITEMS: usize = 6;

fn show<T: std::fmt::Debug>(value: T) -> String {
    format!("{value:?}")
}

// Values made before the loop: every worker and the writer share these nodes.
fn shared(dec: &mut Decider) -> Vec<Real> {
    let mut values = vec![real("2").sqrt(), real("3").sqrt(), real("5").sqrt()];
    // x^5 - x - 1 has one real root, and x^3 - 2 another.
    for coefficients in [&[-1, -1, 0, 0, 0, 1][..], &[-2, 0, 0, 1]] {
        for root in roots(dec, &poly(coefficients)) {
            values.push(Real::from_algebraic(root.value));
        }
    }
    // Certificates the writer holds before the loop: a worker starts from
    // them, and a record that reused one replays as a reuse.
    for width in ["1e-30", "1e-90"] {
        // Under a small allowance this may refuse, like anything later.
        let _ = dec.enclose(&(&values[3] * &values[4]), &r(width));
    }
    values
}

// One loop item: decisions that need refinement, on shared and on fresh nodes.
fn item(shared: &[Real], index: usize, dec: &mut Decider) -> Vec<String> {
    let mut seen = Vec::new();
    let n = Real::from_i64(index as i64 + 2);
    let [two, three, five, alpha, cube] = shared else {
        panic!("five shared values");
    };
    // An exact zero: (sqrt 2 + sqrt 3)^2 = 5 + 2 sqrt 6.
    let six = (two * three).square().sqrt();
    let zero = &(two + three).square() - &(&Real::from_i64(5) + &(&Real::from_i64(2) * &six));
    seen.push(show(dec.sign(&(&zero * &n))));
    // A near tie that needs a few hundred bits.
    let big = Real::from(r("1e40")) * n.clone();
    let tie = &(&big.square() + &Real::one()).sqrt() - &big;
    seen.push(show(
        dec.sign(&(&tie - &(Real::one() / (&big * &Real::from_i64(2))))),
    ));
    seen.push(show(dec.compare(&(alpha * &n), &(cube * five))));
    // The root's own equation, and a product that is exactly rational.
    let equation = &(&alpha.pow_u32(5) - alpha) - &Real::one();
    seen.push(show(dec.sign(&(&equation + &(&zero * &n)))));
    seen.push(show(dec.rational_value(&(&(two * two) * &n))));
    seen.push(show(dec.rational_identity(&(&(cube * cube) * cube))));
    let width = r("1e-60");
    seen.push(show(dec.enclose(&(&(alpha * two) + &n.sqrt()), &width)));
    // The same question to the same width, stopping at the first sign: a
    // different decision, with a different ball and less work.
    seen.push(show(
        dec.enclose_sign(&(&(alpha * two) + &n.sqrt()), &width),
    ));
    seen.push(show(dec.enclose_sign(&(&zero * &n), &width)));
    seen.push(show(dec.enclose_at(&(&(cube * three) - &n.sqrt()), 300)));
    let both = [&(alpha + &n.sqrt()), &(alpha * &n.sqrt())];
    seen.push(show(dec.enclose_all_at(&both, 200)));
    seen.push(show(
        dec.enclose_budgeted_at(&(&(five * alpha) / &n.sqrt()), 400),
    ));
    // Fresh roots, isolated inside the item, then used.
    let shift = index as i64 + 3;
    let found = dec.real_roots(&poly(&[-shift, 1, 0, 0, 1]));
    if let Outcome::Proven(found) = &found {
        for root in found {
            let value = Real::from_algebraic(root.value.clone());
            seen.push(show(dec.sign(&(&(&value * two) - &n))));
            seen.push(show(dec.enclose(&(&value + alpha), &width)));
        }
    }
    seen.push(show(found.map(|roots| roots.len())));
    seen
}

type Run = (Vec<Vec<String>>, DecisionStats, u64, bool);

fn run(threads: usize, wait: bool, limit: u64) -> Run {
    set_threads(NonZeroUsize::new(threads));
    set_run_ahead_loops(true);
    set_run_ahead_wait(wait);
    let ((seen, stats), work) = with_operation_budget(limit, || {
        let mut dec = Decider::default();
        let values = shared(&mut dec);
        static SITE: AheadSite = AheadSite::new();
        let mut ahead = RunAhead::new(&SITE, &dec, ITEMS);
        let mut seen = Vec::new();
        for index in 0..ITEMS {
            ahead.at(index, || {
                let values = values.clone();
                move |index: usize, dec: &mut Decider| {
                    item(&values, index, dec);
                }
            });
            seen.push(item(&values, index, &mut dec));
        }
        drop(ahead);
        (seen, dec.stats().clone())
    });
    set_threads(None);
    set_run_ahead_loops(false);
    set_run_ahead_wait(false);
    (seen, stats, work.used, work.exhausted)
}

#[test]
fn workers_change_no_outcome_statistic_or_work() {
    let _setting = SETTING.lock().unwrap_or_else(|e| e.into_inner());
    let one = run(1, false, u64::MAX);
    assert!(!one.3);
    assert!(one.2 > 100_000, "the loop does real work: {}", one.2);
    let before = run_ahead_adopted();
    assert_eq!(
        run(2, true, u64::MAX),
        one,
        "two threads, every record adopted"
    );
    let adopted = run_ahead_adopted() - before;
    assert!(adopted > 20, "records were adopted: {adopted}");
    assert_eq!(
        run(4, true, u64::MAX),
        one,
        "four threads, every record adopted"
    );
    for threads in [2, 4, 8] {
        for _ in 0..3 {
            assert_eq!(run(threads, false, u64::MAX), one, "{threads} threads");
        }
    }
}

#[test]
fn the_allowance_runs_out_at_the_same_charge() {
    let _setting = SETTING.lock().unwrap_or_else(|e| e.into_inner());
    let total = run(1, false, u64::MAX).2;
    let mut limits = vec![total, total - 1, total + 1];
    limits.extend((1..8).map(|step| total / 8 * step + step * 7));
    for limit in limits {
        let one = run(1, false, limit);
        assert_eq!(one.3, limit < total, "limit {limit} of {total}");
        assert_eq!(
            run(4, true, limit),
            one,
            "limit {limit}, every record adopted"
        );
        assert_eq!(run(4, false, limit), one, "limit {limit}");
    }
}

// A computation outside any single decision: checked arctangents, which reuse
// the request's pi enclosures, and enclosures over shared roots.
fn computation(shared: &[Real], index: usize, dec: &Decider) -> Vec<String> {
    let [two, _, _, alpha, cube] = shared else {
        panic!("five shared values");
    };
    let n = Real::from_i64(index as i64 + 2);
    let steep = Ball::from_rational(&r("17/3"), 160 + 32 * (index as u32 % 3));
    vec![
        show(steep.atan_checked(dec)),
        show(Ball::from_rational(&r("1/7"), 200).atan_checked(dec)),
        show(dec.enclose_at(&(&(alpha * two) + &n.sqrt()), 300)),
        show(dec.enclose_at(&(&(cube * alpha) - &n), 200)),
        show(dec.enclose_all_at(&[alpha, &(cube * &n)], 250)),
    ]
}

#[test]
fn recorded_computations_charge_what_cold_ones_do() {
    let _setting = SETTING.lock().unwrap_or_else(|e| e.into_inner());
    const COUNT: usize = 9;
    let alone = |limit: u64| {
        set_threads(NonZeroUsize::new(1));
        let (seen, work) = with_operation_budget(limit, || {
            let mut dec = Decider::default();
            let values = shared(&mut dec);
            (0..COUNT)
                .map(|index| computation(&values, index, &dec))
                .collect::<Vec<_>>()
        });
        set_threads(None);
        (seen, work.used, work.exhausted)
    };
    let helped = |limit: u64| {
        set_threads(NonZeroUsize::new(3));
        let (seen, work) = with_operation_budget(limit, || {
            let mut dec = Decider::default();
            let values = shared(&mut dec);
            let helpers = Helpers::new(&dec).expect("three threads");
            // Every computation twice: adopting a record again must charge
            // what a second cold run charges.
            let recorded = std::thread::scope(|scope| {
                scope
                    .spawn(|| {
                        helpers.run(|dec| {
                            (0..COUNT)
                                .map(|index| record(|| computation(&values, index, dec)))
                                .collect::<Vec<_>>()
                        })
                    })
                    .join()
                    .unwrap()
            });
            recorded
                .iter()
                .enumerate()
                .map(|(index, (value, charges))| {
                    let charges = charges.as_ref().expect("a usable record");
                    if charges.adopt() {
                        value.clone()
                    } else {
                        computation(&values, index, &dec)
                    }
                })
                .collect::<Vec<_>>()
        });
        set_threads(None);
        (seen, work.used, work.exhausted)
    };
    let cold = alone(u64::MAX);
    assert!(cold.1 > 10_000, "the computations do real work: {}", cold.1);
    let before = run_ahead_adopted();
    assert_eq!(helped(u64::MAX), cold);
    assert_eq!(run_ahead_adopted() - before, COUNT as u64);
    for step in 1..12 {
        let limit = cold.1 / 12 * step + step;
        assert_eq!(helped(limit), alone(limit), "limit {limit} of {}", cold.1);
    }
}
