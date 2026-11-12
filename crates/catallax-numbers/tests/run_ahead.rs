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
                    if charges.adopt(&dec) {
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

// A reserved allowance run inside recorded code starts a counter of its own
// and sets the proof caches aside, which no record can describe. The helper
// stays a helper inside it, and what it was recording is not adopted.
#[test]
fn a_nested_allowance_inside_a_record_stays_isolated_and_is_not_adopted() {
    let _setting = SETTING.lock().unwrap_or_else(|e| e.into_inner());
    set_threads(NonZeroUsize::new(2));
    let ((seen, usable, shared_inside, shared_after), work) =
        with_operation_budget(u64::MAX, || {
            let dec = Decider::default();
            let helpers = Helpers::new(&dec).expect("two threads");
            std::thread::scope(|scope| {
                scope
                    .spawn(|| {
                        helpers.run(|dec| {
                            let mut shared_inside = None;
                            let (seen, charges) = record(|| {
                                let first = show(Ball::pi_checked(96, dec));
                                let nested = OperationAllowance::reserve(100_000).unwrap();
                                let (second, inner) = nested.run(|| {
                                    shared_inside = Some(shared_caches());
                                    show(Ball::pi_checked(160, dec))
                                });
                                assert!(inner.used > 0 && !inner.exhausted);
                                // Charged after the nested counter is gone.
                                let third = show(Ball::pi_checked(224, dec));
                                (first, second, third)
                            });
                            (seen, charges.is_some(), shared_inside, shared_caches())
                        })
                    })
                    .join()
                    .unwrap()
            })
        });
    set_threads(None);
    assert!(seen.0.starts_with("Ok") && seen.1.starts_with("Ok") && seen.2.starts_with("Ok"));
    assert_eq!(
        shared_inside,
        Some(false),
        "a helper inside a nested allowance"
    );
    assert!(!shared_after);
    assert!(!usable, "a record across a nested allowance");
    assert_eq!(work.used, 0, "nothing a helper charges reaches the writer");
}

// A helper works to the writer's deadline, and nothing is adopted once it has
// passed: the same computation run by the writer would refuse.
#[test]
fn helpers_and_adoption_keep_the_writers_deadline() {
    let _setting = SETTING.lock().unwrap_or_else(|e| e.into_inner());
    set_threads(NonZeroUsize::new(2));
    let (seen, work) = with_operation_budget(u64::MAX, || {
        let mut late = Decider::default();
        late.set_deadline(Some(std::time::Instant::now()));
        std::thread::sleep(std::time::Duration::from_millis(5));
        assert!(late.check_deadline().is_err());
        // A helper of a writer that is already late refuses as it would.
        let helpers = Helpers::new(&late).expect("two threads");
        let helper_is_late = std::thread::scope(|scope| {
            scope
                .spawn(|| helpers.run(|dec| dec.check_deadline().is_err()))
                .join()
                .unwrap()
        });
        // A record made in time is adopted in time and not afterwards.
        let timely = Decider::default();
        let helpers = Helpers::new(&timely).expect("two threads");
        let (value, charges) = std::thread::scope(|scope| {
            scope
                .spawn(|| helpers.run(|dec| record(|| show(Ball::pi_checked(200, dec)))))
                .join()
                .unwrap()
        });
        let charges = charges.expect("a usable record");
        // The writer's own deadline check above cost it a unit already.
        let before = current_operation_work().used;
        let refused = !charges.adopt(&late);
        let charged_when_late = current_operation_work().used - before;
        let adopted = charges.adopt(&timely);
        let charged_in_time = current_operation_work().used - before;
        (
            (helper_is_late, refused, charged_when_late, adopted, value),
            (before, charged_in_time),
        )
    });
    set_threads(None);
    let ((helper_is_late, refused, charged_when_late, adopted, value), (before, charged_in_time)) =
        seen;
    assert!(helper_is_late, "a helper has the writer's deadline");
    assert!(refused, "nothing is adopted past the deadline");
    assert_eq!(charged_when_late, 0);
    assert!(adopted && charged_in_time > 0 && value.starts_with("Ok"));
    assert_eq!(work.used, before + charged_in_time);
}

// A whole loop item: a value computed from shared inputs alone. Every third
// one isolates roots, which no record can describe, so the writer computes it.
fn whole(shared: &[Real], index: usize, dec: &mut Decider) -> Vec<String> {
    let [two, three, five, alpha, cube] = shared else {
        panic!("five shared values");
    };
    let n = Real::from_i64(index as i64 + 2);
    let six = (two * three).square().sqrt();
    let zero = &(two + three).square() - &(&Real::from_i64(5) + &(&Real::from_i64(2) * &six));
    let big = Real::from(r("1e40")) * n.clone();
    let tie = &(&big.square() + &Real::one()).sqrt() - &big;
    let width = r("1e-60");
    let mut seen = vec![
        show(dec.sign(&(&zero * &n))),
        show(dec.sign(&(&tie - &(Real::one() / (&big * &Real::from_i64(2)))))),
        show(dec.compare(&(alpha * &n), &(cube * five))),
        show(dec.rational_value(&(&(two * two) * &n))),
        show(dec.enclose(&(&(alpha * two) + &n.sqrt()), &width)),
        show(dec.enclose_sign(&(&(cube * three) - &n.sqrt()), &width)),
        show(dec.enclose_at(&(&(cube * alpha) - &n), 300)),
        show(Ball::from_rational(&r("17/3"), 160 + 32 * (index as u32 % 3)).atan_checked(dec)),
    ];
    dec.note_fact_decision();
    if index % 3 == 2 {
        seen.push(show(
            dec.real_roots(&poly(&[-(index as i64) - 3, 1, 0, 0, 1]))
                .map(|roots| roots.len()),
        ));
    }
    seen
}

fn run_whole(threads: usize, wait: bool, limit: u64) -> Run {
    set_threads(NonZeroUsize::new(threads));
    set_run_ahead_wait(wait);
    let ((seen, stats), work) = with_operation_budget(limit, || {
        let mut dec = Decider::default();
        let values = shared(&mut dec);
        static SITE: AheadSite = AheadSite::new();
        let mut ahead = ValuesAhead::new(&SITE, &dec, ITEMS);
        let mut seen = Vec::new();
        for index in 0..ITEMS {
            let adopted = ahead.take(index, &mut dec, || {
                let values = values.clone();
                move |index: usize, dec: &mut Decider| whole(&values, index, dec)
            });
            seen.push(adopted.unwrap_or_else(|| whole(&values, index, &mut dec)));
        }
        drop(ahead);
        (seen, dec.stats().clone())
    });
    set_threads(None);
    set_run_ahead_wait(false);
    (seen, stats, work.used, work.exhausted)
}

#[test]
fn whole_items_from_other_threads_change_nothing() {
    let _setting = SETTING.lock().unwrap_or_else(|e| e.into_inner());
    let one = run_whole(1, false, u64::MAX);
    assert!(!one.3);
    assert!(one.2 > 50_000, "the loop does real work: {}", one.2);
    let before = run_ahead_adopted();
    assert_eq!(run_whole(2, true, u64::MAX), one, "two threads, waiting");
    // Items 1, 3 and 4 can be adopted; 2 and 5 isolate roots.
    assert_eq!(run_ahead_adopted() - before, 3);
    assert_eq!(run_whole(4, true, u64::MAX), one, "four threads, waiting");
    for threads in [2, 4, 8] {
        for _ in 0..3 {
            assert_eq!(
                run_whole(threads, false, u64::MAX),
                one,
                "{threads} threads"
            );
        }
    }
    let total = one.2;
    let mut limits = vec![total, total - 1, total + 1];
    limits.extend((1..10).map(|step| total / 10 * step + step * 3));
    for limit in limits {
        let one = run_whole(1, false, limit);
        assert_eq!(one.3, limit < total, "limit {limit} of {total}");
        assert_eq!(run_whole(4, true, limit), one, "limit {limit}, waiting");
        assert_eq!(run_whole(3, false, limit), one, "limit {limit}");
    }
}

// The pool keeps the workers it once had. After the thread count drops, only
// the first few may take work, and a loop that waits for its workers must
// still get one of those.
#[test]
fn a_pool_larger_than_the_thread_count_still_serves_a_waiting_loop() {
    use std::sync::mpsc;
    use std::time::Duration;
    let _setting = SETTING.lock().unwrap_or_else(|e| e.into_inner());
    // Grow the pool to seven workers.
    let one = run_whole(1, false, u64::MAX);
    assert_eq!(run_whole(8, true, u64::MAX), one);
    // Two threads: one worker may work. A loop of two items wakes one.
    for attempt in 0..20 {
        let (done, finished) = mpsc::channel();
        std::thread::spawn(move || {
            set_threads(NonZeroUsize::new(2));
            set_run_ahead_wait(true);
            static SITE: AheadSite = AheadSite::new();
            let mut dec = Decider::default();
            let mut ahead = ValuesAhead::<usize>::new(&SITE, &dec, 2);
            let taken: Vec<_> = (0..2)
                .map(|index| ahead.take(index, &mut dec, || |index: usize, _: &mut Decider| index))
                .collect();
            drop(ahead);
            let _ = done.send(taken);
        });
        let taken = finished.recv_timeout(Duration::from_secs(20));
        set_threads(None);
        set_run_ahead_wait(false);
        let taken = taken.unwrap_or_else(|_| {
            panic!("attempt {attempt}: the loop waited for a worker that was never woken")
        });
        assert_eq!(taken[0], None, "the writer computes the first item");
    }
}

// A pool worker works to the writer's deadline, and the writer waits for a
// worker's item no longer than its own deadline.
#[test]
fn pool_workers_keep_the_writers_deadline() {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::time::{Duration, Instant};
    let _setting = SETTING.lock().unwrap_or_else(|e| e.into_inner());
    set_threads(NonZeroUsize::new(2));

    // A writer that is already late: its worker's decider refuses too.
    set_run_ahead_wait(true);
    let worker_was_late = Arc::new(AtomicBool::new(false));
    let mut late = Decider::default();
    late.set_deadline(Some(Instant::now()));
    std::thread::sleep(Duration::from_millis(5));
    {
        static SITE: AheadSite = AheadSite::new();
        let mut ahead = ValuesAhead::<()>::new(&SITE, &late, 3);
        let flag = worker_was_late.clone();
        let adopted = ahead.take(0, &mut late, || {
            move |_: usize, dec: &mut Decider| {
                if dec.check_deadline().is_err() {
                    flag.store(true, Ordering::Relaxed);
                }
            }
        });
        assert!(adopted.is_none(), "the writer computes the first item");
    }
    set_run_ahead_wait(false);
    assert!(
        worker_was_late.load(Ordering::Relaxed),
        "a pool worker has the writer's deadline"
    );

    // A worker has had an item for a second and a half when the writer
    // reaches it, 200 milliseconds before the writer's deadline. The writer
    // gives the worker until the deadline, not another second and a half.
    fn item(index: usize) {
        let on_worker = std::thread::current()
            .name()
            .is_some_and(|name| name.starts_with("catallax-ahead-"));
        std::thread::sleep(Duration::from_millis(if on_worker {
            4000
        } else if index == 1 {
            1500
        } else {
            5
        }));
    }
    static SLOW: AheadSite = AheadSite::new();
    let began = Instant::now();
    let mut dec = Decider::default();
    dec.set_deadline(Some(began + Duration::from_millis(1700)));
    let mut ahead = ValuesAhead::<()>::new(&SLOW, &dec, 3);
    let mut waited = Duration::ZERO;
    for index in 0..3 {
        let asked = Instant::now();
        let adopted = ahead.take(index, &mut dec, || {
            |index: usize, _: &mut Decider| item(index)
        });
        if index == 2 {
            waited = asked.elapsed();
        }
        if adopted.is_none() {
            item(index);
        }
    }
    drop(ahead);
    set_threads(None);
    assert!(
        waited < Duration::from_millis(900),
        "the writer waited {waited:?} for a worker, past its deadline"
    );
}
