//! Independent box-pair proofs run ahead on every core. The request's own
//! thread takes each result in order and adopts the helper's recorded charges,
//! so the answer, a refusal and the operation work are those of one thread.
//! With one thread, or when a record cannot describe the work, the request's
//! thread evaluates the item itself.

use catallax_numbers::{Charges, Decider, Helpers};
use std::ops::ControlFlow;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex, MutexGuard, PoisonError};

enum Slot<T> {
    Open,
    Taken,
    Writer,
    Done(Box<T>, Option<Charges>),
    // A helper unwound before it settled the item; the writer evaluates it.
    Abandoned,
}

struct Board<T> {
    slots: Vec<Mutex<Slot<T>>>,
    changed: Condvar,
    next: AtomicUsize,
    // Helpers still taking items. A waiting writer stops waiting for an open
    // item once none is left to take it.
    exploring: AtomicUsize,
}

fn locked<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(PoisonError::into_inner)
}

/// Evaluates `items` and hands each value to `consume` in index order. Stops at
/// the first `Break`; helpers are told to stop and their remaining work is
/// discarded without being charged.
pub(super) fn in_order<I, T, R>(
    items: &[I],
    evaluate: impl Fn(&I, &mut Decider) -> T + Sync,
    dec: &mut Decider,
    mut consume: impl FnMut(usize, T, &mut Decider) -> ControlFlow<R>,
) -> ControlFlow<R>
where
    I: Sync,
    T: Send,
{
    let helpers = if items.len() > 1 {
        Helpers::new(dec)
    } else {
        None
    };
    let Some(helpers) = helpers else {
        for (index, item) in items.iter().enumerate() {
            let value = evaluate(item, dec);
            consume(index, value, dec)?;
        }
        return ControlFlow::Continue(());
    };
    let board = Board {
        slots: (0..items.len()).map(|_| Mutex::new(Slot::Open)).collect(),
        changed: Condvar::new(),
        next: AtomicUsize::new(0),
        exploring: AtomicUsize::new(0),
    };
    let wait = catallax_numbers::run_ahead_waits();
    let evaluate = &evaluate;
    std::thread::scope(|scope| {
        // The scope joins the helpers on every path out, so they are told to
        // stop on every path out too.
        struct Stop<'h>(&'h Helpers);
        impl Drop for Stop<'_> {
            fn drop(&mut self) {
                self.0.stop();
            }
        }
        let _stop = Stop(&helpers);
        let mut started = 0;
        for _ in 0..helpers.count().min(items.len()) {
            board.exploring.fetch_add(1, Ordering::SeqCst);
            let spawned = std::thread::Builder::new()
                .name("catallax-assembly".into())
                .stack_size(helpers.stack_size())
                .spawn_scoped(scope, || {
                    helpers.run(|worker| explore(items, &board, &helpers, evaluate, worker));
                });
            if spawned.is_err() {
                board.exploring.fetch_sub(1, Ordering::SeqCst);
                break;
            }
            started += 1;
        }
        for index in 0..items.len() {
            let value = take(items, &board, index, wait && started > 0, evaluate, dec);
            consume(index, value, dec)?;
        }
        ControlFlow::Continue(())
    })
}

fn explore<I, T>(
    items: &[I],
    board: &Board<T>,
    helpers: &Helpers,
    evaluate: &(impl Fn(&I, &mut Decider) -> T + Sync),
    worker: &mut Decider,
) {
    // Every path out, a panic included, leaves the board. Taking each slot's
    // lock before the notice reaches a writer that is about to wait on one.
    struct Leave<'b, T>(&'b Board<T>);
    impl<T> Drop for Leave<'_, T> {
        fn drop(&mut self) {
            self.0.exploring.fetch_sub(1, Ordering::SeqCst);
            for slot in &self.0.slots {
                drop(locked(slot));
            }
            self.0.changed.notify_all();
        }
    }
    // The item is settled on every path out, a panic included: a slot left
    // taken would be waited for.
    struct Settle<'b, T>(&'b Board<T>, usize);
    impl<T> Drop for Settle<'_, T> {
        fn drop(&mut self) {
            let mut slot = locked(&self.0.slots[self.1]);
            if matches!(*slot, Slot::Taken) {
                *slot = Slot::Abandoned;
            }
            drop(slot);
            self.0.changed.notify_all();
        }
    }
    let _leave = Leave(board);
    loop {
        if helpers.stopped() {
            return;
        }
        let index = board.next.fetch_add(1, Ordering::Relaxed);
        if index >= items.len() {
            return;
        }
        {
            let mut slot = locked(&board.slots[index]);
            if !matches!(*slot, Slot::Open) {
                continue;
            }
            *slot = Slot::Taken;
        }
        let settle = Settle(board, index);
        let (value, charges) = catallax_numbers::record(|| evaluate(&items[index], worker));
        *locked(&board.slots[index]) = Slot::Done(Box::new(value), charges);
        drop(settle);
    }
}

/// The writer's value for one item: a helper's, adopted with its charges, or
/// its own.
fn take<I, T>(
    items: &[I],
    board: &Board<T>,
    index: usize,
    wait: bool,
    evaluate: &impl Fn(&I, &mut Decider) -> T,
    dec: &mut Decider,
) -> T {
    let mut slot = locked(&board.slots[index]);
    loop {
        match &*slot {
            Slot::Open if !wait || board.exploring.load(Ordering::SeqCst) == 0 => {
                *slot = Slot::Writer;
                drop(slot);
                return evaluate(&items[index], dec);
            }
            Slot::Abandoned => {
                *slot = Slot::Writer;
                drop(slot);
                return evaluate(&items[index], dec);
            }
            Slot::Open | Slot::Taken => {
                slot = board
                    .changed
                    .wait(slot)
                    .unwrap_or_else(PoisonError::into_inner);
            }
            Slot::Writer => unreachable!("the writer takes each item once"),
            Slot::Done(..) => break,
        }
    }
    let Slot::Done(value, charges) = std::mem::replace(&mut *slot, Slot::Writer) else {
        unreachable!("matched above");
    };
    drop(slot);
    if let Some(charges) = charges
        && charges.adopt(dec)
    {
        return *value;
    }
    evaluate(&items[index], dec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::num::NonZeroUsize;
    use std::sync::atomic::AtomicBool;
    use std::sync::{Arc, mpsc};
    use std::time::Duration;

    // Runs sixteen items on four threads whose helpers unwind on every item
    // they take; returns whether the call panicked once the scope joined.
    fn unwinding_helpers() -> bool {
        let unwound = Arc::new(AtomicBool::new(false));
        let consumed = Arc::new(Mutex::new(Vec::new()));
        let (done, finished) = mpsc::channel();
        {
            let (unwound, consumed) = (unwound.clone(), consumed.clone());
            std::thread::spawn(move || {
                let items: Vec<usize> = (0..16).collect();
                let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    in_order(
                        &items,
                        |&item, _| {
                            if std::thread::current().name() == Some("catallax-assembly") {
                                unwound.store(true, Ordering::SeqCst);
                                panic!("a helper unwinds on item {item}");
                            }
                            // The writer holds its first item until a helper has
                            // taken one and unwound, so a later item waits on it.
                            let start = std::time::Instant::now();
                            while item == 0
                                && !unwound.load(Ordering::SeqCst)
                                && start.elapsed() < Duration::from_secs(30)
                            {
                                std::thread::yield_now();
                            }
                            item
                        },
                        &mut Decider::default(),
                        |_, value, _| {
                            locked(&consumed).push(value);
                            ControlFlow::<()>::Continue(())
                        },
                    )
                }));
                let _ = done.send(outcome.is_err());
            });
        }
        let panicked = finished
            .recv_timeout(Duration::from_secs(120))
            .expect("the writer waited for an item no helper would settle");
        assert!(unwound.load(Ordering::SeqCst), "no helper took an item");
        // Every item reaches the writer in order.
        assert_eq!(*locked(&consumed), (0..16).collect::<Vec<_>>());
        panicked
    }

    #[test]
    fn a_helper_that_unwinds_leaves_its_items_to_the_writer() {
        struct Restore;
        impl Drop for Restore {
            fn drop(&mut self) {
                catallax_numbers::set_threads(None);
                catallax_numbers::set_run_ahead_wait(false);
            }
        }
        let _restore = Restore;
        catallax_numbers::set_threads(NonZeroUsize::new(4));
        // The helpers' panics still reach the caller when the scope joins.
        assert!(unwinding_helpers());
        // A writer that waits for its helpers stops waiting once none is left.
        catallax_numbers::set_run_ahead_wait(true);
        assert!(unwinding_helpers());
    }
}
