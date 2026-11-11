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
}

struct Board<T> {
    slots: Vec<Mutex<Slot<T>>>,
    changed: Condvar,
    next: AtomicUsize,
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
            let spawned = std::thread::Builder::new()
                .name("catallax-assembly".into())
                .stack_size(helpers.stack_size())
                .spawn_scoped(scope, || {
                    helpers.run(|worker| explore(items, &board, &helpers, evaluate, worker));
                });
            if spawned.is_err() {
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
        let (value, charges) = catallax_numbers::record(|| evaluate(&items[index], worker));
        *locked(&board.slots[index]) = Slot::Done(Box::new(value), charges);
        board.changed.notify_all();
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
            Slot::Open if !wait => {
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
