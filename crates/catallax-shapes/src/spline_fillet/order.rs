//! Independent tasks proved on every core, adopted in task order.
//!
//! Helper threads take tasks in order and record their charges. The request's
//! own thread adopts each task's record in order, or proves the task itself
//! when no helper has finished it, so results, refusals and work are those of
//! one thread.

use catallax_numbers::{Charges, Decider, Helpers, Refusal};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};

enum Slot<T> {
    Open,
    Taken,
    Writer,
    Done(Box<Result<T, Refusal>>, Option<Charges>),
}

struct Board<T> {
    slots: Vec<Mutex<Slot<T>>>,
    changed: Condvar,
    next: AtomicUsize,
}

fn locked<T>(mutex: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// `task(0..count)`, in order, the first refusal ending the run.
pub(super) fn in_order<T, F>(count: usize, task: F, dec: &mut Decider) -> Result<Vec<T>, Refusal>
where
    T: Send,
    F: Fn(usize, &mut Decider) -> Result<T, Refusal> + Sync,
{
    let helpers = if count > 1 { Helpers::new(dec) } else { None };
    let Some(helpers) = helpers else {
        return (0..count).map(|index| task(index, dec)).collect();
    };
    let board = Board {
        slots: (0..count).map(|_| Mutex::new(Slot::Open)).collect(),
        changed: Condvar::new(),
        next: AtomicUsize::new(0),
    };
    let wait = catallax_numbers::run_ahead_waits();
    std::thread::scope(|scope| {
        struct Stop<'h>(&'h Helpers);
        impl Drop for Stop<'_> {
            fn drop(&mut self) {
                self.0.stop();
            }
        }
        let _stop = Stop(&helpers);
        let mut started = 0;
        for _ in 0..helpers.count() {
            let spawned = std::thread::Builder::new()
                .name("catallax-spline-fillet".into())
                .stack_size(helpers.stack_size())
                .spawn_scoped(scope, || {
                    helpers.run(|worker| explore(&task, &board, &helpers, worker));
                });
            if spawned.is_err() {
                break;
            }
            started += 1;
        }
        let mut out = Vec::with_capacity(count);
        for index in 0..count {
            out.push(take(&task, &board, index, wait && started > 0, dec)?);
        }
        Ok(out)
    })
}

fn explore<T, F>(task: &F, board: &Board<T>, helpers: &Helpers, worker: &mut Decider)
where
    F: Fn(usize, &mut Decider) -> Result<T, Refusal>,
{
    loop {
        if helpers.stopped() {
            return;
        }
        let index = board.next.fetch_add(1, Ordering::Relaxed);
        if index >= board.slots.len() {
            return;
        }
        {
            let mut slot = locked(&board.slots[index]);
            if !matches!(*slot, Slot::Open) {
                continue;
            }
            *slot = Slot::Taken;
        }
        let (value, charges) = catallax_numbers::record(|| task(index, worker));
        *locked(&board.slots[index]) = Slot::Done(Box::new(value), charges);
        board.changed.notify_all();
    }
}

fn take<T, F>(
    task: &F,
    board: &Board<T>,
    index: usize,
    wait: bool,
    dec: &mut Decider,
) -> Result<T, Refusal>
where
    F: Fn(usize, &mut Decider) -> Result<T, Refusal>,
{
    let mut slot = locked(&board.slots[index]);
    loop {
        match &*slot {
            Slot::Open if !wait => {
                *slot = Slot::Writer;
                drop(slot);
                return task(index, dec);
            }
            Slot::Open | Slot::Taken => {
                slot = board
                    .changed
                    .wait(slot)
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
            }
            Slot::Writer => unreachable!("the writer takes each task once"),
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
    task(index, dec)
}
