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

// One lock guards every slot. The writer waits on `changed` for slot after
// slot, and a condition variable must always be used with the same mutex:
// macOS panics on a second one.
struct Board<T> {
    slots: Mutex<Vec<Slot<T>>>,
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
        slots: Mutex::new((0..count).map(|_| Slot::Open).collect()),
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
        {
            let mut slots = locked(&board.slots);
            match slots.get(index) {
                None => return,
                Some(Slot::Open) => slots[index] = Slot::Taken,
                Some(_) => continue,
            }
        }
        let (value, charges) = catallax_numbers::record(|| task(index, worker));
        locked(&board.slots)[index] = Slot::Done(Box::new(value), charges);
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
    let mut slots = locked(&board.slots);
    loop {
        match &slots[index] {
            Slot::Open if !wait => {
                slots[index] = Slot::Writer;
                drop(slots);
                return task(index, dec);
            }
            Slot::Open | Slot::Taken => {
                slots = board
                    .changed
                    .wait(slots)
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
            }
            Slot::Writer => unreachable!("the writer takes each task once"),
            Slot::Done(..) => break,
        }
    }
    let Slot::Done(value, charges) = std::mem::replace(&mut slots[index], Slot::Writer) else {
        unreachable!("matched above");
    };
    drop(slots);
    if let Some(charges) = charges
        && charges.adopt(dec)
    {
        return *value;
    }
    task(index, dec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::num::NonZeroUsize;
    use std::time::Duration;

    // Helpers prove the tasks while the writer waits for each one in turn, so
    // the writer waits on the board for task after task. Every wait must use
    // the one lock the condition variable has seen: macOS panics on another.
    #[test]
    fn the_writer_waits_for_task_after_task_on_one_lock() {
        struct Restore;
        impl Drop for Restore {
            fn drop(&mut self) {
                catallax_numbers::set_threads(None);
                catallax_numbers::set_run_ahead_wait(false);
            }
        }
        let _restore = Restore;
        catallax_numbers::set_threads(NonZeroUsize::new(4));
        catallax_numbers::set_run_ahead_wait(true);
        // Later tasks take longer, so the writer finds each one unfinished.
        let task = |index: usize, _: &mut Decider| {
            std::thread::sleep(Duration::from_millis(10 * (index as u64 + 1)));
            Ok(index)
        };
        let done = in_order(8, task, &mut Decider::default()).unwrap();
        assert_eq!(done, (0..8).collect::<Vec<_>>());
    }
}
