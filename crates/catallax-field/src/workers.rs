//! Shared field workers. Only independent arithmetic runs here; admission is ordered.

use rayon::prelude::*;
use std::{cell::Cell, sync::LazyLock};

static POOL: LazyLock<Option<rayon::ThreadPool>> = LazyLock::new(|| {
    rayon::ThreadPoolBuilder::new()
        .num_threads(std::thread::available_parallelism().map_or(1, usize::from))
        .thread_name(|i| format!("catallax-field-{i}"))
        .build()
        .ok()
});

thread_local! {
    static LIMIT: Cell<Option<usize>> = const { Cell::new(None) };
}

/// Number of available field workers. Defaults to all processors available to
/// this process. A pool that cannot start falls back to synchronous execution.
pub fn worker_count() -> usize {
    let available = POOL
        .as_ref()
        .map_or(1, rayon::ThreadPool::current_num_threads);
    LIMIT.with(|limit| limit.get().unwrap_or(available).min(available))
}

/// Limits this request's independent field jobs without rebuilding the shared
/// pool. The limit applies on the coordinator thread and restores on unwind.
/// Zero is refused. Other requests keep their own limits.
pub fn with_workers<T>(
    workers: usize,
    run: impl FnOnce() -> T,
) -> Result<T, catallax_numbers::Refusal> {
    if workers == 0 {
        return Err(catallax_numbers::Refusal::new(
            "field.workers",
            "at least one field worker is required",
        ));
    }
    struct Restore(Option<usize>);
    impl Drop for Restore {
        fn drop(&mut self) {
            LIMIT.with(|limit| limit.set(self.0));
        }
    }
    let _restore = Restore(LIMIT.with(|limit| limit.replace(Some(workers))));
    Ok(run())
}

pub(crate) fn map<T: Sync, R: Send>(items: &[T], run: impl Fn(&T) -> R + Sync) -> Vec<R> {
    map_init(items, || (), |_, item| run(item))
}

pub(crate) fn map_init<T: Sync, R: Send, S>(
    items: &[T],
    init: impl Fn() -> S + Sync,
    run: impl Fn(&mut S, &T) -> R + Sync,
) -> Vec<R> {
    if items.is_empty() {
        return Vec::new();
    }
    let workers = worker_count().min(items.len());
    if workers <= 1 {
        let mut state = init();
        return items.iter().map(|item| run(&mut state, item)).collect();
    }
    let chunk = items.len().div_ceil(workers);
    POOL.as_ref()
        .expect("multiple workers need a pool")
        .install(|| {
            items
                .par_chunks(chunk)
                .map(|items| {
                    with_workers(1, || {
                        let mut state = init();
                        items
                            .iter()
                            .map(|item| run(&mut state, item))
                            .collect::<Vec<_>>()
                    })
                    .expect("one worker is valid")
                })
                .collect::<Vec<_>>()
        })
        .into_iter()
        .flatten()
        .collect()
}

pub(crate) fn numeric_map<T: Sync, R: Send>(
    items: &[T],
    run: impl Fn(&mut crate::Workspace, &T) -> Result<R, catallax_numbers::Refusal> + Sync,
) -> Result<Vec<R>, catallax_numbers::Refusal> {
    use catallax_numbers::OperationTask;
    use std::{num::NonZeroUsize, sync::Mutex};
    let mut values = Vec::with_capacity(items.len());
    for items in items.chunks(256) {
        let tasks = OperationTask::capture().split(NonZeroUsize::new(items.len()).unwrap());
        let jobs = items
            .iter()
            .zip(tasks.into_iter().map(|task| Mutex::new(Some(task))))
            .collect::<Vec<_>>();
        let results = map_init(
            &jobs,
            crate::Workspace::default,
            |workspace, (item, task)| {
                let task = task.lock().unwrap().take().expect("one share per input");
                task.run(|| run(workspace, item))
            },
        );
        // Pay every completed share before admitting results or allocating retries.
        let admitted = results
            .into_iter()
            .map(|result| result.commit_share())
            .collect::<Result<Vec<_>, _>>()?;
        let mut workspace = crate::Workspace::default();
        for (item, value) in items.iter().zip(admitted) {
            let value = match value {
                Some(value) => value,
                None => OperationTask::capture()
                    .run(|| run(&mut workspace, item))
                    .commit()?,
            };
            values.push(value?);
        }
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computed_numeric_work_never_exceeds_the_one_request_cap() {
        use catallax_numbers::{Decider, with_operation_budget};
        use std::sync::atomic::{AtomicU64, Ordering};
        for workers in [1, 2, 4, 8, worker_count()] {
            with_workers(workers, || {
                for limit in [10, 1000] {
                    let physical = AtomicU64::new(0);
                    let (result, work) = with_operation_budget(limit, || {
                        numeric_map(&[0; 16], |_, &item| {
                            Decider::default().charge_work(6)?;
                            physical.fetch_add(6, Ordering::Relaxed);
                            Ok(item)
                        })
                    });
                    assert!(physical.load(Ordering::Relaxed) <= limit);
                    if limit == 10 {
                        assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
                        assert_eq!(work.used, 7);
                        assert!(work.exhausted);
                        assert_eq!(physical.load(Ordering::Relaxed), 6);
                    } else {
                        assert_eq!(result.unwrap(), vec![0; 16]);
                        assert_eq!(work.used, 112);
                        assert!(!work.exhausted);
                    }
                }
            })
            .unwrap();
        }
    }

    #[test]
    fn limits_restore_and_outputs_have_input_order() {
        let available = worker_count();
        assert!(with_workers(0, || ()).is_err());
        with_workers(1, || {
            assert_eq!(worker_count(), 1);
            let panic = std::panic::catch_unwind(|| with_workers(2, || panic!("unwind limit")));
            assert!(panic.is_err());
            assert_eq!(worker_count(), 1);
        })
        .unwrap();
        assert_eq!(worker_count(), available);
        for count in [1, 2, 4, 8, available] {
            with_workers(count, || {
                let items = (0..257).collect::<Vec<_>>();
                assert_eq!(
                    map(&items, |i| i * i),
                    items.iter().map(|i| i * i).collect::<Vec<_>>()
                );
            })
            .unwrap();
        }
    }
}
