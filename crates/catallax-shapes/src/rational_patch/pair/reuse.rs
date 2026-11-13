//! Exact root requests reuse only committed levels and one cell's own work.

use catallax_numbers::{BudgetPolicy, IsolatedRoot, PolySystem, RootBox};
use std::{
    cell::RefCell,
    collections::{HashMap, hash_map::DefaultHasher},
    hash::{Hash, Hasher},
    sync::Arc,
};

#[allow(clippy::option_option)]
pub(super) type Answer = Option<Option<IsolatedRoot>>;
pub(super) type Snapshot = Arc<HashMap<u64, Vec<Arc<Certificate>>>>;

pub(super) struct Certificate {
    hash: u64,
    system: Arc<PolySystem>,
    bounds: RootBox,
    policy: BudgetPolicy,
    precision: u32,
    answer: Answer,
}

impl Certificate {
    fn matches(
        &self,
        system: &PolySystem,
        bounds: &RootBox,
        policy: BudgetPolicy,
        precision: u32,
    ) -> bool {
        self.policy == policy
            && self.precision == precision
            && self.bounds == *bounds
            && self.system.equations() == system.equations()
    }
}

#[derive(Default)]
pub(super) struct Cache {
    saved: Snapshot,
}

impl Cache {
    pub(super) fn snapshot(&self) -> Snapshot {
        self.saved.clone()
    }

    pub(super) fn commit(&mut self, records: Vec<Arc<Certificate>>) {
        if records.is_empty() {
            return;
        }
        let saved = Arc::make_mut(&mut self.saved);
        for record in records {
            let bucket = saved.entry(record.hash).or_default();
            if !bucket.iter().any(|known| {
                known.matches(
                    &record.system,
                    &record.bounds,
                    record.policy,
                    record.precision,
                )
            }) {
                bucket.push(record);
            }
        }
    }
}

struct Pending {
    base: Snapshot,
    records: Vec<Arc<Certificate>>,
}

thread_local! {
    static PENDING: RefCell<Option<Pending>> = const { RefCell::new(None) };
}

/// A helper's results become shared only after the entire level is adopted.
/// Every cell in that level sees the same immutable prefix, on every thread.
pub(super) fn scope<T>(
    base: Option<&Snapshot>,
    run: impl FnOnce() -> T,
) -> (T, Vec<Arc<Certificate>>) {
    struct Restore(Option<Pending>);
    impl Drop for Restore {
        fn drop(&mut self) {
            PENDING.with(|pending| {
                drop(pending.replace(self.0.take()));
            });
        }
    }
    let old = PENDING.with(|pending| {
        pending.replace(base.map(|base| Pending {
            base: base.clone(),
            records: Vec::new(),
        }))
    });
    let _restore = Restore(old);
    let value = run();
    let records = PENDING.with(|pending| {
        pending
            .borrow_mut()
            .as_mut()
            .map(|pending| std::mem::take(&mut pending.records))
            .unwrap_or_default()
    });
    (value, records)
}

pub(super) struct Request<'a> {
    hash: Option<u64>,
    system: &'a Arc<PolySystem>,
    bounds: &'a RootBox,
    policy: BudgetPolicy,
    precision: u32,
}

impl<'a> Request<'a> {
    pub(super) fn new(
        system: &'a Arc<PolySystem>,
        bounds: &'a RootBox,
        policy: BudgetPolicy,
        precision: u32,
    ) -> Self {
        let hash = PENDING.with(|pending| {
            pending.borrow().as_ref().map(|_| {
                let mut hasher = DefaultHasher::new();
                system.equations().hash(&mut hasher);
                bounds.lo().hash(&mut hasher);
                bounds.hi().hash(&mut hasher);
                policy.max_work.hash(&mut hasher);
                policy.max_precision_bits.hash(&mut hasher);
                precision.hash(&mut hasher);
                hasher.finish()
            })
        });
        Self {
            hash,
            system,
            bounds,
            policy,
            precision,
        }
    }

    pub(super) fn find(&self) -> Option<Answer> {
        let hash = self.hash?;
        PENDING.with(|pending| {
            let pending = pending.borrow();
            let pending = pending.as_ref()?;
            pending
                .base
                .get(&hash)
                .into_iter()
                .flatten()
                .chain(pending.records.iter().filter(|record| record.hash == hash))
                .find(|record| {
                    record.matches(self.system, self.bounds, self.policy, self.precision)
                })
                .map(|record| record.answer.clone())
        })
    }

    pub(super) fn remember(self, answer: &Answer) {
        let Some(hash) = self.hash else {
            return;
        };
        PENDING.with(|pending| {
            if let Some(pending) = pending.borrow_mut().as_mut() {
                pending.records.push(Arc::new(Certificate {
                    hash,
                    system: self.system.clone(),
                    bounds: self.bounds.clone(),
                    policy: self.policy,
                    precision: self.precision,
                    answer: answer.clone(),
                }));
            }
        });
    }
}
