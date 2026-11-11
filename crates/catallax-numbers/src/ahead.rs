//! Run-ahead proofs: every core works on a request, and the answer, the
//! refusals and the work count stay those of one thread.
//!
//! The request's own thread, the writer, runs the kernel exactly as it always
//! has: one model, one work counter, one set of proof caches, everything in
//! order. When it enters a loop whose items are independent proofs, it hands
//! the later items to worker threads. A worker runs an item on its own copy of
//! the model and throws the copy away. What it keeps is a record of each
//! decision it had to refine: the question, the outcome, and a tape of the work
//! the decision charged.
//!
//! When the writer reaches a decision and a record of the same question is
//! ready, it adopts the record instead of redoing the arithmetic. Adopting
//! replays the tape against the writer's own caches, so the writer charges what
//! it would have charged had it done the arithmetic itself:
//!
//! - work outside the proof caches is replayed as the amount recorded;
//! - each algebraic root enclosure goes through the writer's real cache
//!   protocol, which decides between a reuse and a cold proof exactly as it
//!   does without workers. Only the cold arithmetic is taken from the record.
//!
//! A record is used only when the writer's question is the recorded question
//! node for node (`structure.rs`), under the same policy. A decision's outcome
//! depends on nothing else: caches change what a proof costs, never what it
//! proves. A record that would not fit the remaining allowance is ignored and
//! the writer decides cold, so running out of the allowance happens at the
//! same charge as before.
//!
//! Nothing a worker does can reach the writer except through a record. Workers
//! have their own work counter and proof caches, and they neither read nor
//! fill the shared caches whose fills are charged (`shared_caches`).
//! With one thread no worker exists and no record is ever made.
//!
//! A whole computation can be recorded the same way ([`record`]), and the
//! writer then skips the work between its decisions too. The STEP pcurve fit
//! does that on every thread by default. Loops hand their items over only when
//! [`set_run_ahead_loops`] asks, because a record of one decision saves the
//! writer less than the workers cost it on the parts measured so far.

use crate::{
    BudgetPolicy, Decider, DecisionStats, Rational, Real, UPoly,
    algebraic::AlgebraicRoot,
    decider::{Bracket, Work},
    structure::{Pairing, fingerprint},
};
use std::{
    any::Any,
    cell::{Cell, RefCell},
    collections::{BTreeMap, HashMap},
    hash::{Hash, Hasher},
    num::NonZeroUsize,
    sync::{
        Arc, Condvar, Mutex, MutexGuard, OnceLock, PoisonError,
        atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering},
    },
    time::{Duration, Instant},
};

// ---------------------------------------------------------------------------
// Thread count
// ---------------------------------------------------------------------------

static FORCED_THREADS: AtomicUsize = AtomicUsize::new(0);
static LOOPS_TO_WORKERS: AtomicBool = AtomicBool::new(false);
static WAIT_FOR_WORKERS: AtomicBool = AtomicBool::new(false);

fn default_threads() -> usize {
    static DEFAULT: OnceLock<usize> = OnceLock::new();
    *DEFAULT.get_or_init(|| {
        std::env::var("CATALLAX_THREADS")
            .ok()
            .and_then(|value| value.trim().parse::<usize>().ok())
            .filter(|&count| count > 0)
            .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, NonZeroUsize::get))
    })
}

/// Whether the writer waits for its workers before it goes on: see
/// [`set_run_ahead_wait`]. Code with helper threads of its own asks this too.
#[doc(hidden)]
pub fn run_ahead_waits() -> bool {
    wait_for_workers()
}

// The `CATALLAX_RUN_AHEAD` environment variable: `loops` or `wait`.
fn run_ahead_from_environment(mode: &str) -> bool {
    static MODE: OnceLock<String> = OnceLock::new();
    MODE.get_or_init(|| {
        std::env::var("CATALLAX_RUN_AHEAD")
            .map(|value| value.trim().to_owned())
            .unwrap_or_default()
    }) == mode
}

fn wait_for_workers() -> bool {
    WAIT_FOR_WORKERS.load(Ordering::Relaxed) || run_ahead_from_environment("wait")
}

fn loops_go_to_workers() -> bool {
    LOOPS_TO_WORKERS.load(Ordering::Relaxed)
        || run_ahead_from_environment("loops")
        || wait_for_workers()
}

/// The number of threads the exact kernel uses for one request: the request's
/// own thread plus its workers. Every available core by default, or the count
/// in the `CATALLAX_THREADS` environment variable, or the count last passed to
/// [`set_threads`]. The count never changes an answer, a refusal or a work
/// count; it changes how long they take.
pub fn threads() -> usize {
    match FORCED_THREADS.load(Ordering::Relaxed) {
        0 => default_threads(),
        count => count,
    }
}

/// Forces the kernel's thread count for this process, or restores the default
/// with `None`. One thread runs every proof on the request's own thread.
pub fn set_threads(count: Option<NonZeroUsize>) {
    FORCED_THREADS.store(count.map_or(0, NonZeroUsize::get), Ordering::Relaxed);
}

/// Hands the loops marked with [`run_ahead!`](crate::run_ahead) to worker
/// threads, which decide their later proofs ahead of the request's own thread.
/// The `CATALLAX_RUN_AHEAD=loops` environment variable does the same.
///
/// Off by default. It changes no answer, refusal or work count, but on the
/// parts measured so far the workers slow the request's thread more than
/// their records save it (docs/RUN-AHEAD.md, section 5). Whole computations
/// recorded with [`record`] use every thread without this.
pub fn set_run_ahead_loops(on: bool) {
    LOOPS_TO_WORKERS.store(on, Ordering::Relaxed);
}

/// Hands loops to workers as [`set_run_ahead_loops`] does, and makes every
/// loop wait until its workers have finished before the writer goes on, so
/// that each recorded decision is adopted. Tests use this to compare complete
/// replays against one thread without depending on timing. The
/// `CATALLAX_RUN_AHEAD=wait` environment variable does the same.
#[doc(hidden)]
pub fn set_run_ahead_wait(wait: bool) {
    WAIT_FOR_WORKERS.store(wait, Ordering::Relaxed);
}

static ADOPTED: AtomicU64 = AtomicU64::new(0);

/// How many recorded decisions writers in this process have adopted so far.
/// For tests and timing reports; nothing in the kernel reads it.
#[doc(hidden)]
pub fn run_ahead_adopted() -> u64 {
    ADOPTED.load(Ordering::Relaxed)
}

// ---------------------------------------------------------------------------
// Fork count
// ---------------------------------------------------------------------------

static EPOCH: AtomicU64 = AtomicU64::new(1);

/// The number of loops handed to workers so far, plus one. Expression nodes
/// and roots keep the value they were made under: anything made before a
/// loop's fork is shared, unchanged, by the writer and that loop's workers.
pub(crate) fn epoch() -> u64 {
    EPOCH.load(Ordering::Relaxed)
}

// ---------------------------------------------------------------------------
// Questions and records
// ---------------------------------------------------------------------------

/// One decision, as asked.
pub(crate) enum Question {
    Sign(Real),
    /// To a width, and whether to stop at the first ball that excludes zero.
    Enclose(Real, Rational, bool),
    EncloseAt(Real, u32),
    EncloseBudgetedAt(Real, u32),
    EncloseAllAt(Vec<Real>, u32),
    RationalIdentity(Real),
    RationalValue(Real),
    Roots(Arc<UPoly>),
}

pub(crate) struct Key {
    pub question: Question,
    pub policy: BudgetPolicy,
}

impl Key {
    fn tag(&self) -> u8 {
        match &self.question {
            Question::Sign(_) => 0,
            Question::Enclose(..) => 1,
            Question::EncloseAt(..) => 2,
            Question::EncloseBudgetedAt(..) => 3,
            Question::EncloseAllAt(..) => 4,
            Question::RationalIdentity(_) => 5,
            Question::RationalValue(_) => 6,
            Question::Roots(_) => 7,
        }
    }

    fn hash(&self) -> u64 {
        let mut state = std::collections::hash_map::DefaultHasher::new();
        self.tag().hash(&mut state);
        self.policy.max_precision_bits.hash(&mut state);
        self.policy.max_work.hash(&mut state);
        match &self.question {
            Question::Sign(x) | Question::RationalIdentity(x) | Question::RationalValue(x) => {
                fingerprint(x).hash(&mut state);
            }
            Question::Enclose(x, width, stop_at_sign) => {
                fingerprint(x).hash(&mut state);
                width.hash(&mut state);
                stop_at_sign.hash(&mut state);
            }
            Question::EncloseAt(x, precision) | Question::EncloseBudgetedAt(x, precision) => {
                fingerprint(x).hash(&mut state);
                precision.hash(&mut state);
            }
            Question::EncloseAllAt(values, precision) => {
                for x in values {
                    fingerprint(x).hash(&mut state);
                }
                precision.hash(&mut state);
            }
            Question::Roots(polynomial) => polynomial.hash(&mut state),
        }
        state.finish()
    }

    // `self` is the recorded question.
    fn same(&self, mine: &Key, epoch: u64) -> Option<Pairing> {
        if self.policy != mine.policy {
            return None;
        }
        let mut pairing = Pairing::default();
        let same = match (&self.question, &mine.question) {
            (Question::Sign(a), Question::Sign(b))
            | (Question::RationalIdentity(a), Question::RationalIdentity(b))
            | (Question::RationalValue(a), Question::RationalValue(b)) => pairing.same(a, b, epoch),
            (Question::Enclose(a, w, s), Question::Enclose(b, v, t)) => {
                w == v && s == t && pairing.same(a, b, epoch)
            }
            (Question::EncloseAt(a, p), Question::EncloseAt(b, q))
            | (Question::EncloseBudgetedAt(a, p), Question::EncloseBudgetedAt(b, q)) => {
                p == q && pairing.same(a, b, epoch)
            }
            (Question::EncloseAllAt(a, p), Question::EncloseAllAt(b, q)) => {
                p == q
                    && a.len() == b.len()
                    && a.iter().zip(b).all(|(a, b)| pairing.same(a, b, epoch))
            }
            (Question::Roots(a), Question::Roots(b)) => a == b,
            _ => false,
        };
        same.then_some(pairing)
    }
}

/// The cold arithmetic of one root enclosure from the root's own interval.
pub(crate) struct RootProof {
    pub result: Result<(crate::Ball, Bracket), &'static str>,
    /// Request work the cold arithmetic charged.
    pub ops: u64,
    /// Decision work it used: the certificate's cost when it succeeded.
    pub used: u64,
}

/// One enclosure of an algebraic root inside a recorded decision.
struct RootSite {
    root: Arc<AlgebraicRoot>,
    precision: u32,
    policy: BudgetPolicy,
    /// Decision work used when the enclosure was asked for.
    used: u64,
    /// The cold arithmetic, or `None` when the worker reused a certificate
    /// the writer already held when it handed the work over.
    proof: Option<Arc<RootProof>>,
}

/// Cold root proofs shared by the workers of one loop or one computation.
type SharedRoots = Mutex<BTreeMap<crate::operation::RootProofKey, Arc<RootProof>>>;

/// The cold arithmetic of one checked pi enclosure.
struct PiProof {
    result: Result<crate::Ball, crate::Refusal>,
    /// Request work the series charged.
    ops: u64,
}

enum Item {
    /// Request work charged outside the proof caches.
    Work(u64),
    Root(RootSite),
    /// A checked pi enclosure at this precision, under this policy.
    Pi(u32, BudgetPolicy, Arc<PiProof>),
}

struct Record {
    outcome: Box<dyn Any + Send + Sync>,
    stats: DecisionStats,
    tape: Vec<Item>,
    /// No replay of the tape can charge more than this.
    most: u64,
}

enum Slot {
    Pending,
    Ready(Arc<Record>),
    Abandoned,
}

struct Entry {
    key: Key,
    slot: Mutex<Slot>,
    settled: Condvar,
    /// When a worker began deciding it.
    claimed: Instant,
}

/// A decision is worth a record once it has taken a worker this long. Finding
/// a record and proving it answers the same question costs about as much as
/// a quick decision does, so quick ones are decided again instead.
const WORTH_A_RECORD: Duration = Duration::from_micros(250);

/// A worker that has been at one decision this long says so, for other
/// threads to wait for its record instead of starting the same proof.
const WORTH_A_WAIT: Duration = Duration::from_millis(1);

fn lock<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(PoisonError::into_inner)
}

impl Entry {
    /// The record, once its worker has made it. A waiting thread gives the
    /// worker as long again as it has already had, then decides for itself: a
    /// worker the machine has set aside must not hold the request up.
    fn wait(&self) -> Option<Arc<Record>> {
        let mut slot = lock(&self.slot);
        let until = (!wait_for_workers()).then(|| Instant::now() + self.claimed.elapsed());
        loop {
            match &*slot {
                Slot::Ready(record) => return Some(record.clone()),
                Slot::Abandoned => return None,
                Slot::Pending => match until {
                    None => {
                        slot = self
                            .settled
                            .wait(slot)
                            .unwrap_or_else(PoisonError::into_inner);
                    }
                    Some(until) => {
                        let left = until.saturating_duration_since(Instant::now());
                        if left.is_zero() {
                            return None;
                        }
                        slot = self
                            .settled
                            .wait_timeout(slot, left)
                            .unwrap_or_else(PoisonError::into_inner)
                            .0;
                    }
                },
            }
        }
    }

    fn settle(&self, slot: Slot) {
        *lock(&self.slot) = slot;
        self.settled.notify_all();
    }
}

const SHARDS: usize = 16;

/// The records of one loop.
struct Table {
    epoch: u64,
    shards: [Mutex<HashMap<u64, Vec<Arc<Entry>>>>; SHARDS],
}

impl Table {
    fn new(epoch: u64) -> Self {
        Self {
            epoch,
            shards: std::array::from_fn(|_| Mutex::new(HashMap::new())),
        }
    }

    fn shard(&self, hash: u64) -> &Mutex<HashMap<u64, Vec<Arc<Entry>>>> {
        &self.shards[(hash >> 32) as usize % SHARDS]
    }

    fn find(&self, key: &Key, hash: u64) -> Option<(Arc<Entry>, Pairing)> {
        let candidates = lock(self.shard(hash)).get(&hash).cloned()?;
        candidates.into_iter().find_map(|entry| {
            let pairing = entry.key.same(key, self.epoch)?;
            Some((entry, pairing))
        })
    }

    /// Enters `key` with `slot` and returns the entry, unless the question
    /// is already entered.
    fn enter(&self, key: Key, hash: u64, slot: Slot, began: Instant) -> Option<Arc<Entry>> {
        let mut shard = lock(self.shard(hash));
        let entries = shard.entry(hash).or_default();
        for entry in entries.iter() {
            if entry.key.same(&key, self.epoch).is_some() {
                return None;
            }
        }
        let entry = Arc::new(Entry {
            key,
            slot: Mutex::new(slot),
            settled: Condvar::new(),
            claimed: began,
        });
        entries.push(entry.clone());
        Some(entry)
    }

    fn remove(&self, entry: &Arc<Entry>, hash: u64) {
        let mut shard = lock(self.shard(hash));
        if let Some(entries) = shard.get_mut(&hash) {
            entries.retain(|other| !Arc::ptr_eq(other, entry));
            if entries.is_empty() {
                shard.remove(&hash);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Thread state
// ---------------------------------------------------------------------------

struct Recording {
    tape: Vec<Item>,
    /// Request work at the last tape item.
    mark: u64,
    most: u64,
    unusable: bool,
}

struct Worker {
    /// Cold root proofs of this item, by root definition and precision.
    roots: BTreeMap<crate::operation::RootProofKey, Arc<RootProof>>,
    /// Cold root proofs of every worker beside this one.
    shared: Arc<SharedRoots>,
    /// The writer's root proofs when it handed the work over.
    inherited: crate::operation::Proofs,
    /// Cold pi enclosures of this item, by precision.
    pis: BTreeMap<u32, Arc<PiProof>>,
    recording: Option<Recording>,
    /// The decision being recorded, until it ends.
    deciding: Option<Deciding>,
    /// The fork this worker runs under: older roots are shared with the writer.
    epoch: u64,
}

struct Deciding {
    /// The question, until a claim moves it into the table.
    key: Option<Key>,
    hash: u64,
    began: Instant,
    /// The claim, once the decision has run long enough to be waited for.
    claim: Option<Arc<Entry>>,
}

thread_local! {
    // Whether a decision on this thread may use a table at all.
    static ENGAGED: Cell<bool> = const { Cell::new(false) };
    static IN_DECISION: Cell<bool> = const { Cell::new(false) };
    static RECORDING: Cell<bool> = const { Cell::new(false) };
    // The writer's open loops, outermost first.
    static WRITER: RefCell<Vec<Arc<Table>>> = const { RefCell::new(Vec::new()) };
    static WORKER: RefCell<Option<Worker>> = const { RefCell::new(None) };
    // The item this worker thread is running. Kept apart from `WORKER`, which
    // is borrowed while a tape is written.
    static ITEM: RefCell<Option<(Arc<Job>, usize)>> = const { RefCell::new(None) };
    // Set when the writer no longer wants what this helper thread is doing.
    static STOPPED: RefCell<Option<Arc<AtomicBool>>> = const { RefCell::new(None) };
    // Inside a nested allowance on a worker thread: nothing is recorded.
    static SUSPENDED: Cell<bool> = const { Cell::new(false) };
}

/// Whether this thread is a worker running a loop item ahead of the writer.
pub(crate) fn is_worker() -> bool {
    crate::operation::worker()
}

/// Whether this thread may read and fill a shared cache whose fill is charged.
///
/// A worker must not: the writer would find the cache filled, skip the charge
/// and report less work than one thread does. Call this at every cache that
/// lives on a shared object (a curve, a surface, a number) and whose miss
/// spends decision work; on `false`, compute the value without the cache.
pub fn shared_caches() -> bool {
    if !is_worker() {
        return true;
    }
    // The writer might have found the cache filled and charged less than
    // this worker is about to: whatever is being recorded cannot be replayed.
    unrecordable();
    false
}

/// Tells a recording that the code it records read a cache of this thread's
/// own whose hits spend less decision work than its misses, such as the
/// intersection plans. The writer's cache differs, so the record is dropped
/// and the writer runs the code itself. Does nothing outside a recording.
pub fn private_cache() {
    unrecordable();
}

/// A worker's regular look round, every few dozen charges. Returns whether
/// its item is no longer wanted. A decision that has run long enough is
/// entered in the table as pending, so other threads wait for its record.
pub(crate) fn worker_tick() -> bool {
    if item_cancelled() {
        return true;
    }
    if RECORDING.get() {
        WORKER.with(|worker| {
            let Ok(mut worker) = worker.try_borrow_mut() else {
                return;
            };
            let Some(deciding) = worker.as_mut().and_then(|worker| worker.deciding.as_mut()) else {
                return;
            };
            if deciding.claim.is_some()
                || deciding.key.is_none()
                || deciding.began.elapsed() < WORTH_A_WAIT
            {
                return;
            }
            let table = ITEM.with(|item| {
                item.borrow()
                    .as_ref()
                    .and_then(|(job, _)| job.tables.last().cloned())
            });
            if let (Some(table), Some(key)) = (table, deciding.key.take()) {
                deciding.claim = table.enter(key, deciding.hash, Slot::Pending, deciding.began);
            }
        });
    }
    false
}

/// Whether this worker's item is no longer wanted.
fn item_cancelled() -> bool {
    ITEM.with(|item| {
        item.borrow().as_ref().is_some_and(|(job, index)| {
            job.cancelled.load(Ordering::Relaxed) || job.floor.load(Ordering::Relaxed) > *index
        })
    }) || STOPPED.with(|stopped| {
        stopped
            .borrow()
            .as_ref()
            .is_some_and(|stopped| stopped.load(Ordering::Relaxed))
    })
}

/// Whether a recorded decision is in progress on this worker thread.
pub(crate) fn recording() -> bool {
    RECORDING.get()
}

/// Marks the decision being recorded as one that cannot be replayed: it
/// touched a proof cache the tape does not describe.
pub(crate) fn unrecordable() {
    if RECORDING.get() {
        WORKER.with(|worker| {
            if let Some(recording) = worker
                .borrow_mut()
                .as_mut()
                .and_then(|worker| worker.recording.as_mut())
            {
                recording.unusable = true;
            }
        });
    }
}

/// A worker thread's recording, set aside while a nested allowance runs.
pub(crate) struct Suspended {
    worker: bool,
    recording: bool,
    engaged: bool,
    suspended: bool,
}

/// Stops this worker thread recording until the value is dropped, and marks
/// what it was recording as unusable. A nested allowance starts its counter
/// from zero and sets the proof caches aside; a tape describes neither. With
/// `worker` false this does nothing: the writer's own scopes are its business.
pub(crate) fn suspend(worker: bool) -> Suspended {
    if !worker {
        return Suspended {
            worker,
            recording: false,
            engaged: false,
            suspended: false,
        };
    }
    unrecordable();
    Suspended {
        worker,
        recording: RECORDING.replace(false),
        engaged: ENGAGED.replace(false),
        suspended: SUSPENDED.replace(true),
    }
}

impl Drop for Suspended {
    fn drop(&mut self) {
        if self.worker {
            RECORDING.set(self.recording);
            ENGAGED.set(self.engaged);
            SUSPENDED.set(self.suspended);
        }
    }
}

struct InDecision(bool);
impl InDecision {
    fn enter() -> Self {
        Self(IN_DECISION.replace(true))
    }
}
impl Drop for InDecision {
    fn drop(&mut self) {
        IN_DECISION.set(self.0);
    }
}

// ---------------------------------------------------------------------------
// Deciding
// ---------------------------------------------------------------------------

/// Where a decision keeps its statistics, and its deadline.
pub(crate) trait Host {
    fn stats(&mut self) -> Option<&mut DecisionStats>;
    fn expired(&self) -> bool;
}

/// Whether a decision on this thread should be offered to [`decide`].
pub(crate) fn engaged() -> bool {
    ENGAGED.get() && !IN_DECISION.get()
}

/// Decides `key` by `cold`, or adopts a record of it. The flag is true when a
/// record was adopted: `cold` did not run and changed nothing of the host's.
pub(crate) fn decide<H: Host, T: Clone + Send + Sync + 'static>(
    host: &mut H,
    key: Key,
    cold: impl FnOnce(&mut H) -> T,
) -> (T, bool) {
    let hash = key.hash();
    if is_worker() {
        return decide_ahead(host, key, hash, cold);
    }
    let found = WRITER.with(|tables| {
        tables.borrow().iter().rev().find_map(|table| {
            let (entry, pairing) = table.find(&key, hash)?;
            Some((entry, pairing, table.epoch))
        })
    });
    if let Some((entry, pairing, epoch)) = found
        && let Some(record) = entry.wait()
        && !host.expired()
        && crate::operation::affordable(record.most)
        && inherited_proofs_held(&record.tape)
    {
        if let Some(stats) = host.stats() {
            stats.absorb(&record.stats);
        }
        replay(&record.tape, &pairing, epoch);
        ADOPTED.fetch_add(1, Ordering::Relaxed);
        let outcome = record
            .outcome
            .downcast_ref::<T>()
            .expect("a question has one outcome type");
        return (outcome.clone(), true);
    }
    let _decision = InDecision::enter();
    (cold(host), false)
}

fn decide_ahead<H: Host, T: Clone + Send + Sync + 'static>(
    host: &mut H,
    key: Key,
    hash: u64,
    cold: impl FnOnce(&mut H) -> T,
) -> (T, bool) {
    let tables = ITEM.with(|item| {
        item.borrow()
            .as_ref()
            .map(|(job, _)| job.tables.clone())
            .unwrap_or_default()
    });
    let Some(own) = tables.last() else {
        let _decision = InDecision::enter();
        return (cold(host), false);
    };
    // A record from this loop or an enclosing one serves this worker too.
    if let Some(record) = tables
        .iter()
        .rev()
        .find_map(|table| table.find(&key, hash))
        .and_then(|(entry, _)| entry.wait())
    {
        if let Some(stats) = host.stats() {
            stats.absorb(&record.stats);
        }
        let outcome = record
            .outcome
            .downcast_ref::<T>()
            .expect("a question has one outcome type");
        return (outcome.clone(), true);
    }

    // Whatever the decision leaves behind is cleared on every path out,
    // including a panic: a claim nobody settled would be waited for.
    struct Recorded<'a>(&'a Table);
    impl Drop for Recorded<'_> {
        fn drop(&mut self) {
            RECORDING.set(false);
            let deciding = WORKER.with(|worker| {
                let mut worker = worker.try_borrow_mut().ok()?;
                let worker = worker.as_mut()?;
                worker.recording = None;
                worker.deciding.take()
            });
            if let Some(Deciding {
                claim: Some(entry),
                hash,
                ..
            }) = deciding
            {
                self.0.remove(&entry, hash);
                entry.settle(Slot::Abandoned);
            }
        }
    }
    let began = Instant::now();
    let saved = host.stats().map(std::mem::take);
    WORKER.with(|worker| {
        if let Some(worker) = worker.borrow_mut().as_mut() {
            worker.recording = Some(Recording {
                tape: Vec::new(),
                mark: crate::operation::used(),
                most: 0,
                unusable: false,
            });
            worker.deciding = Some(Deciding {
                key: Some(key),
                hash,
                began,
                claim: None,
            });
        }
    });
    RECORDING.set(true);
    let recorded = Recorded(own);
    let outcome = {
        let _decision = InDecision::enter();
        cold(host)
    };
    RECORDING.set(false);
    let (recording, deciding) = WORKER.with(|worker| {
        let mut worker = worker.borrow_mut();
        match worker.as_mut() {
            Some(worker) => (worker.recording.take(), worker.deciding.take()),
            None => (None, None),
        }
    });
    let mut delta = DecisionStats::default();
    if let (Some(stats), Some(saved)) = (host.stats(), saved) {
        delta = std::mem::replace(stats, saved);
        stats.absorb(&delta);
    }
    // A cancelled item latches the work counter as exhausted; whatever it
    // decided after that is not a decision.
    let record = recording
        .filter(|recording| !recording.unusable && !crate::operation::exhausted())
        .map(|mut recording| {
            let rest = crate::operation::used() - recording.mark;
            if rest != 0 {
                recording.tape.push(Item::Work(rest));
                recording.most = recording.most.saturating_add(rest);
            }
            Arc::new(Record {
                outcome: Box::new(outcome.clone()),
                stats: delta,
                tape: recording.tape,
                most: recording.most,
            })
        });
    if let Some(Deciding {
        key,
        hash,
        began,
        claim,
    }) = deciding
    {
        match (claim, record) {
            (Some(entry), Some(record)) => entry.settle(Slot::Ready(record)),
            (Some(entry), None) => {
                own.remove(&entry, hash);
                entry.settle(Slot::Abandoned);
            }
            // A quick decision is not worth the writer's search for it,
            // unless a test wants every record adopted.
            (None, Some(record)) => {
                if let Some(key) = key
                    && (began.elapsed() >= WORTH_A_RECORD || wait_for_workers())
                {
                    own.enter(key, hash, Slot::Ready(record), began);
                }
            }
            (None, None) => {}
        }
    }
    drop(recorded);
    (outcome, false)
}

/// Whether the writer still holds every certificate a worker reused from it.
/// It does, short of a nested allowance that set its caches aside; a tape that
/// reused one the writer lacks cannot say what the proof costs cold.
fn inherited_proofs_held(tape: &[Item]) -> bool {
    tape.iter().all(|item| match item {
        Item::Root(site) if site.proof.is_none() => {
            crate::operation::root_proof(&site.root, site.precision).is_some_and(|certificate| {
                certificate.cost <= site.policy.max_work.saturating_sub(site.used)
            })
        }
        _ => true,
    })
}

/// Charges the writer what a cold decision would charge it now.
fn replay(tape: &[Item], pairing: &Pairing, epoch: u64) {
    // Roots the decision made for itself are made afresh for each replay, as a
    // cold decision makes them.
    let mut made: HashMap<*const AlgebraicRoot, Arc<AlgebraicRoot>> = HashMap::new();
    for item in tape {
        match item {
            Item::Work(units) => {
                crate::operation::charge(*units);
            }
            Item::Pi(precision, policy, proof) => {
                let dec = Decider::new(*policy);
                let _ = crate::ball::pi_enclosure(*precision, Some(&dec), &mut |_| {
                    crate::operation::charge(proof.ops);
                    proof.result.clone()
                });
            }
            Item::Root(site) => {
                let root = if let Some(root) = pairing.root(&site.root) {
                    root
                } else if site.root.epoch < epoch {
                    site.root.clone()
                } else {
                    made.entry(Arc::as_ptr(&site.root))
                        .or_insert_with(|| Arc::new(AlgebraicRoot::clone(&site.root)))
                        .clone()
                };
                let mut work = Work::resume(site.policy, site.used);
                let proof = &site.proof;
                // The decision's outcome is the record's; only the charges and
                // the cache entries of this call matter here.
                let _ = crate::decider::enclose_root(
                    &root,
                    site.precision,
                    &mut work,
                    true,
                    &mut |bracket, scratch, work| match proof {
                        Some(proof) if scratch => {
                            work.adopt(proof.ops, proof.used);
                            proof.result.clone()
                        }
                        _ => crate::decider::bisect_root(&root, site.precision, bracket, work),
                    },
                );
            }
        }
    }
}

/// A worker's enclosure of an algebraic root inside a recorded decision.
/// The decision's own work is what one thread's would be: a proof already
/// known is reused under the same rule as the request's proof cache.
pub(crate) fn record_root(
    root: &Arc<AlgebraicRoot>,
    precision: u32,
    work: &mut Work,
) -> Result<crate::Ball, &'static str> {
    let key = crate::operation::RootProofKey::new(root, precision);
    let allowance = work.policy.max_work.saturating_sub(work.used);
    let (policy, used) = (work.policy, work.used);
    // This item's proofs, then the writer's certificates, then the other
    // workers' proofs.
    enum Known {
        Proof(Arc<RootProof>),
        Certificate(Arc<crate::algebraic::RootEnclosure>),
    }
    let state = WORKER.with(|worker| {
        let mut worker = worker.borrow_mut();
        let worker = worker.as_mut()?;
        if let Some(recording) = worker.recording.as_mut() {
            let before = crate::operation::used() - recording.mark;
            if before != 0 {
                recording.tape.push(Item::Work(before));
                recording.most = recording.most.saturating_add(before);
            }
        }
        let known = worker
            .roots
            .get(&key)
            .filter(|proof| proof.used <= allowance)
            .cloned()
            .map(Known::Proof)
            .or_else(|| {
                worker
                    .inherited
                    .root(root, precision)
                    .filter(|certificate| certificate.cost <= allowance)
                    .map(Known::Certificate)
            });
        Some((known, worker.shared.clone()))
    });
    let Some((known, shared)) = state else {
        // Recording is only ever on with a worker's state in place.
        return crate::decider::bisect_root(root, precision, Bracket::of_root(root), work)
            .map(|(ball, _)| ball);
    };
    let known = known.or_else(|| {
        lock(&shared)
            .get(&key)
            .filter(|proof| proof.used <= allowance)
            .cloned()
            .map(Known::Proof)
    });
    let (proof, result) = match known {
        Some(Known::Proof(proof)) => {
            let result = work
                .reuse_proof(proof.used)
                .and_then(|()| proof.result.clone())
                .map(|(ball, _)| ball);
            (Some(proof), result)
        }
        Some(Known::Certificate(certificate)) => {
            let result = work
                .reuse_proof(certificate.cost)
                .map(|()| certificate.ball.clone());
            (None, result)
        }
        None => {
            let ops = crate::operation::used();
            let result = crate::decider::bisect_root(root, precision, Bracket::of_root(root), work);
            let proof = Arc::new(RootProof {
                result: result.clone(),
                ops: crate::operation::used() - ops,
                used: work.used - used,
            });
            if proof.result.is_ok() {
                lock(&shared)
                    .entry(key.clone())
                    .or_insert_with(|| proof.clone());
            }
            (Some(proof), result.map(|(ball, _)| ball))
        }
    };
    WORKER.with(|worker| {
        let mut worker = worker.borrow_mut();
        let Some(worker) = worker.as_mut() else {
            return;
        };
        if let Some(proof) = &proof
            && proof.result.is_ok()
        {
            worker.roots.entry(key).or_insert_with(|| proof.clone());
        }
        if let Some(recording) = worker.recording.as_mut() {
            // A reuse charges one unit, a cached certificate its cost and one,
            // a cold proof its own work; a proof continued from a shorter
            // certificate charges less than both together. A certificate the
            // writer held is a reuse, which is checked before the replay.
            recording.most = recording.most.saturating_add(match &proof {
                Some(proof) => proof.ops.saturating_add(proof.used).saturating_add(1),
                None => 1,
            });
            recording.tape.push(Item::Root(RootSite {
                root: root.clone(),
                precision,
                policy,
                used,
                proof,
            }));
            recording.mark = crate::operation::used();
        }
    });
    result
}

/// A worker's checked pi enclosure inside recorded code.
pub(crate) fn record_pi(precision: u32, dec: &Decider) -> Result<crate::Ball, crate::Refusal> {
    let known = WORKER.with(|worker| {
        let mut worker = worker.borrow_mut();
        let worker = worker.as_mut()?;
        if let Some(recording) = worker.recording.as_mut() {
            let before = crate::operation::used() - recording.mark;
            if before != 0 {
                recording.tape.push(Item::Work(before));
                recording.most = recording.most.saturating_add(before);
            }
        }
        worker.pis.get(&precision).cloned()
    });
    let (proof, result) = match known {
        Some(proof) => {
            // A reuse still honours cancellation, as the request's cache does.
            let result = dec.check_deadline().and_then(|()| proof.result.clone());
            (proof, result)
        }
        None => {
            let ops = crate::operation::used();
            let result = dec
                .check_deadline()
                .and_then(|()| crate::ball::pi_series(precision, Some(dec)));
            // The entry check belongs to the protocol, not to the series.
            let proof = Arc::new(PiProof {
                result: result.clone(),
                ops: (crate::operation::used() - ops).saturating_sub(1),
            });
            (proof, result)
        }
    };
    WORKER.with(|worker| {
        let mut worker = worker.borrow_mut();
        let Some(worker) = worker.as_mut() else {
            return;
        };
        if proof.result.is_ok() {
            worker.pis.entry(precision).or_insert_with(|| proof.clone());
        }
        if let Some(recording) = worker.recording.as_mut() {
            recording.most = recording.most.saturating_add(proof.ops).saturating_add(1);
            recording
                .tape
                .push(Item::Pi(precision, *dec.policy(), proof));
            recording.mark = crate::operation::used();
        }
    });
    result
}

// ---------------------------------------------------------------------------
// Whole computations recorded on helper threads
// ---------------------------------------------------------------------------

/// The request work one computation charged on a helper thread, in a form
/// the writer can charge itself.
pub struct Charges {
    tape: Vec<Item>,
    most: u64,
    epoch: u64,
}

impl Charges {
    /// Charges this thread what the recorded computation would charge it now,
    /// and leaves its proof caches as that computation would. Returns `false`
    /// and charges nothing when the work might not fit the remaining
    /// allowance; the caller then runs the computation itself.
    ///
    /// Only for a computation whose result depends on nothing but values the
    /// recording thread and this thread share: its charges may depend on the
    /// proof caches, and this accounts for that, but nothing else of this
    /// thread's may have steered it.
    ///
    /// `dec` is the decider the computation would have used here. Past its
    /// deadline nothing is adopted: the computation run here would refuse.
    pub fn adopt(&self, dec: &Decider) -> bool {
        if Host::expired(dec)
            || !crate::operation::affordable(self.most)
            || !inherited_proofs_held(&self.tape)
        {
            return false;
        }
        replay(&self.tape, &Pairing::default(), self.epoch);
        ADOPTED.fetch_add(1, Ordering::Relaxed);
        true
    }
}

/// Runs `code` on a helper thread and records what it charges, for the writer
/// to adopt with the value. The record is `None` when the code touched a
/// proof cache a record cannot describe, or was stopped; the value is then
/// not to be used either.
pub fn record<T>(code: impl FnOnce() -> T) -> (T, Option<Charges>) {
    if !is_worker() || RECORDING.get() || SUSPENDED.get() {
        return (code(), None);
    }
    struct Stop;
    impl Drop for Stop {
        fn drop(&mut self) {
            RECORDING.set(false);
            WORKER.with(|worker| {
                if let Ok(mut worker) = worker.try_borrow_mut()
                    && let Some(worker) = worker.as_mut()
                {
                    worker.recording = None;
                }
            });
        }
    }
    let epoch = WORKER.with(|worker| {
        let mut worker = worker.borrow_mut();
        let worker = worker.as_mut()?;
        worker.recording = Some(Recording {
            tape: Vec::new(),
            mark: crate::operation::used(),
            most: 0,
            unusable: false,
        });
        Some(worker.epoch)
    });
    let Some(epoch) = epoch else {
        return (code(), None);
    };
    RECORDING.set(true);
    let stop = Stop;
    let value = {
        // Decisions inside the code are part of its record.
        let _decision = InDecision::enter();
        code()
    };
    let recording = WORKER.with(|worker| {
        worker
            .borrow_mut()
            .as_mut()
            .and_then(|worker| worker.recording.take())
    });
    drop(stop);
    let charges = recording
        .filter(|recording| !recording.unusable && !crate::operation::exhausted())
        .map(|mut recording| {
            let rest = crate::operation::used() - recording.mark;
            if rest != 0 {
                recording.tape.push(Item::Work(rest));
                recording.most = recording.most.saturating_add(rest);
            }
            Charges {
                tape: recording.tape,
                most: recording.most,
                epoch,
            }
        });
    (value, charges)
}

/// Helper threads for one computation of the writer's, such as one curve fit.
///
/// The writer makes this, starts its own scoped threads and has each call
/// [`Helpers::run`]. Inside, a helper [`record`]s pieces of the computation;
/// the writer adopts their [`Charges`] in its own order.
pub struct Helpers {
    policy: BudgetPolicy,
    /// The writer's deadline: a helper refuses when the writer would.
    deadline: Option<Instant>,
    reduced_angle_series: bool,
    reuse_enclosures: bool,
    epoch: u64,
    stopped: Arc<AtomicBool>,
    proofs: crate::operation::Proofs,
    roots: Arc<SharedRoots>,
}

impl Helpers {
    /// Helpers for the calling writer, or `None` with one thread, on a worker
    /// thread, or where a test clock counts this thread's deadline checks.
    pub fn new(dec: &Decider) -> Option<Self> {
        if threads() < 2
            || is_worker()
            || crate::decider::test_clock_installed()
            || cfg!(feature = "work-profile")
        {
            return None;
        }
        let (reduced_angle_series, reuse_enclosures) = crate::operation::preferences();
        Some(Self {
            policy: *dec.policy(),
            deadline: dec.deadline(),
            reduced_angle_series,
            reuse_enclosures,
            epoch: EPOCH.fetch_add(1, Ordering::Relaxed) + 1,
            stopped: Arc::new(AtomicBool::new(false)),
            proofs: crate::operation::proofs(),
            roots: Arc::default(),
        })
    }

    /// How many helper threads to start: the kernel's threads less the writer.
    pub fn count(&self) -> usize {
        threads().saturating_sub(1)
    }

    /// The stack size a helper thread needs to run kernel code.
    pub fn stack_size(&self) -> usize {
        WORKER_STACK
    }

    /// Tells the helpers to stop: their charges begin to fail, so whatever
    /// they are computing ends early and records nothing.
    pub fn stop(&self) {
        self.stopped.store(true, Ordering::Relaxed);
    }

    /// Whether [`Helpers::stop`] has been called.
    pub fn stopped(&self) -> bool {
        self.stopped.load(Ordering::Relaxed)
    }

    /// Runs `work` on the calling helper thread with a work counter and proof
    /// caches of its own and a decider with the writer's policy. Nothing it
    /// charges reaches the writer except through [`record`].
    pub fn run<T>(&self, work: impl FnOnce(&mut Decider) -> T) -> T {
        struct Leave;
        impl Drop for Leave {
            fn drop(&mut self) {
                IN_DECISION.set(false);
                RECORDING.set(false);
                STOPPED.with(|stopped| {
                    if let Ok(mut stopped) = stopped.try_borrow_mut() {
                        *stopped = None;
                    }
                });
                WORKER.with(|worker| {
                    if let Ok(mut worker) = worker.try_borrow_mut() {
                        *worker = None;
                    }
                });
                crate::operation::leave_worker();
            }
        }
        crate::operation::enter_worker(
            self.reduced_angle_series,
            self.reuse_enclosures,
            &self.proofs,
        );
        STOPPED.with(|stopped| *stopped.borrow_mut() = Some(self.stopped.clone()));
        WORKER.with(|worker| {
            *worker.borrow_mut() = Some(Worker {
                roots: BTreeMap::new(),
                shared: self.roots.clone(),
                inherited: self.proofs.clone(),
                pis: BTreeMap::new(),
                recording: None,
                deciding: None,
                epoch: self.epoch,
            });
        });
        let _leave = Leave;
        let mut dec = Decider::new(self.policy);
        dec.set_deadline(self.deadline);
        work(&mut dec)
    }
}

// ---------------------------------------------------------------------------
// Loops handed to workers
// ---------------------------------------------------------------------------

type Body = Box<dyn Fn(usize, &mut Decider) + Send + Sync>;

struct Job {
    body: Body,
    len: usize,
    /// The next item no worker has taken.
    next: AtomicUsize,
    /// Items below this are done or under way on the writer.
    floor: AtomicUsize,
    cancelled: AtomicBool,
    policy: BudgetPolicy,
    /// The writer's deadline: a worker refuses when the writer would.
    deadline: Option<Instant>,
    reduced_angle_series: bool,
    reuse_enclosures: bool,
    /// The fork this job's workers run under.
    epoch: u64,
    /// The enclosing loops' tables, then this loop's own. A loop of whole
    /// items has none.
    tables: Vec<Arc<Table>>,
    /// The writer's root proofs when it handed the loop over.
    proofs: crate::operation::Proofs,
    /// Cold root proofs of this loop's workers.
    roots: Arc<SharedRoots>,
    /// Items workers have finished or will never start.
    finished: Mutex<usize>,
    all_finished: Condvar,
}

impl Job {
    fn take(&self) -> Option<usize> {
        loop {
            if self.cancelled.load(Ordering::Relaxed) {
                return None;
            }
            let index = self.next.fetch_add(1, Ordering::Relaxed);
            if index >= self.len {
                self.next.store(self.len, Ordering::Relaxed);
                return None;
            }
            if index > self.floor.load(Ordering::Relaxed) {
                return Some(index);
            }
            self.finish();
        }
    }

    fn exhausted(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed) || self.next.load(Ordering::Relaxed) >= self.len
    }

    fn finish(&self) {
        *lock(&self.finished) += 1;
        self.all_finished.notify_all();
    }
}

struct Pool {
    state: Mutex<PoolState>,
    wake: Condvar,
}

struct PoolState {
    /// Open loops with items left, oldest first.
    jobs: Vec<Arc<Job>>,
    workers: usize,
}

static POOL: Pool = Pool {
    state: Mutex::new(PoolState {
        jobs: Vec::new(),
        workers: 0,
    }),
    wake: Condvar::new(),
};

// Kernel recursion follows expression depth. A worker runs the same code as
// the writer, so it needs at least the writer's stack; the reservation is
// address space, not memory.
#[cfg(target_pointer_width = "64")]
const WORKER_STACK: usize = 256 << 20;
#[cfg(not(target_pointer_width = "64"))]
const WORKER_STACK: usize = 8 << 20;

// Returns whether any worker exists to take the job.
fn submit(job: &Arc<Job>, items: usize) -> bool {
    let wanted = threads().saturating_sub(1);
    let mut state = lock(&POOL.state);
    state.jobs.push(job.clone());
    while state.workers < wanted {
        let id = state.workers;
        let spawned = std::thread::Builder::new()
            .name(format!("catallax-ahead-{id}"))
            .stack_size(WORKER_STACK)
            .spawn(move || work(id));
        if spawned.is_err() {
            break;
        }
        state.workers += 1;
    }
    let staffed = state.workers != 0;
    let wake = items.min(state.workers);
    drop(state);
    for _ in 0..wake {
        POOL.wake.notify_one();
    }
    staffed
}

fn withdraw(job: &Arc<Job>) {
    job.cancelled.store(true, Ordering::Relaxed);
    lock(&POOL.state)
        .jobs
        .retain(|other| !Arc::ptr_eq(other, job));
}

fn work(id: usize) {
    loop {
        let (job, index) = {
            let mut state = lock(&POOL.state);
            loop {
                state.jobs.retain(|job| !job.exhausted());
                // The innermost open loop is the one the writer is in now.
                let taken = (id < threads().saturating_sub(1))
                    .then(|| {
                        state
                            .jobs
                            .iter()
                            .rev()
                            .find_map(|job| job.take().map(|index| (job.clone(), index)))
                    })
                    .flatten();
                if let Some(taken) = taken {
                    break taken;
                }
                state = POOL
                    .wake
                    .wait(state)
                    .unwrap_or_else(PoisonError::into_inner);
            }
        };
        run_item(&job, index);
        job.finish();
    }
}

fn run_item(job: &Arc<Job>, index: usize) {
    struct Leave;
    impl Drop for Leave {
        fn drop(&mut self) {
            ENGAGED.set(false);
            IN_DECISION.set(false);
            RECORDING.set(false);
            ITEM.with(|item| {
                if let Ok(mut item) = item.try_borrow_mut() {
                    *item = None;
                }
            });
            WORKER.with(|worker| {
                if let Ok(mut worker) = worker.try_borrow_mut() {
                    *worker = None;
                }
            });
            crate::operation::leave_worker();
        }
    }
    crate::operation::enter_worker(job.reduced_angle_series, job.reuse_enclosures, &job.proofs);
    ITEM.with(|item| *item.borrow_mut() = Some((job.clone(), index)));
    WORKER.with(|worker| {
        *worker.borrow_mut() = Some(Worker {
            roots: BTreeMap::new(),
            shared: job.roots.clone(),
            inherited: job.proofs.clone(),
            pis: BTreeMap::new(),
            recording: None,
            deciding: None,
            epoch: job.epoch,
        });
    });
    ENGAGED.set(true);
    let _leave = Leave;
    let mut dec = Decider::new(job.policy);
    dec.set_deadline(job.deadline);
    // A panic in an item run ahead is the writer's to meet, if it is real.
    let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        (job.body)(index, &mut dec);
    }));
}

/// One loop's place in the source. It remembers how long the loop's recent
/// runs took, so that a loop that is usually quick wakes no worker. The
/// memory only decides when workers start; it cannot change a result.
pub struct AheadSite {
    /// A running average, in nanoseconds.
    recent: AtomicU64,
    /// For a loop of whole items: what an item took with the workers' help,
    /// as a share of what it takes the writer alone, in 1024ths. A running
    /// average; zero while nothing has been measured.
    helped: AtomicU64,
    /// How many times a loop here was worth handing over.
    asked: AtomicU64,
}

impl AheadSite {
    /// A site with no history.
    pub const fn new() -> Self {
        Self {
            recent: AtomicU64::new(0),
            helped: AtomicU64::new(0),
            asked: AtomicU64::new(0),
        }
    }
}

impl Default for AheadSite {
    fn default() -> Self {
        Self::new()
    }
}

/// A loop is worth waking workers for once it, or its recent runs on average,
/// took this long. Below it the workers cost the writer more than they save:
/// waking them is a system call each, they contend with it for the reference
/// counts of the expressions they share, and only a decision that takes a
/// quarter of a millisecond or more is worth a record.
const SLOW_LOOP_NANOS: u64 = 20_000_000;

/// A request under a work allowance wakes no worker before it has done this
/// much work. Most requests are smaller than this and finish in milliseconds.
const WARM_UP_WORK: u64 = 20_000;

/// A loop whose later items workers run ahead of the writer.
///
/// The writer makes one of these before a loop of independent proofs and
/// calls [`RunAhead::at`] as it starts each item. Dropping it withdraws the
/// items no worker has finished. It changes no answer, refusal or work count.
pub struct RunAhead {
    site: &'static AheadSite,
    len: usize,
    policy: BudgetPolicy,
    deadline: Option<Instant>,
    /// When the loop began, on a thread that may use workers.
    began: Option<Instant>,
    job: Option<Arc<Job>>,
}

impl RunAhead {
    /// A loop of `len` items. With one thread, on a worker thread, or unless
    /// [`set_run_ahead_loops`] asked for it, this does nothing and the loop
    /// runs as it always has.
    pub fn new(site: &'static AheadSite, dec: &Decider, len: usize) -> Self {
        let eligible = len >= 2
            && threads() >= 2
            && loops_go_to_workers()
            && !is_worker()
            && !crate::decider::test_clock_installed()
            && !cfg!(feature = "work-profile");
        Self {
            site,
            len,
            policy: *dec.policy(),
            deadline: dec.deadline(),
            began: eligible.then(Instant::now),
            job: None,
        }
    }

    /// The writer is starting item `index` itself. Workers leave earlier items
    /// alone from here on. If the loop has proved expensive and no worker has
    /// been given its items yet, `make` builds the item body and workers start
    /// on the items after `index`.
    ///
    /// The body gets an item's index and a fresh decider with the writer's
    /// policy. It must own what it uses, do the item on its own copy of the
    /// model and keep nothing: only its decisions are recorded.
    pub fn at<F>(&mut self, index: usize, make: impl FnOnce() -> F)
    where
        F: Fn(usize, &mut Decider) + Send + Sync + 'static,
    {
        if let Some(job) = &self.job {
            job.floor.fetch_max(index, Ordering::Relaxed);
            return;
        }
        let Some(began) = self.began else {
            return;
        };
        if index + 1 >= self.len || crate::operation::exhausted() {
            return;
        }
        let wait = wait_for_workers();
        if !wait
            && crate::operation::operation_budget_active()
            && crate::operation::used() < WARM_UP_WORK
        {
            return;
        }
        // Reading the clock costs as much as a cheap item; in a long loop,
        // look at it only now and then.
        if !wait
            && self.site.recent.load(Ordering::Relaxed) < SLOW_LOOP_NANOS
            && ((index >= 8 && !index.is_multiple_of(8))
                || began.elapsed().as_nanos() < u128::from(SLOW_LOOP_NANOS))
        {
            return;
        }
        let epoch = EPOCH.fetch_add(1, Ordering::Relaxed) + 1;
        let table = Arc::new(Table::new(epoch));
        let mut tables = WRITER.with(|tables| tables.borrow().clone());
        tables.push(table.clone());
        let (reduced_angle_series, reuse_enclosures) = crate::operation::preferences();
        let job = Arc::new(Job {
            body: Box::new(make()),
            len: self.len,
            next: AtomicUsize::new(index + 1),
            floor: AtomicUsize::new(index),
            cancelled: AtomicBool::new(false),
            policy: self.policy,
            deadline: self.deadline,
            reduced_angle_series,
            reuse_enclosures,
            epoch,
            tables,
            proofs: crate::operation::proofs(),
            roots: Arc::default(),
            finished: Mutex::new(index + 1),
            all_finished: Condvar::new(),
        });
        WRITER.with(|tables| tables.borrow_mut().push(table));
        ENGAGED.set(true);
        if submit(&job, self.len - index - 1) && wait {
            let mut finished = lock(&job.finished);
            while *finished < self.len {
                finished = job
                    .all_finished
                    .wait(finished)
                    .unwrap_or_else(PoisonError::into_inner);
            }
        }
        self.job = Some(job);
    }
}

impl Drop for RunAhead {
    fn drop(&mut self) {
        if let Some(began) = self.began {
            let spent = u64::try_from(began.elapsed().as_nanos()).unwrap_or(u64::MAX);
            let recent = self.site.recent.load(Ordering::Relaxed);
            self.site
                .recent
                .store(recent / 2 + spent / 2, Ordering::Relaxed);
        }
        let Some(job) = self.job.take() else {
            return;
        };
        withdraw(&job);
        let own = job.tables.last().expect("a job has its own table");
        WRITER.with(|tables| {
            let mut tables = tables.borrow_mut();
            tables.retain(|table| !Arc::ptr_eq(table, own));
            ENGAGED.set(!tables.is_empty());
        });
    }
}

/// A [`RunAhead`] for the loop at this place in the source, with that place's
/// own [`AheadSite`]. Takes the writer's decider and the number of items.
#[macro_export]
macro_rules! run_ahead {
    ($dec:expr, $len:expr) => {{
        static SITE: $crate::AheadSite = $crate::AheadSite::new();
        $crate::RunAhead::new(&SITE, $dec, $len)
    }};
}

// ---------------------------------------------------------------------------
// Loops whose items other threads compute whole
// ---------------------------------------------------------------------------

/// An item is worth another thread once the writer's items take this long.
/// Below it, waking a worker and replaying its tape cost what the item does.
const WORTH_AN_ITEM_NANOS: u64 = 200_000;

/// The items left are worth handing over once they add up to this.
const WORTH_A_HAND_OVER_NANOS: u64 = 2_000_000;

/// A loop keeps handing over while its items have taken at most this share of
/// the writer's own time for them, in 1024ths. On a machine whose cores are
/// all busy the workers add nothing, and the loop stays with the writer.
const WORTH_THE_HELP: u64 = 920;

/// A loop that has stopped handing over tries again once in this many times.
const TRY_AGAIN: u64 = 16;

enum ValueSlot<T> {
    /// Nobody has started the item.
    Free,
    /// A worker began it at this instant.
    Taken(Instant),
    /// A worker finished it. `None` when its record cannot be used.
    Done(Option<(T, Charges, DecisionStats)>),
    /// The writer computes it, or already has its value.
    Writers,
}

struct Values<T> {
    slots: Mutex<Vec<ValueSlot<T>>>,
    settled: Condvar,
}

/// A loop whose items are computations of values from shared inputs. The
/// writer visits the items in order, and other threads compute later ones
/// whole. Each item's value and charges are taken once: from a worker when
/// its record is usable and fits the allowance, otherwise from the writer's
/// own run. The writer waits a while for an item a worker has begun before it
/// runs the item itself. No answer, refusal or work count changes.
///
/// An item must depend on nothing but its index and values every thread
/// shares, and must keep nothing: no model, cache or decider of the writer's
/// may steer it or be changed by it. Its value and what it charges are then
/// all it does. A worker's item that meets a charged cache is not adopted;
/// the writer computes that one itself.
pub struct ValuesAhead<T> {
    site: &'static AheadSite,
    len: usize,
    policy: BudgetPolicy,
    deadline: Option<Instant>,
    /// When the loop began, on a thread that may use workers.
    began: Option<Instant>,
    /// Items the writer has reached.
    reached: usize,
    /// What one item took the writer in this loop, once measured.
    each: Option<u64>,
    /// When the loop was handed over, at which item, and what an item was
    /// taking the writer alone then.
    handed: Option<(Instant, usize, u64)>,
    /// The loop was worth handing over and was kept, because handing over
    /// here has not been paying. Decided once for the loop.
    kept: bool,
    job: Option<(Arc<Job>, Arc<Values<T>>)>,
}

impl<T: Send + Sync + 'static> ValuesAhead<T> {
    /// A loop of `len` items. With one thread or on a worker thread this does
    /// nothing, and the loop runs as it always has.
    pub fn new(site: &'static AheadSite, dec: &Decider, len: usize) -> Self {
        let eligible = len >= 2
            && threads() >= 2
            && !is_worker()
            && !crate::decider::test_clock_installed()
            && !cfg!(feature = "work-profile");
        Self {
            site,
            len,
            policy: *dec.policy(),
            deadline: dec.deadline(),
            began: eligible.then(Instant::now),
            reached: 0,
            each: None,
            handed: None,
            kept: false,
            job: None,
        }
    }

    /// The writer is at item `index`. Returns the item's value when another
    /// thread computed it; its charges and statistics are then the writer's,
    /// as if it had computed the item here. Returns `None` when the writer
    /// must compute the item itself.
    ///
    /// `make` builds the item computation for the other threads, once, when
    /// the loop has proved long enough to hand over. The computation gets an
    /// item's index and a fresh decider with the writer's policy, and must own
    /// what it uses.
    pub fn take<F>(
        &mut self,
        index: usize,
        dec: &mut Decider,
        make: impl FnOnce() -> F,
    ) -> Option<T>
    where
        F: Fn(usize, &mut Decider) -> T + Send + Sync + 'static,
    {
        self.reached = index + 1;
        if self.job.is_none() {
            self.hand_over(index, make);
            return None;
        }
        let deadline = self.deadline;
        let (job, values) = self.job.as_ref()?;
        job.floor.fetch_max(index, Ordering::Relaxed);
        let mut slots = lock(&values.slots);
        let wait = wait_for_workers();
        let mut patience = None;
        loop {
            match &slots[index] {
                ValueSlot::Free | ValueSlot::Writers => {
                    slots[index] = ValueSlot::Writers;
                    return None;
                }
                ValueSlot::Taken(since) => {
                    // The worker gets as long again as it has had, and no
                    // longer than the writer's own deadline. Then the writer
                    // computes the item itself: a worker the machine has set
                    // aside must not hold the request up, and past its
                    // deadline the writer has a refusal to give.
                    let now = Instant::now();
                    let patient = *patience.get_or_insert_with(|| now + since.elapsed());
                    let until = match (wait, deadline) {
                        (true, None) => None,
                        (true, Some(deadline)) => Some(deadline),
                        (false, None) => Some(patient),
                        (false, Some(deadline)) => Some(patient.min(deadline)),
                    };
                    let Some(until) = until else {
                        slots = values
                            .settled
                            .wait(slots)
                            .unwrap_or_else(PoisonError::into_inner);
                        continue;
                    };
                    let left = until.saturating_duration_since(now);
                    if left.is_zero() {
                        slots[index] = ValueSlot::Writers;
                        return None;
                    }
                    slots = values
                        .settled
                        .wait_timeout(slots, left)
                        .unwrap_or_else(PoisonError::into_inner)
                        .0;
                }
                ValueSlot::Done(_) => {
                    let ValueSlot::Done(record) =
                        std::mem::replace(&mut slots[index], ValueSlot::Writers)
                    else {
                        unreachable!("matched above")
                    };
                    drop(slots);
                    let (value, charges, stats) = record?;
                    if !charges.adopt(dec) {
                        return None;
                    }
                    if let Some(own) = Host::stats(dec) {
                        own.absorb(&stats);
                    }
                    return Some(value);
                }
            }
        }
    }

    fn hand_over<F>(&mut self, index: usize, make: impl FnOnce() -> F)
    where
        F: Fn(usize, &mut Decider) -> T + Send + Sync + 'static,
    {
        let Some(began) = self.began else {
            return;
        };
        if index + 1 >= self.len || crate::operation::exhausted() {
            return;
        }
        let wait = wait_for_workers();
        if !wait {
            if self.kept
                || (crate::operation::operation_budget_active()
                    && crate::operation::used() < WARM_UP_WORK)
            {
                return;
            }
            // What the items before this one took the writer each, or what
            // this loop's items took it in its recent runs. Reading the clock
            // costs as much as a cheap item, so look only now and then.
            let each = if index == 0 {
                let recent = self.site.recent.load(Ordering::Relaxed);
                // Old news fades, so that a loop whose items have become
                // quick is measured again before long.
                self.site
                    .recent
                    .store(recent - recent / 8, Ordering::Relaxed);
                recent
            } else if index >= 8 && !index.is_multiple_of(8) {
                return;
            } else {
                let each =
                    u64::try_from(began.elapsed().as_nanos()).unwrap_or(u64::MAX) / index as u64;
                self.each = Some(each);
                each
            };
            let left = (self.len - index - 1) as u64;
            if each < WORTH_AN_ITEM_NANOS || each.saturating_mul(left) < WORTH_A_HAND_OVER_NANOS {
                return;
            }
            // Has handing this loop over been quicker than the writer alone?
            let helped = self.site.helped.load(Ordering::Relaxed);
            let asked = self.site.asked.fetch_add(1, Ordering::Relaxed);
            if helped > WORTH_THE_HELP && !asked.is_multiple_of(TRY_AGAIN) {
                self.kept = true;
                return;
            }
            self.handed = Some((Instant::now(), index, each));
        }
        let item = make();
        let values = Arc::new(Values {
            slots: Mutex::new((0..self.len).map(|_| ValueSlot::Free).collect()),
            settled: Condvar::new(),
        });
        let body = {
            let values = values.clone();
            move |index: usize, dec: &mut Decider| {
                // The slot is settled on every path out, a panic included: a
                // slot left taken would be waited for.
                struct Settle<'a, T>(&'a Values<T>, usize);
                impl<T> Drop for Settle<'_, T> {
                    fn drop(&mut self) {
                        let mut slots = lock(&self.0.slots);
                        if matches!(slots[self.1], ValueSlot::Taken(_)) {
                            slots[self.1] = ValueSlot::Done(None);
                        }
                        drop(slots);
                        self.0.settled.notify_all();
                    }
                }
                {
                    let mut slots = lock(&values.slots);
                    if !matches!(slots[index], ValueSlot::Free) {
                        return;
                    }
                    slots[index] = ValueSlot::Taken(Instant::now());
                }
                let settle = Settle(&values, index);
                let (value, charges) = record(|| item(index, dec));
                let record = charges.map(|charges| (value, charges, Decider::stats(dec).clone()));
                let mut slots = lock(&values.slots);
                if matches!(slots[index], ValueSlot::Taken(_)) {
                    slots[index] = ValueSlot::Done(record);
                }
                drop(slots);
                drop(settle);
            }
        };
        let epoch = EPOCH.fetch_add(1, Ordering::Relaxed) + 1;
        let (reduced_angle_series, reuse_enclosures) = crate::operation::preferences();
        let job = Arc::new(Job {
            body: Box::new(body),
            len: self.len,
            next: AtomicUsize::new(index + 1),
            floor: AtomicUsize::new(index),
            cancelled: AtomicBool::new(false),
            policy: self.policy,
            deadline: self.deadline,
            reduced_angle_series,
            reuse_enclosures,
            epoch,
            tables: Vec::new(),
            proofs: crate::operation::proofs(),
            roots: Arc::default(),
            finished: Mutex::new(index + 1),
            all_finished: Condvar::new(),
        });
        if submit(&job, self.len - index - 1) && wait {
            let mut finished = lock(&job.finished);
            while *finished < self.len {
                finished = job
                    .all_finished
                    .wait(finished)
                    .unwrap_or_else(PoisonError::into_inner);
            }
        }
        self.job = Some((job, values));
    }
}

impl<T> Drop for ValuesAhead<T> {
    fn drop(&mut self) {
        // The site remembers what an item takes the writer alone. A loop the
        // workers shared says nothing about that after the hand-over.
        let each = match (&self.job, self.began) {
            (None, Some(began)) if self.reached != 0 => Some(
                u64::try_from(began.elapsed().as_nanos()).unwrap_or(u64::MAX) / self.reached as u64,
            ),
            _ => self.each,
        };
        if let Some(each) = each {
            let recent = self.site.recent.load(Ordering::Relaxed);
            self.site
                .recent
                .store(recent / 2 + each / 2, Ordering::Relaxed);
        }
        // And how the items went with the workers' help, against the writer
        // alone: the next loop here is handed over only while that pays.
        if let Some((since, index, alone)) = self.handed
            && self.reached > index
            && alone != 0
        {
            let each = u64::try_from(since.elapsed().as_nanos()).unwrap_or(u64::MAX)
                / (self.reached - index) as u64;
            let share = (each.saturating_mul(1024) / alone).clamp(16, 16 * 1024);
            let helped = self.site.helped.load(Ordering::Relaxed);
            let helped = if helped == 0 {
                share
            } else {
                (helped * 3 + share) / 4
            };
            self.site.helped.store(helped, Ordering::Relaxed);
        }
        if let Some((job, _)) = self.job.take() {
            withdraw(&job);
        }
    }
}

/// A [`ValuesAhead`] for the loop at this place in the source, with that
/// place's own [`AheadSite`]. Takes the writer's decider and the number of
/// items.
#[macro_export]
macro_rules! values_ahead {
    ($dec:expr, $len:expr) => {{
        static SITE: $crate::AheadSite = $crate::AheadSite::new();
        $crate::ValuesAhead::new(&SITE, $dec, $len)
    }};
}
