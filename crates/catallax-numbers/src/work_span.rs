//! Development-only work and span measurement.
//!
//! Compiled only with the `work-span` feature. The observer replays one
//! request on an ideal machine with unlimited workers and free scheduling, in
//! request work units. Code runs in its normal order; the observer only keeps
//! clocks. Work charged on the current strand advances its clock. Inside a
//! loop marked with [`parallel`], every task begins at the time the loop body
//! reached it, not after the previous task. A read of a tracked resource moves
//! the reader's clock forward to the time its last writer wrote it, so data
//! dependencies, including hits on proof caches, are kept. The span is the
//! finishing time of the whole request on that machine; work divided by span
//! bounds the speed-up of any schedule of the same tasks.
//!
//! Three clocks run side by side, one per [`Model`]. They differ only in which
//! loops run their tasks side by side and which reads wait.
//!
//! Observation never charges work, asks a numeric question or changes an
//! outcome. Measure one request at a time: hooks from every thread are merged.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::{
    borrow::Cow,
    collections::{BTreeMap, HashMap},
    hash::{Hash, Hasher},
    panic::Location,
    sync::{Arc, Mutex, PoisonError},
};

static ACTIVE: AtomicBool = AtomicBool::new(false);
static STATE: Mutex<Option<State>> = Mutex::new(None);

/// What may run side by side on the ideal machine.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Model {
    /// Only the tasks of [`steps`] loops, the recipe's operations, run side
    /// by side: feature-level parallelism alone.
    Steps,
    /// Tasks of [`steps`] and [`parallel`] loops run side by side; every
    /// tracked read waits for its writer. Runs exactly the work the
    /// sequential build ran.
    Independent,
    /// Also runs the candidates of [`search`] loops side by side. Running a
    /// search this way can do more work than the sequential build, which
    /// stops at the first accepted candidate; that extra work is not counted.
    Speculative,
    /// As [`Model::Speculative`], but reads of shared indexes, fact books and
    /// proof caches do not wait. A bound for a kernel whose tasks no longer
    /// share them, before the cost of not sharing them, which is not counted.
    Unshared,
    /// As [`Model::Unshared`], but each task's [`append`] regions, the merges
    /// of its results into the loop's shared lists, leave the task and join
    /// one ordered commit chain per loop run: plan in parallel, append in
    /// order. The task goes on without waiting for the merge.
    Ordered,
}

/// Every model, in report order.
pub const MODELS: [Model; 5] = [
    Model::Steps,
    Model::Independent,
    Model::Speculative,
    Model::Unshared,
    Model::Ordered,
];
const N: usize = MODELS.len();

/// Whether `kind` names shared bookkeeping rather than a constructed entity:
/// a fact book, an index, an arena's length or a proof or bounds cache.
pub fn shared(kind: &str) -> bool {
    kind.starts_with("fact ")
        || kind.contains("cache")
        || kind == "conic supports"
        || kind == EXTENT
}

/// The resource kind of a store's length, one per store copy.
pub const EXTENT: &str = "arena extent";

/// Work on one critical-path segment, attributed to its innermost region.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathRow {
    /// Region path from the outermost region inward, joined by ` > `.
    pub region: String,
    /// Work units of that region on the critical path.
    pub work: u64,
}

/// Work charged inside one region, excluding its nested regions, and the
/// time the region was open.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionRow {
    /// Region path, joined by ` > `.
    pub region: String,
    /// Work units charged directly in the region.
    pub work: u64,
    /// Per model, the summed time from entering the region to leaving it on
    /// its own strand, nested regions included. Regions that ran side by
    /// side in one loop overlap, so their sum can exceed the loop's span.
    pub span: [u64; N],
}

/// One loop whose tasks may run side by side.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LoopRow {
    /// Region path of the loop, joined by ` > `.
    pub region: String,
    /// Whether the loop stops at its first accepted candidate.
    pub search: bool,
    /// Times the loop ran.
    pub runs: u64,
    /// Tasks over all runs.
    pub tasks: u64,
    /// Work units over all runs, including work between tasks.
    pub work: u64,
    /// Per model, the sum over runs of the time from loop entry to loop exit.
    pub span: [u64; N],
    /// Largest work of one task.
    pub largest_task: u64,
}

/// Reads that moved a strand's clock forward, by resource kind.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DependencyRow {
    /// Resource kind, such as an arena type or a proof cache.
    pub kind: String,
    /// Whether [`Model::Unshared`] ignores this kind.
    pub shared: bool,
    /// Reads that found a later write than the reader's clock, in
    /// [`Model::Speculative`].
    pub waits: u64,
    /// Work units the reader's clock moved over all those reads.
    pub delay: u64,
    /// Waits whose writer ran in a task still open beside the reader's task.
    pub sibling_waits: u64,
}

/// Top-level decider calls from one entry point and calling line. A call
/// outside any other decider call is one indivisible step of its strand.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DecisionRow {
    /// Decider entry point.
    pub label: &'static str,
    /// Calling source location, as `file:line`.
    pub caller: String,
    /// Region path of the largest call.
    pub region: String,
    /// Completed calls.
    pub calls: u64,
    /// Work units of all calls, including nested calls.
    pub work: u64,
    /// Work units of the largest call.
    pub largest: u64,
}

/// One model's span and critical path.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SpanRow {
    /// The model.
    pub model: Option<Model>,
    /// Finishing time on the ideal machine, in work units.
    pub span: u64,
    /// The critical path from first to last segment, adjacent equal regions merged.
    pub path: Vec<PathRow>,
}

/// The result of one observed request.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct WorkSpanReport {
    /// Work units charged while observing.
    pub work: u64,
    /// One row per model, in [`MODELS`] order.
    pub spans: Vec<SpanRow>,
    /// Work per region.
    pub regions: Vec<RegionRow>,
    /// Parallel and search loops.
    pub loops: Vec<LoopRow>,
    /// Dependencies that delayed a strand.
    pub dependencies: Vec<DependencyRow>,
    /// Top-level decider calls by calling line, largest single call first.
    pub decisions: Vec<DecisionRow>,
    /// Tracked reads, for scale.
    pub reads: u64,
    /// Tracked writes.
    pub writes: u64,
}

/// Observes `run` and reports its work and span. Hooks from every thread
/// count, so a request answered on its own thread is observed whole; run no
/// other kernel work in the process meanwhile. Nested observations are rejected.
pub fn with_work_span<T>(run: impl FnOnce() -> T) -> (T, WorkSpanReport) {
    observe(None, run)
}

/// As [`with_work_span`], but counts only hooks on the calling thread, so
/// other threads' work, such as concurrent tests, cannot leak in.
pub fn with_work_span_on_this_thread<T>(run: impl FnOnce() -> T) -> (T, WorkSpanReport) {
    observe(Some(std::thread::current().id()), run)
}

fn observe<T>(only: Option<std::thread::ThreadId>, run: impl FnOnce() -> T) -> (T, WorkSpanReport) {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            ACTIVE.store(false, Ordering::SeqCst);
            *lock() = None;
        }
    }
    {
        let mut state = lock();
        assert!(state.is_none(), "work-span observations cannot nest");
        let mut fresh = State::new();
        fresh.only = only;
        *state = Some(fresh);
    }
    let _reset = Reset;
    ACTIVE.store(true, Ordering::SeqCst);
    let result = run();
    ACTIVE.store(false, Ordering::SeqCst);
    let report = lock().take().map(State::report).unwrap_or_default();
    (result, report)
}

/// A stable 64-bit key for a tracked resource.
pub fn key(value: &impl Hash) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

/// Notes a read of `(kind, key)`: the reader cannot run before its last writer.
#[inline]
pub fn read(kind: &'static str, key: u64) {
    read_at(kind, 0, key);
}

/// Notes a write of `(kind, key)` at the current strand's time.
#[inline]
pub fn write(kind: &'static str, key: u64) {
    write_at(kind, 0, key);
}

/// Notes a read of entry `index` of the store or cache named by `owner`, an
/// [`Instance`] number. Distinct owners never share a resource.
#[inline]
pub fn read_at(kind: &'static str, owner: u64, index: u64) {
    if ACTIVE.load(Ordering::Relaxed) {
        with_state(|state| state.read((kind, owner, index)));
    }
}

/// Notes a write of entry `index` of the store or cache named by `owner`.
#[inline]
pub fn write_at(kind: &'static str, owner: u64, index: u64) {
    if ACTIVE.load(Ordering::Relaxed) {
        with_state(|state| state.write((kind, owner, index)));
    }
}

/// Gives `to` the last writer of `from`: a copy of a store starts with its
/// original's history.
pub fn copy_at(kind: &'static str, from: (u64, u64), to: (u64, u64)) {
    if ACTIVE.load(Ordering::Relaxed) {
        with_state(|state| state.copy((kind, from.0, from.1), (kind, to.0, to.1)));
    }
}

/// A number that names one store, cache or copy of a store for the observer.
/// Every value, including every clone, gets a fresh number, so independent
/// state never shares a resource.
pub struct Instance(u64);

static INSTANCES: AtomicU64 = AtomicU64::new(1);

impl Instance {
    /// A fresh number.
    pub fn new() -> Self {
        Self(INSTANCES.fetch_add(1, Ordering::Relaxed))
    }

    /// The number.
    pub fn id(&self) -> u64 {
        self.0
    }
}

impl Default for Instance {
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for Instance {
    fn clone(&self) -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for Instance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Instance")
    }
}

/// The observer's view of a keyed table that is copied with its store. Every
/// write is a new, never reused generation, and each entry remembers the
/// generation it holds. A copy keeps the generations it inherited, so a later
/// write to either copy, including a rewrite of an inherited entry, is
/// invisible to the other. Use one per kind.
#[derive(Clone, Default, Debug)]
pub struct Entries {
    generations: im::OrdMap<u64, u64>,
}

impl Entries {
    /// Notes a read of entry `index`.
    pub fn read(&self, kind: &'static str, index: u64) {
        if ACTIVE.load(Ordering::Relaxed)
            && let Some(&generation) = self.generations.get(&index)
        {
            read_at(kind, generation, index);
        }
    }

    /// Notes a write of entry `index` by this copy.
    pub fn write(&mut self, kind: &'static str, index: u64) {
        let generation = Instance::new().id();
        self.generations.insert(index, generation);
        write_at(kind, generation, index);
    }
}

/// Writes `(kind, key)` when dropped, at the time the owner finished.
#[must_use = "the write happens when dropped"]
pub struct WriteOnDrop(Option<(&'static str, u64, u64)>);

/// Defers a write to the end of the enclosing scope.
pub fn write_on_drop(kind: &'static str, key: u64) -> WriteOnDrop {
    write_on_drop_at(kind, 0, key)
}

/// Defers a write of an owned entry to the end of the enclosing scope.
pub fn write_on_drop_at(kind: &'static str, owner: u64, index: u64) -> WriteOnDrop {
    WriteOnDrop(
        ACTIVE
            .load(Ordering::Relaxed)
            .then_some((kind, owner, index)),
    )
}

impl Drop for WriteOnDrop {
    fn drop(&mut self) {
        if let Some((kind, owner, index)) = self.0 {
            write_at(kind, owner, index);
        }
    }
}

/// A named series region. Dropping it closes the region.
#[must_use = "a region closes when dropped"]
pub struct Region(bool);

/// Opens a named region; work inside it is attributed to it.
pub fn region(name: impl Into<Cow<'static, str>>) -> Region {
    if !ACTIVE.load(Ordering::Relaxed) {
        return Region(false);
    }
    let name = name.into();
    Region(with_state(|state| state.enter_region(&name)).is_some())
}

impl Drop for Region {
    fn drop(&mut self) {
        if self.0 {
            with_state(State::leave_region);
        }
    }
}

/// An ordered merge of a task's results into its loop's shared lists.
#[must_use = "an append region closes when dropped"]
pub struct Append(bool);

/// Opens an append region. Only [`Model::Ordered`] treats it specially.
pub fn append() -> Append {
    if !ACTIVE.load(Ordering::Relaxed) {
        return Append(false);
    }
    Append(with_state(State::enter_append).is_some())
}

impl Drop for Append {
    fn drop(&mut self) {
        if self.0 {
            with_state(State::leave_append);
        }
    }
}

/// A loop whose tasks may run side by side. Dropping it joins every task.
#[must_use = "a parallel loop joins its tasks when dropped"]
pub struct Parallel(bool);

/// Opens a named loop whose iterations are independent tasks except through
/// tracked resources. The caller must keep every other loop-carried value out
/// of the tasks or track it with [`read`] and [`write`]. Commutative
/// collection of results, such as inserting one entry per task, is allowed.
pub fn parallel(name: &'static str) -> Parallel {
    open(name, Kind::Parallel)
}

/// Opens the loop over a recipe's operations. Its tasks are the only ones
/// [`Model::Steps`] runs side by side; otherwise it is a [`parallel`] loop.
pub fn steps(name: &'static str) -> Parallel {
    open(name, Kind::Steps)
}

/// Opens a named loop that tries candidates in order and stops at the first
/// accepted one. Only [`Model::Speculative`] and [`Model::Unshared`] run its
/// tasks side by side.
pub fn search(name: &'static str) -> Parallel {
    open(name, Kind::Search)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Steps,
    Parallel,
    Search,
}

fn open(name: &'static str, kind: Kind) -> Parallel {
    if !ACTIVE.load(Ordering::Relaxed) {
        return Parallel(false);
    }
    Parallel(with_state(|state| state.enter_loop(name, kind)).is_some())
}

impl Parallel {
    /// Starts one task of this loop. Dropping it ends the task.
    pub fn task(&self) -> Task {
        Task(self.0 && with_state(State::enter_task).is_some())
    }
}

impl Drop for Parallel {
    fn drop(&mut self) {
        if self.0 {
            with_state(State::leave_loop);
        }
    }
}

/// One task of a [`Parallel`] loop.
#[must_use = "a task ends when dropped"]
pub struct Task(bool);

impl Drop for Task {
    fn drop(&mut self) {
        if self.0 {
            with_state(State::leave_task);
        }
    }
}

pub(crate) fn advance(units: u64) {
    if ACTIVE.load(Ordering::Relaxed) {
        with_state(|state| state.advance(units));
    }
}

pub(crate) fn decision_enter() {
    if ACTIVE.load(Ordering::Relaxed) {
        with_state(|state| {
            state.decision_depth += 1;
            if state.decision_depth == 1 {
                state.decision_start = state.work;
            }
        });
    }
}

pub(crate) fn decision_leave(label: &'static str, caller: &'static Location<'static>) {
    if ACTIVE.load(Ordering::Relaxed) {
        with_state(|state| state.decision_leave(label, caller));
    }
}

fn lock() -> std::sync::MutexGuard<'static, Option<State>> {
    STATE.lock().unwrap_or_else(PoisonError::into_inner)
}

fn with_state<T>(f: impl FnOnce(&mut State) -> T) -> Option<T> {
    // A run-ahead worker computes on its own copy and its own counter; only
    // the writer's work and reads belong to the request.
    if crate::ahead::is_worker() {
        return None;
    }
    let mut guard = lock();
    let state = guard.as_mut()?;
    if state
        .only
        .is_some_and(|thread| thread != std::thread::current().id())
    {
        return None;
    }
    Some(f(state))
}

/// A persistent list of critical-path segments, newest first.
struct Segment {
    previous: Option<Arc<Segment>>,
    label: usize,
    work: u64,
}

impl Drop for Segment {
    fn drop(&mut self) {
        let mut previous = self.previous.take();
        while let Some(segment) = previous {
            match Arc::try_unwrap(segment) {
                Ok(mut segment) => previous = segment.previous.take(),
                Err(_) => break,
            }
        }
    }
}

/// A strand's position: its clock, the path that led to it, and work since
/// the last segment boundary.
#[derive(Clone, Default)]
struct Clock {
    time: u64,
    path: Option<Arc<Segment>>,
    pending: u64,
}

impl Clock {
    fn flush(&mut self, label: usize) {
        if self.pending > 0 {
            self.path = Some(Arc::new(Segment {
                previous: self.path.take(),
                label,
                work: self.pending,
            }));
            self.pending = 0;
        }
    }

    fn advance(&mut self, units: u64) {
        self.time += units;
        self.pending += units;
    }
}

type Clocks = [Clock; N];

struct LoopFrame {
    label: usize,
    kind: Kind,
    /// Per model, where the next task starts: the loop body's own clock.
    base: Clocks,
    /// Per model, the latest finishing task so far.
    end: Clocks,
    entered: [u64; N],
    entered_work: u64,
    task_work: Option<u64>,
    tasks: u64,
    largest: u64,
    /// The [`Model::Ordered`] commit chain of this loop run.
    commit: Clock,
}

impl LoopFrame {
    fn side_by_side(&self, model: usize) -> bool {
        match MODELS[model] {
            Model::Steps => self.kind == Kind::Steps,
            Model::Independent => self.kind != Kind::Search,
            Model::Speculative | Model::Unshared | Model::Ordered => true,
        }
    }
}

type Resource = (&'static str, u64, u64);

#[derive(Clone)]
struct Written {
    clocks: Clocks,
    /// Depth of open loops and the task number in each when written.
    tasks: Arc<[(usize, u64)]>,
}

struct Decisions {
    label: &'static str,
    caller: &'static Location<'static>,
    calls: u64,
    work: u64,
    largest: u64,
    region: usize,
}

struct State {
    now: Clocks,
    work: u64,
    labels: Vec<(Option<usize>, String)>,
    label_ids: HashMap<(Option<usize>, String), usize>,
    stack: Vec<usize>,
    /// Per open region, each model's clock when it was entered.
    entered: Vec<[u64; N]>,
    region_span: BTreeMap<usize, [u64; N]>,
    loops: Vec<LoopFrame>,
    task_ids: Vec<(usize, u64)>,
    task_counter: u64,
    written: HashMap<Resource, Written>,
    region_work: BTreeMap<usize, u64>,
    loop_rows: BTreeMap<usize, LoopRow>,
    dependencies: BTreeMap<&'static str, DependencyRow>,
    decision_depth: u32,
    decision_start: u64,
    decisions: HashMap<(&'static str, &'static str, u32), Decisions>,
    reads: u64,
    writes: u64,
    /// Open append regions; the outermost saves the task's ordered clock.
    appends: Vec<Option<Clock>>,
    /// The only thread whose hooks count, when restricted.
    only: Option<std::thread::ThreadId>,
}

const ORDERED: usize = 4;

impl State {
    fn new() -> Self {
        let mut state = Self {
            now: Default::default(),
            work: 0,
            labels: Vec::new(),
            label_ids: HashMap::new(),
            stack: Vec::new(),
            entered: Vec::new(),
            region_span: BTreeMap::new(),
            loops: Vec::new(),
            task_ids: Vec::new(),
            task_counter: 0,
            written: HashMap::new(),
            region_work: BTreeMap::new(),
            loop_rows: BTreeMap::new(),
            dependencies: BTreeMap::new(),
            decision_depth: 0,
            decision_start: 0,
            decisions: HashMap::new(),
            reads: 0,
            writes: 0,
            appends: Vec::new(),
            only: None,
        };
        debug_assert_eq!(MODELS[ORDERED], Model::Ordered);
        let root = state.intern(None, "request");
        state.stack.push(root);
        state.entered.push([0; N]);
        state
    }

    fn label(&self) -> usize {
        *self
            .stack
            .last()
            .expect("the request region is never closed")
    }

    fn intern(&mut self, parent: Option<usize>, name: &str) -> usize {
        let key = (parent, name.to_owned());
        if let Some(&id) = self.label_ids.get(&key) {
            return id;
        }
        let id = self.labels.len();
        self.labels.push(key.clone());
        self.label_ids.insert(key, id);
        id
    }

    fn name(&self, mut label: usize) -> String {
        let mut parts = Vec::new();
        loop {
            let (parent, name) = &self.labels[label];
            parts.push(name.as_str());
            match parent {
                Some(parent) => label = *parent,
                None => break,
            }
        }
        parts.reverse();
        parts.join(" > ")
    }

    fn flush(&mut self) {
        let label = self.label();
        for clock in &mut self.now {
            clock.flush(label);
        }
    }

    /// Work, reads and writes inside a loop but outside its tasks happen
    /// after every task so far: they are a join, never overlapped.
    fn settle(&mut self) {
        let Some(frame) = self.loops.last() else {
            return;
        };
        if frame.task_work.is_some() {
            return;
        }
        let label = self.label();
        for (now, end) in self.now.iter_mut().zip(&frame.end) {
            if end.time > now.time {
                now.flush(label);
                *now = end.clone();
            }
        }
    }

    fn advance(&mut self, units: u64) {
        self.settle();
        self.work += units;
        for clock in &mut self.now {
            clock.advance(units);
        }
        *self.region_work.entry(self.label()).or_default() += units;
    }

    fn read(&mut self, resource: Resource) {
        self.reads += 1;
        self.settle();
        let kind = resource.0;
        let Some(written) = self.written.get(&resource) else {
            return;
        };
        let shared = shared(kind);
        let label = self.label();
        let mut moved = None;
        for (model, (now, clock)) in self.now.iter_mut().zip(&written.clocks).enumerate() {
            let ignored = shared && matches!(MODELS[model], Model::Unshared | Model::Ordered);
            if clock.time <= now.time || ignored {
                continue;
            }
            if MODELS[model] == Model::Speculative {
                moved = Some(clock.time - now.time);
            }
            now.flush(label);
            *now = clock.clone();
        }
        let Some(delay) = moved else {
            return;
        };
        // A writer in a task of an open loop that differs from the reader's
        // own task in that loop ran beside the reader on the ideal machine.
        let sibling = written
            .tasks
            .iter()
            .zip(&self.task_ids)
            .any(|(written, reading)| written.0 == reading.0 && written.1 != reading.1);
        let row = self
            .dependencies
            .entry(kind)
            .or_insert_with(|| DependencyRow {
                kind: kind.to_owned(),
                shared,
                ..DependencyRow::default()
            });
        row.waits += 1;
        row.delay += delay;
        row.sibling_waits += u64::from(sibling);
    }

    // Every resource names one entry of one store, cache or copy, so the
    // last write in program order is the one a later read sees.
    fn write(&mut self, resource: Resource) {
        self.writes += 1;
        self.settle();
        self.flush();
        self.written.insert(
            resource,
            Written {
                clocks: self.now.clone(),
                tasks: self.task_ids.as_slice().into(),
            },
        );
    }

    fn copy(&mut self, from: Resource, to: Resource) {
        if let Some(written) = self.written.get(&from).cloned() {
            self.written.insert(to, written);
        } else {
            self.written.remove(&to);
        }
    }

    fn enter_region(&mut self, name: &str) {
        self.settle();
        self.flush();
        let label = self.intern(Some(self.label()), name);
        self.stack.push(label);
        self.entered
            .push(self.now.each_ref().map(|clock| clock.time));
    }

    fn leave_region(&mut self) {
        self.settle();
        self.flush();
        if self.stack.len() > 1 {
            let label = self.stack.pop().expect("an open region");
            let entered = self.entered.pop().expect("an open region");
            let span = self.region_span.entry(label).or_insert([0; N]);
            for ((total, now), entered) in span.iter_mut().zip(&self.now).zip(entered) {
                *total += now.time.saturating_sub(entered);
            }
        }
    }

    fn enter_loop(&mut self, name: &'static str, kind: Kind) {
        self.settle();
        self.enter_region(name);
        let commit = self.now[ORDERED].clone();
        self.loops.push(LoopFrame {
            commit,
            label: self.label(),
            kind,
            base: self.now.clone(),
            end: self.now.clone(),
            entered: self.now.each_ref().map(|clock| clock.time),
            entered_work: self.work,
            task_work: None,
            tasks: 0,
            largest: 0,
        });
    }

    fn enter_task(&mut self) {
        self.flush();
        self.task_counter += 1;
        let depth = self.loops.len();
        let counter = self.task_counter;
        let (now, work) = (self.now.clone(), self.work);
        let Some(frame) = self.loops.last_mut() else {
            return;
        };
        frame.base = now;
        frame.task_work = Some(work);
        frame.tasks += 1;
        self.task_ids.push((depth, counter));
    }

    fn leave_task(&mut self) {
        self.flush();
        let work = self.work;
        let Some(frame) = self.loops.last_mut() else {
            return;
        };
        if let Some(start) = frame.task_work.take() {
            frame.largest = frame.largest.max(work - start);
            self.task_ids.pop();
        }
        for model in 0..N {
            let now = &mut self.now[model];
            if now.time > frame.end[model].time {
                frame.end[model] = now.clone();
            }
            // A side-by-side task leaves the loop body where it found it;
            // a serial one moves it on to its own finish.
            if frame.side_by_side(model) {
                *now = frame.base[model].clone();
            }
        }
    }

    fn leave_loop(&mut self) {
        self.settle();
        self.flush();
        let Some(frame) = self.loops.pop() else {
            return;
        };
        for (now, end) in self.now.iter_mut().zip(&frame.end) {
            if end.time > now.time {
                *now = end.clone();
            }
        }
        if frame.commit.time > self.now[ORDERED].time {
            self.now[ORDERED] = frame.commit.clone();
        }
        let row = self.loop_rows.entry(frame.label).or_default();
        row.search = frame.kind == Kind::Search;
        row.runs += 1;
        row.tasks += frame.tasks;
        row.work += self.work - frame.entered_work;
        for (span, (now, entered)) in row.span.iter_mut().zip(self.now.iter().zip(frame.entered)) {
            *span += now.time - entered;
        }
        row.largest_task = row.largest_task.max(frame.largest);
        self.leave_region();
    }

    fn enter_append(&mut self) {
        let in_task = self
            .loops
            .last()
            .is_some_and(|frame| frame.task_work.is_some());
        if !self.appends.is_empty() || !in_task {
            self.appends.push(None);
            return;
        }
        self.flush();
        let frame = self.loops.last().expect("an open task has a loop");
        let task = self.now[ORDERED].clone();
        if frame.commit.time > task.time {
            self.now[ORDERED] = frame.commit.clone();
        }
        self.appends.push(Some(task));
    }

    fn leave_append(&mut self) {
        let Some(saved) = self.appends.pop() else {
            return;
        };
        let Some(task) = saved else {
            return;
        };
        self.flush();
        if let Some(frame) = self.loops.last_mut() {
            frame.commit = self.now[ORDERED].clone();
        }
        self.now[ORDERED] = task;
    }

    fn decision_leave(&mut self, label: &'static str, caller: &'static Location<'static>) {
        if self.decision_depth == 0 {
            return;
        }
        self.decision_depth -= 1;
        if self.decision_depth > 0 {
            return;
        }
        let work = self.work - self.decision_start;
        let region = self.label();
        let row = self
            .decisions
            .entry((label, caller.file(), caller.line()))
            .or_insert(Decisions {
                label,
                caller,
                calls: 0,
                work: 0,
                largest: 0,
                region,
            });
        row.calls += 1;
        row.work += work;
        if work > row.largest {
            row.largest = work;
            row.region = region;
        }
    }

    fn path(&self, clock: &Clock) -> Vec<PathRow> {
        let mut path = Vec::new();
        let mut segment = clock.path.clone();
        while let Some(current) = segment {
            path.push((current.label, current.work));
            segment = current.previous.clone();
        }
        path.reverse();
        let mut merged: Vec<PathRow> = Vec::new();
        for (label, work) in path {
            let region = self.name(label);
            match merged.last_mut() {
                Some(last) if last.region == region => last.work += work,
                _ => merged.push(PathRow { region, work }),
            }
        }
        merged
    }

    fn report(mut self) -> WorkSpanReport {
        // Close anything left open by unwinding, innermost first.
        while !self.loops.is_empty() {
            if self
                .loops
                .last()
                .is_some_and(|frame| frame.task_work.is_some())
            {
                self.leave_task();
            }
            self.leave_loop();
        }
        self.flush();
        let spans = MODELS
            .iter()
            .zip(&self.now)
            .map(|(&model, clock)| SpanRow {
                model: Some(model),
                span: clock.time,
                path: self.path(clock),
            })
            .collect();
        let root = self.stack[0];
        self.region_span
            .insert(root, self.now.each_ref().map(|clock| clock.time));
        let mut labels: Vec<_> = self.region_work.keys().copied().collect();
        labels.extend(self.region_span.keys().copied());
        labels.sort_unstable();
        labels.dedup();
        let regions = labels
            .into_iter()
            .map(|label| RegionRow {
                region: self.name(label),
                work: self.region_work.get(&label).copied().unwrap_or(0),
                span: self.region_span.get(&label).copied().unwrap_or([0; N]),
            })
            .collect();
        let loops = self
            .loop_rows
            .iter()
            .map(|(&label, row)| LoopRow {
                region: self.name(label),
                ..row.clone()
            })
            .collect();
        let mut decisions: Vec<_> = self
            .decisions
            .values()
            .map(|row| DecisionRow {
                label: row.label,
                caller: format!("{}:{}", row.caller.file(), row.caller.line()),
                region: self.name(row.region),
                calls: row.calls,
                work: row.work,
                largest: row.largest,
            })
            .collect();
        decisions.sort_by(|a, b| {
            b.largest
                .cmp(&a.largest)
                .then_with(|| a.caller.cmp(&b.caller))
                .then_with(|| a.label.cmp(b.label))
        });
        WorkSpanReport {
            work: self.work,
            spans,
            regions,
            loops,
            dependencies: self.dependencies.into_values().collect(),
            decisions,
            reads: self.reads,
            writes: self.writes,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn charge(units: u64) {
        advance(units);
    }

    // The observer is process-wide, so its tests take turns.
    static SERIAL: Mutex<()> = Mutex::new(());

    #[test]
    fn long_paths_drop_iteratively_and_keep_shared_prefixes_alive() {
        std::thread::Builder::new()
            .stack_size(128 * 1024)
            .spawn(|| {
                let mut clock = Clock::default();
                clock.advance(1);
                clock.flush(0);
                let shared = clock.clone();
                for _ in 0..100_000 {
                    clock.advance(1);
                    clock.flush(1);
                }
                drop(clock);
                assert_eq!(shared.path.as_ref().unwrap().work, 1);
                assert!(shared.path.as_ref().unwrap().previous.is_none());
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn operation_charges_reach_the_observer_without_reborrowing_the_budget() {
        let _serial = SERIAL.lock().unwrap_or_else(PoisonError::into_inner);
        let (((), work), report) = with_work_span_on_this_thread(|| {
            crate::with_operation_budget(10, || {
                assert!(!crate::operation::charge(4));
                assert!(crate::operation::charge(7));
            })
        });
        assert!(work.exhausted);
        assert_eq!(work.used, 4);
        assert_eq!(report.work, work.used);
        assert!(report.spans.iter().all(|row| row.span == work.used));
    }

    #[test]
    fn span_follows_tasks_and_their_tracked_dependencies() {
        let _serial = SERIAL.lock().unwrap_or_else(PoisonError::into_inner);
        // Concurrent tests charge work on their own threads.
        let ((), report) = with_work_span_on_this_thread(|| {
            charge(5);
            let _region = region("stage");
            let tasks = steps("loop");
            for (index, units) in [10u64, 30, 40].into_iter().enumerate() {
                let _task = tasks.task();
                match index {
                    // Waits for the first task's write at time 15.
                    1 => read("value", 1),
                    // Shared bookkeeping waits in every model but the last.
                    2 => read("fact class", 1),
                    _ => {}
                }
                charge(units);
                if index == 0 {
                    write("value", 1);
                    write("fact class", 1);
                }
            }
            // Work inside the loop after its tasks joins them first.
            charge(2);
            drop(tasks);
            // A search runs its candidates in turn unless speculating.
            let candidates = search("search");
            for units in [4u64, 6] {
                let _task = candidates.task();
                charge(units);
            }
            drop(candidates);
            charge(1);
            // Ordered appends leave their tasks for one commit chain. Here the
            // second merge reads the first, so only the ordered model lets
            // the second task's own work overlap the first.
            let merges = parallel("merges");
            for index in 0..2u64 {
                let _task = merges.task();
                charge(5);
                {
                    let _append = append();
                    if index == 1 {
                        read("list", 0);
                    }
                    charge(3);
                    write("list", 0);
                }
                charge(1);
            }
        });
        assert_eq!(report.work, 116);
        // 5, then max(10, 10 + 30, 10 + 40) = 50 or, without the shared wait,
        // max(10, 40, 40) = 40, then 2, the search's 10 or 6, then 1. The
        // merges take 8 + 3 + 1 = 12 when the second waits inline and
        // 5 + 3 + 3 = 11 on the commit chain.
        let spans: Vec<_> = report.spans.iter().map(|row| row.span).collect();
        assert_eq!(spans, [86, 80, 76, 66, 65]);
        for row in &report.spans {
            assert_eq!(row.path.iter().map(|row| row.work).sum::<u64>(), row.span);
        }
        let row = &report.loops[0];
        assert_eq!((row.runs, row.tasks, row.work), (1, 3, 82));
        assert_eq!(row.span, [52, 52, 52, 42, 42]);
        assert_eq!(row.largest_task, 40);
        assert_eq!(report.loops[1].span, [10, 10, 6, 6, 6]);
        // The request region spans everything.
        let request = &report.regions[0];
        assert_eq!(request.region, "request");
        assert_eq!(request.span, [86, 80, 76, 66, 65]);
        let waits: Vec<_> = report
            .dependencies
            .iter()
            .map(|row| {
                (
                    row.kind.as_str(),
                    row.shared,
                    row.waits,
                    row.delay,
                    row.sibling_waits,
                )
            })
            .collect();
        assert_eq!(
            waits,
            [
                ("fact class", true, 1, 10, 1),
                ("list", false, 1, 3, 1),
                ("value", false, 1, 10, 1)
            ]
        );
        // Outside an observation every hook is inert.
        read("value", 1);
        write("value", 1);
        let _ignored = region("outside");
        assert!(lock().is_none());
    }

    // A copy of a table must not see a later rewrite of an entry it
    // inherited: the original's rewrite is the copy's past no longer.
    #[test]
    fn a_copy_does_not_see_its_originals_later_rewrite() {
        let _serial = SERIAL.lock().unwrap_or_else(PoisonError::into_inner);
        let ((), report) = with_work_span_on_this_thread(|| {
            let mut original = Entries::default();
            original.write("edge coedges", 0);
            let copy = original.clone();
            let tasks = parallel("rewrite");
            {
                let _task = tasks.task();
                charge(10);
                original.write("edge coedges", 0);
            }
            {
                let _task = tasks.task();
                copy.read("edge coedges", 0);
                charge(5);
            }
        });
        // max(10, 5) = 10; with the rewrite visible to the copy it was 15.
        // Only the feature-level model runs this loop's tasks in turn.
        let spans: Vec<_> = report.spans.iter().map(|row| row.span).collect();
        assert_eq!(spans, [15, 10, 10, 10, 10]);
        assert!(report.dependencies.is_empty());
    }
}
