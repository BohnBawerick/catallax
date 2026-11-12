//! Recorded construction facts.

use im::{OrdMap, OrdSet, Vector};
use std::collections::{BTreeSet, VecDeque};

use crate::{Watch, ZeroProof};

/// Kinds of recorded relation.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum FactKind {
    /// Same point set.
    Coincident,
    /// Sharing an axis line.
    Coaxial,
    /// Sharing a centre.
    Concentric,
    /// Parallel normals, directions or axes.
    Parallel,
    /// Tangent, optionally along or at the `at` key.
    Tangent,
    /// The first subject lies on the second.
    Incident,
}

/// Where a fact came from.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum FactOrigin {
    /// Guaranteed by the named constructing routine.
    Construction(&'static str),
    /// Proven by an exact method.
    Proof(ZeroProof),
}

/// One recorded relation between two keyed objects.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Fact<K> {
    /// The relation.
    pub kind: FactKind,
    /// First subject.
    pub a: K,
    /// Second subject.
    pub b: K,
    /// Where the relation holds, when that is a separate object.
    pub at: Option<K>,
    /// Where the fact came from.
    pub origin: FactOrigin,
}

/// One thing a book was asked that a later record can change the answer to.
///
/// A book only grows, but its answers are not all about what it holds for
/// good: whether a fact was recorded, which class a subject is in, what a
/// class is tied to. [`FactBook::watch`] keeps a list of these while work
/// runs on a copy, and [`FactBook::unchanged_where_read`] says whether the
/// book the copy came from still answers each the same.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum FactRead<K> {
    /// Whether this fact is in the book.
    Recorded(FactKind, K, K, Option<K>),
    /// What a subject is tangent to.
    Tangents(K),
    /// Which class of a relation a subject is in, if any.
    Class(FactKind, K),
    /// The class a class was merged into.
    Parent(usize),
    /// How many subjects a class holds.
    Size(usize),
    /// What a class lies on.
    Outgoing(usize),
    /// What lies on a class.
    Incoming(usize),
    /// The facts that name a subject.
    About(K),
    /// The whole list of facts.
    All,
}

/// Append-only, deterministic store of facts.
#[derive(Clone, Debug)]
pub struct FactBook<K: Clone + Ord> {
    facts: Vector<Fact<K>>,
    keys: OrdSet<(FactKind, K, K, Option<K>)>,
    edges: OrdMap<(FactKind, K), OrdSet<K>>,
    classes: OrdMap<(FactKind, K), usize>,
    components: Vector<Component>,
    // What this book and its copies have been asked since `watch`.
    watch: Watch<FactRead<K>>,
    // Work-span observation only: a dense number for each recorded subject,
    // so tangent and per-subject reads have their own dependency keys.
    #[cfg(feature = "work-span")]
    subjects: OrdMap<K, u64>,
    #[cfg(feature = "work-span")]
    span_nodes: crate::work_span::Entries,
    #[cfg(feature = "work-span")]
    span_subjects: crate::work_span::Entries,
}

// Union by size bounds the read-only root walk. Incidence edges name roots,
// not all historical copies of an equivalent geometric object.
#[derive(Clone, Debug)]
struct Component {
    parent: usize,
    size: usize,
    outgoing: OrdSet<usize>,
    incoming: OrdSet<usize>,
}

impl<K: Copy + Ord> Default for FactBook<K> {
    fn default() -> Self {
        Self::new()
    }
}

// Dependency keys for the work-span observer: a component's parent, size,
// and outgoing and incoming incidence sets, four entries per component; and
// per subject its facts and its tangents. Each copy of a book keeps its own
// writers, so copies never see each other's later changes.
#[cfg(feature = "work-span")]
impl<K: Copy + Ord> FactBook<K> {
    fn span_class(&self, node: usize) {
        self.span_nodes.read("fact class", 4 * node as u64);
    }
    fn span_set_class(&mut self, node: usize) {
        self.span_nodes.write("fact class", 4 * node as u64);
    }
    fn span_size(&self, node: usize) {
        self.span_nodes.read("fact class size", 4 * node as u64 + 1);
    }
    fn span_set_size(&mut self, node: usize) {
        self.span_nodes
            .write("fact class size", 4 * node as u64 + 1);
    }
    fn span_links(&self, node: usize, incoming: bool) {
        self.span_nodes
            .read("fact links", 4 * node as u64 + 2 + u64::from(incoming));
    }
    fn span_set_links(&mut self, node: usize, incoming: bool) {
        self.span_nodes
            .write("fact links", 4 * node as u64 + 2 + u64::from(incoming));
    }
}

impl<K: Copy + Ord> FactBook<K> {
    /// An empty book.
    pub fn new() -> Self {
        Self {
            facts: Vector::new(),
            keys: OrdSet::new(),
            edges: OrdMap::new(),
            classes: OrdMap::new(),
            components: Vector::new(),
            watch: Watch::default(),
            #[cfg(feature = "work-span")]
            subjects: OrdMap::new(),
            #[cfg(feature = "work-span")]
            span_nodes: Default::default(),
            #[cfg(feature = "work-span")]
            span_subjects: Default::default(),
        }
    }

    /// Records a fact; `false` if an equal fact is already present.
    pub fn record(&mut self, fact: Fact<K>) -> bool {
        if fact.a == fact.b {
            return false;
        }
        let symmetric = fact.kind != FactKind::Incident;
        let (a, b) = if symmetric && fact.b < fact.a {
            (fact.b, fact.a)
        } else {
            (fact.a, fact.b)
        };
        self.note(|| FactRead::Recorded(fact.kind, a, b, fact.at));
        if self.keys.insert((fact.kind, a, b, fact.at)).is_some() {
            return false;
        }
        #[cfg(feature = "work-span")]
        for subject in [Some(fact.a), Some(fact.b), fact.at].into_iter().flatten() {
            let next = self.subjects.len() as u64;
            let id = *self.subjects.entry(subject).or_insert(next);
            self.span_subjects.write("fact list", 2 * id);
            if fact.kind == FactKind::Tangent {
                self.span_subjects.write("fact tangents", 2 * id + 1);
            }
        }
        if fact.kind == FactKind::Tangent {
            self.edges
                .entry((fact.kind, fact.a))
                .or_default()
                .insert(fact.b);
            self.edges
                .entry((fact.kind, fact.b))
                .or_default()
                .insert(fact.a);
        }
        match fact.kind {
            FactKind::Tangent => {}
            FactKind::Incident => {
                let a = self.component(FactKind::Coincident, fact.a);
                let b = self.component(FactKind::Coincident, fact.b);
                self.components[a].outgoing.insert(b);
                self.components[b].incoming.insert(a);
                #[cfg(feature = "work-span")]
                {
                    self.span_set_links(a, false);
                    self.span_set_links(b, true);
                }
            }
            kind => {
                let a = self.component(kind, fact.a);
                let b = self.component(kind, fact.b);
                self.unite(a, b);
            }
        }
        self.facts.push_back(fact);
        #[cfg(feature = "work-span")]
        self.span_subjects.write("fact list", u64::MAX);
        true
    }

    /// Whether the relation holds, with the closure rules in the piece description.
    pub fn holds(&self, kind: FactKind, a: K, b: K) -> bool {
        if kind == FactKind::Tangent {
            #[cfg(feature = "work-span")]
            if let Some(&id) = self.subjects.get(&a) {
                self.span_subjects.read("fact tangents", 2 * id + 1);
            }
            self.note(|| FactRead::Tangents(a));
            return self
                .edges
                .get(&(kind, a))
                .is_some_and(|neighbors| neighbors.contains(&b));
        }
        if kind != FactKind::Incident && a == b {
            return true;
        }
        let relation = if kind == FactKind::Incident {
            FactKind::Coincident
        } else {
            kind
        };
        self.note(|| FactRead::Class(relation, a));
        self.note(|| FactRead::Class(relation, b));
        let (Some(&a), Some(&b)) = (
            self.classes.get(&(relation, a)),
            self.classes.get(&(relation, b)),
        ) else {
            return false;
        };
        let (a, b) = (self.root(a), self.root(b));
        if kind != FactKind::Incident {
            return a == b;
        }
        #[cfg(feature = "work-span")]
        self.span_links(a, false);
        self.note(|| FactRead::Outgoing(a));
        if self.components[a].outgoing.contains(&b) {
            return true;
        }
        if self.components[a].outgoing.is_empty() {
            return false;
        }
        let mut seen = BTreeSet::from([a]);
        let mut queue = VecDeque::from([a]);
        while let Some(key) = queue.pop_front() {
            #[cfg(feature = "work-span")]
            self.span_links(key, false);
            self.note(|| FactRead::Outgoing(key));
            for &next in &self.components[key].outgoing {
                if next == b {
                    return true;
                }
                if seen.insert(next) {
                    queue.push_back(next);
                }
            }
        }
        false
    }

    fn root(&self, mut node: usize) -> usize {
        #[cfg(feature = "work-span")]
        self.span_class(node);
        self.note(|| FactRead::Parent(node));
        while self.components[node].parent != node {
            node = self.components[node].parent;
            #[cfg(feature = "work-span")]
            self.span_class(node);
            self.note(|| FactRead::Parent(node));
        }
        node
    }

    fn component(&mut self, kind: FactKind, key: K) -> usize {
        self.note(|| FactRead::Class(kind, key));
        if let Some(&node) = self.classes.get(&(kind, key)) {
            return self.root(node);
        }
        let node = self.components.len();
        self.components.push_back(Component {
            parent: node,
            size: 1,
            outgoing: OrdSet::new(),
            incoming: OrdSet::new(),
        });
        self.classes.insert((kind, key), node);
        #[cfg(feature = "work-span")]
        self.span_set_class(node);
        node
    }

    fn unite(&mut self, mut a: usize, mut b: usize) {
        if a == b {
            return;
        }
        self.note(|| FactRead::Size(a));
        self.note(|| FactRead::Size(b));
        if self.components[a].size < self.components[b].size {
            std::mem::swap(&mut a, &mut b);
        }
        // Which classes the smaller one was tied to decides what is rewired.
        self.note(|| FactRead::Outgoing(b));
        self.note(|| FactRead::Incoming(b));
        self.components[b].parent = a;
        self.components[a].size += self.components[b].size;
        #[cfg(feature = "work-span")]
        {
            self.span_size(a);
            self.span_size(b);
            self.span_set_size(a);
            self.span_set_class(b);
            self.span_links(b, false);
            self.span_links(b, true);
            self.span_set_links(a, false);
            self.span_set_links(a, true);
            self.span_set_links(b, false);
            self.span_set_links(b, true);
        }
        // Rewire both directions so later queries never scan dead members.
        // An incidence within a merged class remains a real directed self edge.
        let outgoing = std::mem::take(&mut self.components[b].outgoing);
        let incoming = std::mem::take(&mut self.components[b].incoming);
        for next in outgoing {
            let next = self.root(next);
            #[cfg(feature = "work-span")]
            self.span_set_links(next, true);
            self.components[next].incoming.remove(&b);
            self.components[next].incoming.insert(a);
            self.components[a].outgoing.insert(next);
        }
        for previous in incoming {
            let previous = self.root(previous);
            #[cfg(feature = "work-span")]
            self.span_set_links(previous, false);
            self.components[previous].outgoing.remove(&b);
            self.components[previous].outgoing.insert(a);
            self.components[a].incoming.insert(previous);
        }
        self.components[a].outgoing.remove(&b);
        self.components[a].incoming.remove(&b);
    }

    /// Facts that name `key` as `a`, `b` or `at`, in recording order.
    pub fn facts_about(&self, key: K) -> Vec<&Fact<K>> {
        #[cfg(feature = "work-span")]
        if let Some(&id) = self.subjects.get(&key) {
            self.span_subjects.read("fact list", 2 * id);
        }
        self.note(|| FactRead::About(key));
        self.facts
            .iter()
            .filter(|fact| fact.a == key || fact.b == key || fact.at == Some(key))
            .collect()
    }

    /// All facts in recording order.
    pub fn iter(&self) -> im::vector::Iter<'_, Fact<K>> {
        #[cfg(feature = "work-span")]
        self.span_subjects.read("fact list", u64::MAX);
        self.note(|| FactRead::All);
        self.facts.iter()
    }

    /// Number of facts.
    pub fn len(&self) -> usize {
        #[cfg(feature = "work-span")]
        self.span_subjects.read("fact list", u64::MAX);
        self.note(|| FactRead::All);
        self.facts.len()
    }

    /// Whether the book is empty.
    pub fn is_empty(&self) -> bool {
        #[cfg(feature = "work-span")]
        self.span_subjects.read("fact list", u64::MAX);
        self.note(|| FactRead::All);
        self.facts.is_empty()
    }

    fn note(&self, read: impl FnOnce() -> FactRead<K>) {
        self.watch.note(read);
    }

    /// Starts a list of what this book is asked from here on. Take a copy of
    /// the book first: the list is for checking the original.
    pub fn watch(&mut self) {
        self.watch.begin();
    }

    /// Ends the list and returns it, each read once. Panics when none was
    /// begun.
    pub fn watched(&mut self) -> Vec<FactRead<K>> {
        self.watch.end()
    }

    /// Whether this book answers every read in `reads` as `base` does. This
    /// book must be `base` or have grown from it, and `reads` must come from
    /// work on a copy of `base`; `old` says whether a subject was in `base`'s
    /// store. A read about a subject or class the work made itself is the
    /// work's own business and always passes.
    ///
    /// When this holds, the work would have found the same answers here, so
    /// it would have done the same and recorded the same.
    pub fn unchanged_where_read(
        &self,
        base: &Self,
        reads: &[FactRead<K>],
        old: impl Fn(K) -> bool,
    ) -> bool {
        let classes = base.components.len();
        let field = |node: usize, same: fn(&Component, &Component) -> bool| {
            node >= classes || same(&base.components[node], &self.components[node])
        };
        reads.iter().all(|read| match *read {
            FactRead::Recorded(kind, a, b, at) => {
                !(old(a) && old(b) && at.is_none_or(&old)) || {
                    let key = (kind, a, b, at);
                    base.keys.contains(&key) == self.keys.contains(&key)
                }
            }
            FactRead::Tangents(key) => {
                !old(key) || {
                    let key = (FactKind::Tangent, key);
                    base.edges.get(&key) == self.edges.get(&key)
                }
            }
            FactRead::Class(kind, key) => {
                !old(key) || base.classes.get(&(kind, key)) == self.classes.get(&(kind, key))
            }
            FactRead::Parent(node) => field(node, |a, b| a.parent == b.parent),
            FactRead::Size(node) => field(node, |a, b| a.size == b.size),
            FactRead::Outgoing(node) => field(node, |a, b| a.outgoing == b.outgoing),
            FactRead::Incoming(node) => field(node, |a, b| a.incoming == b.incoming),
            FactRead::About(key) => {
                !old(key)
                    || self
                        .facts
                        .iter()
                        .skip(base.facts.len())
                        .all(|fact| fact.a != key && fact.b != key && fact.at != Some(key))
            }
            FactRead::All => self.facts.len() == base.facts.len(),
        })
    }
}

#[cfg(test)]
mod watch_tests {
    use super::*;

    fn fact(kind: FactKind, a: u32, b: u32) -> Fact<u32> {
        Fact {
            kind,
            a,
            b,
            at: None,
            origin: FactOrigin::Construction("test"),
        }
    }

    // Subjects 0 to 9 are in the base book's store; work makes 10 and up.
    fn old(key: u32) -> bool {
        key < 10
    }

    fn base() -> FactBook<u32> {
        let mut book = FactBook::new();
        book.record(fact(FactKind::Incident, 1, 2));
        book.record(fact(FactKind::Incident, 2, 3));
        book.record(fact(FactKind::Coincident, 4, 5));
        book.record(fact(FactKind::Tangent, 6, 7));
        book
    }

    /// Runs `work` on a watched copy of `base` and returns what it read.
    fn reads(base: &FactBook<u32>, work: impl FnOnce(&mut FactBook<u32>)) -> Vec<FactRead<u32>> {
        let mut copy = base.clone();
        copy.watch();
        work(&mut copy);
        copy.watched()
    }

    #[test]
    fn a_book_that_has_not_changed_answers_the_same() {
        let base = base();
        let reads = reads(&base, |book| {
            assert!(book.holds(FactKind::Incident, 1, 3));
            assert!(!book.holds(FactKind::Incident, 3, 1));
            assert!(book.holds(FactKind::Coincident, 5, 4));
            assert!(book.holds(FactKind::Tangent, 7, 6));
            assert_eq!(book.facts_about(2).len(), 2);
            assert_eq!(book.len(), 4);
            book.record(fact(FactKind::Incident, 10, 1));
            assert!(book.holds(FactKind::Incident, 10, 3));
        });
        assert!(!reads.is_empty());
        assert!(base.unchanged_where_read(&base, &reads, old));
    }

    #[test]
    fn records_about_other_subjects_change_no_answer() {
        let base = base();
        let reads = reads(&base, |book| {
            assert!(book.holds(FactKind::Incident, 1, 3));
            assert!(book.holds(FactKind::Coincident, 4, 5));
            book.record(fact(FactKind::Incident, 10, 1));
            assert!(book.holds(FactKind::Incident, 10, 3));
        });
        // Other work put its own point on 3 and tied two subjects this work
        // never asked about.
        let mut grown = base.clone();
        grown.record(fact(FactKind::Incident, 11, 3));
        grown.record(fact(FactKind::Coincident, 8, 9));
        grown.record(fact(FactKind::Tangent, 8, 9));
        assert!(grown.unchanged_where_read(&base, &reads, old));
    }

    #[test]
    fn a_record_the_work_could_have_used_is_seen() {
        let base = base();
        // The work found no path from 1 to 8.
        let asked = reads(&base, |book| {
            assert!(!book.holds(FactKind::Incident, 1, 8));
        });
        let mut grown = base.clone();
        grown.record(fact(FactKind::Incident, 3, 8));
        assert!(grown.holds(FactKind::Incident, 1, 8));
        assert!(!grown.unchanged_where_read(&base, &asked, old));

        // The work recorded a fact that was new to its copy.
        let recorded = reads(&base, |book| {
            assert!(book.record(fact(FactKind::Coincident, 8, 9)));
        });
        let mut grown = base.clone();
        assert!(grown.record(fact(FactKind::Coincident, 8, 9)));
        assert!(!grown.unchanged_where_read(&base, &recorded, old));

        // The work found a subject in no class.
        let asked = reads(&base, |book| {
            assert!(!book.holds(FactKind::Coincident, 8, 4));
        });
        let mut grown = base.clone();
        grown.record(fact(FactKind::Coincident, 8, 5));
        assert!(!grown.unchanged_where_read(&base, &asked, old));

        // The work asked what touches 6.
        let asked = reads(&base, |book| {
            assert!(!book.holds(FactKind::Tangent, 6, 8));
        });
        let mut grown = base.clone();
        grown.record(fact(FactKind::Tangent, 8, 6));
        assert!(!grown.unchanged_where_read(&base, &asked, old));

        // The work listed the facts about 2, and counted them all.
        let listed = reads(&base, |book| {
            let _ = book.facts_about(2);
        });
        let counted = reads(&base, |book| {
            let _ = book.len();
        });
        let mut grown = base.clone();
        grown.record(fact(FactKind::Incident, 11, 2));
        assert!(!grown.unchanged_where_read(&base, &listed, old));
        assert!(!grown.unchanged_where_read(&base, &counted, old));
        let mut elsewhere = base.clone();
        elsewhere.record(fact(FactKind::Incident, 11, 3));
        assert!(elsewhere.unchanged_where_read(&base, &listed, old));
        assert!(!elsewhere.unchanged_where_read(&base, &counted, old));
    }

    #[test]
    fn work_on_a_copy_of_a_watched_copy_is_listed_too() {
        let base = base();
        let reads = reads(&base, |book| {
            let scratch = book.clone();
            assert!(!scratch.holds(FactKind::Incident, 1, 8));
        });
        let mut grown = base.clone();
        grown.record(fact(FactKind::Incident, 3, 8));
        assert!(!grown.unchanged_where_read(&base, &reads, old));
    }
}
