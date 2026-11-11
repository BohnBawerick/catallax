//! Recorded construction facts.

use im::{OrdMap, OrdSet, Vector};
use std::collections::{BTreeSet, VecDeque};

use crate::ZeroProof;

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

/// Append-only, deterministic store of facts.
#[derive(Clone, Debug)]
pub struct FactBook<K: Clone + Ord> {
    facts: Vector<Fact<K>>,
    keys: OrdSet<(FactKind, K, K, Option<K>)>,
    edges: OrdMap<(FactKind, K), OrdSet<K>>,
    classes: OrdMap<(FactKind, K), usize>,
    components: Vector<Component>,
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
        while self.components[node].parent != node {
            node = self.components[node].parent;
            #[cfg(feature = "work-span")]
            self.span_class(node);
        }
        node
    }

    fn component(&mut self, kind: FactKind, key: K) -> usize {
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
        if self.components[a].size < self.components[b].size {
            std::mem::swap(&mut a, &mut b);
        }
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
        self.facts
            .iter()
            .filter(|fact| fact.a == key || fact.b == key || fact.at == Some(key))
            .collect()
    }

    /// All facts in recording order.
    pub fn iter(&self) -> im::vector::Iter<'_, Fact<K>> {
        #[cfg(feature = "work-span")]
        self.span_subjects.read("fact list", u64::MAX);
        self.facts.iter()
    }

    /// Number of facts.
    pub fn len(&self) -> usize {
        #[cfg(feature = "work-span")]
        self.span_subjects.read("fact list", u64::MAX);
        self.facts.len()
    }

    /// Whether the book is empty.
    pub fn is_empty(&self) -> bool {
        #[cfg(feature = "work-span")]
        self.span_subjects.read("fact list", u64::MAX);
        self.facts.is_empty()
    }
}
