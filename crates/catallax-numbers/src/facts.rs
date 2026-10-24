//! Recorded construction facts.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

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
pub struct FactBook<K> {
    facts: Vec<Fact<K>>,
    keys: BTreeSet<(FactKind, K, K, Option<K>)>,
    edges: BTreeMap<(FactKind, K), BTreeSet<K>>,
    classes: BTreeMap<(FactKind, K), usize>,
    components: Vec<Component>,
}

// Union by size bounds the read-only root walk. Incidence edges name roots,
// not all historical copies of an equivalent geometric object.
#[derive(Clone, Debug)]
struct Component {
    parent: usize,
    size: usize,
    outgoing: BTreeSet<usize>,
    incoming: BTreeSet<usize>,
}

impl<K: Copy + Ord> Default for FactBook<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: Copy + Ord> FactBook<K> {
    /// An empty book.
    pub fn new() -> Self {
        Self {
            facts: Vec::new(),
            keys: BTreeSet::new(),
            edges: BTreeMap::new(),
            classes: BTreeMap::new(),
            components: Vec::new(),
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
        if !self.keys.insert((fact.kind, a, b, fact.at)) {
            return false;
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
            }
            kind => {
                let a = self.component(kind, fact.a);
                let b = self.component(kind, fact.b);
                self.unite(a, b);
            }
        }
        self.facts.push(fact);
        true
    }

    /// Whether the relation holds, with the closure rules in the piece description.
    pub fn holds(&self, kind: FactKind, a: K, b: K) -> bool {
        if kind == FactKind::Tangent {
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
        if self.components[a].outgoing.contains(&b) {
            return true;
        }
        if self.components[a].outgoing.is_empty() {
            return false;
        }
        let mut seen = BTreeSet::from([a]);
        let mut queue = VecDeque::from([a]);
        while let Some(key) = queue.pop_front() {
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
        while self.components[node].parent != node {
            node = self.components[node].parent;
        }
        node
    }

    fn component(&mut self, kind: FactKind, key: K) -> usize {
        if let Some(&node) = self.classes.get(&(kind, key)) {
            return self.root(node);
        }
        let node = self.components.len();
        self.components.push(Component {
            parent: node,
            size: 1,
            outgoing: BTreeSet::new(),
            incoming: BTreeSet::new(),
        });
        self.classes.insert((kind, key), node);
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
        // Rewire both directions so later queries never scan dead members.
        // An incidence within a merged class remains a real directed self edge.
        let outgoing = std::mem::take(&mut self.components[b].outgoing);
        let incoming = std::mem::take(&mut self.components[b].incoming);
        for next in outgoing {
            let next = self.root(next);
            self.components[next].incoming.remove(&b);
            self.components[next].incoming.insert(a);
            self.components[a].outgoing.insert(next);
        }
        for previous in incoming {
            let previous = self.root(previous);
            self.components[previous].outgoing.remove(&b);
            self.components[previous].outgoing.insert(a);
            self.components[a].incoming.insert(previous);
        }
        self.components[a].outgoing.remove(&b);
        self.components[a].incoming.remove(&b);
    }

    /// Facts that name `key` as `a`, `b` or `at`, in recording order.
    pub fn facts_about(&self, key: K) -> Vec<&Fact<K>> {
        self.facts
            .iter()
            .filter(|fact| fact.a == key || fact.b == key || fact.at == Some(key))
            .collect()
    }

    /// All facts in recording order.
    pub fn iter(&self) -> std::slice::Iter<'_, Fact<K>> {
        self.facts.iter()
    }

    /// Number of facts.
    pub fn len(&self) -> usize {
        self.facts.len()
    }

    /// Whether the book is empty.
    pub fn is_empty(&self) -> bool {
        self.facts.is_empty()
    }
}
