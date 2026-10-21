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
        self.edges
            .entry((fact.kind, fact.a))
            .or_default()
            .insert(fact.b);
        if symmetric {
            self.edges
                .entry((fact.kind, fact.b))
                .or_default()
                .insert(fact.a);
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
        let mut seen = BTreeSet::new();
        let mut queue = VecDeque::from([(a, false)]);
        seen.insert((a, false));
        while let Some((key, incident)) = queue.pop_front() {
            for edge_kind in [kind, FactKind::Coincident] {
                if edge_kind != kind && kind != FactKind::Incident {
                    continue;
                }
                if let Some(neighbors) = self.edges.get(&(edge_kind, key)) {
                    for next in neighbors {
                        let used = incident || edge_kind == FactKind::Incident;
                        if *next == b && (kind != FactKind::Incident || used) {
                            return true;
                        }
                        if seen.insert((*next, used)) {
                            queue.push_back((*next, used));
                        }
                    }
                }
            }
        }
        false
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
