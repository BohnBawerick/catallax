//! Exact structural identity of expression graphs, for run-ahead proofs.
//!
//! A decision's work depends on which expression nodes are shared, not only on
//! the expression as a tree: one evaluation charges each distinct node once.
//! Two questions therefore count as the same question only when their graphs
//! are isomorphic node for node. The fingerprint here only finds candidates;
//! `Pairing` proves the match and never relies on the fingerprint.

use crate::{
    Real,
    algebraic::{AlgebraicRepr, AlgebraicRoot},
    real::{Expr, Node},
};
use std::{
    collections::HashMap,
    hash::{Hash, Hasher},
    sync::Arc,
};

fn hasher(tag: u8) -> std::collections::hash_map::DefaultHasher {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    tag.hash(&mut hasher);
    hasher
}

fn children(node: &Node) -> (Option<&Real>, Option<&Real>) {
    match &node.expr {
        Expr::Rational(_) | Expr::Algebraic(_) | Expr::Coordinate(..) => (None, None),
        Expr::Neg(a) | Expr::Sqrt(a) | Expr::Square(a) | Expr::Pow(a, _) => (Some(a), None),
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => (Some(a), Some(b)),
    }
}

fn tag(expr: &Expr) -> u8 {
    match expr {
        Expr::Rational(_) => 0,
        Expr::Algebraic(_) => 1,
        Expr::Coordinate(..) => 2,
        Expr::Add(..) => 3,
        Expr::Sub(..) => 4,
        Expr::Mul(..) => 5,
        Expr::Div(..) => 6,
        Expr::Neg(_) => 7,
        Expr::Sqrt(_) => 8,
        Expr::Square(_) => 9,
        Expr::Pow(..) => 10,
    }
}

/// A cheap sorting key for finding a recorded question: the kinds and the
/// hardware bounds of the top few nodes. It is the same for two graphs that
/// match node for node, and takes the same time whatever their size. It
/// proves nothing; [`Pairing`] does.
pub(crate) fn fingerprint(x: &Real) -> u64 {
    fn mix(node: &Node, depth: u32, state: &mut impl Hasher) {
        tag(&node.expr).hash(state);
        node.approx.lo().to_bits().hash(state);
        node.approx.hi().to_bits().hash(state);
        node.fast_defined.hash(state);
        if depth != 0 {
            let (first, second) = children(node);
            for child in [first, second].into_iter().flatten() {
                mix(&child.node, depth - 1, state);
            }
        }
    }
    let mut state = hasher(0);
    mix(&x.node, 2, &mut state);
    state.finish()
}

/// A proved one-to-one match between the nodes of a recorded question and the
/// nodes of the question being asked.
///
/// Nodes and roots made before the record's fork (`epoch`) are shared by both
/// sides, so they may only match themselves; whatever they refer to is older
/// still and needs no visit. Everything newer is visited and entered here.
#[derive(Default)]
pub(crate) struct Pairing {
    forward: HashMap<*const Node, *const Node>,
    backward: HashMap<*const Node, *const Node>,
    roots: HashMap<*const AlgebraicRoot, Arc<AlgebraicRoot>>,
    roots_backward: HashMap<*const AlgebraicRoot, *const AlgebraicRoot>,
}

impl Pairing {
    /// The asking side's root for a recorded root, when the question holds it.
    pub(crate) fn root(&self, recorded: &Arc<AlgebraicRoot>) -> Option<Arc<AlgebraicRoot>> {
        self.roots.get(&Arc::as_ptr(recorded)).cloned()
    }

    fn pair_roots(
        &mut self,
        recorded: &Arc<AlgebraicRoot>,
        mine: &Arc<AlgebraicRoot>,
        epoch: u64,
    ) -> bool {
        let (a, b) = (Arc::as_ptr(recorded), Arc::as_ptr(mine));
        if a != b && (recorded.epoch < epoch || mine.epoch < epoch) {
            return false;
        }
        match self.roots.get(&a) {
            Some(known) => return Arc::as_ptr(known) == b,
            None => {
                if self.roots_backward.contains_key(&b) {
                    return false;
                }
            }
        }
        if a != b
            && (recorded.polynomial != mine.polynomial
                || recorded.lo != mine.lo
                || recorded.hi != mine.hi)
        {
            return false;
        }
        self.roots.insert(a, mine.clone());
        self.roots_backward.insert(b, a);
        true
    }

    /// Whether `mine` is the recorded question `recorded`, node for node.
    /// Several values of one question share this pairing.
    pub(crate) fn same(&mut self, recorded: &Real, mine: &Real, epoch: u64) -> bool {
        let mut pending: Vec<(&Node, &Node)> = vec![(&recorded.node, &mine.node)];
        while let Some((a, b)) = pending.pop() {
            let (pa, pb) = (std::ptr::from_ref(a), std::ptr::from_ref(b));
            if pa == pb {
                if a.epoch < epoch {
                    continue;
                }
            } else if a.epoch < epoch || b.epoch < epoch {
                return false;
            }
            match self.forward.get(&pa) {
                Some(&known) => {
                    if known != pb {
                        return false;
                    }
                    continue;
                }
                None => {
                    if self.backward.contains_key(&pb) {
                        return false;
                    }
                }
            }
            if a.fast_defined != b.fast_defined
                || a.approx.lo().to_bits() != b.approx.lo().to_bits()
                || a.approx.hi().to_bits() != b.approx.hi().to_bits()
            {
                return false;
            }
            let same_node = match (&a.expr, &b.expr) {
                (Expr::Rational(x), Expr::Rational(y)) => x == y,
                (Expr::Algebraic(x), Expr::Algebraic(y)) => match (&x.repr, &y.repr) {
                    (AlgebraicRepr::Rational(x), AlgebraicRepr::Rational(y)) => x == y,
                    (AlgebraicRepr::Root(x), AlgebraicRepr::Root(y)) => {
                        self.pair_roots(x, y, epoch)
                    }
                    _ => false,
                },
                (Expr::Coordinate(x, i), Expr::Coordinate(y, j)) => {
                    i == j
                        && (Arc::ptr_eq(x, y)
                            || (x.system().equations() == y.system().equations()
                                && x.root_box() == y.root_box()))
                }
                (Expr::Pow(_, n), Expr::Pow(_, m)) => n == m,
                (x, y) => tag(x) == tag(y),
            };
            if !same_node {
                return false;
            }
            self.forward.insert(pa, pb);
            self.backward.insert(pb, pa);
            let (a_first, a_second) = children(a);
            let (b_first, b_second) = children(b);
            for (x, y) in [(a_first, b_first), (a_second, b_second)] {
                match (x, y) {
                    (Some(x), Some(y)) => pending.push((&x.node, &y.node)),
                    (None, None) => {}
                    _ => return false,
                }
            }
        }
        true
    }
}
