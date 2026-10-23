//! Reduced ordered Boolean diagrams for regional material-side predicates.
//! Variables are ordered by operand, then by whether the two germs are shared.
use catallax_numbers::{Decider, Refusal};
use catallax_recipe::BooleanKind;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct State(usize);
impl State {
    pub(super) const OUTSIDE: Self = Self(0);
    pub(super) const INSIDE: Self = Self(1);
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum Germ {
    Shared,
    Positive,
    Negative,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Variable {
    operand: usize,
    germ: Germ,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Node {
    variable: Variable,
    outside: State,
    inside: State,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Op {
    And,
    Or,
}

/// One request owns this arena. Neighbouring cells reuse canonical nodes and
/// apply results; no node or memo table is allocated anew for each cell.
pub(super) struct Boolean {
    nodes: Vec<Node>,
    variables: Vec<[Option<State>; 3]>,
    unique: BTreeMap<Node, State>,
    applied: BTreeMap<(Op, State, State), State>,
    negated: BTreeMap<State, State>,
}
impl Boolean {
    pub(super) fn new(operands: usize) -> Self {
        Self {
            nodes: Vec::new(),
            variables: vec![[None; 3]; operands],
            unique: BTreeMap::new(),
            applied: BTreeMap::new(),
            negated: BTreeMap::new(),
        }
    }
    fn node(&self, state: State) -> Option<Node> {
        state.0.checked_sub(2).map(|i| self.nodes[i])
    }
    fn make(&mut self, node: Node, dec: &mut Decider) -> Result<State, Refusal> {
        dec.check_deadline()?;
        if node.outside == node.inside {
            return Ok(node.outside);
        }
        if let Some(&state) = self.unique.get(&node) {
            return Ok(state);
        }
        let state = State(self.nodes.len() + 2);
        self.nodes.push(node);
        self.unique.insert(node, state);
        Ok(state)
    }
    pub(super) fn variable(
        &mut self,
        operand: usize,
        germ: Germ,
        dec: &mut Decider,
    ) -> Result<State, Refusal> {
        let index = match germ {
            Germ::Shared => 0,
            Germ::Positive => 1,
            Germ::Negative => 2,
        };
        if let Some(state) = self.variables[operand][index] {
            return Ok(state);
        }
        let state = self.make(
            Node {
                variable: Variable { operand, germ },
                outside: State::OUTSIDE,
                inside: State::INSIDE,
            },
            dec,
        )?;
        self.variables[operand][index] = Some(state);
        Ok(state)
    }
    fn negate(&mut self, state: State, dec: &mut Decider) -> Result<State, Refusal> {
        dec.check_deadline()?;
        if state == State::OUTSIDE {
            return Ok(State::INSIDE);
        }
        if state == State::INSIDE {
            return Ok(State::OUTSIDE);
        }
        if let Some(&other) = self.negated.get(&state) {
            return Ok(other);
        }
        let node = self.node(state).expect("nonconstant Boolean node");
        let outside = self.negate(node.outside, dec)?;
        let inside = self.negate(node.inside, dec)?;
        let other = self.make(
            Node {
                variable: node.variable,
                outside,
                inside,
            },
            dec,
        )?;
        self.negated.insert(state, other);
        self.negated.insert(other, state);
        Ok(other)
    }
    fn apply(&mut self, op: Op, a: State, b: State, dec: &mut Decider) -> Result<State, Refusal> {
        dec.check_deadline()?;
        if a == b {
            return Ok(a);
        }
        let (a, b) = if a < b { (a, b) } else { (b, a) };
        match op {
            Op::And if a == State::OUTSIDE => return Ok(State::OUTSIDE),
            Op::And if a == State::INSIDE => return Ok(b),
            Op::Or if a == State::OUTSIDE => return Ok(b),
            Op::Or if a == State::INSIDE => return Ok(State::INSIDE),
            _ => (),
        }
        if let Some(&state) = self.applied.get(&(op, a, b)) {
            return Ok(state);
        }
        let an = self.node(a).expect("constant cases settled");
        let bn = self.node(b).expect("constant cases settled");
        let variable = an.variable.min(bn.variable);
        let (ao, ai) = if an.variable == variable {
            (an.outside, an.inside)
        } else {
            (a, a)
        };
        let (bo, bi) = if bn.variable == variable {
            (bn.outside, bn.inside)
        } else {
            (b, b)
        };
        let outside = self.apply(op, ao, bo, dec)?;
        let inside = self.apply(op, ai, bi, dec)?;
        let state = self.make(
            Node {
                variable,
                outside,
                inside,
            },
            dec,
        )?;
        self.applied.insert((op, a, b), state);
        Ok(state)
    }
    pub(super) fn combine(
        &mut self,
        kind: BooleanKind,
        a: State,
        b: State,
        dec: &mut Decider,
    ) -> Result<State, Refusal> {
        match kind {
            BooleanKind::Union => self.apply(Op::Or, a, b, dec),
            BooleanKind::Intersect => self.apply(Op::And, a, b, dec),
            BooleanKind::Subtract => {
                let b = self.negate(b, dec)?;
                self.apply(Op::And, a, b, dec)
            }
        }
    }
    pub(super) fn unresolved(&self, roots: impl IntoIterator<Item = State>) -> Vec<usize> {
        let mut stack: Vec<_> = roots.into_iter().collect();
        let mut visited = BTreeSet::new();
        let mut operands = BTreeSet::new();
        while let Some(state) = stack.pop() {
            if !visited.insert(state) {
                continue;
            }
            if let Some(node) = self.node(state) {
                operands.insert(node.variable.operand);
                stack.extend([node.outside, node.inside]);
            }
        }
        operands.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn value(b: &Boolean, mut state: State, input: impl Fn(Variable) -> bool) -> bool {
        while let Some(node) = b.node(state) {
            state = if input(node.variable) {
                node.inside
            } else {
                node.outside
            };
        }
        state == State::INSIDE
    }
    #[test]
    fn distributive_recipes_have_one_canonical_state_and_exact_truth_table() {
        let mut b = Boolean::new(3);
        let mut d = Decider::default();
        let [x, y, z] = std::array::from_fn(|i| b.variable(i, Germ::Shared, &mut d).unwrap());
        let xy = b.combine(BooleanKind::Intersect, x, y, &mut d).unwrap();
        let xz = b.combine(BooleanKind::Intersect, z, x, &mut d).unwrap();
        let left = b.combine(BooleanKind::Union, xy, xz, &mut d).unwrap();
        let yz = b.combine(BooleanKind::Union, z, y, &mut d).unwrap();
        let right = b.combine(BooleanKind::Intersect, yz, x, &mut d).unwrap();
        assert_eq!(left, right);
        let cut = b.combine(BooleanKind::Subtract, left, y, &mut d).unwrap();
        for bits in 0..8 {
            let inputs = std::array::from_fn::<_, 3, _>(|i| bits & (1 << i) != 0);
            assert_eq!(
                value(&b, cut, |v| inputs[v.operand]),
                inputs[0] && inputs[2] && !inputs[1]
            );
        }
        for node in &b.nodes {
            for child in [node.outside, node.inside] {
                if let Some(child) = b.node(child) {
                    assert!(node.variable < child.variable);
                }
            }
        }
    }
    #[test]
    fn a_shared_unresolved_cut_cancels_without_becoming_constant_material() {
        let mut b = Boolean::new(5);
        let mut d = Decider::default();
        let cut = b.variable(4, Germ::Shared, &mut d).unwrap();
        let plus = b
            .combine(BooleanKind::Subtract, State::INSIDE, cut, &mut d)
            .unwrap();
        let minus = b
            .combine(BooleanKind::Subtract, State::INSIDE, cut, &mut d)
            .unwrap();
        assert_eq!(plus, minus);
        assert_eq!(b.unresolved([plus, minus]), vec![4]);
        assert!(!value(&b, plus, |_| true));
        assert!(value(&b, plus, |_| false));
        let positive = b.variable(4, Germ::Positive, &mut d).unwrap();
        let negative = b.variable(4, Germ::Negative, &mut d).unwrap();
        let exposed = b
            .combine(BooleanKind::Subtract, positive, negative, &mut d)
            .unwrap();
        assert!(value(&b, exposed, |v| v.germ == Germ::Positive));
        assert!(!value(&b, exposed, |_| true));
    }
}
