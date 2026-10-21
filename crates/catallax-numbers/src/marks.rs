//! Marks: the record of every not-exact answer.

use crate::{NotExact, Rational};

/// What a mark's bound measures.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum BoundKind {
    /// A length in model units.
    Distance,
    /// An angle in radians.
    Angle,
    /// The decided quantity itself.
    Value,
}

/// The record of one not-exact answer.
#[derive(Clone, PartialEq, Debug)]
pub struct Mark {
    /// Stable reason code.
    pub code: &'static str,
    /// Free text for people; not part of the contract.
    pub detail: String,
    /// Proven bound; never negative.
    pub bound: Rational,
    /// What the bound measures.
    pub kind: BoundKind,
}

impl Mark {
    /// A mark with these parts.
    pub fn new(
        code: &'static str,
        kind: BoundKind,
        bound: Rational,
        detail: impl Into<String>,
    ) -> Self {
        assert!(bound.sign() != crate::Sign::Negative, "negative mark bound");
        Self {
            code,
            kind,
            bound,
            detail: detail.into(),
        }
    }

    /// A mark that keeps the code of `not_exact`, with the bound converted by the caller.
    pub fn from_not_exact(
        not_exact: &NotExact,
        kind: BoundKind,
        bound: Rational,
        detail: impl Into<String>,
    ) -> Self {
        Self::new(not_exact.code, kind, bound, detail)
    }
}

/// Ordered list of marks.
#[derive(Clone, PartialEq, Debug, Default)]
pub struct MarkLog {
    marks: Vec<Mark>,
}

impl MarkLog {
    /// An empty log.
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends a mark.
    pub fn push(&mut self, mark: Mark) {
        self.marks.push(mark);
    }

    /// Appends every mark of `other`, in order.
    pub fn extend(&mut self, other: &MarkLog) {
        self.marks.extend(other.marks.iter().cloned());
    }

    /// Marks in order.
    pub fn iter(&self) -> std::slice::Iter<'_, Mark> {
        self.marks.iter()
    }

    /// Number of marks.
    pub fn len(&self) -> usize {
        self.marks.len()
    }

    /// Whether the log is empty.
    pub fn is_empty(&self) -> bool {
        self.marks.is_empty()
    }

    /// Largest bound of this kind, if any.
    pub fn worst(&self, kind: BoundKind) -> Option<&Rational> {
        self.marks
            .iter()
            .filter(|mark| mark.kind == kind)
            .map(|mark| &mark.bound)
            .max()
    }
}
