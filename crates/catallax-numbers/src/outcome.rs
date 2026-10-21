//! The decision outcome vocabulary and operation-level answers.

use std::fmt;

use crate::{MarkLog, Rational};

/// The result of one decision. See `docs/ARCHITECTURE.md`.
#[must_use]
#[derive(Clone, Debug, PartialEq)]
pub enum Outcome<T> {
    /// The strict answer is certified; the quantity is proven non-zero.
    Proven(T),
    /// The quantity is proven exactly zero.
    ProvenCoincident(Coincidence),
    /// The budget ran out; the quantity is proven to lie within the bound.
    NotExact(NotExact),
    /// The decision cannot be answered.
    Refused(Refusal),
}

impl<T> Outcome<T> {
    /// Maps the proven answer, keeping the other variants.
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> Outcome<U> {
        match self {
            Self::Proven(value) => Outcome::Proven(f(value)),
            Self::ProvenCoincident(c) => Outcome::ProvenCoincident(c),
            Self::NotExact(n) => Outcome::NotExact(n),
            Self::Refused(r) => Outcome::Refused(r),
        }
    }

    /// The proven answer, if any.
    pub fn proven(&self) -> Option<&T> {
        if let Self::Proven(value) = self {
            Some(value)
        } else {
            None
        }
    }

    /// Whether the outcome is `Proven`.
    pub fn is_proven(&self) -> bool {
        matches!(self, Self::Proven(_))
    }

    /// Whether the outcome is `ProvenCoincident`.
    pub fn is_coincident(&self) -> bool {
        matches!(self, Self::ProvenCoincident(_))
    }

    /// Whether the outcome is `NotExact`.
    pub fn is_not_exact(&self) -> bool {
        matches!(self, Self::NotExact(_))
    }

    /// Whether the outcome is `Refused`.
    pub fn is_refused(&self) -> bool {
        matches!(self, Self::Refused(_))
    }
}

/// A sign that is proven non-zero.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum StrictSign {
    /// Below zero.
    Negative,
    /// Above zero.
    Positive,
}

/// An order between two values that are proven different.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum StrictOrder {
    /// The first is smaller.
    Less,
    /// The first is larger.
    Greater,
}

/// How a zero was proven.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum ZeroProof {
    /// Exact rational evaluation gave zero.
    ExactArithmetic,
    /// An exact algebraic method showed it.
    Algebraic,
    /// An enclosure narrower than a proven root-separation bound.
    SeparationBound,
    /// A recorded construction fact (produced by piece 2, not by this piece).
    ConstructionFact,
}

/// The proof attached to a `ProvenCoincident` outcome.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Coincidence {
    /// How the zero was proven.
    pub proof: ZeroProof,
}

/// A decision whose budget ran out. The true quantity lies in `[-bound, bound]`.
#[derive(Clone, PartialEq, Debug)]
pub struct NotExact {
    /// Proven bound on the magnitude of the quantity; never negative.
    pub bound: Rational,
    /// Stable reason code, such as `numbers.budget-exhausted`.
    pub code: &'static str,
    /// Free text for people; not part of the contract.
    pub detail: String,
}

/// A named refusal.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Refusal {
    /// Stable code of the form `piece.reason`.
    pub code: &'static str,
    /// Free text for people; not part of the contract.
    pub detail: String,
}

impl Refusal {
    /// A refusal with this code and detail.
    pub fn new(code: &'static str, detail: impl Into<String>) -> Self {
        Self {
            code,
            detail: detail.into(),
        }
    }
}

impl fmt::Display for Refusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.code, self.detail)
    }
}

impl std::error::Error for Refusal {}

/// The value of an operation together with the marks of every not-exact
/// decision met on the way.
#[derive(Clone, Debug)]
pub struct Answer<T> {
    /// The value.
    pub value: T,
    /// Marks of not-exact decisions; empty when the answer is exact.
    pub marks: MarkLog,
}

impl<T> Answer<T> {
    /// An exact answer with an empty mark log.
    pub fn exact(value: T) -> Self {
        Self {
            value,
            marks: MarkLog::new(),
        }
    }

    /// An answer with these marks.
    pub fn with_marks(value: T, marks: MarkLog) -> Self {
        Self { value, marks }
    }

    /// Whether the mark log is empty.
    pub fn is_exact(&self) -> bool {
        self.marks.is_empty()
    }

    /// Maps the value, keeping the marks.
    pub fn map<U>(self, f: impl FnOnce(T) -> U) -> Answer<U> {
        Answer {
            value: f(self.value),
            marks: self.marks,
        }
    }
}

/// The result of an operation: an answer with marks, or a refusal.
pub type Answered<T> = Result<Answer<T>, Refusal>;
