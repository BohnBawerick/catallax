//! Piece 1 of catallax: numbers and decisions.
//!
//! Exact integers, rationals and real algebraic numbers; proven enclosures
//! (`Interval`, `Ball`); polynomials; certified roots; the decision interface
//! with its four outcomes; time budgets; construction facts; and marks.
//!
//! The contract is `docs/pieces/01-numbers-and-decisions.md`.

#[macro_use]
mod macros;

mod algebraic;
mod angle;
mod ball;
mod decider;
mod facts;
mod generator;
mod height;
mod integer;
mod interval;
mod marks;
mod operation;
mod outcome;
pub use operation::{OperationWork, current_operation_work, with_operation_budget};
mod poly;
mod quadratic;
mod rational;
mod real;
mod real_poly;
mod real_poly_roots;
mod system;
mod tower;
mod trace;
pub use trace::{OperationEntry, OperationSpan, with_operation_log};

#[cfg(test)]
mod tests;

pub use algebraic::{Algebraic, RealRoot};
pub use angle::sin_cos_degrees;
pub use ball::Ball;
pub use decider::{BudgetPolicy, Decider, DecisionStats};
pub use facts::{Fact, FactBook, FactKind, FactOrigin};
pub use integer::{Integer, ParseErrorKind, ParseNumberError, Sign};
pub use interval::Interval;
pub use marks::{BoundKind, Mark, MarkLog};
pub use outcome::{
    Answer, Answered, Coincidence, NotExact, Outcome, Refusal, StrictOrder, StrictSign, ZeroProof,
};
pub use poly::{MPoly, Monomial, UPoly};
pub use rational::Rational;
pub use real::Real;
pub use real_poly::RealPoly;
pub use system::{BoxRoots, IsolatedRoot, PolySystem, RootBox};
