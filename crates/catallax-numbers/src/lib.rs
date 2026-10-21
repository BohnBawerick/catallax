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
mod ball;
mod decider;
mod facts;
mod generator;
mod height;
mod integer;
mod interval;
mod marks;
mod outcome;
mod poly;
mod quadratic;
mod rational;
mod real;
mod system;
mod tower;

#[cfg(test)]
mod tests;

pub use algebraic::{Algebraic, RealRoot};
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
pub use system::{BoxRoots, IsolatedRoot, PolySystem, RootBox};
