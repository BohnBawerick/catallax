//! Piece 1 of catallax: numbers and decisions.
//!
//! Exact integers, rationals and real algebraic numbers; proven enclosures
//! (`Interval`, `Ball`); polynomials; certified roots; the decision interface
//! with its four outcomes; time budgets; construction facts; and marks.
//!
//! The contract is `docs/pieces/01-numbers-and-decisions.md`.

#[macro_use]
mod macros;

mod ahead;
pub use ahead::{
    AheadSite, Charges, Helpers, RunAhead, private_cache, record, run_ahead_adopted,
    run_ahead_waits, set_run_ahead_loops, set_run_ahead_wait, set_threads, shared_caches, threads,
};
mod algebraic;
mod angle;
mod ball;
mod certificate;
mod decider;
mod facts;
mod generator;
mod height;
mod identity;
mod integer;
mod interval;
mod marks;
mod operation;
mod outcome;
pub use operation::{
    OperationAllowance, OperationWork, current_operation_work, operation_budget_active,
    prefer_reduced_angle_series, reuse_operation_enclosures, with_operation_budget,
};
mod poly;
#[cfg(feature = "work-profile")]
mod profile;
#[cfg(feature = "work-profile")]
pub use profile::{WorkProfileRow, with_work_profile};
mod quadratic;
mod rational;
mod real;
mod real_poly;
mod real_poly_roots;
mod structure;
mod system;
mod tower;
mod trace;
pub use trace::{OperationEntry, OperationSpan, with_operation_log};

#[cfg(test)]
mod tests;

pub use algebraic::{Algebraic, RealRoot};
pub use angle::sin_cos_degrees;
pub use ball::Ball;
pub use certificate::{
    CertificateIdentity, CertificateScope, DecisionBook, DecisionEvidence, EvidenceBook,
    EvidenceBound, EvidenceUnit, ExceptionSource, FeatureSizeLowerBound, GEOMETRIC_GAP_EXCEEDED,
    GEOMETRIC_GAP_UNPROVEN, GapCertificate, GeometricCertificate, GeometricEvidence,
    GeometricProof, ImportToleranceAcceptance, LocalException, LocalToleranceSettings,
    ModellingTolerance, ProofContext, RecordedChoice, RegionCertificate,
};
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
