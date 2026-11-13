//! Rule-bound answers without an exact proof. No operation selects a rule yet.

use std::fmt;

use crate::{
    Answer, Answered, CertificateIdentity, CertificateScope, EvidenceBound, EvidenceUnit, Rational,
    Refusal, Sign,
};

/// A versioned rule fixed before evaluating a named kind of hard case.
/// Rules belong in code, not in requests or result-dependent recovery logic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BoundedRule {
    /// Stable name of the kind of case.
    pub case_kind: &'static str,
    /// Stable, versioned rule name.
    pub name: &'static str,
    /// Exact rational upper limit, in `unit`. Acceptance includes equality.
    pub level: &'static str,
    /// Units shared by the level and the measured bound.
    pub unit: EvidenceUnit,
}

/// Checked bound against an immutable kept definition, without an exact proof.
/// This checks the record's structure and level, not the measurement method.
/// The producer must establish the bound over the entire named scope.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundedEvidence {
    rule: &'static BoundedRule,
    level: Rational,
    measured: Rational,
    definition: CertificateIdentity<String>,
}

impl BoundedEvidence {
    /// Records the bound under a static rule selected before doing the work.
    /// A residual with unknown units or a bound above the rule is not accepted.
    pub fn new(
        rule: &'static BoundedRule,
        definition: CertificateIdentity<String>,
        measured: EvidenceBound,
    ) -> Result<Self, Refusal> {
        let invalid = |detail| Refusal::new("numbers.invalid-bounded-answer", detail);
        let level = Rational::parse(rule.level)
            .map_err(|_| invalid("bounded rule level is not an exact rational"))?;
        if rule.case_kind.trim().is_empty()
            || rule.name.trim().is_empty()
            || level.sign() == Sign::Negative
            || measured.upper.sign() == Sign::Negative
            || rule.unit == EvidenceUnit::Raw
            || measured.unit != rule.unit
        {
            return Err(invalid(
                "bounded answer needs a named rule and nonnegative bounds in known matching units",
            ));
        }
        if measured.upper > level {
            return Err(Refusal::new(
                "numbers.bounded-level-exceeded",
                "measured bound exceeds the preset case level",
            ));
        }
        if [
            &definition.revision,
            &definition.representation,
            &definition.target,
        ]
        .iter()
        .any(|key| key.trim().is_empty())
        {
            return Err(invalid(
                "bounded answer must identify its kept definition and revisions",
            ));
        }
        definition.scope.check()?;
        let empty_scope_key = match &definition.scope {
            CertificateScope::WholeSolid => false,
            CertificateScope::Entities(keys) => keys.iter().any(|key| key.trim().is_empty()),
            CertificateScope::ParameterBox { chart, .. } => chart.trim().is_empty(),
            CertificateScope::Query(key) => key.trim().is_empty(),
        };
        if empty_scope_key {
            return Err(invalid(
                "bounded answer scope must identify its kept domain",
            ));
        }
        Ok(Self {
            rule,
            level,
            measured: measured.upper,
            definition,
        })
    }

    /// Preset case rule, unchanged by the measurement.
    pub fn rule(&self) -> &'static BoundedRule {
        self.rule
    }
    /// Exact preset level.
    pub fn level(&self) -> &Rational {
        &self.level
    }
    /// Established upper bound against the kept definition.
    pub fn measured(&self) -> &Rational {
        &self.measured
    }
    /// Kept definition, representation revision and complete measured scope.
    pub fn definition(&self) -> &CertificateIdentity<String> {
        &self.definition
    }
}

impl fmt::Display for BoundedEvidence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "bounded: {} under {}: {} <= {} {:?} against {} at {}",
            self.rule.case_kind,
            self.rule.name,
            self.measured,
            self.level,
            self.rule.unit,
            self.definition.target,
            self.definition.revision
        )
    }
}

/// An answer that meets a preset case level without an exact proof.
/// It has no conversion to [`Answer`], even when its measured bound is zero.
///
/// ```compile_fail
/// use catallax_numbers::{Answer, BoundedAnswer};
/// fn promote(value: BoundedAnswer<u32>) -> Answer<u32> { value.into() }
/// ```
#[must_use]
#[derive(Clone, Debug)]
pub struct BoundedAnswer<T> {
    value: T,
    evidence: Box<BoundedEvidence>,
}

impl<T> BoundedAnswer<T> {
    /// Attaches checked whole-result evidence to the measured value.
    pub fn new(value: T, evidence: BoundedEvidence) -> Self {
        Self {
            value,
            evidence: Box::new(evidence),
        }
    }
    /// The bounded value. Reading it does not grant proof status.
    pub fn value(&self) -> &T {
        &self.value
    }
    /// The rule, level and bound against the kept definition.
    pub fn evidence(&self) -> &BoundedEvidence {
        &self.evidence
    }
}

/// Operation delivery, with bounded answers separate from proven answers.
/// Existing operations keep [`Answered`] until they implement a named rule.
#[must_use]
#[derive(Clone, Debug)]
pub enum OperationResult<T> {
    /// Existing proof-bearing answer, including its diagnostic marks.
    Proven(Answer<T>),
    /// Meets its preset case level without an exact proof.
    Bounded(BoundedAnswer<T>),
    /// No acceptable whole-result bound was established.
    Refused(Refusal),
}

impl<T> From<Answered<T>> for OperationResult<T> {
    fn from(result: Answered<T>) -> Self {
        match result {
            Ok(answer) => Self::Proven(answer),
            Err(refusal) => Self::Refused(refusal),
        }
    }
}

impl<T> OperationResult<T> {
    /// Only a proven variant supplies an existing proof-bearing answer.
    pub fn proven(&self) -> Option<&Answer<T>> {
        match self {
            Self::Proven(answer) => Some(answer),
            _ => None,
        }
    }
}

impl<T: fmt::Display> fmt::Display for OperationResult<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Proven(answer) => write!(f, "proven: {}", answer.value),
            Self::Bounded(answer) => write!(f, "{}; value: {}", answer.evidence, answer.value),
            Self::Refused(refusal) => write!(f, "refused: {refusal}"),
        }
    }
}
