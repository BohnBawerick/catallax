//! Records of numerical checks. These records carry no complete geometric bound.

use crate::{Rational, Refusal, Sign};

/// One evaluated parameter location, not an enclosure of the surrounding domain.
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub struct NumericalSample {
    /// Curve parameter, optionally followed by surface u and v.
    pub parameters: Vec<Rational>,
    /// Upper distance at this location only, in millimetres.
    pub distance_upper_mm: Rational,
}

/// A check performed by a numerical producer. Unperformed checks have no entry.
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub enum NumericalCheck {
    /// A discrete check, such as closed edge uses or a simple parameter polygon.
    Passed {
        /// Stable check name.
        name: String,
        /// Entity or body to which the check applies.
        scope: String,
    },
    /// Distances at all the listed locations. No assertion between samples.
    Samples {
        /// Stable check name.
        name: String,
        /// Entity or incidence to which these locations belong.
        scope: String,
        /// Actual evaluated locations and their distance enclosures.
        samples: Vec<NumericalSample>,
    },
}

/// An actual repair, retained as history rather than a bound on a later body.
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub enum NumericalRepair {
    /// A vertex relocation. The distance applies only to this point movement.
    VertexMoved {
        /// Stable vertex key.
        vertex: String,
        /// Position before the repair, in millimetres.
        before: [Rational; 3],
        /// Position after the repair, in millimetres.
        after: [Rational; 3],
        /// Enclosed distance between these two positions.
        distance_upper_mm: Rational,
    },
}

/// Tolerance, actual checks and repair history for one immutable numerical result.
///
/// Unlike [`crate::BoundedEvidence`], this makes no whole-face or whole-body
/// deviation claim. Even zero at every sample does not grant proof status.
/// This constructor validates a record; the geometry producer must run the checks.
///
/// ```compile_fail
/// use catallax_numbers::{BoundedEvidence, NumericalEvidence};
/// fn promote(evidence: NumericalEvidence) -> BoundedEvidence { evidence.into() }
/// ```
#[derive(Clone, Debug, PartialEq, Eq, catallax_numbers::Transfer)]
pub struct NumericalEvidence {
    tolerance_mm: Rational,
    checks: Vec<NumericalCheck>,
    repairs: Vec<NumericalRepair>,
}

impl NumericalEvidence {
    /// Validates nonempty named checks and their sampled distances under a stated
    /// positive tolerance. A failed or missing required check belongs in a refusal.
    pub fn new(
        tolerance_mm: Rational,
        checks: Vec<NumericalCheck>,
        repairs: Vec<NumericalRepair>,
    ) -> Result<Self, Refusal> {
        let invalid = || {
            Refusal::new(
                "numbers.invalid-numerical-evidence",
                "numerical evidence needs a positive tolerance and nonempty named checks",
            )
        };
        if tolerance_mm.sign() != Sign::Positive || checks.is_empty() {
            return Err(invalid());
        }
        for check in &checks {
            let (name, scope) = match check {
                NumericalCheck::Passed { name, scope } => (name, scope),
                NumericalCheck::Samples {
                    name,
                    scope,
                    samples,
                } => {
                    if samples.is_empty() {
                        return Err(invalid());
                    }
                    for sample in samples {
                        if sample.parameters.is_empty()
                            || sample.distance_upper_mm.sign() == Sign::Negative
                        {
                            return Err(invalid());
                        }
                        if sample.distance_upper_mm > tolerance_mm {
                            return Err(Refusal::new(
                                "numbers.numerical-tolerance-exceeded",
                                format!(
                                    "{name} on {scope}: sampled distance {} mm exceeds {} mm",
                                    sample.distance_upper_mm, tolerance_mm
                                ),
                            ));
                        }
                    }
                    (name, scope)
                }
            };
            if name.trim().is_empty() || scope.trim().is_empty() {
                return Err(invalid());
            }
        }
        for repair in &repairs {
            let NumericalRepair::VertexMoved {
                vertex,
                before,
                after,
                distance_upper_mm,
            } = repair;
            let squared = (0..3).fold(Rational::zero(), |sum, i| {
                let d = &after[i] - &before[i];
                sum + &d * &d
            });
            if vertex.trim().is_empty()
                || distance_upper_mm.sign() == Sign::Negative
                || distance_upper_mm * distance_upper_mm < squared
            {
                return Err(invalid());
            }
        }
        Ok(Self {
            tolerance_mm,
            checks,
            repairs,
        })
    }

    /// The tolerance used for this result's checks, in millimetres.
    pub fn tolerance_mm(&self) -> &Rational {
        &self.tolerance_mm
    }
    /// The checks actually completed for this result.
    pub fn checks(&self) -> &[NumericalCheck] {
        &self.checks
    }
    /// Actions on this result's history. Empty means no repair occurred.
    /// Historical distances are not composed into a bound on the current body.
    pub fn repairs(&self) -> &[NumericalRepair] {
        &self.repairs
    }
}
