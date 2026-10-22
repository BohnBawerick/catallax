//! Geometric proof records, separate from diagnostic marks.
//!
//! These are proof claims, not authority to change geometry. Structural checking
//! verifies units, bounds, identities and dependencies. The geometric consumer
//! must also verify each named theorem against the immutable definitions.

use std::collections::{BTreeMap, BTreeSet, btree_map::Entry};

use crate::{Decider, Rational, Real, Refusal, Sign};

/// A geometric conversion, scope or whole-result proof is missing.
pub const GEOMETRIC_GAP_UNPROVEN: &str = "numbers.geometric-gap-unproven";
/// The available geometric upper bound does not meet the strict modelling cap.
pub const GEOMETRIC_GAP_EXCEEDED: &str = "numbers.geometric-gap-exceeded";

fn unproven(detail: impl Into<String>) -> Refusal {
    Refusal::new(GEOMETRIC_GAP_UNPROVEN, detail)
}

/// The modelling sameness setting, fixed before an operation starts.
/// Output accuracy and entity-local approximation caps are separate settings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModellingTolerance {
    millimetres: Rational,
}

impl Default for ModellingTolerance {
    fn default() -> Self {
        Self {
            millimetres: Rational::from_ratio(1, 1_000_000).unwrap(),
        }
    }
}

impl ModellingTolerance {
    /// The recorded 1 nm value, in millimetres. There is no growth setter.
    pub fn millimetres(&self) -> &Rational {
        &self.millimetres
    }

    /// Checks a nonnegative geometric bound strictly below the recorded cap.
    /// This checks arithmetic only, not the theorem that produced the bound.
    pub fn check(&self, bound: &EvidenceBound) -> Result<(), Refusal> {
        bound.check_millimetres()?;
        if bound.upper >= self.millimetres {
            return Err(Refusal::new(
                GEOMETRIC_GAP_EXCEEDED,
                format!(
                    "available geometric upper bound {} mm must be strictly below {} mm",
                    bound.upper, self.millimetres
                ),
            ));
        }
        Ok(())
    }
}

/// Units declared by a geometric caller, never inferred from a `Real` expression.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EvidenceUnit {
    /// Distance or coordinate in millimetres.
    Millimetres,
    /// Angle in radians.
    Radians,
    /// Dimensionless scalar parameter.
    Scalar,
    /// Area in square millimetres.
    SquareMillimetres,
    /// Volume in cubic millimetres.
    CubicMillimetres,
    /// A residual whose geometric units have not been established.
    Raw,
}

/// An outward rational magnitude bound and its declared units.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EvidenceBound {
    /// Nonnegative upper bound. Imported records must be checked before use.
    pub upper: Rational,
    /// Meaning of the quantity, supplied by the producing piece.
    pub unit: EvidenceUnit,
}

impl EvidenceBound {
    /// Encloses the magnitude of a shared expression under the caller's budget.
    /// Declaring units does not prove a residual-to-geometry conversion.
    pub fn enclose_magnitude(
        value: &Real,
        unit: EvidenceUnit,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if width.sign() != Sign::Positive {
            return Err(unproven("evidence enclosure width must be positive"));
        }
        let enclosed = dec.enclose(value, width)?;
        Ok(Self {
            upper: enclosed.lo().abs().max(enclosed.hi().abs()),
            unit,
        })
    }

    /// Adds consecutive bounds in the same units without rounding or clamping.
    /// The caller must separately prove that addition describes its geometry.
    pub fn add(&self, other: &Self) -> Result<Self, Refusal> {
        self.check_nonnegative()?;
        other.check_nonnegative()?;
        if self.unit != other.unit {
            return Err(unproven("cannot add evidence bounds with different units"));
        }
        Ok(Self {
            upper: &self.upper + &other.upper,
            unit: self.unit,
        })
    }

    fn check_nonnegative(&self) -> Result<(), Refusal> {
        if self.upper.sign() == Sign::Negative {
            return Err(unproven("negative evidence upper bound"));
        }
        Ok(())
    }

    fn check_millimetres(&self) -> Result<(), Refusal> {
        self.check_nonnegative()?;
        if self.unit != EvidenceUnit::Millimetres {
            return Err(unproven("geometric gap evidence must be in millimetres"));
        }
        Ok(())
    }
}

/// The complete domain of a claim. A witness alone cannot stand for a face.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CertificateScope<K> {
    /// The oriented whole-solid signed-distance metric over Euclidean space.
    WholeSolid,
    /// The complete named entities, including their finite trims.
    Entities(BTreeSet<K>),
    /// A finite closed parameter box in the named exact chart.
    ParameterBox {
        /// Immutable chart definition identity.
        chart: K,
        /// Ordered lower/upper parameter bounds.
        bounds: Vec<(Rational, Rational)>,
    },
    /// One query, with its exact domain retained under this identity.
    Query(K),
}

impl<K> CertificateScope<K> {
    fn check(&self) -> Result<(), Refusal> {
        match self {
            Self::Entities(entities) if entities.is_empty() => {
                Err(unproven("entity evidence has no covered entities"))
            }
            Self::ParameterBox { bounds, .. }
                if bounds.is_empty() || bounds.iter().any(|(lo, hi)| lo > hi) =>
            {
                Err(unproven("invalid finite evidence parameter domain"))
            }
            _ => Ok(()),
        }
    }
}

/// Immutable identities of both the exact target and its representation.
/// Keys must name recipe/definition revisions, not pointers or completion order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CertificateIdentity<K> {
    /// Exact recipe revision.
    pub revision: K,
    /// Represented body, face, curve or query revision.
    pub representation: K,
    /// Exact recipe-prefix or geometric definition being represented.
    pub target: K,
    /// Entire domain covered by the claim.
    pub scope: CertificateScope<K>,
}

/// Producer location retained independently of the original numeric reason.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProofContext<K> {
    /// Producing piece or checked construction routine.
    pub producer: &'static str,
    /// Recipe step, when applicable.
    pub step: Option<K>,
    /// Boolean stage, when applicable.
    pub stage: Option<&'static str>,
}

/// Original unresolved decision preserved beside a geometric conversion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DecisionEvidence<K> {
    /// Original stable decision code, not replaced by the geometric refusal code.
    pub code: &'static str,
    /// Outward enclosure of the original quantity, in its original units.
    pub enclosure: [Rational; 2],
    /// Original quantity's declared units.
    pub unit: EvidenceUnit,
    /// All affected entities.
    pub entities: BTreeSet<K>,
    /// Chosen branch in the immutable definition.
    pub branch: K,
}

/// Geometric obligation named by a record. No variant is a theorem by assertion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GeometricProof<K> {
    /// Exact construction identity, to be checked against the named definition.
    ExactDefinition(K),
    /// Scoped conversion with conditioning and branch-correspondence evidence.
    Conversion {
        /// Original residual, affected entities and selected branch.
        decision: DecisionEvidence<K>,
        /// The producer's residual-to-geometry theorem identity.
        theorem: K,
    },
    /// Operation-specific transfer from its input certificates.
    Transfer(K),
    /// Complete oriented coverage against the original recipe target.
    WholeResult(K),
}

/// One shared proof dependency. It remains untrusted until its theorem is checked.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeometricEvidence<K> {
    /// Exact target, represented object and covered domain.
    pub identity: CertificateIdentity<K>,
    /// Available upper bound in millimetres.
    pub bound: EvidenceBound,
    /// Producer-owned geometric proof obligation.
    pub proof: GeometricProof<K>,
    /// Shared dependencies in canonical identity order.
    pub dependencies: BTreeSet<K>,
    /// Producing piece, step and stage.
    pub context: ProofContext<K>,
}

/// A certificate claim, separate from marks and from proof verification status.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeometricCertificate<K> {
    /// Immutable target and representation binding.
    pub identity: CertificateIdentity<K>,
    /// Nonnegative geometric upper bound in millimetres.
    pub bound: EvidenceBound,
    /// Root evidence identity, never an empty-log exactness heuristic.
    pub evidence: K,
}

/// Append-only shared proof records, ordered by stable definition identity.
#[derive(Clone, Debug)]
pub struct EvidenceBook<K> {
    tolerance: ModellingTolerance,
    records: BTreeMap<K, GeometricEvidence<K>>,
}

impl<K: Ord> EvidenceBook<K> {
    /// Starts a record with the modelling setting fixed before any proof work.
    pub fn new(tolerance: ModellingTolerance) -> Self {
        Self {
            tolerance,
            records: BTreeMap::new(),
        }
    }

    /// The immutable modelling setting attached to this record.
    pub fn tolerance(&self) -> &ModellingTolerance {
        &self.tolerance
    }

    /// Adds evidence without changing an existing shared proof identity.
    /// Equal repeated records deduplicate; conflicting records refuse.
    pub fn record(
        &mut self,
        id: K,
        evidence: GeometricEvidence<K>,
        dec: &Decider,
    ) -> Result<bool, Refusal> {
        dec.check_deadline()?;
        match self.records.entry(id) {
            Entry::Vacant(entry) => {
                entry.insert(evidence);
                Ok(true)
            }
            Entry::Occupied(entry) if entry.get() == &evidence => Ok(false),
            Entry::Occupied(_) => Err(unproven("conflicting shared geometric proof identity")),
        }
    }

    /// A dependency by its stable identity. Absence is not exactness.
    pub fn get(&self, id: &K) -> Option<&GeometricEvidence<K>> {
        self.records.get(id)
    }

    /// Shared evidence in canonical identity order, independent of insertion order.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &GeometricEvidence<K>)> {
        self.records.iter()
    }

    /// Checks binding, units, the strict cap and the reachable dependency graph.
    /// This does not verify geometric theorems, topology, coverage or orientation,
    /// and must not be used alone to publish a successful geometric result.
    pub fn check_structure(
        &self,
        certificate: &GeometricCertificate<K>,
        dec: &Decider,
    ) -> Result<(), Refusal> {
        dec.check_deadline()?;
        self.tolerance.check(&certificate.bound)?;
        certificate.identity.scope.check()?;
        let root = self
            .records
            .get(&certificate.evidence)
            .ok_or_else(|| unproven("missing root geometric evidence"))?;
        if root.identity != certificate.identity {
            return Err(unproven(
                "certificate is bound to a different target or representation",
            ));
        }
        root.bound.check_millimetres()?;
        if certificate.bound.upper < root.bound.upper {
            return Err(unproven("certificate understates its root evidence bound"));
        }
        let mut active = BTreeSet::new();
        let mut checked = BTreeSet::new();
        let mut pending = vec![(&certificate.evidence, false)];
        while let Some((id, finish)) = pending.pop() {
            dec.check_deadline()?;
            if finish {
                active.remove(id);
                checked.insert(id);
                continue;
            }
            if checked.contains(id) {
                continue;
            }
            if !active.insert(id) {
                return Err(unproven("cyclic geometric proof dependencies"));
            }
            let record = self
                .records
                .get(id)
                .ok_or_else(|| unproven("missing geometric proof dependency"))?;
            record.bound.check_millimetres()?;
            record.identity.scope.check()?;
            if record.context.producer.is_empty() {
                return Err(unproven("geometric evidence has no producing piece"));
            }
            if let GeometricProof::Conversion { decision, .. } = &record.proof
                && (decision.code.is_empty()
                    || decision.enclosure[0] > decision.enclosure[1]
                    || decision.entities.is_empty())
            {
                return Err(unproven("incomplete original decision evidence"));
            }
            pending.push((id, true));
            pending.extend(
                record
                    .dependencies
                    .iter()
                    .rev()
                    .map(|dependency| (dependency, false)),
            );
        }
        Ok(())
    }
}

/// Settled branch and supporting evidence. This is not an exact construction fact.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecordedChoice<K> {
    /// Selected branch in the definition.
    pub branch: K,
    /// Supporting scoped evidence identity.
    pub evidence: K,
}

/// Decide-once record for one immutable recipe revision and work request.
/// Callers include scope, operands and branch selectors in each decision key.
#[derive(Clone, Debug)]
pub struct DecisionBook<K> {
    choices: BTreeMap<K, RecordedChoice<K>>,
}

impl<K: Ord> Default for DecisionBook<K> {
    fn default() -> Self {
        Self {
            choices: BTreeMap::new(),
        }
    }
}

impl<K: Ord> DecisionBook<K> {
    /// Reads a settled choice without repeating its numeric or geometric decision.
    pub fn get(&self, id: &K) -> Option<&RecordedChoice<K>> {
        self.choices.get(id)
    }

    /// Runs the deciding closure only for a new key. A refusal records no choice.
    /// Cancellation is checked even when an earlier choice is reused.
    pub fn decide_once(
        &mut self,
        id: K,
        dec: &mut Decider,
        decide: impl FnOnce(&mut Decider) -> Result<RecordedChoice<K>, Refusal>,
    ) -> Result<&RecordedChoice<K>, Refusal> {
        dec.check_deadline()?;
        match self.choices.entry(id) {
            Entry::Occupied(entry) => Ok(entry.into_mut()),
            Entry::Vacant(entry) => {
                let choice = decide(dec)?;
                dec.check_deadline()?;
                Ok(entry.insert(choice))
            }
        }
    }

    /// Imports a settled choice without allowing a contradictory later branch.
    /// Evidence must be verified by the importing geometric consumer.
    pub fn record(
        &mut self,
        id: K,
        choice: RecordedChoice<K>,
        dec: &Decider,
    ) -> Result<bool, Refusal> {
        dec.check_deadline()?;
        match self.choices.entry(id) {
            Entry::Vacant(entry) => {
                entry.insert(choice);
                Ok(true)
            }
            Entry::Occupied(entry) if entry.get() == &choice => Ok(false),
            Entry::Occupied(_) => Err(unproven("conflicting decide-once branch or evidence")),
        }
    }

    /// Settled choices in canonical decision-key order.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &RecordedChoice<K>)> {
        self.choices.iter()
    }
}
