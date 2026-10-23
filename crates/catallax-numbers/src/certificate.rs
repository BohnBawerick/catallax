//! Geometric proof records, separate from diagnostic marks.
//!
//! These are proof claims, not authority to change geometry. Structural checking
//! verifies units, bounds, identities and dependencies. The geometric consumer
//! must also verify each named theorem against the immutable definitions.

use std::collections::{BTreeMap, BTreeSet, btree_map::Entry};

use crate::{Decider, Rational, Real, Refusal, Sign};

/// A geometric conversion, scope or whole-result proof is missing.
pub const GEOMETRIC_GAP_UNPROVEN: &str = "numbers.geometric-gap-unproven";
/// Available decision or local geometric evidence does not meet its named limit.
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
    /// Complete oriented whole-result cover, optionally with a distance summary.
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
    /// Explicit representation beside an immutable definition, with a checked
    /// whole-domain approximation or imported-geometry correspondence theorem.
    Approximation {
        /// Exact definition being represented, not yesterday's fitted controls.
        definition: K,
        /// Producer-owned correspondence, branch and neighbourhood theorem.
        theorem: K,
    },
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
    /// Optional geometric upper bound in millimetres. A complete result cover
    /// need not assert a uniform metric; absence never means zero displacement.
    pub bound: Option<EvidenceBound>,
    /// Producer-owned geometric proof obligation.
    pub proof: GeometricProof<K>,
    /// Shared dependencies in canonical identity order.
    pub dependencies: BTreeSet<K>,
    /// Producing piece, step and stage.
    pub context: ProofContext<K>,
}

/// A scoped gap claim under the strict 1 nm decision/legacy-output rule.
/// This is separate from a whole-result certificate containing local exceptions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GapCertificate<K> {
    /// Immutable target and representation binding.
    pub identity: CertificateIdentity<K>,
    /// Nonnegative geometric upper bound in millimetres.
    pub bound: EvidenceBound,
    /// Root evidence identity, never an empty-log exactness heuristic.
    pub evidence: K,
}

/// Fixed local approximation policy recorded before proof work starts.
/// Only a file-specific explicit import acceptance can relax its import limit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalToleranceSettings {
    kernel_target_mm: Rational,
    kernel_cap_mm: Rational,
    import_cap_mm: Rational,
}

impl Default for LocalToleranceSettings {
    fn default() -> Self {
        Self {
            kernel_target_mm: Rational::from_ratio(1, 10_000).unwrap(),
            kernel_cap_mm: Rational::from_ratio(1, 1_000).unwrap(),
            import_cap_mm: Rational::from_ratio(1, 100).unwrap(),
        }
    }
}

impl LocalToleranceSettings {
    /// Engineering fit target, 100 nm; exceeding it is not itself a cap failure.
    pub fn kernel_target_mm(&self) -> &Rational {
        &self.kernel_target_mm
    }
    /// Inclusive kernel-made approximation cap, 1 micrometre.
    pub fn kernel_cap_mm(&self) -> &Rational {
        &self.kernel_cap_mm
    }
    /// Inclusive default import cap, 10 micrometres.
    pub fn import_cap_mm(&self) -> &Rational {
        &self.import_cap_mm
    }
}

/// Explicit acceptance belonging to one immutable imported file revision.
/// The consumer must verify that `record` is the recipe's actual authorisation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ImportToleranceAcceptance<K> {
    /// Accepted file revision, never a global or later-file override.
    pub file: K,
    /// Accepted geometric limit with explicit millimetre units.
    pub limit: EvidenceBound,
    /// Immutable authorisation record.
    pub record: K,
}

/// Why an entity is an explicit exception to exact agreement.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExceptionSource<K> {
    /// A fitted representation of a retained exact construction law.
    KernelApproximation,
    /// Geometry declared by an imported file, not an invented intended design.
    Import {
        /// Immutable source file revision.
        file: K,
        /// Source uncertainty, distinct from deviation from declared geometry.
        source_uncertainty: Option<EvidenceBound>,
        /// Explicit file-specific acceptance of a looser limit.
        acceptance: Option<ImportToleranceAcceptance<K>>,
    },
    /// Effect of one bounded choice; it cannot borrow the fitting allowance.
    BoundedDecision {
        /// Original converted decision evidence, also listed in the result.
        decision: K,
    },
}

/// A claimed neighbouring-size lower bound, covered by the region's theorem.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FeatureSizeLowerBound {
    /// Explicit claim that this kind of neighbour is absent. The geometric
    /// consumer must check adjacency; this never stands for an unknown size.
    NoNeighbour,
    /// Positive proved lower bound, in millimetres.
    Millimetres(Rational),
}

/// Local deviation and limits for a complete represented trimmed domain.
/// These public records are claims; arithmetic checks are not geometry proofs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalException<K> {
    /// Lower endpoint of the attained deviation enclosure, in millimetres.
    pub lower_mm: Rational,
    /// Upper endpoint over both corresponding domains, including boundaries.
    pub bound: EvidenceBound,
    /// Applicable source hard limit as recorded when the work was requested.
    pub limit: EvidenceBound,
    /// Approximation, import or bounded-choice provenance.
    pub source: ExceptionSource<K>,
    /// Lower bound on the shortest neighbouring edge.
    pub shortest_edge: FeatureSizeLowerBound,
    /// Lower bound on the narrowest neighbouring face.
    pub narrowest_face: FeatureSizeLowerBound,
}

impl<K: Ord> LocalException<K> {
    /// Checks units, interval order, source cap and both tenth-of-feature limits.
    /// A fit above 100 nm but within its hard limits passes; inspect
    /// `kernel_target_met` separately. No value or setting is ever clamped.
    pub fn check(
        &self,
        tolerance: &ModellingTolerance,
        settings: &LocalToleranceSettings,
    ) -> Result<(), Refusal> {
        self.bound.check_millimetres()?;
        self.limit.check_millimetres()?;
        if self.lower_mm < Rational::zero() || self.lower_mm > self.bound.upper {
            return Err(unproven("invalid local deviation enclosure"));
        }
        let cap = match &self.source {
            ExceptionSource::KernelApproximation => settings.kernel_cap_mm(),
            ExceptionSource::BoundedDecision { .. } => tolerance.millimetres(),
            ExceptionSource::Import {
                file,
                source_uncertainty,
                acceptance,
            } => {
                if let Some(uncertainty) = source_uncertainty {
                    uncertainty.check_millimetres()?;
                }
                if let Some(accepted) = acceptance {
                    accepted.limit.check_millimetres()?;
                    if &accepted.file != file || accepted.limit.upper < *settings.import_cap_mm() {
                        return Err(unproven(
                            "import acceptance is not a looser limit for this file",
                        ));
                    }
                    &accepted.limit.upper
                } else {
                    settings.import_cap_mm()
                }
            }
        };
        if self.limit.upper != *cap {
            return Err(unproven(
                "recorded local source limit differs from its applicable policy",
            ));
        }
        if matches!(self.source, ExceptionSource::BoundedDecision { .. }) {
            tolerance.check(&self.bound)?;
        } else if self.bound.upper > *cap {
            return Err(Refusal::new(
                GEOMETRIC_GAP_EXCEEDED,
                format!(
                    "available local upper bound {} mm exceeds inclusive source cap {} mm",
                    self.bound.upper, cap
                ),
            ));
        }
        let ten_upper = &self.bound.upper * &Rational::from_i64(10);
        for (name, size) in [
            ("edge", &self.shortest_edge),
            ("face", &self.narrowest_face),
        ] {
            if let FeatureSizeLowerBound::Millimetres(lower) = size {
                if lower <= &Rational::zero() {
                    return Err(unproven(
                        "neighbouring feature has no positive size lower bound",
                    ));
                }
                if &ten_upper > lower {
                    return Err(Refusal::new(
                        GEOMETRIC_GAP_EXCEEDED,
                        format!(
                            "available local upper bound {} mm exceeds one tenth of the neighbouring {name} lower bound {lower} mm",
                            self.bound.upper
                        ),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Numeric engineering-target status; `None` for non-kernel exceptions.
    /// This does not replace cap or geometric proof checking.
    pub fn kernel_target_met(&self, settings: &LocalToleranceSettings) -> Option<bool> {
        matches!(self.source, ExceptionSource::KernelApproximation)
            .then(|| self.bound.upper <= *settings.kernel_target_mm())
    }
}

/// One exact region or explicit local exception in a complete result cover.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionCertificate<K> {
    /// Exact definition, immutable representation and complete trimmed domain.
    pub identity: CertificateIdentity<K>,
    /// Checked geometry theorem and shared dependencies.
    pub evidence: K,
    /// `None` requires a zero-displacement proof, not an empty mark log.
    pub exception: Option<LocalException<K>>,
}

/// Complete exact/exception cover and oriented material correspondence claim.
/// The geometric consumer must verify topology and the named theorems.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeometricCertificate<K> {
    /// Whole represented body and immutable recipe target.
    pub identity: CertificateIdentity<K>,
    /// Complete oriented correspondence and validity evidence.
    pub evidence: K,
    /// Every actual region exactly once, keyed by its representation identity.
    pub regions: BTreeMap<K, RegionCertificate<K>>,
    /// Original converted bounded choices; each retains the strict 1 nm rule.
    pub decisions: BTreeSet<K>,
    /// Optional proved uniform signed-distance summary, not a sub-1-nm gate.
    /// A maximum of local bounds alone cannot establish this claim.
    pub uniform_bound: Option<EvidenceBound>,
}

/// Append-only shared proof records, ordered by stable definition identity.
#[derive(Clone, Debug)]
pub struct EvidenceBook<K> {
    tolerance: ModellingTolerance,
    local_settings: LocalToleranceSettings,
    records: BTreeMap<K, GeometricEvidence<K>>,
}

impl<K: Ord> EvidenceBook<K> {
    /// Starts a record with the modelling setting fixed before any proof work.
    pub fn new(tolerance: ModellingTolerance) -> Self {
        Self {
            tolerance,
            local_settings: LocalToleranceSettings::default(),
            records: BTreeMap::new(),
        }
    }

    /// The immutable modelling setting attached to this record.
    pub fn tolerance(&self) -> &ModellingTolerance {
        &self.tolerance
    }

    /// Immutable local approximation policy recorded with this evidence book.
    pub fn local_settings(&self) -> &LocalToleranceSettings {
        &self.local_settings
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

    /// Checks a scoped gap's binding, strict cap and dependency graph.
    /// This does not verify its geometric theorem or authorize publication.
    pub fn check_gap_structure(
        &self,
        certificate: &GapCertificate<K>,
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
        let bound = root
            .bound
            .as_ref()
            .ok_or_else(|| unproven("gap evidence has no geometric bound"))?;
        bound.check_millimetres()?;
        if certificate.bound.upper < bound.upper {
            return Err(unproven("certificate understates its root evidence bound"));
        }
        self.check_graph(&certificate.evidence, dec).map(|_| ())
    }

    fn check_graph<'a>(&'a self, root: &'a K, dec: &Decider) -> Result<BTreeSet<&'a K>, Refusal> {
        let mut active = BTreeSet::new();
        let mut checked = BTreeSet::new();
        let mut pending = vec![(root, false)];
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
            if let Some(bound) = &record.bound {
                bound.check_millimetres()?;
            }
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
        Ok(checked)
    }

    /// Checks a whole-result claim against the actual representation inventory.
    /// Source caps are local and inclusive; decisions remain strictly below 1 nm.
    /// This checks structure and arithmetic only. The caller still verifies
    /// geometry, adjacency, exact definitions, complete coverage and orientation.
    pub fn check_structure(
        &self,
        certificate: &GeometricCertificate<K>,
        expected_regions: &BTreeSet<K>,
        dec: &Decider,
    ) -> Result<(), Refusal> {
        dec.check_deadline()?;
        if !matches!(certificate.identity.scope, CertificateScope::WholeSolid) {
            return Err(unproven(
                "result certificate does not cover the whole solid",
            ));
        }
        if certificate.regions.len() != expected_regions.len()
            || certificate
                .regions
                .keys()
                .any(|key| !expected_regions.contains(key))
        {
            return Err(unproven(
                "result certificate omits or invents a represented region",
            ));
        }
        let checked = self.check_graph(&certificate.evidence, dec)?;
        let root = self.records.get(&certificate.evidence).unwrap();
        if root.identity != certificate.identity
            || !matches!(root.proof, GeometricProof::WholeResult(_))
        {
            return Err(unproven(
                "result certificate is not bound to its whole correspondence proof",
            ));
        }
        if certificate.uniform_bound != root.bound {
            return Err(unproven(
                "uniform metric summary does not match its complete proof",
            ));
        }
        for decision in &certificate.decisions {
            if !checked.contains(decision) {
                return Err(unproven(
                    "result lists a decision absent from its proof graph",
                ));
            }
            let record = self.records.get(decision).unwrap();
            if !matches!(record.proof, GeometricProof::Conversion { .. }) {
                return Err(unproven(
                    "bounded choice is not original converted decision evidence",
                ));
            }
            self.tolerance.check(
                record
                    .bound
                    .as_ref()
                    .ok_or_else(|| unproven("bounded choice has no geometric gap proof"))?,
            )?;
        }
        for (key, region) in &certificate.regions {
            dec.check_deadline()?;
            if !checked.contains(&region.evidence)
                || region.identity.revision != certificate.identity.revision
                || &region.identity.representation != key
                || !matches!(&region.identity.scope,
                    CertificateScope::Entities(entities) if entities.len() == 1 && entities.contains(key))
            {
                return Err(unproven("region binding or complete domain is missing"));
            }
            let record = self.records.get(&region.evidence).unwrap();
            if record.identity != region.identity {
                return Err(unproven(
                    "region evidence names a different definition or representation",
                ));
            }
            let bound = record
                .bound
                .as_ref()
                .ok_or_else(|| unproven("region correspondence has no geometric bound"))?;
            match &region.exception {
                None if !bound.upper.is_zero() => {
                    return Err(unproven(
                        "changed region is not listed as a local exception",
                    ));
                }
                None => {}
                Some(exception) => {
                    exception.check(&self.tolerance, &self.local_settings)?;
                    if exception.bound.upper < bound.upper {
                        return Err(unproven("local record understates its geometric evidence"));
                    }
                    match &exception.source {
                        ExceptionSource::KernelApproximation | ExceptionSource::Import { .. } => {
                            if !matches!(&record.proof, GeometricProof::Approximation { definition, .. }
                                if definition == &region.identity.target)
                            {
                                return Err(unproven(
                                    "fit/import exception lacks its retained-definition correspondence",
                                ));
                            }
                        }
                        ExceptionSource::BoundedDecision { decision } => {
                            if !certificate.decisions.contains(decision)
                                || (&region.evidence != decision
                                    && !record.dependencies.contains(decision))
                            {
                                return Err(unproven(
                                    "decision-only exception lacks its original bounded choice",
                                ));
                            }
                        }
                    }
                }
            }
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
