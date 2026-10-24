//! Opt-in bounded geometry beside immutable recipe definitions.
//! Numeric evidence records remain claims; the complete prism gate supplies
//! their geometric authority before a certified build can leave this module.

mod geometry;
mod source;

use crate::{Built, EllipticOffsetFit, Profile, fit_elliptic_offset};
use catallax_numbers::{
    BoundKind, CertificateIdentity, CertificateScope, Decider, DecisionBook, EvidenceBook,
    EvidenceBound, EvidenceUnit, ExceptionSource, FeatureSizeLowerBound, GEOMETRIC_GAP_EXCEEDED,
    GEOMETRIC_GAP_UNPROVEN, GeometricCertificate, GeometricEvidence, GeometricProof,
    LocalException, Mark, MarkLog, ModellingTolerance, Outcome, ProofContext, Rational, Real,
    RecordedChoice, Refusal, RegionCertificate, Sign, StrictSign,
};
use catallax_recipe::{BooleanKind, Origin, Recipe, StepId};
use catallax_shapes::{
    EllipticOffsetDefinition, SplineSurface, Surface, certify_elliptic_offset_sweep,
};
use catallax_topology::{Body, BodyId, BodyKind, Model};
use geometry::{PrismTopology, ProfileSupport};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    io,
    sync::Arc,
};

/// Content-bound evidence identity; allocation-free keys have a fixed byte order.
pub type CertificateKey = [u8; 32];

pub(super) fn unproven(detail: impl Into<String>) -> Refusal {
    Refusal::new(GEOMETRIC_GAP_UNPROVEN, detail)
}
fn unsupported(detail: impl Into<String>) -> Refusal {
    Refusal::new("features.unsupported-offset-family", detail)
}
fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(reason) => Err(Refusal::new(reason.code, reason.detail)),
        Outcome::Refused(reason) => Err(reason),
    }
}

struct DigestWriter(Sha256);
impl io::Write for DigestWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn revision(recipe: &Recipe) -> Result<CertificateKey, Refusal> {
    let mut digest = DigestWriter(Sha256::new());
    digest.0.update(b"catallax.certified-recipe-definition/1\0");
    serde_json::to_writer(&mut digest, recipe)
        .map_err(|error| unproven(format!("cannot retain canonical recipe identity: {error}")))?;
    Ok(digest.0.finalize().into())
}
fn key(revision: &CertificateKey, step: &StepId, label: &str, index: u32) -> CertificateKey {
    let mut digest = Sha256::new();
    digest.update(b"catallax.elliptic-prism-certificate/1\0");
    digest.update(revision);
    for text in [step.as_str(), label] {
        digest.update((text.len() as u64).to_le_bytes());
        digest.update(text.as_bytes());
    }
    digest.update(index.to_le_bytes());
    digest.finalize().into()
}

#[derive(Clone, Debug)]
struct SourceProof {
    geometry: source::PrismSource,
    evidence: CertificateKey,
    body: BodyId,
    radius: Rational,
}

/// Exact non-growing construction target, independent of the fitted controls.
#[derive(Clone, Debug)]
pub struct PrismTarget {
    source: Arc<SourceProof>,
    definition: EllipticOffsetDefinition,
    axial: [Real; 2],
    empty: bool,
}
impl PrismTarget {
    /// Original analytic support, accumulated signed normal side and exact distance.
    pub fn definition(&self) -> &EllipticOffsetDefinition {
        &self.definition
    }
    /// Exact cap coordinates relative to the original support's unit axis.
    pub fn axial(&self) -> [&Real; 2] {
        self.axial.each_ref()
    }
    /// Exact unit axis, retained from the original circular sweep.
    pub fn axis(&self) -> &catallax_shapes::Vector3 {
        &self.source.geometry.axis
    }
    /// Accumulated inward erosion in millimetres, never fitted controls as input.
    pub fn distance_mm(&self) -> &Rational {
        self.definition
            .distance()
            .rational_ref()
            .expect("recipe offsets retain rational lengths")
    }
    /// Whether the exact regularised target has no three-dimensional material.
    pub fn is_empty(&self) -> bool {
        self.empty
    }
}

#[derive(Clone, Debug)]
enum Representation {
    Analytic,
    Spline(SplineSurface),
    Empty,
}

/// A complete body record with retained definitions and a geometry-owned proof.
/// Public numeric records can be inspected or replayed, but cannot mint this type.
#[derive(Clone, Debug)]
pub struct BodyCertificate {
    body: BodyId,
    step: StepId,
    target: PrismTarget,
    representation: Representation,
    topology: Option<PrismTopology>,
    origins: Arc<[Vec<Origin>; 3]>,
    record: GeometricCertificate<CertificateKey>,
    choice: (CertificateKey, RecordedChoice<CertificateKey>),
    fit_target: Option<Rational>,
}
impl BodyCertificate {
    /// Actual immutable represented body.
    pub fn body(&self) -> BodyId {
        self.body
    }
    /// Producing recipe node.
    pub fn step(&self) -> &StepId {
        &self.step
    }
    /// Immutable exact target, including accumulated erosion and subsequent caps.
    pub fn target(&self) -> &PrismTarget {
        &self.target
    }
    /// Complete face/trim and whole-result evidence claim.
    pub fn record(&self) -> &GeometricCertificate<CertificateKey> {
        &self.record
    }
    /// Genuine checked spline support when this body is approximated.
    pub fn spline(&self) -> Option<&SplineSurface> {
        match &self.representation {
            Representation::Spline(surface) => Some(surface),
            _ => None,
        }
    }
    /// Lower cap, upper cap and side handles; empty results have no faces.
    pub fn faces(&self) -> Option<[catallax_topology::FaceId; 3]> {
        self.topology.as_ref().map(|t| t.faces)
    }
    /// Numeric region bound to this actual face, not an assumed map ordering.
    /// Handles from outside this certificate's complete body cover return `None`.
    pub fn region_for_face(
        &self,
        face: catallax_topology::FaceId,
    ) -> Option<&RegionCertificate<CertificateKey>> {
        if !self.topology.as_ref()?.faces.contains(&face) {
            return None;
        }
        let region = key(
            &self.record.identity.revision,
            &self.step,
            "face",
            face.index(),
        );
        self.record.regions.get(&region)
    }
    /// Shared physical lower/upper rings, not artificial knot or chart edges.
    pub fn rings(&self) -> Option<[catallax_topology::EdgeId; 2]> {
        self.topology.as_ref().map(|t| t.rings)
    }
    /// Producer-selected fitting target, fixed before fitting and never enlarged.
    pub fn fit_target_mm(&self) -> Option<&Rational> {
        self.fit_target.as_ref()
    }

    /// Rechecks an untrusted numeric claim against actual body geometry and the
    /// retained exact law. Understated marks, incomplete scopes and stale choices
    /// do not gain authority from a valid numeric dependency graph.
    pub fn verify_record(
        &self,
        model: &Model,
        candidate: &GeometricCertificate<CertificateKey>,
        evidence: &EvidenceBook<CertificateKey>,
        decisions: &DecisionBook<CertificateKey>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        dec.check_deadline()?;
        if candidate.identity != self.record.identity
            || decisions.get(&self.choice.0) != Some(&self.choice.1)
            || evidence.tolerance() != &ModellingTolerance::default()
        {
            return Err(unproven(
                "body proof has a stale target, representation or decide-once scope",
            ));
        }
        let expected: BTreeSet<_> = self.record.regions.keys().copied().collect();
        evidence.check_structure(candidate, &expected, dec)?;
        let Surface::EllipticCylinder(original) = self.target.definition.original() else {
            unreachable!("retained strict elliptic source");
        };
        match &self.representation {
            Representation::Empty => {
                if !self.target.empty
                    || model.body(self.body).kind != BodyKind::Empty
                    || !model.body(self.body).lumps.is_empty()
                    || !model.body(self.body).marks.is_empty()
                    || !candidate.regions.is_empty()
                    || candidate.uniform_bound.is_some()
                {
                    return Err(unproven(
                        "empty certificate does not describe an exact empty body",
                    ));
                }
            }
            representation => {
                let (profile, bound) = match representation {
                    Representation::Analytic => (
                        ProfileSupport::Analytic,
                        EvidenceBound {
                            upper: Rational::zero(),
                            unit: EvidenceUnit::Millimetres,
                        },
                    ),
                    Representation::Spline(surface) => {
                        if !surface
                            .definition()
                            .original()
                            .same_definition(self.target.definition.original())
                            || surface.definition().side() != self.target.definition.side()
                            || !surface
                                .definition()
                                .distance()
                                .same_definition(self.target.definition.distance())
                            || !surface.axial()[0].same_definition(&self.target.axial[0])
                            || !surface.axial()[1].same_definition(&self.target.axial[1])
                        {
                            return Err(unproven(
                                "checked spline proof belongs to a different retained target or finite domain",
                            ));
                        }
                        // This immutable support was independently checked against
                        // this retained target before the body existed. The body
                        // gate binds its exact controls and complete trims below.
                        (ProfileSupport::Spline(surface), surface.bound().clone())
                    }
                    Representation::Empty => unreachable!("handled empty result"),
                };
                let topology = geometry::check_prism(
                    model,
                    self.body,
                    original,
                    self.target.axial(),
                    profile,
                    &self.origins,
                    dec,
                )?;
                if self.topology.as_ref().is_none_or(|expected| {
                    expected.faces != topology.faces || expected.rings != topology.rings
                }) {
                    return Err(unproven(
                        "body proof is detached from its complete actual topology",
                    ));
                }
                let uniform = candidate
                    .uniform_bound
                    .as_ref()
                    .ok_or_else(|| unproven("prism has no whole-result correspondence bound"))?;
                if uniform.unit != EvidenceUnit::Millimetres || uniform.upper < bound.upper {
                    return Err(unproven(
                        "uniform bound understates the independent whole-body theorem",
                    ));
                }
                let (edge_lower, face_lower) = size_bounds(&self.target, &bound, dec)?;
                for (id, region) in &candidate.regions {
                    if region.identity != self.record.regions[id].identity {
                        return Err(unproven("a face or complete trim domain was rebound"));
                    }
                    if bound.upper.is_zero() {
                        if region.exception.is_some() {
                            return Err(unproven(
                                "exact source acquired an unrelated local exception",
                            ));
                        }
                    } else {
                        let exception = region.exception.as_ref().ok_or_else(|| {
                            unproven("changed cap or side is omitted from the exception cover")
                        })?;
                        if exception.bound.upper < bound.upper
                            || !matches!(exception.source, ExceptionSource::KernelApproximation)
                            || !matches!(&exception.shortest_edge, FeatureSizeLowerBound::Millimetres(value) if value <= &edge_lower)
                            || !matches!(&exception.narrowest_face, FeatureSizeLowerBound::Millimetres(value) if value <= &face_lower)
                        {
                            return Err(unproven(
                                "local metric, fitting provenance or actual neighbouring-size proof is false",
                            ));
                        }
                    }
                }
                let marks = &model.body(self.body).marks;
                if bound.upper.is_zero() {
                    if !marks.is_empty() {
                        return Err(unproven("exact source has an unscoped diagnostic mark"));
                    }
                } else if marks.len() != 1
                    || !marks.iter().all(|mark| {
                        mark.code == "features.certified-offset"
                            && mark.kind == BoundKind::Distance
                            && mark.bound == bound.upper
                    })
                {
                    return Err(unproven(
                        "body mark is missing, forged or detached from its full certificate",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Opt-in recipe build whose final result has a checked complete certificate.
#[derive(Clone, Debug)]
pub struct CertifiedBuilt {
    built: Built,
    evidence: EvidenceBook<CertificateKey>,
    decisions: DecisionBook<CertificateKey>,
    certificates: BTreeMap<StepId, Arc<BodyCertificate>>,
}
impl CertifiedBuilt {
    /// Actual B-rep and unchanged immutable recipe snapshot.
    pub fn built(&self) -> &Built {
        &self.built
    }
    /// Fixed modelling setting and append-only shared evidence graph.
    pub fn evidence(&self) -> &EvidenceBook<CertificateKey> {
        &self.evidence
    }
    /// Settled choices in canonical scope order.
    pub fn decisions(&self) -> &DecisionBook<CertificateKey> {
        &self.decisions
    }
    /// Complete certificate for a bounded/derived body-producing node.
    pub fn certificate(&self, step: &StepId) -> Option<&BodyCertificate> {
        self.certificates.get(step).map(Arc::as_ref)
    }
    /// Final result certificate; construction refuses rather than returning a missing record.
    pub fn result_certificate(&self) -> &BodyCertificate {
        self.certificate(&self.built.recipe().result)
            .expect("certified final result invariant")
    }
    /// Revalidates actual geometry and numeric evidence before a consumer publishes it.
    pub fn verify(&self, dec: &mut Decider) -> Result<(), Refusal> {
        self.result_certificate().verify_record(
            &self.built.model,
            self.result_certificate().record(),
            &self.evidence,
            &self.decisions,
            dec,
        )
    }
}

struct PreparedBody {
    target: PrismTarget,
    representation: Representation,
    body: BodyId,
    origins: Arc<[Vec<Origin>; 3]>,
    choice: (CertificateKey, RecordedChoice<CertificateKey>),
    fit_target: Option<Rational>,
    parent: Option<CertificateKey>,
}

pub(crate) struct CertificationContext<'a> {
    recipe: &'a Recipe,
    snapshot: Option<Arc<Recipe>>,
    revision: Option<CertificateKey>,
    evidence: EvidenceBook<CertificateKey>,
    decisions: DecisionBook<CertificateKey>,
    sources: BTreeMap<StepId, Arc<SourceProof>>,
    certificates: BTreeMap<StepId, Arc<BodyCertificate>>,
}
impl<'a> CertificationContext<'a> {
    pub fn new(recipe: &'a Recipe) -> Self {
        Self {
            recipe,
            snapshot: None,
            revision: None,
            evidence: EvidenceBook::new(ModellingTolerance::default()),
            decisions: DecisionBook::default(),
            sources: BTreeMap::new(),
            certificates: BTreeMap::new(),
        }
    }
    pub fn snapshot(&self) -> Option<Arc<Recipe>> {
        self.snapshot.clone()
    }
    fn retain(&mut self) -> Result<CertificateKey, Refusal> {
        if let Some(revision) = self.revision {
            return Ok(revision);
        }
        let revision = revision(self.recipe)?;
        self.snapshot = Some(Arc::new(self.recipe.clone()));
        self.revision = Some(revision);
        Ok(revision)
    }
    pub fn finish(self, built: Built, dec: &mut Decider) -> Result<CertifiedBuilt, Refusal> {
        if !self.certificates.contains_key(&built.recipe().result) {
            return Err(unsupported(
                "the opt-in final result has no supported complete prism certificate",
            ));
        }
        let result = CertifiedBuilt {
            built,
            evidence: self.evidence,
            decisions: self.decisions,
            certificates: self.certificates,
        };
        result.verify(dec)?;
        Ok(result)
    }

    pub fn offset(
        &mut self,
        model: &mut Model,
        profiles: &BTreeMap<StepId, Profile>,
        step: &StepId,
        input: (&StepId, BodyId, &Rational),
        signed_distance: &Rational,
        dec: &mut Decider,
    ) -> Result<(BodyId, Rational), Refusal> {
        if signed_distance > &Rational::zero() {
            return Err(unsupported(
                "outward solid offsets need rounded edge/corner pieces, not a prism sweep",
            ));
        }
        let revision = self.retain()?;
        let parent = self.certificates.get(input.0).cloned();
        let (source, previous_distance, previous_axial, origins) = if let Some(parent) = &parent {
            if parent.body != input.1 {
                return Err(unproven("offset input body and certificate disagree"));
            }
            (
                parent.target.source.clone(),
                parent.target.distance_mm().clone(),
                parent.target.axial.clone(),
                parent.origins.clone(),
            )
        } else {
            let source = if let Some(source) = self.sources.get(input.0) {
                if source.body != input.1 {
                    return Err(unproven("cached exact source changed body identity"));
                }
                source.clone()
            } else {
                let geometry =
                    source::recognize(self.recipe, input.0, profiles, model, input.1, dec)?;
                let evidence = key(&revision, input.0, "exact-source", input.1.index());
                self.evidence.record(
                    evidence,
                    GeometricEvidence {
                        identity: CertificateIdentity {
                            revision,
                            representation: key(&revision, input.0, "body", input.1.index()),
                            target: key(&revision, input.0, "target", 0),
                            scope: CertificateScope::WholeSolid,
                        },
                        bound: Some(EvidenceBound {
                            upper: Rational::zero(),
                            unit: EvidenceUnit::Millimetres,
                        }),
                        proof: GeometricProof::ExactDefinition(key(
                            &revision, input.0, "target", 0,
                        )),
                        dependencies: BTreeSet::new(),
                        context: ProofContext {
                            producer: "features.checked-circular-sweep-prism/1",
                            step: Some(key(&revision, input.0, "step", 0)),
                            stage: Some("source"),
                        },
                    },
                    dec,
                )?;
                let source = Arc::new(SourceProof {
                    geometry,
                    evidence,
                    body: input.1,
                    radius: input.2.clone(),
                });
                self.sources.insert(input.0.clone(), source.clone());
                source
            };
            let origins = source.geometry.origins.clone();
            (
                source.clone(),
                Rational::zero(),
                source.geometry.axial.clone(),
                origins,
            )
        };
        let magnitude = -signed_distance;
        let distance = previous_distance + magnitude.clone();
        let axial = [
            previous_axial[0].clone() + Real::from_rational(magnitude.clone()),
            previous_axial[1].clone() - Real::from_rational(magnitude),
        ];
        let Surface::EllipticCylinder(original) = &source.geometry.original else {
            unreachable!("checked elliptic source");
        };
        let definition =
            EllipticOffsetDefinition::new(original.clone(), Real::from_rational(-&distance), dec)?;
        let decision_key = key(&revision, step, "branch-decision", 0);
        let decision_evidence = key(&revision, step, "branch-proof", 0);
        let identity_branch = key(&revision, step, "identity", 0);
        let empty_branch = key(&revision, step, "empty", 0);
        let regular_branch = key(&revision, step, "regular-inward", 0);
        let inherited_empty = parent.as_ref().is_some_and(|p| p.target.empty);
        let choice = self.decisions.decide_once(decision_key, dec, |dec| {
            let branch = if inherited_empty || sign(&(&axial[1] - &axial[0]), dec)? != Sign::Positive
                || sign(&(&original.minor_sq.sqrt() - &Real::from_rational(distance.clone())), dec)? != Sign::Positive
            { empty_branch } else if signed_distance.is_zero() { identity_branch } else {
                let reach = &original.minor_sq / &original.major_sq.sqrt();
                if sign(&(reach - Real::from_rational(distance.clone())), dec)? != Sign::Positive {
                    return Err(Refusal::new("features.offset-nonregular-branch", "target is not empty; the inward normal branch needs a different contact/branch construction"));
                }
                regular_branch
            };
            Ok(RecordedChoice { branch, evidence: decision_evidence })
        })?.clone();
        let target = PrismTarget {
            source,
            definition,
            axial,
            empty: choice.branch == empty_branch,
        };
        let (body, representation, fit_target) = if target.empty {
            (
                model.add_body(Body {
                    kind: BodyKind::Empty,
                    lumps: Vec::new(),
                    marks: MarkLog::new(),
                }),
                Representation::Empty,
                None,
            )
        } else if choice.branch == identity_branch {
            let representation = parent.as_ref().map_or(Representation::Analytic, |parent| {
                parent.representation.clone()
            });
            (
                input.1,
                representation,
                parent.as_ref().and_then(|p| p.fit_target.clone()),
            )
        } else {
            let (body, representation, fit_target) = self.fitted(model, &target, &origins, dec)?;
            (body, representation, Some(fit_target))
        };
        let parent_evidence = parent.as_ref().map(|p| p.record.evidence);
        self.record_body(
            model,
            step,
            PreparedBody {
                target,
                representation,
                body,
                origins,
                choice: (decision_key, choice),
                fit_target,
                parent: parent_evidence,
            },
            dec,
        )
    }

    fn fitted(
        &self,
        model: &mut Model,
        target: &PrismTarget,
        origins: &[Vec<Origin>; 3],
        dec: &mut Decider,
    ) -> Result<(BodyId, Representation, Rational), Refusal> {
        let Surface::EllipticCylinder(original) = target.definition.original() else {
            unreachable!("checked elliptic target");
        };
        let inradius = positive_lower(
            &(&original.minor_sq.sqrt() - target.definition.distance()),
            dec,
        )?;
        let height = positive_lower(&(&target.axial[1] - &target.axial[0]), dec)?;
        let fit_target = self
            .evidence
            .local_settings()
            .kernel_target_mm()
            .clone()
            .min(&inradius / &Rational::from_i64(16))
            .min(&height / &Rational::from_i64(16));
        let EllipticOffsetFit { surface, .. } =
            fit_elliptic_offset(&target.definition, target.axial.clone(), &fit_target, dec)?;
        let proof =
            certify_elliptic_offset_sweep(&surface, &target.definition, target.axial(), dec)?;
        let spline = SplineSurface::from_offset(proof);
        let mut marks = MarkLog::new();
        marks.push(Mark::new(
            "features.certified-offset",
            BoundKind::Distance,
            spline.bound().upper.clone(),
            "complete checked elliptic-prism correspondence; the mark alone is not authority",
        ));
        let body = geometry::build_prism(model, &spline, origins, marks);
        Ok((body, Representation::Spline(spline), fit_target))
    }

    pub fn has_certificate(&self, step: &StepId) -> bool {
        self.certificates.contains_key(step)
    }

    pub fn boolean(
        &mut self,
        model: &mut Model,
        step: &StepId,
        kind: &BooleanKind,
        inputs: [(&StepId, BodyId); 2],
        dec: &mut Decider,
    ) -> Result<(BodyId, Rational), Refusal> {
        let parents = inputs.map(|(id, _)| self.certificates.get(id).cloned());
        for (parent, (_, body)) in parents.iter().zip(inputs) {
            if parent.as_ref().is_some_and(|parent| parent.body != body) {
                return Err(unproven("boolean input and retained certificate disagree"));
            }
        }
        let slot = parents
            .iter()
            .position(Option::is_some)
            .ok_or_else(|| unsupported("bounded composition has no certified input"))?;
        let parent = parents[slot].as_ref().expect("selected certified operand");
        let revision = self.retain()?;
        let decision_key = key(&revision, step, "branch-decision", 0);
        let decision_evidence = key(&revision, step, "branch-proof", 0);
        if inputs[0].0 == inputs[1].0 {
            let empty = *kind == BooleanKind::Subtract || parent.target.empty;
            let branch = key(
                &revision,
                step,
                if empty {
                    "empty-set-identity"
                } else {
                    "same-set-identity"
                },
                0,
            );
            let choice = self
                .decisions
                .decide_once(decision_key, dec, |_| {
                    Ok(RecordedChoice {
                        branch,
                        evidence: decision_evidence,
                    })
                })?
                .clone();
            let mut target = parent.target.clone();
            target.empty = empty;
            let (body, representation) = if empty {
                (
                    model.add_body(Body {
                        kind: BodyKind::Empty,
                        lumps: Vec::new(),
                        marks: MarkLog::new(),
                    }),
                    Representation::Empty,
                )
            } else {
                (parent.body, parent.representation.clone())
            };
            return self.record_body(
                model,
                step,
                PreparedBody {
                    target,
                    representation,
                    body,
                    origins: parent.origins.clone(),
                    choice: (decision_key, choice),
                    fit_target: parent.fit_target.clone(),
                    parent: Some(parent.record.evidence),
                },
                dec,
            );
        }
        if *kind != BooleanKind::Intersect || parents[1 - slot].is_some() {
            return Err(unsupported(
                "this bounded piece composes through exact axial-box intersections only",
            ));
        }
        let tool = inputs[1 - slot].0;
        let op = self
            .recipe
            .steps
            .iter()
            .find(|node| &node.id == tool)
            .map(|node| &node.op)
            .ok_or_else(|| unproven("retained trim definition is missing"))?;
        let Surface::EllipticCylinder(original) = parent.target.definition.original() else {
            unreachable!("checked elliptic target");
        };
        let (clipping, cap_origins) =
            source::axial_box(original, parent.target.axis(), op, tool, dec)?;
        let branches: [CertificateKey; 8] =
            std::array::from_fn(|mask| key(&revision, step, "axial-trim-selection", mask as u32));
        let choice = self
            .decisions
            .decide_once(decision_key, dec, |dec| {
                let mut mask = 0;
                if sign(&(&clipping[0] - &parent.target.axial[0]), dec)? == Sign::Positive {
                    mask |= 1;
                }
                if sign(&(&clipping[1] - &parent.target.axial[1]), dec)? == Sign::Negative {
                    mask |= 2;
                }
                let lo = if mask & 1 != 0 {
                    &clipping[0]
                } else {
                    &parent.target.axial[0]
                };
                let hi = if mask & 2 != 0 {
                    &clipping[1]
                } else {
                    &parent.target.axial[1]
                };
                if parent.target.empty || sign(&(hi - lo), dec)? != Sign::Positive {
                    mask |= 4;
                }
                Ok(RecordedChoice {
                    branch: branches[mask],
                    evidence: decision_evidence,
                })
            })?
            .clone();
        let mask = branches
            .iter()
            .position(|branch| branch == &choice.branch)
            .ok_or_else(|| unproven("axial trim reused an unrelated settled selection"))?;
        let mut target = parent.target.clone();
        let mut origins = parent.origins.clone();
        if mask & 1 != 0 {
            target.axial[0] = clipping[0].clone();
        }
        if mask & 2 != 0 {
            target.axial[1] = clipping[1].clone();
        }
        if mask & 3 != 0 {
            let mut changed = (*origins).clone();
            for index in 0..2 {
                if mask & (1 << index) != 0 {
                    changed[index] = cap_origins[index].clone();
                }
            }
            origins = Arc::new(changed);
        }
        target.empty |= mask & 4 != 0;
        let (body, representation, fit_target) = if target.empty {
            (
                model.add_body(Body {
                    kind: BodyKind::Empty,
                    lumps: Vec::new(),
                    marks: MarkLog::new(),
                }),
                Representation::Empty,
                None,
            )
        } else if mask == 0 {
            (
                parent.body,
                parent.representation.clone(),
                parent.fit_target.clone(),
            )
        } else if matches!(parent.representation, Representation::Analytic) {
            let body = catallax_boolean::boolean(
                model,
                inputs[0].1,
                inputs[1].1,
                catallax_boolean::BoolOp::Intersect,
                dec,
            )
            .map_err(|error| error.refusal)?
            .body;
            (body, Representation::Analytic, None)
        } else {
            let (body, representation, fit_target) = self.fitted(model, &target, &origins, dec)?;
            (body, representation, Some(fit_target))
        };
        self.record_body(
            model,
            step,
            PreparedBody {
                target,
                representation,
                body,
                origins,
                choice: (decision_key, choice),
                fit_target,
                parent: Some(parent.record.evidence),
            },
            dec,
        )
    }

    fn record_body(
        &mut self,
        model: &Model,
        step: &StepId,
        prepared: PreparedBody,
        dec: &mut Decider,
    ) -> Result<(BodyId, Rational), Refusal> {
        let PreparedBody {
            target,
            representation,
            body,
            origins,
            choice,
            fit_target,
            parent,
        } = prepared;
        let revision = self.revision.expect("retained certified revision");
        let target_key = key(&revision, step, "target", 0);
        let context = ProofContext {
            producer: "features.checked-elliptic-prism/1",
            step: Some(key(&revision, step, "step", 0)),
            stage: Some("complete-correspondence"),
        };
        let mut dependencies = BTreeSet::from([target.source.evidence, choice.1.evidence]);
        if let Some(parent) = parent {
            dependencies.insert(parent);
        }
        self.evidence.record(
            choice.1.evidence,
            GeometricEvidence {
                identity: CertificateIdentity {
                    revision,
                    representation: choice.0,
                    target: target_key,
                    scope: CertificateScope::Query(choice.0),
                },
                bound: None,
                proof: GeometricProof::ExactDefinition(choice.1.branch),
                dependencies: BTreeSet::from([target.source.evidence]),
                context: context.clone(),
            },
            dec,
        )?;
        let Surface::EllipticCylinder(original) = target.definition.original() else {
            unreachable!("checked elliptic target");
        };
        let (topology, bound) = match &representation {
            Representation::Empty => (None, None),
            Representation::Analytic => (
                Some(geometry::check_prism(
                    model,
                    body,
                    original,
                    target.axial(),
                    ProfileSupport::Analytic,
                    &origins,
                    dec,
                )?),
                Some(EvidenceBound {
                    upper: Rational::zero(),
                    unit: EvidenceUnit::Millimetres,
                }),
            ),
            Representation::Spline(surface) => (
                Some(geometry::check_prism(
                    model,
                    body,
                    original,
                    target.axial(),
                    ProfileSupport::Spline(surface),
                    &origins,
                    dec,
                )?),
                Some(surface.bound().clone()),
            ),
        };
        let mut regions = BTreeMap::new();
        if let (Some(topology), Some(bound)) = (&topology, &bound) {
            let (edge, face) = size_bounds(&target, bound, dec)?;
            for (index, actual) in topology.faces.iter().enumerate() {
                let region = key(&revision, step, "face", actual.index());
                let evidence = key(&revision, step, "face-correspondence", index as u32);
                let identity = CertificateIdentity {
                    revision,
                    representation: region,
                    target: key(&revision, step, "target-face", index as u32),
                    scope: CertificateScope::Entities(BTreeSet::from([region])),
                };
                self.evidence.record(
                    evidence,
                    GeometricEvidence {
                        identity: identity.clone(),
                        bound: Some(bound.clone()),
                        proof: if bound.upper.is_zero() {
                            GeometricProof::ExactDefinition(identity.target)
                        } else {
                            GeometricProof::Approximation {
                                definition: identity.target,
                                theorem: key(
                                    &revision,
                                    step,
                                    if index == 2 {
                                        "periodic-convex-sweep"
                                    } else {
                                        "complete-cap-fan"
                                    },
                                    index as u32,
                                ),
                            }
                        },
                        dependencies: BTreeSet::from([target.source.evidence, choice.1.evidence]),
                        context: context.clone(),
                    },
                    dec,
                )?;
                dependencies.insert(evidence);
                let exception = (!bound.upper.is_zero()).then(|| LocalException {
                    lower_mm: Rational::zero(),
                    bound: bound.clone(),
                    limit: EvidenceBound {
                        upper: self.evidence.local_settings().kernel_cap_mm().clone(),
                        unit: EvidenceUnit::Millimetres,
                    },
                    source: ExceptionSource::KernelApproximation,
                    shortest_edge: FeatureSizeLowerBound::Millimetres(edge.clone()),
                    narrowest_face: FeatureSizeLowerBound::Millimetres(face.clone()),
                });
                regions.insert(
                    region,
                    RegionCertificate {
                        identity,
                        evidence,
                        exception,
                    },
                );
            }
        }
        let identity = CertificateIdentity {
            revision,
            representation: key(&revision, step, "body", body.index()),
            target: target_key,
            scope: CertificateScope::WholeSolid,
        };
        let root = key(&revision, step, "whole-body-correspondence", 0);
        self.evidence.record(
            root,
            GeometricEvidence {
                identity: identity.clone(),
                bound: bound.clone(),
                proof: GeometricProof::WholeResult(key(
                    &revision,
                    step,
                    "convex-prism-homeomorphism",
                    0,
                )),
                dependencies,
                context,
            },
            dec,
        )?;
        let record = GeometricCertificate {
            identity,
            evidence: root,
            regions,
            decisions: BTreeSet::new(),
            uniform_bound: bound.clone(),
        };
        let certificate = BodyCertificate {
            body,
            step: step.clone(),
            target,
            representation,
            topology,
            origins,
            record,
            choice,
            fit_target,
        };
        certificate.verify_record(
            model,
            &certificate.record,
            &self.evidence,
            &self.decisions,
            dec,
        )?;
        let radius = match bound {
            None => Rational::zero(),
            Some(bound) => &certificate.target.source.radius + &bound.upper,
        };
        self.certificates
            .insert(step.clone(), Arc::new(certificate));
        Ok((body, radius))
    }
}

fn positive_lower(value: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    let max = dec.policy().max_precision_bits;
    let mut bits = 64.min(max);
    if bits == 0 {
        return Err(Refusal::new(
            "numbers.precision-limit",
            "no geometric-size working precision",
        ));
    }
    loop {
        let enclosure = dec.enclose_budgeted_at(value, bits)?;
        if enclosure.lo() > Rational::zero() {
            return Ok(enclosure.lo());
        }
        if bits == max {
            return Err(unproven(
                "positive neighbouring-size lower bound was not separated",
            ));
        }
        bits = bits.saturating_mul(2).min(max);
    }
}
fn size_bounds(
    target: &PrismTarget,
    error: &EvidenceBound,
    dec: &mut Decider,
) -> Result<(Rational, Rational), Refusal> {
    let Surface::EllipticCylinder(original) = target.definition.original() else {
        unreachable!("checked elliptic target");
    };
    let radius = positive_lower(
        &(&original.minor_sq.sqrt()
            - target.definition.distance()
            - Real::from_rational(error.upper.clone())),
        dec,
    )?;
    let height = positive_lower(&(&target.axial[1] - &target.axial[0]), dec)?;
    // The uniform convex support-function bound retains a disk of this radius.
    // Its ring length is >=4r; all adjacent cap/side widths are >=min(2r,h).
    let edge = &radius * &Rational::from_i64(4);
    let face = (&radius * &Rational::from_i64(2)).min(height);
    if error.upper > &face / &Rational::from_i64(10) {
        return Err(Refusal::new(
            GEOMETRIC_GAP_EXCEEDED,
            "whole-profile fit exceeds its proved neighbouring-face tenth",
        ));
    }
    Ok((edge, face))
}
