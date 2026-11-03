//! Complete certificate of a body with fitted rolling-ball bands: a
//! plane/cylinder blend continued over countersunk mouths, or rings between
//! nonparallel cylinders. Each model face has one region. Analytic faces keep
//! their exact definition, and their boundaries with a carrier follow that
//! carrier's rows, so the only fitted surfaces are the named bands.

use super::{CertificateKey, Refusal, key, positive_lower, unproven};
use catallax_numbers::{
    BoundKind, CertificateIdentity, CertificateScope, Decider, DecisionBook, EvidenceBook,
    EvidenceBound, EvidenceUnit, ExceptionSource, FeatureSizeLowerBound, GEOMETRIC_GAP_EXCEEDED,
    GeometricCertificate, GeometricEvidence, GeometricProof, LocalException, ModellingTolerance,
    ProofContext, Rational, Real, RecordedChoice, RegionCertificate, Sign,
};
use catallax_recipe::StepId;
use catallax_shapes::{
    BSplineCurve, BoundaryRow, Curve, Point3, RollingEnvelope, RollingEnvelopeDefinition,
    RollingSaddle, RollingSaddleDefinition, SaddleLaw, SaddleRing, SplineBoundary, SplineSupport,
    Surface, Vector3,
};
use catallax_topology::{BodyId, EdgeId, FaceId, Model, envelope_rows, saddle_rows};
use std::collections::{BTreeMap, BTreeSet};

/// The code of the body mark every fitted junction carries.
pub const FIT_MARK: &str = "features.rolling-envelope-fit";

/// One fitted face with its exact retained definition.
#[derive(Clone, Debug)]
pub struct EnvelopeRegion {
    face: FaceId,
    definition: RollingEnvelopeDefinition,
    bound: EvidenceBound,
    edges: Vec<(EdgeId, BoundaryRow)>,
}
impl EnvelopeRegion {
    /// Actual fitted face.
    pub fn face(&self) -> FaceId {
        self.face
    }
    /// Exact rolling-ball definition the carrier follows.
    pub fn definition(&self) -> &RollingEnvelopeDefinition {
        &self.definition
    }
    /// Checked two-way carrier bound.
    pub fn bound(&self) -> &EvidenceBound {
        &self.bound
    }
    /// The face's four boundary edges and their chart rows.
    pub fn edges(&self) -> &[(EdgeId, BoundaryRow)] {
        &self.edges
    }
}

/// One fitted ring between two cylinders with its exact retained definition.
#[derive(Clone, Debug)]
pub struct SaddleRegion {
    face: FaceId,
    definition: RollingSaddleDefinition,
    bound: EvidenceBound,
    rails: [EdgeId; 2],
}
impl SaddleRegion {
    /// Actual fitted face.
    pub fn face(&self) -> FaceId {
        self.face
    }
    /// Exact rolling-ball definition the carrier follows.
    pub fn definition(&self) -> &RollingSaddleDefinition {
        &self.definition
    }
    /// Checked two-way carrier bound.
    pub fn bound(&self) -> &EvidenceBound {
        &self.bound
    }
    /// The wide and narrow rails.
    pub fn rails(&self) -> [EdgeId; 2] {
        self.rails
    }
}

/// A complete mixed-body record with its retained definitions. Public numeric
/// records can be inspected, but only a checked model can mint this type.
#[derive(Clone, Debug)]
pub struct JunctionCertificate {
    body: BodyId,
    base: BodyId,
    step: StepId,
    faces: Vec<FaceId>,
    envelopes: Vec<EnvelopeRegion>,
    saddles: Vec<SaddleRegion>,
    sizes: (Rational, Rational),
    trims: BTreeMap<FaceId, Rational>,
    record: GeometricCertificate<CertificateKey>,
    choice: (CertificateKey, RecordedChoice<CertificateKey>),
    fit_target: Rational,
}

fn envelope_of(model: &Model, face: FaceId) -> Option<&RollingEnvelope> {
    match model.geom().surface(model.face(face).surface) {
        Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => Some(envelope),
        _ => None,
    }
}

// Lower bounds on the shortest edge and narrowest face beside this envelope. A
// row is at least its chord. The plane neck between mouth and wall is the
// narrowest width of any neighbour, and the envelope spans at least that gap.
fn neighbouring_sizes(
    definition: &RollingEnvelopeDefinition,
    dec: &mut Decider,
) -> Result<(Rational, Rational), Refusal> {
    let negative = definition.corners(Sign::Negative);
    let positive = definition.corners(Sign::Positive);
    let radius = definition.radius();
    let two = Real::from_i64(2);
    let chord =
        |a: &catallax_shapes::Point3, c: &catallax_shapes::Point3| a.minus(c).norm_squared().sqrt();
    let lengths = [
        chord(&negative[0], &positive[0]),
        chord(&negative[1], &positive[1]),
        radius * &two.sqrt(),
    ];
    let mut edge = None::<Rational>;
    for length in &lengths {
        let lower = positive_lower(length, dec)?;
        edge = Some(edge.map_or(lower.clone(), |value| value.min(lower)));
    }
    let gap = &(definition.offset() - definition.mouth()) - definition.body();
    let face = positive_lower(&gap, dec)?;
    Ok((edge.expect("three rows"), face))
}

fn saddle_of(model: &Model, face: FaceId) -> Option<&RollingSaddle> {
    match model.geom().surface(model.face(face).surface) {
        Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => Some(saddle),
        _ => None,
    }
}

// The rail's other edges on its wall, and the rail's nearest approach to each.
fn wall_necks(
    model: &Model,
    face: FaceId,
    rail: EdgeId,
    law: &SaddleLaw,
    origin: &Point3,
    dec: &mut Decider,
) -> Result<Vec<Real>, Refusal> {
    const BITS: u32 = 96;
    let Curve::BSpline(curve) = model.geom().curve(model.edge(rail).curve) else {
        return Err(unproven("a saddle rail is not a spline"));
    };
    let controls: Vec<_> = curve.controls().cloned().collect();
    // The hull's nearest approach to a plane, signed by the side it lies on.
    let side_of = |at: &Point3, normal: &Vector3, dec: &mut Decider| -> Result<Rational, Refusal> {
        let unit = normal.scaled(&(&Real::one() / &normal.norm_squared().sqrt()));
        let mut nearest = None::<Rational>;
        for control in &controls {
            dec.check_deadline()?;
            let height = dec.enclose_budgeted_at(&control.minus(at).dot(&unit), BITS)?;
            let value = if height.lo() > Rational::zero() {
                height.lo()
            } else if height.hi() < Rational::zero() {
                height.hi()
            } else {
                return Err(unproven("a saddle rail is not clear of its wall's edge"));
            };
            if nearest.as_ref().is_some_and(|old| {
                (old.sign() == Sign::Positive) != (value.sign() == Sign::Positive)
            }) {
                return Err(unproven("a saddle rail is not clear of its wall's edge"));
            }
            nearest = Some(match nearest {
                Some(old) if old.abs() <= value.abs() => old,
                _ => value,
            });
        }
        nearest.ok_or_else(|| unproven("a saddle rail has no controls"))
    };
    let mut necks = Vec::new();
    for coedge in model.coedges_of_edge(rail) {
        let Some((_, wall)) = model.owner_of_coedge(coedge) else {
            continue;
        };
        if wall == face {
            continue;
        }
        for lp in &model.face(wall).loops {
            let catallax_topology::Loop::Coedges(coedges) = model.loop_(*lp) else {
                continue;
            };
            for other in coedges {
                let other = model.coedge(*other).edge;
                if other == rail {
                    continue;
                }
                let neck = match model.geom().curve(model.edge(other).curve) {
                    Curve::Circle(circle) => side_of(&circle.center, &circle.normal, dec)?.abs(),
                    Curve::Ellipse(ellipse) => {
                        side_of(&ellipse.center, &ellipse.normal, dec)?.abs()
                    }
                    // The far ring of a through bore lies across the plane that
                    // holds the wide axis and `across`.
                    Curve::BSpline(far) => {
                        let near = side_of(origin, &law.along, dec)?;
                        let start = far.control(0).minus(origin).dot(&law.along);
                        let mut clear = None::<Rational>;
                        for control in far.controls() {
                            let height = dec.enclose_budgeted_at(
                                &control.minus(origin).dot(&law.along),
                                BITS,
                            )?;
                            let value = height.lo().abs().min(height.hi().abs());
                            clear = Some(clear.map_or(value.clone(), |old| old.min(value)));
                        }
                        if super::sign(&start, dec)? == near.sign() {
                            return Err(unproven("two saddle rails share a side of their wall"));
                        }
                        near.abs() + clear.expect("a rail has controls")
                    }
                    _ => return Err(unproven("a saddle's wall has an edge of no known width")),
                };
                necks.push(Real::from_rational(neck));
            }
        }
    }
    Ok(necks)
}

// Lower bounds on the shortest edge and narrowest face beside a ring.
//
// The narrow rail winds once round its cylinder, so it is at least that
// circumference. The wide rail's shadow on `across` is `2 p R / P` wide and it
// goes there and back. The band is at least as wide as the narrow rail's
// distance from the wide cylinder. With unit contact rays A and B from the ball
// centre, `D = A.B` is at most `sqrt(1 - S^2 rho)` for `rho = 1 - x^2 / P^2` at
// the largest `x`, by Cauchy. Writing a narrow contact's distance from the
// wide axis as `w`, `|w^2 - R^2| >= r min(R, P) l` for `l = 2 (1 - D)`, and
// `w + R <= 2 (R + r)`. A carrier's rails sit within its bound of the law's.
fn saddle_sizes(
    model: &Model,
    face: FaceId,
    saddle: &RollingSaddle,
    dec: &mut Decider,
) -> Result<(Rational, Rational), Refusal> {
    let definition = saddle.definition();
    let law = definition.law();
    // These bounds rely on a ring that winds round the narrow cylinder.
    if law.ring != SaddleRing::Through {
        return Err(unproven("a two-cylinder loop has no size bounds yet"));
    }
    let rails = saddle_rows(model, face).ok_or_else(|| unproven("a ring lost a rail"))?;
    let big = &law.wide + &(&law.radius * &Real::from_i64(law.sides[0]));
    let small = &law.narrow + &(&law.radius * &Real::from_i64(law.sides[1]));
    let reach = &(law.gap.square().sqrt() + small.clone()) / &big;
    let open = &(&Real::one() - &law.sine.square()) + &(&law.sine.square() * &reach.square());
    let spread = (&Real::one() - &open.sqrt()) * Real::from_i64(2);
    let near = if law.sides[0] > 0 { &law.wide } else { &big };
    let lift = &(&law.radius * near) * &spread;
    let width = &lift / &((&law.wide + &law.radius) * Real::from_i64(2));
    let slack = Real::from_rational(saddle.bound().upper.clone()) * Real::from_i64(2);
    let turn = &law.narrow * &Real::from_i64(6);
    let shadow = &(&small * &law.wide) * &(&Real::from_i64(4) / &big);
    let edge = positive_lower(&turn, dec)?.min(positive_lower(&shadow, dec)?);
    let mut narrowest = positive_lower(&(&width - &slack), dec)?;
    for (rail, _) in rails {
        for neck in wall_necks(model, face, rail, &law, definition.origin(), dec)? {
            narrowest = narrowest.min(positive_lower(&neck, dec)?);
        }
    }
    Ok((edge, narrowest))
}

fn same_real(a: &Real, c: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(super::sign(&(a - c), dec)? == Sign::Zero)
}

fn same_vector(a: &Vector3, c: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(same_real(&a.x, &c.x, dec)? && same_real(&a.y, &c.y, dec)? && same_real(&a.z, &c.z, dec)?)
}

fn same_law(
    a: &RollingSaddleDefinition,
    c: &RollingSaddleDefinition,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let (first, second) = (a.law(), c.law());
    if first.sides != second.sides || first.branch != second.branch {
        return Ok(false);
    }
    for (x, y) in [
        (&first.wide, &second.wide),
        (&first.narrow, &second.narrow),
        (&first.radius, &second.radius),
        (&first.gap, &second.gap),
        (&first.cosine, &second.cosine),
        (&first.sine, &second.sine),
    ] {
        if !same_real(x, y, dec)? {
            return Ok(false);
        }
    }
    for (x, y) in [
        (a.origin().minus(c.origin()), Vector3::zero()),
        (first.foot.minus(&second.foot), Vector3::zero()),
        (first.across, second.across),
        (first.along, second.along),
        (first.axis, second.axis),
    ] {
        if !same_vector(&x, &y, dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}

// Each rail is the carrier's own row and stays on its exact cylinder, so a
// boundary cannot move without both its faces moving with it.
fn check_rails(
    model: &Model,
    region: &SaddleRegion,
    saddle: &RollingSaddle,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let rows =
        saddle_rows(model, region.face).ok_or_else(|| unproven("a ring lost a retained rail"))?;
    let law = saddle.definition().law();
    let narrow_axis = law.narrow_axis();
    let axes = [
        (
            &law.wide,
            saddle.definition().origin().clone(),
            law.axis.clone(),
        ),
        (&law.narrow, law.foot.clone(), narrow_axis),
    ];
    for (index, ((rail, _), (radius, at, axis))) in rows.iter().zip(axes).enumerate() {
        if *rail != region.rails[index] {
            return Err(unproven("a ring lost a retained rail"));
        }
        let mut neighbours = model
            .coedges_of_edge(*rail)
            .into_iter()
            .filter_map(|coedge| model.owner_of_coedge(coedge).map(|(_, face)| face))
            .filter(|face| *face != region.face);
        let wall = neighbours
            .next()
            .ok_or_else(|| unproven("a rail has no neighbouring face"))?;
        let Surface::Cylinder(cylinder) = model.geom().surface(model.face(wall).surface) else {
            return Err(unproven("a rail left its exact neighbouring surface"));
        };
        let offset = at.minus(&cylinder.origin).cross(&cylinder.axis);
        if !same_real(&cylinder.radius_sq, &radius.square(), dec)?
            || !same_vector(&axis.cross(&cylinder.axis), &Vector3::zero(), dec)?
            || !same_vector(&offset, &Vector3::zero(), dec)?
        {
            return Err(unproven("a rail left its exact neighbouring surface"));
        }
    }
    Ok(())
}

impl JunctionCertificate {
    /// Actual immutable represented body.
    pub fn body(&self) -> BodyId {
        self.body
    }
    /// Unfilleted body the blend continued from.
    pub fn base(&self) -> BodyId {
        self.base
    }
    /// Producing recipe node.
    pub fn step(&self) -> &StepId {
        &self.step
    }
    /// Complete face/trim and whole-result evidence claim.
    pub fn record(&self) -> &GeometricCertificate<CertificateKey> {
        &self.record
    }
    /// Every represented face, each with one region.
    pub fn faces(&self) -> &[FaceId] {
        &self.faces
    }
    /// The named fitted faces.
    pub fn envelopes(&self) -> &[EnvelopeRegion] {
        &self.envelopes
    }
    /// The named fitted rings.
    pub fn saddles(&self) -> &[SaddleRegion] {
        &self.saddles
    }
    /// Numeric region bound to this actual face.
    pub fn region_for_face(&self, face: FaceId) -> Option<&RegionCertificate<CertificateKey>> {
        if !self.faces.contains(&face) {
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
    /// A face's own bound to its definition: a carrier's fit, or the fit of any
    /// carrier row that trims an exact surface, or zero.
    pub fn face_bound(&self, face: FaceId) -> Rational {
        self.trims
            .get(&face)
            .cloned()
            .unwrap_or_else(Rational::zero)
    }
    /// Producer-selected fitting target, the local tenth or the kernel target.
    pub fn fit_target_mm(&self) -> &Rational {
        &self.fit_target
    }

    // Checks one represented body against the retained definitions, independent
    // of any numeric claim.
    pub(super) fn check_body(
        &self,
        model: &Model,
        body: BodyId,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if !model.check(body).is_empty() {
            return Err(unproven("the represented body fails its topology check"));
        }
        let faces = model.faces_of_body(body);
        if faces != self.faces || model.blend_base(body) != Some(self.base) {
            return Err(unproven(
                "junction proof is detached from its complete actual topology",
            ));
        }
        for region in &self.envelopes {
            let envelope = envelope_of(model, region.face)
                .ok_or_else(|| unproven("a named envelope region lost its carrier"))?;
            let definition = envelope.definition();
            let retained = &region.definition;
            for (a, c) in [
                (definition.radius(), retained.radius()),
                (definition.body(), retained.body()),
                (definition.mouth(), retained.mouth()),
                (definition.offset(), retained.offset()),
                (&definition.origin().x, &retained.origin().x),
                (&definition.origin().y, &retained.origin().y),
                (&definition.origin().z, &retained.origin().z),
                (&definition.toward().x, &retained.toward().x),
                (&definition.toward().y, &retained.toward().y),
                (&definition.toward().z, &retained.toward().z),
                (&definition.axis().x, &retained.axis().x),
                (&definition.axis().y, &retained.axis().y),
                (&definition.axis().z, &retained.axis().z),
            ] {
                if super::sign(&(a - c), dec)? != Sign::Zero {
                    return Err(unproven(
                        "a carrier follows a different retained definition",
                    ));
                }
            }
            if envelope.bound().upper > region.bound.upper {
                return Err(unproven("a carrier's checked bound grew past its record"));
            }
            check_rows(model, region, envelope, dec)?;
        }
        for region in &self.saddles {
            let saddle = saddle_of(model, region.face)
                .ok_or_else(|| unproven("a named ring lost its carrier"))?;
            if !same_law(saddle.definition(), &region.definition, dec)? {
                return Err(unproven(
                    "a carrier follows a different retained definition",
                ));
            }
            if saddle.bound().upper > region.bound.upper {
                return Err(unproven("a carrier's checked bound grew past its record"));
            }
            check_rails(model, region, saddle, dec)?;
        }
        Ok(())
    }

    /// Rechecks an untrusted numeric claim against actual body geometry and the
    /// retained definitions.
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
        self.check_body(model, self.body, dec)?;
        let mut worst = Rational::zero();
        let fitted = self
            .envelopes
            .iter()
            .map(|region| (region.face, &region.bound))
            .chain(
                self.saddles
                    .iter()
                    .map(|region| (region.face, &region.bound)),
            );
        for (face, upper) in fitted {
            worst = worst.max(upper.upper.clone());
            let id = key(
                &self.record.identity.revision,
                &self.step,
                "face",
                face.index(),
            );
            let exception = candidate
                .regions
                .get(&id)
                .and_then(|actual| actual.exception.as_ref())
                .ok_or_else(|| unproven("a fitted face is omitted from the exception cover"))?;
            let (edge, face) = &self.sizes;
            if exception.bound.upper < upper.upper
                || !matches!(exception.source, ExceptionSource::KernelApproximation)
                || !matches!(&exception.shortest_edge, FeatureSizeLowerBound::Millimetres(value) if value <= edge)
                || !matches!(&exception.narrowest_face, FeatureSizeLowerBound::Millimetres(value) if value <= face)
            {
                return Err(unproven(
                    "local metric, fitting provenance or actual neighbouring-size proof is false",
                ));
            }
        }
        for (id, region) in &candidate.regions {
            if region.identity != self.record.regions[id].identity {
                return Err(unproven("a face or complete trim domain was rebound"));
            }
        }
        let uniform = candidate
            .uniform_bound
            .as_ref()
            .ok_or_else(|| unproven("junction has no whole-result correspondence bound"))?;
        if uniform.unit != EvidenceUnit::Millimetres || uniform.upper < worst {
            return Err(unproven(
                "uniform bound understates the independent whole-body theorem",
            ));
        }
        let marks = &model.body(self.body).marks;
        if !marks.iter().any(|mark| {
            mark.code == FIT_MARK && mark.kind == BoundKind::Distance && mark.bound == worst
        }) {
            return Err(unproven(
                "body mark is missing, forged or detached from its full certificate",
            ));
        }
        Ok(())
    }
}

// Each row stays on the carrier and on its exact neighbour, so a boundary cannot
// move without both its faces moving with it.
fn check_rows(
    model: &Model,
    region: &EnvelopeRegion,
    envelope: &RollingEnvelope,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let definition = envelope.definition();
    let actual = envelope_rows(model, region.face, dec)?;
    if actual.len() != region.edges.len()
        || actual
            .iter()
            .zip(&region.edges)
            .any(|((edge, row, _), expected)| (*edge, *row) != *expected)
    {
        return Err(unproven("an envelope lost a retained boundary row"));
    }
    for (edge, row, _) in &actual {
        let shared = model.edge(*edge);
        let neighbour = model
            .coedges_of_edge(*edge)
            .into_iter()
            .filter_map(|coedge| model.owner_of_coedge(coedge).map(|(_, face)| face))
            .find(|face| *face != region.face)
            .ok_or_else(|| unproven("a boundary row has no neighbouring face"))?;
        let surface = model.geom().surface(model.face(neighbour).surface);
        let exact = match (row, surface) {
            (BoundaryRow::Rail, Surface::Cylinder(cylinder)) => {
                let curve = model.geom().curve(shared.curve);
                let rail = BSplineCurve::boundary(envelope.carrier(), SplineBoundary::VMin);
                matches!(curve, Curve::BSpline(actual) if actual.same_definition(&rail))
                    && super::sign(&(&cylinder.radius_sq - &definition.body().square()), dec)?
                        == Sign::Zero
            }
            (BoundaryRow::Mouth, Surface::Cone(_)) => {
                matches!(model.geom().curve(shared.curve), Curve::Circle(_))
            }
            (BoundaryRow::Negative | BoundaryRow::Positive, Surface::Torus(torus)) => {
                matches!(model.geom().curve(shared.curve), Curve::Circle(_))
                    && super::sign(&(&torus.minor - definition.radius()), dec)? == Sign::Zero
            }
            _ => false,
        };
        if !exact {
            return Err(unproven(
                "a boundary row left its carrier or exact neighbouring surface",
            ));
        }
        for vertex in [shared.start, shared.end].into_iter().flatten() {
            let point = model.geom().position(model.vertex(vertex).point).clone();
            if !envelope.on_boundary(&point, dec)? {
                return Err(unproven("a row end left the carrier's boundary"));
            }
        }
    }
    Ok(())
}

/// Certifies an already built result, or returns `None` when no face is fitted.
pub(super) fn certify(
    model: &Model,
    body: BodyId,
    step: &StepId,
    revision: CertificateKey,
    evidence: &mut EvidenceBook<CertificateKey>,
    decisions: &mut DecisionBook<CertificateKey>,
    dec: &mut Decider,
) -> Result<Option<JunctionCertificate>, Refusal> {
    let faces = model.faces_of_body(body);
    if !faces
        .iter()
        .any(|face| envelope_of(model, *face).is_some() || saddle_of(model, *face).is_some())
    {
        return Ok(None);
    }
    let base = model
        .blend_base(body)
        .ok_or_else(|| unproven("a fitted blend has no recorded unfilleted base"))?;
    let target_key = key(&revision, step, "target", 0);
    let decision_key = key(&revision, step, "branch-decision", 0);
    let decision_evidence = key(&revision, step, "branch-proof", 0);
    let branch = key(&revision, step, "rolling-envelope-continuation", 0);
    let choice = decisions
        .decide_once(decision_key, dec, |_| {
            Ok(RecordedChoice {
                branch,
                evidence: decision_evidence,
            })
        })?
        .clone();
    let context = ProofContext {
        producer: "features.checked-rolling-junction/1",
        step: Some(key(&revision, step, "step", 0)),
        stage: Some("complete-correspondence"),
    };
    let source = key(&revision, step, "exact-source", base.index());
    evidence.record(
        source,
        GeometricEvidence {
            identity: CertificateIdentity {
                revision,
                representation: key(&revision, step, "source-body", base.index()),
                target: target_key,
                scope: CertificateScope::WholeSolid,
            },
            bound: Some(EvidenceBound {
                upper: Rational::zero(),
                unit: EvidenceUnit::Millimetres,
            }),
            proof: GeometricProof::ExactDefinition(target_key),
            dependencies: BTreeSet::new(),
            context: ProofContext {
                stage: Some("source"),
                ..context.clone()
            },
        },
        dec,
    )?;
    evidence.record(
        choice.evidence,
        GeometricEvidence {
            identity: CertificateIdentity {
                revision,
                representation: decision_key,
                target: target_key,
                scope: CertificateScope::Query(decision_key),
            },
            bound: None,
            proof: GeometricProof::ExactDefinition(choice.branch),
            dependencies: BTreeSet::from([source]),
            context: context.clone(),
        },
        dec,
    )?;
    let mut envelopes = Vec::new();
    let mut saddles = Vec::new();
    let mut trims = BTreeMap::<FaceId, Rational>::new();
    let mut worst = Rational::zero();
    let mut sizes: Option<(Rational, Rational)> = None;
    for &face in &faces {
        if let Some(envelope) = envelope_of(model, face) {
            let bound = envelope.bound().clone();
            let (edge, width) = neighbouring_sizes(envelope.definition(), dec)?;
            sizes = Some(match sizes {
                Some((a, c)) => (a.min(edge), c.min(width)),
                None => (edge, width),
            });
            worst = worst.max(bound.upper.clone());
            let edges: Vec<_> = envelope_rows(model, face, dec)?
                .into_iter()
                .map(|(edge, row, _)| (edge, row))
                .collect();
            for coedge in edges
                .iter()
                .flat_map(|(edge, _)| model.coedges_of_edge(*edge))
            {
                if let Some((_, owner)) = model.owner_of_coedge(coedge) {
                    let trim = trims.entry(owner).or_insert_with(Rational::zero);
                    *trim = trim.clone().max(bound.upper.clone());
                }
            }
            envelopes.push(EnvelopeRegion {
                face,
                definition: envelope.definition().clone(),
                bound,
                edges,
            });
        } else if let Some(saddle) = saddle_of(model, face) {
            let bound = saddle.bound().clone();
            let (edge, width) = saddle_sizes(model, face, saddle, dec)?;
            sizes = Some(match sizes {
                Some((a, c)) => (a.min(edge), c.min(width)),
                None => (edge, width),
            });
            worst = worst.max(bound.upper.clone());
            let rails = saddle_rows(model, face).ok_or_else(|| unproven("a ring lost a rail"))?;
            for coedge in rails
                .iter()
                .flat_map(|(edge, _)| model.coedges_of_edge(*edge))
            {
                if let Some((_, owner)) = model.owner_of_coedge(coedge) {
                    let trim = trims.entry(owner).or_insert_with(Rational::zero);
                    *trim = trim.clone().max(bound.upper.clone());
                }
            }
            saddles.push(SaddleRegion {
                face,
                definition: saddle.definition().clone(),
                bound,
                rails: rails.map(|(edge, _)| edge),
            });
        }
    }
    let sizes = sizes.expect("at least one fitted band");
    let settings = evidence.local_settings();
    let fit_target = settings
        .kernel_target_mm()
        .clone()
        .min(&sizes.0 / &Rational::from_i64(10))
        .min(&sizes.1 / &Rational::from_i64(10));
    let kernel_cap = settings.kernel_cap_mm().clone();
    if worst > fit_target {
        return Err(Refusal::new(
            GEOMETRIC_GAP_EXCEEDED,
            "a carrier exceeds its local target or neighbouring tenth",
        ));
    }
    let mut regions = BTreeMap::new();
    let mut dependencies = BTreeSet::from([source, choice.evidence]);
    for (index, &face) in faces.iter().enumerate() {
        let region = key(&revision, step, "face", face.index());
        let record = key(&revision, step, "face-correspondence", index as u32);
        let identity = CertificateIdentity {
            revision,
            representation: region,
            target: key(&revision, step, "target-face", index as u32),
            scope: CertificateScope::Entities(BTreeSet::from([region])),
        };
        let fitted = envelopes
            .iter()
            .find(|region| region.face == face)
            .map(|region| &region.bound)
            .or(saddles
                .iter()
                .find(|region| region.face == face)
                .map(|region| &region.bound));
        let bound = fitted.cloned().unwrap_or(EvidenceBound {
            upper: Rational::zero(),
            unit: EvidenceUnit::Millimetres,
        });
        evidence.record(
            record,
            GeometricEvidence {
                identity: identity.clone(),
                bound: Some(bound.clone()),
                proof: if fitted.is_some() {
                    GeometricProof::Approximation {
                        definition: identity.target,
                        theorem: key(
                            &revision,
                            step,
                            "checked-rolling-ball-envelope",
                            index as u32,
                        ),
                    }
                } else {
                    GeometricProof::ExactDefinition(identity.target)
                },
                dependencies: BTreeSet::from([source, choice.evidence]),
                context: context.clone(),
            },
            dec,
        )?;
        dependencies.insert(record);
        let exception = fitted.map(|_| LocalException {
            lower_mm: Rational::zero(),
            bound,
            limit: EvidenceBound {
                upper: kernel_cap.clone(),
                unit: EvidenceUnit::Millimetres,
            },
            source: ExceptionSource::KernelApproximation,
            shortest_edge: FeatureSizeLowerBound::Millimetres(sizes.0.clone()),
            narrowest_face: FeatureSizeLowerBound::Millimetres(sizes.1.clone()),
        });
        regions.insert(
            region,
            RegionCertificate {
                identity,
                evidence: record,
                exception,
            },
        );
    }
    let identity = CertificateIdentity {
        revision,
        representation: key(&revision, step, "body", body.index()),
        target: target_key,
        scope: CertificateScope::WholeSolid,
    };
    let uniform = EvidenceBound {
        upper: worst,
        unit: EvidenceUnit::Millimetres,
    };
    let root = key(&revision, step, "whole-body-correspondence", 0);
    evidence.record(
        root,
        GeometricEvidence {
            identity: identity.clone(),
            bound: Some(uniform.clone()),
            proof: GeometricProof::WholeResult(key(
                &revision,
                step,
                "rolling-junction-homeomorphism",
                0,
            )),
            dependencies,
            context,
        },
        dec,
    )?;
    let certificate = JunctionCertificate {
        body,
        base,
        step: step.clone(),
        faces,
        envelopes,
        saddles,
        sizes,
        trims,
        record: GeometricCertificate {
            identity,
            evidence: root,
            regions,
            decisions: BTreeSet::new(),
            uniform_bound: Some(uniform),
        },
        choice: (decision_key, choice),
        fit_target,
    };
    certificate.verify_record(model, &certificate.record, evidence, decisions, dec)?;
    Ok(Some(certificate))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CertifiedBuilt, evaluate_certified};
    use catallax_numbers::{BudgetPolicy, GEOMETRIC_GAP_UNPROVEN};
    use catallax_recipe::Recipe;
    use catallax_shapes::{Pose, Vector3};
    use catallax_topology::{Coedge, Edge, Face, Loop};
    use serde_json::{Value, json};

    fn decider() -> Decider {
        Decider::new(BudgetPolicy {
            max_work: 20_000_000,
            ..Default::default()
        })
    }
    fn built(radius: &str) -> CertifiedBuilt {
        let mut recipe: Value =
            serde_json::from_str(include_str!("../../../../parts/P5-flanged-housing.json"))
                .unwrap();
        recipe["steps"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|step| step["id"] == "filleted")
            .unwrap()["radius"] = json!(radius);
        evaluate_certified(&Recipe::parse(&recipe.to_string()).unwrap(), &mut decider()).unwrap()
    }
    fn rejected(result: Result<(), Refusal>) {
        assert_eq!(result.unwrap_err().code, GEOMETRIC_GAP_UNPROVEN);
    }

    // A rail's new edge in the envelope face alone, or in both its faces.
    fn moved_rail(certified: &CertifiedBuilt, both: bool) -> (Model, EnvelopeRegion) {
        let certificate = certified.junction().unwrap();
        let mut model = certified.built().model.clone();
        let region = &certificate.envelopes()[0];
        let (rail, _) = *region
            .edges()
            .iter()
            .find(|(_, row)| *row == BoundaryRow::Rail)
            .unwrap();
        let shared = model.edge(rail).clone();
        let moved = model
            .geom()
            .curve(shared.curve)
            .transformed(&Pose {
                translation: Vector3 {
                    x: Real::zero(),
                    y: Real::from_rational(Rational::from_ratio(1, 100).unwrap()),
                    z: Real::zero(),
                },
                ..Pose::identity()
            })
            .unwrap();
        let curve = model.geom_mut().add_curve(moved);
        let edge = model.add_edge(Edge { curve, ..shared });
        let replaced = |model: &mut Model, face: FaceId| -> FaceId {
            let old = model.face(face).clone();
            let mut loops = Vec::new();
            for lp in &old.loops {
                let Loop::Coedges(coedges) = model.loop_(*lp).clone() else {
                    loops.push(*lp);
                    continue;
                };
                let fresh: Vec<_> = coedges
                    .iter()
                    .map(|coedge| {
                        let use_ = model.coedge(*coedge).clone();
                        model.add_coedge(Coedge {
                            edge: if use_.edge == rail { edge } else { use_.edge },
                            ..use_
                        })
                    })
                    .collect();
                loops.push(model.add_loop(Loop::Coedges(fresh)));
            }
            model.add_face(Face { loops, ..old })
        };
        let neighbour = model
            .coedges_of_edge(rail)
            .into_iter()
            .filter_map(|coedge| model.owner_of_coedge(coedge).map(|(_, face)| face))
            .find(|face| *face != region.face)
            .unwrap();
        let face = replaced(&mut model, region.face);
        if both {
            replaced(&mut model, neighbour);
        }
        let mut forged = region.clone();
        forged.face = face;
        for entry in &mut forged.edges {
            if entry.0 == rail {
                entry.0 = edge;
            }
        }
        (model, forged)
    }

    fn ring(radius: &str) -> CertifiedBuilt {
        let cylinder = |id: &str, radius: &str, height: &str, rotation, at| {
            json!({"id":id,"op":"cylinder","radius":radius,"height":height,
                "pose":{"rotation":rotation,"translation":at}})
        };
        let side = json!({"side":{"segment":"wall"}});
        let recipe = json!({"format":"catallax-recipe/3","name":"ring","result":"round","steps":[
            cylinder("wide", "10", "60", ["1", "1", "0", "0"], ["0", "30", "0"]),
            cylinder("wall", "4", "30", ["1", "0", "0", "0"], ["1", "0", "0"]),
            {"id":"joined","op":"boolean","kind":"union","target":"wide","tool":"wall"},
            {"id":"round","op":"fillet","body":"joined","radius":radius,"edges":[
                {"between":[{"step":"wide","role":side},{"step":"wall","role":side}]}]}]});
        let recipe = Recipe::parse(&recipe.to_string()).unwrap();
        let built = crate::evaluate(&recipe, &mut decider()).unwrap();
        crate::certify_junction(built, &mut decider()).unwrap()
    }

    #[test]
    fn forged_ring_records_are_rejected() {
        let certified = ring("1");
        let certificate = certified.junction().unwrap();
        let model = &certified.built().model;
        let check = |certificate: &JunctionCertificate,
                     candidate: &GeometricCertificate<CertificateKey>,
                     model: &Model| {
            certificate.verify_record(
                model,
                candidate,
                certified.evidence(),
                certified.decisions(),
                &mut decider(),
            )
        };
        check(certificate, certificate.record(), model).unwrap();

        // The fitted face left out of the exception cover.
        let mut missing = certificate.record().clone();
        let band = certificate.saddles()[0].face();
        let id = certificate
            .region_for_face(band)
            .unwrap()
            .identity
            .representation;
        missing.regions.get_mut(&id).unwrap().exception = None;
        rejected(check(certificate, &missing, model));

        // A model of the same shape with another radius keeps every handle.
        let other = ring("5/4");
        assert_eq!(other.built().result, certified.built().result);
        rejected(check(
            certificate,
            certificate.record(),
            &other.built().model,
        ));

        // Another law, swapped rails, and a bound below the carrier's own.
        let mut forged = certificate.clone();
        forged.saddles[0].definition = other.junction().unwrap().saddles[0].definition.clone();
        rejected(check(&forged, certificate.record(), model));
        let mut forged = certificate.clone();
        forged.saddles[0].rails.swap(0, 1);
        rejected(check(&forged, certificate.record(), model));
        let mut forged = certificate.clone();
        forged.saddles[0].bound.upper = Rational::zero();
        rejected(check(&forged, certificate.record(), model));
    }

    #[test]
    fn forged_junction_records_are_rejected() {
        let certified = built("2");
        let certificate = certified.junction().unwrap();
        let model = &certified.built().model;
        let check = |candidate: &GeometricCertificate<CertificateKey>, model: &Model| {
            certificate.verify_record(
                model,
                candidate,
                certified.evidence(),
                certified.decisions(),
                &mut decider(),
            )
        };
        check(certificate.record(), model).unwrap();
        assert_eq!(certificate.faces().len(), 22);
        assert_eq!(certificate.envelopes().len(), 4);
        assert_eq!(
            certificate
                .record()
                .regions
                .values()
                .filter(|region| region.exception.is_some())
                .count(),
            4
        );

        let mut missing = certificate.record().clone();
        let fitted = certificate.region_for_face(certificate.envelopes()[0].face());
        let id = fitted.unwrap().identity.representation;
        missing.regions.remove(&id);
        rejected(check(&missing, model));

        let mut stale = certificate.record().clone();
        stale.identity.revision[0] ^= 1;
        rejected(check(&stale, model));

        // A model of the same shape with another radius keeps every handle.
        let other = built("5/2");
        assert_eq!(other.built().result, certified.built().result);
        rejected(check(certificate.record(), &other.built().model));

        for both in [false, true] {
            let (forged, region) = moved_rail(&certified, both);
            let envelope = envelope_of(&forged, region.face).unwrap();
            rejected(check_rows(&forged, &region, envelope, &mut decider()));
        }
    }
}
