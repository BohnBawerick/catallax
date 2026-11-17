//! ShapeFix-style source-topology repairs, published as one candidate transaction.
//! Geometry measurements establish each requested increase; the tolerance
//! ledger records the increases actually committed. It is not their proof.

use std::collections::BTreeMap;

use catallax_numbers::{
    Ball, Decider, Rational, Refusal, current_operation_work, operation_budget_active,
    with_operation_budget,
};
use serde::{Deserialize, Serialize};

use super::{Error, ToleranceEntity, ToleranceRecord, TolerantBody, TolerantFile, TolerantUse};
use crate::{Curve, Placement, StepFile, Surface, Vec3};

mod checks;
pub(crate) mod constructed;
mod native_sewing;
mod sewing;
mod wire;

/// A completed healing stage, with deterministic work charged to this request.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HealingStage {
    /// Source body name.
    pub body: String,
    /// Sewing, wire, face, shell or solid.
    pub stage: String,
    /// Work used by this stage, including unsuccessful candidate checks.
    pub work_used: u64,
}

/// A repair performed on a private candidate; never on the caller's source.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum HealingChange {
    /// Two finite free edges were matched by geometry and endpoint identity.
    EdgesSewn {
        /// Source body.
        body: String,
        /// Retained EDGE_CURVE.
        retained: String,
        /// Replaced EDGE_CURVE; original geometry remains in the source.
        replaced: String,
        /// Whether the replaced edge's uses change direction.
        reversed: bool,
        /// Source endpoint correspondence distances, in millimetres.
        endpoint_distance_upper_mm: [String; 2],
        /// Which finite-geometry checks established this correspondence.
        checks: Vec<String>,
    },
    /// A source endpoint now uses a shared retained vertex.
    VerticesJoined {
        /// Source body.
        body: String,
        /// Replaced vertex.
        source: String,
        /// Retained vertex.
        retained: String,
        /// Actual point displacement enclosure, in millimetres.
        distance_upper_mm: String,
    },
    /// A loop was reordered without removing a source edge use.
    WireReordered {
        /// Source body.
        body: String,
        /// Source face.
        face: String,
        /// Bound index in source order.
        boundary: usize,
        /// Original use indices in their new order.
        order: Vec<usize>,
        /// Original use indices whose traversal was reversed.
        reversed_use_indices: Vec<usize>,
    },
    /// An oriented boundary was reversed.
    BoundaryReversed {
        /// Source body.
        body: String,
        /// Source face.
        face: String,
        /// Bound index.
        boundary: usize,
    },
    /// A face and all its bounds were reversed together.
    FaceReversed {
        /// Source body.
        body: String,
        /// Source face.
        face: String,
        /// Shell consistency or whole-solid outward orientation.
        cause: String,
    },
    /// A finite source line segment below the initial uncertainty was removed.
    SmallLineRemoved {
        /// Source body.
        body: String,
        /// Removed source edge; it is not rebound to a different edge.
        edge: String,
        /// Retained endpoint representing the collapsed segment.
        retained_vertex: String,
        /// Measured source segment length enclosure, in millimetres.
        source_length_upper_mm: String,
    },
    /// Source surface shells joined by sampled finite-edge correspondences.
    ShellsSewn {
        /// Source body.
        body: String,
        /// Original shell names in this connected component.
        source_shells: Vec<String>,
        /// Shell name retained by the rebuilt candidate.
        retained_shell: String,
    },
    /// A face trace was reconstructed from its finite 3D curve and surface.
    TraceRebuilt {
        /// Source body.
        body: String,
        /// Source face.
        face: String,
        /// Source edge.
        edge: String,
        /// Replaced trace, or no original trace.
        source_pcurve: Option<u64>,
        /// Evaluated reconstruction stations.
        samples: usize,
        /// Largest measured 3D correspondence residual in millimetres.
        maximum_residual_mm: String,
    },
    /// Adjacent boundary curves were cut at their measured crossing.
    WireCrossingCut {
        /// Source body.
        body: String,
        /// Source face.
        face: String,
        /// Source boundary.
        boundary: usize,
        /// The two affected source edges.
        edges: [String; 2],
        /// New common point, in millimetres.
        point_mm: [String; 3],
        /// Natural curve parameters at the cut.
        parameters: [String; 2],
    },
    /// A collapsed surface boundary received a missing degenerate edge.
    DegenerateEdgeInserted {
        /// Source body.
        body: String,
        /// Source face.
        face: String,
        /// Source boundary.
        boundary: usize,
        /// New edge identity.
        edge: String,
        /// Existing physical pole vertex.
        vertex: String,
        /// Largest measured pole-collapse residual in millimetres.
        maximum_residual_mm: String,
    },
}

impl HealingChange {
    /// Source body owning this repair, including derived boundary entities.
    pub fn body(&self) -> &str {
        match self {
            Self::EdgesSewn { body, .. }
            | Self::VerticesJoined { body, .. }
            | Self::WireReordered { body, .. }
            | Self::BoundaryReversed { body, .. }
            | Self::FaceReversed { body, .. }
            | Self::SmallLineRemoved { body, .. }
            | Self::ShellsSewn { body, .. }
            | Self::TraceRebuilt { body, .. }
            | Self::WireCrossingCut { body, .. }
            | Self::DegenerateEdgeInserted { body, .. } => body,
        }
    }
}

/// A discrete check actually completed by the healer.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HealingCheck {
    /// Source body.
    pub body: String,
    /// Entity or shell scope.
    pub scope: String,
    /// Exact description of the local check; no implicit global validity.
    pub check: String,
}

/// A remaining condition to be handled by body-tier assignment.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HealingIssue {
    /// Source body.
    pub body: String,
    /// Entity or shell scope.
    pub scope: String,
    /// Stable local diagnostic.
    pub code: String,
}

/// A geometric quantity evaluated before proposing tolerance growth.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HealingMeasurement {
    /// Source body.
    pub body: String,
    /// Source entity or incidence, including endpoint index where relevant.
    pub scope: String,
    /// The measured geometric quantity and its domain.
    pub check: String,
    /// Actual attained upper distance, in millimetres.
    pub distance_upper_mm: String,
}

/// Complete repair and check history of a successfully published candidate.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HealingReport {
    /// Immutable source to which all entity names refer.
    pub source_sha256: String,
    /// Every committed topology change, in order.
    pub changes: Vec<HealingChange>,
    /// Every completed discrete check.
    pub checks: Vec<HealingCheck>,
    /// Geometry measurements behind requested local tolerances.
    pub measurements: Vec<HealingMeasurement>,
    /// Conditions which remain unresolved; these never grant solid validity.
    pub issues: Vec<HealingIssue>,
    /// Actual work for all five stages, in body/stage order.
    pub stages: Vec<HealingStage>,
    /// Work prepaid before cloning the immutable source and repair candidate.
    pub preparation_work_used: u64,
    /// Work checking and committing the gathered tolerance plan.
    pub publication_work_used: u64,
    /// Any topology repair or tolerance growth forces an approximate result.
    pub approximate: bool,
}

/// Typed refusal and work already spent; no candidate geometry is published.
#[derive(Debug)]
pub struct HealingFailure {
    /// Original refusal, including the numerical budget/deadline code.
    pub error: Error,
    /// Preparation, sewing, wire, face, shell, solid or publication.
    pub failed_stage: &'static str,
    /// Source body being processed, absent for file-wide stages.
    pub body: Option<String>,
    /// Completed stages and, if applicable, the partial failed stage.
    pub stages: Vec<HealingStage>,
    /// Work spent preparing the immutable source and candidate.
    pub preparation_work_used: u64,
    /// Work spent validating source geometry and publishing tolerance changes.
    pub publication_work_used: u64,
}

/// A healed candidate with immutable input, derived geometry and actual checks.
/// Body tiers and solid admission remain separate consumers of these checks.
#[derive(Clone, Debug)]
pub struct HealedFile {
    source: TolerantFile,
    candidate: TolerantFile,
    report: HealingReport,
}
impl HealedFile {
    /// Original source topology and tolerances, unchanged by healing.
    pub fn source(&self) -> &TolerantFile {
        &self.source
    }
    /// Derived source bodies; these do not grant solid admission.
    pub fn bodies(&self) -> &[TolerantBody] {
        &self.candidate.bodies
    }
    /// Derived definition, including explicitly recorded topology repairs.
    pub fn definition(&self, dec: &mut Decider) -> Result<StepFile, Error> {
        self.candidate.definition(dec)
    }
    /// Current entity tolerance values and growth history, scoped to the source.
    /// Reconstructing derived geometry also requires the healing report.
    pub fn tolerances(&self) -> ToleranceRecord {
        self.candidate.record()
    }
    /// Repairs, attained checks and unresolved conditions.
    pub fn report(&self) -> &HealingReport {
        &self.report
    }
}

#[derive(Clone)]
struct Growth {
    body: usize,
    entity: ToleranceEntity,
    required: Rational,
    cause: String,
}

#[derive(Default)]
struct Plan {
    growth: Vec<Growth>,
    changes: Vec<HealingChange>,
    checks: Vec<HealingCheck>,
    measurements: Vec<HealingMeasurement>,
    issues: Vec<HealingIssue>,
}
impl Plan {
    fn need(&mut self, body: usize, entity: ToleranceEntity, required: Rational, cause: &str) {
        self.growth.push(Growth {
            body,
            entity,
            required,
            cause: cause.into(),
        });
    }
    fn passed(&mut self, body: &TolerantBody, scope: &str, check: &str) {
        self.checks.push(HealingCheck {
            body: body.name.clone(),
            scope: scope.into(),
            check: check.into(),
        });
    }
    fn unresolved(&mut self, body: &TolerantBody, scope: &str, code: &str) {
        self.issues.push(HealingIssue {
            body: body.name.clone(),
            scope: scope.into(),
            code: code.into(),
        });
    }
    fn measured(&mut self, body: &TolerantBody, scope: &str, check: &str, upper: &Rational) {
        self.measurements.push(HealingMeasurement {
            body: body.name.clone(),
            scope: scope.into(),
            check: check.into(),
            distance_upper_mm: upper.to_string(),
        });
    }
}

fn failure(code: &'static str, detail: &str) -> Error {
    Refusal::new(code, detail).into()
}

pub(crate) struct NativeTopologyRepair {
    pub candidate: TolerantBody,
    pub changes: Vec<HealingChange>,
    pub checks: Vec<HealingCheck>,
    pub measurements: Vec<HealingMeasurement>,
    pub issues: Vec<HealingIssue>,
}

/// Repair source walks and shell orientations without asking for a global
/// geometry or solid certificate. The native sampled tier performs its own fits.
pub(crate) fn repair_native_topology(
    source: &TolerantBody,
    parsed: &crate::geometry_source::ParsedStep,
    scale: &Rational,
    angles: &crate::geometry_source::AngleUnit,
    curves: &mut BTreeMap<u64, std::sync::Arc<super::native::NativeCurve>>,
    stage: &mut super::admission::NativeStage,
    dec: &mut Decider,
) -> Result<NativeTopologyRepair, Error> {
    *stage = super::admission::NativeStage::BodyBuild;
    let mut decoded = source.materialize(dec)?;
    let source = decoded.to_mut();
    for edge in &mut source.edges {
        if let Curve::Retained(geometry) = &edge.curve
            && let Some(curve) = geometry.recognized_line(dec)?
        {
            edge.curve = curve;
        }
    }
    for face in source.shells.iter_mut().flat_map(|shell| &mut shell.faces) {
        if let Surface::Retained(geometry) = &face.surface
            && let Some(surface) = geometry.recognized_plane(dec)?
        {
            face.surface = surface;
        }
    }
    let count = source.vertices.len()
        + source.edges.len()
        + source
            .shells
            .iter()
            .map(|shell| shell.faces.len())
            .sum::<usize>();
    dec.charge_work(count as u64 + 1)?;
    let mut candidate = source.clone();
    let mut plan = Plan::default();
    let mut span = Rational::zero();
    for axis in 0..3 {
        let mut values = candidate
            .vertices
            .iter()
            .map(|vertex| &vertex.source.point[axis]);
        if let Some(first) = values.next() {
            let (mut low, mut high) = (first, first);
            for value in values {
                low = low.min(value);
                high = high.max(value);
            }
            span = &span + &(high - low);
        }
    }
    candidate.options.maximum_mm =
        (&span + &candidate.options.initial_mm).max(source.options.maximum_mm.clone());
    *stage = super::admission::NativeStage::WireConnection;
    wire::fix(0, source, &mut candidate, &mut plan, dec)?;
    if let Some(issue) = plan
        .issues
        .iter()
        .find(|issue| issue.code.starts_with("wire.open") || issue.code == "wire.empty-boundary")
    {
        return Err(Refusal::new(
            "import.tolerant-non-manifold",
            format!("{}: {}", issue.scope, issue.code),
        )
        .into());
    }
    candidate.options.maximum_mm = source.options.maximum_mm.clone();
    *stage = super::admission::NativeStage::Sewing;
    native_sewing::sew(
        0,
        &mut candidate,
        native_sewing::Source {
            parsed,
            scale,
            angles,
        },
        curves,
        &mut plan,
        dec,
    )?;
    *stage = super::admission::NativeStage::ShellIncidence;
    checks::shells(&mut candidate, &mut plan, dec)?;
    if let Some(issue) = plan
        .issues
        .iter()
        .find(|issue| issue.code.starts_with("shell."))
    {
        return Err(Refusal::new(
            "import.tolerant-non-manifold",
            format!("{}: {}", issue.scope, issue.code),
        )
        .into());
    }
    *stage = super::admission::NativeStage::ToleranceGrowth;
    for need in &plan.growth {
        candidate.options.maximum_mm = candidate
            .options
            .maximum_mm
            .clone()
            .max(need.required.clone());
        match &need.entity {
            ToleranceEntity::Vertex(name) => {
                candidate.grow_vertex_tolerance(name, need.required.clone(), &need.cause, dec)?;
            }
            ToleranceEntity::Edge(name) => {
                candidate.grow_edge_tolerance(name, need.required.clone(), &need.cause, dec)?;
            }
        }
    }
    Ok(NativeTopologyRepair {
        candidate,
        changes: plan.changes,
        checks: plan.checks,
        measurements: plan.measurements,
        issues: plan.issues,
    })
}

/// Repair a private candidate, then publish its geometry and tolerance changes
/// together. Cap, topology, deadline or work failure leaves `source` untouched.
/// Unresolved checks remain explicit for later surface/solid tier assignment.
pub fn heal(source: &TolerantFile, dec: &mut Decider) -> Result<HealedFile, Error> {
    attempt(source, dec).map_err(|failure| failure.error)
}

/// As [`heal`], retaining stage measurements even on cap or work refusal.
/// Failure diagnostics never contain a partially repaired candidate.
pub fn attempt(
    source: &TolerantFile,
    dec: &mut Decider,
) -> Result<HealedFile, Box<HealingFailure>> {
    if !operation_budget_active() {
        return with_operation_budget(catallax_topology::numerical::DEFAULT_WORK, || {
            attempt(source, dec)
        })
        .0;
    }
    let start = current_operation_work().used;
    let mut phase_start = start;
    let mut failed_stage = "preparation";
    let mut failed_body = None;
    let mut preparation_work_used = 0;
    let mut stages = Vec::new();
    let result = super::budget(|| {
        // Prepay the source inventory before cloning the candidate transaction.
        charge_source(source, dec)?;
        charge_source(source, dec)?;
        let original = source.clone();
        let mut candidate = source.clone();
        preparation_work_used = current_operation_work().used.saturating_sub(start);
        let mut plan = Plan::default();
        for (index, body) in candidate.bodies.iter_mut().enumerate() {
            for stage in ["sewing", "wire", "face", "shell", "solid"] {
                phase_start = current_operation_work().used;
                failed_stage = stage;
                failed_body = Some(index);
                let stage_result = match stage {
                    "sewing" => sewing::sew(index, body, &mut plan, dec),
                    "wire" => wire::fix(index, &source.bodies[index], body, &mut plan, dec),
                    "face" => checks::faces(index, body, &mut plan, dec),
                    "shell" => checks::shells(body, &mut plan, dec),
                    "solid" => checks::solid(body, &mut plan, dec),
                    _ => unreachable!(),
                };
                stages.push(HealingStage {
                    body: body.name.clone(),
                    stage: stage.into(),
                    work_used: current_operation_work().used.saturating_sub(phase_start),
                });
                stage_result?;
            }
        }
        let before = current_operation_work().used;
        phase_start = before;
        failed_stage = "publication";
        failed_body = None;
        for (index, body) in candidate.bodies.iter().enumerate() {
            checks::source_geometry(index, &source.bodies[index], body, &mut plan, dec)?;
        }
        // Gather maxima and endpoint propagation before making any ledger
        // change. Even a late refusal discards this whole private candidate.
        let mut growth: BTreeMap<(usize, bool, String), (Rational, String)> = BTreeMap::new();
        for need in &plan.growth {
            dec.charge_work(1)?;
            if need.required.is_zero() {
                continue;
            }
            let body = &candidate.bodies[need.body];
            if compare(&need.required, &body.options.maximum_mm, dec)?.is_gt() {
                return Err(failure(
                    "import.healing-cap",
                    "a measured repair exceeds its body's cap",
                ));
            }
            let (edge, name) = match &need.entity {
                ToleranceEntity::Vertex(name) => (false, name),
                ToleranceEntity::Edge(name) => (true, name),
            };
            let mut add = |is_edge: bool, name: &str| -> Result<(), Error> {
                let value = growth
                    .entry((need.body, is_edge, name.into()))
                    .or_insert_with(|| (Rational::zero(), String::new()));
                if compare(&need.required, &value.0, dec)?.is_gt() {
                    *value = (need.required.clone(), need.cause.clone());
                }
                Ok(())
            };
            add(edge, name)?;
            if edge {
                let edge = &body.edges[*body.edge_indices.get(name).ok_or_else(|| {
                    failure("import.healing-topology", "growth names an unknown edge")
                })?];
                for vertex in edge.vertices {
                    add(false, &body.vertices[vertex].source.name)?;
                }
            }
        }
        for edges_first in [true, false] {
            for ((body, edge, name), (required, cause)) in &growth {
                if *edge != edges_first {
                    continue;
                }
                if *edge {
                    candidate.bodies[*body].grow_edge_tolerance(
                        name,
                        required.clone(),
                        cause,
                        dec,
                    )?;
                } else {
                    candidate.bodies[*body].grow_vertex_tolerance(
                        name,
                        required.clone(),
                        cause,
                        dec,
                    )?;
                }
            }
        }
        dec.charge_work(candidate.bodies.len() as u64)?;
        super::charge_text(&source.source_sha256, dec)?;
        let approximate = !plan.changes.is_empty()
            || candidate
                .bodies
                .iter()
                .any(TolerantBody::has_tolerance_repairs);
        let report = HealingReport {
            source_sha256: source.source_sha256.clone(),
            changes: plan.changes,
            checks: plan.checks,
            measurements: plan.measurements,
            issues: plan.issues,
            stages: std::mem::take(&mut stages),
            preparation_work_used,
            publication_work_used: current_operation_work().used.saturating_sub(before),
            approximate,
        };
        let result = HealedFile {
            source: original,
            candidate,
            report,
        };
        if let Err(error) = dec.check_deadline() {
            stages = result.report.stages;
            return Err(error.into());
        }
        Ok(result)
    });
    result.map_err(|error| {
        Box::new(HealingFailure {
            error,
            failed_stage,
            body: failed_body.map(|index| source.bodies[index].name.clone()),
            stages,
            preparation_work_used: if failed_stage == "preparation" {
                current_operation_work().used.saturating_sub(start)
            } else {
                preparation_work_used
            },
            publication_work_used: if failed_stage == "publication" {
                current_operation_work().used.saturating_sub(phase_start)
            } else {
                0
            },
        })
    })
}

fn arithmetic(a: &Rational, b: &Rational, op: char, dec: &mut Decider) -> Result<Rational, Error> {
    dec.reserve_exact_arithmetic(a, b, matches!(op, '*' | '/'))?;
    Ok(match op {
        '+' => a + b,
        '-' => a - b,
        '*' => a * b,
        '/' => a / b,
        _ => unreachable!(),
    })
}

fn compare(a: &Rational, b: &Rational, dec: &mut Decider) -> Result<std::cmp::Ordering, Error> {
    dec.reserve_exact_arithmetic(a, b, true)?;
    dec.reserve_exact_arithmetic(a, b, true)?;
    Ok(a.cmp(b))
}

fn sub(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Vec3, Error> {
    Ok([
        arithmetic(&a[0], &b[0], '-', dec)?,
        arithmetic(&a[1], &b[1], '-', dec)?,
        arithmetic(&a[2], &b[2], '-', dec)?,
    ])
}
fn dot(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Rational, Error> {
    let mut value = Rational::zero();
    for i in 0..3 {
        let product = arithmetic(&a[i], &b[i], '*', dec)?;
        value = arithmetic(&value, &product, '+', dec)?;
    }
    Ok(value)
}

fn distance_squared(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Rational, Error> {
    let d = sub(a, b, dec)?;
    dot(&d, &d, dec)
}

fn distance_upper(
    squared: &Rational,
    cap: &Rational,
    dec: &mut Decider,
) -> Result<Rational, Error> {
    let squared_cap = arithmetic(cap, cap, '*', dec)?;
    if compare(squared, &squared_cap, dec)?.is_gt() {
        return Err(failure(
            "import.healing-cap",
            "measured source correspondence exceeds the healing cap",
        ));
    }
    charge_value(squared, dec)?;
    let bits = squared
        .numerator()
        .bit_length()
        .saturating_add(squared.denominator().bit_length());
    dec.charge_work(
        bits.saturating_mul(bits)
            .saturating_mul(4)
            .saturating_add(1),
    )?;
    // The cheap exact squared comparison also proves cap as an upper bound.
    // It can tighten an outward-rounded square root at the cap itself.
    Ok(Ball::from_rational(squared, 128)
        .sqrt()
        .ok_or_else(|| {
            failure(
                "import.healing-measurement",
                "distance could not be enclosed",
            )
        })?
        .hi()
        .min(cap.clone()))
}

pub(super) fn charge_value(value: &Rational, dec: &mut Decider) -> Result<(), Error> {
    let bits = value
        .numerator()
        .bit_length()
        .saturating_add(value.denominator().bit_length());
    dec.charge_work(1_u64.saturating_add(bits))?;
    Ok(())
}

fn charge_vector(value: &Vec3, dec: &mut Decider) -> Result<(), Error> {
    for value in value {
        charge_value(value, dec)?;
    }
    Ok(())
}

fn charge_placement(value: &Placement, dec: &mut Decider) -> Result<(), Error> {
    for value in [&value.origin, &value.axis, &value.ref_dir] {
        charge_vector(value, dec)?;
    }
    Ok(())
}

fn charge_curve(value: &Curve, dec: &mut Decider) -> Result<(), Error> {
    match value {
        Curve::Retained(_) => dec.charge_work(1)?,
        Curve::Line(a, b) => {
            charge_vector(a, dec)?;
            charge_vector(b, dec)?;
        }
        Curve::Circle(p, r) => {
            charge_placement(p, dec)?;
            charge_value(r, dec)?;
        }
        Curve::Ellipse(p, a, b) => {
            charge_placement(p, dec)?;
            charge_value(a, dec)?;
            charge_value(b, dec)?;
        }
        Curve::BSpline(spline) => {
            for value in spline
                .geometry
                .knots()
                .iter()
                .chain(spline.geometry.weights())
            {
                charge_value(value, dec)?;
            }
            for point in spline.geometry.controls() {
                charge_vector(point, dec)?;
            }
        }
    }
    Ok(())
}

fn charge_surface(value: &Surface, dec: &mut Decider) -> Result<(), Error> {
    match value {
        Surface::Retained(_) => dec.charge_work(1)?,
        Surface::Plane(p) => charge_placement(p, dec)?,
        Surface::Cylinder(p, r) | Surface::Cone45(p, r) | Surface::Sphere(p, r) => {
            charge_placement(p, dec)?;
            charge_value(r, dec)?;
        }
        Surface::Torus(p, a, b) => {
            charge_placement(p, dec)?;
            charge_value(a, dec)?;
            charge_value(b, dec)?;
        }
        Surface::Cone(p, r, angle) => {
            charge_placement(p, dec)?;
            charge_value(r, dec)?;
            charge_value(&angle.value, dec)?;
        }
        Surface::BSpline(spline) => {
            for knots in spline.geometry.knots() {
                for value in knots {
                    charge_value(value, dec)?;
                }
            }
            for value in spline.geometry.weights() {
                charge_value(value, dec)?;
            }
            for point in spline.geometry.controls() {
                charge_vector(point, dec)?;
            }
        }
    }
    Ok(())
}

fn charge_source(source: &TolerantFile, dec: &mut Decider) -> Result<(), Error> {
    super::charge_text(&source.source_sha256, dec)?;
    charge_value(&source.mm_per_unit, dec)?;
    for body in &source.bodies {
        charge_body(body, dec)?;
    }
    Ok(())
}

pub(super) fn charge_body(body: &TolerantBody, dec: &mut Decider) -> Result<(), Error> {
    super::charge_text(&body.name, dec)?;
    charge_value(&body.options.initial_mm, dec)?;
    charge_value(&body.options.maximum_mm, dec)?;
    charge_value(&body.maximum_attained_mm, dec)?;
    for vertex in &body.vertices {
        dec.charge_work(
            3_u64
                .saturating_mul(vertex.source.name.len() as u64)
                .saturating_add(1),
        )?;
        charge_vector(&vertex.source.point, dec)?;
        charge_value(&vertex.tolerance_mm, dec)?;
    }
    for edge in &body.edges {
        dec.charge_work(
            3_u64
                .saturating_mul(edge.name.len() as u64)
                .saturating_add(1),
        )?;
        charge_curve(&edge.curve, dec)?;
        charge_value(&edge.tolerance_mm, dec)?;
    }
    for shell in &body.shells {
        super::charge_text(&shell.name, dec)?;
        for face in &shell.faces {
            super::charge_text(&face.name, dec)?;
            charge_surface(&face.surface, dec)?;
            for boundary in &face.loops {
                dec.charge_work(1 + boundary.edges.len() as u64)?;
            }
        }
    }
    for change in &body.changes {
        super::charge_change(change, dec)?;
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct EdgeUse {
    shell: usize,
    face: usize,
    forward: bool,
}

fn usage(body: &TolerantBody, dec: &mut Decider) -> Result<Vec<Vec<EdgeUse>>, Error> {
    let mut result = vec![Vec::new(); body.edges.len()];
    for (s, shell) in body.shells.iter().enumerate() {
        for (f, face) in shell.faces.iter().enumerate() {
            for boundary in &face.loops {
                for coedge in &boundary.edges {
                    dec.charge_work(1)?;
                    result[coedge.edge].push(EdgeUse {
                        shell: s,
                        face: f,
                        forward: coedge.forward == boundary.orientation,
                    });
                }
            }
        }
    }
    Ok(result)
}

fn endpoints(body: &TolerantBody, coedge: &TolerantUse) -> [usize; 2] {
    let ends = body.edges[coedge.edge].vertices;
    if coedge.forward {
        ends
    } else {
        [ends[1], ends[0]]
    }
}
