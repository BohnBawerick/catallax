//! Default tolerant B-reps: measured local tolerances and sampled native fits.
//! Certified construction remains a separate, higher-authority result.

use std::sync::Arc;

use catallax_numbers::{Decider, Rational, Refusal};
use serde::Serialize;

use crate::ImportRefusal;
use crate::geometry_source::{AngleUnit, ParsedStep};

use super::{Error, ToleranceOptions, TolerantBody, TolerantFile, healing};

mod evidence;
mod fit;
mod step;

/// A local tolerance, in millimetres, backed by the stated sampled measurement.
#[derive(Clone, Debug, Serialize)]
pub struct EntityTolerance {
    /// Original source body identity.
    pub body: String,
    /// Vertex, edge or face.
    pub kind: &'static str,
    /// Original entity identity within the immutable source file.
    pub entity: String,
    /// Current local tolerance, including evaluation uncertainty.
    pub tolerance_mm: f64,
    /// Largest observed fit residual, not a whole-domain deviation bound.
    pub maximum_sample_residual_mm: f64,
    /// Number of measured sample incidences.
    pub samples: usize,
}

/// An actual sampled check. It supplies no unsampled geometric theorem.
#[derive(Clone, Debug, Serialize)]
pub struct SampledCheck {
    /// Source body.
    pub body: String,
    /// Original edge or face scope.
    pub scope: String,
    /// Check method and its sampled authority.
    pub check: &'static str,
    /// Number of actual samples.
    pub samples: usize,
    /// Largest observed residual including evaluation uncertainty.
    pub maximum_residual_mm: f64,
    /// Fits which reached the numerical stopping width.
    pub converged: usize,
}

/// One finite edge/surface correspondence sample, in original source parameters.
#[derive(Clone, Debug, Serialize)]
pub struct UvSample {
    /// Fraction along the finite source edge's chosen parameter interval.
    pub fraction: f64,
    /// Native source curve chart and parameter.
    pub curve: super::native::NativeLocation,
    /// Native surface chart and UV.
    pub surface: super::native::NativeLocation,
    /// Final physical curve correspondence residual in millimetres. In a
    /// body taken with the tolerance its file recorded, that tolerance:
    /// nothing was measured.
    pub residual_mm: f64,
    /// Final topology-vertex incidence at a boundary endpoint, in millimetres.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vertex_residual_mm: Option<f64>,
    /// Numerical convergence status, not a proof of a nearest point.
    pub converged: bool,
}

/// A finite native trim, retaining original source geometry and its UV samples.
#[derive(Clone, Debug, Serialize)]
pub struct NativeTrim {
    /// Source body containing this edge use.
    pub body: String,
    /// Source face.
    pub face: String,
    /// Source boundary index.
    pub boundary: usize,
    /// Position in the repaired boundary, including inserted pole edges.
    pub use_index: usize,
    /// Original shared edge name.
    pub edge: String,
    /// Traversal against the source EDGE_CURVE.
    pub forward: bool,
    /// Original PCURVE used as a seed, if available.
    pub source_pcurve: Option<u64>,
    /// Actual native correspondence samples.
    pub samples: Vec<UvSample>,
    /// Rebuilt boundary trace, rather than a retained source PCURVE.
    pub reconstructed: bool,
    /// Plane projection keeps the source spline's exact knots and weights.
    pub projected: Option<ProjectedTrace>,
    /// The source PCURVE at the parameter of its 3D curve, where the gap
    /// rule for imports gave the use one (decision Q9). It replaces the
    /// source PCURVE as the use's face copy.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reparametrised: Option<ReparametrisedTrace>,
}

/// A face copy in the parameter of its 3D curve: a polynomial B-spline in the
/// surface's original UV units, with knots in the 3D curve's native
/// parameter.
#[derive(Clone, Debug, Serialize)]
pub struct ReparametrisedTrace {
    /// Polynomial degree.
    pub degree: usize,
    /// Full knot vector, ascending.
    pub knots: Vec<f64>,
    /// Controls in the surface's original UV units.
    pub controls: Vec<[f64; 2]>,
    /// Spans the change of parameter took.
    pub spans: usize,
}

impl ReparametrisedTrace {
    /// The UV at a parameter of the 3D curve, by de Boor's recurrence
    /// (Piegl and Tiller, "The NURBS Book", algorithm A3.1). Past either end
    /// the end span's polynomial continues.
    pub fn uv(&self, parameter: f64) -> [f64; 2] {
        let (degree, knots) = (self.degree, &self.knots);
        let last = self.controls.len() - 1;
        let mut span = degree;
        while span < last && parameter >= knots[span + 1] {
            span += 1;
        }
        let mut points: Vec<[f64; 2]> = self.controls[span - degree..=span].to_vec();
        for round in 1..=degree {
            for index in (round..=degree).rev() {
                let low = knots[span + index - degree];
                let high = knots[span + 1 + index - round];
                let share = (parameter - low) / (high - low);
                points[index] = std::array::from_fn(|axis| {
                    (1.0 - share) * points[index - 1][axis] + share * points[index][axis]
                });
            }
        }
        points[degree]
    }
}

/// What the gap rule for imports read of one edge use that carries a source
/// PCURVE (decision Q9, and section 4.2 step 2 of the Checkpoint 2
/// decisions): "measure the geometric gap both ways, reparametrise, then
/// continue with the same-parameter gap".
#[derive(Clone, Debug, Serialize)]
pub struct NativeImportGap {
    /// Source body.
    pub body: String,
    /// Source face.
    pub face: String,
    /// Source edge.
    pub edge: String,
    /// Source boundary index.
    pub boundary: usize,
    /// Position in the boundary.
    pub use_index: usize,
    /// The source PCURVE that was read.
    pub source_pcurve: u64,
    /// Rule T4's tolerance of the same-parameter gap of the source PCURVE,
    /// as the A+ checker read it, in millimetres.
    pub parametric_mm: f64,
    /// The geometric gap, the larger of its two directions with its
    /// rounding error, in millimetres. Absent where the rule leaves the use
    /// as it is: a tolerance at rule T1's least floor, or a straight
    /// residual.
    pub geometric_mm: Option<f64>,
    /// Rule T4's tolerance of the same-parameter gap of the copy in the 3D
    /// curve's parameter, where one was made, in millimetres.
    pub reparametrised_mm: Option<f64>,
    /// `as-imported`, `geometric`, `reparametrised` or `not-lowered`.
    pub outcome: &'static str,
    /// Whether the copy the rule read is the use's face copy. `false` where
    /// this admission then rebuilt the trace from the 3D curve, so the
    /// reading describes a copy the use no longer has.
    pub kept: bool,
    /// Work this use counted against the one allowance of its edge.
    pub work: u64,
}

/// A plane trace that preserves the source curve's parameterization.
#[derive(Clone, Debug, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ProjectedTrace {
    /// Projected source spline with unchanged exact knots and weights.
    Spline {
        /// Source 3D B-spline basis.
        source_curve: u64,
        /// Measured native plane coordinates of its controls.
        controls: Vec<[f64; 2]>,
    },
    /// Analytic conic in the plane's original UV units.
    Conic {
        /// UV centre.
        centre: [f64; 2],
        /// Unit UV directions for cosine and sine, including reflected frames.
        axes: [[f64; 2]; 2],
        /// Radii along those directions.
        radii: [f64; 2],
    },
}

/// Sampled natural closure of an edgeless native face.
#[derive(Clone, Debug, Serialize)]
pub struct NativeDomain {
    /// Source body.
    pub body: String,
    /// Source face.
    pub face: String,
    /// Original source parameter limits, one [minimum, maximum] per axis.
    pub bounds: [[f64; 2]; 2],
}

/// A numerical warning retained on the result rather than promoted to a theorem.
#[derive(Clone, Debug, Serialize)]
pub struct NativeWarning {
    /// Owning source body; absent only for document-wide authority limitations.
    pub body: Option<String>,
    /// Source scope.
    pub scope: String,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Actual limitation of the measurement or source interpretation.
    pub detail: String,
}

/// Complete lower-tier provenance. Geometry and source hashes stay immutable.
#[derive(Clone, Debug, Serialize)]
pub struct TolerantEvidence {
    /// Versioned numerical method.
    pub method: &'static str,
    /// Immutable source identity.
    pub source_sha256: String,
    /// The tier never grants certified Model authority.
    pub authority: &'static str,
    /// All measured local entity tolerances.
    pub tolerances: Vec<EntityTolerance>,
    /// Completed numerical checks.
    pub sampled_checks: Vec<SampledCheck>,
    /// Actual topology repairs.
    pub repairs: Vec<healing::HealingChange>,
    /// Completed discrete topology checks.
    pub topology_checks: Vec<healing::HealingCheck>,
    /// Original measured topology-repair displacements.
    pub topology_measurements: Vec<healing::HealingMeasurement>,
    /// Explicit remaining limitations, not suppressed admission errors.
    pub warnings: Vec<NativeWarning>,
    /// Finite source-native edge/face correspondences.
    pub native_trims: Vec<NativeTrim>,
    /// Natural domains whose opposing rims or poles were sampled for closure.
    pub native_domains: Vec<NativeDomain>,
    /// The gap rule's reading of every edge use that carries a source
    /// PCURVE the checker reads.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub import_gaps: Vec<NativeImportGap>,
}

impl TolerantEvidence {
    fn checkpoint(&self) -> [usize; 9] {
        [
            self.tolerances.len(),
            self.sampled_checks.len(),
            self.repairs.len(),
            self.topology_checks.len(),
            self.topology_measurements.len(),
            self.warnings.len(),
            self.native_trims.len(),
            self.native_domains.len(),
            self.import_gaps.len(),
        ]
    }
    fn rollback(&mut self, lengths: [usize; 9]) {
        self.tolerances.truncate(lengths[0]);
        self.sampled_checks.truncate(lengths[1]);
        self.repairs.truncate(lengths[2]);
        self.topology_checks.truncate(lengths[3]);
        self.topology_measurements.truncate(lengths[4]);
        self.warnings.truncate(lengths[5]);
        self.native_trims.truncate(lengths[6]);
        self.native_domains.truncate(lengths[7]);
        self.import_gaps.truncate(lengths[8]);
    }
}

/// Counts of a closed, oriented source-topology candidate.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct TolerantCounts {
    /// Physical faces.
    pub faces: usize,
    /// Active physical edges.
    pub edges: usize,
    /// Active physical vertices.
    pub vertices: usize,
    /// Closed shell count.
    pub shells: usize,
}

/// A sampled tolerant solid. No exact Model extraction is available.
#[derive(Debug)]
pub struct TolerantSolid {
    pub(crate) body: TolerantBody,
    counts: TolerantCounts,
    source_index: usize,
    recorded_mm: Option<Rational>,
}
impl TolerantSolid {
    /// The tolerance the file recorded at export, where this body was taken
    /// with it and nothing was measured (check once, roadmap item 3.0A
    /// step 5). `None` for a body that was measured.
    pub fn recorded_tolerance_mm(&self) -> Option<&Rational> {
        self.recorded_mm.as_ref()
    }
    /// Source identity of this body's definition.
    pub fn name(&self) -> &str {
        self.body.name()
    }
    /// Discrete topology counts, not an embedding certificate.
    pub fn counts(&self) -> TolerantCounts {
        self.counts
    }
    /// Native numerical source topology for display and downstream tolerant work.
    pub fn definition(&self) -> &TolerantBody {
        &self.body
    }
    /// Position in the immutable source inventory, independent of refused siblings.
    pub fn source_index(&self) -> usize {
        self.source_index
    }
}

/// Rebuilt native STEP bytes and their exact source-body root bindings.
/// The file has sampled tolerant authority, not certified Model authority.
#[derive(Debug)]
pub struct NativeStepFile {
    bytes: Vec<u8>,
    body_roots: Vec<Vec<u64>>,
    checks: Vec<SampledCheck>,
}
impl NativeStepFile {
    /// Actual serialized bytes.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
    /// STEP instance IDs of every material root per source body, in input order.
    pub fn source_body_roots(&self) -> &[Vec<u64>] {
        &self.body_roots
    }
    /// Export-time native UV correspondence measurements.
    pub fn checks(&self) -> &[SampledCheck] {
        &self.checks
    }
    /// Consume the file without copying its bytes.
    pub fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

#[derive(Serialize)]
struct ExportEvidence<'a> {
    export_scope: &'static str,
    source_bodies: Option<&'a [usize]>,
    source_body_roots: &'a [Vec<u64>],
    source_document_evidence: evidence::Scope<'a>,
    export_checks: &'a [SampledCheck],
}

/// The first native construction stage that prevents a body from being used.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum NativeStage {
    /// The shared request allowance expired before this body could start.
    RequestLimit,
    /// Decode this body's retained source topology and carriers.
    BodyBuild,
    /// Order and reconnect its source wires.
    WireConnection,
    /// Sew mutually matching free boundaries within the sewing cap.
    Sewing,
    /// Check connected shell incidence and orientation.
    ShellIncidence,
    /// Check the body's closed manifold incidence.
    ManifoldIncidence,
    /// Fit finite 3D edges and their endpoints.
    CurveFit,
    /// Cut a measured wire crossing.
    CrossingCut,
    /// Fit face coordinates to the source carriers.
    SurfaceFit,
    /// Reconstruct a missing or inconsistent face trace.
    TraceRebuild,
    /// Insert a measured collapsed pole boundary.
    PoleBoundary,
    /// Commit measured local entity tolerances.
    ToleranceGrowth,
}

/// A body-local first blocker. Resource refusals remain resource observations.
#[derive(Clone, Debug)]
pub struct NativeBodyFailure {
    /// Actual repair or check active at the first refusal.
    pub stage: NativeStage,
    /// Original typed failure, not inferred from a file-level outcome.
    pub error: Error,
    /// Work spent on this body before the refusal.
    pub work_used: u64,
}
impl NativeBodyFailure {
    /// A work, deadline or numerical-resolution limit, not a geometry verdict.
    pub fn resource_limited(&self) -> bool {
        super::native::control(&self.error)
    }
}

/// A tolerant source document, including its original occurrence graph.
#[derive(Debug)]
pub struct TolerantImport {
    parsed: Arc<ParsedStep>,
    source: Arc<TolerantFile>,
    bodies: Vec<TolerantSolid>,
    evidence: TolerantEvidence,
    carriers: fit::Carriers,
    body_errors: Vec<Option<NativeBodyFailure>>,
}
impl TolerantImport {
    /// Immutable supported source definitions. Occurrence links remain in the document.
    pub fn source(&self) -> &TolerantFile {
        &self.source
    }
    /// Tolerant body definitions, in source-model order.
    pub fn bodies(&self) -> &[TolerantSolid] {
        &self.bodies
    }
    /// Typed body-local failures. Successful siblings remain usable.
    pub fn body_failures(&self) -> impl Iterator<Item = (usize, &NativeBodyFailure)> {
        self.body_errors
            .iter()
            .enumerate()
            .filter_map(|(index, error)| error.as_ref().map(|error| (index, error)))
    }
    /// First blocker for one source-model slot, without renumbering siblings.
    pub fn body_failure(&self, source_index: usize) -> Option<&NativeBodyFailure> {
        self.body_errors.get(source_index).and_then(Option::as_ref)
    }

    pub(crate) fn checked_definition(
        &self,
        index: usize,
        dec: &mut Decider,
    ) -> Result<crate::StepFile, Error> {
        if self.carriers.angles != AngleUnit::Radians {
            return Err(invalid(
                "import.checked-angular-context",
                "checked carrier conversion requires radians",
            ));
        }
        let mut file = Self::recognize_definition(self.body_definition(index, dec)?, dec)?;
        for model in file.surface_models.drain(..) {
            if model.shells.len() != 1 {
                return Err(invalid(
                    "import.checked-body-shells",
                    "checked upgrade requires one closed material shell",
                ));
            }
            file.solids.push(crate::Solid {
                name: model.name,
                faces: model.shells.into_iter().next().expect("one shell").faces,
            });
        }
        Ok(file)
    }
    pub(crate) fn checked_source_definition(
        &self,
        index: usize,
        dec: &mut Decider,
    ) -> Result<crate::StepFile, Error> {
        if self.carriers.angles != AngleUnit::Radians {
            return Err(invalid(
                "import.checked-angular-context",
                "checked carrier conversion requires radians",
            ));
        }
        Self::recognize_definition(
            self.source
                .definition_of(&self.source.bodies[index..=index], dec)?,
            dec,
        )
    }

    fn recognize_definition(
        mut file: crate::StepFile,
        dec: &mut Decider,
    ) -> Result<crate::StepFile, Error> {
        for face in file
            .solids
            .iter_mut()
            .flat_map(|body| &mut body.faces)
            .chain(
                file.surface_models
                    .iter_mut()
                    .flat_map(|body| &mut body.shells)
                    .flat_map(|shell| &mut shell.faces),
            )
        {
            for edge in face
                .loops
                .iter_mut()
                .flat_map(|boundary| &mut boundary.edges)
            {
                if let crate::Curve::Retained(source) = &edge.curve {
                    edge.curve = source.recognized_curve(dec)?;
                }
            }
            if let crate::Surface::Retained(source) = &face.surface {
                face.surface = source.recognized_surface(&face.loops, dec)?;
            }
        }
        Ok(file)
    }
    /// Reconstruct one admitted native candidate as source geometry, without
    /// converting it into a certified kernel model.
    pub fn body_definition(
        &self,
        index: usize,
        dec: &mut Decider,
    ) -> Result<crate::StepFile, Error> {
        let body = self
            .bodies
            .get(index)
            .ok_or_else(|| invalid("import.tolerant-body-index", "body index is out of range"))?;
        self.source
            .definition_of(std::slice::from_ref(&body.body), dec)
    }
    /// All measured tolerances, samples, repairs and limitations.
    pub fn evidence(&self) -> &TolerantEvidence {
        &self.evidence
    }
    /// Actual immutable source identity.
    pub fn source_sha256(&self) -> &str {
        &self.evidence.source_sha256
    }
    /// Rebuild the source document's B-reps and preserve occurrence placement.
    /// The file explicitly carries the tolerant label and per-entity evidence.
    pub fn export_step(&self, dec: &mut Decider) -> Result<NativeStepFile, Error> {
        if self.bodies.len() != self.source.bodies.len() {
            return Err(invalid(
                "import.tolerant-incomplete-document",
                "export usable bodies separately while siblings remain source surfaces",
            ));
        }
        let built = step::rebuild(self, &self.bodies, dec)?;
        let metadata = serde_json::to_string(&ExportEvidence {
            export_scope: "source_document",
            source_bodies: None,
            source_body_roots: &built.body_roots,
            source_document_evidence: evidence::Scope::new(&self.evidence, None),
            export_checks: &built.checks,
        })
        .map_err(|error| Refusal::new("import.tolerant-evidence", error.to_string()))?;
        let bytes = self.parsed.write_step_graph(
            &built.changes,
            self.parsed.schemas()?,
            &metadata,
            None,
            crate::geometry_source::StepGraphKind::Native,
            dec,
        )?;
        Ok(NativeStepFile {
            bytes,
            body_roots: built.body_roots,
            checks: built.checks,
        })
    }
    /// Export selected body-local geometry. Product occurrences and unrelated
    /// bodies are excluded; use `export_step` to preserve the source document.
    pub fn export_body(&self, index: usize, dec: &mut Decider) -> Result<NativeStepFile, Error> {
        let body = self
            .bodies
            .get(index)
            .ok_or_else(|| invalid("import.tolerant-body-index", "body index is out of range"))?;
        let mut built = step::rebuild(self, std::slice::from_ref(body), dec)?;
        let root = step::standalone_frame(
            self,
            &mut built.changes,
            body.name(),
            &built.body_roots[0],
            dec,
        )?;
        let metadata = serde_json::to_string(&ExportEvidence {
            export_scope: "body_local_geometry",
            source_bodies: Some(&[body.source_index]),
            source_body_roots: &built.body_roots,
            source_document_evidence: evidence::Scope::new(&self.evidence, Some(body.name())),
            export_checks: &built.checks,
        })
        .map_err(|error| Refusal::new("import.tolerant-evidence", error.to_string()))?;
        let bytes = self.parsed.write_step_graph(
            &built.changes,
            self.parsed.schemas()?,
            &metadata,
            Some(&[root]),
            crate::geometry_source::StepGraphKind::Native,
            dec,
        )?;
        Ok(NativeStepFile {
            bytes,
            body_roots: built.body_roots,
            checks: built.checks,
        })
    }
    /// Source-to-OCCT dimensional factors for independent native export checks.
    /// These describe source units, not an adjustment made by a repair.
    pub fn source_unit_factors(&self) -> [f64; 3] {
        let angle = match &self.carriers.angles {
            AngleUnit::Radians => 1.0,
            AngleUnit::Degrees => std::f64::consts::PI / 180.0,
            AngleUnit::ScaledRadians(scale) => scale.to_f64_nearest(),
        };
        [self.source.mm_per_unit.to_f64_nearest(), angle, 1.0]
    }
    /// Original parsed geometry and occurrence declarations for native display.
    pub fn parsed_source(&self) -> &ParsedStep {
        &self.parsed
    }
    /// Already compiled physical source surface, in its original native UV.
    pub fn native_surface(&self, entity: u64) -> Option<&super::native::NativeSurface> {
        self.carriers.surfaces.get(&entity).map(Arc::as_ref)
    }
    /// Already compiled source curve, without a finite-edge proof requirement.
    pub fn native_curve(&self, entity: u64) -> Option<&super::native::NativeCurve> {
        self.carriers.curves.get(&entity).map(Arc::as_ref)
    }
}

fn invalid(code: &'static str, detail: impl Into<String>) -> Error {
    Refusal::new(code, detail).into()
}

fn source_angles(contexts: &[crate::geometry_source::ContextUnits]) -> Result<AngleUnit, Error> {
    let mut selected = None;
    for units in contexts {
        let angles = match units.plane_angle.clone().unwrap_or(AngleUnit::Radians) {
            AngleUnit::ScaledRadians(scale) if scale == Rational::one() => AngleUnit::Radians,
            angles => angles,
        };
        if selected.as_ref().is_some_and(|prior| prior != &angles) {
            return Err(invalid(
                "import.tolerant-ambiguous-units",
                "3D representation contexts assign different plane-angle units",
            ));
        }
        selected = Some(angles);
    }
    Ok(selected.unwrap_or(AngleUnit::Radians))
}

/// Admit a separately labelled tolerant document without whole-domain proof gates.
/// Request deadlines and externally installed budgets remain effective.
///
/// A file that records its tolerance is taken with it (check once); see
/// [`crate::recorded`] and [`TolerantSolid::recorded_tolerance_mm`].
pub fn import_step(
    bytes: &[u8],
    expected_sha256: &str,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<TolerantImport, Error> {
    let (parsed, source, angles) = prepare(bytes, expected_sha256, options, dec)?;
    // The record is in the file's length unit, which this reading has.
    let recorded = crate::recorded::written_tolerance(bytes)
        .and_then(|written| written.in_millimetres(&source.mm_per_unit));
    let imported = admit(
        parsed,
        source,
        angles,
        AdmissionPolicy::Repair,
        recorded.as_ref().map(|record| &record.tolerance_mm),
        dec,
    )?;
    if imported.bodies.is_empty() {
        return Err(imported
            .body_errors
            .into_iter()
            .flatten()
            .next()
            .map(|failure| failure.error)
            .unwrap_or_else(|| {
                invalid(
                    "import.tolerant-broken-topology",
                    "source contains no usable body",
                )
            }));
    }
    Ok(imported)
}

pub(crate) fn prepare(
    bytes: &[u8],
    expected_sha256: &str,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<(Arc<ParsedStep>, Arc<TolerantFile>, AngleUnit), Error> {
    if options.initial_mm < Rational::zero()
        || options.maximum_mm <= Rational::zero()
        || options.initial_mm > options.maximum_mm
    {
        return Err(invalid(
            "import.invalid-tolerances",
            "invalid source tolerance controls",
        ));
    }
    let parsed = Arc::new(ParsedStep::read(bytes, dec)?);
    if parsed.source_sha256() != expected_sha256 {
        return Err(ImportRefusal::HashMismatch {
            expected: expected_sha256.into(),
            actual: parsed.source_sha256().into(),
        }
        .into());
    }
    parsed.schemas()?;
    let (scale, declarations, contexts) = parsed.native_inventory(dec)?;
    let angles = source_angles(&contexts)?;
    let mut source = TolerantFile {
        source_sha256: expected_sha256.into(),
        mm_per_unit: scale,
        options: options.clone(),
        bodies: Vec::with_capacity(declarations.len()),
    };
    for declaration in declarations {
        source.bodies.push(TolerantBody::retained(
            &format!("#{}", declaration.entity),
            declaration.kind == crate::geometry_source::ModelDeclarationKind::Solid,
            parsed.retained_body(declaration.entity, &source.mm_per_unit),
            &options,
        ));
    }
    Ok((parsed, Arc::new(source), angles))
}

pub(crate) enum AdmissionPolicy<'a> {
    Repair,
    MainFirst(&'a [bool]),
}

/// One pass over the bodies still to be admitted.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Round {
    /// Check once: the file records its tolerance, so each use's face copy
    /// is read as written and nothing is measured. A body this pass cannot
    /// take is measured by the passes after it.
    Recorded,
    /// Every entity is measured (rule T3), with boundary repair or without.
    Measured { repair_boundaries: bool },
}

/// Whether an edge's hard work cap refused this body (Checkpoint 2
/// decisions, section 4.2 step 8: "A hard per-edge work cap refuses by
/// name").
fn edge_work_cap(error: &Error) -> bool {
    matches!(error, Error::Numerical(refusal)
        if refusal.code == catallax_topology::numerical::gap::WORK_CAP)
}

/// Runs the rounds over the source bodies. `attempt` builds and checks one
/// body in one round. `needed`, where given, says which bodies are tried at
/// all. A round is not started once `exhausted` says the request's allowance
/// is spent.
///
/// A body one round admits is not tried again, and a body one round refuses
/// is tried by the next, with one exception. Every round gives each edge of
/// the body a fresh work allowance, so a body that an edge refused by its
/// work cap would get a second allowance for that edge in the next round, and
/// the repair round may read the edge through another carrier. The cap is
/// hard, so that refusal ends the body: it is not tried again and its
/// refusal stands. Other bodies are not affected.
///
/// Returns the admitted bodies in the order they were admitted and, by
/// source index, the last refusal of each body that was not.
fn over_rounds<T>(
    bodies: usize,
    rounds: &[Round],
    needed: Option<&[bool]>,
    mut attempt: impl FnMut(Round, usize) -> Result<T, NativeBodyFailure>,
    exhausted: impl Fn() -> bool,
) -> (Vec<T>, Vec<Option<NativeBodyFailure>>) {
    let mut admitted = Vec::with_capacity(bodies);
    let mut refusals: Vec<Option<NativeBodyFailure>> = (0..bodies).map(|_| None).collect();
    // Admitted, or refused for good.
    let mut ended = vec![false; bodies];
    for round in rounds {
        for index in 0..bodies {
            if ended[index] || needed.is_some_and(|needed| !needed[index]) {
                continue;
            }
            match attempt(*round, index) {
                Ok(body) => {
                    ended[index] = true;
                    refusals[index] = None;
                    admitted.push(body);
                }
                Err(failure) => {
                    ended[index] = edge_work_cap(&failure.error);
                    refusals[index] = Some(failure);
                }
            }
        }
        if exhausted() {
            break;
        }
    }
    (admitted, refusals)
}

/// `recorded_mm` is the tolerance the file recorded at export, if any
/// ([`crate::recorded`]).
pub(crate) fn admit(
    parsed: Arc<ParsedStep>,
    source: Arc<TolerantFile>,
    angles: AngleUnit,
    policy: AdmissionPolicy<'_>,
    recorded_mm: Option<&Rational>,
    dec: &mut Decider,
) -> Result<TolerantImport, Error> {
    let edge_work = catallax_topology::numerical::gap::EDGE_WORK;
    admit_under(edge_work, parsed, source, angles, policy, recorded_mm, dec)
}

/// [`admit`] with `edge_work` units for each edge's gap readings in place of
/// the checker's cap. Only the test of the cap gives less.
fn admit_under(
    edge_work: u64,
    parsed: Arc<ParsedStep>,
    source: Arc<TolerantFile>,
    angles: AngleUnit,
    policy: AdmissionPolicy<'_>,
    recorded_mm: Option<&Rational>,
    dec: &mut Decider,
) -> Result<TolerantImport, Error> {
    let mut fits = fit::Context::new(
        &parsed,
        &source.mm_per_unit,
        angles.clone(),
        Default::default(),
        source.options.maximum_mm.to_f64_nearest(),
    );
    fits.edge_work = edge_work;
    let mut evidence = TolerantEvidence {
        method: "catallax-tolerant-native-v1",
        source_sha256: source.source_sha256.clone(),
        authority: "tolerant: discrete topology and sampled native fits; no whole-domain certificate",
        tolerances: Vec::new(),
        sampled_checks: Vec::new(),
        repairs: Vec::new(),
        topology_checks: Vec::new(),
        topology_measurements: Vec::new(),
        warnings: Vec::new(),
        native_trims: Vec::new(),
        native_domains: Vec::new(),
        import_gaps: Vec::new(),
    };
    let (measured, incomplete): (&[bool], Option<&[bool]>) = match policy {
        AdmissionPolicy::Repair => (&[true], None),
        AdmissionPolicy::MainFirst(incomplete) => (&[false, true], Some(incomplete)),
    };
    let rounds: Vec<Round> = recorded_mm
        .map(|_| Round::Recorded)
        .into_iter()
        .chain(
            measured
                .iter()
                .map(|&repair_boundaries| Round::Measured { repair_boundaries }),
        )
        .collect();
    // One source body in one round: its topology repair, then its check.
    let attempt = |round: Round, source_index: usize| -> Result<TolerantSolid, NativeBodyFailure> {
        let original = &source.bodies[source_index];
        // A recorded pass reads the written source parameters, as the first
        // measured pass of the main-first policy does.
        fits.repair_boundaries = matches!(
            round,
            Round::Measured {
                repair_boundaries: true
            }
        );
        let recorded = match round {
            Round::Recorded => recorded_mm,
            Round::Measured { .. } => None,
        };
        let checkpoint = evidence.checkpoint();
        let before = catallax_numbers::current_operation_work().used;
        fits.stage = NativeStage::BodyBuild;
        if catallax_numbers::current_operation_work().exhausted {
            return Err(NativeBodyFailure {
                stage: NativeStage::RequestLimit,
                error: ImportRefusal::ResourceLimit.into(),
                work_used: 0,
            });
        }
        let candidate = (|| {
            let repaired = healing::repair_native_topology(
                original,
                &parsed,
                &source.mm_per_unit,
                &angles,
                &mut fits.curves,
                &mut fits.stage,
                dec,
            )?;
            if recorded.is_some() && !repaired.changes.is_empty() {
                // The record is of the file's own entities. A repaired
                // topology is not theirs, so the body is measured.
                return Err(invalid(
                    "import.tolerant-broken-topology",
                    "the file's topology needed a repair, so its recorded tolerance does not describe the repaired body",
                ));
            }
            evidence.repairs.extend(repaired.changes);
            evidence.topology_checks.extend(repaired.checks);
            evidence.topology_measurements.extend(repaired.measurements);
            for issue in repaired.issues {
                if issue.code.starts_with("shell.")
                    || issue.code.starts_with("wire.open")
                    || issue.code == "wire.empty-boundary"
                {
                    return Err(invalid(
                        "import.tolerant-non-manifold",
                        format!("{} {}: {}", issue.body, issue.scope, issue.code),
                    ));
                }
                evidence.warnings.push(NativeWarning {
                    body: Some(issue.body),
                    scope: issue.scope,
                    code: "topology-repair-limitation",
                    detail: issue.code,
                });
            }
            let mut body = repaired.candidate;
            let counts = match recorded {
                Some(recorded) => {
                    fit::check_recorded_body(&mut fits, &mut body, &mut evidence, recorded, dec)?
                }
                None => fit::check_body(&mut fits, &mut body, &mut evidence, dec)?,
            };
            dec.check_deadline()?;
            if catallax_numbers::current_operation_work().exhausted {
                return Err(ImportRefusal::ResourceLimit.into());
            }
            Ok(TolerantSolid {
                body,
                counts,
                source_index,
                recorded_mm: recorded.cloned(),
            })
        })();
        candidate.map_err(|error| {
            evidence.rollback(checkpoint);
            if recorded.is_some() && !super::native::control(&error) {
                // The measured passes follow. This says why the file's
                // record was not used.
                evidence.warnings.push(NativeWarning {
                    body: Some(original.name.clone()),
                    scope: "body".into(),
                    code: "recorded-tolerance-unused",
                    detail: format!(
                        "the file records a tolerance, but this reader could not take the body as written at {:?} ({error:?}); it is measured instead",
                        fits.stage
                    ),
                });
            }
            NativeBodyFailure {
                stage: fits.stage,
                error,
                work_used: catallax_numbers::current_operation_work()
                    .used
                    .saturating_sub(before),
            }
        })
    };
    let (mut bodies, body_errors) =
        over_rounds(source.bodies.len(), &rounds, incomplete, attempt, || {
            catallax_numbers::current_operation_work().exhausted
        });
    bodies.sort_unstable_by_key(TolerantSolid::source_index);
    evidence.warnings.push(NativeWarning {
        body: None,
        scope: "document".into(),
        code: "sampled-authority",
        detail: "Native edge/surface fits are sampled. Global embedding, unsampled self-intersection and a complete source-to-result deviation bound are not claimed.".into(),
    });
    let carriers = fits.into_carriers();
    Ok(TolerantImport {
        parsed,
        source,
        bodies,
        evidence,
        carriers,
        body_errors,
    })
}

#[cfg(test)]
#[path = "../../tests/common/topped_tube.rs"]
mod topped_tube;

#[cfg(test)]
mod round_tests {
    use super::*;
    use catallax_topology::numerical::gap::{EDGE_WORK, WORK_CAP};

    const AS_WRITTEN: Round = Round::Measured {
        repair_boundaries: false,
    };
    const REPAIRED: Round = Round::Measured {
        repair_boundaries: true,
    };

    fn refused(code: &'static str) -> NativeBodyFailure {
        NativeBodyFailure {
            stage: NativeStage::CurveFit,
            error: Error::Numerical(Refusal::new(code, "a test refusal")),
            work_used: 0,
        }
    }

    #[test]
    fn a_body_an_edge_cap_refused_is_not_measured_again_and_the_others_go_on() {
        // Three bodies, as written and then with their boundaries repaired.
        // An edge of the first spends its work allowance while the request
        // has allowance to spare. The second fails on its geometry as written
        // and passes once repaired. The third passes as written.
        let mut tried = Vec::new();
        let (admitted, refusals) = over_rounds(
            3,
            &[AS_WRITTEN, REPAIRED],
            None,
            |round, body| {
                let repaired = round == REPAIRED;
                tried.push((repaired, body));
                match (body, repaired) {
                    (0, _) => Err(refused(WORK_CAP)),
                    (1, false) => Err(refused("import.tolerant-fit")),
                    _ => Ok(body),
                }
            },
            || false,
        );
        // The repair round does not give the capped body's edges a second
        // allowance. A geometry failure is still repaired.
        assert_eq!(tried, [(false, 0), (false, 1), (false, 2), (true, 1)]);
        assert_eq!(admitted, [2, 1]);
        let codes: Vec<_> = refusals
            .iter()
            .map(
                |refusal| match refusal.as_ref().map(|failure| &failure.error) {
                    Some(Error::Numerical(refusal)) => Some(refusal.code),
                    _ => None,
                },
            )
            .collect();
        assert_eq!(codes, [Some(WORK_CAP), None, None]);
    }

    #[test]
    fn a_spent_request_and_a_body_that_is_not_needed_are_as_before() {
        // A body that is not needed is never tried, and no round starts once
        // the request's allowance is spent: a failure of the first round then
        // stands, as it did.
        let mut tried = Vec::new();
        let (admitted, refusals) = over_rounds(
            2,
            &[AS_WRITTEN, REPAIRED],
            Some(&[false, true]),
            |round, body| {
                tried.push((round == REPAIRED, body));
                Err::<usize, _>(refused("import.tolerant-fit"))
            },
            || true,
        );
        assert_eq!(tried, [(false, 1)]);
        assert!(admitted.is_empty());
        assert!(refusals[0].is_none() && refusals[1].is_some());
    }

    /// `second` appended to `first` as one file of two bodies: the entities
    /// of `second` follow those of `first`, numbered after them.
    fn one_file(first: &str, second: &str) -> Vec<u8> {
        let data = |text: &str| -> (usize, usize) {
            let start = text.find("DATA;").expect("a data section") + "DATA;".len();
            let end = text.rfind("ENDSEC;").expect("the end of the data section");
            (start, end)
        };
        let offset = first
            .lines()
            .filter_map(|line| {
                line.strip_prefix('#')?
                    .split('=')
                    .next()?
                    .parse::<u64>()
                    .ok()
            })
            .max()
            .expect("an entity");
        let (start, end) = data(second);
        let mut renumbered = String::new();
        let mut rest = &second[start..end];
        while let Some(at) = rest.find('#') {
            renumbered.push_str(&rest[..=at]);
            rest = &rest[at + 1..];
            let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
            if digits > 0 {
                let id: u64 = rest[..digits].parse().expect("an entity number");
                renumbered.push_str(&(id + offset).to_string());
                rest = &rest[digits..];
            }
        }
        renumbered.push_str(rest);
        let (_, end) = data(first);
        format!(
            "{}{}{}",
            &first[..end],
            renumbered.trim_start(),
            &first[end..]
        )
        .into_bytes()
    }

    /// Section 4.2 step 8: "A hard per-edge work cap refuses by name." A file
    /// of two bodies, read as an import through the measured rounds a reopen
    /// runs. Both are the square tube with a spline top, whose u/v curves the
    /// gap rule reads. The first has one u/v curve at another speed, so that
    /// edge counts more than any edge of the second.
    ///
    /// Under the checker's cap both bodies are measured and admitted, which
    /// gives each edge's count. With one unit less than the costliest edge
    /// counts, the first body is refused by the cap's name, in the stage that
    /// measures its faces, after its edges were read under an allowance they
    /// could spend. The second body is measured and admitted as before, and
    /// carries no refusal.
    #[test]
    fn a_capped_edge_refuses_its_own_body_after_it_was_measured_and_no_other() {
        use std::collections::BTreeMap;
        let plain = String::from_utf8(topped_tube::unstamped(
            &topped_tube::written(&topped_tube::spline_topped_tube()).0,
        ))
        .unwrap();
        let stretched = topped_tube::with_one_uv_curve_at_another_speed(&plain);
        let bytes = one_file(&stretched, &plain);
        let sha = crate::sha256_hex(&bytes);
        let admit_at = |edge_work: u64| {
            catallax_numbers::with_operation_budget(100_000_000_000, || {
                let mut dec = Decider::default();
                let (parsed, source, angles) =
                    prepare(&bytes, &sha, ToleranceOptions::default(), &mut dec).unwrap();
                let needed = vec![true; source.bodies.len()];
                admit_under(
                    edge_work,
                    parsed,
                    source,
                    angles,
                    AdmissionPolicy::MainFirst(&needed),
                    None,
                    &mut dec,
                )
                .unwrap()
            })
        };

        let (whole, _) = admit_at(EDGE_WORK);
        let names: Vec<&str> = whole
            .source
            .bodies
            .iter()
            .map(|body| body.name.as_str())
            .collect();
        assert_eq!(names.len(), 2, "{names:?}");
        assert_eq!(whole.bodies.len(), 2, "{:?}", whole.body_errors);
        let mut counted: BTreeMap<(&str, &str), u64> = BTreeMap::new();
        for row in &whole.evidence.import_gaps {
            *counted
                .entry((row.body.as_str(), row.edge.as_str()))
                .or_default() += row.work;
        }
        let costliest = |body: &str| {
            counted
                .iter()
                .filter(|((of, _), _)| *of == body)
                .map(|(_, work)| *work)
                .max()
                .expect("an edge the gap rule read")
        };
        let (first, second) = (costliest(names[0]), costliest(names[1]));
        assert!(second > 0 && first > second, "{first} and {second}");

        let (capped, request) = admit_at(first - 1);
        assert!(!request.exhausted, "{request:?}");
        let failure = capped.body_errors[0]
            .as_ref()
            .expect("the body with the capped edge");
        assert!(
            matches!(&failure.error, Error::Numerical(refusal) if refusal.code == WORK_CAP),
            "{failure:?}"
        );
        // The stage that reads each face's edge uses through the gap rule,
        // past the body's decoding and its topology, with work spent on the
        // body. Every edge had an allowance it could spend, `first - 1`, and
        // only a reading under an edge's allowance gives the cap's name.
        assert_eq!(failure.stage, NativeStage::SurfaceFit);
        assert!(failure.work_used > 0, "{failure:?}");
        assert!(capped.body_errors[1].is_none(), "{:?}", capped.body_errors);
        let admitted: Vec<usize> = capped
            .bodies
            .iter()
            .map(TolerantSolid::source_index)
            .collect();
        assert_eq!(admitted, [1]);
        // The admitted body's edges were read as under the checker's cap.
        let read: Vec<_> = capped
            .evidence
            .import_gaps
            .iter()
            .map(|row| (row.body.as_str(), row.edge.as_str(), row.work))
            .collect();
        let before: Vec<_> = whole
            .evidence
            .import_gaps
            .iter()
            .filter(|row| row.body == names[1])
            .map(|row| (row.body.as_str(), row.edge.as_str(), row.work))
            .collect();
        assert!(!read.is_empty());
        assert_eq!(read, before);
        println!(
            "two tubes: the costliest edge counts {first} units in the first body and {second} in the second; under {} the first body is refused after {} units of its own",
            first - 1,
            failure.work_used
        );
    }
}
