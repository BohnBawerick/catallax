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
    /// Measured physical fit residual in millimetres.
    pub residual_mm: f64,
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
}

impl TolerantEvidence {
    fn checkpoint(&self) -> [usize; 8] {
        [
            self.tolerances.len(),
            self.sampled_checks.len(),
            self.repairs.len(),
            self.topology_checks.len(),
            self.topology_measurements.len(),
            self.warnings.len(),
            self.native_trims.len(),
            self.native_domains.len(),
        ]
    }
    fn rollback(&mut self, lengths: [usize; 8]) {
        self.tolerances.truncate(lengths[0]);
        self.sampled_checks.truncate(lengths[1]);
        self.repairs.truncate(lengths[2]);
        self.topology_checks.truncate(lengths[3]);
        self.topology_measurements.truncate(lengths[4]);
        self.warnings.truncate(lengths[5]);
        self.native_trims.truncate(lengths[6]);
        self.native_domains.truncate(lengths[7]);
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
}
impl TolerantSolid {
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
        let bytes = self.parsed.write_tolerant_step(
            &built.changes,
            self.parsed.schemas()?,
            &metadata,
            None,
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
            &self.parsed,
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
        let bytes = self.parsed.write_tolerant_step(
            &built.changes,
            self.parsed.schemas()?,
            &metadata,
            Some(&[root]),
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

fn source_angles(parsed: &ParsedStep, dec: &mut Decider) -> Result<AngleUnit, Error> {
    let mut selected = None;
    for entity in parsed.entities() {
        if parsed
            .entity_types(entity)?
            .contains(&"GLOBAL_UNIT_ASSIGNED_CONTEXT")
            && parsed.representation_context(entity)?.dimension == Some(3)
        {
            let units = parsed.context_units(entity, dec)?;
            let angles = match units.plane_angle.unwrap_or(AngleUnit::Radians) {
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
    }
    Ok(selected.unwrap_or(AngleUnit::Radians))
}

/// Admit a separately labelled tolerant document without whole-domain proof gates.
/// Request deadlines and externally installed budgets remain effective.
pub fn import_step(
    bytes: &[u8],
    expected_sha256: &str,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<TolerantImport, Error> {
    let (parsed, source, angles) = prepare(bytes, expected_sha256, options, dec)?;
    let imported = admit(parsed, source, angles, dec)?;
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
    let angles = source_angles(&parsed, dec)?;
    let (scale, declarations) = parsed.native_inventory(dec)?;
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

pub(crate) fn admit(
    parsed: Arc<ParsedStep>,
    source: Arc<TolerantFile>,
    angles: AngleUnit,
    dec: &mut Decider,
) -> Result<TolerantImport, Error> {
    let mut fits = fit::Context::new(
        &parsed,
        &source.mm_per_unit,
        angles.clone(),
        Default::default(),
        source.options.maximum_mm.to_f64_nearest(),
    );
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
    };
    let mut bodies = Vec::with_capacity(source.bodies.len());
    let mut body_errors = vec![None; source.bodies.len()];
    for (source_index, original) in source.bodies.iter().enumerate() {
        let checkpoint = evidence.checkpoint();
        let before = catallax_numbers::current_operation_work().used;
        fits.stage = NativeStage::BodyBuild;
        if catallax_numbers::current_operation_work().exhausted {
            body_errors[source_index] = Some(NativeBodyFailure {
                stage: NativeStage::RequestLimit,
                error: ImportRefusal::ResourceLimit.into(),
                work_used: 0,
            });
            continue;
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
            let counts = fit::check_body(&mut fits, &mut body, &mut evidence, dec)?;
            dec.check_deadline()?;
            if catallax_numbers::current_operation_work().exhausted {
                return Err(ImportRefusal::ResourceLimit.into());
            }
            Ok(TolerantSolid {
                body,
                counts,
                source_index,
            })
        })();
        match candidate {
            Ok(body) => bodies.push(body),
            Err(error) => {
                evidence.rollback(checkpoint);
                body_errors[source_index] = Some(NativeBodyFailure {
                    stage: fits.stage,
                    error,
                    work_used: catallax_numbers::current_operation_work()
                        .used
                        .saturating_sub(before),
                });
            }
        }
    }
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
