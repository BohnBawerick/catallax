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
    /// Original shared edge name.
    pub edge: String,
    /// Traversal against the source EDGE_CURVE.
    pub forward: bool,
    /// Original PCURVE used as a seed, if available.
    pub source_pcurve: Option<u64>,
    /// Actual native correspondence samples.
    pub samples: Vec<UvSample>,
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

/// A tolerant source document, including its original occurrence graph.
#[derive(Debug)]
pub struct TolerantImport {
    parsed: Arc<ParsedStep>,
    source: Arc<TolerantFile>,
    bodies: Vec<TolerantSolid>,
    evidence: TolerantEvidence,
    carriers: fit::Carriers,
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
    pub fn export_body(&self, index: usize, dec: &mut Decider) -> Result<Vec<u8>, Error> {
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
            source_bodies: Some(&[index]),
            source_body_roots: &built.body_roots,
            source_document_evidence: evidence::Scope::new(&self.evidence, Some(body.name())),
            export_checks: &built.checks,
        })
        .map_err(|error| Refusal::new("import.tolerant-evidence", error.to_string()))?;
        Ok(self.parsed.write_tolerant_step(
            &built.changes,
            self.parsed.schemas()?,
            &metadata,
            Some(&[root]),
            dec,
        )?)
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

fn source_angles(
    parsed: &ParsedStep,
    scale: &Rational,
    dec: &mut Decider,
) -> Result<AngleUnit, Error> {
    for entity in parsed.entities() {
        if parsed
            .entity_types(entity)?
            .contains(&"GLOBAL_UNIT_ASSIGNED_CONTEXT")
        {
            let context = parsed.representation_context(entity)?;
            if context.dimension == Some(3) {
                let units = parsed.context_units(entity, dec)?;
                if units.length_mm.as_ref().is_none_or(|value| value == scale)
                    && let Some(angles) = units.plane_angle
                {
                    return Ok(angles);
                }
            }
        }
    }
    Ok(AngleUnit::Radians)
}

/// Admit a separately labelled tolerant document without whole-domain proof gates.
/// Request deadlines and externally installed budgets remain effective.
pub fn import_step(
    bytes: &[u8],
    expected_sha256: &str,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<TolerantImport, Error> {
    let (parsed, source) = prepare(bytes, expected_sha256, options, dec)?;
    admit(parsed, source, dec)
}

pub(crate) fn prepare(
    bytes: &[u8],
    expected_sha256: &str,
    options: ToleranceOptions,
    dec: &mut Decider,
) -> Result<(Arc<ParsedStep>, Arc<TolerantFile>), Error> {
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
    let definition = parsed.native_definition(dec)?;
    let mut source = TolerantFile {
        source_sha256: expected_sha256.into(),
        mm_per_unit: definition.scale,
        options: options.clone(),
        bodies: Vec::with_capacity(definition.bodies.len()),
    };
    for definition in definition.bodies {
        let mut body = TolerantBody::new(&definition.name, definition.solid, &options);
        for shell in definition.shells {
            body.shell(&shell.name, shell.closed, &shell.faces, dec)?;
        }
        source.bodies.push(body);
    }
    Ok((parsed, Arc::new(source)))
}

pub(crate) fn admit(
    parsed: Arc<ParsedStep>,
    source: Arc<TolerantFile>,
    dec: &mut Decider,
) -> Result<TolerantImport, Error> {
    let angles = source_angles(&parsed, &source.mm_per_unit, dec)?;
    let repaired = healing::repair_native_topology(&source, &parsed, &angles, dec)?;
    let mut evidence = TolerantEvidence {
        method: "catallax-tolerant-native-v1",
        source_sha256: source.source_sha256.clone(),
        authority: "tolerant: discrete topology and sampled native fits; no whole-domain certificate",
        tolerances: Vec::new(),
        sampled_checks: Vec::new(),
        repairs: repaired.changes,
        topology_checks: repaired.checks,
        topology_measurements: repaired.measurements,
        warnings: Vec::new(),
        native_trims: Vec::new(),
        native_domains: Vec::new(),
    };
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
    let mut bodies = Vec::with_capacity(repaired.candidate.bodies.len());
    let mut fits = fit::Context::new(&parsed, &source.mm_per_unit, angles, repaired.curves);
    for mut body in repaired.candidate.bodies {
        let counts = fit::check_body(&mut fits, &mut body, &mut evidence, dec)?;
        bodies.push(TolerantSolid { body, counts });
    }
    if bodies.is_empty() {
        return Err(invalid(
            "import.tolerant-broken-topology",
            "source document contains no finite B-rep body",
        ));
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
    })
}
