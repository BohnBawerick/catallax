//! The manifest, formats `exact-brep-manifest/2` and `exact-brep-manifest/3`.

use std::collections::BTreeMap;

/// Answer status.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, serde::Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Status {
    /// The kernel answered (exactly or with marks).
    Answered,
    /// A named case met its preset level without an exact proof.
    Bounded,
    /// The kernel refused with a named reason.
    Refused,
    /// The question file was invalid.
    InvalidQuestion,
    /// The kernel failed internally.
    InternalError,
}

/// Kernel identity.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct KernelInfo {
    /// Always `exact-brep`.
    pub name: String,
    /// Crate version.
    pub version: String,
    /// Source commit.
    pub commit: String,
}

/// A proven bound as outward-rounded decimal strings.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct Bound {
    /// Lower end.
    pub lo: String,
    /// Upper end.
    pub hi: String,
}

/// Counts of the answer body.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct ResultSummary {
    /// `solid` or `empty`.
    pub kind: String,
    /// Lumps.
    pub lumps: usize,
    /// Shells.
    pub shells: usize,
    /// Faces.
    pub faces: usize,
    /// Edges.
    pub edges: usize,
    /// Ring edges.
    pub ring_edges: usize,
    /// Vertices.
    pub vertices: usize,
    /// Non-zero face counts by surface kind name.
    pub faces_by_kind: BTreeMap<String, usize>,
    /// Non-zero edge counts by curve kind name.
    pub edges_by_kind: BTreeMap<String, usize>,
}

/// Requested properties with proven bounds.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct Properties {
    /// Volume, if requested.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub volume: Option<Bound>,
    /// Area, if requested.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub area: Option<Bound>,
    /// Centre of mass, if requested: `Some(None)` for an empty body.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub centroid: Option<Option<[Bound; 3]>>,
    /// Centroidal unit-density inertia, if requested; null for an empty body.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inertia: Option<Option<[[Bound; 3]; 3]>>,
}

/// What the STEP export wrote.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct StepSummary {
    /// Present only for a rule-bound export without an exact proof.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bounded: Option<crate::BoundedSummary>,
    /// File name inside the output directory.
    pub file: String,
    /// Proven upper bound on deviation, in millimetres.
    pub max_deviation: String,
    /// Intersection branches written as splines.
    pub approximated_curves: usize,
    /// Vertices added for STEP.
    pub added_vertices: usize,
    /// Edges added for STEP.
    pub added_edges: usize,
    /// `"required"` when the file meets the seam contract of a
    /// `exact-brep-question/3` question; absent otherwise.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub seams: Option<&'static str>,
}

impl StepSummary {
    /// Summarizes a STEP result while retaining bounded evidence as such.
    /// A refusal produces no successful export summary.
    pub fn from_result(
        result: &catallax_step::StepResult,
        file: impl Into<String>,
        seams: Option<&'static str>,
    ) -> Result<Self, catallax_numbers::Refusal> {
        use catallax_numbers::OperationResult;
        let (report, bounded) = match result {
            OperationResult::Proven(answer) => (&answer.value, None),
            OperationResult::Bounded(answer) => (
                answer.value(),
                Some(crate::BoundedSummary(answer.evidence().clone())),
            ),
            OperationResult::Refused(refusal) => return Err(refusal.clone()),
        };
        Ok(Self {
            bounded,
            file: file.into(),
            max_deviation: report.max_deviation.to_decimal_bounds(40).1,
            approximated_curves: report.approximated_curves,
            added_vertices: report.added_vertices,
            added_edges: report.added_edges,
            seams,
        })
    }
}

/// One mark on the answer.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct MarkEntry {
    /// Stable code.
    pub code: String,
    /// `distance`, `angle` or `value`.
    pub kind: String,
    /// Proven bound, as a decimal string rounded up.
    pub bound: String,
    /// Step that produced it, if any.
    pub step: Option<String>,
    /// Free text.
    pub detail: String,
}

/// The refusal, for refused, invalid or failed answers.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct RefusalEntry {
    /// Stable code.
    pub code: String,
    /// Step that refused, if any.
    pub step: Option<String>,
    /// Boolean stage that broke, if any.
    pub stage: Option<String>,
    /// Free text.
    pub detail: String,
}

/// Per-step measurements.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct StepStats {
    /// Step id.
    pub id: String,
    /// Wall time in milliseconds.
    pub wall_ms: u64,
    /// Decisions made.
    pub decisions: u64,
    /// Decisions that ended not exact.
    pub not_exact: u64,
}

/// Measurements of the whole answer; may differ between runs.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct Stats {
    /// Wall time in milliseconds.
    pub wall_ms: u64,
    /// Decisions by outcome, as in `catallax_numbers::DecisionStats`.
    pub decisions: catallax_numbers::DecisionStats,
    /// Per-step measurements.
    pub per_step: Vec<StepStats>,
}

/// One answer's manifest.
#[derive(Clone, PartialEq, Eq, Debug, serde::Serialize)]
pub struct Manifest {
    /// Present only for a bounded answer, separately from proof-bearing answers.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bounded: Option<crate::BoundedSummary>,
    /// Format tag, written first: [`crate::MANIFEST_FORMAT`], or
    /// [`crate::SEAMED_MANIFEST_FORMAT`] for a seamed question.
    #[serde(skip)]
    pub format: &'static str,
    /// The question id, if it could be read.
    pub question_id: Option<String>,
    /// Kernel identity.
    pub kernel: KernelInfo,
    /// Status.
    pub status: Status,
    /// Answered with no marks.
    pub exact: bool,
    /// Counts of the answer body.
    pub result: Option<ResultSummary>,
    /// Measured properties.
    pub properties: Option<Properties>,
    /// STEP export summary.
    pub step: Option<StepSummary>,
    /// Every mark on the answer.
    pub marks: Vec<MarkEntry>,
    /// The refusal, if any.
    pub refusal: Option<RefusalEntry>,
    /// Measurements.
    pub stats: Stats,
}

impl Manifest {
    /// Labels a whole-result bound separately from an existing proven answer.
    /// The producer must bind the evidence to this manifest's complete result.
    pub fn with_bounded_evidence(mut self, evidence: catallax_numbers::BoundedEvidence) -> Self {
        self.status = Status::Bounded;
        self.exact = false;
        self.refusal = None;
        self.bounded = Some(crate::BoundedSummary(evidence));
        self
    }
}

/// The manifest as JSON text, with keys in the documented order.
pub fn render_manifest(manifest: &Manifest) -> String {
    #[derive(serde::Serialize)]
    struct Wire<'a> {
        format: &'static str,
        #[serde(flatten)]
        manifest: &'a Manifest,
    }
    let mut text = serde_json::to_string_pretty(&Wire {
        format: manifest.format,
        manifest,
    })
    .expect("manifest fields serialize as JSON");
    text.push('\n');
    text
}

impl serde::Serialize for Stats {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut s = serializer.serialize_struct("Stats", 13)?;
        s.serialize_field("wall_ms", &self.wall_ms)?;
        s.serialize_field("decisions", &self.decisions.decisions)?;
        s.serialize_field("filter", &self.decisions.filter)?;
        s.serialize_field("refined", &self.decisions.refined)?;
        s.serialize_field("exact", &self.decisions.exact)?;
        s.serialize_field("coincident_exact", &self.decisions.coincident_exact)?;
        s.serialize_field(
            "coincident_separation",
            &self.decisions.coincident_separation,
        )?;
        s.serialize_field("coincident_fact", &self.decisions.coincident_fact)?;
        s.serialize_field("not_exact", &self.decisions.not_exact)?;
        s.serialize_field("refused", &self.decisions.refused)?;
        s.serialize_field("max_bits", &self.decisions.max_bits)?;
        s.serialize_field("work", &self.decisions.work)?;
        s.serialize_field("per_step", &self.per_step)?;
        s.end()
    }
}
