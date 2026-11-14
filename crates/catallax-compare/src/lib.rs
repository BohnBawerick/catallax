//! Two-way sampled display comparison. Agreement is not a global Hausdorff
//! certificate or a pre-export accuracy certificate. Imported definitions also
//! compare exact containment samples and independent volume enclosures.
//! See `docs/archive/PLAN-2026-11.md` for coverage and measured costs.
#![forbid(unsafe_code)]

use catallax_numbers::{Interval, Rational, Refusal};
use catallax_recipe::{Origin, Recipe, StepId};
use catallax_shapes::RPoint3;
use std::sync::{Arc, Mutex};

mod check;
mod import_checks;
mod steps;
pub use steps::{Construction, IndependentStep, KernelBuild, independent_steps};

/// Closed model-space edit region, in millimetres.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    /// Inclusive lower corner.
    pub min: [Rational; 3],
    /// Inclusive upper corner.
    pub max: [Rational; 3],
}
impl Region {
    pub(crate) fn contains(&self, p: &RPoint3, padding: &Rational) -> bool {
        [&p.x, &p.y, &p.z]
            .into_iter()
            .enumerate()
            .all(|(i, x)| *x >= &self.min[i] - padding && *x <= &self.max[i] + padding)
    }
}

/// Deterministic sampling and proof controls.
#[derive(Clone, Debug)]
pub struct Options {
    /// Relative and absolute volume difference limits for imported parts.
    pub import_volume_accuracy: catallax_measure::Accuracy,
    /// Allowed display discrepancy. Default is 0.01 mm.
    pub accuracy: Rational,
    /// Display mesh controls. Default chord error is 0.005 mm.
    pub mesh: catallax_mesh::MeshOptions,
    /// Rational intervals per quarter circle and primitive parameter interval.
    /// This is a sampling density, not a global error certificate.
    pub divisions: u32,
    /// Minimum dot product of unit normals. Default is 999/1000.
    pub normal_cosine: Rational,
    /// Restricts queries to this closed region, padded by mesh error forward.
    /// `None` checks the complete sampled boundary.
    pub region: Option<Region>,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            import_volume_accuracy: Default::default(),
            accuracy: Rational::from_ratio(1, 100).unwrap(),
            mesh: catallax_mesh::MeshOptions {
                chord_tolerance: Rational::from_ratio(1, 200).unwrap(),
                ..Default::default()
            },
            divisions: 4,
            normal_cosine: Rational::from_ratio(999, 1000).unwrap(),
            region: None,
        }
    }
}

/// Which direction supplied an observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pass {
    /// Every selected display mesh vertex, with per-face deviation bounds.
    Forward,
    /// Independent exact field boundary samples projected onto trimmed faces.
    Reverse,
}
/// Why the kernels disagree, or why a sample could not be settled.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Reason {
    /// A certified distance lower bound exceeds the display accuracy.
    Gap,
    /// A certified field origin is absent from the matching B-rep face.
    Origins,
    /// Certified normal enclosures violate the allowed cosine.
    Direction,
    /// A query, projection, build, mesh, or coverage proof is unfinished.
    Unresolved(String),
}
/// A located comparison observation. Empty names mean there was no matching face.
#[derive(Clone, Debug, PartialEq)]
pub struct Finding {
    /// Direction of the check.
    pub pass: Pass,
    /// Exact sample coordinate, not a rounded diagnostic position.
    pub point: RPoint3,
    /// Canonical origin sets of the B-rep candidates.
    pub brep_faces: Vec<Vec<Origin>>,
    /// Canonical field origin set at the witness.
    pub field_faces: Vec<Origin>,
    /// Disagreement or unresolved reason.
    pub reason: Reason,
    /// Proven unsigned gap enclosure, after accounting for mesh deviation.
    /// `None` for naming/direction failures or unfinished distance proofs.
    pub gap: Option<Interval>,
}
/// An operation without implementations in both kernels.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NotComparable {
    /// Recipe step.
    pub step: StepId,
    /// Operation's schema name.
    pub operation: String,
    /// Specific missing implementation.
    pub reason: String,
}
/// Aggregate outcome. Unresolved observations never become agreement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    /// All selected samples and face associations were settled within bounds.
    Agree,
    /// At least one proven discrepancy. Other observations may be unresolved.
    Disagree,
    /// No proven discrepancy, but some required proof or coverage is missing.
    Unresolved,
}
/// Complete report, published only for its still-current revision.
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    /// Requested membership samples and volume bounds for imported definitions.
    pub import_agreement: Option<ImportAgreement>,
    /// Revision supplied by the controller.
    pub revision: u64,
    /// Scope of this report; regional agreement says nothing outside this box.
    pub region: Option<Region>,
    /// Aggregate outcome.
    pub verdict: Verdict,
    /// Number of distinct selected mesh vertices queried.
    pub forward_vertices: usize,
    /// Number of certified surviving field samples checked.
    pub reverse_samples: usize,
    /// Maximum certified distance upper bound over checked samples, including
    /// forward mesh error. `None` means at least one distance is unbounded.
    /// This is not a whole-boundary Hausdorff certificate.
    pub sampled_distance_upper: Option<Rational>,
    /// Proven disagreements and explicit unresolved observations, in fixed order.
    pub findings: Vec<Finding>,
    /// Every unsupported operation, in recipe order.
    pub not_comparable: Vec<NotComparable>,
    /// Original kernel refusal, when a query prevented completion.
    pub refusal: Option<Refusal>,
}

/// One independently classified imported-part sample.
#[derive(Clone, Debug, PartialEq)]
pub struct ImportPoint {
    /// Exact model-space position.
    pub point: RPoint3,
    /// B-rep containment, with boundary coincidence retained.
    pub brep: catallax_field::Location,
    /// Containment from retained field definitions.
    pub field: catallax_field::Location,
}

/// Sampled containment and certified volume enclosures, not a global certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct ImportAgreement {
    /// The fixed 125-point grid, including points outside the source bounds.
    pub points: Vec<ImportPoint>,
    /// B-rep volume enclosure in cubic millimetres.
    pub brep_volume: catallax_measure::Enclosure,
    /// Independently integrated field volume in cubic millimetres.
    pub field_volume: catallax_measure::Enclosure,
    /// Worst difference over both complete volume enclosures.
    pub volume_difference_upper: Rational,
    /// Requested relative/absolute difference bound at the certified magnitude.
    pub volume_difference_allowed: Rational,
    /// Whether the complete difference enclosure fits the requested bound.
    pub volume_agrees: bool,
}

impl ImportAgreement {
    /// Whether all containment samples and the volume difference agree.
    pub fn agrees(&self) -> bool {
        self.volume_agrees && self.points.iter().all(|p| p.brep == p.field)
    }
}

#[derive(Default)]
struct State {
    revision: u64,
    latest: Option<Report>,
}
/// Owns revision invalidation and atomic publication. It does not track app edits
/// or start threads. Start a revision before dispatching work to a background worker.
#[derive(Clone, Default)]
pub struct Controller(Arc<Mutex<State>>);
impl Controller {
    /// Invalidates the previous revision and its published report.
    pub fn begin(&self) -> Revision {
        let mut state = self.0.lock().unwrap();
        state.revision = state
            .revision
            .checked_add(1)
            .expect("revision space exhausted");
        state.latest = None;
        Revision {
            controller: self.clone(),
            id: state.revision,
            #[cfg(test)]
            checkpoint: None,
        }
    }
    /// A complete snapshot for the current revision, or no published result.
    pub fn latest(&self) -> Option<Report> {
        self.0.lock().unwrap().latest.clone()
    }
}
/// Cooperative cancellation token and publication authority for one revision.
#[derive(Clone)]
pub struct Revision {
    controller: Controller,
    id: u64,
    #[cfg(test)]
    checkpoint: Option<Arc<(std::sync::Barrier, std::sync::Barrier)>>,
}
/// Completion status. No partial or stale report is returned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Completion {
    /// A complete report was published while its revision was current.
    Published,
    /// The revision was invalidated; no report was published.
    Cancelled,
}
impl Revision {
    /// Revision identifier, independent of worker scheduling.
    pub fn id(&self) -> u64 {
        self.id
    }
    /// Whether a newer edit has invalidated this work.
    pub fn is_cancelled(&self) -> bool {
        self.controller.0.lock().unwrap().revision != self.id
    }
    /// Invalidates this revision without cancelling a newer one.
    pub fn cancel(&self) {
        let mut state = self.controller.0.lock().unwrap();
        if state.revision == self.id {
            state.revision = state
                .revision
                .checked_add(1)
                .expect("revision space exhausted");
            state.latest = None;
        }
    }
    fn publish(&self, mut report: Report) -> Completion {
        let mut state = self.controller.0.lock().unwrap();
        if state.revision != self.id {
            return Completion::Cancelled;
        }
        // Sampling can retain a budget refusal as an unresolved finding.
        // Keep its request-wide code alongside the partial counts and findings.
        if report.refusal.is_none() && catallax_numbers::current_operation_work().exhausted {
            report.refusal = Some(Refusal::new(
                "numbers.operation-budget",
                "comparison exhausted the enclosing request's work allowance",
            ));
        }
        report.revision = self.id;
        state.latest = Some(report);
        Completion::Published
    }
    /// Builds both kernels and runs the complete two-way display check.
    /// No thread or wall-clock deadline influences a geometric verdict.
    pub fn check(&self, recipe: &Recipe, options: &Options) -> Completion {
        self.check_with_files(recipe, &catallax_features::ImportFiles::new(), options)
    }
    /// Checks both kernels against the same immutable STEP sources. Imports add
    /// exact containment samples and volume enclosures to the boundary check.
    pub fn check_with_files(
        &self,
        recipe: &Recipe,
        files: &catallax_features::ImportFiles<'_>,
        options: &Options,
    ) -> Completion {
        check::run(self, recipe, None, files, options)
    }
    /// Checks an independently built B-rep against a field recipe. Useful for
    /// cached builds and for reproducing kernel discrepancies without altering kernels.
    pub fn check_built(
        &self,
        built: &catallax_features::Built,
        field: &Recipe,
        options: &Options,
    ) -> Completion {
        check::run(
            self,
            field,
            Some(built),
            &catallax_features::ImportFiles::new(),
            options,
        )
    }
}

#[cfg(test)]
mod tests;
