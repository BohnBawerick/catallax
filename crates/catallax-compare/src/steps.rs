//! Independent prefix builds, without display or a claim of geometric agreement.

use catallax_features::{BuildError, ImportFiles};
use catallax_numbers::{Decider, Refusal, with_operation_budget, with_operation_log};
use catallax_recipe::{Recipe, StepId, ValueKind};

/// One kernel's construction outcome for a complete recipe prefix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Construction {
    /// A sketch, path or scan is a reference definition rather than a solid.
    Reference,
    /// The kernel constructed the requested body independently.
    Built {
        /// Number of result face candidates, before pointwise trim queries.
        faces: usize,
        /// Whether the kernel proved the empty solid.
        empty: bool,
    },
    /// The first failure, which can precede the requested step.
    Refused {
        /// The first operation that prevented this prefix from building.
        blocked_at: Option<String>,
        /// Original typed kernel refusal.
        refusal: Refusal,
    },
}

/// Construction and deterministic work for one kernel and prefix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KernelBuild {
    /// Construction outcome, independent of the other kernel.
    pub construction: Construction,
    /// Work used by this complete prefix build.
    pub work: u64,
    /// Unchanged per-build allowance supplied by the caller.
    pub limit: u64,
}

/// One row per authored recipe step, in recipe order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IndependentStep {
    /// Authored step id.
    pub step: StepId,
    /// Operation schema name.
    pub operation: String,
    /// B-rep construction, never borrowed by the field.
    pub brep: KernelBuild,
    /// Field construction from recipe definitions, never a B-rep conversion.
    pub field: KernelBuild,
}

/// Builds every body-valued prefix in both kernels with a fresh thread, decider
/// and the same per-build work allowance. Call outside an operation budget.
/// Reference rows report no body; they are validated by the recipe and by each
/// subsequent body build. This is construction coverage, not sampled agreement,
/// display, measurement, export or a whole-boundary distance certificate.
pub fn independent_steps(
    recipe: &Recipe,
    files: &ImportFiles<'_>,
    max_work: u64,
) -> Result<Vec<IndependentStep>, Refusal> {
    recipe
        .check()
        .map_err(|e| Refusal::new("compare.recipe", e.to_string()))?;
    let mut rows = Vec::with_capacity(recipe.steps.len());
    for (index, step) in recipe.steps.iter().enumerate() {
        let reference = || KernelBuild {
            construction: Construction::Reference,
            work: 0,
            limit: max_work,
        };
        let (brep, field) = if step.op.kind() != ValueKind::Body {
            (reference(), reference())
        } else {
            let mut prefix = recipe.clone();
            prefix.steps.truncate(index + 1);
            prefix.result = step.id.clone();
            prefix
                .materials
                .retain(|m| prefix.steps.iter().any(|s| s.id == m.region));
            let build = |field: bool| {
                std::thread::scope(|scope| {
                    std::thread::Builder::new().stack_size(8 * 1024 * 1024).spawn_scoped(scope, || {
                        let ((construction, work), log) = with_operation_log(|| with_operation_budget(max_work, || {
                            let mut dec = Decider::default();
                            if field {
                                catallax_field::Plan::compile_with_files_checked(&prefix, files, &mut dec)
                                    .map(|p| Construction::Built {
                                        faces: p.faces().iter().enumerate().filter(|(i,_)| p.result_face(*i).is_some()).count(),
                                        empty: p.is_proven_empty(),
                                    }).map_err(|r| (None, r))
                            } else {
                                catallax_features::evaluate_with_files(&prefix, files, &mut dec)
                                    .and_then(|b| {
                                        if let Some(s) = b.stand_ins.first() {
                                            return Err(BuildError::Geometry { step: s.step.clone(), refusal: Refusal::new("compare.envelope-only", "region stand-in is not an independent construction of the field feature") });
                                        }
                                        Ok(Construction::Built { faces: b.model.faces_of_body(b.result).len(), empty: b.model.body(b.result).kind == catallax_topology::BodyKind::Empty })
                                    }).map_err(|e| match e {
                                        BuildError::Geometry { step, refusal } => (Some(step.to_string()), refusal),
                                        BuildError::Recipe(e) => (None, Refusal::new("compare.recipe", e.to_string())),
                                    })
                            }
                        }));
                        let construction = match construction {
                            Ok(c) if !work.exhausted => c,
                            result => {
                                let (blocked_at, refusal) = result.err().unwrap_or_else(|| (None, Refusal::new("compare.work-exhausted", format!("prefix exhausted {max_work} work units"))));
                                Construction::Refused { blocked_at: blocked_at.or_else(|| log.iter().find(|e| e.refusal.is_some()).map(|e| e.id.clone())), refusal }
                            }
                        };
                        KernelBuild { construction, work: work.used, limit: max_work }
                    }).expect("start independent construction thread").join().expect("independent construction panicked")
                })
            };
            (build(false), build(true))
        };
        rows.push(IndependentStep {
            step: step.id.clone(),
            operation: step.op.name().into(),
            brep,
            field,
        });
    }
    Ok(rows)
}
