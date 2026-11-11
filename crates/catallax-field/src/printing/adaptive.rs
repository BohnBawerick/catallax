//! Explicit adaptive print requests, cell certificates and actionable failures.
//!
//! The existing two-way interval certificate is authoritative. Grid-plane
//! propagation keeps unequal cells conforming; no hanging nodes are stitched by
//! a tolerance. Neither a failed candidate nor an alternative is an export.
use super::{Mesh, Options, certified, refusal};
use crate::Plan;
use catallax_numbers::{Decider, Rational, Refusal};

/// Stage active when a request stops.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stage {
    /// Validate request and prepare exact definitions and bounds.
    Preparation,
    /// Construct the conforming candidate and check its topology.
    Candidate,
    /// Bound every mesh patch's distance to the true boundary.
    MeshToBoundary,
    /// Cover the complete true boundary, including unresolved boxes.
    BoundaryToMesh,
    /// Refine the cells implicated by a failed certificate.
    Refinement,
    /// Check output facet representation on processor workers.
    Representation,
    /// The approximation's material topology must match the source definition.
    MaterialTopology,
}
/// Bound attached to an emitted candidate cell after BOTH global checks pass.
/// The conservative global bound is retained per cell, not inferred from pitch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CellCertificate {
    /// Exact candidate-cell extents in millimetres.
    pub bounds: [[Rational; 2]; 3],
    /// Contiguous emitted facets owned by this cell.
    pub triangles: std::ops::Range<usize>,
    /// Mesh patches in this cell are within this distance of the true boundary;
    /// true boundary points in this cell are within it of the complete mesh.
    pub max_deviation: Rational,
}
/// A separate, unexecuted request suggestion. It is NOT estimated to fit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoarserChoice {
    /// User-visible label; never present this choice as a certified output.
    pub label: &'static str,
    /// Twice the original accuracy, retaining both capacity limits.
    pub options: Options,
    /// Why no capacity or calibrated time prediction accompanies the choice.
    pub estimate: &'static str,
}
/// Original refusal and requested accuracy are preserved, including undecided evidence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Failure {
    /// Stage that stopped.
    pub stage: Stage,
    /// Original typed outcome, not a guessed classification.
    pub cause: Refusal,
    /// Unchanged requested controls.
    pub requested: Options,
    /// Decider work performed during this call (not a throughput estimate).
    pub work: u64,
    /// Charged request units including checkpoints, if an allowance is active.
    /// Otherwise explicitly unavailable; no worker creates a fresh allowance.
    pub request_work: Option<u64>,
    /// Separate choice only for a valid positive-accuracy request.
    pub coarser: Option<CoarserChoice>,
}
/// Immutable successful geometry and retained cell evidence.
#[derive(Debug)]
pub struct Output {
    mesh: Mesh,
    cells: Vec<CellCertificate>,
    requested: Options,
}
impl Output {
    /// Certified mesh. Borrowing cannot alter the geometry bound to this report.
    pub fn mesh(&self) -> &Mesh {
        &self.mesh
    }
    /// Emitted cells; strictly classified interior/exterior cells emit no facets.
    pub fn cells(&self) -> &[CellCertificate] {
        &self.cells
    }
    /// Original request; no implicit coarsening occurs.
    pub fn requested(&self) -> &Options {
        &self.requested
    }
}
pub(super) struct Trace {
    pub stage: Stage,
    pub cells: Vec<CellCertificate>,
}
/// Processor count selected by default, respecting the OS availability query.
pub fn default_workers() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}
/// Adaptive certified candidate and complete two-way coverage, with all available
/// cores used for independent output facet validation. Refinement/certification
/// remain ordered and serial; this is not a claim of parallel speed-up for them.
pub fn mesh(plan: &Plan, options: &Options, dec: &mut Decider) -> Result<Output, Box<Failure>> {
    mesh_with_workers(plan, options, default_workers(), dec)
}
/// Same semantic output and deterministic accounting for every positive worker
/// count. Workers perform only prepaid fixed-cost representation checks; no
/// numeric Decider or thread-local operation allowance escapes to a worker.
pub fn mesh_with_workers(
    plan: &Plan,
    options: &Options,
    workers: usize,
    dec: &mut Decider,
) -> Result<Output, Box<Failure>> {
    let before = dec.stats().work;
    let request_before = catallax_numbers::operation_budget_active()
        .then(|| catallax_numbers::current_operation_work().used);
    let mut trace = Trace {
        stage: Stage::Preparation,
        cells: Vec::new(),
    };
    let result = (|| {
        dec.check_deadline()?;
        if options.tolerance <= Rational::zero()
            || options.max_cells == 0
            || options.max_vertices < 4
            || workers == 0
        {
            return Err(refusal(
                "field.print-options",
                "positive accuracy, cell cap and worker count, and at least four vertices required",
            ));
        }
        let mesh = certified::run(plan, options, &Rational::zero(), Some(&mut trace), dec)?;
        trace.stage = Stage::MaterialTopology;
        // Two-sided Hausdorff coverage alone allows a tiny hole to disappear
        // inside a mixed cell. Until a general cell-isotopy proof exists, admit
        // only a single convex block (or proven empty), whose surface is one
        // genus-zero component. Other candidates stay private and refuse.
        let block = matches!(plan.instructions[plan.root], crate::Instruction::Leaf(i)
            if matches!(&plan.leaves[i], crate::geometry::Leaf::Primitive(p) if p.kind == crate::metric::Kind::Block));
        if !plan.is_proven_empty()
            && (!block
                || mesh.component_count() != 1
                || mesh.positions.len() as i64 - (3 * mesh.triangles.len() / 2) as i64
                    + mesh.triangles.len() as i64
                    != 2)
        {
            return Err(refusal(
                "field.print-detail-topology",
                "general mixed-cell material topology is not certified; distance coverage alone cannot authorize losing a small void",
            ));
        }
        trace.stage = Stage::Representation;
        facets(&mesh, workers, dec)?;
        Ok(mesh)
    })();
    match result {
        Ok(mesh) => Ok(Output { mesh, cells: trace.cells, requested: options.clone() }),
        Err(cause) => Err(Box::new(Failure {
            stage: trace.stage, cause, requested: options.clone(), work: dec.stats().work - before,
            request_work: request_before.map(|before| catallax_numbers::current_operation_work().used - before),
            coarser: (options.tolerance > Rational::zero()).then(|| CoarserChoice {
                label: "UNVALIDATED coarser request; requires explicit selection and a fresh proof",
                options: Options { tolerance: options.tolerance.mul_pow2(1), ..options.clone() },
                estimate: "Fit, cells, triangles, peak memory and calibrated time unavailable; coarsening may not resolve this refusal",
            }),
        })),
    }
}

fn facets(mesh: &Mesh, workers: usize, dec: &mut Decider) -> Result<(), Refusal> {
    // Fixed batches and reservations, independent of scheduling and worker count.
    // Every worker is joined before checking cancellation or publishing anything.
    for batch in mesh.triangles.chunks(1024) {
        dec.charge_work(batch.len() as u64)?;
        let width = batch.len().div_ceil(workers.min(batch.len()));
        let result = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            let mut spawn_failure = None;
            for chunk in batch.chunks(width) {
                match std::thread::Builder::new().spawn_scoped(scope, move || {
                    chunk.iter().all(|t| {
                        let p = t.map(|v| mesh.positions[v].map(f64::from));
                        let a: [f64; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
                        let b: [f64; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
                        let normal: [f64; 3] = std::array::from_fn(|i| {
                            a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3]
                        });
                        normal.iter().all(|x| x.is_finite()) && normal.iter().any(|x| *x != 0.)
                    })
                }) {
                    Ok(handle) => handles.push(handle),
                    Err(_) => {
                        spawn_failure = Some(refusal(
                            "field.print-worker",
                            "could not start facet validation worker",
                        ));
                        break;
                    }
                }
            }
            let mut valid = true;
            let mut panic = false;
            for handle in handles {
                match handle.join() {
                    Ok(ok) => valid &= ok,
                    Err(_) => panic = true,
                }
            }
            if let Some(error) = spawn_failure {
                Err(error)
            } else if panic {
                Err(refusal(
                    "field.print-worker",
                    "facet validation worker panicked",
                ))
            } else if !valid {
                Err(refusal(
                    "field.print-coordinate",
                    "a facet has a nonfinite or zero output normal",
                ))
            } else {
                Ok(())
            }
        });
        dec.check_deadline()?;
        result?;
    }
    Ok(())
}
