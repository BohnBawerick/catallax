//! Explicitly unproved field preview, beside the certified display request.
//!
//! The preview reads the same operation definitions and uses the same charts,
//! trim predicates and cell classification as [`super::mesh_plan`]. It differs
//! in two places. Where the certified request would refuse a cell, that cell is
//! counted and left empty. A buried carrier that can never resolve refuses
//! before any cell. No B-rep geometry, sampling or hull fills a gap. And where
//! classification has proved which side of a single-face support is material
//! on an open cell, the preview takes the face and its orientation from that
//! proof and skips the exact witness checks the certified request repeats there.

use super::{Mesh, Options, Triangle, mesh_plan_inner};
use crate::Plan;
use catallax_numbers::{Decider, Rational, Refusal};

/// Visual sampling targets and resource limits for an unproved preview.
#[derive(Clone, Debug)]
pub struct PreviewOptions {
    /// Visual sampling target in millimetres, not a guaranteed deviation bound.
    pub tolerance: Rational,
    /// Maximum distinct binary64 display vertices. Exceeding it refuses.
    pub max_vertices: usize,
    /// Maximum visited parameter cells. Cells beyond it become gaps.
    pub max_cells: usize,
}
impl From<&Options> for PreviewOptions {
    fn from(options: &Options) -> Self {
        Self {
            tolerance: options.tolerance.clone(),
            max_vertices: options.max_vertices,
            max_cells: options.max_cells,
        }
    }
}
impl Default for PreviewOptions {
    fn default() -> Self {
        Self::from(&Options::default())
    }
}

/// Unproved display triangles tagged with field registry handles.
///
/// It carries no certificate and no deviation bound, and has no writer or
/// conversion to the certified [`Mesh`]. It may have gaps and need not sew.
///
/// ```compile_fail
/// use catallax_field::display::{Mesh, PreviewMesh};
/// fn certify(preview: PreviewMesh) -> Mesh { preview.into() }
/// ```
///
/// ```compile_fail
/// use catallax_field::display::PreviewMesh;
/// fn export(preview: PreviewMesh) { preview.write_obj(&mut Vec::new()).unwrap(); }
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct PreviewMesh {
    /// Binary64 millimetre coordinates.
    pub positions: Vec<[f64; 3]>,
    /// Triangles tagged with the field face that owned the emitting cell.
    pub triangles: Vec<Triangle>,
    /// Visited parameter cells, including rejected and deferred cells.
    pub cells: usize,
    /// Cells the certified request would have refused, left empty.
    pub unresolved_cells: usize,
    /// Face queries answered from their cell's proved material sides instead
    /// of exact witness checks. The certified request answers every one exactly.
    pub skipped_checks: usize,
}

#[derive(Default)]
pub(super) struct Tally {
    pub(super) gaps: usize,
    pub(super) skipped_checks: usize,
}

/// Refusals that end a preview: the caller's allowance, deadline or vertex limit.
pub(super) fn is_allowance(refusal: &Refusal) -> bool {
    matches!(
        refusal.code,
        "numbers.operation-budget" | "numbers.deadline"
    ) || (refusal.code == "field.display-budget" && refusal.detail.starts_with("display vertex"))
}
pub(super) fn is_cell_limit(refusal: &Refusal) -> bool {
    refusal.code == "field.display-budget" && refusal.detail.starts_with("refinement limit")
}

/// Builds the field's unproved preview from its own definitions.
///
/// Charge this call to the caller's allowance; it neither resets nor enlarges an
/// active one. An exhausted allowance or deadline, the vertex limit, an
/// unsupported chart, a proved trim obstruction or a lattice refuses without a
/// partial preview.
pub fn preview_plan(
    plan: &Plan,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh, Refusal> {
    let mut tally = Tally::default();
    let mesh: Mesh = mesh_plan_inner(
        plan,
        &Options {
            tolerance: options.tolerance.clone(),
            max_vertices: options.max_vertices,
            max_cells: options.max_cells,
        },
        dec,
        None,
        Some(&mut tally),
    )?;
    Ok(PreviewMesh {
        positions: mesh.positions,
        triangles: mesh.triangles,
        cells: mesh.cells,
        unresolved_cells: tally.gaps,
        skipped_checks: tally.skipped_checks,
    })
}
