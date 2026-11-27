//! Unproved display meshes. None of these records carries a geometric certificate.

use std::ops::Range;

use catallax_numbers::{Decider, Rational, Refusal, refund_on_refusal};
use catallax_topology::{BodyId, EdgeId, FaceId, Model};

use crate::{Mesh, MeshOptions, Triangle};

pub(crate) mod edges;
mod imported;
mod tolerant;
pub use imported::{StepPreview, preview_step};
pub use tolerant::preview_tolerant_body;
pub(crate) mod strips;

/// Visual sampling targets, not promised deviation bounds.
#[derive(Clone, Debug)]
pub struct PreviewOptions {
    /// Target chord distance in millimetres. No whole-face bound is claimed.
    pub chord_tolerance: Rational,
    /// Target normal angle in radians. No whole-face bound is claimed.
    pub angular_deviation: Option<Rational>,
    /// Limit on generated vertices, including temporary periodic chart copies.
    /// Also caps candidate visits and added edge references during vertex-weld repair.
    pub max_vertices: usize,
}

impl Default for PreviewOptions {
    fn default() -> Self {
        Self::from(&MeshOptions::default())
    }
}

impl From<&MeshOptions> for PreviewOptions {
    fn from(options: &MeshOptions) -> Self {
        Self {
            chord_tolerance: options.chord_tolerance.clone(),
            angular_deviation: options.angular_deviation.clone(),
            max_vertices: options.max_vertices,
        }
    }
}

/// A face's pickable triangles, without accuracy or quality claims.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PreviewFace<F = FaceId> {
    /// Input-model face handle or imported STEP entity name.
    pub face: F,
    /// Contiguous range in [`PreviewMesh::triangles`].
    pub triangles: Range<usize>,
}

/// A pickable B-rep edge, without a deviation claim.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PreviewEdge<E = EdgeId> {
    /// Input-model edge handle or imported STEP entity name.
    pub edge: E,
    /// Polyline in edge direction. Closed edges repeat their first vertex.
    pub vertices: Vec<usize>,
}

/// Explicitly unproved display data. Not an `Answer<Mesh>` or an export mesh.
///
/// B-rep handles resolve through the retained model; source keys through its file. No
/// watertightness, deviation, angular, aspect, or export guarantee is made.
/// There is deliberately no conversion to [`crate::Mesh`] and no export writer.
#[derive(Clone, Debug, PartialEq)]
pub struct PreviewMesh<F = FaceId, E = EdgeId> {
    /// Binary64 millimetre coordinates.
    pub positions: Vec<[f64; 3]>,
    /// Face-labelled triangles following the input's stated surface sense.
    pub triangles: Vec<Triangle<F>>,
    /// Per-face picking ranges, with no bounds.
    pub faces: Vec<PreviewFace<F>>,
    /// Shared topological edge polylines, with no bounds.
    pub edges: Vec<PreviewEdge<E>>,
}

impl<F, E> Default for PreviewMesh<F, E> {
    fn default() -> Self {
        Self {
            positions: Vec::new(),
            triangles: Vec::new(),
            faces: Vec::new(),
            edges: Vec::new(),
        }
    }
}

impl PreviewMesh {
    pub(crate) fn from_mesh(mesh: Mesh) -> Self {
        Self {
            positions: mesh.positions,
            triangles: mesh.triangles,
            faces: mesh
                .faces
                .into_iter()
                .map(|face| PreviewFace {
                    face: face.face,
                    triangles: face.triangles,
                })
                .collect(),
            edges: mesh
                .edges
                .into_iter()
                .map(|edge| PreviewEdge {
                    edge: edge.edge,
                    vertices: edge.vertices,
                })
                .collect(),
        }
    }
}

fn validate(options: &PreviewOptions) -> Result<(), Refusal> {
    if options.chord_tolerance <= Rational::zero()
        || options.max_vertices < 4
        || options
            .angular_deviation
            .as_ref()
            .is_some_and(|v| *v <= Rational::zero())
    {
        return Err(Refusal::new(
            "mesh.preview-invalid-options",
            "visual targets must be positive and max_vertices at least four",
        ));
    }
    Ok(())
}

/// Samples the retained B-rep for display, without certifying its triangles.
/// Its topological edge polylines are sampled once for all incident faces.
///
/// A geometry or resource refusal at the requested detail selects a coarser
/// sampled grid. The call preserves cancellation and the caller's resource
/// limits, and refuses if even the coarsest preview will not fit.
pub fn preview_mesh_body(
    model: &Model,
    body: BodyId,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh, Refusal> {
    #[cfg(feature = "display-capture")]
    let call =
        crate::capture::Call::start("preview_mesh_body", crate::capture::carriers(model, body));
    let result = preview_checked(model, body, options, dec);
    #[cfg(feature = "display-capture")]
    call.finish(options, &result);
    result
}

fn preview_checked(
    model: &Model,
    body: BodyId,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh, Refusal> {
    validate(options)?;
    let direct = refund_on_refusal(|| {
        if let Some(preview) = crate::triangulate::shared_spline_preview(model, body, options, dec)
        {
            return preview;
        }
        if crate::junction::present(model, body) {
            return crate::junction::preview(model, body, options, dec);
        }
        crate::triangulate::build_preview(model, body, options, None, dec)
    });
    direct.or_else(|direct_refusal| {
        coarsen(options, dec, |divisions, tolerance, dec| {
            crate::float::grid::build(model, body, options, divisions, tolerance, dec)
        })
        .map_err(|coarse_refusal| {
            // A fallback that cannot draw the body is not a new public
            // outcome. Cancellation or real exhaustion still takes priority.
            if coarse_refusal.code == "numbers.deadline"
                || catallax_numbers::current_operation_work().exhausted
            {
                coarse_refusal
            } else {
                direct_refusal
            }
        })
    })
}

/// Draws the retained numerical fit without claiming a bound to its source.
/// Face and edge keys belong to the numerical definition.
pub fn preview_numerical_body(
    body: &catallax_topology::numerical::NumericalBody,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh<String, String>, Refusal> {
    validate(options)?;
    coarsen(options, dec, |divisions, tolerance, dec| {
        let edges = body.definition().edges.len();
        dec.charge_traversal(edges as u64)?;
        let description = crate::float::describe_numerical_body(
            body,
            &crate::float::DisplayOptions {
                chord_tolerance: tolerance,
                max_edge_points: options.max_vertices.saturating_add(edges.saturating_mul(2)),
            },
        )?;
        crate::float::grid::from_description(description, options, divisions, dec)
    })
}

fn coarsen<T>(
    options: &PreviewOptions,
    dec: &mut Decider,
    mut build: impl FnMut(u32, f64, &mut Decider) -> Result<T, Refusal>,
) -> Result<T, Refusal> {
    let chord = options
        .chord_tolerance
        .to_f64_nearest()
        .clamp(f64::MIN_POSITIVE, f64::MAX);
    for (divisions, scale) in [(16, 1.0), (8, 4.0), (4, 16.0), (2, 64.0)] {
        dec.check_deadline()?;
        if let Ok(preview) =
            refund_on_refusal(|| build(divisions, (chord * scale).min(f64::MAX), dec))
        {
            return Ok(preview);
        }
    }
    dec.check_deadline()?;
    build(1, f64::MAX, dec)
}
