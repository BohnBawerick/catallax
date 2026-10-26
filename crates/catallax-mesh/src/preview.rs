//! Unproved display meshes. None of these records carries a geometric certificate.

use std::ops::Range;

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_topology::{BodyId, EdgeId, FaceId, Model};

use crate::{MeshOptions, Triangle};

pub(crate) mod edges;
mod spline;
pub(crate) mod strips;

/// Visual sampling targets, not promised deviation bounds.
#[derive(Clone, Debug)]
pub struct PreviewOptions {
    /// Target chord distance in millimetres. No whole-face bound is claimed.
    pub chord_tolerance: Rational,
    /// Target normal angle in radians. No whole-face bound is claimed.
    pub angular_deviation: Option<Rational>,
    /// Limit on generated vertices, including temporary periodic chart copies.
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
pub struct PreviewFace {
    /// Face handle in the retained input model.
    pub face: FaceId,
    /// Contiguous range in [`PreviewMesh::triangles`].
    pub triangles: Range<usize>,
}

/// A pickable B-rep edge, without a deviation claim.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PreviewEdge {
    /// Edge handle in the retained input model.
    pub edge: EdgeId,
    /// Polyline in edge direction. Closed edges repeat their first vertex.
    pub vertices: Vec<usize>,
}

/// Explicitly unproved display data. Not an `Answer<Mesh>` or an export mesh.
///
/// Face and edge handles resolve origins through the retained model. No
/// watertightness, deviation, angular, aspect, or export guarantee is made.
/// There is deliberately no conversion to [`crate::Mesh`] and no export writer.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PreviewMesh {
    /// Binary64 millimetre coordinates.
    pub positions: Vec<[f64; 3]>,
    /// Face-labelled triangles, wound toward the material-outward normal.
    pub triangles: Vec<Triangle>,
    /// Per-face picking ranges, with no bounds.
    pub faces: Vec<PreviewFace>,
    /// Shared topological edge polylines, with no bounds.
    pub edges: Vec<PreviewEdge>,
}

/// Samples the retained B-rep for display, without certifying its triangles.
///
/// This synchronous call uses the caller's decider and active work allowance.
/// It never starts, replaces, or enlarges an allowance and never calls
/// [`crate::mesh_body`]. Failure returns no partial preview.
pub fn preview_mesh_body(
    model: &Model,
    body: BodyId,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh, Refusal> {
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
    if model.faces_of_body(body).iter().any(|&face| {
        matches!(
            model.geom().surface(model.face(face).surface),
            catallax_shapes::Surface::BSpline(_)
        )
    }) {
        return spline::build(model, body, options, dec);
    }
    crate::triangulate::build_preview(model, body, options, dec)
}
