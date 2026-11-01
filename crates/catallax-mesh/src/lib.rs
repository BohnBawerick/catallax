//! Display meshes derived from exact B-rep definitions.
//!
//! [`mesh_body`] samples each topological edge once. Faces reference the same
//! vertex indices on either side of that edge. Periodic chart copies close
//! seams without adding topological edges. All traversal and refinement orders
//! are fixed; the implementation has no thread-local or process-global state.
//!
//! Checked spline prisms use one adaptive paired ring grid for both planar caps
//! and the ruled side. Whole-span chord bounds include coordinate rounding, and
//! interval orientation checks preserve outward triangles. These bounds concern
//! the represented spline model; its error to an original feature definition is
//! a separate body certificate, not an enlarged display tolerance.
//!
//! Spline faces report certified rounded-facet aspect bounds and unmet requested
//! aspect controls explicitly. A requested edge-length limit is checked against
//! every returned triangle and ring segment; a grid requiring further interior
//! refinement refuses. Whole-spline angular controls remain an explicit coverage
//! refusal rather than an ignored option or an unproved normal bound.

use std::{io::Write, ops::Range};

use catallax_numbers::{Answered, Decider, Rational, Refusal};
use catallax_topology::{BodyId, EdgeId, FaceId, Model};

mod bounded;
mod chart;
mod controls;
mod edges;
mod evaluate;
mod junction;
mod periodic;
mod planar;
/// Unproved preview sampling, separate from certified display and export meshes.
pub mod preview;
pub use preview::{PreviewEdge, PreviewFace, PreviewMesh, PreviewOptions, preview_mesh_body};
mod ruled;
mod spline_prism;
mod triangulate;

/// Display meshing controls. Lengths use the model's millimetres.
#[derive(Clone, Debug)]
pub struct MeshOptions {
    /// Maximum surface-to-chord distance, including coordinate rounding.
    pub chord_tolerance: Rational,
    /// Maximum angle in radians between a facet normal and every analytic
    /// normal over its parameter triangle. Includes coordinate rounding.
    pub angular_deviation: Option<Rational>,
    /// Maximum length of every triangle and boundary-polyline segment.
    pub max_edge_length: Option<Rational>,
    /// Requested longest-edge/shortest-altitude ratio. Reported as unavailable
    /// when triangulation cannot certify it; achieved ratios are always reported.
    pub max_aspect_ratio: Option<Rational>,
    /// Deterministic limit on generated vertices, including temporary chart copies.
    pub max_vertices: usize,
}

impl Default for MeshOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: Rational::from_ratio(1, 100).unwrap(),
            angular_deviation: None,
            max_edge_length: None,
            max_aspect_ratio: None,
            max_vertices: 1_000_000,
        }
    }
}

/// One triangle and the B-rep face to use for picking.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Triangle {
    /// Indices into [`Mesh::positions`], wound toward the face's material-outward normal.
    pub vertices: [usize; 3],
    /// Stable within the input model, without needing a recipe naming layer.
    pub face: FaceId,
}

/// One face's contiguous part of a body mesh.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FaceMesh {
    /// Input B-rep face.
    pub face: FaceId,
    /// Range in [`Mesh::triangles`].
    pub triangles: Range<usize>,
    /// Proven maximum deviation from the represented surface over each triangle's
    /// complete parameter domain. Includes endpoint replacement and binary64 rounding.
    pub max_chord_error: Rational,
    /// Proven upper bound on longest edge divided by shortest altitude.
    /// None when coordinate rounding collapses a triangle to zero area.
    pub max_aspect_ratio: Option<Rational>,
    /// Proven angular bound in radians, when requested and supported.
    pub max_angular_deviation: Option<Rational>,
}

/// The single sampled polyline used by both faces of a B-rep edge.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EdgeMesh {
    /// Input B-rep edge.
    pub edge: EdgeId,
    /// Mesh vertex indices in the edge's direction. A closed edge repeats its first index.
    pub vertices: Vec<usize>,
    /// Proven uniform distance from the exact edge to its polyline.
    pub max_chord_error: Rational,
}

/// A requested display control that was not certified.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnavailableControl {
    /// Face on which the control is unavailable.
    pub face: FaceId,
    /// Name of the requested control.
    pub control: MeshControl,
    /// Why the control could not be certified.
    pub reason: &'static str,
}

/// Optional display controls, distinct from the mandatory chord bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MeshControl {
    /// Surface-to-facet normal deviation.
    AngularDeviation,
    /// Triangle and polyline segment length.
    EdgeLength,
    /// Longest edge divided by shortest altitude.
    AspectRatio,
}
/// A display approximation. It never replaces geometry in the model.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    /// Shared binary64 coordinates, in millimetres.
    pub positions: Vec<[f64; 3]>,
    /// Triangles in body-face order, with canonical cyclic index order.
    pub triangles: Vec<Triangle>,
    /// Per-face picking and accuracy reports.
    pub faces: Vec<FaceMesh>,
    /// Shared B-rep boundary polylines in body-edge order.
    pub edges: Vec<EdgeMesh>,
    /// Requested controls that were not certified. No control is silently ignored.
    pub unavailable_controls: Vec<UnavailableControl>,
}

impl Mesh {
    /// Writes Wavefront OBJ with shared vertices and one `face_<id>` group per face.
    /// Decimal coordinates round-trip to the binary64 values used by the proof.
    pub fn write_obj(&self, out: &mut impl Write) -> std::io::Result<()> {
        writeln!(out, "# catallax display mesh; millimetres")?;
        for p in &self.positions {
            writeln!(out, "v {} {} {}", p[0], p[1], p[2])?;
        }
        for face in &self.faces {
            writeln!(out, "g face_{}", face.face.index())?;
            writeln!(out, "# proven_chord_error {}", face.max_chord_error)?;
            for triangle in &self.triangles[face.triangles.clone()] {
                let [a, b, c] = triangle.vertices.map(|v| v + 1);
                writeln!(out, "f {a} {b} {c}")?;
            }
        }
        Ok(())
    }
}

/// Builds a display mesh without changing the B-rep or the judge/STEP contracts.
/// Analytic faces and certified traced edges carry proven bounds, not sampled
/// estimates. A resource/precision failure returns an error, never a mesh whose
/// advertised tolerance was exceeded. Marks on the input body or its traces are
/// preserved in the answer.
pub fn mesh_body(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Answered<Mesh> {
    if options
        .angular_deviation
        .as_ref()
        .is_some_and(|v| *v <= Rational::zero())
        || options
            .max_edge_length
            .as_ref()
            .is_some_and(|v| *v <= Rational::zero())
        || options
            .max_aspect_ratio
            .as_ref()
            .is_some_and(|v| *v <= Rational::zero())
    {
        return Err(Refusal::new(
            "mesh.invalid-options",
            "optional controls must be positive",
        ));
    }
    if options.chord_tolerance <= Rational::zero() || options.max_vertices < 4 {
        return Err(Refusal::new(
            "mesh.invalid-options",
            "chord tolerance must be positive and max_vertices at least four",
        ));
    }
    triangulate::build(model, body, options, dec)
}

fn failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("mesh.precision-or-budget", detail)
}

#[cfg(test)]
mod tests;
