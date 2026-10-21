//! Display meshes derived from exact B-rep definitions.
//!
//! [`mesh_body`] samples each topological edge once. Faces reference the same
//! vertex indices on either side of that edge. Periodic chart copies close
//! seams without adding topological edges. All traversal and refinement orders
//! are fixed; the implementation has no thread-local or process-global state.

use std::{io::Write, ops::Range};

use catallax_numbers::{Answered, Decider, Rational, Refusal};
use catallax_topology::{BodyId, EdgeId, FaceId, Model};

mod chart;
mod edges;
mod triangulate;

/// Display meshing controls. Lengths use the model's millimetres.
#[derive(Clone, Debug)]
pub struct MeshOptions {
    /// Maximum surface-to-chord distance, including coordinate rounding.
    pub chord_tolerance: Rational,
    /// Deterministic limit on generated vertices, including temporary chart copies.
    pub max_vertices: usize,
}

impl Default for MeshOptions {
    fn default() -> Self {
        Self {
            chord_tolerance: Rational::from_ratio(1, 100).unwrap(),
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
    /// Proven maximum deviation from the analytic surface over each triangle's
    /// parameter triangle. Includes endpoint replacement and binary64 rounding.
    pub max_chord_error: Rational,
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
