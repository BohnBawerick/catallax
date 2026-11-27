//! Display meshes derived from exact B-rep definitions.
//!
//! [`mesh_body`] samples each topological edge once. Faces reference the same
//! vertex indices on either side of that edge. Periodic chart copies close
//! seams without adding topological edges. Each face is then built as its own
//! result against those samples, and the body places the results in face
//! order. A face reused from an earlier mesh passes the same placement checks
//! as a fresh one. All traversal and refinement orders are fixed; the
//! implementation has no thread-local or process-global state.
//!
//! Spline faces use rational charts split at their knots and checked periodic
//! joins. Whole-span chord bounds include coordinate rounding, and interval
//! orientation checks preserve outward triangles. These bounds concern the
//! represented spline model. Any fit error to the feature definition has its
//! own body certificate.
//!
//! Spline faces report certified rounded-facet aspect bounds and unmet requested
//! aspect controls explicitly. A requested edge-length limit is checked against
//! every returned triangle and ring segment. Angular controls compare rounded
//! facets with normal enclosures over their whole parameter triangles.

use std::{collections::BTreeMap, io::Write, ops::Range};

use catallax_numbers::{Answered, Decider, Rational, Refusal};
use catallax_topology::{BodyId, EdgeId, FaceId, Model};

mod assembly;
mod bounded;
#[cfg(feature = "display-capture")]
pub mod capture;
mod chart;
mod controls;
mod display;
pub use display::{
    DisplayMeasurement, DisplayMesh, SPLINE_DISPLAY_RULE, SplineCarrierDisplay, display_mesh_body,
    display_spline_fillet_carrier, display_spline_fillet_face,
};
mod edges;
mod evaluate;
pub mod float;
mod junction;
mod periodic;
mod planar;
/// Unproved preview sampling, separate from certified display and export meshes.
pub mod preview;
pub use preview::{
    PreviewEdge, PreviewFace, PreviewMesh, PreviewOptions, StepPreview, preview_mesh_body,
    preview_step, preview_tolerant_body,
};
mod rational_atlas;
mod reuse;
pub use reuse::FaceCache;
mod ruled;
mod triangulate;

/// Evidence used when building a face display.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceLevel {
    /// Visual sampling without a certified deviation.
    Preview,
    /// A complete positional bound for display, without an export certificate.
    Bounded,
    /// Certified display geometry built from the exact model.
    Proven,
}

impl FaceLevel {
    const ALL: [Self; 3] = [Self::Preview, Self::Bounded, Self::Proven];
}

/// Display meshing controls. Lengths use the model's millimetres.
#[derive(Clone, Debug, PartialEq, Eq)]
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

/// One triangle with its input face or source-entity picking key.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Triangle<F = FaceId> {
    /// Indices into the owning mesh's positions, following its stated face sense.
    pub vertices: [usize; 3],
    /// Stable input picking key. Certified meshes retain B-rep face handles.
    pub face: F,
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
    mesh_body_impl(model, body, options, None, dec)
}

/// Reuses display triangles from an unmodified previous [`mesh_body`] result.
///
/// Retain the previous model and its meshing options alongside the mesh. Handles
/// may change across edits: faces match by exact definitions and ordered trims,
/// not by handle values. Shared edges are sampled anew before matching, so a
/// changed neighbour's sampling invalidates both incident face meshes.
///
/// Different controls disable reuse. Coupled whole-body spline-prism, drafted
/// cylinder and rational-atlas routes retain their full-build behavior. The
/// previous mesh must belong to `previous_model` and `previous_options`, and
/// must not have been modified. No geometry, topology or tolerance is changed.
pub fn mesh_body_reusing(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    previous_model: &Model,
    previous_options: &MeshOptions,
    previous: &Mesh,
    dec: &mut Decider,
) -> Answered<Mesh> {
    let edit = (options == previous_options).then_some(reuse::Edit {
        previous: Some(reuse::Input {
            model: previous_model,
            geometry: previous,
            source: previous,
            records: None,
        }),
        capture: None,
    });
    mesh_body_impl(model, body, options, edit, dec)
}

/// Builds a checked display and retains per-face provenance for the next edit.
///
/// Supply the previous model, unmodified mesh and its matching cache together.
/// A change in options starts a new cache. A hit still requires the same face
/// definition, effective chord, sampled boundary and earlier shared chords.
/// Expanded rational bodies match their physical cells, not merged picking faces.
pub fn mesh_body_cached(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    previous: Option<(&Model, &Mesh, &FaceCache)>,
    dec: &mut Decider,
) -> Result<(catallax_numbers::Answer<Mesh>, FaceCache), Refusal> {
    let mut capture = reuse::Capture::default();
    let edit = reuse::Edit {
        previous: previous.and_then(|(model, mesh, cache)| cache.input(model, mesh, options)),
        capture: Some(&mut capture),
    };
    let answer = mesh_body_impl(model, body, options, Some(edit), dec)?;
    Ok((answer, FaceCache::new(options.clone(), capture)))
}

/// The point identities and certified topology of a retained B-rep.
///
/// Count queries do not emit render buffers. Emission consumes the schedule
/// instead of tracing edges or triangulating faces again.
pub struct FaceLevelPlan<'a> {
    model: &'a Model,
    expansion: Option<&'a catallax_topology::ChartExpansion>,
    schedule: triangulate::MeshSchedule,
    counts: [BTreeMap<FaceId, usize>; 3],
}

impl FaceLevelPlan<'_> {
    /// Number of distinct indices referenced by this logical face at `level`.
    /// Shared edge samples count once within the face, even across chart cells.
    pub fn point_count(&self, face: FaceId, level: FaceLevel) -> Option<usize> {
        self.counts[level as usize].get(&face).copied()
    }

    /// Assembles the certified mesh from the retained samples and topology.
    pub fn emit(self) -> Answered<Mesh> {
        self.emit_at(FaceLevel::Proven)
    }

    /// Emits the earliest retained display schedule, without a proof claim.
    pub fn emit_preview(self) -> Result<PreviewMesh, Refusal> {
        Ok(PreviewMesh::from_mesh(
            self.emit_at(FaceLevel::Preview)?.value,
        ))
    }

    /// Emits the retained positional schedule and its complete deviation bound.
    pub fn emit_bounded(self) -> Result<DisplayMesh, Refusal> {
        let answer = self.emit_at(FaceLevel::Bounded)?;
        Ok(display::rendering(
            answer.value,
            answer.marks,
            "complete-cell chord bounds on the retained face schedule",
        ))
    }

    /// Emits independently selected logical face levels. All incident faces
    /// use the same shared edge samples, regardless of the selection.
    pub fn emit_selected(
        self,
        level: impl FnMut(FaceId) -> FaceLevel,
    ) -> Result<PreviewMesh, Refusal> {
        Ok(PreviewMesh::from_mesh(self.emit_selected_at(level)?.value))
    }

    fn emit_at(self, level: FaceLevel) -> Answered<Mesh> {
        self.emit_selected_at(|_| level)
    }

    fn emit_selected_at(self, mut level: impl FnMut(FaceId) -> FaceLevel) -> Answered<Mesh> {
        let Some(expansion) = self.expansion else {
            return self.schedule.emit_selected(self.model, level);
        };
        let choices: BTreeMap<_, _> = self.counts[FaceLevel::Proven as usize]
            .keys()
            .map(|&face| (face, level(face)))
            .collect();
        let mut answer = self
            .schedule
            .emit_selected(self.model, |face| choices[&expansion.face(face)])?;
        rational_atlas::remap(self.model, expansion, &mut answer.value)?;
        Ok(answer)
    }
}

/// Runs shared-edge sampling and face refinement once before render emission.
///
/// Earlier admissible rounds reuse the same samples as the final certified
/// round. A face that needs its final topology at every level reports the
/// same count at each level.
/// The preview and positional targets are ten and two times the requested
/// chord, respectively. Optional facet controls apply to the certified mesh.
pub fn plan_face_levels<'a>(
    model: &'a Model,
    body: BodyId,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<FaceLevelPlan<'a>, Refusal> {
    prepare_face_levels(model, body, options, None, dec)
}

/// Cost of planning one physical face, excluding shared edge preparation.
#[derive(Clone, Debug)]
pub struct FacePreflight {
    /// The face in the planned topology, which may be an expanded chart cell.
    pub face: FaceId,
    /// Physical carrier of the planned face.
    pub kind: catallax_shapes::SurfaceKind,
    /// Time from preparing face-local inputs through retaining its topology.
    pub elapsed: std::time::Duration,
    /// Charged work during that interval, when an operation budget is active.
    pub work: u64,
}

/// Like [`plan_face_levels`], with per-carrier measurements for a benchmark.
pub fn plan_face_levels_profiled<'a>(
    model: &'a Model,
    body: BodyId,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<(FaceLevelPlan<'a>, Vec<FacePreflight>), Refusal> {
    let mut profile = Vec::new();
    let plan = prepare_face_levels(model, body, options, Some(&mut profile), dec)?;
    Ok((plan, profile))
}

fn prepare_face_levels<'a>(
    model: &'a Model,
    body: BodyId,
    options: &MeshOptions,
    profile: Option<&mut Vec<FacePreflight>>,
    dec: &mut Decider,
) -> Result<FaceLevelPlan<'a>, Refusal> {
    validate_mesh_options(options)?;
    let expansion = model
        .chart_expansion(body)
        .filter(|cover| cover.needs_cells());
    let prepared = expansion.and_then(|cover| {
        catallax_numbers::reuse_operation_enclosures();
        rational_atlas::expanded_surface(model, body, cover)
    });
    let (mesh_model, mesh_body) = prepared.as_ref().map_or_else(
        || (model, expansion.map_or(body, |cover| cover.body())),
        |(model, body)| (model, *body),
    );
    let schedule = triangulate::prepare(mesh_model, mesh_body, options, None, true, profile, dec)?;
    let counts = FaceLevel::ALL.map(|level| {
        schedule.point_counts(level, |face| {
            expansion.map_or(face, |cover| cover.face(face))
        })
    });
    Ok(FaceLevelPlan {
        model,
        expansion,
        schedule,
        counts,
    })
}

fn mesh_body_impl(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    edit: Option<reuse::Edit<'_, '_>>,
    dec: &mut Decider,
) -> Answered<Mesh> {
    #[cfg(feature = "display-capture")]
    let call = capture::Call::start("mesh_body_impl", capture::carriers(model, body));
    let result = mesh_body_checked(model, body, options, edit, dec);
    #[cfg(feature = "display-capture")]
    call.finish(options, &result);
    result
}

fn mesh_body_checked(
    model: &Model,
    body: BodyId,
    options: &MeshOptions,
    edit: Option<reuse::Edit<'_, '_>>,
    dec: &mut Decider,
) -> Answered<Mesh> {
    if let Some(expansion) = model
        .chart_expansion(body)
        .filter(|cover| cover.needs_cells())
    {
        return rational_atlas::build(model, body, expansion, options, edit, dec);
    }
    validate_mesh_options(options)?;
    triangulate::build(model, body, options, edit, dec)
}

fn validate_mesh_options(options: &MeshOptions) -> Result<(), Refusal> {
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
    Ok(())
}

fn failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("mesh.precision-or-budget", detail)
}

#[cfg(test)]
mod tests;
