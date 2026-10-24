//! Application entries for named B-rep builds and their display approximations.
//!
//! [`build`] accepts Catallax recipes. [`build_soliddraft`] first converts the
//! application's JSON format. No SolidDraft backend or external checkout is needed.
//! [`build_certified`] adds geometry-bound inward offset prisms, optional display
//! and original-definition properties without changing the exact build entry.
//! Handles belong to the returned model, not to a later rebuild. Store origins and
//! names for persistent references; do not persist numeric handles.

mod certified;
/// Saved feature-pattern editing with explicit per-instance bindings.
pub mod pattern;
pub mod soliddraft;
pub use certified::{CertifiedApplicationBody, build_certified};

use std::{collections::BTreeMap, fmt};

use catallax_features::{BuildError, Built, EvaluationCheckpoint, evaluate, evaluate_incremental};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{Answer, Decider, Refusal, operation_budget_active, with_operation_budget};
use catallax_recipe::{EdgeName, Origin, Recipe, VertexName};
use catallax_topology::{BodyId, EdgeId, FaceId, Model, VertexId};

/// Default aggregate arithmetic allowance for one application build and display.
/// An existing caller-owned operation allowance takes precedence.
pub const DEFAULT_BUILD_WORK: u64 = 20_000_000;

fn bounded_build<T>(run: impl FnOnce() -> Result<T, Error>) -> Result<T, Error> {
    if operation_budget_active() {
        return run();
    }
    let (result, work) = with_operation_budget(DEFAULT_BUILD_WORK, run);
    if work.exhausted {
        Err(Error::ResourceLimit)
    } else {
        result
    }
}

/// Names of every entity in the result, keyed by that build's picking handles.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Names {
    /// Complete origin sets, including the union on merged faces.
    pub faces: BTreeMap<FaceId, Vec<Origin>>,
    /// An edge can lack a unique origin-pair name. Its handle remains pickable.
    pub edges: BTreeMap<EdgeId, Option<EdgeName>>,
    /// Vertex names include all incident faces' first origins.
    pub vertices: BTreeMap<VertexId, Option<VertexName>>,
}

/// An exact body and its separate, certified display approximation.
#[derive(Clone, Debug)]
pub struct ApplicationBody {
    /// Model, result body, intermediate bodies and immutable recipe definition.
    pub built: Built,
    /// Shared positions, face-labelled triangles and edge polylines. Marks are retained.
    pub display: Answer<Mesh>,
    /// The names behind every result picking handle.
    pub names: Names,
}

/// An application body whose exact evaluation can be consumed by the next edit.
#[derive(Clone, Debug)]
pub struct IncrementalApplicationBody {
    /// Owned prefix. Read its exact result through [`EvaluationCheckpoint::built`].
    pub checkpoint: EvaluationCheckpoint,
    /// Display approximation of that result.
    pub display: Answer<Mesh>,
    /// Picking names of that result.
    pub names: Names,
}

impl ApplicationBody {
    /// Requests proven mass properties, including centroidal unit-density inertia.
    /// A refusal leaves the exact body, picking names and display untouched.
    pub fn mass_properties(
        &self,
        accuracy: &catallax_measure::Accuracy,
        decider: &mut Decider,
    ) -> catallax_numbers::Answered<catallax_measure::MassProperties> {
        catallax_measure::mass_properties(&self.built.model, self.built.result, accuracy, decider)
    }

    /// Requests only centroidal unit-density inertia.
    pub fn inertia(
        &self,
        accuracy: &catallax_measure::Accuracy,
        decider: &mut Decider,
    ) -> catallax_numbers::Answered<Option<catallax_measure::InertiaTensor>> {
        catallax_measure::inertia(&self.built.model, self.built.result, accuracy, decider)
    }

    /// Requests the exact body's enclosing axis-aligned box, without display padding.
    pub fn bounding_box(
        &self,
        accuracy: &catallax_measure::Accuracy,
        decider: &mut Decider,
    ) -> catallax_numbers::Answered<Option<catallax_measure::BoundingBox>> {
        catallax_measure::bounding_box(&self.built.model, self.built.result, accuracy, decider)
    }

    /// Requests an area-weighted centroid of one exact face.
    pub fn face_centroid(
        &self,
        face: FaceId,
        accuracy: &catallax_measure::Accuracy,
        decider: &mut Decider,
    ) -> catallax_numbers::Answered<[catallax_measure::Enclosure; 3]> {
        catallax_measure::face_centroid(&self.built.model, face, accuracy, decider)
    }

    /// Requests a proven enclosing box of one exact face.
    pub fn face_bounding_box(
        &self,
        face: FaceId,
        accuracy: &catallax_measure::Accuracy,
        decider: &mut Decider,
    ) -> catallax_numbers::Answered<catallax_measure::BoundingBox> {
        catallax_measure::face_bounding_box(&self.built.model, face, accuracy, decider)
    }

    /// Requests the arc length of an exact trimmed edge.
    pub fn edge_length(
        &self,
        edge: EdgeId,
        accuracy: &catallax_measure::Accuracy,
        decider: &mut Decider,
    ) -> catallax_numbers::Answered<catallax_measure::Enclosure> {
        catallax_measure::edge_length(&self.built.model, edge, accuracy, decider)
    }

    /// Requests the half-arc-length point, not the mean of the edge's endpoints.
    pub fn edge_midpoint(
        &self,
        edge: EdgeId,
        accuracy: &catallax_measure::Accuracy,
        decider: &mut Decider,
    ) -> catallax_numbers::Answered<[catallax_measure::Enclosure; 3]> {
        catallax_measure::edge_midpoint(&self.built.model, edge, accuracy, decider)
    }
}

/// Failure stage. No variant returns a partial body or substitutes approximate geometry.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum Error {
    /// The application input could not be converted.
    Conversion(soliddraft::ConversionError),
    /// Structural recipe error or step-local geometric refusal, unchanged from features.
    Build(BuildError),
    /// The exact body failed the topology check.
    InvalidBody(Vec<catallax_topology::TopologyFault>),
    /// Display failed, with its original stable refusal code.
    Display(Refusal),
    /// The default whole-build work allowance expired; no body is published.
    ResourceLimit,
}

impl Error {
    /// Stable code for UI error classification, without parsing human explanations.
    pub fn code(&self) -> &'static str {
        match self {
            Self::Conversion(soliddraft::ConversionError::Invalid { .. }) => "soliddraft.invalid",
            Self::Conversion(soliddraft::ConversionError::NotYet { code, .. }) => code,
            Self::Build(BuildError::Recipe(error)) => error.code(),
            Self::Build(BuildError::Geometry { refusal, .. }) | Self::Display(refusal) => {
                refusal.code
            }
            Self::InvalidBody(_) => "application.invalid-body",
            Self::ResourceLimit => "application.work-limit",
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Conversion(e) => e.fmt(f),
            Self::Build(e) => e.fmt(f),
            Self::InvalidBody(faults) => write!(f, "application.invalid-body: {faults:?}"),
            Self::Display(e) => write!(f, "{}: {}", e.code, e.detail),
            Self::ResourceLimit => {
                write!(f, "application.work-limit: aggregate build work exhausted")
            }
        }
    }
}
impl std::error::Error for Error {}

/// Builds and checks the exact body before deriving display triangles and names.
/// Display controls never change the recipe, topology or origin sets.
pub fn build(
    recipe: &Recipe,
    options: &MeshOptions,
    decider: &mut Decider,
) -> Result<ApplicationBody, Error> {
    bounded_build(|| build_inner(recipe, options, decider))
}

fn build_inner(
    recipe: &Recipe,
    options: &MeshOptions,
    decider: &mut Decider,
) -> Result<ApplicationBody, Error> {
    let built = evaluate(recipe, decider).map_err(Error::Build)?;
    let (display, names) = display_and_names(&built, options, decider)?;
    Ok(ApplicationBody {
        built,
        display,
        names,
    })
}

/// Resumes an unchanged prefix, then checks and meshes only the requested result.
/// A prefix mismatch rebuilds from the start, as in [`evaluate_incremental`].
pub fn build_incremental(
    recipe: &Recipe,
    prefix: Option<EvaluationCheckpoint>,
    options: &MeshOptions,
    decider: &mut Decider,
) -> Result<IncrementalApplicationBody, Error> {
    let checkpoint = evaluate_incremental(recipe, prefix, decider).map_err(Error::Build)?;
    let (display, names) = display_and_names(checkpoint.built(), options, decider)?;
    Ok(IncrementalApplicationBody {
        checkpoint,
        display,
        names,
    })
}

fn display_and_names(
    built: &Built,
    options: &MeshOptions,
    decider: &mut Decider,
) -> Result<(Answer<Mesh>, Names), Error> {
    let faults = built.model.check(built.result);
    if !faults.is_empty() {
        return Err(Error::InvalidBody(faults));
    }
    let display =
        mesh_body(&built.model, built.result, options, decider).map_err(Error::Display)?;
    let names = result_names(&built.model, built.result);
    Ok((display, names))
}

fn result_names(model: &Model, body: BodyId) -> Names {
    Names {
        faces: model
            .faces_of_body(body)
            .into_iter()
            .map(|face| (face, model.face_origins(face).to_vec()))
            .collect(),
        edges: model
            .edges_of_body(body)
            .into_iter()
            .map(|edge| (edge, model.edge_name(body, edge)))
            .collect(),
        vertices: model
            .vertices_of_body(body)
            .into_iter()
            .map(|vertex| (vertex, model.vertex_name(body, vertex)))
            .collect(),
    }
}

/// A SolidDraft build, with the conversion's original-tag to step-id table.
#[derive(Clone, Debug)]
pub struct SolidDraftBody {
    /// The body, mesh and names.
    pub body: ApplicationBody,
    /// Original tags, including those whose spelling needed normalisation.
    pub steps: BTreeMap<String, catallax_recipe::StepId>,
}

/// Converts `soliddraft-recipe/1` JSON and builds it through [`build`].
pub fn build_soliddraft(
    json: &str,
    options: &MeshOptions,
    decider: &mut Decider,
) -> Result<SolidDraftBody, Error> {
    let converted = soliddraft::convert(json).map_err(Error::Conversion)?;
    let body = build(&converted.recipe, options, decider)?;
    Ok(SolidDraftBody {
        body,
        steps: converted.steps,
    })
}
