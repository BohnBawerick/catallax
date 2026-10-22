//! Application entry for named, exact B-rep builds and their display approximations.
//!
//! [`build`] accepts Catallax recipes. [`build_soliddraft`] first converts the
//! application's JSON format. No SolidDraft backend or external checkout is needed.
//! Handles belong to the returned model, not to a later rebuild. Store origins and
//! names for persistent references; do not persist numeric handles.

pub mod soliddraft;

use std::{collections::BTreeMap, fmt};

use catallax_features::{BuildError, Built, evaluate};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{Answer, Decider, Refusal};
use catallax_recipe::{EdgeName, Origin, Recipe, VertexName};
use catallax_topology::{EdgeId, FaceId, VertexId};

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
    let built = evaluate(recipe, decider).map_err(Error::Build)?;
    let faults = built.model.check(built.result);
    if !faults.is_empty() {
        return Err(Error::InvalidBody(faults));
    }
    let display =
        mesh_body(&built.model, built.result, options, decider).map_err(Error::Display)?;
    let model = &built.model;
    let body = built.result;
    let names = Names {
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
    };
    Ok(ApplicationBody {
        built,
        display,
        names,
    })
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
