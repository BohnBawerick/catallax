//! Mesh Booleans and decimation as recipe steps, RECIPE.md sections 4.21 and
//! 4.22.
//!
//! A mesh step keeps its exact triangles beside its B-rep body, so the next
//! mesh step reads triangles, not faces. Each result's B-rep is built as an
//! imported mesh's is (section 4.19): coplanar neighbours from the same step
//! join one planar face.

use std::sync::Arc;

use catallax_import::mesh::{ClosedMesh, boolean, decimate};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{BooleanKind, Origin, Role, StepId};
use catallax_shapes::Isometry;
use catallax_topology::{BodyId, Model};

/// The exact triangles of a mesh step.
#[derive(Clone, Debug)]
pub(crate) struct MeshValue {
    /// The triangles in their own frame.
    mesh: Arc<ClosedMesh>,
    /// How they sit in the model; `None` when they are in model coordinates.
    placement: Option<Isometry>,
    /// Face origins, and each triangle's index into them; no labels means
    /// every triangle takes the first.
    origins: Arc<[Origin]>,
    labels: Arc<[u32]>,
}

impl MeshValue {
    /// An imported mesh as read, before its placement.
    pub(crate) fn imported(mesh: ClosedMesh, placement: Isometry, step: &StepId) -> Self {
        Self {
            mesh: Arc::new(mesh),
            placement: Some(placement),
            origins: Arc::from([Origin {
                step: step.clone(),
                role: Role::Mesh,
            }]),
            labels: Arc::from([]),
        }
    }

    fn label(&self, t: u32) -> u32 {
        self.labels.get(t as usize).copied().unwrap_or(0)
    }

    /// The triangles in model coordinates.
    fn placed(&self, dec: &mut Decider) -> Result<Arc<ClosedMesh>, Refusal> {
        let Some(placement) = &self.placement else {
            return Ok(Arc::clone(&self.mesh));
        };
        let motion = placement
            .rational_matrix()
            .zip(placement.translation.rational_components());
        let Some((rotation, translation)) = motion else {
            return Err(Refusal::new(
                "import.mesh-pose",
                "a mesh pose must be an exact rational motion",
            ));
        };
        Ok(Arc::new(self.mesh.placed(rotation, translation, dec)?))
    }

    /// The result value and its B-rep body with its bounding radius.
    fn build(
        model: &mut Model,
        mesh: ClosedMesh,
        origins: Vec<Origin>,
        labels: Vec<u32>,
        dec: &mut Decider,
    ) -> Result<(BodyId, Rational, Self), Refusal> {
        let value = Self {
            mesh: Arc::new(mesh),
            placement: None,
            origins: origins.into(),
            labels: labels.into(),
        };
        let (body, radius) = crate::mesh_import::build_labelled(
            model,
            &value.mesh,
            &Isometry::identity(),
            &value.origins,
            &value.labels,
            dec,
        )?;
        Ok((body, radius, value))
    }
}

/// `mesh_boolean`: the exact Boolean of two mesh steps.
pub(crate) fn boolean(
    model: &mut Model,
    kind: &BooleanKind,
    target: &MeshValue,
    tool: &MeshValue,
    dec: &mut Decider,
) -> Result<(BodyId, Rational, MeshValue), Refusal> {
    let (a, b) = (target.placed(dec)?, tool.placed(dec)?);
    let kind = match kind {
        BooleanKind::Union => boolean::Kind::Union,
        BooleanKind::Subtract => boolean::Kind::Subtract,
        BooleanKind::Intersect => boolean::Kind::Intersect,
    };
    let outcome = boolean::boolean(&a, &b, kind, dec)?;
    // One origin list for both operands, each origin once, in operand order.
    let mut origins: Vec<Origin> = Vec::new();
    let mut slot = |origin: &Origin| -> u32 {
        match origins.iter().position(|o| o == origin) {
            Some(i) => i as u32,
            None => {
                origins.push(origin.clone());
                origins.len() as u32 - 1
            }
        }
    };
    let maps: [Vec<u32>; 2] = [target, tool].map(|v| v.origins.iter().map(&mut slot).collect());
    let labels = outcome
        .sources
        .iter()
        .map(|&(side, t)| {
            let value = if side == 0 { target } else { tool };
            maps[side as usize][value.label(t) as usize]
        })
        .collect();
    MeshValue::build(model, outcome.mesh, origins, labels, dec)
}

/// `mesh_decimate`: the mesh step with edges collapsed towards `triangles`.
pub(crate) fn decimate(
    model: &mut Model,
    input: &MeshValue,
    triangles: u32,
    dec: &mut Decider,
) -> Result<(BodyId, Rational, MeshValue), Refusal> {
    let placed = input.placed(dec)?;
    let outcome = decimate::decimate(&placed, triangles as usize, dec)?;
    let labels = outcome.sources.iter().map(|&t| input.label(t)).collect();
    MeshValue::build(model, outcome.mesh, input.origins.to_vec(), labels, dec)
}
