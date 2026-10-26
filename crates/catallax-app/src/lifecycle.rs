//! Caller-scheduled preview and certified display for one immutable edit.
//!
//! Construction keeps its existing proofs. Only certified display is deferred;
//! this entry does not promise a whole-model or second-kernel certificate.

use std::sync::{Arc, Mutex};

use catallax_features::Built;
use catallax_mesh::{MeshOptions, PreviewOptions};
use catallax_numbers::{Decider, OperationAllowance, OperationWork, Refusal};
use catallax_recipe::{EdgeName, Origin, Recipe, VertexName};
use catallax_topology::{EdgeId, FaceId, VertexId};

use crate::{ApplicationBody, Error, PreviewDisplay, build_preview, display_and_names};

/// Explicit caps for construction plus preview, and for subsequent proof.
/// Unused units never move between phases. The default total remains 20M;
/// callers choose its split rather than receiving two default allowances.
#[derive(Clone, Copy, Debug)]
pub struct WorkAllocation {
    /// Construction, unproved mesh, edge polylines and selection names.
    pub preview: u64,
    /// Topology validation and certified display of the retained exact build.
    pub proof: u64,
}

impl WorkAllocation {
    /// Splits the existing application default. Returns `None` for a preview
    /// cap above the total. Explicit caller budgets may use the public fields.
    pub fn application(preview: u64) -> Option<Self> {
        Some(Self {
            preview,
            proof: crate::DEFAULT_BUILD_WORK.checked_sub(preview)?,
        })
    }

    fn reserve(self) -> Result<(OperationAllowance, OperationAllowance), Error> {
        let total = self.preview.checked_add(self.proof).ok_or_else(|| {
            Error::Allocation(Refusal {
                code: "application.invalid-allocation",
                detail: "preview and proof caps overflow the total work allowance".into(),
            })
        })?;
        let mut preview = OperationAllowance::reserve(total).map_err(Error::Allocation)?;
        let proof = preview.split_off(self.proof).expect("validated allocation");
        Ok((preview, proof))
    }
}

/// Actual arithmetic performed in each phase, distinct from reserved caps.
#[derive(Clone, Copy, Debug)]
pub struct WorkReport {
    /// Reserved work, including unused or cancelled units.
    pub allocation: WorkAllocation,
    /// Actual construction and preview work.
    pub preview: OperationWork,
    /// Actual proof work; absent until the caller runs proof.
    pub proof: Option<OperationWork>,
}

/// Delivery state. A confirmed result has the same guarantees as [`crate::build`],
/// not a whole-boundary certificate. Requested measurements and exports must use
/// the confirmed exact result and perform their own scoped proof.
#[derive(Clone, Debug)]
pub enum Status {
    /// Visible geometry and names only, without a certified display.
    Unproved,
    /// The retained exact build passed topology and certified display.
    Confirmed(Arc<ApplicationBody>),
    /// Proof refused. Keep this refusal with the visible failed edit.
    Failed(Error),
}

/// Immutable delivery. Its fields cannot be paired with another revision.
#[derive(Clone, Debug)]
pub struct Snapshot {
    revision: u64,
    identity: Arc<()>,
    preview: Arc<PreviewDisplay>,
    status: Status,
    work: WorkReport,
}

impl Snapshot {
    /// Monotonic identifier within this controller.
    pub fn revision(&self) -> u64 {
        self.revision
    }
    /// Unproved mesh and revision-local numeric picking handles.
    pub fn preview(&self) -> &PreviewDisplay {
        &self.preview
    }
    /// Whether this edit is awaiting proof, confirmed or failed.
    pub fn status(&self) -> &Status {
        &self.status
    }
    /// Reserved caps and actual per-phase work.
    pub fn work(&self) -> WorkReport {
        self.work
    }
    /// Saves origins instead of carrying a numeric face handle into proof.
    pub fn select_face(&self, face: FaceId) -> Option<Selection> {
        Some(self.selection(SelectionName::Face(
            self.preview.names.faces.get(&face)?.clone(),
        )))
    }
    /// Saves an edge's incident-origin name, if it has one.
    pub fn select_edge(&self, edge: EdgeId) -> Option<Selection> {
        Some(self.selection(SelectionName::Edge(
            self.preview.names.edges.get(&edge)?.clone()?,
        )))
    }
    /// Saves a vertex's incident-origin name, if it has one.
    pub fn select_vertex(&self, vertex: VertexId) -> Option<Selection> {
        Some(self.selection(SelectionName::Vertex(
            self.preview.names.vertices.get(&vertex)?.clone()?,
        )))
    }
    fn selection(&self, name: SelectionName) -> Selection {
        Selection {
            identity: self.identity.clone(),
            name,
        }
    }
    /// Resolves only against this revision's confirmed result, never against
    /// preview IDs or another controller whose revision number happens to match.
    /// Returned handles belong only to the confirmed body in this snapshot.
    pub fn resolve(&self, selection: &Selection) -> Option<ResolvedSelection> {
        if !Arc::ptr_eq(&self.identity, &selection.identity) {
            return None;
        }
        let Status::Confirmed(body) = &self.status else {
            return None;
        };
        let model = &body.built.model;
        let result = body.built.result;
        Some(match &selection.name {
            SelectionName::Face(origins) => ResolvedSelection::Faces(
                body.names
                    .faces
                    .iter()
                    .filter_map(|(&face, names)| {
                        origins
                            .iter()
                            .any(|origin| names.contains(origin))
                            .then_some(face)
                    })
                    .collect(),
            ),
            SelectionName::Edge(name) => {
                ResolvedSelection::Edges(model.edges_with_name(result, name))
            }
            SelectionName::Vertex(name) => ResolvedSelection::Vertices(
                body.names
                    .vertices
                    .iter()
                    .filter_map(|(&vertex, candidate)| {
                        (candidate.as_ref() == Some(name)).then_some(vertex)
                    })
                    .collect(),
            ),
        })
    }
}

#[derive(Clone, Debug)]
enum SelectionName {
    Face(Vec<Origin>),
    Edge(EdgeName),
    Vertex(VertexName),
}

/// A selection name bound to one edit, not a persistent numeric handle.
#[derive(Clone, Debug)]
pub struct Selection {
    identity: Arc<()>,
    name: SelectionName,
}

/// Named pieces of the confirmed result. Empty lists mean no matching piece.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ResolvedSelection {
    /// All split or merged faces containing a selected origin.
    Faces(Vec<FaceId>),
    /// Edges with the selected incident origins.
    Edges(Vec<EdgeId>),
    /// Vertices with the selected incident origins.
    Vertices(Vec<VertexId>),
}

#[derive(Default)]
struct State {
    revision: u64,
    identity: Option<Arc<()>>,
    latest: Option<Snapshot>,
}

/// Atomic edit invalidation and publication. Never starts worker threads.
#[derive(Clone, Default)]
pub struct Controller(Arc<Mutex<State>>);

impl Controller {
    /// Invalidates old preview/proof work immediately, before constructing the
    /// new edit. Previously retained snapshots still describe their old edit.
    pub fn begin(&self) -> Revision {
        let mut state = self.0.lock().unwrap();
        state.revision = state
            .revision
            .checked_add(1)
            .expect("revision space exhausted");
        let identity = Arc::new(());
        state.identity = Some(identity.clone());
        state.latest = None;
        Revision {
            controller: self.clone(),
            id: state.revision,
            identity,
        }
    }
    /// Current delivery, without copying geometry or mesh buffers.
    pub fn latest(&self) -> Option<Snapshot> {
        self.0.lock().unwrap().latest.clone()
    }
}

/// Single-use authority to construct and deliver one edit.
pub struct Revision {
    controller: Controller,
    id: u64,
    identity: Arc<()>,
}

impl Revision {
    /// Constructs a complete selectable preview without certified meshing.
    /// Both phase caps are reserved now from an active caller allowance, or
    /// from the explicit standalone total. The reservation is never refunded.
    /// A stale revision returns `Ok(None)` and publishes nothing.
    pub fn preview(
        self,
        recipe: &Recipe,
        preview_options: &PreviewOptions,
        proof_options: &MeshOptions,
        allocation: WorkAllocation,
        decider: &mut Decider,
    ) -> Result<Option<ProofJob>, Error> {
        if !self.is_current() {
            return Ok(None);
        }
        let (preview_allowance, proof_allowance) = allocation.reserve()?;
        let (answer, work) =
            preview_allowance.run(|| build_preview(recipe, preview_options, decider));
        let body = answer?;
        if work.exhausted {
            return Err(Error::ResourceLimit);
        }
        let snapshot = Snapshot {
            revision: self.id,
            identity: self.identity.clone(),
            preview: Arc::new(body.preview),
            status: Status::Unproved,
            work: WorkReport {
                allocation,
                preview: work,
                proof: None,
            },
        };
        if !self.publish(snapshot.clone()) {
            return Ok(None);
        }
        Ok(Some(ProofJob {
            revision: self,
            snapshot,
            built: body.built,
            options: proof_options.clone(),
            allowance: proof_allowance,
        }))
    }
    fn is_current(&self) -> bool {
        self.controller
            .0
            .lock()
            .unwrap()
            .identity
            .as_ref()
            .is_some_and(|identity| Arc::ptr_eq(identity, &self.identity))
    }
    fn publish(&self, snapshot: Snapshot) -> bool {
        let mut state = self.controller.0.lock().unwrap();
        if !state
            .identity
            .as_ref()
            .is_some_and(|identity| Arc::ptr_eq(identity, &self.identity))
        {
            return false;
        }
        state.latest = Some(snapshot);
        true
    }
}

/// Owned retained geometry and reserved proof work. Dropping it postpones proof
/// indefinitely; it never changes an unproved preview into success.
pub struct ProofJob {
    revision: Revision,
    snapshot: Snapshot,
    built: Built,
    options: MeshOptions,
    allowance: OperationAllowance,
}

impl ProofJob {
    /// Runs proof when the caller schedules it. No recipe rebuild and no new
    /// allowance. Publication is separate so a queued completion can be stale.
    /// The caller owns the decider's policy and wall-clock deadline.
    pub fn prove(self, decider: &mut Decider) -> ProofCompletion {
        let Self {
            revision,
            mut snapshot,
            built,
            options,
            allowance,
        } = self;
        if !revision.is_current() {
            return ProofCompletion {
                revision,
                snapshot: None,
            };
        }
        let (answer, work) = allowance.run(|| display_and_names(&built, &options, decider));
        snapshot.work.proof = Some(work);
        snapshot.status = match answer {
            Err(error) => Status::Failed(error),
            Ok(_) if work.exhausted => Status::Failed(Error::ResourceLimit),
            Ok((display, names)) => Status::Confirmed(Arc::new(ApplicationBody {
                built,
                display,
                names,
            })),
        };
        ProofCompletion {
            revision,
            snapshot: Some(snapshot),
        }
    }
}

/// A completed job awaiting atomic publication. Cannot be forged or replayed.
pub struct ProofCompletion {
    revision: Revision,
    snapshot: Option<Snapshot>,
}

impl ProofCompletion {
    /// Confirms or marks the same preview failed. Returns false for stale work;
    /// neither a late success nor a late refusal can replace a newer edit.
    pub fn publish(self) -> bool {
        self.snapshot
            .is_some_and(|snapshot| self.revision.publish(snapshot))
    }
}
