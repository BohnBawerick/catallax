//! Exact ordered face boundaries, independent of the display mesh.
//!
//! [`ApplicationBody::face_loops`] and [`IncrementalApplicationBody::face_loops`]
//! return a result face's stored boundaries exactly as the topology holds them:
//! loops in stored order, coedges in stored cycle order, each with its travel
//! relative to the edge and the edge's exact [`EdgeCurve`]. No loop is labelled
//! outer or inner, and the order is that of one model, not a promise across
//! rebuilds. A face with zero loops is a whole closed surface. A vertex loop
//! is a single point boundary such as a cone apex.
//!
//! Coedge `sense` is the direction of travel relative to the edge direction, and
//! the edge's own `sense` is relative to its carrier. The face `sense` is the
//! face normal relative to its surface normal. [`LoopCoedge::start`] and
//! [`LoopCoedge::end`] are the vertices in travel order, `None` for a ring edge.
//! Per-edge `span` stays an explicit capability result, so an unsupported
//! carrier does not hide the rest of the face. An underlying refusal propagates
//! and leaves the body untouched. A face outside the result returns `app.loop-face`.
//!
//! Handles belong to the queried model. Persist origin names, not handles.

use std::collections::{BTreeMap, btree_map::Entry};

use catallax_numbers::{Decider, Refusal};
use catallax_recipe::{Origin, VertexName};
use catallax_shapes::SurfaceId;
use catallax_topology::{BodyId, CoedgeId, EdgeId, FaceId, Loop, LoopId, Model, Sense, VertexId};

use crate::curves::{CurveEndpoint, EdgeCurve, edge_curve};
use crate::{ApplicationBody, IncrementalApplicationBody};

/// One use of an edge in a loop, in stored cycle order.
#[derive(Clone, Debug)]
pub struct LoopCoedge {
    /// Model-local coedge handle.
    pub coedge: CoedgeId,
    /// Direction of travel relative to the edge direction.
    pub sense: Sense,
    /// The edge's exact query, shared in value with any other use of the edge.
    pub edge: EdgeCurve,
    /// Vertex where travel along this coedge starts; `None` for a ring edge.
    pub start: Option<VertexId>,
    /// Vertex where travel along this coedge ends; `None` for a ring edge.
    pub end: Option<VertexId>,
}

/// A face boundary.
#[derive(Clone, Debug)]
pub enum FaceLoop {
    /// A closed cycle of coedges in stored order.
    Coedges {
        /// Model-local loop handle.
        id: LoopId,
        /// The cycle. Consecutive travel ends meet the next travel start.
        coedges: Vec<LoopCoedge>,
    },
    /// A single point boundary.
    Vertex {
        /// Model-local loop handle.
        id: LoopId,
        /// The point, with its exact definition and coordinates.
        endpoint: Box<CurveEndpoint>,
        /// Stable vertex name when the result admits one.
        name: Option<VertexName>,
    },
}

/// A face's surface, sense, origins and ordered boundaries.
#[derive(Clone, Debug)]
pub struct FaceLoops {
    /// Model-local face handle.
    pub face: FaceId,
    /// Model-local surface handle.
    pub surface: SurfaceId,
    /// Face normal relative to the surface normal.
    pub sense: Sense,
    /// Complete canonical origin set, including the union on merged faces.
    pub origins: Vec<Origin>,
    /// Loops in stored order, empty for a whole closed surface.
    pub loops: Vec<FaceLoop>,
}

impl ApplicationBody {
    /// Queries one result face's ordered boundaries without reading the display mesh.
    /// Uses the caller's decider and any active work allowance.
    pub fn face_loops(&self, face: FaceId, decider: &mut Decider) -> Result<FaceLoops, Refusal> {
        face_loops(&self.built.model, self.built.result, face, decider)
    }
}

impl IncrementalApplicationBody {
    /// Queries a face of this checkpoint's result, with the same contract as
    /// [`ApplicationBody::face_loops`]. Handles from another edit are not stable.
    pub fn face_loops(&self, face: FaceId, decider: &mut Decider) -> Result<FaceLoops, Refusal> {
        let built = self.checkpoint.built();
        face_loops(&built.model, built.result, face, decider)
    }
}

/// Queries one face of `body`, composing [`edge_curve`] for every distinct edge.
///
/// A face outside the body returns `app.loop-face`. Edge-query refusals
/// propagate. No healing, ordering or outer-loop inference is applied.
pub fn face_loops(
    model: &Model,
    body: BodyId,
    face: FaceId,
    decider: &mut Decider,
) -> Result<FaceLoops, Refusal> {
    decider.check_deadline()?;
    if !model.faces_of_body(body).contains(&face) {
        return Err(Refusal::new(
            "app.loop-face",
            "face is not in the queried body",
        ));
    }
    let stored = model.face(face);
    let mut edges: BTreeMap<EdgeId, EdgeCurve> = BTreeMap::new();
    let mut loops = Vec::with_capacity(stored.loops.len());
    for &id in &stored.loops {
        decider.check_deadline()?;
        loops.push(match model.loop_(id) {
            Loop::Vertex(vertex) => {
                let point = model.vertex(*vertex).point;
                FaceLoop::Vertex {
                    id,
                    endpoint: Box::new(CurveEndpoint {
                        vertex: *vertex,
                        point,
                        definition: model.geom().point(point).clone(),
                        position: model.geom().position(point),
                    }),
                    name: model.vertex_name(body, *vertex),
                }
            }
            Loop::Coedges(cycle) => {
                let mut coedges = Vec::with_capacity(cycle.len());
                for &coedge in cycle {
                    let use_ = model.coedge(coedge);
                    if let Entry::Vacant(slot) = edges.entry(use_.edge) {
                        slot.insert(edge_curve(model, body, use_.edge, decider)?);
                    }
                    let source = model.edge(use_.edge);
                    let (start, end) = match use_.sense {
                        Sense::Same => (source.start, source.end),
                        Sense::Reversed => (source.end, source.start),
                    };
                    coedges.push(LoopCoedge {
                        coedge,
                        sense: use_.sense,
                        edge: edges[&use_.edge].clone(),
                        start,
                        end,
                    });
                }
                FaceLoop::Coedges { id, coedges }
            }
        });
    }
    Ok(FaceLoops {
        face,
        surface: stored.surface,
        sense: stored.sense,
        origins: model.face_origins(face).to_vec(),
        loops,
    })
}

#[cfg(test)]
mod tests;
