//! Shared boundary replacement for elementary edge blends and bevels.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::Refusal;
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Lump, Model, Sense, Shell, VertexId};

/// A gap closes with one edge or a connected oriented path. Applying a path
/// may reverse the whole chain, never its individual connections.
pub(crate) enum Gap {
    Edge(EdgeId),
    Path(Vec<(EdgeId, Sense)>),
}

impl Gap {
    fn len(&self) -> usize {
        match self {
            Self::Edge(_) => 1,
            Self::Path(path) => path.len(),
        }
    }

    fn append(
        &self,
        model: &Model,
        end: Option<VertexId>,
        start: Option<VertexId>,
        closed: &mut Vec<(EdgeId, Sense)>,
    ) -> Result<(), Refusal> {
        let wrong = || {
            Refusal::new(
                "edge-blend.wrong-connection",
                "connector path does not join the two replacement rails",
            )
        };
        match self {
            Self::Edge(edge) => {
                let sense = if ends(model, *edge, Sense::Same) == (end, start) {
                    Sense::Same
                } else if ends(model, *edge, Sense::Reversed) == (end, start) {
                    Sense::Reversed
                } else {
                    return Err(wrong());
                };
                closed.push((*edge, sense));
            }
            Self::Path(path) => {
                let (Some(first), Some(last)) = (path.first(), path.last()) else {
                    return Err(wrong());
                };
                let first = ends(model, first.0, first.1).0;
                let last = ends(model, last.0, last.1).1;
                let forward = if (first, last) == (end, start) {
                    true
                } else if (last, first) == (end, start) {
                    false
                } else {
                    return Err(wrong());
                };
                let mut at = end;
                for offset in 0..path.len() {
                    let index = if forward {
                        offset
                    } else {
                        path.len() - 1 - offset
                    };
                    let (edge, sense) = path[index];
                    let sense = if forward {
                        sense
                    } else {
                        crate::extrude::reverse(sense)
                    };
                    let (from, to) = ends(model, edge, sense);
                    if from != at {
                        return Err(wrong());
                    }
                    closed.push((edge, sense));
                    at = to;
                }
                if at != start {
                    return Err(wrong());
                }
            }
        }
        Ok(())
    }
}

/// Replacements refer to the immutable input skin. A boundary shared by two
/// output faces has one edge handle; only its two coedges differ. A removed
/// edge collapses to a point of its face, where its neighbours in the loop
/// already meet. A consumed face leaves its shell with no replacement of
/// its own.
#[derive(Default)]
pub(crate) struct Replacement {
    pub edges: BTreeMap<(FaceId, EdgeId), EdgeId>,
    pub consumed: BTreeSet<FaceId>,
    pub removed: BTreeSet<(FaceId, EdgeId)>,
    pub gaps: BTreeMap<(FaceId, VertexId), Gap>,
    pub faces: Vec<(FaceId, FaceId)>,
}

pub(crate) fn ends(
    model: &Model,
    edge: EdgeId,
    sense: Sense,
) -> (Option<VertexId>, Option<VertexId>) {
    let e = model.edge(edge);
    if sense == Sense::Same {
        (e.start, e.end)
    } else {
        (e.end, e.start)
    }
}

pub(crate) fn apply(
    model: &mut Model,
    body: BodyId,
    edit: &Replacement,
) -> Result<BodyId, Refusal> {
    let input = model.body(body).clone();
    let gap_edges: usize = edit.gaps.values().map(Gap::len).sum();
    let mut lumps = Vec::with_capacity(input.lumps.len());
    for lump in &input.lumps {
        let input_shells = model.lump(*lump).shells.clone();
        let mut shells = Vec::with_capacity(input_shells.len());
        for shell in input_shells {
            let input_faces = model.shell(shell).faces.clone();
            let mut faces = Vec::with_capacity(input_faces.len() + edit.faces.len());
            for old in &input_faces {
                if edit.consumed.contains(old) {
                    continue;
                }
                if !edit.edges.keys().any(|(face, _)| face == old)
                    && !edit.removed.iter().any(|(face, _)| face == old)
                    && !edit.gaps.keys().any(|(face, _)| face == old)
                {
                    faces.push(*old);
                    continue;
                }
                let face = model.face(*old).clone();
                let mut paths = Vec::with_capacity(face.loops.len());
                let mut apex = None;
                for lp in &face.loops {
                    let Loop::Coedges(cs) = model.loop_(*lp) else {
                        if let Loop::Vertex(v) = model.loop_(*lp)
                            && apex.replace(*v).is_some()
                        {
                            return Err(Refusal::new(
                                "edge-blend.singular-loops",
                                "a changed face with multiple point loops needs explicit singular trimming",
                            ));
                        }
                        continue;
                    };
                    // Each kept coedge with the input coedge it replaces.
                    let mut path = Vec::with_capacity(cs.len() + gap_edges);
                    for c in cs {
                        let original = model.coedge(*c);
                        if edit.removed.contains(&(*old, original.edge)) {
                            continue;
                        }
                        let edge = edit
                            .edges
                            .get(&(*old, original.edge))
                            .copied()
                            .unwrap_or(original.edge);
                        path.push(((edge, original.sense), *c));
                    }
                    let mut closed = Vec::with_capacity(path.len() + gap_edges);
                    for i in 0..path.len() {
                        let (current, input) = path[i];
                        closed.push(current);
                        let (next, _) = path[(i + 1) % path.len()];
                        let end = ends(model, current.0, current.1).1;
                        let start = ends(model, next.0, next.1).0;
                        if end == start {
                            continue;
                        }
                        let original = model.coedge(input);
                        let vertex =
                            ends(model, original.edge, original.sense)
                                .1
                                .ok_or_else(|| {
                                    Refusal::new(
                                        "edge-blend.open-ring",
                                        "a ring replacement cannot have a gap",
                                    )
                                })?;
                        let connector = edit.gaps.get(&(*old, vertex)).ok_or_else(|| {
                            Refusal::new(
                                "edge-blend.missing-connection",
                                "replacement rails do not share a sewn boundary",
                            )
                        })?;
                        connector.append(model, end, start, &mut closed)?;
                    }
                    paths.push(closed);
                }
                let surface = model.geom().surface(face.surface).clone();
                let origins = model.face_origins(*old).to_vec();
                faces.push(crate::extrude::face(
                    model, surface, face.sense, paths, apex, origins,
                ));
            }
            faces.extend(
                edit.faces
                    .iter()
                    .filter_map(|(owner, face)| input_faces.contains(owner).then_some(*face)),
            );
            shells.push(model.add_shell(Shell { faces }));
        }
        lumps.push(model.add_lump(Lump { shells }));
    }
    let mut output = input;
    output.lumps = lumps;
    let output = model.add_body(output);
    let faults = model.check(output);
    if faults.is_empty() {
        Ok(output)
    } else {
        Err(Refusal::new(
            "edge-blend.invalid-topology",
            format!("{faults:?}"),
        ))
    }
}
