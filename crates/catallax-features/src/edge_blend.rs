//! Shared boundary replacement for elementary edge blends and bevels.

use std::collections::BTreeMap;

use catallax_numbers::Refusal;
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Lump, Model, Sense, Shell, VertexId};

/// Replacements refer to the immutable input skin. A boundary shared by two
/// output faces has one edge handle; only its two coedges differ.
#[derive(Default)]
pub(crate) struct Replacement {
    pub edges: BTreeMap<(FaceId, EdgeId), EdgeId>,
    pub gaps: BTreeMap<(FaceId, VertexId), EdgeId>,
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
    let mut lumps = Vec::with_capacity(input.lumps.len());
    for lump in &input.lumps {
        let input_shells = model.lump(*lump).shells.clone();
        let mut shells = Vec::with_capacity(input_shells.len());
        for shell in input_shells {
            let input_faces = model.shell(shell).faces.clone();
            let mut faces = Vec::with_capacity(input_faces.len() + edit.faces.len());
            for old in &input_faces {
                if !edit.edges.keys().any(|(face, _)| face == old)
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
                    let mut path = Vec::with_capacity(cs.len() + edit.gaps.len());
                    for c in cs {
                        let c = model.coedge(*c);
                        let edge = edit.edges.get(&(*old, c.edge)).copied().unwrap_or(c.edge);
                        path.push((edge, c.sense));
                    }
                    let mut closed = Vec::with_capacity(path.len() + edit.gaps.len());
                    for i in 0..path.len() {
                        let current = path[i];
                        closed.push(current);
                        let next = path[(i + 1) % path.len()];
                        let end = ends(model, current.0, current.1).1;
                        let start = ends(model, next.0, next.1).0;
                        if end == start {
                            continue;
                        }
                        let original = model.coedge(cs[i]);
                        let vertex =
                            ends(model, original.edge, original.sense)
                                .1
                                .ok_or_else(|| {
                                    Refusal::new(
                                        "edge-blend.open-ring",
                                        "a ring replacement cannot have a gap",
                                    )
                                })?;
                        let connector = *edit.gaps.get(&(*old, vertex)).ok_or_else(|| {
                            Refusal::new(
                                "edge-blend.missing-connection",
                                "replacement rails do not share a sewn boundary",
                            )
                        })?;
                        let sense = if ends(model, connector, Sense::Same) == (end, start) {
                            Sense::Same
                        } else if ends(model, connector, Sense::Reversed) == (end, start) {
                            Sense::Reversed
                        } else {
                            return Err(Refusal::new(
                                "edge-blend.wrong-connection",
                                "connector does not join the two replacement rails",
                            ));
                        };
                        closed.push((connector, sense));
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
