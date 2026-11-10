//! Export exact sewn cells while reporting the published logical model keys.

use super::*;
use std::collections::BTreeMap;

pub(super) fn write<W: Write>(
    model: &Model,
    bodies: &[BodyId],
    options: &StepOptions,
    dec: &mut Decider,
    out: &mut W,
) -> Answered<StepReport> {
    let mut expanded = Vec::with_capacity(bodies.len());
    let mut faces = BTreeMap::new();
    let mut edges = BTreeMap::new();
    let mut added_edges = 0;
    let mut added_vertices = 0;
    for body in bodies {
        if let Some(cover) = model.chart_expansion(*body) {
            expanded.push(cover.body());
            for face in model.faces_of_body(cover.body()) {
                faces.insert(face, cover.face(face));
            }
            for (logical, pieces) in cover.edges() {
                for (piece, _) in pieces {
                    edges.insert(*piece, *logical);
                }
            }
            let source = model.counts(*body);
            let cells = model.counts(cover.body());
            added_edges += cells.edges - source.edges;
            added_vertices += cells.vertices - source.vertices;
        } else {
            expanded.push(*body);
        }
    }
    let mut answer = write_prepared(model, &expanded, options, dec, out, true)?;
    answer.value.added_edges += added_edges;
    answer.value.added_vertices += added_vertices;
    for region in &mut answer.value.face_map {
        if let Some(logical) = faces.get(&region.face) {
            region.face = *logical;
        }
    }
    for repair in &mut answer.value.regularizations {
        if let Some(logical) = faces.get(&repair.face) {
            repair.face = *logical;
        }
        for face in &mut repair.merged_faces {
            if let Some(logical) = faces.get(face) {
                *face = *logical;
            }
        }
        for edge in &mut repair.edges {
            if let Some(logical) = edges.get(edge) {
                *edge = *logical;
            }
        }
    }
    Ok(answer)
}
