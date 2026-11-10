//! Exact chart expansion with logical picking keys and physical edge polylines.

use super::*;
use catallax_topology::ChartExpansion;
use std::collections::BTreeMap;

pub(super) fn build(
    model: &Model,
    expansion: &ChartExpansion,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Answered<Mesh> {
    catallax_numbers::reuse_operation_enclosures();
    let mut answer = mesh_body(model, expansion.body(), options, dec)?;
    let mesh = &mut answer.value;
    let mut faces: BTreeMap<FaceId, FaceMesh> = BTreeMap::new();
    for mut face in std::mem::take(&mut mesh.faces) {
        face.face = expansion.face(face.face);
        let count = face.triangles.len();
        match faces.entry(face.face) {
            std::collections::btree_map::Entry::Vacant(entry) => {
                face.triangles = 0..count;
                entry.insert(face);
            }
            std::collections::btree_map::Entry::Occupied(mut entry) => {
                let old = entry.get_mut();
                old.triangles.end += count;
                old.max_chord_error = old.max_chord_error.clone().max(face.max_chord_error);
                old.max_aspect_ratio = old
                    .max_aspect_ratio
                    .take()
                    .zip(face.max_aspect_ratio)
                    .map(|(a, b)| a.max(b));
                old.max_angular_deviation = old
                    .max_angular_deviation
                    .take()
                    .zip(face.max_angular_deviation)
                    .map(|(a, b)| a.max(b));
            }
        }
    }
    for triangle in &mut mesh.triangles {
        triangle.face = expansion.face(triangle.face);
    }
    mesh.triangles
        .sort_unstable_by_key(|triangle| (triangle.face, triangle.vertices));
    let mut offset = 0;
    for mut face in faces.into_values() {
        let count = face.triangles.len();
        face.triangles = offset..offset + count;
        offset += count;
        mesh.faces.push(face);
    }
    let mut pieces: BTreeMap<_, _> = std::mem::take(&mut mesh.edges)
        .into_iter()
        .map(|edge| (edge.edge, edge))
        .collect();
    for (logical, chain) in expansion.edges() {
        let count = chain
            .iter()
            .map(|(edge, _)| {
                pieces
                    .get(edge)
                    .map(|p| p.vertices.len())
                    .ok_or_else(|| failure("a physical chart boundary has no certified polyline"))
            })
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .sum::<usize>();
        let mut vertices = Vec::with_capacity(count.saturating_sub(chain.len().saturating_sub(1)));
        let mut error = Rational::zero();
        for (edge, forward) in chain {
            let mut piece = pieces
                .remove(edge)
                .ok_or_else(|| failure("a physical chart boundary is consumed twice"))?;
            if !forward {
                piece.vertices.reverse();
            }
            if !vertices.is_empty() && vertices.last() != piece.vertices.first() {
                return Err(failure(
                    "logical edge pieces do not share their exact sewn mesh vertex",
                ));
            }
            let skip = usize::from(!vertices.is_empty());
            vertices.extend(piece.vertices.into_iter().skip(skip));
            error = error.max(piece.max_chord_error);
        }
        if model.edge(*logical).start.is_none() && vertices.first() != vertices.last() {
            return Err(failure(
                "a logical ring edge has an open certified polyline",
            ));
        }
        mesh.edges.push(EdgeMesh {
            edge: *logical,
            vertices,
            max_chord_error: error,
        });
    }
    for control in &mut mesh.unavailable_controls {
        control.face = expansion.face(control.face);
    }
    let key = |control: &UnavailableControl| {
        (
            control.face,
            match control.control {
                MeshControl::AngularDeviation => 0,
                MeshControl::EdgeLength => 1,
                MeshControl::AspectRatio => 2,
            },
            control.reason,
        )
    };
    mesh.unavailable_controls.sort_unstable_by_key(key);
    mesh.unavailable_controls
        .dedup_by_key(|control| key(control));
    Ok(answer)
}
