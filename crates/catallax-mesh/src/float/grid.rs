use std::{collections::BTreeMap, ops::Range};

use catallax_numbers::{Decider, Refusal};
use catallax_topology::{BodyId, Model};

use super::{
    DisplayDescription, DisplayFace, DisplayOptions, DisplaySurface, Uv, Xyz, describe_body,
};
use super::{evaluate::Projector, trim};
use crate::{PreviewEdge, PreviewFace, PreviewMesh, PreviewOptions, Triangle};

fn limit() -> Refusal {
    Refusal::new(
        "mesh.preview-resource-limit",
        "the sampled preview exceeds max_vertices",
    )
}

fn trim_error(error: trim::TrimRefusal) -> Refusal {
    match error {
        trim::TrimRefusal::MeshSize => limit(),
        _ => Refusal::new(
            "mesh.preview-grid",
            format!("cannot triangulate the sampled trim: {error:?}"),
        ),
    }
}

fn append<F, E>(
    mesh: &mut PreviewMesh<F, E>,
    position: Xyz,
    max_vertices: usize,
) -> Result<usize, Refusal> {
    if mesh.positions.len() == max_vertices {
        return Err(limit());
    }
    if !position.iter().all(|value| value.is_finite()) {
        return Err(Refusal::new(
            "mesh.preview-grid",
            "sample grid exceeds finite coordinates",
        ));
    }
    let index = mesh.positions.len();
    mesh.positions.push(position);
    Ok(index)
}

pub(crate) fn build(
    model: &Model,
    body: BodyId,
    options: &PreviewOptions,
    divisions: u32,
    tolerance: f64,
    dec: &mut Decider,
) -> Result<PreviewMesh, Refusal> {
    dec.check_deadline()?;
    let edge_count = model.edges_of_body(body).len();
    dec.charge_traversal(edge_count as u64)?;
    let description = describe_body(
        model,
        body,
        &DisplayOptions {
            chord_tolerance: tolerance,
            max_edge_points: options
                .max_vertices
                .saturating_add(edge_count.saturating_mul(2)),
        },
    )?;
    from_description(description, options, divisions, dec)
}

pub(crate) fn from_description<F: Clone, E>(
    description: DisplayDescription<F, E>,
    options: &PreviewOptions,
    divisions: u32,
    dec: &mut Decider,
) -> Result<PreviewMesh<F, E>, Refusal> {
    let world = |p: Xyz| std::array::from_fn(|axis| p[axis] + description.origin[axis]);
    let mut mesh = PreviewMesh::default();
    let mut vertices = BTreeMap::new();
    let mut edge_samples = Vec::new();
    let mut ranges = Vec::with_capacity(description.edges.len());
    for edge in description.edges {
        dec.charge_traversal(edge.points.len() as u64)?;
        let first = u32::try_from(edge_samples.len()).map_err(|_| limit())?;
        let closed = edge
            .points
            .first()
            .zip(edge.points.last())
            .is_some_and(|(a, b)| a.map(f64::to_bits) == b.map(f64::to_bits));
        let count = edge.points.len();
        let mut indices = Vec::with_capacity(count);
        for (i, point) in edge.points.into_iter().enumerate() {
            let endpoint = !closed && (i == 0 || i + 1 == count);
            let position = world(point);
            let index = if endpoint {
                match vertices.entry(point.map(f64::to_bits)) {
                    std::collections::btree_map::Entry::Occupied(e) => *e.get(),
                    std::collections::btree_map::Entry::Vacant(e) => {
                        *e.insert(append(&mut mesh, position, options.max_vertices)?)
                    }
                }
            } else if i + 1 == count && closed {
                indices[0]
            } else {
                append(&mut mesh, position, options.max_vertices)?
            };
            indices.push(index);
            edge_samples.push(index);
        }
        ranges.push(first..u32::try_from(edge_samples.len()).map_err(|_| limit())?);
        mesh.edges.push(PreviewEdge {
            edge: edge.edge,
            vertices: indices,
        });
    }
    for face in &description.faces {
        dec.check_deadline()?;
        let mut projector = Projector::new(&face.surface, None)
            .ok_or_else(|| Refusal::new("mesh.preview-grid", "sampled face has no evaluator"))?;
        let start = mesh.triangles.len();
        if face.loops.iter().all(|lp| lp.uses.is_empty()) {
            untrimmed(
                face,
                &mut projector,
                divisions,
                &world,
                &mut mesh,
                options,
                dec,
            )?;
        } else {
            trimmed(
                face,
                &mut projector,
                &ranges,
                &edge_samples,
                divisions,
                &world,
                &mut mesh,
                options,
                dec,
            )?;
        }
        if start == mesh.triangles.len() {
            return Err(Refusal::new(
                "mesh.preview-grid",
                "sampled face has no triangles",
            ));
        }
        mesh.triangles[start..].sort_unstable_by_key(|triangle| triangle.vertices);
        mesh.faces.push(PreviewFace {
            face: face.face.clone(),
            triangles: start..mesh.triangles.len(),
        });
    }
    Ok(mesh)
}

#[allow(clippy::too_many_arguments)]
fn trimmed<F: Clone, E>(
    face: &DisplayFace<F>,
    projector: &mut Projector<'_>,
    ranges: &[Range<u32>],
    edge_samples: &[usize],
    divisions: u32,
    world: &impl Fn(Xyz) -> Xyz,
    mesh: &mut PreviewMesh<F, E>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let nodes = face
        .loops
        .iter()
        .flat_map(|lp| &lp.uses)
        .map(|edge| edge.samples.len())
        .sum::<usize>();
    dec.charge_traversal((nodes as u64).saturating_add(u64::from(divisions).pow(2)))?;
    let trim = trim::triangulate(
        face,
        &|i| ranges.get(i).cloned(),
        &|a, b| edge_samples[a as usize] == edge_samples[b as usize],
        [divisions; 2],
        options.max_vertices.saturating_add(nodes),
        usize::MAX,
        true,
    )
    .map_err(trim_error)?;
    let mut indices = Vec::with_capacity(trim.vertices.len());
    for node in trim.vertices {
        let index = if node.boundary != trim::NO_BOUNDARY {
            edge_samples[node.boundary as usize]
        } else {
            let uv: Uv = std::array::from_fn(|axis| {
                let [lo, hi] = trim.bounds[axis];
                lo + node.uv[axis] * (hi - lo)
            });
            let p = projector.jet(canonical(face, uv)).0;
            append(mesh, world(p), options.max_vertices)?
        };
        indices.push(index);
    }
    for &triangle in trim.indices.as_chunks::<3>().0 {
        let vertices = triangle.map(|i| indices[i as usize]);
        emit(mesh, face, vertices);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn untrimmed<F: Clone, E>(
    face: &DisplayFace<F>,
    projector: &mut Projector<'_>,
    divisions: u32,
    world: &impl Fn(Xyz) -> Xyz,
    mesh: &mut PreviewMesh<F, E>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let divisions = divisions.max(3) as usize;
    let width = divisions + 1;
    dec.charge_traversal((width * width) as u64)?;
    let mut grid = Vec::with_capacity(width * width);
    for v in 0..width {
        for u in 0..width {
            let at = [u, v];
            let uv = std::array::from_fn(|axis| {
                let [lo, hi] = face.bounds[axis];
                lo + at[axis] as f64 / divisions as f64 * (hi - lo)
            });
            let vertex = if u == divisions && face.periods[0].is_some() {
                grid[v * width]
            } else if v == divisions && face.periods[1].is_some() {
                grid[u]
            } else if u > 0
                && (v == 0 || v == divisions)
                && matches!(face.surface, DisplaySurface::Sphere { .. })
            {
                grid[v * width]
            } else {
                let p = projector.jet(uv).0;
                append(mesh, world(p), options.max_vertices)?
            };
            grid.push(vertex);
        }
    }
    for v in 0..divisions {
        for u in 0..divisions {
            let a = v * width + u;
            for corners in [[a, a + 1, a + width + 1], [a, a + width + 1, a + width]] {
                let mut vertices = corners.map(|i| grid[i]);
                if face.reversed {
                    vertices.swap(1, 2);
                }
                emit(mesh, face, vertices);
            }
        }
    }
    Ok(())
}

fn emit<F: Clone, E>(
    mesh: &mut PreviewMesh<F, E>,
    face: &DisplayFace<F>,
    mut vertices: [usize; 3],
) {
    if vertices[0] == vertices[1] || vertices[1] == vertices[2] || vertices[2] == vertices[0] {
        return;
    }
    let first = (0..3).min_by_key(|&i| vertices[i]).unwrap();
    vertices.rotate_left(first);
    mesh.triangles.push(Triangle {
        vertices,
        face: face.face.clone(),
    });
}

fn canonical<F>(face: &DisplayFace<F>, mut uv: Uv) -> Uv {
    for (axis, value) in uv.iter_mut().enumerate() {
        if let Some(period) = face.periods[axis] {
            let lo = face.bounds[axis][0];
            *value = lo + (*value - lo).rem_euclid(period);
        }
    }
    uv
}
