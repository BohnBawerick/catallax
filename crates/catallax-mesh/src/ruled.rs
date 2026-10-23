//! Certified axial strips for elliptic cylinders and finite cylinder fillets.

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Loop, Model, Sense};

use crate::{
    FaceMesh, MeshControl, MeshOptions, Triangle, UnavailableControl,
    chart::{self, Chart, Uv, align, chord_bound},
    edges::Boundaries,
};

#[derive(Clone)]
struct Sample {
    vertex: usize,
    uv: Uv,
    error: Rational,
}

/// Returns false when the sampled trim is not an axial strip, or needs interior
/// refinement. The general periodic triangulator then handles that face.
pub(crate) fn build(
    model: &Model,
    face_id: FaceId,
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<bool, Refusal> {
    let face = model.face(face_id);
    if face.loops.len() == 1 {
        return open_strip(model, face_id, chart, boundaries, options, dec);
    }
    if matches!(chart.surface, Surface::Cylinder(_)) {
        return Ok(false);
    }
    if face.loops.len() != 2 {
        return Ok(false);
    }
    let mut rings = [Vec::<Sample>::new(), Vec::new()];
    let mut directions = [false; 2];
    for (index, &boundary) in face.loops.iter().enumerate() {
        let Loop::Coedges(coedges) = model.loop_(boundary) else {
            return Ok(false);
        };
        let ring = &mut rings[index];
        for &coedge in coedges {
            let coedge = model.coedge(coedge);
            let edge = &boundaries.mesh.edges[boundaries.edges[&coedge.edge]];
            let count = edge.vertices.len();
            if count < 2 {
                return Ok(false);
            }
            let vertex_at = |i| {
                edge.vertices[if coedge.sense == Sense::Same {
                    i
                } else {
                    count - 1 - i
                }]
            };
            let joined = !ring.is_empty();
            if joined && ring.last().unwrap().vertex != vertex_at(0) {
                return Ok(false);
            }
            ring.reserve(count - usize::from(joined));
            for i in usize::from(joined)..count {
                let vertex = vertex_at(i);
                ring.push(Sample {
                    vertex,
                    uv: chart.canonical(chart.project(boundaries.mesh.positions[vertex])),
                    error: Rational::zero(),
                });
            }
        }
        if ring.len() < 4 || ring.first().unwrap().vertex != ring.last().unwrap().vertex {
            return Ok(false);
        }
        ring.pop();
        let mut winding = 0.0;
        let mut direction = None;
        for i in 0..ring.len() {
            let a = ring[i].uv;
            let b = align(ring[(i + 1) % ring.len()].uv, a, chart.periodic);
            let delta = b[0] - a[0];
            if delta == 0.0 || direction.is_some_and(|forward| forward != (delta > 0.0)) {
                return Ok(false);
            }
            direction = Some(delta > 0.0);
            winding += delta;
        }
        if winding.abs().to_bits() != 1.0f64.to_bits() {
            return Ok(false);
        }
        if !direction.unwrap() {
            ring.reverse();
        }
        let first = (0..ring.len())
            .min_by(|&a, &b| ring[a].uv[0].total_cmp(&ring[b].uv[0]))
            .unwrap();
        ring.rotate_left(first);
        directions[index] = (face.sense == Sense::Same) == direction.unwrap();
    }
    if directions[0] == directions[1] {
        return Ok(false);
    }
    if !directions[0] {
        rings.swap(0, 1);
    }
    let [mut lower, mut upper] = rings;
    let first = upper
        .iter()
        .rposition(|p| p.uv[0] <= lower[0].uv[0])
        .unwrap_or(upper.len() - 1);
    upper.rotate_left(first);
    if upper[0].uv[0] > lower[0].uv[0] {
        upper[0].uv[0] -= 1.0;
    }
    for i in 1..upper.len() {
        if upper[i].uv[0] <= upper[i - 1].uv[0] {
            upper[i].uv[0] += 1.0;
        }
    }
    for ring in [&mut lower, &mut upper] {
        let mut last = ring[0].clone();
        last.uv[0] += 1.0;
        ring.push(last);
    }
    stitch(
        lower,
        upper,
        (face_id, face.sense),
        chart,
        boundaries,
        options,
        dec,
    )
}

fn stitch(
    mut lower: Vec<Sample>,
    mut upper: Vec<Sample>,
    face: (FaceId, Sense),
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<bool, Refusal> {
    let (face_id, sense) = face;
    for ring in [&mut lower, &mut upper] {
        for sample in ring.iter_mut() {
            let values = chart.evaluate(chart.canonical(sample.uv), dec)?;
            sample.error = chart::residual(&boundaries.mesh.positions[sample.vertex], &values);
        }
    }
    let hessian = chart.hessian(0.0, dec)?;
    let mut triangles = Vec::with_capacity(lower.len() + upper.len() - 2);
    let mut maximum = Rational::zero();
    let mut aspect = Rational::zero();
    let (mut a, mut b) = (0, 0);
    while a + 1 < lower.len() || b + 1 < upper.len() {
        dec.check_deadline()?;
        let samples = if b + 1 == upper.len()
            || (a + 1 < lower.len() && lower[a + 1].uv[0] <= upper[b + 1].uv[0])
        {
            let samples = [&lower[a], &lower[a + 1], &upper[b]];
            a += 1;
            samples
        } else {
            let samples = [&lower[a], &upper[b + 1], &upper[b]];
            b += 1;
            samples
        };
        let uvs = samples.map(|p| p.uv);
        let [p, q, r] = uvs.map(|uv| uv.map(chart::q));
        let orientation = (&q[0] - &p[0]) * (&r[1] - &p[1]) - (&q[1] - &p[1]) * (&r[0] - &p[0]);
        if orientation <= Rational::zero() {
            return Ok(false);
        }
        let error =
            chord_bound(&uvs, &hessian) + samples.iter().map(|p| &p.error).max().unwrap().clone();
        if error > options.chord_tolerance {
            return Ok(false);
        }
        let mut vertices = samples.map(|p| p.vertex);
        let positions = vertices.map(|i| boundaries.mesh.positions[i]);
        if let Some(limit) = &options.max_edge_length {
            for i in 0..3 {
                if crate::controls::length_squared(positions[i], positions[(i + 1) % 3])
                    > limit * limit
                {
                    return Ok(false);
                }
            }
        }
        let (ratio, angular_reached) = crate::controls::quality(
            chart,
            &uvs,
            positions,
            options.angular_deviation.as_ref(),
            dec,
        )?;
        let Some(ratio) = ratio else {
            return Ok(false);
        };
        if !angular_reached {
            return Ok(false);
        }
        maximum = maximum.max(error);
        aspect = aspect.max(ratio);
        if sense == Sense::Reversed {
            vertices.swap(1, 2);
        }
        let first = (0..3).min_by_key(|&i| vertices[i]).unwrap();
        vertices.rotate_left(first);
        triangles.push(Triangle {
            vertices,
            face: face_id,
        });
    }
    triangles.sort_unstable_by_key(|triangle| triangle.vertices);
    let start = boundaries.mesh.triangles.len();
    boundaries.mesh.triangles.extend(triangles);
    boundaries.mesh.faces.push(FaceMesh {
        face: face_id,
        triangles: start..boundaries.mesh.triangles.len(),
        max_chord_error: maximum,
        max_aspect_ratio: Some(aspect.clone()),
        max_angular_deviation: options.angular_deviation.clone(),
    });
    if options
        .max_aspect_ratio
        .as_ref()
        .is_some_and(|limit| aspect > *limit)
    {
        boundaries
            .mesh
            .unavailable_controls
            .push(UnavailableControl {
                face: face_id,
                control: MeshControl::AspectRatio,
                reason: "the certified aspect-ratio bound exceeds the requested limit",
            });
    }
    Ok(true)
}

fn open_strip(
    model: &Model,
    face_id: FaceId,
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<bool, Refusal> {
    let face = model.face(face_id);
    let Loop::Coedges(coedges) = model.loop_(face.loops[0]) else {
        return Ok(false);
    };
    if coedges.len() != 4 {
        return Ok(false);
    }
    let mut chains: [Vec<Sample>; 4] = std::array::from_fn(|_| Vec::new());
    let mut axial = [false; 4];
    for (i, &coedge) in coedges.iter().enumerate() {
        let coedge = model.coedge(coedge);
        let vertices = &boundaries.mesh.edges[boundaries.edges[&coedge.edge]].vertices;
        if vertices.len() < 2 {
            return Ok(false);
        }
        for j in 0..vertices.len() {
            let j = if coedge.sense == Sense::Same {
                j
            } else {
                vertices.len() - 1 - j
            };
            let vertex = vertices[j];
            let uv = chart.canonical(chart.project(boundaries.mesh.positions[vertex]));
            chains[i].push(Sample {
                vertex,
                uv,
                error: Rational::zero(),
            });
        }
        axial[i] = chains[i]
            .iter()
            .all(|p| p.uv[0].to_bits() == chains[i][0].uv[0].to_bits());
        if axial[i] && chains[i].len() != 2 {
            return Ok(false);
        }
    }
    let caps: Vec<_> = (0..4).filter(|&i| !axial[i]).collect();
    if caps.len() != 2 || (caps[0] + 2) % 4 != caps[1] {
        return Ok(false);
    }
    let mut lower = None;
    let mut upper = None;
    for index in caps {
        let mut chain = std::mem::take(&mut chains[index]);
        let mut direction = None;
        for i in 1..chain.len() {
            chain[i].uv = align(chain[i].uv, chain[i - 1].uv, chart.periodic);
            let delta = chain[i].uv[0] - chain[i - 1].uv[0];
            if delta == 0.0 || direction.is_some_and(|forward| forward != (delta > 0.0)) {
                return Ok(false);
            }
            direction = Some(delta > 0.0);
        }
        let forward = direction.unwrap();
        if !forward {
            chain.reverse();
        }
        if forward == (face.sense == Sense::Same) {
            lower = Some(chain);
        } else {
            upper = Some(chain);
        }
    }
    let (Some(lower), Some(mut upper)) = (lower, upper) else {
        return Ok(false);
    };
    let shift = libm::round(lower[0].uv[0] - upper[0].uv[0]);
    for sample in &mut upper {
        sample.uv[0] += shift;
    }
    if lower[0].uv[0].to_bits() != upper[0].uv[0].to_bits()
        || lower.last().unwrap().uv[0].to_bits() != upper.last().unwrap().uv[0].to_bits()
    {
        return Ok(false);
    }
    stitch(
        lower,
        upper,
        (face_id, face.sense),
        chart,
        boundaries,
        options,
        dec,
    )
}
