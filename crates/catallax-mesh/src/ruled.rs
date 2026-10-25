//! Certified strips and layered cone/toroidal bands sharing their original rails.

use catallax_numbers::{Decider, Interval, Rational, Refusal};
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
    filtered_error: f64,
}

fn layers(
    lower: &[Sample],
    upper: &[Sample],
    face: (FaceId, Sense),
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<bool, Refusal> {
    let [lo, hi] = [lower[0].uv[1], upper[0].uv[1]];
    // Horizontal, disjoint rails bound one periodic rectangle. Nonconstant
    // trims retain the general triangulator rather than inferring a band.
    if lo >= hi
        || lower.iter().any(|p| p.uv[1].to_bits() != lo.to_bits())
        || upper.iter().any(|p| p.uv[1].to_bits() != hi.to_bits())
    {
        return Ok(false);
    }
    let mut columns: Vec<_> = lower[..lower.len() - 1]
        .iter()
        .chain(&upper[..upper.len() - 1])
        .map(|p| p.uv[0] - libm::floor(p.uv[0]))
        .collect();
    columns.sort_unstable_by(f64::total_cmp);
    columns.dedup_by(|a, b| (*a - *b).abs() <= 16.0 / 281_474_976_710_656.0);
    let du = columns
        .windows(2)
        .map(|p| p[1] - p[0])
        .chain(std::iter::once(
            columns[0] + 1.0 - columns[columns.len() - 1],
        ))
        .fold(0.0_f64, f64::max)
        .next_up();
    let hessian = chart.hessian(lo.abs().max(hi.abs()), dec)?;
    if chord_bound(&[[0.0, 0.0], [du, 0.0]], &hessian) > options.chord_tolerance.mul_pow2(-2) {
        return Ok(false);
    }
    let angular_step = options
        .angular_deviation
        .as_ref()
        .map(|limit| {
            catallax_numbers::Ball::pi_checked(chart::BITS, dec)
                .map(|pi| limit / &pi.hi().mul_pow2(1))
        })
        .transpose()?;
    let chord_limit = &options.chord_tolerance * &Rational::from_ratio(3, 4).unwrap();
    let mut rows = 1usize;
    loop {
        let dv = (hi - lo) / rows as f64;
        let span = if matches!(chart.surface, Surface::Torus(_)) {
            &chart::q(du) + &chart::q(dv)
        } else {
            chart::q(du)
        };
        if chord_bound(&[[0.0, 0.0], [du, dv]], &hessian) <= chord_limit
            && angular_step.as_ref().is_none_or(|step| span <= *step)
        {
            break;
        }
        let Some(next) = rows.checked_mul(2) else {
            return Ok(false);
        };
        rows = next;
        if columns.len().saturating_mul(rows) > options.max_vertices {
            return Ok(false);
        }
    }
    if columns
        .len()
        .saturating_mul(rows - 1)
        .saturating_add(boundaries.mesh.positions.len())
        > options.max_vertices
    {
        return Ok(false);
    }
    let mut previous: Option<Vec<Sample>> = None;
    for row in 1..=rows {
        dec.check_deadline()?;
        let before = previous.as_deref().unwrap_or(lower);
        if row == rows {
            return sampled_stitch(before, upper, face, chart, boundaries, options, dec);
        }
        let v = lo + (hi - lo) * row as f64 / rows as f64;
        let mut next = Vec::with_capacity(columns.len() + 1);
        for &u in &columns {
            let uv = [u, v];
            let (position, error) = chart.sample(uv, None, &options.chord_tolerance, dec)?;
            let vertex = boundaries.mesh.positions.len();
            boundaries.mesh.positions.push(position);
            next.push(Sample {
                vertex,
                uv,
                filtered_error: Interval::from_rational(&error).hi(),
                error,
            });
        }
        let mut seam = next[0].clone();
        seam.uv[0] += 1.0;
        next.push(seam);
        if !sampled_stitch(before, &next, face, chart, boundaries, options, dec)? {
            return Ok(false);
        }
        previous = Some(next);
    }
    unreachable!("a positive row count always includes the upper rail")
}

/// Returns false when the sampled trim is not a certified strip or band, or
/// needs general refinement. The periodic triangulator then handles that face.
pub(crate) fn build(
    model: &Model,
    face_id: FaceId,
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<bool, Refusal> {
    let face = model.face(face_id);
    let band = matches!(chart.surface, Surface::Cone(_) | Surface::Torus(_));
    if band && (face.loops.len() != 2 || !chart.periodic[0] || chart.periodic[1]) {
        return Ok(false);
    }
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
                    filtered_error: 0.0,
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
    for ring in [&mut lower, &mut upper] {
        for sample in ring.iter_mut() {
            let (_, error) = chart.sample(
                chart.canonical(sample.uv),
                Some(boundaries.mesh.positions[sample.vertex]),
                &options.chord_tolerance,
                dec,
            )?;
            sample.error = error;
            sample.filtered_error = Interval::from_rational(&sample.error).hi();
        }
    }
    let start = (
        boundaries.mesh.positions.len(),
        boundaries.mesh.triangles.len(),
        boundaries.mesh.faces.len(),
        boundaries.mesh.unavailable_controls.len(),
    );
    let band = matches!(chart.surface, Surface::Cone(_) | Surface::Torus(_));
    let complete = if band {
        layers(&lower, &upper, face, chart, boundaries, options, dec)?
    } else {
        sampled_stitch(&lower, &upper, face, chart, boundaries, options, dec)?
    };
    if !complete {
        boundaries.mesh.positions.truncate(start.0);
        boundaries.mesh.triangles.truncate(start.1);
        boundaries.mesh.faces.truncate(start.2);
        boundaries.mesh.unavailable_controls.truncate(start.3);
        return Ok(false);
    }
    if band {
        let faces = &mut boundaries.mesh.faces;
        let maximum = faces[start.2..]
            .iter()
            .map(|f| &f.max_chord_error)
            .max()
            .unwrap()
            .clone();
        let aspect = faces[start.2..]
            .iter()
            .filter_map(|f| f.max_aspect_ratio.as_ref())
            .max()
            .cloned();
        faces.truncate(start.2 + 1);
        let face = &mut faces[start.2];
        face.triangles = start.1..boundaries.mesh.triangles.len();
        face.max_chord_error = maximum;
        face.max_aspect_ratio = aspect;
        let unavailable = boundaries.mesh.unavailable_controls.len() > start.3;
        boundaries
            .mesh
            .unavailable_controls
            .truncate(start.3 + usize::from(unavailable));
    }
    boundaries.mesh.triangles[start.1..].sort_unstable_by_key(|triangle| triangle.vertices);
    Ok(true)
}

fn sampled_stitch(
    lower: &[Sample],
    upper: &[Sample],
    face: (FaceId, Sense),
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<bool, Refusal> {
    let (face_id, sense) = face;
    let vmax = lower
        .iter()
        .chain(upper)
        .map(|p| p.uv[1].abs())
        .fold(0.0_f64, f64::max);
    let hessian = chart.hessian(vmax, dec)?;
    let filtered_hessian = hessian.each_ref().map(Interval::from_rational);
    let filtered_tolerance = Interval::from_rational(&options.chord_tolerance);
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
        if !chart::positive_turn(uvs[0], uvs[1], uvs[2]) {
            return Ok(false);
        }
        let filtered_error = samples
            .iter()
            .map(|p| p.filtered_error)
            .fold(0.0_f64, f64::max);
        let filtered = Interval::point(filtered_error)
            .map(|error| chart::filtered_chord_bound(&uvs, &filtered_hessian) + error);
        let error = match filtered {
            Some(bound) if bound.is_finite() && bound.hi() <= filtered_tolerance.lo() => {
                chart::q(bound.hi())
            }
            _ => {
                chord_bound(&uvs, &hessian)
                    + samples.iter().map(|p| &p.error).max().unwrap().clone()
            }
        };
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
                filtered_error: 0.0,
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
