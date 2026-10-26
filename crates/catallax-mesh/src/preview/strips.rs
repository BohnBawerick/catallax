//! Unproved rectangular analytic bands; other trims use the shared chart CDT.
use std::collections::BTreeMap;

use catallax_numbers::{Decider, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Loop, Model, Sense};

use super::edges::Boundaries;
use crate::{
    PreviewFace, PreviewOptions, Triangle,
    chart::{self, Chart, Uv, align},
    failure,
};

type Sample = (usize, Uv);

fn chain(
    model: &Model,
    coedges: &[catallax_topology::CoedgeId],
    chart: &Chart,
    boundaries: &Boundaries,
) -> Vec<Sample> {
    let mut result: Vec<Sample> = Vec::new();
    for &coedge in coedges {
        let coedge = model.coedge(coedge);
        let vertices = &boundaries.mesh.edges[boundaries.edges[&coedge.edge]].vertices;
        for j in usize::from(!result.is_empty())..vertices.len() {
            let index = vertices[if coedge.sense == Sense::Same {
                j
            } else {
                vertices.len() - 1 - j
            }];
            let uv = chart.project(boundaries.mesh.positions[index]);
            let uv = result.last().map_or(uv, |p| align(uv, p.1, chart.periodic));
            result.push((index, uv));
        }
    }
    result
}

pub(crate) fn restrict(model: &Model, face: FaceId, chart: &mut Chart, boundaries: &Boundaries) {
    if !matches!(chart.surface, Surface::Torus(_)) || model.face(face).loops.len() != 2 {
        return;
    }
    let mut means = [0.0; 2];
    let mut increasing = [false; 2];
    for (i, &boundary) in model.face(face).loops.iter().enumerate() {
        let Loop::Coedges(coedges) = model.loop_(boundary) else {
            return;
        };
        let ring = chain(model, coedges, chart, boundaries);
        let (Some(first), Some(last)) = (ring.first(), ring.last()) else {
            return;
        };
        if ((last.1[0] - first.1[0]).abs() - 1.0).abs() > 1e-10
            || (last.1[1] - first.1[1]).abs() > 1e-10
        {
            return;
        }
        means[i] = ring.iter().map(|p| p.1[1]).sum::<f64>() / ring.len() as f64;
        increasing[i] = (last.1[0] > first.1[0]) == (model.face(face).sense == Sense::Same);
    }
    if increasing[0] == increasing[1] {
        return;
    }
    if !increasing[0] {
        means.swap(0, 1);
    }
    means[1] += libm::floor(means[0] - means[1]) + 1.0;
    chart.restrict(1, (means[0] + means[1]) * 0.5);
}

pub(crate) fn build(
    model: &Model,
    face: FaceId,
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &PreviewOptions,
    dec: &Decider,
) -> Result<bool, Refusal> {
    if chart.periodic[1]
        || !matches!(
            chart.surface,
            Surface::Cylinder(_)
                | Surface::EllipticCylinder(_)
                | Surface::Cone(_)
                | Surface::Torus(_)
        )
    {
        return Ok(false);
    }
    let definition = model.face(face);
    let mut rails = Vec::new();
    if definition.loops.len() == 2 {
        for &boundary in &definition.loops {
            let Loop::Coedges(coedges) = model.loop_(boundary) else {
                return Ok(false);
            };
            rails.push(chain(model, coedges, chart, boundaries));
        }
    } else if definition.loops.len() == 1 {
        let Loop::Coedges(coedges) = model.loop_(definition.loops[0]) else {
            return Ok(false);
        };
        if coedges.len() != 4 {
            return Ok(false);
        }
        for &coedge in coedges {
            let rail = chain(model, &[coedge], chart, boundaries);
            if rail.len() >= 2 && (rail.last().unwrap().1[0] - rail[0].1[0]).abs() > 1e-12 {
                rails.push(rail);
            }
        }
    }
    if rails.len() != 2 || rails.iter().any(|r| r.len() < 2) {
        return Ok(false);
    }
    for rail in &mut rails {
        if rail.last().unwrap().1[0] < rail[0].1[0] {
            rail.reverse();
        }
        if rail.windows(2).any(|p| p[1].1[0] <= p[0].1[0]) {
            return Ok(false);
        }
        let v = rail[0].1[1];
        if matches!(chart.surface, Surface::Torus(_))
            && rail
                .iter()
                .any(|p| (p.1[1] - v).abs() > 1e-12 * v.abs().max(1.0))
        {
            return Ok(false);
        }
    }
    if rails[0][0].1[1] > rails[1][0].1[1] {
        rails.swap(0, 1);
    }
    let lo = rails[0][0].1[1];
    let hi = rails[1][0].1[1];
    if lo >= hi {
        return Ok(false);
    }
    let closed = rails.iter().all(|r| r[0].0 == r.last().unwrap().0);
    let mut columns = Vec::new();
    if closed {
        for rail in &rails {
            columns.extend(rail.iter().map(|p| p.1[0] - libm::floor(p.1[0])));
        }
        columns.extend([0.0, 1.0]);
    } else {
        let shift = libm::round(rails[0][0].1[0] - rails[1][0].1[0]);
        for p in &mut rails[1] {
            p.1[0] += shift;
        }
        if (rails[0][0].1[0] - rails[1][0].1[0]).abs() > 1e-12
            || (rails[0].last().unwrap().1[0] - rails[1].last().unwrap().1[0]).abs() > 1e-12
        {
            return Ok(false);
        }
        for rail in &rails {
            columns.extend(rail.iter().map(|p| p.1[0]));
        }
    }
    columns.sort_unstable_by(f64::total_cmp);
    columns.dedup_by(|a, b| (*a - *b).abs() <= 16.0 / 281_474_976_710_656.0);
    let rows = if matches!(chart.surface, Surface::Torus(_)) {
        let h = chart.hessian(lo.abs().max(hi.abs()), dec)?[2].to_f64_nearest();
        let chord_step = (2.0 * options.chord_tolerance.to_f64_nearest() / h).sqrt();
        let angle_step = options
            .angular_deviation
            .as_ref()
            .map_or(f64::INFINITY, |v| {
                v.to_f64_nearest() / std::f64::consts::TAU * 0.5
            });
        let count = ((hi - lo) / chord_step.min(angle_step)).ceil().max(1.0);
        if !count.is_finite() || count > options.max_vertices as f64 {
            return Err(failure("preview band rows exceed max_vertices"));
        }
        count as usize
    } else {
        1
    };
    if columns
        .len()
        .saturating_mul(rows + 1)
        .saturating_add(boundaries.mesh.positions.len())
        > options.max_vertices
    {
        return Err(failure("preview band exceeds max_vertices"));
    }
    let mut indices = BTreeMap::new();
    for rail in &rails {
        for &(index, uv) in rail {
            indices.insert(chart::key(chart.canonical(uv)), index);
        }
    }
    let start = boundaries.mesh.triangles.len();
    let mut previous: Vec<usize> = Vec::new();
    let mut current = Vec::with_capacity(columns.len());
    for row in 0..=rows {
        current.clear();
        for &u in &columns {
            dec.check_deadline()?;
            let v_at = |rail: &[Sample]| {
                let mut u = u;
                if closed {
                    u += libm::round((rail[0].1[0] + rail.last().unwrap().1[0]) * 0.5 - u);
                }
                let next = rail
                    .partition_point(|p| p.1[0] < u)
                    .clamp(1, rail.len() - 1);
                let a = rail[next - 1].1;
                let b = rail[next].1;
                let t = ((u - a[0]) / (b[0] - a[0])).clamp(0.0, 1.0);
                a[1] + (b[1] - a[1]) * t
            };
            let lower = v_at(&rails[0]);
            let upper = v_at(&rails[1]);
            let v = lower + (upper - lower) * row as f64 / rows as f64;
            let uv = chart.canonical([u, v]);
            let key = chart::key(uv);
            let index = if let Some(&index) = indices.get(&key) {
                index
            } else {
                let point = chart.preview_position(uv, dec)?;
                let index = boundaries.mesh.positions.len();
                boundaries.mesh.positions.push(point);
                indices.insert(key, index);
                index
            };
            current.push(index);
        }
        if row > 0 {
            for (a, b) in previous.windows(2).zip(current.windows(2)) {
                for mut vertices in [[a[0], a[1], b[1]], [a[0], b[1], b[0]]] {
                    if vertices[0] == vertices[1]
                        || vertices[1] == vertices[2]
                        || vertices[0] == vertices[2]
                    {
                        continue;
                    }
                    if definition.sense == Sense::Reversed {
                        vertices.swap(1, 2);
                    }
                    let first = (0..3).min_by_key(|&i| vertices[i]).unwrap();
                    vertices.rotate_left(first);
                    boundaries.mesh.triangles.push(Triangle { vertices, face });
                }
            }
        }
        std::mem::swap(&mut previous, &mut current);
    }
    boundaries.mesh.triangles[start..].sort_unstable_by_key(|t| t.vertices);
    boundaries.mesh.faces.push(PreviewFace {
        face,
        triangles: start..boundaries.mesh.triangles.len(),
    });
    Ok(true)
}
