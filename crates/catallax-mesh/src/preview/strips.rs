//! Unproved rectangular analytic bands; other trims use the shared chart CDT.
use catallax_numbers::{Decider, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Loop, Model, Sense};

use super::edges::Boundaries;
use crate::{
    PreviewFace, PreviewOptions, Triangle,
    chart::{Chart, Uv, align},
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

// A consistently oriented, closed lift of a single trim selects a disk,
// not the complementary region. Don't tile the empty part of the carrier.
fn restrict_disk(chart: &mut Chart, ring: &[Sample], sense: Sense) {
    let Some(&(first, origin)) = ring.first() else {
        return;
    };
    if ring.len() < 4
        || ring.last() != Some(&(first, origin))
        || ring.iter().any(|(_, uv)| uv.iter().any(|v| !v.is_finite()))
    {
        return;
    }
    let mut bounds = origin.map(|v| [v, v]);
    let mut area = catallax_numbers::Interval::point(0.0).unwrap();
    let origin = origin.map(|v| catallax_numbers::Interval::point(v).unwrap());
    for pair in ring.windows(2) {
        let [a, b] = [pair[0], pair[1]].map(|(_, uv)| {
            std::array::from_fn::<_, 2, _>(|i| {
                catallax_numbers::Interval::point(uv[i]).unwrap() - origin[i]
            })
        });
        area = area + a[0] * b[1] - a[1] * b[0];
        for (i, [lo, hi]) in bounds.iter_mut().enumerate() {
            *lo = lo.min(pair[1].1[i]);
            *hi = hi.max(pair[1].1[i]);
        }
    }
    if (sense == Sense::Same && area.lo() <= 0.0)
        || (sense == Sense::Reversed && area.hi() >= 0.0)
        || bounds.iter().any(|[lo, hi]| hi - lo >= 1.0 || hi <= lo)
    {
        return;
    }
    for (i, [lo, hi]) in bounds.into_iter().enumerate() {
        if chart.periodic[i] {
            chart.restrict(i, lo.midpoint(hi));
        }
    }
}

pub(crate) fn restrict(model: &Model, face: FaceId, chart: &mut Chart, boundaries: &Boundaries) {
    let definition = model.face(face);
    if matches!(chart.surface, Surface::Sphere(_) | Surface::Torus(_))
        && let [boundary] = definition.loops.as_slice()
        && let Loop::Coedges(coedges) = model.loop_(*boundary)
    {
        restrict_disk(
            chart,
            &chain(model, coedges, chart, boundaries),
            definition.sense,
        );
        return;
    }
    if !matches!(chart.surface, Surface::Torus(_)) || definition.loops.len() != 2 {
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

// Advance one rail at a time: a station on one rail never inserts a
// vertex into a shared boundary on its neighbour.
fn stitch(
    lower: &[Sample],
    upper: &[Sample],
    face: FaceId,
    sense: Sense,
    triangles: &mut Vec<Triangle>,
    dec: &Decider,
) -> Result<(), Refusal> {
    let (mut i, mut j) = (0, 0);
    while i + 1 < lower.len() || j + 1 < upper.len() {
        dec.check_deadline()?;
        let mut vertices = if j + 1 == upper.len()
            || i + 1 < lower.len() && lower[i + 1].1[0] < upper[j + 1].1[0]
        {
            i += 1;
            [lower[i - 1].0, lower[i].0, upper[j].0]
        } else {
            j += 1;
            [lower[i].0, upper[j].0, upper[j - 1].0]
        };
        if vertices[0] == vertices[1] || vertices[1] == vertices[2] || vertices[0] == vertices[2] {
            continue;
        }
        if sense == Sense::Reversed {
            vertices.swap(1, 2);
        }
        let first = (0..3).min_by_key(|&k| vertices[k]).unwrap();
        vertices.rotate_left(first);
        triangles.push(Triangle { vertices, face });
    }
    Ok(())
}

pub(crate) fn build(
    model: &Model,
    face: FaceId,
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &PreviewOptions,
    dec: &Decider,
) -> Result<bool, Refusal> {
    // A disk's short sides are sampled edges too; the ring strip only
    // takes boundary rows, not side columns.
    if chart.periodic[1]
        || matches!(chart.surface, Surface::Torus(_)) && !chart.periodic[0]
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
        if coedges.len() < 4 {
            return Ok(false);
        }
        // Rails are the runs of curved coedges between rulings. A run may join
        // several edges, such as a slot end wall's top arc and the blend that cuts it.
        let curved: Vec<_> = coedges
            .iter()
            .map(|&coedge| {
                let rail = chain(model, &[coedge], chart, boundaries);
                rail.len() >= 2 && (rail.last().unwrap().1[0] - rail[0].1[0]).abs() > 1e-12
            })
            .collect();
        let Some(ruling) = curved.iter().position(|curved| !curved) else {
            return Ok(false);
        };
        let mut sides = coedges.iter().zip(&curved).filter(|(_, curved)| !**curved);
        for _ in 0..2 {
            let Some((&side, _)) = sides.next() else {
                return Ok(false);
            };
            let edge = model.coedge(side).edge;
            if boundaries.mesh.edges[boundaries.edges[&edge]]
                .vertices
                .len()
                != 2
            {
                return Ok(false);
            }
        }
        if sides.next().is_some() {
            return Ok(false);
        }
        let mut run = Vec::new();
        for i in 1..=coedges.len() {
            let index = (ruling + i) % coedges.len();
            if curved[index] {
                run.push(coedges[index]);
            } else if !run.is_empty() {
                rails.push(chain(model, &run, chart, boundaries));
                run.clear();
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
    if closed {
        for rail in &mut rails {
            if (rail.last().unwrap().1[0] - rail[0].1[0] - 1.0).abs() > 1e-12 {
                return Ok(false);
            }
            rail.pop();
            let first_turn = libm::floor(rail[0].1[0]);
            let seam = rail.partition_point(|p| p.1[0] < first_turn + 1.0);
            for (_, uv) in rail.iter_mut() {
                uv[0] -= libm::floor(uv[0]);
            }
            rail.rotate_left(seam);
            let (index, [u, v]) = rail[0];
            rail.push((index, [u + 1.0, v]));
        }
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
    }
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
    let start = boundaries.mesh.triangles.len();
    let mut previous = Vec::new();
    if rows > 1 {
        let mut columns = Vec::new();
        for rail in &rails {
            columns.extend(
                rail[..rail.len() - usize::from(closed)]
                    .iter()
                    .map(|p| p.1[0]),
            );
        }
        columns.sort_unstable_by(f64::total_cmp);
        columns.dedup_by(|a, b| (*a - *b).abs() <= 16.0 / 281_474_976_710_656.0);
        if columns
            .len()
            .saturating_mul(rows - 1)
            .saturating_add(boundaries.mesh.positions.len())
            > options.max_vertices
        {
            return Err(failure("preview band exceeds max_vertices"));
        }
        let mut current = Vec::with_capacity(columns.len() + usize::from(closed));
        for row in 1..rows {
            let v = lo + (hi - lo) * row as f64 / rows as f64;
            for &u in &columns {
                let point = chart.preview_position(chart.canonical([u, v]), dec)?;
                let index = boundaries.mesh.positions.len();
                boundaries.mesh.positions.push(point);
                current.push((index, [u, v]));
            }
            if closed {
                let (index, [u, v]) = current[0];
                current.push((index, [u + 1.0, v]));
            }
            let lower = if row == 1 { &rails[0] } else { &previous };
            stitch(
                lower,
                &current,
                face,
                definition.sense,
                &mut boundaries.mesh.triangles,
                dec,
            )?;
            std::mem::swap(&mut previous, &mut current);
            current.clear();
        }
    }
    let lower = if rows == 1 { &rails[0] } else { &previous };
    stitch(
        lower,
        &rails[1],
        face,
        definition.sense,
        &mut boundaries.mesh.triangles,
        dec,
    )?;
    boundaries.mesh.triangles[start..].sort_unstable_by_key(|t| t.vertices);
    boundaries.mesh.faces.push(PreviewFace {
        face,
        triangles: start..boundaries.mesh.triangles.len(),
    });
    Ok(true)
}
