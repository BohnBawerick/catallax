//! Structured fallback for loopless analytic faces that exhaust the periodic CDT.

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Model, Sense};

use crate::{
    FaceMesh, MeshControl, MeshOptions, Triangle, UnavailableControl,
    chart::{Chart, chord_bound},
    edges::Boundaries,
    failure,
};

pub(crate) fn eligible(model: &Model, face: FaceId, chart: &Chart, options: &MeshOptions) -> bool {
    model.face(face).loops.is_empty()
        && options.max_edge_length.is_none()
        && matches!(chart.surface, Surface::Sphere(_) | Surface::Torus(_))
}

pub(crate) fn build(
    model: &Model,
    face: FaceId,
    chart: &Chart,
    boundaries: &mut Boundaries,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<(), Refusal> {
    let sphere = matches!(chart.surface, Surface::Sphere(_));
    let width = if sphere { 0.5 } else { 1.0 };
    let origin = if sphere { -0.25 } else { 0.0 };
    let hessian = chart.hessian(0.25, dec)?;
    // Floating estimates choose the lattice only. Every facet is certified below.
    let h = hessian.each_ref().map(|v| v.to_f64_nearest());
    let ratio = libm::sqrt(h[2] / h[0]);
    let mut dv = libm::sqrt(
        7.0 * options.chord_tolerance.to_f64_nearest()
            / (h[0] * ratio * ratio + 2.0 * h[1] * ratio + h[2]),
    );
    if let Some(angle) = &options.angular_deviation {
        dv = dv.min(angle.to_f64_nearest() / (std::f64::consts::TAU * (1.0 + ratio)));
    }
    let counts = [libm::ceil(1.0 / (ratio * dv)), libm::ceil(width / dv)];
    if counts
        .iter()
        .any(|n| !n.is_finite() || *n > options.max_vertices as f64)
    {
        return Err(failure("periodic lattice exceeds max_vertices"));
    }
    let [mut columns, mut rows] = counts.map(|n| (n as usize).max(8));
    let start = boundaries.mesh.positions.len();
    loop {
        dec.check_deadline()?;
        let count = columns
            .checked_mul(rows - usize::from(sphere))
            .and_then(|n| n.checked_add(if sphere { 2 } else { 0 }))
            .ok_or_else(|| failure("periodic lattice size overflow"))?;
        if count > options.max_vertices.saturating_sub(start) {
            return Err(failure("periodic lattice exceeds max_vertices"));
        }
        let u_at = |column| column as f64 / columns as f64;
        let v_at = |row| origin + width * row as f64 / rows as f64;
        let du = (0..columns)
            .map(|i| (u_at(i + 1) - u_at(i)).next_up())
            .fold(0.0_f64, f64::max);
        let dv = (0..rows)
            .map(|i| (v_at(i + 1) - v_at(i)).next_up())
            .fold(0.0_f64, f64::max);
        let mut residual = Rational::zero();
        let mut sample = |uv| -> Result<(), Refusal> {
            let (point, error) = chart.sample(uv, None, &options.chord_tolerance, dec)?;
            boundaries.mesh.positions.push(point);
            if error > residual {
                residual = error;
            }
            Ok(())
        };
        if sphere {
            sample([0.0, -0.25])?;
        }
        for row in usize::from(sphere)..rows {
            for column in 0..columns {
                sample([u_at(column), v_at(row)])?;
            }
        }
        if sphere {
            sample([0.0, 0.25])?;
        }
        let error = chord_bound(&[[0.0, 0.0], [du, dv]], &hessian) + residual;
        let mut triangles = Vec::with_capacity(2 * columns * (rows - usize::from(sphere)));
        let mut aspect = Rational::zero();
        let mut reached = error <= options.chord_tolerance;
        let vertex = |column: usize, row: usize| {
            if sphere && row == 0 {
                start
            } else if sphere && row == rows {
                start + count - 1
            } else {
                start
                    + usize::from(sphere)
                    + (row - usize::from(sphere)) % rows * columns
                    + column % columns
            }
        };
        'facets: for row in 0..rows {
            dec.check_deadline()?;
            let v = v_at(row);
            let next_v = v_at(row + 1);
            for column in 0..columns {
                let u = u_at(column);
                let next_u = u_at(column + 1);
                let corners = [[u, v], [next_u, v], [next_u, next_v], [u, next_v]];
                let indices = [
                    vertex(column, row),
                    vertex(column + 1, row),
                    vertex(column + 1, row + 1),
                    vertex(column, row + 1),
                ];
                for selection in [[0, 1, 2], [0, 2, 3]] {
                    let mut vertices = selection.map(|i| indices[i]);
                    if vertices[0] == vertices[1]
                        || vertices[1] == vertices[2]
                        || vertices[0] == vertices[2]
                    {
                        continue;
                    }
                    let (quality, angular) = crate::controls::quality(
                        chart,
                        &selection.map(|i| corners[i]),
                        vertices.map(|i| boundaries.mesh.positions[i]),
                        options.angular_deviation.as_ref(),
                        dec,
                    )?;
                    let Some(quality) = quality else {
                        reached = false;
                        break 'facets;
                    };
                    if !angular {
                        reached = false;
                        break 'facets;
                    }
                    aspect = aspect.max(quality);
                    if model.face(face).sense == Sense::Reversed {
                        vertices.swap(1, 2);
                    }
                    let first = (0..3).min_by_key(|&i| vertices[i]).unwrap();
                    vertices.rotate_left(first);
                    triangles.push(Triangle { vertices, face });
                }
            }
        }
        if reached {
            triangles.sort_unstable_by_key(|t| t.vertices);
            let begin = boundaries.mesh.triangles.len();
            boundaries.mesh.triangles.extend(triangles);
            boundaries.mesh.faces.push(FaceMesh {
                face,
                triangles: begin..boundaries.mesh.triangles.len(),
                max_chord_error: error,
                max_angular_deviation: options.angular_deviation.clone(),
                max_aspect_ratio: Some(aspect.clone()),
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
                        face,
                        control: MeshControl::AspectRatio,
                        reason: "the certified aspect-ratio bound exceeds the requested limit",
                    });
            }
            return Ok(());
        }
        boundaries.mesh.positions.truncate(start);
        columns = columns
            .checked_mul(2)
            .ok_or_else(|| failure("periodic lattice size overflow"))?;
        rows = rows
            .checked_mul(2)
            .ok_or_else(|| failure("periodic lattice size overflow"))?;
    }
}
