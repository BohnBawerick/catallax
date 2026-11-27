//! Structured fallback for loopless analytic faces that exhaust the periodic CDT.

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Model, Sense};

use crate::{
    FaceMesh, MeshControl, MeshOptions, Triangle, UnavailableControl,
    assembly::FaceDraft,
    chart::{Chart, chord_bound},
    failure,
};

pub(crate) fn eligible(model: &Model, face: FaceId, chart: &Chart, options: &MeshOptions) -> bool {
    model.face(face).loops.is_empty()
        && options.max_edge_length.is_none()
        && matches!(chart.surface, Surface::Sphere(_) | Surface::Torus(_))
}

#[derive(Clone, Copy)]
struct Lattice {
    columns: usize,
    rows: usize,
    sphere: bool,
}

impl Lattice {
    fn point_count(self) -> Result<usize, Refusal> {
        self.columns
            .checked_mul(self.rows - usize::from(self.sphere))
            .and_then(|n| n.checked_add(if self.sphere { 2 } else { 0 }))
            .ok_or_else(|| failure("periodic lattice size overflow"))
    }

    fn refined(self) -> Result<Self, Refusal> {
        let double = |n: usize| {
            n.checked_mul(2)
                .ok_or_else(|| failure("periodic lattice size overflow"))
        };
        Ok(Self {
            columns: double(self.columns)?,
            rows: double(self.rows)?,
            ..self
        })
    }

    fn u(self, column: usize) -> f64 {
        column as f64 / self.columns as f64
    }

    fn v(self, row: usize) -> f64 {
        let t = row as f64 / self.rows as f64;
        if self.sphere { 0.5 * t - 0.25 } else { t }
    }

    fn index(self, column: usize, row: usize, count: usize) -> usize {
        if self.sphere && row == 0 {
            0
        } else if self.sphere && row == self.rows {
            count - 1
        } else {
            let first_row = usize::from(self.sphere);
            first_row + (row - first_row) % self.rows * self.columns + column % self.columns
        }
    }

    // Move old points back to front to retain their storage and row-major order.
    // A doubled grid visits every old station at an even index; the poles count once.
    fn sample(
        self,
        previous: Option<Self>,
        points: &mut [[f64; 3]],
        mut evaluate: impl FnMut([f64; 2]) -> Result<[f64; 3], Refusal>,
    ) -> Result<(), Refusal> {
        let count = points.len();
        if self.sphere {
            points[count - 1] = if let Some(old) = previous {
                points[old.point_count()? - 1]
            } else {
                evaluate([0.0, 0.25])?
            };
        }
        for row in (usize::from(self.sphere)..self.rows).rev() {
            for column in (0..self.columns).rev() {
                let index = self.index(column, row, count);
                points[index] = if let Some(old) = previous
                    && row.is_multiple_of(2)
                    && column.is_multiple_of(2)
                {
                    points[old.index(column / 2, row / 2, count)]
                } else {
                    evaluate([self.u(column), self.v(row)])?
                };
            }
        }
        if self.sphere && previous.is_none() {
            points[0] = evaluate([0.0, -0.25])?;
        }
        Ok(())
    }
}

pub(crate) fn build(
    model: &Model,
    face: FaceId,
    chart: &Chart,
    draft: &mut FaceDraft<'_>,
    options: &MeshOptions,
    dec: &Decider,
) -> Result<(), Refusal> {
    let sphere = matches!(chart.surface, Surface::Sphere(_));
    let width = if sphere { 0.5 } else { 1.0 };
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
    let [columns, rows] = counts.map(|n| (n as usize).max(8));
    let mut lattice = Lattice {
        columns,
        rows,
        sphere,
    };
    let start = draft.mesh.positions.len();
    let mut previous = None;
    let mut residual = Rational::zero();
    loop {
        dec.check_deadline()?;
        let count = lattice.point_count()?;
        if count > options.max_vertices.saturating_sub(start) {
            return Err(failure("periodic lattice exceeds max_vertices"));
        }
        let Lattice { columns, rows, .. } = lattice;
        let u_at = |column| lattice.u(column);
        let v_at = |row| lattice.v(row);
        let du = (0..columns)
            .map(|i| (u_at(i + 1) - u_at(i)).next_up())
            .fold(0.0_f64, f64::max);
        let dv = (0..rows)
            .map(|i| (v_at(i + 1) - v_at(i)).next_up())
            .fold(0.0_f64, f64::max);
        lattice.sample(
            previous,
            draft.mesh.positions.grow_tail(start, count),
            |uv| {
                let (point, error) = chart.sample(uv, None, &options.chord_tolerance, dec)?;
                if error > residual {
                    residual = error;
                }
                Ok(point)
            },
        )?;
        let error = &chord_bound(&[[0.0, 0.0], [du, dv]], &hessian) + &residual;
        let mut triangles = Vec::with_capacity(2 * columns * (rows - usize::from(sphere)));
        let mut aspect = Rational::zero();
        let mut reached = error <= options.chord_tolerance;
        let vertex = |column, row| start + lattice.index(column, row, count);
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
                        vertices.map(|i| draft.mesh.positions[i]),
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
            let begin = draft.mesh.triangles.len();
            draft.mesh.triangles.extend(triangles);
            draft.mesh.faces.push(FaceMesh {
                face,
                triangles: begin..draft.mesh.triangles.len(),
                max_chord_error: error,
                max_angular_deviation: options.angular_deviation.clone(),
                max_aspect_ratio: Some(aspect.clone()),
            });
            if options
                .max_aspect_ratio
                .as_ref()
                .is_some_and(|limit| aspect > *limit)
            {
                draft.mesh.unavailable_controls.push(UnavailableControl {
                    face,
                    control: MeshControl::AspectRatio,
                    reason: "the certified aspect-ratio bound exceeds the requested limit",
                });
            }
            return Ok(());
        }
        previous = Some(lattice);
        lattice = lattice.refined()?;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refined_periodic_lattice_matches_fresh_samples_without_evaluating_old_sites() {
        for sphere in [false, true] {
            let coarse = Lattice {
                columns: 9,
                rows: 11,
                sphere,
            };
            let mut points = vec![[0.0; 3]; coarse.point_count().unwrap()];
            let evaluate = |[u, v]: [f64; 2]| Ok([u, v, u - v]);
            coarse.sample(None, &mut points, evaluate).unwrap();
            let mut old = coarse;
            for _ in 0..2 {
                let fine = old.refined().unwrap();
                let new_count = fine.point_count().unwrap();
                let old_count = points.len();
                points.resize(new_count, [f64::NAN; 3]);
                let mut evaluated = 0;
                fine.sample(Some(old), &mut points, |uv| {
                    evaluated += 1;
                    evaluate(uv)
                })
                .unwrap();

                let mut fresh = vec![[f64::NAN; 3]; new_count];
                fine.sample(None, &mut fresh, evaluate).unwrap();
                assert_eq!(evaluated, new_count - old_count);
                let bits = |p: &[f64; 3]| p.map(f64::to_bits);
                assert!(points.iter().map(bits).eq(fresh.iter().map(bits)));
                old = fine;
            }
        }
    }
}
