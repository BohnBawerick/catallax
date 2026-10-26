//! Sampled shared rings for the existing retained spline-prism family.
use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Ball, Decider, Interval, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, Model, Sense, SplinePrism};

use crate::{PreviewEdge, PreviewFace, PreviewMesh, PreviewOptions, Triangle, chart, failure};

type Pair = [[f64; 3]; 2];

struct Sampler<'a> {
    options: &'a PreviewOptions,
    dec: &'a Decider,
    mesh: PreviewMesh,
    rings: [Vec<usize>; 2],
}

impl Sampler<'_> {
    fn pair(&mut self, points: Pair) -> Result<[usize; 2], Refusal> {
        if self.mesh.positions.len().saturating_add(2) > self.options.max_vertices {
            return Err(failure("spline preview reached max_vertices"));
        }
        let first = self.mesh.positions.len();
        self.mesh.positions.extend(points);
        Ok([first, first + 1])
    }

    fn segment(
        &mut self,
        pieces: [&TracePiece; 2],
        interval: [f64; 2],
        ends: [[usize; 2]; 2],
        depth: u32,
    ) -> Result<(), Refusal> {
        self.dec.check_deadline()?;
        let [lo, hi] = interval;
        let [first, last] = ends;
        let mid = (lo + hi) * 0.5;
        let [a, b, c] = [lo + (hi - lo) * 0.25, mid, lo + (hi - lo) * 0.75]
            .map(|t| sample(pieces, t, self.dec));
        let points = [a?, b?, c?];
        let mut small = true;
        let tolerance = self.options.chord_tolerance.to_f64_nearest();
        for ring in 0..2 {
            let a = self.mesh.positions[first[ring]];
            let b = self.mesh.positions[last[ring]];
            for (i, p) in points.iter().enumerate() {
                let t = (i + 1) as f64 * 0.25;
                let distance = (0..3)
                    .map(|k| (p[ring][k] - (a[k] + (b[k] - a[k]) * t)).powi(2))
                    .sum::<f64>()
                    .sqrt();
                small &= distance <= tolerance * 0.5;
            }
            if let Some(angle) = &self.options.angular_deviation {
                let a: [f64; 3] = std::array::from_fn(|i| points[0][ring][i] - a[i]);
                let b: [f64; 3] = std::array::from_fn(|i| b[i] - points[2][ring][i]);
                let dot = (0..3).map(|i| a[i] * b[i]).sum::<f64>();
                let norm = (a.iter().map(|x| x * x).sum::<f64>()
                    * b.iter().map(|x| x * x).sum::<f64>())
                .sqrt();
                small &=
                    norm > 0.0 && (dot / norm).clamp(-1.0, 1.0).acos() <= angle.to_f64_nearest();
            }
        }
        if small {
            for (ring, vertex) in last.into_iter().enumerate() {
                self.rings[ring].push(vertex);
            }
            return Ok(());
        }
        if depth >= 32 {
            return Err(failure(
                "spline preview sampling did not resolve its visual targets",
            ));
        }
        let middle = self.pair(points[1])?;
        self.segment(pieces, [lo, mid], [first, middle], depth + 1)?;
        self.segment(pieces, [mid, hi], [middle, last], depth + 1)
    }
}

fn sample(pieces: [&TracePiece; 2], t: f64, dec: &Decider) -> Result<Pair, Refusal> {
    let point = |piece: &TracePiece| {
        if let Some(jets) = piece.taylor_interval(Interval::point(t).unwrap(), 0, dec)? {
            let values = jets[0].map(|x| x.mid());
            if values.iter().all(|x| x.is_finite()) {
                return Ok(values);
            }
        }
        chart::position(&piece.evaluate(&Ball::from_rational(&chart::q(t), chart::BITS), dec)?)
    };
    Ok([point(pieces[0])?, point(pieces[1])?])
}

pub(super) fn build(
    model: &Model,
    body: BodyId,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<PreviewMesh, Refusal> {
    let checked = SplinePrism::check(model, body, dec)?;
    let spline = checked.spline();
    let Surface::EllipticCylinder(original) = spline.definition().original() else {
        unreachable!("checked spline prism")
    };
    let rings = checked.rings();
    let lower = trace(model.geom(), model.edge(rings[0]).curve, None, None, dec)?;
    let upper = trace(model.geom(), model.edge(rings[1]).curve, None, None, dec)?;
    if lower.value.pieces().len() != upper.value.pieces().len() || lower.value.pieces().is_empty() {
        return Err(failure("spline preview boundary spans do not correspond"));
    }
    let mut sampler = Sampler {
        options,
        dec,
        mesh: PreviewMesh::default(),
        rings: [Vec::new(), Vec::new()],
    };
    let first = sampler.pair(sample(
        [&lower.value.pieces()[0], &upper.value.pieces()[0]],
        0.0,
        dec,
    )?)?;
    for (ring, vertex) in first.into_iter().enumerate() {
        sampler.rings[ring].push(vertex);
    }
    let mut before = first;
    for (index, (a, b)) in lower
        .value
        .pieces()
        .iter()
        .zip(upper.value.pieces())
        .enumerate()
    {
        let last = if index + 1 == lower.value.pieces().len() {
            first
        } else {
            sampler.pair(sample([a, b], 1.0, dec)?)?
        };
        sampler.segment([a, b], [0.0, 1.0], [before, last], 0)?;
        before = last;
    }
    let mut centers = [0; 2];
    for (index, center_index) in centers.iter_mut().enumerate() {
        if sampler.mesh.positions.len() >= options.max_vertices {
            return Err(failure("spline preview cap reached max_vertices"));
        }
        let center = original
            .origin
            .plus(&spline.axis().scaled(spline.axial()[index]));
        let position = chart::position(&[
            dec.enclose_at(&center.x, chart::BITS)?,
            dec.enclose_at(&center.y, chart::BITS)?,
            dec.enclose_at(&center.z, chart::BITS)?,
        ])?;
        *center_index = sampler.mesh.positions.len();
        sampler.mesh.positions.push(position);
    }
    let faces = checked.faces();
    for face in model.faces_of_body(body) {
        let slot = faces
            .iter()
            .position(|&id| id == face)
            .ok_or_else(|| failure("spline preview face inventory differs"))?;
        let start = sampler.mesh.triangles.len();
        if slot < 2 {
            for pair in sampler.rings[slot].windows(2) {
                let vertices = if slot == 0 {
                    [centers[slot], pair[1], pair[0]]
                } else {
                    [centers[slot], pair[0], pair[1]]
                };
                sampler.mesh.triangles.push(Triangle { vertices, face });
            }
        } else {
            for (a, b) in sampler.rings[0].windows(2).zip(sampler.rings[1].windows(2)) {
                sampler.mesh.triangles.extend([
                    Triangle {
                        vertices: [a[0], a[1], b[1]],
                        face,
                    },
                    Triangle {
                        vertices: [a[0], b[1], b[0]],
                        face,
                    },
                ]);
            }
        }
        sampler.mesh.faces.push(PreviewFace {
            face,
            triangles: start..sampler.mesh.triangles.len(),
        });
    }
    for (index, mut vertices) in sampler.rings.into_iter().enumerate() {
        if model.edge(rings[index]).sense == Sense::Reversed {
            vertices.reverse();
        }
        sampler.mesh.edges.push(PreviewEdge {
            edge: rings[index],
            vertices,
        });
    }
    Ok(sampler.mesh)
}
