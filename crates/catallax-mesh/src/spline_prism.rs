//! A ruled strip and two convex cap fans, sampled on one common parameter grid.

use crate::{EdgeMesh, FaceMesh, Mesh, MeshOptions, Triangle, chart, edges::trace_chord, failure};
use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Answer, Answered, Ball, Decider, Rational, Real, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, FaceId, Model, Sense, SplinePrism};

fn cross(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
fn oriented(
    positions: &[[f64; 3]],
    triangle: [usize; 3],
    normal: &[Ball; 3],
) -> Result<(), Refusal> {
    let points = triangle
        .map(|index| positions[index].map(|x| Ball::from_rational(&chart::q(x), chart::BITS)));
    let a = std::array::from_fn(|i| &points[1][i] - &points[0][i]);
    let b = std::array::from_fn(|i| &points[2][i] - &points[0][i]);
    let actual = cross(&a, &b);
    let dot = &(&actual[0] * &normal[0]) + &(&actual[1] * &normal[1]) + &actual[2] * &normal[2];
    if dot.lo() <= Rational::zero() {
        return Err(failure(
            "binary64 triangle material orientation is not separated from zero",
        ));
    }
    Ok(())
}

struct Sampler<'a> {
    options: &'a MeshOptions,
    dec: &'a Decider,
    limit: Rational,
    axis: [Ball; 3],
    mesh: Mesh,
    rings: [Vec<usize>; 2],
    errors: [Rational; 2],
    side: Vec<[usize; 3]>,
}
impl Sampler<'_> {
    fn pair(
        &mut self,
        pieces: [&TracePiece; 2],
        parameter: &Rational,
    ) -> Result<[usize; 2], Refusal> {
        if self.mesh.positions.len().saturating_add(4) > self.options.max_vertices {
            return Err(failure("paired spline sampling exceeds the vertex limit"));
        }
        let mut indices = [0; 2];
        for index in 0..2 {
            let position =
                pieces[index].evaluate(&Ball::from_rational(parameter, chart::BITS), self.dec)?;
            indices[index] = self.mesh.positions.len();
            self.mesh.positions.push(chart::position(&position)?);
        }
        Ok(indices)
    }
    fn segment(
        &mut self,
        pieces: [&TracePiece; 2],
        interval: [Rational; 2],
        endpoints: [[usize; 2]; 2],
        depth: u32,
    ) -> Result<(), Refusal> {
        self.dec.check_deadline()?;
        let [lo, hi] = &interval;
        let [first, last] = endpoints;
        let proofs = [
            trace_chord(
                pieces[0],
                [lo, hi],
                [
                    &self.mesh.positions[first[0]],
                    &self.mesh.positions[last[0]],
                ],
                self.dec,
            )?,
            trace_chord(
                pieces[1],
                [lo, hi],
                [
                    &self.mesh.positions[first[1]],
                    &self.mesh.positions[last[1]],
                ],
                self.dec,
            )?,
        ];
        if proofs.iter().any(|proof| proof.replacement >= self.limit) {
            return Err(failure(
                "binary64 endpoint rounding exceeds the fixed mesh/feature margin",
            ));
        }
        if proofs.iter().all(|proof| proof.error <= self.limit) {
            let chord =
                std::array::from_fn(|i| &proofs[0].endpoints[1][i] - &proofs[0].endpoints[0][i]);
            let normal = cross(&chord, &self.axis);
            let triangles = [[first[0], last[0], last[1]], [first[0], last[1], first[1]]];
            for triangle in triangles {
                oriented(&self.mesh.positions, triangle, &normal)?;
                self.side.push(triangle);
            }
            for index in 0..2 {
                self.rings[index].push(last[index]);
                self.errors[index] = self.errors[index].clone().max(proofs[index].error.clone());
            }
            return Ok(());
        }
        if depth >= 48 {
            return Err(failure(
                "spline display subdivision reached its deterministic precision limit",
            ));
        }
        let mid = (lo + hi).mul_pow2(-1);
        let middle = self.pair(pieces, &mid)?;
        self.segment(
            pieces,
            [lo.clone(), mid.clone()],
            [first, middle],
            depth + 1,
        )?;
        self.segment(pieces, [mid, hi.clone()], [middle, last], depth + 1)
    }
}

pub(crate) fn build(
    model: &Model,
    body: BodyId,
    face_order: &[FaceId],
    options: &MeshOptions,
    dec: &mut Decider,
) -> Answered<Mesh> {
    let checked = SplinePrism::check(model, body, dec)?;
    let spline = checked.spline();
    let Surface::EllipticCylinder(original) = spline.definition().original() else {
        unreachable!("checked elliptic sweep");
    };
    let radius = dec
        .enclose_at(
            &(&original.minor_sq.sqrt()
                - spline.definition().distance()
                - Real::from_rational(spline.bound().upper.clone())),
            chart::BITS,
        )?
        .lo();
    let height = dec.enclose_at(spline.height(), chart::BITS)?.lo();
    if radius <= Rational::zero() || height <= Rational::zero() {
        return Err(failure(
            "display needs a positive proved interior/height margin",
        ));
    }
    let feature_margin = radius.min(height) / Rational::from_i64(16);
    let limit = options.chord_tolerance.mul_pow2(-2).min(feature_margin);
    let axis = [
        dec.enclose_at(&spline.axis().x, chart::BITS)?,
        dec.enclose_at(&spline.axis().y, chart::BITS)?,
        dec.enclose_at(&spline.axis().z, chart::BITS)?,
    ];
    let rings = checked.rings();
    let lower = trace(model.geom(), model.edge(rings[0]).curve, None, None, dec)?;
    let upper = trace(model.geom(), model.edge(rings[1]).curve, None, None, dec)?;
    if lower.value.pieces().len() != upper.value.pieces().len() || lower.value.pieces().is_empty() {
        return Err(failure(
            "shared spline boundaries have different finite span covers",
        ));
    }
    let mut marks = model.body(body).marks.clone();
    marks.extend(&lower.marks);
    marks.extend(&upper.marks);
    let mut sampler = Sampler {
        options,
        dec,
        limit,
        axis,
        mesh: Mesh::default(),
        rings: [Vec::new(), Vec::new()],
        errors: [Rational::zero(), Rational::zero()],
        side: Vec::new(),
    };
    let first = sampler.pair(
        [&lower.value.pieces()[0], &upper.value.pieces()[0]],
        &Rational::zero(),
    )?;
    for index in 0..2 {
        sampler.rings[index].push(first[index]);
    }
    let mut previous = first;
    let span_count = lower.value.pieces().len();
    for (index, (lower, upper)) in lower
        .value
        .pieces()
        .iter()
        .zip(upper.value.pieces())
        .enumerate()
    {
        let pieces = [lower, upper];
        let last = if index + 1 == span_count {
            first
        } else {
            sampler.pair(pieces, &Rational::one())?
        };
        sampler.segment(
            pieces,
            [Rational::zero(), Rational::one()],
            [previous, last],
            0,
        )?;
        previous = last;
    }
    let mut centers = [0; 2];
    let mut center_errors = [Rational::zero(), Rational::zero()];
    for index in 0..2 {
        let center = original
            .origin
            .plus(&spline.axis().scaled(spline.axial()[index]));
        let enclosure = [
            sampler.dec.enclose_at(&center.x, chart::BITS)?,
            sampler.dec.enclose_at(&center.y, chart::BITS)?,
            sampler.dec.enclose_at(&center.z, chart::BITS)?,
        ];
        let position = chart::position(&enclosure)?;
        let error = chart::residual(&position, &enclosure);
        if error >= sampler.limit {
            return Err(failure(
                "cap center rounding exceeds the mesh/feature margin",
            ));
        }
        centers[index] = sampler.mesh.positions.len();
        sampler.mesh.positions.push(position);
        center_errors[index] = error;
    }
    let faces = checked.faces();
    for &face in face_order {
        let slot = faces
            .iter()
            .position(|actual| *actual == face)
            .ok_or_else(|| failure("display face inventory differs from checked prism"))?;
        let start = sampler.mesh.triangles.len();
        let error = if slot < 2 {
            let normal = if slot == 0 {
                sampler.axis.each_ref().map(|x| -x)
            } else {
                sampler.axis.clone()
            };
            for pair in sampler.rings[slot].windows(2) {
                let triangle = if slot == 0 {
                    [centers[slot], pair[1], pair[0]]
                } else {
                    [centers[slot], pair[0], pair[1]]
                };
                oriented(&sampler.mesh.positions, triangle, &normal)?;
                sampler.mesh.triangles.push(Triangle {
                    vertices: triangle,
                    face,
                });
            }
            sampler.errors[slot]
                .clone()
                .max(center_errors[slot].clone())
        } else {
            for &vertices in &sampler.side {
                sampler.mesh.triangles.push(Triangle { vertices, face });
            }
            sampler.errors[0].clone().max(sampler.errors[1].clone())
        };
        if error > options.chord_tolerance {
            return Err(failure(
                "proved spline mesh error exceeds the requested output accuracy",
            ));
        }
        sampler.mesh.faces.push(FaceMesh {
            face,
            triangles: start..sampler.mesh.triangles.len(),
            max_chord_error: error,
        });
    }
    for (index, vertices) in sampler.rings.into_iter().enumerate() {
        let mut vertices = vertices;
        if model.edge(rings[index]).sense == Sense::Reversed {
            vertices.reverse();
        }
        sampler.mesh.edges.push(EdgeMesh {
            edge: rings[index],
            vertices,
            max_chord_error: sampler.errors[index].clone(),
        });
    }
    crate::triangulate::check_closed(&sampler.mesh)?;
    Ok(Answer {
        value: sampler.mesh,
        marks,
    })
}
