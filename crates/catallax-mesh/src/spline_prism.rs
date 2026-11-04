//! A ruled strip and two convex cap fans, sampled on one common parameter grid.

use crate::{EdgeMesh, FaceMesh, Mesh, MeshOptions, Triangle, chart, edges::trace_chord, failure};
use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Answer, Answered, Ball, Decider, Interval, Rational, Real, Refusal};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, FaceId, Model, Sense, SplinePrism};

fn cross(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
pub(super) fn oriented(
    positions: &[[f64; 3]],
    triangle: [usize; 3],
    normal: &[Ball; 3],
) -> Result<(Rational, Rational), Refusal> {
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
    let area2 = (&(&actual[0].square() + &actual[1].square()) + &actual[2].square())
        .sqrt()
        .ok_or_else(|| failure("rounded triangle area is not enclosed"))?;
    if area2.lo() <= Rational::zero() {
        return Err(failure(
            "rounded triangle quality exceeds certificate precision",
        ));
    }
    let longest = (0..3)
        .map(|i| {
            crate::controls::length_squared(
                positions[triangle[i]],
                positions[triangle[(i + 1) % 3]],
            )
        })
        .max()
        .unwrap();
    Ok((&longest / &area2.lo(), longest))
}

fn balls(values: [Interval; 3]) -> [Ball; 3] {
    values.map(|value| {
        Ball::from_rational(&chart::q(value.lo()), chart::BITS)
            .hull(&Ball::from_rational(&chart::q(value.hi()), chart::BITS))
    })
}

fn point(piece: &TracePiece, parameter: &Rational, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    if let Some(jets) = piece.taylor_interval(Interval::from_rational(parameter), 0, dec)?
        && jets[0].iter().all(Interval::is_finite)
    {
        return Ok(balls(jets[0]));
    }
    piece.evaluate(&Ball::from_rational(parameter, chart::BITS), dec)
}

fn span(
    piece: &TracePiece,
    interval: [&Rational; 2],
    positions: [&[f64; 3]; 2],
    limit: &Rational,
    need_tangent: bool,
    dec: &Decider,
) -> Result<(crate::edges::ChordProof, Option<[Ball; 3]>), Refusal> {
    let [lo, hi] = interval;
    let domain = Interval::from_rational(lo).hull(&Interval::from_rational(hi));
    if let Some(jets) = piece.taylor_interval(domain, 2, dec)?
        && jets.iter().flatten().all(Interval::is_finite)
    {
        let endpoints = [point(piece, lo, dec)?, point(piece, hi, dec)?];
        let replacement = chart::residual(positions[0], &endpoints[0])
            .max(chart::residual(positions[1], &endpoints[1]));
        let acceleration = jets[2]
            .iter()
            .fold(Interval::point(0.0).unwrap(), |sum, x| {
                sum + Interval::point(x.lo().abs().max(x.hi().abs())).unwrap()
            });
        if acceleration.is_finite() && replacement < *limit {
            let width = hi - lo;
            let curvature = &(&chart::q(acceleration.hi()) * &width)
                * &width
                * Rational::from_ratio(1, 4).unwrap();
            let error = &curvature + &replacement;
            return Ok((
                crate::edges::ChordProof {
                    error,
                    replacement,
                    endpoints,
                },
                need_tangent.then(|| balls(jets[1])),
            ));
        }
    }
    let proof = trace_chord(piece, interval, positions, dec)?;
    let tangent = if need_tangent {
        let domain =
            Ball::from_rational(lo, chart::BITS).hull(&Ball::from_rational(hi, chart::BITS));
        Some(piece.derivative(&domain, dec)?)
    } else {
        None
    };
    Ok((proof, tangent))
}

// The positive dot selects the acute branch. The squared cross product then
// bounds sin(angle) over the entire normal box, including facet rounding.
fn angular_within(
    positions: &[[f64; 3]],
    triangle: [usize; 3],
    normal: &[Ball; 3],
    sine_squared: &Rational,
) -> bool {
    let points = triangle
        .map(|index| positions[index].map(|x| Ball::from_rational(&chart::q(x), chart::BITS)));
    let a = std::array::from_fn(|i| &points[1][i] - &points[0][i]);
    let b = std::array::from_fn(|i| &points[2][i] - &points[0][i]);
    let facet = cross(&a, &b);
    let zero = Ball::from_rational(&Rational::zero(), chart::BITS);
    let dot = (0..3).fold(zero.clone(), |sum, i| sum + &facet[i] * &normal[i]);
    if dot.lo() <= Rational::zero() {
        return false;
    }
    let norm_squared =
        |vector: &[Ball; 3]| vector.iter().fold(zero.clone(), |sum, x| sum + x.square());
    let transverse = norm_squared(&cross(&facet, normal));
    let scale = norm_squared(&facet).lo() * norm_squared(normal).lo();
    scale > Rational::zero() && transverse.hi() <= &scale * sine_squared
}

struct Sampler<'a> {
    options: &'a MeshOptions,
    dec: &'a Decider,
    limit: Rational,
    axis: [Ball; 3],
    sine_squared: Option<Rational>,
    mesh: Mesh,
    rings: [Vec<usize>; 2],
    errors: [Rational; 2],
    side: Vec<[usize; 3]>,
    side_aspect: Rational,
    max_length_squared: Rational,
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
            let position = point(pieces[index], parameter, self.dec)?;
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
            span(
                pieces[0],
                [lo, hi],
                [
                    &self.mesh.positions[first[0]],
                    &self.mesh.positions[last[0]],
                ],
                &self.limit,
                self.sine_squared.is_some(),
                self.dec,
            )?,
            span(
                pieces[1],
                [lo, hi],
                [
                    &self.mesh.positions[first[1]],
                    &self.mesh.positions[last[1]],
                ],
                &self.limit,
                false,
                self.dec,
            )?,
        ];
        if proofs
            .iter()
            .any(|(proof, _)| proof.replacement >= self.limit)
        {
            return Err(failure(
                "binary64 endpoint rounding exceeds the fixed mesh/feature margin",
            ));
        }
        let chord =
            std::array::from_fn(|i| &proofs[0].0.endpoints[1][i] - &proofs[0].0.endpoints[0][i]);
        let normal = cross(&chord, &self.axis);
        let triangles = [[first[0], last[0], last[1]], [first[0], last[1], first[1]]];
        let angular_reached = proofs.iter().all(|(proof, _)| proof.error <= self.limit)
            && self.sine_squared.as_ref().is_none_or(|limit| {
                // A checked prism is an exact ruled translation, so its normal is
                // independent of height and the lower ring covers the whole side.
                let normals = cross(proofs[0].1.as_ref().unwrap(), &self.axis);
                triangles.iter().all(|&triangle| {
                    angular_within(&self.mesh.positions, triangle, &normals, limit)
                })
            });
        if angular_reached {
            for triangle in triangles {
                let (aspect, length_squared) = oriented(&self.mesh.positions, triangle, &normal)?;
                if aspect > self.side_aspect {
                    self.side_aspect = aspect;
                }
                if self.options.max_edge_length.is_some()
                    && length_squared > self.max_length_squared
                {
                    self.max_length_squared = length_squared;
                }
                self.side.push(triangle);
            }
            for (index, (proof, _)) in proofs.into_iter().enumerate() {
                self.rings[index].push(last[index]);
                if proof.error > self.errors[index] {
                    self.errors[index] = proof.error;
                }
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
    // Certifying a smaller acute angle also satisfies larger requests.
    let sine_squared = options.angular_deviation.as_ref().map(|angle| {
        let sine = Ball::from_rational(&angle.clone().min(Rational::one()), chart::BITS)
            .sin()
            .lo()
            .max(Rational::zero());
        &sine * &sine
    });
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
        sine_squared,
        mesh: Mesh::default(),
        rings: [Vec::new(), Vec::new()],
        errors: [Rational::zero(), Rational::zero()],
        side: Vec::new(),
        side_aspect: Rational::zero(),
        max_length_squared: Rational::zero(),
    };
    let first = sampler.pair(
        [&lower.value.pieces()[0], &upper.value.pieces()[0]],
        &Rational::zero(),
    )?;
    for (index, &vertex) in first.iter().enumerate() {
        sampler.rings[index].push(vertex);
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
        let mut max_aspect_ratio = if slot == 2 {
            std::mem::take(&mut sampler.side_aspect)
        } else {
            Rational::zero()
        };
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
                let (aspect, length_squared) =
                    oriented(&sampler.mesh.positions, triangle, &normal)?;
                if let Some(limit) = &sampler.sine_squared
                    && !angular_within(&sampler.mesh.positions, triangle, &normal, limit)
                {
                    return Err(failure(
                        "cap rounding exceeds the requested angular deviation",
                    ));
                }
                max_aspect_ratio = max_aspect_ratio.max(aspect);
                if options.max_edge_length.is_some() && length_squared > sampler.max_length_squared
                {
                    sampler.max_length_squared = length_squared;
                }
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
        if options
            .max_aspect_ratio
            .as_ref()
            .is_some_and(|limit| max_aspect_ratio > *limit)
        {
            sampler
                .mesh
                .unavailable_controls
                .push(crate::UnavailableControl {
                    face,
                    control: crate::MeshControl::AspectRatio,
                    reason: "the certified spline-facet aspect bound exceeds the requested limit",
                });
        }
        sampler.mesh.faces.push(FaceMesh {
            face,
            triangles: start..sampler.mesh.triangles.len(),
            max_chord_error: error,
            max_aspect_ratio: Some(max_aspect_ratio),
            max_angular_deviation: options.angular_deviation.clone(),
        });
    }
    if options
        .max_edge_length
        .as_ref()
        .is_some_and(|limit| sampler.max_length_squared > limit * limit)
    {
        return Err(Refusal::new(
            "mesh.spline-edge-length",
            "the shared cap/strip grid needs additional interior refinement to meet the requested edge length",
        ));
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn angular_certificate_covers_normal_variation_and_material_sense() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let ball = |value| Ball::from_rational(&Rational::parse(value).unwrap(), chart::BITS);
        let normal = [ball("-1/4").hull(&ball("1/4")), ball("0"), ball("1")];
        let tight = Rational::parse("1/100").unwrap();
        let loose = Rational::parse("1/4").unwrap();
        assert!(!angular_within(&positions, [0, 1, 2], &normal, &tight));
        assert!(angular_within(&positions, [0, 1, 2], &normal, &loose));
        assert!(!angular_within(&positions, [0, 2, 1], &normal, &loose));
        assert!(!angular_within(&positions, [0, 0, 1], &normal, &loose));
    }
}
