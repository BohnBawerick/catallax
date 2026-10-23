use std::collections::BTreeMap;

use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Ball, Decider, MarkLog, Rational, Refusal};
use catallax_shapes::Point3;
use catallax_topology::{BodyId, EdgeId, Loop, Model, Sense, VertexId};

use crate::{
    EdgeMesh, Mesh, MeshOptions,
    chart::{self, BITS, Chart, align, chord_bound},
    failure,
};

pub(crate) struct Boundaries {
    pub(crate) mesh: Mesh,
    pub(crate) edges: BTreeMap<EdgeId, usize>,
    pub(crate) vertices: BTreeMap<VertexId, usize>,
    pub(crate) exact_positions: BTreeMap<usize, Point3>,
    pub(crate) marks: MarkLog,
}

impl Boundaries {
    pub(crate) fn vertex(
        &mut self,
        model: &Model,
        id: VertexId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if let Some(&v) = self.vertices.get(&id) {
            return Ok(v);
        }
        let point = model.geom().position(model.vertex(id).point);
        let xyz = [&point.x, &point.y, &point.z]
            .map(|x| dec.enclose_at(x, BITS))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        let p = chart::position(&xyz.try_into().unwrap())?;
        let index = self.mesh.positions.len();
        self.mesh.positions.push(p);
        self.vertices.insert(id, index);
        Ok(index)
    }
}

pub(crate) struct ChordProof {
    pub(crate) error: Rational,
    pub(crate) replacement: Rational,
    pub(crate) endpoints: [[Ball; 3]; 2],
}

/// Taylor's full-interval chord bound, including actual output-point replacement.
pub(crate) fn trace_chord(
    piece: &TracePiece,
    interval: [&Rational; 2],
    positions: [&[f64; 3]; 2],
    dec: &Decider,
) -> Result<ChordProof, Refusal> {
    let [lo, hi] = interval;
    let a = piece.evaluate(&Ball::from_rational(lo, BITS), dec)?;
    let b = piece.evaluate(&Ball::from_rational(hi, BITS), dec)?;
    let interval = Ball::from_rational(lo, BITS).hull(&Ball::from_rational(hi, BITS));
    let second = piece.taylor(&interval, 2, dec)?.pop().unwrap();
    let acceleration = second
        .iter()
        .fold(Rational::zero(), |sum, x| sum + x.lo().abs().max(x.hi().abs()));
    let width = hi - lo;
    let replacement = chart::residual(positions[0], &a).max(chart::residual(positions[1], &b));
    let curvature = &(&acceleration * &width) * &width * Rational::from_ratio(1, 4).unwrap();
    let error = &curvature + &replacement;
    Ok(ChordProof { error, replacement, endpoints: [a, b] })
}

struct Sampler<'a> {
    piece: &'a TracePiece,
    charts: &'a [&'a Chart],
    options: &'a MeshOptions,
    dec: &'a Decider,
    positions: &'a mut Vec<[f64; 3]>,
    indices: &'a mut Vec<usize>,
    error: &'a mut Rational,
    exact_positions: Option<&'a mut BTreeMap<usize, Point3>>,
}

impl Sampler<'_> {
    fn segment(
        &mut self,
        lo: Rational,
        hi: Rational,
        first: usize,
        last: usize,
        depth: u32,
        endpoints: [&[Ball; 3]; 2],
    ) -> Result<(), Refusal> {
        self.dec.check_deadline()?;
        let [a, b] = endpoints;
        let interval = Ball::from_rational(&lo, BITS).hull(&Ball::from_rational(&hi, BITS));
        let second = self.piece.taylor(&interval, 2, self.dec)?.pop().unwrap();
        let acceleration = second
            .iter()
            .fold(Rational::zero(), |s, x| s + x.lo().abs().max(x.hi().abs()));
        let width = &hi - &lo;
        let bound = &(&acceleration * &width) * &width * Rational::from_ratio(1, 4).unwrap()
            + chart::residual(&self.positions[first], a)
                .max(chart::residual(&self.positions[last], b));
        let mid = (&lo + &hi).mul_pow2(-1);
        let middle_values = self
            .piece
            .evaluate(&Ball::from_rational(&mid, BITS), self.dec)?;
        let middle = chart::position(&middle_values)?;
        let mut small = bound <= self.options.chord_tolerance.mul_pow2(-3);
        if let Some(limit) = &self.options.max_edge_length {
            small &= crate::controls::length_squared(self.positions[first], self.positions[last])
                <= limit * limit;
        }
        for chart in self.charts {
            let a = chart.project(self.positions[first]);
            let b = align(chart.project(self.positions[last]), a, chart.periodic);
            let m = align(chart.project(middle), a, chart.periodic);
            let hessian = chart.hessian(a[1].abs().max(b[1].abs()).max(m[1].abs()), self.dec)?;
            small &= (0..2).all(|i| !chart.periodic[i] || (b[i] - a[i]).abs() <= 0.125);
            if let Some(angle) = &self.options.angular_deviation {
                let mut span = (crate::chart::q(b[0]) - crate::chart::q(a[0])).abs();
                if matches!(
                    chart.surface,
                    catallax_shapes::Surface::Sphere(_) | catallax_shapes::Surface::Torus(_)
                ) {
                    span = span + (crate::chart::q(b[1]) - crate::chart::q(a[1])).abs();
                }
                if !matches!(chart.surface, catallax_shapes::Surface::Plane(_)) {
                    small &= span * Rational::from_i64(16) <= *angle;
                }
            }
            small &= chord_bound(&[a, b, m], &hessian) <= self.options.chord_tolerance.mul_pow2(-3);
        }
        if small {
            *self.error = self.error.clone().max(bound);
            self.indices.push(last);
            return Ok(());
        }
        if depth >= 48 || self.positions.len() >= self.options.max_vertices {
            return Err(failure(
                "edge refinement exhausted the deterministic vertex/precision limit",
            ));
        }
        let index = self.positions.len();
        self.positions.push(middle);
        if let Some(points) = &mut self.exact_positions {
            points.insert(index, self.piece.evaluate_exact(&mid, self.dec)?);
        }
        self.segment(
            lo,
            mid.clone(),
            first,
            index,
            depth + 1,
            [a, &middle_values],
        )?;
        self.segment(mid, hi, index, last, depth + 1, [&middle_values, b])
    }
}

pub(crate) fn sample(
    model: &Model,
    body: BodyId,
    charts: &BTreeMap<catallax_topology::FaceId, Chart>,
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<Boundaries, Refusal> {
    let mut result = Boundaries {
        mesh: Mesh::default(),
        edges: BTreeMap::new(),
        vertices: BTreeMap::new(),
        exact_positions: BTreeMap::new(),
        marks: model.body(body).marks.clone(),
    };
    let mut adjacent = BTreeMap::<EdgeId, Vec<&Chart>>::new();
    for (&face, chart) in charts {
        for &l in &model.face(face).loops {
            match model.loop_(l) {
                Loop::Coedges(cs) => {
                    for &c in cs {
                        adjacent
                            .entry(model.coedge(c).edge)
                            .or_default()
                            .push(chart);
                    }
                }
                Loop::Vertex(v) => {
                    result.vertex(model, *v, dec)?;
                }
            }
        }
    }
    for edge_id in model.edges_of_body(body) {
        let edge = model.edge(edge_id);
        let (from, to) = if edge.sense == Sense::Same {
            (edge.start, edge.end)
        } else {
            (edge.end, edge.start)
        };
        let traced = trace(
            model.geom(),
            edge.curve,
            from.map(|v| model.vertex(v).point),
            to.map(|v| model.vertex(v).point),
            dec,
        )?;
        result.marks.extend(&traced.marks);
        let pieces = traced.value.pieces();
        let initial_values = pieces
            .first()
            .map(|piece| piece.evaluate(&Ball::from_rational(&Rational::zero(), BITS), dec))
            .transpose()?;
        let planar = adjacent[&edge_id]
            .iter()
            .any(|chart| matches!(chart.surface, catallax_shapes::Surface::Plane(_)));
        let first = if let Some(v) = from {
            result.vertex(model, v, dec)?
        } else {
            let p = chart::position(initial_values.as_ref().unwrap())?;
            let index = result.mesh.positions.len();
            result.mesh.positions.push(p);
            if planar {
                result
                    .exact_positions
                    .insert(index, pieces[0].evaluate_exact(&Rational::zero(), dec)?);
            }
            index
        };
        let end = if let Some(v) = to {
            result.vertex(model, v, dec)?
        } else {
            first
        };
        let mut indices = vec![first];
        let mut error = Rational::zero();
        for (i, piece) in pieces.iter().enumerate() {
            let next_values;
            let first_values = if i == 0 {
                initial_values.as_ref().unwrap()
            } else {
                next_values = piece.evaluate(&Ball::from_rational(&Rational::zero(), BITS), dec)?;
                &next_values
            };
            let last_values = piece.evaluate(&Ball::from_rational(&Rational::one(), BITS), dec)?;
            let last = if i + 1 == pieces.len() {
                end
            } else {
                let p = chart::position(&last_values)?;
                let index = result.mesh.positions.len();
                result.mesh.positions.push(p);
                if planar {
                    result
                        .exact_positions
                        .insert(index, piece.evaluate_exact(&Rational::one(), dec)?);
                }
                index
            };
            let start = *indices.last().unwrap();
            Sampler {
                piece,
                charts: &adjacent[&edge_id],
                options,
                dec,
                positions: &mut result.mesh.positions,
                indices: &mut indices,
                error: &mut error,
                exact_positions: planar.then_some(&mut result.exact_positions),
            }
            .segment(
                Rational::zero(),
                Rational::one(),
                start,
                last,
                0,
                [first_values, &last_values],
            )?;
        }
        // Chart seams can add a trace sample closer to an endpoint than the
        // periodic display grid can distinguish. Remove only interior mesh
        // samples, retaining the exact edge and both topological endpoints.
        let mut i = 1;
        while i < indices.len() {
            let (a, b) = (indices[i - 1], indices[i]);
            let collapses = adjacent[&edge_id].iter().any(|chart| {
                chart::key(chart.canonical(chart.project(result.mesh.positions[a])))
                    == chart::key(chart.canonical(chart.project(result.mesh.positions[b])))
            });
            if collapses && !(a == first && b == end) {
                let displacement = (0..3).fold(Rational::zero(), |sum, k| {
                    sum + (chart::q(result.mesh.positions[a][k])
                        - chart::q(result.mesh.positions[b][k]))
                    .abs()
                });
                let bound = &error + &displacement;
                if bound <= options.chord_tolerance.mul_pow2(-3) {
                    error = bound;
                    if b == end {
                        indices.remove(i - 1);
                        i = i.saturating_sub(1).max(1);
                    } else {
                        indices.remove(i);
                    }
                    continue;
                }
            }
            i += 1;
        }
        if edge.sense == Sense::Reversed {
            indices.reverse();
        }
        result.edges.insert(edge_id, result.mesh.edges.len());
        result.mesh.edges.push(EdgeMesh {
            edge: edge_id,
            vertices: indices,
            max_chord_error: error,
        });
    }
    Ok(result)
}
