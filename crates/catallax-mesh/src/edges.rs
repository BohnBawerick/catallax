use std::collections::BTreeMap;

use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Ball, Decider, MarkLog, Rational, Refusal};
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

struct Sampler<'a> {
    piece: &'a TracePiece,
    charts: &'a [&'a Chart],
    options: &'a MeshOptions,
    dec: &'a Decider,
    positions: &'a mut Vec<[f64; 3]>,
    indices: &'a mut Vec<usize>,
    error: &'a mut Rational,
}

impl Sampler<'_> {
    fn segment(
        &mut self,
        lo: Rational,
        hi: Rational,
        first: usize,
        last: usize,
        depth: u32,
    ) -> Result<(), Refusal> {
        self.dec.check_deadline()?;
        let a = self
            .piece
            .evaluate(&Ball::from_rational(&lo, BITS), self.dec)?;
        let b = self
            .piece
            .evaluate(&Ball::from_rational(&hi, BITS), self.dec)?;
        let interval = Ball::from_rational(&lo, BITS).hull(&Ball::from_rational(&hi, BITS));
        let second = self.piece.taylor(&interval, 2, self.dec)?.pop().unwrap();
        let acceleration = second
            .iter()
            .fold(Rational::zero(), |s, x| s + x.lo().abs().max(x.hi().abs()));
        let width = &hi - &lo;
        let bound = &(&acceleration * &width) * &width * Rational::from_ratio(1, 4).unwrap()
            + chart::residual(&self.positions[first], &a)
                .max(chart::residual(&self.positions[last], &b));
        let mid = (&lo + &hi).mul_pow2(-1);
        let middle = chart::position(
            &self
                .piece
                .evaluate(&Ball::from_rational(&mid, BITS), self.dec)?,
        )?;
        let mut small = bound <= self.options.chord_tolerance.mul_pow2(-3);
        for chart in self.charts {
            let a = chart.project(self.positions[first]);
            let b = align(chart.project(self.positions[last]), a, chart.periodic);
            let m = align(chart.project(middle), a, chart.periodic);
            let hessian = chart.hessian(a[1].abs().max(b[1].abs()).max(m[1].abs()));
            small &= (0..2).all(|i| !chart.periodic[i] || (b[i] - a[i]).abs() <= 0.125);
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
        self.segment(lo, mid.clone(), first, index, depth + 1)?;
        self.segment(mid, hi, index, last, depth + 1)
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
        let first = if let Some(v) = from {
            result.vertex(model, v, dec)?
        } else {
            let p = chart::position(
                &pieces[0].evaluate(&Ball::from_rational(&Rational::zero(), BITS), dec)?,
            )?;
            let index = result.mesh.positions.len();
            result.mesh.positions.push(p);
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
            let last = if i + 1 == pieces.len() {
                end
            } else {
                let p = chart::position(
                    &piece.evaluate(&Ball::from_rational(&Rational::one(), BITS), dec)?,
                )?;
                let index = result.mesh.positions.len();
                result.mesh.positions.push(p);
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
            }
            .segment(Rational::zero(), Rational::one(), start, last, 0)?;
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
