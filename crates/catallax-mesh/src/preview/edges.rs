//! Shared, unproved edge polylines. Sampling never constructs chord certificates.
use std::collections::BTreeMap;

use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Ball, Decider, Interval, Refusal};
use catallax_topology::{BodyId, EdgeId, FaceId, Model, Sense, VertexId};

use crate::{
    PreviewEdge, PreviewMesh, PreviewOptions,
    chart::{self, Chart, align},
    failure,
};

pub(crate) struct Boundaries {
    pub(crate) mesh: PreviewMesh,
    pub(crate) edges: BTreeMap<EdgeId, usize>,
    pub(crate) vertices: BTreeMap<VertexId, usize>,
}

impl Boundaries {
    fn push(&mut self, point: [f64; 3], options: &PreviewOptions) -> Result<usize, Refusal> {
        if self.mesh.positions.len() >= options.max_vertices {
            return Err(failure("preview boundary reached max_vertices"));
        }
        let index = self.mesh.positions.len();
        self.mesh.positions.push(point);
        Ok(index)
    }
    fn vertex(
        &mut self,
        model: &Model,
        vertex: VertexId,
        options: &PreviewOptions,
        dec: &Decider,
    ) -> Result<usize, Refusal> {
        if let Some(&index) = self.vertices.get(&vertex) {
            return Ok(index);
        }
        let p = model.geom().position(model.vertex(vertex).point);
        let point = chart::position(&[
            dec.enclose_at(&p.x, chart::BITS)?,
            dec.enclose_at(&p.y, chart::BITS)?,
            dec.enclose_at(&p.z, chart::BITS)?,
        ])?;
        let index = self.push(point, options)?;
        self.vertices.insert(vertex, index);
        Ok(index)
    }
}

fn sample(piece: &TracePiece, t: f64, dec: &Decider) -> Result<[f64; 3], Refusal> {
    if let Some(jets) = piece.taylor_interval(Interval::point(t).unwrap(), 0, dec)? {
        let values = jets[0].map(|x| x.mid());
        if values.iter().all(|x| x.is_finite()) {
            return Ok(values);
        }
    }
    chart::position(&piece.evaluate(&Ball::from_rational(&chart::q(t), chart::BITS), dec)?)
}

struct Sampler<'a> {
    piece: &'a TracePiece,
    charts: &'a [&'a Chart],
    options: &'a PreviewOptions,
    dec: &'a Decider,
    boundaries: &'a mut Boundaries,
    indices: &'a mut Vec<usize>,
}

impl Sampler<'_> {
    fn segment(&mut self, interval: [f64; 2], ends: [usize; 2], depth: u32) -> Result<(), Refusal> {
        self.dec.check_deadline()?;
        let [lo, hi] = interval;
        let [first, last] = ends;
        let a = self.boundaries.mesh.positions[first];
        let b = self.boundaries.mesh.positions[last];
        let mid = (lo + hi) * 0.5;
        let [p, q, r] = [lo + (hi - lo) * 0.25, mid, lo + (hi - lo) * 0.75]
            .map(|t| sample(self.piece, t, self.dec));
        let points = [p?, q?, r?];
        let mut small = points.iter().enumerate().all(|(i, p)| {
            let t = (i + 1) as f64 * 0.25;
            (0..3)
                .map(|k| (p[k] - (a[k] + (b[k] - a[k]) * t)).powi(2))
                .sum::<f64>()
                .sqrt()
                <= self.options.chord_tolerance.to_f64_nearest() * 0.25
        });
        for chart in self.charts {
            let uv_a = chart.project(a);
            let uv_b = align(chart.project(b), uv_a, chart.periodic);
            let uv_m = align(chart.project(points[1]), uv_a, chart.periodic);
            let limit = self
                .options
                .angular_deviation
                .as_ref()
                .map_or(0.125, |angle| {
                    (angle.to_f64_nearest() / std::f64::consts::TAU * 0.5).min(0.125)
                });
            for i in 0..2 {
                let angular = chart.periodic[i]
                    || (i == 1 && matches!(chart.surface, catallax_shapes::Surface::Sphere(_)));
                if angular {
                    small &= (uv_b[i] - uv_a[i]).abs() <= limit;
                    small &= (uv_m[i] - (uv_a[i] + uv_b[i]) * 0.5).abs() <= limit * 0.25;
                }
            }
        }
        if small {
            self.indices.push(last);
            return Ok(());
        }
        if depth >= 32 {
            return Err(failure(
                "preview edge sampling did not resolve visual targets",
            ));
        }
        let middle = self.boundaries.push(points[1], self.options)?;
        self.segment([lo, mid], [first, middle], depth + 1)?;
        self.segment([mid, hi], [middle, last], depth + 1)
    }
}

pub(crate) fn build(
    model: &Model,
    body: BodyId,
    charts: &BTreeMap<FaceId, Chart>,
    options: &PreviewOptions,
    dec: &mut Decider,
) -> Result<Boundaries, Refusal> {
    let mut result = Boundaries {
        mesh: PreviewMesh::default(),
        edges: BTreeMap::new(),
        vertices: BTreeMap::new(),
    };
    for vertex in model.vertices_of_body(body) {
        result.vertex(model, vertex, options, dec)?;
    }
    let mut adjacent = BTreeMap::<EdgeId, Vec<&Chart>>::new();
    for (&face, chart) in charts {
        for &boundary in &model.face(face).loops {
            if let catallax_topology::Loop::Coedges(coedges) = model.loop_(boundary) {
                for &coedge in coedges {
                    adjacent
                        .entry(model.coedge(coedge).edge)
                        .or_default()
                        .push(chart);
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
        let pieces = traced.value.pieces();
        let piece = pieces
            .first()
            .ok_or_else(|| failure("preview edge has no finite trace"))?;
        let first = if let Some(vertex) = from {
            result.vertices[&vertex]
        } else {
            result.push(sample(piece, 0.0, dec)?, options)?
        };
        let end = to.map_or(first, |vertex| result.vertices[&vertex]);
        let mut indices = vec![first];
        for (i, piece) in pieces.iter().enumerate() {
            let last = if i + 1 == pieces.len() {
                end
            } else {
                result.push(sample(piece, 1.0, dec)?, options)?
            };
            let start = *indices.last().unwrap();
            Sampler {
                piece,
                charts: &adjacent[&edge_id],
                options,
                dec,
                boundaries: &mut result,
                indices: &mut indices,
            }
            .segment([0.0, 1.0], [start, last], 0)?;
        }
        if edge.sense == Sense::Reversed {
            indices.reverse();
        }
        result.edges.insert(edge_id, result.mesh.edges.len());
        result.mesh.edges.push(PreviewEdge {
            edge: edge_id,
            vertices: indices,
        });
    }
    Ok(result)
}
