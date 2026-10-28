//! Shared physical edge subdivision for export-only chart cuts.

use crate::{
    chart::Meridian,
    number,
    trim::{self, Geometry},
};
use catallax_intersect::Cyclic;
use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Outcome, Refusal, StrictOrder};
use catallax_shapes::{Curve, CurveId, GeomKey, PointId};
use catallax_topology::{
    Coedge, CoedgeId, Edge, EdgeId, FaceId, Loop, Model, Sense, Vertex, VertexId,
};
use std::{cmp::Ordering, collections::BTreeMap};

#[derive(Default)]
pub(crate) struct Cuts {
    points: BTreeMap<EdgeId, Vec<PointId>>,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Endpoint {
    Original(VertexId),
    Cut(EdgeId, PointId),
}

struct Segment {
    from: Endpoint,
    to: Endpoint,
}
struct Split {
    source: EdgeId,
    curve: CurveId,
    reversed: bool,
    segments: Vec<Segment>,
}
pub(crate) struct SplitPlan {
    splits: Vec<Split>,
}

pub(crate) struct SplitEdges {
    /// Walk each original edge in its original start-to-end direction.
    pub replacements: BTreeMap<EdgeId, Vec<Coedge>>,
    pub vertices: BTreeMap<(EdgeId, PointId), VertexId>,
    pub added_edges: usize,
    pub added_vertices: usize,
}

impl Cuts {
    /// Adjacent charts can encounter the same cut through distinct exact queries.
    /// Canonicalize only within an original topological edge; never sew bodies.
    pub fn add(
        &mut self,
        geometry: &Geometry,
        meridian: &mut Meridian,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        for crossing in &mut meridian.crossings {
            let points = self.points.entry(crossing.edge).or_default();
            let mut canonical = None;
            for &point in points.iter() {
                if geometry.same_point(point, crossing.point, dec)? {
                    canonical = Some(point);
                    break;
                }
            }
            if let Some(point) = canonical {
                crossing.point = point;
            } else {
                points.push(crossing.point);
            }
        }
        Ok(())
    }

    pub fn rings(
        &mut self,
        geometry: &mut Geometry,
        model: &Model,
        faces: impl Iterator<Item = FaceId>,
    ) {
        for face in faces {
            for loop_id in &model.face(face).loops {
                let Loop::Coedges(coedges) = model.loop_(*loop_id) else {
                    continue;
                };
                for coedge in coedges {
                    let edge_id = model.coedge(*coedge).edge;
                    let edge = model.edge(edge_id);
                    if edge.start.is_some()
                        || edge.end.is_some()
                        || self.points.contains_key(&edge_id)
                    {
                        continue;
                    }
                    let point = if let Curve::Intersection(curve) = geometry.store.curve(edge.curve)
                    {
                        curve.witness
                    } else {
                        let point = geometry
                            .store
                            .curve(edge.curve)
                            .parameter_origin(&geometry.store);
                        geometry.point(point)
                    };
                    geometry.store.record_fact(Fact {
                        kind: FactKind::Incident,
                        a: GeomKey::Point(point),
                        b: GeomKey::Curve(edge.curve),
                        at: None,
                        origin: FactOrigin::Construction("step.chart-ring-origin"),
                    });
                    self.points.insert(edge_id, vec![point]);
                }
            }
        }
    }

    pub fn prepare(
        &self,
        geometry: &mut Geometry,
        model: &Model,
        dec: &mut Decider,
    ) -> Result<SplitPlan, Refusal> {
        let mut splits = Vec::with_capacity(self.points.len());
        for (&source, points) in &self.points {
            let edge = model.edge(source);
            let reversed = edge.sense == Sense::Reversed;
            let (from, to) = if reversed {
                (edge.end, edge.start)
            } else {
                (edge.start, edge.end)
            };
            let closed = geometry.closed(edge.curve, dec)?;
            let anchor = from.map(|v| model.vertex(v).point).unwrap_or(points[0]);
            let mut points = points.clone();
            trim::try_sort(&mut points, |a, b| {
                if a == b {
                    return Ok(Ordering::Equal);
                }
                if closed {
                    if *a == anchor {
                        return Ok(Ordering::Less);
                    }
                    if *b == anchor {
                        return Ok(Ordering::Greater);
                    }
                    match catallax_intersect::cyclic_order(
                        &geometry.store,
                        edge.curve,
                        [anchor, *a, *b],
                        dec,
                    ) {
                        Outcome::Proven(Cyclic::Forward) => Ok(Ordering::Less),
                        Outcome::Proven(Cyclic::Backward) => Ok(Ordering::Greater),
                        Outcome::ProvenCoincident(_) => {
                            Err(number::budget("cut point order is not strict"))
                        }
                        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
                        Outcome::Refused(r) => Err(r),
                    }
                } else {
                    match catallax_intersect::order_on_curve(
                        &geometry.store,
                        edge.curve,
                        *a,
                        *b,
                        dec,
                    ) {
                        Outcome::Proven(StrictOrder::Less) => Ok(Ordering::Less),
                        Outcome::Proven(StrictOrder::Greater) => Ok(Ordering::Greater),
                        Outcome::ProvenCoincident(_) => {
                            Err(number::budget("cut point order is not strict"))
                        }
                        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
                        Outcome::Refused(r) => Err(r),
                    }
                }
            })?;
            let mut endpoints = Vec::with_capacity(points.len() + 2);
            match (from, to) {
                (Some(from), Some(to)) => {
                    endpoints.push(Endpoint::Original(from));
                    endpoints.extend(points.into_iter().map(|p| Endpoint::Cut(source, p)));
                    endpoints.push(Endpoint::Original(to));
                }
                (None, None) if closed => {
                    endpoints.extend(points.into_iter().map(|p| Endpoint::Cut(source, p)));
                    endpoints.push(endpoints[0]);
                }
                _ => {
                    return Err(Refusal::new(
                        "step.not-solid",
                        "cut edge has no finite boundary",
                    ));
                }
            }
            let segments = endpoints
                .windows(2)
                .map(|p| Segment {
                    from: p[0],
                    to: p[1],
                })
                .collect();
            splits.push(Split {
                source,
                curve: edge.curve,
                reversed,
                segments,
            });
        }
        Ok(SplitPlan { splits })
    }
}

impl SplitPlan {
    /// `model` already contains the exact scratch geometry used to prepare cuts.
    /// Its original topology and all original vertex identities remain untouched.
    pub fn apply(self, model: &mut Model) -> SplitEdges {
        let mut result = SplitEdges {
            replacements: BTreeMap::new(),
            vertices: BTreeMap::new(),
            added_edges: 0,
            added_vertices: 0,
        };
        for split in self.splits {
            result.added_edges += split.segments.len() - 1;
            let mut replacement = Vec::with_capacity(split.segments.len());
            for segment in split.segments {
                let mut vertex = |endpoint| match endpoint {
                    Endpoint::Original(v) => v,
                    Endpoint::Cut(edge, point) => {
                        *result.vertices.entry((edge, point)).or_insert_with(|| {
                            result.added_vertices += 1;
                            model.add_vertex(Vertex { point })
                        })
                    }
                };
                let start = vertex(segment.from);
                let end = vertex(segment.to);
                let edge = model.add_edge(Edge {
                    curve: split.curve,
                    start: Some(start),
                    end: Some(end),
                    sense: Sense::Same,
                });
                replacement.push(Coedge {
                    edge,
                    sense: Sense::Same,
                });
            }
            if split.reversed {
                reverse(&mut replacement);
            }
            result.replacements.insert(split.source, replacement);
        }
        result
    }
}

impl SplitEdges {
    pub fn expand(&self, model: &Model, coedges: &[CoedgeId]) -> Vec<Coedge> {
        let capacity = coedges
            .iter()
            .map(|id| {
                self.replacements
                    .get(&model.coedge(*id).edge)
                    .map_or(1, Vec::len)
            })
            .sum();
        let mut result = Vec::with_capacity(capacity);
        for id in coedges {
            let coedge = model.coedge(*id);
            if let Some(replacement) = self.replacements.get(&coedge.edge) {
                if coedge.sense == Sense::Same {
                    result.extend(replacement.iter().cloned());
                } else {
                    result.extend(replacement.iter().rev().map(|c| Coedge {
                        edge: c.edge,
                        sense: if c.sense == Sense::Same {
                            Sense::Reversed
                        } else {
                            Sense::Same
                        },
                    }));
                }
            } else {
                result.push(coedge.clone());
            }
        }
        result
    }
}

pub(crate) fn reverse(coedges: &mut [Coedge]) {
    coedges.reverse();
    for coedge in coedges {
        coedge.sense = if coedge.sense == Sense::Same {
            Sense::Reversed
        } else {
            Sense::Same
        };
    }
}
