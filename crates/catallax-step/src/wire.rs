//! Surface-CCW chart walks over exact shared physical edges.

use crate::{
    chart::Meridian,
    cut::{self, SplitEdges},
    number,
};
use catallax_numbers::{Fact, FactKind, FactOrigin, Refusal};
use catallax_shapes::{Circle, Curve, GeomKey, Line, Surface, SurfaceId};
use catallax_topology::{
    Coedge, CoedgeId, Edge, EdgeId, Face, FaceId, Loop, Model, Sense, VertexId,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Side {
    Low,
    Interior,
    High,
}

pub(crate) struct Use {
    pub coedge: Coedge,
    /// Endpoints in the underlying curve's natural direction, not wire order.
    pub sides: [Side; 2],
}

pub(crate) struct FaceCut {
    pub surface: SurfaceId,
    pub loops: Vec<Vec<Use>>,
    pub seams: BTreeSet<EdgeId>,
}

pub(crate) struct ChartFace {
    pub touching: bool,
    pub sides: BTreeMap<CoedgeId, [Side; 2]>,
    pub seams: BTreeSet<EdgeId>,
}

impl FaceCut {
    pub fn split_touching(mut self, model: &Model) -> Result<Self, Refusal> {
        self.loops = split_touching_walks(model, self.loops, &self.seams)?;
        Ok(self)
    }

    pub fn store(self, model: &mut Model, source: FaceId) -> (FaceId, ChartFace) {
        let sense = model.face(source).sense;
        let mut loops = Vec::with_capacity(self.loops.len());
        let mut sides = BTreeMap::new();
        for mut wire in self.loops {
            if sense == Sense::Reversed {
                wire.reverse();
                for edge_use in &mut wire {
                    edge_use.coedge.sense = if edge_use.coedge.sense == Sense::Same {
                        Sense::Reversed
                    } else {
                        Sense::Same
                    };
                }
            }
            let mut coedges = Vec::with_capacity(wire.len());
            for edge_use in wire {
                let id = model.add_coedge(edge_use.coedge);
                sides.insert(id, edge_use.sides);
                coedges.push(id);
            }
            loops.push(model.add_loop(Loop::Coedges(coedges)));
        }
        let face = model.add_face(Face {
            surface: self.surface,
            sense,
            loops,
        });
        (
            face,
            ChartFace {
                touching: false,
                sides,
                seams: self.seams,
            },
        )
    }
}
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Node {
    vertex: VertexId,
    side: Side,
    corner: usize,
}

#[derive(Default)]
struct Walks {
    outgoing: BTreeMap<Node, (Node, Use)>,
}

impl Walks {
    fn add(&mut self, model: &Model, edge_use: Use, corners: [usize; 2]) -> Result<(), Refusal> {
        let edge = model.edge(edge_use.coedge.edge);
        let (Some(a), Some(b)) = (edge.start, edge.end) else {
            return Err(Refusal::new(
                "step.not-solid",
                "chart edge has no finite boundary",
            ));
        };
        let (a, b) = if edge.sense == Sense::Same {
            (a, b)
        } else {
            (b, a)
        };
        let from = Node {
            vertex: a,
            side: edge_use.sides[0],
            corner: 0,
        };
        let to = Node {
            vertex: b,
            side: edge_use.sides[1],
            corner: 0,
        };
        let (mut from, mut to) = if edge_use.coedge.sense == edge.sense {
            (from, to)
        } else {
            (to, from)
        };
        // Preserve source-loop incidence at touching trim vertices. Only a
        // meridian cut reconnects walks by physical vertex and chart side.
        if from.side == Side::Interior {
            from.corner = corners[0];
        }
        if to.side == Side::Interior {
            to.corner = corners[1];
        }
        if self.outgoing.insert(from, (to, edge_use)).is_some() {
            return Err(number::budget(
                "periodic cut has no unique outgoing chart edge",
            ));
        }
        Ok(())
    }

    fn finish(mut self, seams: &BTreeSet<EdgeId>) -> Result<Vec<Vec<Use>>, Refusal> {
        let mut loops = Vec::new();
        while let Some((&start, _)) = self.outgoing.first_key_value() {
            let mut at = start;
            let mut wire = Vec::new();
            loop {
                let (next, edge_use) = self.outgoing.remove(&at).ok_or_else(|| {
                    number::budget("periodic chart edges do not form closed cycles")
                })?;
                wire.push(edge_use);
                if next == start {
                    break;
                }
                at = next;
            }
            seam_first(&mut wire, seams)?;
            loops.push(wire);
        }
        Ok(loops)
    }
}

/// Raw STEP translation chooses seam representations at the first occurrence.
/// Keep the native-proven reversed occurrence first for every seam used twice
/// by one cycle. Distinct cycles will receive distinct surface identities.
fn seam_first(wire: &mut [Use], seams: &BTreeSet<EdgeId>) -> Result<(), Refusal> {
    if seams.is_empty() {
        return Ok(());
    }
    let mut occurrences = BTreeMap::<EdgeId, Vec<(usize, bool)>>::new();
    for (i, u) in wire.iter().enumerate() {
        if seams.contains(&u.coedge.edge) {
            occurrences
                .entry(u.coedge.edge)
                .or_default()
                .push((i, u.coedge.sense == Sense::Same));
        }
    }
    if occurrences.is_empty() {
        return Ok(());
    }
    let mut delta = vec![0i64; wire.len() + 1];
    for uses in occurrences.values() {
        if uses.len() != 2 {
            continue;
        }
        let (forward, reversed) = match (uses[0].1, uses[1].1) {
            (true, false) => (uses[0].0, uses[1].0),
            (false, true) => (uses[1].0, uses[0].0),
            _ => {
                return Err(number::budget(
                    "seam occurrences do not have opposite directions",
                ));
            }
        };
        if reversed < forward {
            delta[reversed + 1] += 1;
            delta[forward + 1] -= 1;
        } else {
            delta[reversed + 1] += 1;
            delta[wire.len()] -= 1;
            delta[0] += 1;
            delta[forward + 1] -= 1;
        }
    }
    let mut forbidden = 0;
    for (index, change) in delta.into_iter().take(wire.len()).enumerate() {
        forbidden += change;
        if forbidden == 0 {
            wire.rotate_left(index);
            return Ok(());
        }
    }
    Err(number::budget("no native-compatible seam-first chart walk"))
}

pub(crate) fn cut_face(
    model: &mut Model,
    face: FaceId,
    meridian: Meridian,
    splits: &SplitEdges,
) -> Result<FaceCut, Refusal> {
    let definition = model.face(face).clone();
    let mut surface = model.geom().surface(definition.surface).clone();
    match &mut surface {
        Surface::Cylinder(s) => s.ref_dir = meridian.reference.clone(),
        Surface::Sphere(s) => s.ref_dir = meridian.reference.clone(),
        Surface::Cone(s) => s.ref_dir = meridian.reference.clone(),
        Surface::Torus(s) => s.ref_dir = meridian.reference.clone(),
        Surface::TorusSheet(s) => s.ref_dir = meridian.reference.clone(),
        _ => unreachable!("chart walk on a non-periodic carrier"),
    }
    let surface = model.geom_mut().add_surface(surface);
    let mut cuts = BTreeMap::new();
    for crossing in &meridian.crossings {
        let vertex = splits.vertices[&(crossing.edge, crossing.point)];
        cuts.insert(vertex, crossing.increasing);
    }
    let mut walks = Walks::default();
    let mut corner = 0;
    for loop_id in definition.loops {
        let Loop::Coedges(coedges) = model.loop_(loop_id) else {
            return Err(Refusal::new(
                "step.not-solid",
                "regular periodic chart has a vertex-only boundary",
            ));
        };
        let mut expanded = splits.expand(model, coedges);
        if definition.sense == Sense::Reversed {
            cut::reverse(&mut expanded);
        }
        let length = expanded.len();
        for (index, coedge) in expanded.into_iter().enumerate() {
            let edge = model.edge(coedge.edge);
            let curve = edge.curve;
            let (Some(a), Some(b)) = (edge.start, edge.end) else {
                return Err(Refusal::new(
                    "step.not-solid",
                    "chart edge has no finite boundary",
                ));
            };
            let (a, b) = if edge.sense == Sense::Same {
                (a, b)
            } else {
                (b, a)
            };
            let curve_key = GeomKey::Curve(curve);
            if !model
                .geom()
                .facts()
                .holds(FactKind::Incident, curve_key, GeomKey::Surface(surface))
            {
                if !model.geom().facts().holds(
                    FactKind::Incident,
                    curve_key,
                    GeomKey::Surface(definition.surface),
                ) {
                    return Err(Refusal::new(
                        "step.not-solid",
                        "source trim has no carrier-incidence proof",
                    ));
                }
                // Only the reference direction changed: the implicit carrier
                // equation and the source curve's incidence are unchanged.
                model.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a: curve_key,
                    b: GeomKey::Surface(surface),
                    at: None,
                    origin: FactOrigin::Construction("step.chart-reparameterization"),
                });
            }
            let side = |vertex, at_start| {
                cuts.get(&vertex).map_or(Side::Interior, |increasing| {
                    if at_start == *increasing {
                        Side::Low
                    } else {
                        Side::High
                    }
                })
            };
            walks.add(
                model,
                Use {
                    coedge,
                    sides: [side(a, true), side(b, false)],
                },
                [corner + index, corner + (index + 1) % length],
            )?;
        }
        corner += length;
    }
    let mut seams = BTreeSet::new();
    if !meridian.crossings.is_empty() {
        if !meridian.crossings.len().is_multiple_of(2) {
            return Err(number::budget(
                "periodic chart does not have paired meridian crossings",
            ));
        }
        let reference = number::unit(&meridian.reference);
        let curve = match model.geom().surface(surface) {
            Surface::Cylinder(s) => Curve::Line(Line {
                origin: s
                    .origin
                    .plus(&reference.scaled(&number::sqrt_real(&s.radius_sq))),
                // Reuse the carrier, never normalize an algebraic endpoint difference.
                dir: s.axis.clone(),
            }),
            Surface::Sphere(s) => Curve::Circle(Circle {
                center: s.center.clone(),
                normal: meridian.reference.cross(&s.axis),
                radius_sq: s.radius.square(),
                ref_dir: meridian.reference.clone(),
            }),
            Surface::Cone(s) => Curve::Line(Line {
                origin: s.apex.clone(),
                dir: number::unit(&s.axis).plus(&reference.scaled(&s.slope)),
            }),
            // The tube circle's angle is the torus v parameter: it starts on
            // the outer equator and increases toward the axis direction.
            Surface::Torus(s) => Curve::Circle(Circle {
                center: s.center.plus(&reference.scaled(&s.major)),
                normal: meridian.reference.cross(&s.axis),
                radius_sq: s.minor.square(),
                ref_dir: meridian.reference.clone(),
            }),
            Surface::TorusSheet(s) => Curve::Circle(Circle {
                center: s.center.plus(&reference.scaled(&s.major)),
                normal: meridian.reference.cross(&s.axis),
                radius_sq: s.minor.square(),
                ref_dir: meridian.reference.clone(),
            }),
            _ => unreachable!(),
        };
        let curve = model.geom_mut().add_curve(curve);
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Curve(curve),
            b: GeomKey::Surface(surface),
            at: None,
            origin: FactOrigin::Construction("step.chart-meridian"),
        });
        for pair in meridian.crossings.as_chunks::<2>().0 {
            if !pair[0].enters || pair[1].enters {
                return Err(number::budget(
                    "periodic chart crossing directions do not bound an interior interval",
                ));
            }
            let start = splits.vertices[&(pair[0].edge, pair[0].point)];
            let end = splits.vertices[&(pair[1].edge, pair[1].point)];
            for point in [model.vertex(start).point, model.vertex(end).point] {
                model.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a: GeomKey::Point(point),
                    b: GeomKey::Curve(curve),
                    at: None,
                    origin: FactOrigin::Construction("step.chart-meridian"),
                });
            }
            let edge = model.add_edge(Edge {
                curve,
                start: Some(start),
                end: Some(end),
                sense: Sense::Same,
            });
            seams.insert(edge);
            walks.add(
                model,
                Use {
                    coedge: Coedge {
                        edge,
                        sense: Sense::Same,
                    },
                    sides: [Side::High; 2],
                },
                [0; 2],
            )?;
            walks.add(
                model,
                Use {
                    coedge: Coedge {
                        edge,
                        sense: Sense::Reversed,
                    },
                    sides: [Side::Low; 2],
                },
                [0; 2],
            )?;
        }
    }
    let loops = walks.finish(&seams)?;
    Ok(FaceCut {
        surface,
        loops,
        seams,
    })
}

/// Separate lobes that meet only at a trim vertex before UV nesting.
fn split_touching_walks(
    model: &Model,
    mut pending: Vec<Vec<Use>>,
    seams: &BTreeSet<EdgeId>,
) -> Result<Vec<Vec<Use>>, Refusal> {
    let mut loops = Vec::new();
    pending.reverse();
    while let Some(mut wire) = pending.pop() {
        let mut starts = BTreeMap::new();
        let mut repeated = None;
        for (index, edge_use) in wire.iter().enumerate() {
            let edge = model.edge(edge_use.coedge.edge);
            let (vertex, side) = if edge_use.coedge.sense == Sense::Same {
                (
                    edge.start,
                    edge_use.sides[usize::from(edge.sense != Sense::Same)],
                )
            } else {
                (
                    edge.end,
                    edge_use.sides[usize::from(edge.sense == Sense::Same)],
                )
            };
            let vertex = vertex.expect("finite chart edge");
            if let Some(first) = starts.insert((vertex, side), index) {
                repeated = Some((first, index));
                break;
            }
        }
        if let Some((first, second)) = repeated {
            let tail = wire.split_off(second);
            let cycle = wire.split_off(first);
            wire.extend(tail);
            pending.push(wire);
            pending.push(cycle);
        } else {
            seam_first(&mut wire, seams)?;
            loops.push(wire);
        }
    }
    Ok(loops)
}
