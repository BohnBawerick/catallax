//! Copies topology and remaps all defining geometry through an exact isometry.

use crate::*;
use catallax_numbers::{Fact, FactKind, FactOrigin};
use catallax_recipe::{EdgeName, Origin};
use catallax_shapes::{
    Curve, CurveId, GeomKey, Isometry, Meet, MeetOf, PointDef, PointId, Position, SurfaceId,
};
use std::collections::BTreeMap;

struct Copy<'a, F> {
    model: &'a mut Model,
    map: &'a Isometry,
    origin: F,
    surfaces: BTreeMap<SurfaceId, SurfaceId>,
    curves: BTreeMap<CurveId, CurveId>,
    points: BTreeMap<PointId, PointId>,
    vertices: BTreeMap<VertexId, VertexId>,
    edges: BTreeMap<EdgeId, EdgeId>,
}

fn reverse(sense: Sense) -> Sense {
    match sense {
        Sense::Same => Sense::Reversed,
        Sense::Reversed => Sense::Same,
    }
}

impl<F: FnMut(&Origin) -> Origin> Copy<'_, F> {
    fn incidence(&mut self, a: GeomKey, b: GeomKey, mapped_a: GeomKey, mapped_b: GeomKey) {
        if self.model.geom().facts().holds(FactKind::Incident, a, b) {
            self.model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: mapped_a,
                b: mapped_b,
                at: None,
                origin: FactOrigin::Construction("topology.isometry"),
            });
        }
    }
    fn surface(&mut self, id: SurfaceId) -> SurfaceId {
        if let Some(mapped) = self.surfaces.get(&id) {
            return *mapped;
        }
        let surface = self.map.surface(self.model.geom().surface(id));
        let mapped = self.model.geom_mut().add_surface(surface);
        self.surfaces.insert(id, mapped);
        mapped
    }
    fn point(&mut self, id: PointId) -> PointId {
        if let Some(mapped) = self.points.get(&id) {
            return *mapped;
        }
        let definition = self.model.geom().point(id).clone();
        let point = self.map.apply_point(&self.model.geom().position(id));
        let definition = match definition {
            PointDef::Explicit(_) => PointDef::Explicit(point),
            PointDef::Meet(meet) => {
                let of = match meet.of {
                    MeetOf::Surfaces(ids) => {
                        MeetOf::Surfaces(ids.into_iter().map(|s| self.surface(s)).collect())
                    }
                    MeetOf::CurveSurface(c, s) => {
                        MeetOf::CurveSurface(self.curve(c), self.surface(s))
                    }
                    MeetOf::Curves(a, b) => MeetOf::Curves(self.curve(a), self.curve(b)),
                };
                PointDef::Meet(Meet {
                    of,
                    position: Position::Closed(point),
                })
            }
        };
        let mapped = self.model.geom_mut().add_point(definition);
        self.points.insert(id, mapped);
        mapped
    }
    fn curve(&mut self, id: CurveId) -> CurveId {
        if let Some(mapped) = self.curves.get(&id) {
            return *mapped;
        }
        let curve = self.model.geom().curve(id).clone();
        let curve = if let Curve::Intersection(mut c) = curve {
            c.surfaces = c.surfaces.map(|s| self.surface(s));
            if self.map.reversing() {
                c.surfaces.swap(0, 1);
            }
            c.witness = self.point(c.witness);
            Curve::Intersection(c)
        } else {
            self.map.curve(&curve).expect("analytic curve")
        };
        let mapped = self.model.geom_mut().add_curve(curve);
        self.curves.insert(id, mapped);
        mapped
    }
    fn vertex(&mut self, id: VertexId) -> VertexId {
        if let Some(mapped) = self.vertices.get(&id) {
            return *mapped;
        }
        let point = self.point(self.model.vertex(id).point);
        let mapped = self.model.add_vertex(Vertex { point });
        self.vertices.insert(id, mapped);
        mapped
    }
    fn edge(&mut self, id: EdgeId) -> EdgeId {
        if let Some(mapped) = self.edges.get(&id) {
            return *mapped;
        }
        let edge = self.model.edge(id).clone();
        let curve = self.curve(edge.curve);
        let start = edge.start.map(|v| self.vertex(v));
        let end = edge.end.map(|v| self.vertex(v));
        for (source, mapped) in edge.start.zip(start).into_iter().chain(edge.end.zip(end)) {
            self.incidence(
                GeomKey::Point(self.model.vertex(source).point),
                GeomKey::Curve(edge.curve),
                GeomKey::Point(self.model.vertex(mapped).point),
                GeomKey::Curve(curve),
            );
        }
        let mapped = self.model.add_edge(Edge {
            curve,
            sense: edge.sense,
            start,
            end,
        });
        let constructions = self
            .model
            .edge_constructions(id)
            .iter()
            .filter_map(|name| {
                let [a, b] = name.between();
                EdgeName::new((self.origin)(a), (self.origin)(b)).ok()
            })
            .collect();
        self.model.record_edge_constructions(mapped, constructions);
        self.edges.insert(id, mapped);
        mapped
    }
    fn loop_(&mut self, id: LoopId) -> LoopId {
        let lp = self.model.loop_(id).clone();
        let lp = match lp {
            Loop::Vertex(v) => Loop::Vertex(self.vertex(v)),
            Loop::Coedges(mut ids) => {
                if self.map.reversing() {
                    ids.reverse();
                }
                Loop::Coedges(
                    ids.into_iter()
                        .map(|id| {
                            let coedge = self.model.coedge(id).clone();
                            let edge = self.edge(coedge.edge);
                            let sense = if self.map.reversing() {
                                reverse(coedge.sense)
                            } else {
                                coedge.sense
                            };
                            self.model.add_coedge(Coedge { edge, sense })
                        })
                        .collect(),
                )
            }
        };
        self.model.add_loop(lp)
    }
    fn face(&mut self, id: FaceId) -> FaceId {
        let face = self.model.face(id).clone();
        let origins = self
            .model
            .face_origins(id)
            .iter()
            .map(&mut self.origin)
            .collect();
        let surface = self.surface(face.surface);
        let loops = face.loops.iter().map(|l| self.loop_(*l)).collect();
        for lp in &face.loops {
            match self.model.loop_(*lp).clone() {
                Loop::Vertex(v) => {
                    let point = self.model.vertex(v).point;
                    let mapped = self.point(point);
                    self.incidence(
                        GeomKey::Point(point),
                        GeomKey::Surface(face.surface),
                        GeomKey::Point(mapped),
                        GeomKey::Surface(surface),
                    );
                }
                Loop::Coedges(coedges) => {
                    for coedge in coedges {
                        let curve = self.model.edge(self.model.coedge(coedge).edge).curve;
                        let mapped = self.curve(curve);
                        self.incidence(
                            GeomKey::Curve(curve),
                            GeomKey::Surface(face.surface),
                            GeomKey::Curve(mapped),
                            GeomKey::Surface(surface),
                        );
                    }
                }
            }
        }
        self.model.add_named_face(
            Face {
                surface,
                sense: face.sense,
                loops,
            },
            origins,
        )
    }
}

/// Copies a body by an exact isometry and maps every face origin.
/// Source entities remain immutable; shared edges and defining objects stay shared
/// within the copy. Reflection reverses loop traversal but preserves outward normals.
pub fn transformed_body(
    model: &mut Model,
    body: BodyId,
    map: &Isometry,
    origin: impl FnMut(&Origin) -> Origin,
) -> BodyId {
    let source = model.body(body).clone();
    let mut copy = Copy {
        model,
        map,
        origin,
        surfaces: BTreeMap::new(),
        curves: BTreeMap::new(),
        points: BTreeMap::new(),
        vertices: BTreeMap::new(),
        edges: BTreeMap::new(),
    };
    let lumps = source
        .lumps
        .into_iter()
        .map(|id| {
            let shells = copy
                .model
                .lump(id)
                .shells
                .clone()
                .into_iter()
                .map(|id| {
                    let faces = copy
                        .model
                        .shell(id)
                        .faces
                        .clone()
                        .into_iter()
                        .map(|f| copy.face(f))
                        .collect();
                    copy.model.add_shell(Shell { faces })
                })
                .collect();
            copy.model.add_lump(Lump { shells })
        })
        .collect();
    copy.model.add_body(Body {
        kind: source.kind,
        lumps,
        marks: source.marks,
    })
}
