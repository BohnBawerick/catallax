//! Topology construction from recognised faces.
//!
//! Faces are built exactly as the file states them. Every edge must lie on both
//! its faces and every vertex on its edges, checked in exact arithmetic; a file
//! that misses by any amount is refused here and left to the gap measurement.
//! Vertices join by exact coordinate equality and straight edges by their two
//! joined vertices. Nothing joins by proximity.

use std::collections::BTreeMap;

use catallax_numbers::{Fact, FactKind, FactOrigin, MarkLog, Rational};
use catallax_recipe::{Origin, Role, StepId};
use catallax_shapes::{
    Circle, Cone, Curve as GeomCurve, Cylinder, Ellipse, GeomKey, Line, Plane, Point3, PointDef,
    RPoint3, RVector3, Sphere, Surface as GeomSurface, Torus, Vector3,
};
use catallax_topology::{
    Body, BodyId, BodyKind, Coedge, Edge, EdgeId, Face, Loop, Lump, Model, Sense, Shell, Vertex,
    VertexId,
};

use crate::exact::{cross, dot, is_zero, sub, unit};
use crate::{
    Curve, Face as StepFace, ImportRefusal, OrientedEdge, Placement, Solid, StepFile, Surface, Vec3,
};

fn not_solid(detail: impl Into<String>) -> ImportRefusal {
    ImportRefusal::NotSolid {
        detail: detail.into(),
    }
}

fn unsupported(entity: &str, type_name: &str) -> ImportRefusal {
    ImportRefusal::Unsupported {
        entity: entity.to_string(),
        type_name: type_name.to_string(),
    }
}

fn point(p: &Vec3) -> Point3 {
    RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()).to_real()
}

fn vector(v: &Vec3) -> Vector3 {
    RVector3::new(v[0].clone(), v[1].clone(), v[2].clone()).to_real()
}

/// Builds every solid of a file into `model`, one body each, in file order.
/// Faces carry the origin `{step, imported: {key: "#n"}}` when `step` is given.
pub fn build_file(
    model: &mut Model,
    file: &StepFile,
    step: Option<&StepId>,
) -> Result<Vec<BodyId>, ImportRefusal> {
    file.solids
        .iter()
        .map(|solid| build_solid(model, solid, step))
        .collect()
}

/// Builds one solid into `model` and checks its topology.
pub fn build_solid(
    model: &mut Model,
    solid: &Solid,
    step: Option<&StepId>,
) -> Result<BodyId, ImportRefusal> {
    crate::intersect::check_solid(solid)?;
    let mut b = Builder::new(model, false);
    let mut faces = Vec::new();
    for face in &solid.faces {
        let origins = step
            .map(|step| {
                vec![Origin {
                    step: step.clone(),
                    role: Role::Imported {
                        key: face.name.clone(),
                    },
                }]
            })
            .unwrap_or_default();
        faces.push(b.face(face, origins)?);
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    if let Some(fault) = model.check(body).first() {
        return Err(not_solid(format!(
            "{} in {}: {}",
            fault.code, solid.name, fault.detail
        )));
    }
    Ok(body)
}

type Joined = (EdgeId, [Option<VertexId>; 2]);

struct Builder<'a> {
    model: &'a mut Model,
    /// Scratch models carry every recognised surface and edge for geometric
    /// queries. They skip the incidence checks and the topology rules.
    scratch: bool,
    vertices: BTreeMap<Vec3, VertexId>,
    edges: BTreeMap<String, Joined>,
    /// Straight edges by their joined vertex pair, lower id first.
    segments: BTreeMap<(VertexId, VertexId), Joined>,
}

/// One use of an edge in a loop, before seams are removed.
struct Use {
    edge: String,
    forward: bool,
}

impl<'a> Builder<'a> {
    fn new(model: &'a mut Model, scratch: bool) -> Self {
        Self {
            model,
            scratch,
            vertices: BTreeMap::new(),
            edges: BTreeMap::new(),
            segments: BTreeMap::new(),
        }
    }
}

/// Builds the faces of a solid into a scratch model for geometric queries only.
/// Any recognised surface and edge is carried as written, with no incidence or
/// topology checks and no body.
pub(crate) fn scratch_faces(
    model: &mut Model,
    solid: &Solid,
) -> Result<Vec<catallax_topology::FaceId>, ImportRefusal> {
    let mut b = Builder::new(model, true);
    solid.faces.iter().map(|f| b.face(f, Vec::new())).collect()
}

impl Builder<'_> {
    fn vertex(&mut self, p: &Vec3) -> VertexId {
        if let Some(v) = self.vertices.get(p) {
            return *v;
        }
        let id = self
            .model
            .geom_mut()
            .add_point(PointDef::Explicit(point(p)));
        let v = self.model.add_vertex(Vertex { point: id });
        self.vertices.insert(p.clone(), v);
        v
    }

    fn fact(&mut self, a: GeomKey, b: GeomKey) {
        self.model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a,
            b,
            at: None,
            origin: FactOrigin::Construction("import.step"),
        });
    }

    fn surface(&mut self, face: &StepFace) -> Result<GeomSurface, ImportRefusal> {
        Ok(match &face.surface {
            Surface::Plane(p) => GeomSurface::Plane(Plane {
                origin: point(&p.origin),
                normal: vector(&p.axis),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Cylinder(p, r) => GeomSurface::Cylinder(Cylinder {
                origin: point(&p.origin),
                axis: vector(&p.axis),
                radius_sq: (r * r).into(),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Cone45(p, r) => {
                let axis = unit(&p.axis).ok_or_else(|| {
                    unsupported(&face.name, "CONICAL_SURFACE (axis length is not rational)")
                })?;
                let apex: Vec3 = std::array::from_fn(|i| &p.origin[i] - &(&axis[i] * r));
                GeomSurface::Cone(Cone {
                    apex: point(&apex),
                    axis: vector(&axis),
                    slope: Rational::one().into(),
                    ref_dir: vector(&p.ref_dir),
                })
            }
            Surface::Sphere(p, r) if self.scratch => GeomSurface::Sphere(Sphere {
                center: point(&p.origin),
                radius: r.clone().into(),
                axis: vector(&p.axis),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Torus(p, major, minor) if self.scratch => GeomSurface::Torus(Torus {
                center: point(&p.origin),
                axis: vector(&p.axis),
                major: major.clone().into(),
                minor: minor.clone().into(),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Sphere(..) | Surface::Torus(..) => {
                return Err(unsupported(
                    &face.name,
                    "ADVANCED_FACE (spherical or toroidal face)",
                ));
            }
        })
    }

    fn edge(
        &mut self,
        e: &OrientedEdge,
    ) -> Result<(EdgeId, [Option<VertexId>; 2], bool), ImportRefusal> {
        // Returns the edge, its ends, and whether this entity runs against it.
        let (start, end) = (self.vertex(&e.start.point), self.vertex(&e.end.point));
        let ring = start == end;
        let sense = if e.curve_same_sense {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let line = matches!(e.curve, Curve::Line(..));
        if ring && line {
            return Err(not_solid(format!(
                "edge {} is a closed straight line",
                e.edge
            )));
        }
        let closed_form = matches!(e.curve, Curve::Circle(..))
            || (self.scratch && matches!(e.curve, Curve::Ellipse(..)));
        if ring && !closed_form {
            return Err(unsupported(
                &e.edge,
                "EDGE_CURVE (closed edge on a non-circle)",
            ));
        }
        let key = (start.min(end), start.max(end));
        if line {
            if let Some((id, ends)) = self.segments.get(&key) {
                let against = ends[0] != Some(start);
                return Ok((*id, *ends, against));
            }
        } else if let Some((id, ends)) = self.edges.get(&e.edge) {
            return Ok((*id, *ends, false));
        }
        let curve = match &e.curve {
            Curve::Line(o, d) => GeomCurve::Line(Line {
                origin: point(o),
                dir: vector(d),
            }),
            Curve::Circle(p, r) => GeomCurve::Circle(Circle {
                center: point(&p.origin),
                normal: vector(&p.axis),
                radius_sq: (r * r).into(),
                ref_dir: vector(&p.ref_dir),
            }),
            Curve::Ellipse(p, a, b) if self.scratch => {
                // The longer semi-axis is the major one; a circle is not an ellipse here.
                if a == b {
                    GeomCurve::Circle(Circle {
                        center: point(&p.origin),
                        normal: vector(&p.axis),
                        radius_sq: (a * a).into(),
                        ref_dir: vector(&p.ref_dir),
                    })
                } else {
                    let (major_dir, major, minor) = if a > b {
                        (p.ref_dir.clone(), a, b)
                    } else {
                        (cross(&p.axis, &p.ref_dir), b, a)
                    };
                    GeomCurve::Ellipse(Ellipse {
                        center: point(&p.origin),
                        normal: vector(&p.axis),
                        major_dir: vector(&major_dir),
                        major_sq: (major * major).into(),
                        minor_sq: (minor * minor).into(),
                    })
                }
            }
            Curve::Ellipse(..) => return Err(unsupported(&e.edge, "ELLIPSE")),
        };
        let curve = self.model.geom_mut().add_curve(curve);
        let ends = if ring {
            [None, None]
        } else {
            [Some(start), Some(end)]
        };
        let id = self.model.add_edge(Edge {
            curve,
            sense,
            start: ends[0],
            end: ends[1],
        });
        if line {
            self.segments.insert(key, (id, ends));
        } else {
            self.edges.insert(e.edge.clone(), (id, ends));
        }
        Ok((id, ends, false))
    }

    fn face(
        &mut self,
        face: &StepFace,
        origins: Vec<Origin>,
    ) -> Result<catallax_topology::FaceId, ImportRefusal> {
        let geom = self.surface(face)?;
        let surface = self.model.geom_mut().add_surface(geom);
        let mut loops = Vec::new();
        for bound in &face.loops {
            for chain in Self::chains(&face.name, bound.orientation, &bound.edges)? {
                let mut coedges = Vec::new();
                for (e, forward) in chain {
                    let (edge, _, against) = self.edge(e)?;
                    let sense = if forward != against {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    };
                    let data = self.model.edge(edge).clone();
                    if !self.scratch {
                        check_on_surface(face, e)?;
                        check_on_curve(e)?;
                        self.fact(GeomKey::Curve(data.curve), GeomKey::Surface(surface));
                        for v in [data.start, data.end].into_iter().flatten() {
                            let p = self.model.vertex(v).point;
                            self.fact(GeomKey::Point(p), GeomKey::Curve(data.curve));
                        }
                    }
                    coedges.push(self.model.add_coedge(Coedge { edge, sense }));
                }
                loops.push(self.model.add_loop(Loop::Coedges(coedges)));
            }
        }
        if loops.is_empty() && !self.scratch {
            return Err(unsupported(&face.name, "ADVANCED_FACE (no bounds)"));
        }
        Ok(self.model.add_named_face(
            Face {
                surface,
                sense: if face.same_sense {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                loops,
            },
            origins,
        ))
    }

    /// Splits a bound into closed chains after removing seams, edges used once
    /// each way in the same bound. The orientation flag reverses the bound.
    fn chains<'e>(
        face: &str,
        orientation: bool,
        edges: &'e [OrientedEdge],
    ) -> Result<Vec<Vec<(&'e OrientedEdge, bool)>>, ImportRefusal> {
        let uses: Vec<Use> = edges
            .iter()
            .map(|e| Use {
                edge: e.edge.clone(),
                forward: e.forward,
            })
            .collect();
        let mut kept: Vec<(&OrientedEdge, bool)> = Vec::new();
        for (i, e) in edges.iter().enumerate() {
            let seam = uses
                .iter()
                .enumerate()
                .any(|(j, u)| j != i && u.edge == uses[i].edge && u.forward != uses[i].forward);
            if !seam {
                kept.push((e, e.forward));
            }
        }
        if !orientation {
            kept.reverse();
            for k in &mut kept {
                k.1 = !k.1;
            }
        }
        let ends = |e: &OrientedEdge, forward: bool| {
            if forward {
                (e.start.point.clone(), e.end.point.clone())
            } else {
                (e.end.point.clone(), e.start.point.clone())
            }
        };
        let mut chains: Vec<Vec<(&OrientedEdge, bool)>> = Vec::new();
        let mut current: Vec<(&OrientedEdge, bool)> = Vec::new();
        for (e, forward) in kept {
            if let Some((last, lf)) = current.last()
                && ends(last, *lf).1 != ends(e, forward).0
            {
                return Err(not_solid(format!(
                    "boundary of face {face} is not a connected cycle at {}",
                    e.edge
                )));
            }
            current.push((e, forward));
            let (first, ff) = current[0];
            if ends(first, ff).0 == ends(e, forward).1 {
                chains.push(std::mem::take(&mut current));
            }
        }
        if !current.is_empty() {
            return Err(not_solid(format!("boundary of face {face} does not close")));
        }
        Ok(chains)
    }
}

fn on_plane(origin: &Vec3, normal: &Vec3, p: &Vec3) -> bool {
    dot(&sub(p, origin), normal).is_zero()
}

fn on_axis(placement: &Placement, p: &Vec3) -> bool {
    is_zero(&cross(&sub(p, &placement.origin), &placement.axis))
}

/// Exact test that the edge's curve lies on the face's surface.
fn check_on_surface(face: &StepFace, e: &OrientedEdge) -> Result<(), ImportRefusal> {
    let ok = match (&face.surface, &e.curve) {
        (Surface::Plane(p), Curve::Line(o, d)) => {
            on_plane(&p.origin, &p.axis, o) && dot(d, &p.axis).is_zero()
        }
        (Surface::Plane(p), Curve::Circle(c, _)) => {
            on_plane(&p.origin, &p.axis, &c.origin) && is_zero(&cross(&c.axis, &p.axis))
        }
        (Surface::Cylinder(p, r), Curve::Line(o, d)) => {
            let rel = cross(&sub(o, &p.origin), &p.axis);
            is_zero(&cross(d, &p.axis)) && dot(&rel, &rel) == r * r * dot(&p.axis, &p.axis)
        }
        (Surface::Cylinder(p, r), Curve::Circle(c, cr)) => {
            cr == r && is_zero(&cross(&c.axis, &p.axis)) && on_axis(p, &c.origin)
        }
        (Surface::Cone45(p, r), Curve::Circle(c, cr)) => {
            // Radius equals the signed axial distance from the apex.
            let Some(axis) = unit(&p.axis) else {
                return Err(unsupported(&face.name, "CONICAL_SURFACE"));
            };
            let t = dot(&sub(&c.origin, &p.origin), &axis);
            is_zero(&cross(&c.axis, &p.axis)) && on_axis(p, &c.origin) && *cr == r + &t
        }
        _ => {
            return Err(unsupported(
                &e.edge,
                "EDGE_CURVE (curve and surface pair outside the supported set)",
            ));
        }
    };
    if ok {
        Ok(())
    } else {
        Err(not_solid(format!(
            "edge {} does not lie exactly on face {}",
            e.edge, face.name
        )))
    }
}

/// Exact test that both ends of the edge lie on its curve.
fn check_on_curve(e: &OrientedEdge) -> Result<(), ImportRefusal> {
    let on = |p: &Vec3| match &e.curve {
        Curve::Line(o, d) => is_zero(&cross(&sub(p, o), d)),
        Curve::Circle(c, r) => {
            let rel = sub(p, &c.origin);
            dot(&rel, &c.axis).is_zero() && dot(&rel, &rel) == r * r
        }
        Curve::Ellipse(..) => false,
    };
    if on(&e.start.point) && on(&e.end.point) {
        Ok(())
    } else {
        Err(not_solid(format!(
            "a vertex of edge {} does not lie exactly on its curve",
            e.edge
        )))
    }
}
