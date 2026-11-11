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
    let has_spline = solid.faces.iter().any(|face| {
        matches!(face.surface, Surface::BSpline(_))
            || face
                .loops
                .iter()
                .flat_map(|bound| &bound.edges)
                .any(|edge| matches!(edge.curve, Curve::BSpline(_)))
    });
    if !has_spline {
        crate::intersect::check_solid(solid)?;
    }
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
    if has_spline {
        crate::intersect::check_solid_with_scratch(
            solid,
            Some(crate::curved::Scratch {
                model: model.clone(),
                faces: faces.clone(),
            }),
        )?;
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
    /// Analytic query geometry skips incidence checks. Spline geometry still
    /// needs its support, finite-edge and boundary certificates.
    scratch: bool,
    dec: catallax_numbers::Decider,
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
            dec: catallax_numbers::Decider::default(),
            vertices: BTreeMap::new(),
            edges: BTreeMap::new(),
            segments: BTreeMap::new(),
        }
    }
}

/// Builds the faces of a solid into a scratch model for geometric queries only.
/// Analytic geometry skips incidence checks. Spline geometry retains its
/// admission and incidence checks. No body is built.
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
            Surface::Cone(..) => {
                let cone = crate::cone::rim_cone(face)?.expect("conical face");
                let p = match &face.surface {
                    Surface::Cone(p, ..) => p,
                    _ => unreachable!(),
                };
                GeomSurface::Cone(Cone {
                    apex: point(&cone.apex),
                    axis: vector(&cone.axis),
                    slope: cone.slope.into(),
                    ref_dir: vector(&p.ref_dir),
                })
            }
            Surface::Sphere(p, r) => GeomSurface::Sphere(Sphere {
                center: point(&p.origin),
                radius: r.clone().into(),
                axis: vector(&p.axis),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::Torus(p, major, minor) => GeomSurface::Torus(Torus {
                center: point(&p.origin),
                axis: vector(&p.axis),
                major: major.clone().into(),
                minor: minor.clone().into(),
                ref_dir: vector(&p.ref_dir),
            }),
            Surface::BSpline(source) => {
                crate::spline_build::surface(&face.name, source, &mut self.dec)?
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
            Curve::BSpline(_) => crate::spline_build::curve(e, &mut self.dec)?,
        };
        let line = matches!(curve, GeomCurve::Line(_));
        if line && let Some((id, ends)) = self.segments.get(&key) {
            return Ok((*id, *ends, ends[0] != Some(start)));
        }
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
        // Spheres and tori are closed. A face on one whose bounds reduce to
        // seams, or to one vertex loop, is the complete surface: no loops.
        let closed = matches!(face.surface, Surface::Sphere(..) | Surface::Torus(..));
        let apex = cone_apex(face)?;
        let mut loops = Vec::new();
        let mut points = 0;
        for bound in &face.loops {
            if let Some(v) = &bound.vertex {
                if let Some(apex) = &apex {
                    if v.point != *apex {
                        return Err(not_solid(format!(
                            "vertex {} of a vertex loop is not the apex of face {}",
                            v.name, face.name
                        )));
                    }
                    let vertex = self.vertex(&v.point);
                    if !self.scratch {
                        let p = self.model.vertex(vertex).point;
                        self.fact(GeomKey::Point(p), GeomKey::Surface(surface));
                    }
                    loops.push(self.model.add_loop(Loop::Vertex(vertex)));
                    continue;
                }
                if !closed {
                    return Err(unsupported(
                        &face.name,
                        "VERTEX_LOOP (bound of a face on an open surface)",
                    ));
                }
                if !self.scratch && !on_closed_surface(&face.surface, &v.point) {
                    return Err(not_solid(format!(
                        "vertex {} of a vertex loop does not lie exactly on face {}",
                        v.name, face.name
                    )));
                }
                points += 1;
                continue;
            }
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
                    let spline_surface = matches!(face.surface, Surface::BSpline(_));
                    let spline_edge = matches!(e.curve, Curve::BSpline(_));
                    if !self.scratch || spline_surface || spline_edge {
                        if !spline_surface {
                            let whole_curve = if spline_edge {
                                let incidence = crate::spline_build::on_surface(
                                    self.model.geom().curve(data.curve),
                                    self.model.geom().surface(surface),
                                    e,
                                    &mut self.dec,
                                )
                                .map_err(|error| {
                                    crate::spline_build::refused(
                                        &e.edge,
                                        "EDGE_CURVE",
                                        "curve-surface incidence",
                                        error,
                                    )
                                })?;
                                if incidence == crate::spline_build::Incidence::OffSurface {
                                    return Err(not_solid(format!(
                                        "spline curve-surface incidence: edge {} does not lie exactly on face {}",
                                        e.edge, face.name
                                    )));
                                }
                                incidence == crate::spline_build::Incidence::WholeCurve
                            } else {
                                check_on_surface(face, e)?;
                                true
                            };
                            // Incident facts quantify over the full retained carrier.
                            if whole_curve {
                                self.fact(GeomKey::Curve(data.curve), GeomKey::Surface(surface));
                            }
                        }
                        if !spline_edge {
                            check_on_curve(e)?;
                        }
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
        // A pointed cone may close at its apex through a seam instead of a
        // vertex loop: one rim, and a seam ending at the apex.
        if let Some(apex) = &apex
            && loops.len() == 1
            && !face.loops.iter().any(|b| b.vertex.is_some())
            && seam_points(face).contains(apex)
        {
            let vertex = self.vertex(apex);
            if !self.scratch {
                let p = self.model.vertex(vertex).point;
                self.fact(GeomKey::Point(p), GeomKey::Surface(surface));
            }
            loops.push(self.model.add_loop(Loop::Vertex(vertex)));
        }
        if !self.scratch && points > 0 && (points > 1 || !loops.is_empty()) {
            return Err(unsupported(
                &face.name,
                "VERTEX_LOOP (beside other bounds of the face)",
            ));
        }
        if loops.is_empty() && !self.scratch && !(closed && !face.loops.is_empty()) {
            return Err(unsupported(&face.name, "ADVANCED_FACE (no bounds)"));
        }
        let id = self.model.add_named_face(
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
        );
        if matches!(face.surface, Surface::BSpline(_)) {
            let _ = catallax_topology::RationalPatchFace::of_face(self.model, id, &mut self.dec)
                .map_err(|error| {
                    crate::spline_build::refused(
                        &face.name,
                        "ADVANCED_FACE",
                        "boundary certification",
                        error,
                    )
                })?;
            let curves: Vec<_> = self
                .model
                .face(id)
                .loops
                .iter()
                .flat_map(|loop_id| match self.model.loop_(*loop_id) {
                    Loop::Coedges(coedges) => coedges
                        .iter()
                        .map(|coedge| self.model.edge(self.model.coedge(*coedge).edge).curve)
                        .collect::<Vec<_>>(),
                    Loop::Vertex(_) => Vec::new(),
                })
                .collect();
            for curve in curves {
                self.fact(GeomKey::Curve(curve), GeomKey::Surface(surface));
            }
        }
        Ok(id)
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

/// Seam-free closed chains of one bound, for callers outside the builder.
pub(crate) fn chains<'e>(
    face: &str,
    orientation: bool,
    edges: &'e [OrientedEdge],
) -> Result<Vec<Vec<(&'e OrientedEdge, bool)>>, ImportRefusal> {
    Builder::chains(face, orientation, edges)
}

/// End points of the seams of a face: edges used once each way in one bound.
pub(crate) fn seam_points(face: &StepFace) -> Vec<Vec3> {
    let mut out: Vec<Vec3> = Vec::new();
    for bound in &face.loops {
        for (i, e) in bound.edges.iter().enumerate() {
            let seam = bound
                .edges
                .iter()
                .enumerate()
                .any(|(j, u)| j != i && u.edge == e.edge && u.forward != e.forward);
            if seam {
                for p in [&e.start.point, &e.end.point] {
                    if !out.contains(p) {
                        out.push(p.clone());
                    }
                }
            }
        }
    }
    out
}

/// The apex of a conical face, where a vertex loop may bound it.
fn cone_apex(face: &StepFace) -> Result<Option<Vec3>, ImportRefusal> {
    Ok(match &face.surface {
        Surface::Cone45(p, r) => {
            let axis = unit(&p.axis).ok_or_else(|| {
                unsupported(&face.name, "CONICAL_SURFACE (axis length is not rational)")
            })?;
            Some(std::array::from_fn(|i| &p.origin[i] - &(&axis[i] * r)))
        }
        Surface::Cone(..) => crate::cone::rim_cone(face)?.map(|cone| cone.apex),
        _ => None,
    })
}

fn on_plane(origin: &Vec3, normal: &Vec3, p: &Vec3) -> bool {
    dot(&sub(p, origin), normal).is_zero()
}

fn on_axis(placement: &Placement, p: &Vec3) -> bool {
    is_zero(&cross(&sub(p, &placement.origin), &placement.axis))
}

/// Squared axial height of `p` above a torus or sphere placement, and its
/// squared distance from the centre.
fn height_and_distance(placement: &Placement, p: &Vec3) -> (Rational, Rational) {
    let d = sub(p, &placement.origin);
    let along = dot(&d, &placement.axis);
    (
        &(&along * &along) / &dot(&placement.axis, &placement.axis),
        dot(&d, &d),
    )
}

/// Exact test that a point lies on a toroidal surface with major radius `major`
/// and minor radius `minor`: (rho - major)^2 + h^2 = minor^2, squared out.
fn on_torus(placement: &Placement, major: &Rational, minor: &Rational, p: &Vec3) -> bool {
    let (h2, d2) = height_and_distance(placement, p);
    let rho2 = &d2 - &h2;
    let k = &(&(&rho2 + &(major * major)) + &h2) - &(minor * minor);
    k.sign() != catallax_numbers::Sign::Negative
        && &(&Rational::from_i64(4) * &(major * major)) * &rho2 == &k * &k
}

/// Exact test that a point lies on a spherical or toroidal surface.
fn on_closed_surface(surface: &Surface, p: &Vec3) -> bool {
    match surface {
        Surface::Sphere(c, r) => {
            let d = sub(p, &c.origin);
            dot(&d, &d) == r * r
        }
        Surface::Torus(c, major, minor) => on_torus(c, major, minor, p),
        _ => false,
    }
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
        (Surface::Cone(..), Curve::Circle(c, cr)) => {
            // The built cone passes through the face's rims; check exactly.
            let cone = crate::cone::rim_cone(face)?.expect("conical face");
            let h = dot(&sub(&c.origin, &cone.apex), &cone.axis);
            is_zero(&cross(&c.axis, &cone.axis))
                && is_zero(&cross(&sub(&c.origin, &cone.apex), &cone.axis))
                && h.sign() == catallax_numbers::Sign::Positive
                && *cr == &cone.slope * &h
        }
        (Surface::Sphere(p, r), Curve::Circle(c, cr)) => {
            // The circle's plane cuts the sphere in a circle about the foot of
            // the centre; its radius closes the right triangle.
            let d = sub(&c.origin, &p.origin);
            is_zero(&cross(&d, &c.axis)) && &dot(&d, &d) + &(cr * cr) == r * r
        }
        (Surface::Torus(p, major, minor), Curve::Circle(c, cr))
            if is_zero(&cross(&c.axis, &p.axis)) =>
        {
            // A parallel: centred on the axis, at height h, radius major +- sqrt(minor^2 - h^2).
            let (h2, _) = height_and_distance(p, &c.origin);
            let offset = cr - major;
            on_axis(p, &c.origin) && &(&offset * &offset) + &h2 == minor * minor
        }
        (Surface::Torus(p, major, minor), Curve::Circle(c, cr))
            if dot(&c.axis, &p.axis).is_zero()
                && dot(&c.axis, &sub(&c.origin, &p.origin)).is_zero() =>
        {
            // A meridian: centred on the spine circle, in a plane through the axis.
            let (h2, d2) = height_and_distance(p, &c.origin);
            h2.is_zero() && d2 == major * major && cr == minor
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
        Curve::Ellipse(..) | Curve::BSpline(_) => false,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_spline_incidence_does_not_certify_the_off_trim_span() {
        use catallax_intersect::{CurveSurface, intersect_curve_surface};
        use catallax_numbers::{Decider, Outcome, Real};
        use catallax_shapes::spline_tools::RationalSplineCurve;

        let mut dec = Decider::default();
        let controls = [v(["0", "0", "0"]), v(["1", "0", "0"]), v(["2", "1", "0"])];
        let geometry = RationalSplineCurve::new(
            1,
            ["0", "0", "1", "2", "2"].map(q).to_vec(),
            controls.to_vec(),
            vec![Rational::one(); 3],
            &dec,
        )
        .unwrap();
        let corners = [
            controls[0].clone(),
            controls[1].clone(),
            v(["1", "0", "1"]),
            v(["0", "0", "1"]),
        ];
        let edges = (0..4)
            .map(|i| {
                let start = corners[i].clone();
                let end = corners[(i + 1) % 4].clone();
                OrientedEdge {
                    edge: format!("edge-{i}"),
                    start: crate::Vertex {
                        name: format!("v-{i}"),
                        point: start.clone(),
                    },
                    end: crate::Vertex {
                        name: format!("v-{}", (i + 1) % 4),
                        point: end.clone(),
                    },
                    curve_same_sense: true,
                    forward: true,
                    curve: if i == 0 {
                        Curve::BSpline(crate::StepSplineCurve {
                            geometry: geometry.clone(),
                            closed: Some(false),
                            self_intersect: Some(false),
                        })
                    } else {
                        Curve::Line(start.clone(), sub(&end, &start))
                    },
                }
            })
            .collect();
        let face = StepFace {
            name: "finite-trim".into(),
            surface: Surface::Plane(placement(
                ["0", "0", "0"],
                ["0", "-1", "0"],
                ["1", "0", "0"],
            )),
            same_sense: true,
            loops: vec![crate::FaceLoop {
                outer: true,
                orientation: true,
                vertex: None,
                edges,
            }],
        };
        let mut model = Model::new();
        let face_id = Builder::new(&mut model, false)
            .face(&face, Vec::new())
            .unwrap();
        let surface = model.face(face_id).surface;
        let Loop::Coedges(coedges) = model.loop_(model.face(face_id).loops[0]) else {
            panic!("edge loop")
        };
        let curve = model.edge(model.coedge(coedges[0]).edge).curve;
        let GeomCurve::BSpline(retained) = model.geom().curve(curve) else {
            panic!("retained spline")
        };
        let off_trim = retained
            .point_at(&Real::from_rational(q("1.5")), &mut dec)
            .unwrap();
        let point = model.geom_mut().add_point(PointDef::Explicit(off_trim));
        assert!(matches!(
            model.geom().side(surface, point, &mut dec),
            Outcome::Proven(_)
        ));
        let query = intersect_curve_surface(model.geom_mut(), curve, surface, &mut dec);
        assert!(
            !matches!(&query, Ok(answer) if answer.is_exact() && matches!(&answer.value, CurveSurface::Contained)),
            "the full retained curve leaves the plane: {query:?}"
        );
        assert!(!model.geom().facts().holds(
            FactKind::Incident,
            GeomKey::Curve(curve),
            GeomKey::Surface(surface)
        ));
        assert!(
            query.is_err(),
            "unsupported full-curve query must refuse: {query:?}"
        );
    }

    fn q(s: &str) -> Rational {
        Rational::parse(s).unwrap()
    }
    fn v(c: [&str; 3]) -> Vec3 {
        c.map(q)
    }
    fn placement(origin: [&str; 3], axis: [&str; 3], ref_dir: [&str; 3]) -> Placement {
        Placement {
            origin: v(origin),
            axis: v(axis),
            ref_dir: v(ref_dir),
        }
    }
    fn circle(
        origin: [&str; 3],
        axis: [&str; 3],
        ref_dir: [&str; 3],
        radius: &str,
    ) -> OrientedEdge {
        let p = placement(origin, axis, ref_dir);
        let start: Vec3 = std::array::from_fn(|i| &p.origin[i] + &(&p.ref_dir[i] * &q(radius)));
        let vertex = crate::Vertex {
            name: "#v".into(),
            point: start,
        };
        OrientedEdge {
            edge: "#e".into(),
            start: vertex.clone(),
            end: vertex,
            curve_same_sense: true,
            curve: Curve::Circle(p, q(radius)),
            forward: true,
        }
    }
    fn on(surface: Surface, edge: &OrientedEdge) -> Option<bool> {
        let face = StepFace {
            name: "#f".into(),
            surface,
            same_sense: true,
            loops: Vec::new(),
        };
        match check_on_surface(&face, edge) {
            Ok(()) => Some(true),
            Err(ImportRefusal::NotSolid { .. }) => Some(false),
            Err(_) => None,
        }
    }

    #[test]
    fn circles_on_spheres_close_the_right_triangle_exactly() {
        // Axis (0,0,2) is not unit; the circle at height 3 has radius 4 on a sphere of 5.
        let sphere = || {
            Surface::Sphere(
                placement(["1", "1", "1"], ["0", "0", "2"], ["1", "0", "0"]),
                q("5"),
            )
        };
        assert_eq!(
            on(
                sphere(),
                &circle(["1", "1", "4"], ["0", "0", "1"], ["1", "0", "0"], "4")
            ),
            Some(true)
        );
        assert_eq!(
            on(
                sphere(),
                &circle(["1", "1", "4"], ["0", "0", "1"], ["1", "0", "0"], "4.001")
            ),
            Some(false)
        );
        // Tilted circle planes through an off-centre foot are off the sphere.
        assert_eq!(
            on(
                sphere(),
                &circle(["1", "4", "1"], ["0", "0", "1"], ["1", "0", "0"], "4")
            ),
            Some(false)
        );
        assert_eq!(
            on(
                sphere(),
                &circle(["1", "4", "1"], ["0", "3", "0"], ["1", "0", "0"], "4")
            ),
            Some(true)
        );
    }

    #[test]
    fn torus_parallels_and_meridians_are_checked_exactly() {
        let torus = || {
            Surface::Torus(
                placement(["0", "0", "0"], ["0", "0", "3"], ["1", "0", "0"]),
                q("5"),
                q("5/2"),
            )
        };
        // Parallels at height 3/2: radius 5 +- 2.
        for (h, radius, ok) in [
            ("3/2", "7", true),
            ("3/2", "3", true),
            ("-3/2", "7", true),
            ("3/2", "7.001", false),
            ("0", "15/2", true),
            ("0", "5/2", true),
        ] {
            assert_eq!(
                on(
                    torus(),
                    &circle(["0", "0", h], ["0", "0", "1"], ["1", "0", "0"], radius)
                ),
                Some(ok),
                "{h} {radius}"
            );
        }
        // A parallel-oriented circle off the axis is not on the torus.
        assert_eq!(
            on(
                torus(),
                &circle(["1", "0", "0"], ["0", "0", "1"], ["1", "0", "0"], "7")
            ),
            Some(false)
        );
        // Meridians: on the spine circle, in a plane through the axis.
        assert_eq!(
            on(
                torus(),
                &circle(["3", "4", "0"], ["-4", "3", "0"], ["3", "4", "0"], "5/2")
            ),
            Some(true)
        );
        assert_eq!(
            on(
                torus(),
                &circle(["3", "4", "0"], ["-4", "3", "0"], ["3", "4", "0"], "2")
            ),
            Some(false)
        );
        assert_eq!(
            on(
                torus(),
                &circle(["3", "4", "1"], ["-4", "3", "0"], ["3", "4", "0"], "5/2")
            ),
            Some(false)
        );
        // Any other circle orientation is outside the supported pair set.
        assert_eq!(
            on(
                torus(),
                &circle(["3", "4", "0"], ["1", "0", "1"], ["0", "1", "0"], "5/2")
            ),
            None
        );
    }

    #[test]
    fn points_on_closed_surfaces_are_exact() {
        let torus = Surface::Torus(
            placement(["0", "0", "0"], ["0", "0", "1"], ["1", "0", "0"]),
            q("5"),
            q("5/2"),
        );
        assert!(on_closed_surface(&torus, &v(["7", "0", "3/2"])));
        assert!(on_closed_surface(&torus, &v(["0", "-15/2", "0"])));
        assert!(!on_closed_surface(&torus, &v(["0", "0", "0"])));
        assert!(!on_closed_surface(&torus, &v(["5", "0", "0"])));
        let sphere = Surface::Sphere(
            placement(["0", "0", "0"], ["0", "0", "1"], ["1", "0", "0"]),
            q("5"),
        );
        assert!(on_closed_surface(&sphere, &v(["3", "0", "-4"])));
        assert!(!on_closed_surface(&sphere, &v(["3", "0", "4.001"])));
    }
}
