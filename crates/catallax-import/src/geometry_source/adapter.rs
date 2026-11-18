//! Typed source topology mapped to the existing per-body admission input.
use super::*;
use crate::{
    Curve as C, Face, FaceLoop, OrientedEdge, Placement, StepSplineCurve, Surface as S,
    SurfaceShell, Vertex,
};
use catallax_numbers::Rational;
use std::collections::BTreeSet;

pub(crate) struct Definition {
    pub scale: Rational,
    pub bodies: Vec<BodyDefinition>,
}
pub(crate) struct BodyDefinition {
    pub name: String,
    pub solid: bool,
    pub shells: Vec<SurfaceShell>,
}
struct Adapter<'a> {
    source: &'a ParsedStep,
    scale: Rational,
    dec: &'a mut Decider,
    active: BTreeSet<EntityId>,
    native: bool,
}
fn invalid(id: EntityId, detail: &str) -> ImportRefusal {
    malformed(id, detail)
}
impl ParsedStep {
    pub(crate) fn native_inventory(
        &self,
        dec: &mut Decider,
    ) -> Result<(Rational, Vec<ModelDeclaration>, Vec<ContextUnits>), ImportRefusal> {
        let mut declarations = self.model_declarations(dec)?;
        declarations.sort_unstable_by_key(|model| {
            (model.kind != ModelDeclarationKind::Solid, model.entity)
        });
        if declarations.is_empty() {
            return Err(invalid(0, "no source body declaration"));
        }
        let mut units = Vec::new();
        let scale = if let Some(contexts) = self.model_unit_contexts(&declarations, dec)? {
            let mut selected = None;
            for context in contexts {
                let context = self.context_units(context, dec)?;
                let scale =
                    context
                        .length_mm
                        .as_ref()
                        .ok_or_else(|| ImportRefusal::UnknownUnit {
                            found: "prototype context has no length unit".into(),
                        })?;
                if selected.as_ref().is_some_and(|prior| prior != scale) {
                    return Err(ImportRefusal::UnknownUnit {
                        found: "conflicting length units".into(),
                    });
                }
                if selected.is_none() {
                    selected = Some(scale.clone());
                }
                units.push(context);
            }
            selected.ok_or_else(|| ImportRefusal::UnknownUnit {
                found: "length unit".into(),
            })?
        } else {
            let scale = crate::recognize::source_length_scale(&self.graph, dec)?;
            for entity in self.entities() {
                if self
                    .entity_types(entity)?
                    .contains(&"GLOBAL_UNIT_ASSIGNED_CONTEXT")
                    && self.representation_context(entity)?.dimension == Some(3)
                {
                    units.push(self.context_units(entity, dec)?);
                }
            }
            scale
        };
        Ok((scale, declarations, units))
    }

    pub(crate) fn retained_body(&self, entity: EntityId, scale: &Rational) -> RetainedGeometry {
        RetainedGeometry {
            source: self.clone(),
            entity,
            length_mm: scale.clone(),
        }
    }

    pub(crate) fn body_definition(
        &self,
        entity: EntityId,
        solid: bool,
        scale: &Rational,
        native: bool,
        dec: &mut Decider,
    ) -> Result<BodyDefinition, ImportRefusal> {
        Adapter {
            source: self,
            scale: scale.clone(),
            dec,
            active: BTreeSet::new(),
            native,
        }
        .model(entity, solid)
    }

    pub(crate) fn retained_definition(
        &self,
        dec: &mut Decider,
    ) -> Result<Definition, ImportRefusal> {
        let scale = crate::recognize::source_length_scale(&self.graph, dec)?;
        // A prototype is not an occurrence. Do not publish unplaced duplicates.
        for id in self.entities() {
            let e = self.entity(id)?;
            dec.charge_work(1)
                .map_err(|_| ImportRefusal::ResourceLimit)?;
            if e.parts.iter().any(|(kind, _)| {
                kind == "NEXT_ASSEMBLY_USAGE_OCCURRENCE"
                    || (crate::recognize::unsupported_solid_model(kind)
                        && !matches!(kind.as_str(), "BREP_WITH_VOIDS" | "FACETED_BREP"))
            }) {
                return Err(ImportRefusal::NotSolid {
                    detail: format!(
                        "#{id}: complete source-model interpretation is required before admission"
                    ),
                });
            }
            if e.parts.iter().any(|(kind, _)| {
                matches!(
                    kind.as_str(),
                    "MAPPED_ITEM"
                        | "REPRESENTATION_MAP"
                        | "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"
                        | "ITEM_DEFINED_TRANSFORMATION"
                )
            }) {
                return Err(ImportRefusal::NotSolid {
                    detail: format!(
                        "#{id}: occurrence placement must be resolved before body admission"
                    ),
                });
            }
        }
        let declarations = self.model_declarations(dec)?;
        let mut adapter = Adapter {
            source: self,
            scale: scale.clone(),
            dec,
            active: BTreeSet::new(),
            native: false,
        };
        let mut bodies = Vec::new();
        // Keep the established solid-then-surface inventory ordering.
        for solid in [true, false] {
            for model in &declarations {
                if (model.kind == ModelDeclarationKind::Solid) != solid {
                    continue;
                }
                bodies.push(adapter.model(model.entity, solid)?);
            }
        }
        if bodies.is_empty() {
            return Err(invalid(0, "no finite body declaration"));
        }
        Ok(Definition { scale, bodies })
    }
}
impl Adapter<'_> {
    fn model(&mut self, entity: EntityId, solid: bool) -> Result<BodyDefinition, ImportRefusal> {
        let shells = match self.source.topology(entity, self.dec)? {
            Topology::Solid { outer, voids } => {
                let mut shells = vec![self.shell(outer)?];
                for inner in voids {
                    shells.push(self.shell(inner)?);
                }
                shells
            }
            Topology::SurfaceModel(shells) | Topology::FaceModel(shells) => shells
                .into_iter()
                .map(|id| self.shell(id))
                .collect::<Result<_, _>>()?,
            Topology::GeometricSet { .. } => {
                return Err(ImportRefusal::NotSolid {
                    detail: format!("#{entity}: geometric set has no finite face topology"),
                });
            }
            _ => return Err(invalid(entity, "expected a body declaration")),
        };
        if shells.is_empty() {
            return Err(invalid(entity, "body has no shells"));
        }
        Ok(BodyDefinition {
            name: format!("#{entity}"),
            solid,
            shells,
        })
    }
    fn enter(&mut self, id: EntityId) -> Result<(), ImportRefusal> {
        self.dec
            .charge_work(1)
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        if self.active.len() >= 128 || !self.active.insert(id) {
            return Err(invalid(
                id,
                "cyclic or excessively deep topology or geometry",
            ));
        }
        Ok(())
    }
    fn point(&mut self, id: EntityId) -> Result<crate::Vec3, ImportRefusal> {
        self.source
            .cartesian_point(id, &std::array::from_fn(|_| self.scale.clone()), self.dec)
            .map_err(mapping_error)
    }
    fn vertex(&mut self, id: EntityId) -> Result<Vertex, ImportRefusal> {
        let Topology::Vertex(point) = self.source.topology(id, self.dec)? else {
            return Err(invalid(id, "expected a vertex point"));
        };
        Ok(Vertex {
            name: format!("#{id}"),
            point: self.point(point)?,
        })
    }
    fn shell(&mut self, id: EntityId) -> Result<SurfaceShell, ImportRefusal> {
        self.enter(id)?;
        let result = match self.source.topology(id, self.dec)? {
            Topology::OrientedShell(parent, orientation) => {
                let mut shell = self.shell(parent)?;
                shell.name = format!("#{id}");
                if !orientation {
                    for face in &mut shell.faces {
                        reverse_face(face);
                    }
                }
                Ok(shell)
            }
            Topology::Shell { faces, closed } => Ok(SurfaceShell {
                name: format!("#{id}"),
                closed: closed.unwrap_or(false),
                faces: faces
                    .into_iter()
                    .map(|id| self.face(id))
                    .collect::<Result<_, _>>()?,
            }),
            _ => Err(invalid(id, "expected a shell or connected face set")),
        };
        self.active.remove(&id);
        result
    }
    fn face(&mut self, id: EntityId) -> Result<Face, ImportRefusal> {
        self.enter(id)?;
        let result = match self.source.topology(id, self.dec)? {
            Topology::OrientedFace(parent, orientation) => {
                let mut face = self.face(parent)?;
                face.name = format!("#{id}");
                if !orientation {
                    reverse_face(&mut face);
                }
                Ok(face)
            }
            Topology::Face {
                bounds,
                surface,
                same_sense,
            } => {
                let loops = bounds
                    .into_iter()
                    .map(|id| self.bound(id))
                    .collect::<Result<Vec<_>, _>>()?;
                let surface = match surface {
                    Some(id) => self.surface(id, &loops)?,
                    None => polygon_plane(id, &loops)?,
                };
                Ok(Face {
                    name: format!("#{id}"),
                    surface,
                    same_sense: same_sense.unwrap_or(true),
                    loops,
                })
            }
            _ => Err(invalid(id, "expected a face")),
        };
        self.active.remove(&id);
        result
    }
    fn bound(&mut self, id: EntityId) -> Result<FaceLoop, ImportRefusal> {
        let Topology::FaceBound {
            loop_id,
            orientation,
            outer,
        } = self.source.topology(id, self.dec)?
        else {
            return Err(invalid(id, "expected a face bound"));
        };
        let mut boundary = FaceLoop {
            outer,
            orientation,
            edges: Vec::new(),
            vertex: None,
        };
        match self.source.topology(loop_id, self.dec)? {
            Topology::EdgeLoop(edges) => {
                for edge in edges {
                    boundary.edges.push(self.edge(edge)?);
                }
            }
            Topology::VertexLoop(vertex) => boundary.vertex = Some(self.vertex(vertex)?),
            Topology::PolyLoop(points) => {
                if points.len() < 3 {
                    return Err(invalid(loop_id, "polygon requires three points"));
                }
                let vertices = points
                    .into_iter()
                    .map(|id| {
                        Ok(Vertex {
                            name: format!("#{id}"),
                            point: self.point(id)?,
                        })
                    })
                    .collect::<Result<Vec<_>, ImportRefusal>>()?;
                for i in 0..vertices.len() {
                    let start = vertices[i].clone();
                    let end = vertices[(i + 1) % vertices.len()].clone();
                    let forward = start.name < end.name;
                    let (start, end) = if forward { (start, end) } else { (end, start) };
                    let direction = crate::exact::sub(&end.point, &start.point);
                    boundary.edges.push(OrientedEdge {
                        edge: format!("{}:{}", start.name, end.name),
                        curve: C::Line(start.point.clone(), direction),
                        start,
                        end,
                        curve_same_sense: true,
                        forward,
                    });
                }
            }
            _ => return Err(invalid(loop_id, "expected a finite face loop")),
        }
        Ok(boundary)
    }
    fn edge(&mut self, id: EntityId) -> Result<OrientedEdge, ImportRefusal> {
        self.enter(id)?;
        let result = match self.source.topology(id, self.dec)? {
            Topology::OrientedEdge(parent, orientation) => {
                let mut edge = self.edge(parent)?;
                edge.forward = edge.forward == orientation;
                Ok(edge)
            }
            Topology::Edge {
                vertices,
                curve,
                same_sense,
            } => {
                let start = self.vertex(vertices[0])?;
                let end = self.vertex(vertices[1])?;
                let curve = match curve {
                    Some(id) => self.curve(id)?,
                    None => C::Retained(self.retained(id)),
                };
                Ok(OrientedEdge {
                    edge: format!("#{id}"),
                    start,
                    end,
                    curve,
                    curve_same_sense: same_sense.unwrap_or(true),
                    forward: true,
                })
            }
            _ => Err(invalid(id, "expected an edge")),
        };
        self.active.remove(&id);
        result
    }
    fn retained(&self, id: EntityId) -> RetainedGeometry {
        RetainedGeometry {
            source: self.source.clone(),
            entity: id,
            length_mm: self.scale.clone(),
        }
    }
    fn curve(&mut self, id: EntityId) -> Result<C, ImportRefusal> {
        if self.native {
            self.source.entity(id)?;
            return Ok(C::Retained(self.retained(id)));
        }
        self.enter(id)?;
        let result = match self.source.geometry(id, self.dec)? {
            Geometry::Curve(Curve::SurfaceCurve(parent, _, _, _)) => self.curve(parent),
            Geometry::Curve(Curve::Polyline(_)) => Ok(C::BSpline(StepSplineCurve {
                geometry: self
                    .source
                    .spline_curve(id, &std::array::from_fn(|_| self.scale.clone()), self.dec)
                    .map_err(mapping_error)?,
                closed: None,
                self_intersect: None,
            })),
            Geometry::Curve(_) => {
                match crate::recognize::source_curve(&self.source.graph, id, &self.scale, self.dec)
                {
                    Ok(curve) => Ok(curve),
                    Err(ImportRefusal::Unsupported { .. }) => Ok(C::Retained(self.retained(id))),
                    Err(error) => Err(error),
                }
            }
            _ => Err(invalid(id, "expected a curve carrier")),
        };
        self.active.remove(&id);
        result
    }
    fn surface(&mut self, id: EntityId, loops: &[FaceLoop]) -> Result<S, ImportRefusal> {
        if self.native {
            self.source.entity(id)?;
            return Ok(S::Retained(self.retained(id)));
        }
        self.enter(id)?;
        let result = match self.source.geometry(id, self.dec)? {
            Geometry::Surface(Surface::Torus(_, major, minor)) if minor >= major => {
                Ok(S::Retained(self.retained(id)))
            }
            Geometry::Surface(Surface::Torus(placement, major, minor)) => Ok(S::Torus(
                crate::recognize::source_placement(
                    &self.source.graph,
                    placement,
                    &self.scale,
                    self.dec,
                )?,
                &major * &self.scale,
                &minor * &self.scale,
            )),
            Geometry::Surface(_) => {
                match crate::recognize::source_surface(
                    &self.source.graph,
                    id,
                    loops,
                    &self.scale,
                    self.dec,
                ) {
                    Ok(surface) => Ok(surface),
                    Err(ImportRefusal::Unsupported { .. }) => Ok(S::Retained(self.retained(id))),
                    Err(error) => Err(error),
                }
            }
            _ => Err(invalid(id, "expected a surface carrier")),
        };
        self.active.remove(&id);
        result
    }
}
fn mapping_error(error: MappingError) -> ImportRefusal {
    match error {
        MappingError::Source(error) => error,
        MappingError::Geometry(error) => ImportRefusal::NotSolid {
            detail: format!("{}: {}", error.code, error.detail),
        },
    }
}
fn reverse_face(face: &mut Face) {
    face.same_sense = !face.same_sense;
    for boundary in &mut face.loops {
        boundary.orientation = !boundary.orientation;
    }
}
fn polygon_plane(id: EntityId, loops: &[FaceLoop]) -> Result<S, ImportRefusal> {
    let boundary = loops
        .iter()
        .find(|b| b.outer)
        .or_else(|| loops.first())
        .ok_or_else(|| invalid(id, "faceted face has no boundary"))?;
    let points: Vec<_> = boundary
        .edges
        .iter()
        .map(|edge| {
            if edge.forward == boundary.orientation {
                &edge.start.point
            } else {
                &edge.end.point
            }
        })
        .collect();
    let origin = points
        .first()
        .ok_or_else(|| invalid(id, "faceted face has no polygon"))?;
    for i in 1..points.len().saturating_sub(1) {
        let x = crate::exact::sub(points[i], origin);
        let normal = crate::exact::cross(&x, &crate::exact::sub(points[i + 1], origin));
        if !crate::exact::is_zero(&normal) {
            if points
                .iter()
                .any(|p| !crate::exact::dot(&normal, &crate::exact::sub(p, origin)).is_zero())
            {
                return Err(invalid(id, "faceted polygon is not planar"));
            }
            return Ok(S::Plane(Placement {
                origin: (*origin).clone(),
                axis: normal,
                ref_dir: x,
            }));
        }
    }
    Err(invalid(id, "faceted polygon is collinear"))
}
