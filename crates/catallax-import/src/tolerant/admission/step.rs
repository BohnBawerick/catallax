use std::collections::BTreeMap;

use catallax_numbers::{Decider, Rational};

use crate::geometry_source::{ParsedStep, Topology};
use crate::part21::{Entity, Graph, Param};
use crate::{Curve, Placement, Surface, Vec3};

use super::{Error, TolerantSolid, invalid};

struct Builder<'a> {
    parsed: &'a ParsedStep,
    changes: Graph,
    next: u64,
}
impl Builder<'_> {
    fn add(&mut self, kind: &str, parameters: Vec<Param>) -> Result<u64, Error> {
        let id = self.next;
        self.next = self.next.checked_add(1).ok_or_else(|| invalid("import.tolerant-broken-reference", "STEP entity index overflow"))?;
        self.set(id, kind, parameters);
        Ok(id)
    }
    fn set(&mut self, id: u64, kind: &str, parameters: Vec<Param>) {
        self.changes.insert(id, Entity { parts: vec![(kind.into(), parameters)] });
    }
    fn existing(&self, name: &str, expected: &str) -> Option<u64> {
        let id = name.strip_prefix('#')?.parse().ok()?;
        self.parsed.source_graph().get(&id)?.part(expected)?;
        Some(id)
    }
    fn point(&mut self, position: &Vec3, scale: &Rational) -> Result<u64, Error> {
        let coordinates = position.iter().map(|coordinate| decimal(&(coordinate / scale))).collect::<Result<Vec<_>, _>>()?;
        self.add("CARTESIAN_POINT", vec![text(""), Param::List(coordinates)])
    }
    fn direction(&mut self, direction: &Vec3) -> Result<u64, Error> {
        self.add("DIRECTION", vec![text(""), Param::List(direction.iter().map(decimal).collect::<Result<_, _>>()?)])
    }
    fn placement(&mut self, placement: &Placement, scale: &Rational) -> Result<u64, Error> {
        let origin = self.point(&placement.origin, scale)?;
        let axis = self.direction(&placement.axis)?;
        let reference = self.direction(&placement.ref_dir)?;
        self.add("AXIS2_PLACEMENT_3D", vec![text(""), reference_to(origin), reference_to(axis), reference_to(reference)])
    }
    fn source_curve(&self, name: &str, dec: &mut Decider) -> Result<Option<u64>, Error> {
        let Some(entity) = name.strip_prefix('#').and_then(|value| value.parse().ok()) else { return Ok(None); };
        Ok(match self.parsed.topology(entity, dec)? { Topology::Edge { curve, .. } => curve, _ => None })
    }
    fn source_surface(&self, name: &str, dec: &mut Decider) -> Result<Option<u64>, Error> {
        let Some(entity) = name.strip_prefix('#').and_then(|value| value.parse().ok()) else { return Ok(None); };
        Ok(match self.parsed.topology(entity, dec)? { Topology::Face { surface, .. } => surface, _ => None })
    }
    fn shell_id(&self, mut id: u64, dec: &mut Decider) -> Result<u64, Error> {
        for _ in 0..32 {
            match self.parsed.topology(id, dec)? {
                Topology::OrientedShell(parent, _) => id = parent,
                Topology::Shell { .. } => return Ok(id),
                _ => return Err(invalid("import.tolerant-broken-reference", "source shell reference is not a shell")),
            }
        }
        Err(invalid("import.tolerant-broken-reference", "cyclic source shell reference"))
    }
}

fn text(value: &str) -> Param { Param::Str(value.into()) }
fn reference_to(id: u64) -> Param { Param::Ref(id) }
fn boolean(value: bool) -> Param { Param::Enum(if value { "T" } else { "F" }.into()) }
fn references(values: impl IntoIterator<Item = u64>) -> Param { Param::List(values.into_iter().map(reference_to).collect()) }
fn decimal(value: &Rational) -> Result<Param, Error> {
    let number = value.to_f64_nearest();
    if !number.is_finite() { return Err(invalid("import.tolerant-broken-geometry", "non-finite tolerant STEP coordinate")); }
    Ok(Param::Number(format!("{number:.17e}")))
}

pub(super) fn rebuild(
    parsed: &ParsedStep,
    bodies: &[TolerantSolid],
    scale: &Rational,
    dec: &mut Decider,
) -> Result<Graph, Error> {
    let next = parsed.entities().last().unwrap_or(0).checked_add(1)
        .ok_or_else(|| invalid("import.tolerant-broken-reference", "STEP entity index overflow"))?;
    let mut builder = Builder { parsed, changes: Graph::new(), next };
    let mut expanded_roots = BTreeMap::<u64, Vec<u64>>::new();
    for solid in bodies {
        let body = &solid.body;
        let root: u64 = body.name.strip_prefix('#').and_then(|value| value.parse().ok())
            .ok_or_else(|| invalid("import.tolerant-broken-reference", "native body lacks a source root"))?;
        let old_shells = match parsed.topology(root, dec)? {
            Topology::Solid { outer, voids } => std::iter::once(outer).chain(voids).collect(),
            Topology::SurfaceModel(shells) | Topology::FaceModel(shells) => shells,
            _ => return Err(invalid("import.tolerant-broken-topology", "native body root has no finite shells")),
        };
        let mut vertices = Vec::with_capacity(body.vertices.len());
        for vertex in &body.vertices {
            dec.charge_work(1)?;
            let point = builder.point(&vertex.source.point, scale)?;
            let id = match builder.existing(&vertex.source.name, "VERTEX_POINT") {
                Some(id) => { builder.set(id, "VERTEX_POINT", vec![text(&vertex.source.name), reference_to(point)]); id }
                None => builder.add("VERTEX_POINT", vec![text(&vertex.source.name), reference_to(point)])?,
            };
            vertices.push(id);
        }
        let mut active = vec![false; body.edges.len()];
        for shell in &body.shells { for face in &shell.faces { for boundary in &face.loops { for usage in &boundary.edges { active[usage.edge] = true; } } } }
        let mut edges = vec![None; body.edges.len()];
        for (index, edge) in body.edges.iter().enumerate() {
            if !active[index] { continue; }
            dec.charge_work(1)?;
            let curve = match builder.source_curve(&edge.name, dec)? {
                Some(curve) => curve,
                None => match &edge.curve {
                    Curve::Line(origin, direction) => {
                        let origin = builder.point(origin, scale)?;
                        let direction = builder.direction(direction)?;
                        let vector = builder.add("VECTOR", vec![text(""), reference_to(direction), Param::Number("1.".into())])?;
                        builder.add("LINE", vec![text(""), reference_to(origin), reference_to(vector)])?
                    }
                    _ => return Err(invalid("import.tolerant-broken-reference", "native edge has no exportable source carrier")),
                },
            };
            let parameters = vec![text(&edge.name), reference_to(vertices[edge.vertices[0]]), reference_to(vertices[edge.vertices[1]]), reference_to(curve), boolean(edge.curve_same_sense)];
            let entity = match builder.existing(&edge.name, "EDGE_CURVE") {
                Some(entity) => { builder.set(entity, "EDGE_CURVE", parameters); entity }
                None => builder.add("EDGE_CURVE", parameters)?,
            };
            edges[index] = Some(entity);
        }
        let mut shells = Vec::with_capacity(body.shells.len());
        for (shell_index, shell) in body.shells.iter().enumerate() {
            let mut faces = Vec::with_capacity(shell.faces.len());
            for face in &shell.faces {
                dec.charge_work(1)?;
                let surface = match builder.source_surface(&face.name, dec)? {
                    Some(surface) => surface,
                    None => match &face.surface {
                        Surface::Plane(placement) => {
                            let placement = builder.placement(placement, scale)?;
                            builder.add("PLANE", vec![text(""), reference_to(placement)])?
                        }
                        _ => return Err(invalid("import.tolerant-broken-reference", "native face has no exportable source surface")),
                    },
                };
                let mut bounds = Vec::with_capacity(face.loops.len());
                for boundary in &face.loops {
                    let loop_id = if let Some(vertex) = boundary.vertex {
                        builder.add("VERTEX_LOOP", vec![text(""), reference_to(vertices[vertex])])?
                    } else {
                        let mut uses = Vec::with_capacity(boundary.edges.len());
                        for usage in &boundary.edges {
                            let edge = edges[usage.edge].ok_or_else(|| invalid("import.tolerant-broken-reference", "native trim references an inactive edge"))?;
                            uses.push(builder.add("ORIENTED_EDGE", vec![text(""), Param::Derived, Param::Derived, reference_to(edge), boolean(usage.forward)])?);
                        }
                        builder.add("EDGE_LOOP", vec![text(""), references(uses)])?
                    };
                    bounds.push(builder.add(if boundary.outer { "FACE_OUTER_BOUND" } else { "FACE_BOUND" }, vec![text(""), reference_to(loop_id), boolean(boundary.orientation)])?);
                }
                let parameters = vec![text(&face.name), references(bounds), reference_to(surface), boolean(face.same_sense)];
                let entity = match builder.existing(&face.name, "ADVANCED_FACE").or_else(|| builder.existing(&face.name, "FACE_SURFACE")).or_else(|| builder.existing(&face.name, "FACE")) {
                    Some(entity) => { builder.set(entity, "ADVANCED_FACE", parameters); entity }
                    None => builder.add("ADVANCED_FACE", parameters)?,
                };
                faces.push(entity);
            }
            let parameters = vec![text(&shell.name), references(faces)];
            let entity = match old_shells.get(shell_index) {
                Some(&old) => { let entity = builder.shell_id(old, dec)?; builder.set(entity, "CLOSED_SHELL", parameters); entity }
                None => builder.add("CLOSED_SHELL", parameters)?,
            };
            shells.push(entity);
        }
        if body.declared_solid && shells.len() > 1 {
            let mut voids = Vec::with_capacity(shells.len() - 1);
            for &shell in &shells[1..] {
                voids.push(builder.add("ORIENTED_CLOSED_SHELL", vec![text(""), Param::Derived, reference_to(shell), boolean(true)])?);
            }
            builder.set(root, "BREP_WITH_VOIDS", vec![text(&body.name), reference_to(shells[0]), references(voids)]);
        } else {
            let mut roots = Vec::with_capacity(shells.len());
            for (index, shell) in shells.into_iter().enumerate() {
                let parameters = vec![text(&body.name), reference_to(shell)];
                roots.push(if index == 0 { builder.set(root, "MANIFOLD_SOLID_BREP", parameters); root } else { builder.add("MANIFOLD_SOLID_BREP", parameters)? });
            }
            if roots.len() > 1 { expanded_roots.insert(root, roots); }
        }
    }
    // Preserve the original occurrence/assembly graph. Only representation items
    // which held a promoted multi-shell surface set need item-list expansion.
    for (&entity, source) in parsed.source_graph() {
        for (kind, parameters) in &source.parts {
            if !kind.ends_with("REPRESENTATION") || kind == "DEFINITIONAL_REPRESENTATION" { continue; }
            let Some(Param::List(items)) = parameters.get(1) else { continue; };
            let mut changed = false;
            let mut replacement = Vec::new();
            for item in items {
                if let Param::Ref(root) = item
                    && let Some(roots) = expanded_roots.get(root)
                {
                    replacement.extend(roots.iter().copied().map(reference_to));
                    changed = true;
                } else { replacement.push(item.clone()); }
            }
            let promote = matches!(kind.as_str(), "MANIFOLD_SURFACE_SHAPE_REPRESENTATION" | "GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION")
                && items.iter().any(|item| matches!(item, Param::Ref(root) if bodies.iter().any(|body| body.name().strip_prefix('#').and_then(|name|name.parse::<u64>().ok()) == Some(*root))));
            if changed || promote {
                let mut value = source.clone();
                for (name, parameters) in &mut value.parts {
                    if name == kind {
                        parameters[1] = Param::List(replacement.clone());
                        if promote { *name = "ADVANCED_BREP_SHAPE_REPRESENTATION".into(); }
                    }
                }
                builder.changes.insert(entity, value);
            }
        }
    }
    Ok(builder.changes)
}
