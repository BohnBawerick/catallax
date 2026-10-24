//! Exact selected normal fans of bounded convex planar skins.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, Origin, Role, StepId, ThickenSide, VertexName};
use catallax_shapes::{Circle, Curve, Cylinder, Line, Plane, PointDef, Sphere, Surface};
use catallax_topology::{BodyId, Edge, EdgeId, Model, Sense, Vertex, VertexId};

use crate::{
    cone_skin::reduced,
    extrude::{face, reverse},
    polyhedron::ConvexCell,
    profile::sign,
    sharp_thicken::interior_rim,
};

type EdgeUses = BTreeMap<(usize, usize), Vec<(usize, Sense)>>;
type ArcUse = (usize, usize, EdgeId, Sense, catallax_shapes::Vector3);

/// A bounded analytic normal sector, with names only on physical boundaries.
#[derive(Clone, Debug)]
pub struct NormalSector {
    /// Oriented nonpositive constraints and their physical origins.
    pub supports: Vec<(Surface, Vec<Origin>)>,
}

/// Selected faces and the complete incidence of their convex source skin.
pub struct PlanarThicken {
    cell: ConvexCell,
    selected: Vec<bool>,
    edges: EdgeUses,
    vertices: BTreeMap<usize, Vec<usize>>,
}
fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.convex-normal-fan", detail)
}
impl PlanarThicken {
    /// Retains exact input incidence before constructing either kernel's result.
    pub fn new(cell: ConvexCell, selected: Vec<bool>) -> Result<Self, Refusal> {
        if selected.len() != cell.faces.len() {
            return Err(bad("selection and planar skin differ"));
        }
        let mut edges: EdgeUses = BTreeMap::new();
        let mut vertices: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, polygon) in cell.polygons.iter().enumerate() {
            for j in 0..polygon.len() {
                let a = polygon[j];
                let b = polygon[(j + 1) % polygon.len()];
                edges
                    .entry((a.min(b), a.max(b)))
                    .or_default()
                    .push((i, if a < b { Sense::Same } else { Sense::Reversed }));
                vertices.entry(a).or_default().push(i);
            }
        }
        for uses in edges.values() {
            let [a, b] = uses.as_slice() else {
                return Err(bad("a convex source edge needs two incident faces"));
            };
            if a.1 == b.1 {
                return Err(bad("input edge uses have inconsistent orientation"));
            }
        }
        Ok(Self {
            cell,
            selected,
            edges,
            vertices,
        })
    }
    fn edge_name(&self, uses: &[(usize, Sense)]) -> Result<EdgeName, Refusal> {
        let origin = |i: usize| {
            self.cell.faces[i]
                .origins
                .first()
                .cloned()
                .ok_or_else(|| bad("a normal fan requires named input faces"))
        };
        EdgeName::new(origin(uses[0].0)?, origin(uses[1].0)?).map_err(|e| bad(&e.to_string()))
    }
    fn boundary_names(&self, vertex: usize, step: &StepId) -> Result<Vec<Origin>, Refusal> {
        let mut names = Vec::new();
        for ((a, b), uses) in &self.edges {
            if (*a == vertex || *b == vertex)
                && self.selected[uses[0].0] != self.selected[uses[1].0]
            {
                names.push(Origin {
                    step: step.clone(),
                    role: Role::Rim {
                        edge: self.edge_name(uses)?,
                    },
                });
            }
        }
        names.sort();
        names.dedup();
        Ok(names)
    }
    fn check_front(&self) -> Result<(), Refusal> {
        let Some(start) = self.selected.iter().position(|yes| *yes) else {
            return Err(bad("a normal fan needs a selected source face"));
        };
        let mut seen = vec![false; self.selected.len()];
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(i) = stack.pop() {
            for uses in self.edges.values() {
                if !uses.iter().any(|(f, _)| *f == i) {
                    continue;
                }
                for (j, _) in uses {
                    if self.selected[*j] && !seen[*j] {
                        seen[*j] = true;
                        stack.push(*j);
                    }
                }
            }
        }
        if self
            .selected
            .iter()
            .zip(seen)
            .any(|(yes, reached)| *yes && !reached)
        {
            return Err(bad(
                "disconnected selected normal fans require separate lump ownership",
            ));
        }
        for (v, incident) in &self.vertices {
            if incident.iter().all(|i| self.selected[*i]) {
                continue;
            }
            let sectors = self
                .edges
                .iter()
                .filter(|((a, b), uses)| {
                    (*a == *v || *b == *v) && uses.iter().all(|(i, _)| self.selected[*i])
                })
                .count();
            if sectors > 1 {
                return Err(bad(
                    "multiple edge fans at a boundary vertex need their exact end-fan arrangement",
                ));
            }
        }
        Ok(())
    }
    /// Finite face images and clipped edge/vertex sectors for independent field compilation.
    pub fn pieces(
        &self,
        thickness: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<(Vec<ConvexCell>, Vec<NormalSector>), Refusal> {
        if matches!(side, ThickenSide::Front) {
            self.check_front()?;
        }
        let selected_origins: Vec<_> = self
            .cell
            .faces
            .iter()
            .zip(&self.selected)
            .filter_map(|(f, yes)| yes.then_some(f.origins.as_slice()))
            .collect();
        let mut prisms = Vec::new();
        for (i, yes) in self.selected.iter().enumerate() {
            if *yes {
                let mut prism = self.cell.sharp_sweep(i, thickness, side, step, dec)?;
                for f in &mut prism.faces {
                    f.origins.retain(|o| !interior_rim(o, &selected_origins));
                }
                prisms.push(prism);
            }
        }
        let mut sectors = Vec::new();
        if matches!(side, ThickenSide::Back) {
            return Ok((prisms, sectors));
        }
        for ((a, b), uses) in &self.edges {
            if !uses.iter().all(|(i, _)| self.selected[*i]) {
                continue;
            }
            let axis = self.cell.vertices()[*b].minus(&self.cell.vertices()[*a]);
            let nf = &self.cell.faces[uses[0].0].plane.normal;
            let ng = &self.cell.faces[uses[1].0].plane.normal;
            let cylinder = Cylinder::new(
                self.cell.vertices()[*a].clone(),
                axis.clone(),
                thickness.clone(),
                nf.clone(),
                dec,
            )?;
            let mut supports = vec![(
                Surface::Cylinder(cylinder),
                vec![Origin {
                    step: step.clone(),
                    role: Role::Round {
                        edge: self.edge_name(uses)?,
                    },
                }],
            )];
            for (v, direction) in [(*a, -1), (*b, 1)] {
                supports.push((
                    Surface::Plane(Plane {
                        origin: self.cell.vertices()[v].clone(),
                        normal: axis.scaled(&Real::from_i64(direction)),
                        ref_dir: nf.clone(),
                    }),
                    self.boundary_names(v, step)?,
                ));
            }
            for (normal, neighbour) in [(nf, ng), (ng, nf)] {
                let tangent = neighbour
                    .minus(&normal.scaled(&(neighbour.dot(normal) / normal.norm_squared())));
                supports.push((
                    Surface::Plane(Plane {
                        origin: self.cell.vertices()[*a].clone(),
                        normal: tangent.scaled(&Real::from_i64(-1)),
                        ref_dir: axis.clone(),
                    }),
                    Vec::new(),
                ));
            }
            sectors.push(NormalSector { supports });
        }
        for (v, incident) in &self.vertices {
            if !incident.iter().all(|i| self.selected[*i]) {
                continue;
            }
            let source = &self.cell.faces[incident[0]].plane;
            let sphere = Sphere::new(
                self.cell.vertices()[*v].clone(),
                thickness.clone(),
                source.normal.clone(),
                source.ref_dir.clone(),
                dec,
            )?;
            let names = incident
                .iter()
                .filter_map(|i| self.cell.faces[*i].origins.first().cloned())
                .collect();
            let vertex = VertexName::new(names).map_err(|e| bad(&e.to_string()))?;
            let mut supports = vec![(
                Surface::Sphere(sphere),
                vec![Origin {
                    step: step.clone(),
                    role: Role::RoundCorner { vertex },
                }],
            )];
            for (a, b) in self.edges.keys() {
                let other = if a == v {
                    *b
                } else if b == v {
                    *a
                } else {
                    continue;
                };
                supports.push((
                    Surface::Plane(Plane {
                        origin: self.cell.vertices()[*v].clone(),
                        normal: self.cell.vertices()[other].minus(&self.cell.vertices()[*v]),
                        ref_dir: source.normal.clone(),
                    }),
                    Vec::new(),
                ));
            }
            sectors.push(NormalSector { supports });
        }
        Ok((prisms, sectors))
    }
    pub(crate) fn build(
        &self,
        model: &mut Model,
        thickness: &Real,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<BodyId, Refusal> {
        if matches!(side, ThickenSide::Back) {
            let (prisms, _) = self.pieces(thickness, side, step, dec)?;
            let mut result = catallax_topology::empty(model);
            for prism in prisms {
                let piece = prism.build(model);
                result = catallax_boolean::boolean(
                    model,
                    result,
                    piece,
                    catallax_boolean::BoolOp::Unite,
                    dec,
                )
                .map_err(|e| e.refusal)?
                .body;
            }
            return Ok(result);
        }
        self.check_front()?;
        let mut original = BTreeMap::new();
        let mut displaced = BTreeMap::new();
        let mut rails = BTreeMap::new();
        let mut base_edges = BTreeMap::new();
        let mut far_edges = BTreeMap::new();
        let mut faces = Vec::new();
        for (i, polygon) in self.cell.polygons.iter().enumerate() {
            if !self.selected[i] || polygon.is_empty() {
                continue;
            }
            let source = &self.cell.faces[i];
            let length = reduced(source.plane.normal.norm_squared().sqrt(), dec)?;
            let shift = source.plane.normal.scaled(&(thickness / &length));
            for v in polygon {
                original.entry(*v).or_insert_with(|| {
                    let point = model
                        .geom_mut()
                        .add_point(PointDef::Explicit(self.cell.vertices()[*v].clone()));
                    model.add_vertex(Vertex { point })
                });
                let point = model
                    .geom_mut()
                    .add_point(PointDef::Explicit(self.cell.vertices()[*v].plus(&shift)));
                let far = model.add_vertex(Vertex { point });
                displaced.insert((*v, i), far);
            }
            let mut base_path = Vec::new();
            let mut far_path = Vec::new();
            for j in 0..polygon.len() {
                let a = polygon[j];
                let b = polygon[(j + 1) % polygon.len()];
                let (lo, hi) = (a.min(b), a.max(b));
                let base = *base_edges
                    .entry((lo, hi))
                    .or_insert_with(|| line(model, original[&lo], original[&hi]));
                let far = line(model, displaced[&(lo, i)], displaced[&(hi, i)]);
                far_edges.insert((lo, hi, i), far);
                let sense = if a < b { Sense::Same } else { Sense::Reversed };
                base_path.push((base, sense));
                far_path.push((far, sense));
            }
            base_path.reverse();
            for (_, sense) in &mut base_path {
                *sense = reverse(*sense);
            }
            faces.push(face(
                model,
                Surface::Plane(source.plane.clone()),
                Sense::Reversed,
                vec![base_path],
                None,
                source.origins.clone(),
            ));
            let mut plane = source.plane.clone();
            plane.origin = plane.origin.plus(&shift);
            let names = source
                .origins
                .iter()
                .map(|of| Origin {
                    step: step.clone(),
                    role: Role::Offset {
                        of: Box::new(of.clone()),
                    },
                })
                .collect();
            faces.push(face(
                model,
                Surface::Plane(plane),
                Sense::Same,
                vec![far_path],
                None,
                names,
            ));
        }
        let mut arcs: BTreeMap<usize, Vec<ArcUse>> = BTreeMap::new();
        for ((a, b), uses) in &self.edges {
            let chosen = [
                self.selected[uses[0].0].then_some(uses[0]),
                self.selected[uses[1].0].then_some(uses[1]),
            ];
            match &chosen {
                [None, None] => (),
                [Some((f, sense)), None] | [None, Some((f, sense))] => {
                    let (start, end) = if *sense == Sense::Same {
                        (*a, *b)
                    } else {
                        (*b, *a)
                    };
                    for v in [start, end] {
                        rails
                            .entry((v, *f))
                            .or_insert_with(|| line(model, original[&v], displaced[&(v, *f)]));
                    }
                    let path = vec![
                        (base_edges[&(*a, *b)], *sense),
                        (rails[&(end, *f)], Sense::Same),
                        (far_edges[&(*a, *b, *f)], reverse(*sense)),
                        (rails[&(start, *f)], Sense::Reversed),
                    ];
                    let normal = self.cell.vertices()[end]
                        .minus(&self.cell.vertices()[start])
                        .cross(&self.cell.faces[*f].plane.normal);
                    let plane = Plane {
                        origin: self.cell.vertices()[start].clone(),
                        normal,
                        ref_dir: self.cell.faces[*f].plane.normal.clone(),
                    };
                    faces.push(face(
                        model,
                        Surface::Plane(plane),
                        Sense::Same,
                        vec![path],
                        None,
                        vec![Origin {
                            step: step.clone(),
                            role: Role::Rim {
                                edge: self.edge_name(uses)?,
                            },
                        }],
                    ));
                }
                [Some((f, sf)), Some((g, _))] => {
                    let nf = &self.cell.faces[*f].plane.normal;
                    let ng = &self.cell.faces[*g].plane.normal;
                    let axis = self.cell.vertices()[*b].minus(&self.cell.vertices()[*a]);
                    let arc_sense = match sign(&axis.dot(&nf.cross(ng)), dec)? {
                        Sign::Positive => Sense::Same,
                        Sign::Negative => Sense::Reversed,
                        Sign::Zero => {
                            return Err(bad("a selected edge has no regular angular normal sector"));
                        }
                    };
                    let caps = [*a, *b].map(|v| {
                        let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                            center: self.cell.vertices()[v].clone(),
                            normal: axis.clone(),
                            radius_sq: thickness.square(),
                            ref_dir: nf.clone(),
                        }));
                        model.add_edge(Edge {
                            curve,
                            sense: arc_sense,
                            start: Some(displaced[&(v, *f)]),
                            end: Some(displaced[&(v, *g)]),
                        })
                    });
                    let (at_a, at_b, path) = if *sf == Sense::Same {
                        (
                            Sense::Same,
                            Sense::Reversed,
                            vec![
                                (far_edges[&(*a, *b, *f)], Sense::Reversed),
                                (caps[0], Sense::Same),
                                (far_edges[&(*a, *b, *g)], Sense::Same),
                                (caps[1], Sense::Reversed),
                            ],
                        )
                    } else {
                        (
                            Sense::Reversed,
                            Sense::Same,
                            vec![
                                (far_edges[&(*a, *b, *f)], Sense::Same),
                                (caps[1], Sense::Same),
                                (far_edges[&(*a, *b, *g)], Sense::Reversed),
                                (caps[0], Sense::Reversed),
                            ],
                        )
                    };
                    for (v, edge, sense, direction) in
                        [(*a, caps[0], at_a, -1), (*b, caps[1], at_b, 1)]
                    {
                        let sense = reverse(sense);
                        let (start, end) = if sense == Sense::Same {
                            (*f, *g)
                        } else {
                            (*g, *f)
                        };
                        arcs.entry(v).or_default().push((
                            start,
                            end,
                            edge,
                            sense,
                            axis.scaled(&Real::from_i64(direction)),
                        ));
                    }
                    let cylinder = Cylinder::new(
                        self.cell.vertices()[*a].clone(),
                        axis,
                        thickness.clone(),
                        nf.clone(),
                        dec,
                    )?;
                    faces.push(face(
                        model,
                        Surface::Cylinder(cylinder),
                        Sense::Same,
                        vec![path],
                        None,
                        vec![Origin {
                            step: step.clone(),
                            role: Role::Round {
                                edge: self.edge_name(uses)?,
                            },
                        }],
                    ));
                }
            }
        }
        for (v, uses) in arcs {
            let incident = &self.vertices[&v];
            if !incident.iter().all(|i| self.selected[*i]) {
                let [(start, end, arc, sense, normal)] = uses.as_slice() else {
                    return Err(bad("a boundary end fan branches"));
                };
                let path = vec![
                    (*arc, *sense),
                    (rails[&(v, *end)], Sense::Reversed),
                    (rails[&(v, *start)], Sense::Same),
                ];
                let plane = Plane {
                    origin: self.cell.vertices()[v].clone(),
                    normal: normal.clone(),
                    ref_dir: self.cell.faces[*start].plane.normal.clone(),
                };
                faces.push(face(
                    model,
                    Surface::Plane(plane),
                    Sense::Same,
                    vec![path],
                    None,
                    self.boundary_names(v, step)?,
                ));
                continue;
            }
            let mut next = BTreeMap::new();
            for (start, end, edge, sense, _) in uses {
                if next.insert(start, (end, edge, sense)).is_some() {
                    return Err(bad("the selected spherical normal fan branches"));
                }
            }
            let start = *next
                .keys()
                .next()
                .ok_or_else(|| bad("a selected vertex has no normal fan"))?;
            let mut at = start;
            let mut path = Vec::new();
            loop {
                let (end, edge, sense) = next
                    .remove(&at)
                    .ok_or_else(|| bad("the selected spherical normal fan does not close"))?;
                path.push((edge, sense));
                at = end;
                if at == start {
                    break;
                }
            }
            if !next.is_empty() {
                return Err(bad("the selected spherical normal fan is disconnected"));
            }
            let source = &self.cell.faces[incident[0]].plane;
            let sphere = Sphere::new(
                self.cell.vertices()[v].clone(),
                thickness.clone(),
                source.normal.clone(),
                source.ref_dir.clone(),
                dec,
            )?;
            let names = incident
                .iter()
                .filter_map(|i| self.cell.faces[*i].origins.first().cloned())
                .collect();
            let vertex = VertexName::new(names).map_err(|e| bad(&e.to_string()))?;
            faces.push(face(
                model,
                Surface::Sphere(sphere),
                Sense::Same,
                vec![path],
                None,
                vec![Origin {
                    step: step.clone(),
                    role: Role::RoundCorner { vertex },
                }],
            ));
        }
        crate::sharp_thicken::assemble(model, faces, step)
    }
}
fn line(model: &mut Model, start: VertexId, end: VertexId) -> EdgeId {
    let a = model.geom().position(model.vertex(start).point);
    let b = model.geom().position(model.vertex(end).point);
    let dir = b.minus(&a);
    let curve = model
        .geom_mut()
        .add_curve(Curve::Line(Line { origin: a, dir }));
    model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: Some(start),
        end: Some(end),
    })
}
