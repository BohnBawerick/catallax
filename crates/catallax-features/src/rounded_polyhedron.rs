//! Direct convex offset assembly. Convex normal fans cannot intersect each other.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, MarkLog, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, Origin, Role, StepId, VertexName};
use catallax_shapes::{Circle, Curve, Cylinder, Line, PointDef, Sphere, Surface};
use catallax_topology::{
    Body, BodyId, BodyKind, Edge, EdgeId, Lump, Model, Sense, Shell, Vertex, VertexId,
};

use crate::{
    extrude::{face, reverse},
    polyhedron::ConvexCell,
    profile::sign,
};

type ArcUse = (usize, usize, EdgeId, Sense);
fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.convex-normal-fan", detail)
}

pub(crate) fn build(
    model: &mut Model,
    cell: &ConvexCell,
    distance: &Real,
    step: &StepId,
    marks: MarkLog,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let mut displaced: BTreeMap<(usize, usize), VertexId> = BTreeMap::new();
    let mut lines = BTreeMap::new();
    let mut incidence: BTreeMap<(usize, usize), Vec<(usize, Sense)>> = BTreeMap::new();
    let mut faces = Vec::new();
    let mut vertex_faces: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (i, polygon) in cell.polygons.iter().enumerate() {
        if polygon.is_empty() {
            continue;
        }
        let source = &cell.faces[i];
        let direction = source
            .plane
            .normal
            .scaled(&(distance / &source.plane.normal.norm_squared().sqrt()));
        for v in polygon {
            let point = model
                .geom_mut()
                .add_point(PointDef::Explicit(cell.vertices()[*v].plus(&direction)));
            displaced.insert((*v, i), model.add_vertex(Vertex { point }));
            vertex_faces.entry(*v).or_default().push(i);
        }
        let mut path = Vec::with_capacity(polygon.len());
        for j in 0..polygon.len() {
            let a = polygon[j];
            let b = polygon[(j + 1) % polygon.len()];
            let (lo, hi) = (a.min(b), a.max(b));
            let curve = model.geom_mut().add_curve(Curve::Line(Line {
                origin: cell.vertices()[lo].plus(&direction),
                dir: cell.vertices()[hi].minus(&cell.vertices()[lo]),
            }));
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(displaced[&(lo, i)]),
                end: Some(displaced[&(hi, i)]),
            });
            lines.insert((lo, hi, i), edge);
            let sense = if a < b { Sense::Same } else { Sense::Reversed };
            incidence.entry((lo, hi)).or_default().push((i, sense));
            path.push((edge, sense));
        }
        let mut plane = source.plane.clone();
        plane.origin = plane.origin.plus(&direction);
        faces.push(face(
            model,
            Surface::Plane(plane),
            Sense::Same,
            vec![path],
            None,
            source.origins.clone(),
        ));
    }
    let mut arcs: BTreeMap<usize, Vec<ArcUse>> = BTreeMap::new();
    for ((a, b), uses) in incidence {
        let [(f, sf), (g, sg)] = uses.as_slice() else {
            return Err(bad("a convex edge must have two incident faces"));
        };
        if *sf == *sg {
            return Err(bad("input edge uses have inconsistent orientation"));
        }
        let nf = &cell.faces[*f].plane.normal;
        let ng = &cell.faces[*g].plane.normal;
        let axis = cell.vertices()[b].minus(&cell.vertices()[a]);
        let orientation = sign(&axis.dot(&nf.cross(ng)), dec)?;
        let arc_sense = match orientation {
            Sign::Positive => Sense::Same,
            Sign::Negative => Sense::Reversed,
            Sign::Zero => return Err(bad("a retained edge has no regular angular normal sector")),
        };
        let caps = [a, b].map(|vertex| {
            let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                center: cell.vertices()[vertex].clone(),
                normal: axis.clone(),
                radius_sq: distance.square(),
                ref_dir: nf.clone(),
            }));
            model.add_edge(Edge {
                curve,
                sense: arc_sense,
                start: Some(displaced[&(vertex, *f)]),
                end: Some(displaced[&(vertex, *g)]),
            })
        });
        let (at_a, at_b, path) = if *sf == Sense::Same {
            (
                Sense::Same,
                Sense::Reversed,
                vec![
                    (lines[&(a, b, *f)], Sense::Reversed),
                    (caps[0], Sense::Same),
                    (lines[&(a, b, *g)], Sense::Same),
                    (caps[1], Sense::Reversed),
                ],
            )
        } else {
            (
                Sense::Reversed,
                Sense::Same,
                vec![
                    (lines[&(a, b, *f)], Sense::Same),
                    (caps[1], Sense::Same),
                    (lines[&(a, b, *g)], Sense::Reversed),
                    (caps[0], Sense::Reversed),
                ],
            )
        };
        for (vertex, edge, sense) in [(a, caps[0], at_a), (b, caps[1], at_b)] {
            let sense = reverse(sense);
            let (start, end) = if sense == Sense::Same {
                (*f, *g)
            } else {
                (*g, *f)
            };
            arcs.entry(vertex)
                .or_default()
                .push((start, end, edge, sense));
        }
        let origin = |i: usize| {
            cell.faces[i]
                .origins
                .first()
                .cloned()
                .ok_or_else(|| bad("an offset edge requires named input faces"))
        };
        let name = EdgeName::new(origin(*f)?, origin(*g)?).map_err(|e| bad(&e.to_string()))?;
        let cylinder = Cylinder::new(
            cell.vertices()[a].clone(),
            axis,
            distance.clone(),
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
                role: Role::Round { edge: name },
            }],
        ));
    }
    for (vertex, uses) in arcs {
        let mut next = BTreeMap::new();
        for (start, end, edge, sense) in uses {
            if next.insert(start, (end, edge, sense)).is_some() {
                return Err(bad("the spherical normal fan branches"));
            }
        }
        let Some(start) = next.keys().next().copied() else {
            return Err(bad("an offset vertex has no spherical fan"));
        };
        let mut at = start;
        let mut path = Vec::with_capacity(next.len());
        loop {
            let Some((end, edge, sense)) = next.remove(&at) else {
                return Err(bad("the spherical normal fan does not close"));
            };
            path.push((edge, sense));
            at = end;
            if at == start {
                break;
            }
        }
        if !next.is_empty() {
            return Err(bad("an input vertex has disconnected normal fans"));
        }
        let incident = &vertex_faces[&vertex];
        let source = &cell.faces[incident[0]].plane;
        let mut origins: Vec<_> = incident
            .iter()
            .filter_map(|i| cell.faces[*i].origins.first().cloned())
            .collect();
        origins.sort();
        origins.dedup();
        let name = VertexName::new(origins).map_err(|e| bad(&e.to_string()))?;
        let sphere = Sphere::new(
            cell.vertices()[vertex].clone(),
            distance.clone(),
            source.normal.clone(),
            source.ref_dir.clone(),
            dec,
        )?;
        faces.push(face(
            model,
            Surface::Sphere(sphere),
            Sense::Same,
            vec![path],
            None,
            vec![Origin {
                step: step.clone(),
                role: Role::RoundCorner { vertex: name },
            }],
        ));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    Ok(model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks,
    }))
}

/// Exact distance tube of one closed convex polygon, including its semicircular
/// edge rolls and spherical vertex lunes. The source solid supplies the names.
pub(crate) fn plate(
    model: &mut Model,
    cell: &ConvexCell,
    selected: usize,
    distance: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let polygon = &cell.polygons[selected];
    if polygon.len() < 3 {
        return Ok(catallax_topology::empty(model));
    }
    let source = &cell.faces[selected];
    let direction = source
        .plane
        .normal
        .scaled(&(distance / &source.plane.normal.norm_squared().sqrt()));
    let mut vertices = Vec::with_capacity(polygon.len());
    for v in polygon {
        let pair: [VertexId; 2] = std::array::from_fn(|side| {
            let sign = Real::from_i64(if side == 0 { -1 } else { 1 });
            let point = model.geom_mut().add_point(PointDef::Explicit(
                cell.vertices()[*v].plus(&direction.scaled(&sign)),
            ));
            model.add_vertex(Vertex { point })
        });
        vertices.push(pair);
    }
    let mut lines = Vec::with_capacity(polygon.len());
    for i in 0..polygon.len() {
        let j = (i + 1) % polygon.len();
        let pair: [EdgeId; 2] = std::array::from_fn(|side| {
            let sign = Real::from_i64(if side == 0 { -1 } else { 1 });
            let curve = model.geom_mut().add_curve(Curve::Line(Line {
                origin: cell.vertices()[polygon[i]].plus(&direction.scaled(&sign)),
                dir: cell.vertices()[polygon[j]].minus(&cell.vertices()[polygon[i]]),
            }));
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(vertices[i][side]),
                end: Some(vertices[j][side]),
            })
        });
        lines.push(pair);
    }
    let mut faces = Vec::new();
    for side in 0..2 {
        let signum = Real::from_i64(if side == 0 { -1 } else { 1 });
        let mut plane = source.plane.clone();
        plane.origin = plane.origin.plus(&direction.scaled(&signum));
        plane.normal = plane.normal.scaled(&signum);
        let path = if side == 1 {
            lines.iter().map(|pair| (pair[side], Sense::Same)).collect()
        } else {
            lines
                .iter()
                .rev()
                .map(|pair| (pair[side], Sense::Reversed))
                .collect()
        };
        let origins = source
            .origins
            .iter()
            .map(|of| Origin {
                step: step.clone(),
                role: Role::Inner {
                    of: Box::new(of.clone()),
                },
            })
            .collect();
        faces.push(face(
            model,
            Surface::Plane(plane),
            Sense::Same,
            vec![path],
            None,
            origins,
        ));
    }
    let mut vertex_arcs = vec![[None; 2]; polygon.len()];
    for i in 0..polygon.len() {
        let j = (i + 1) % polygon.len();
        let a = polygon[i];
        let b = polygon[j];
        let axis = cell.vertices()[b].minus(&cell.vertices()[a]);
        let caps = [(i, a), (j, b)].map(|(local, global)| {
            let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                center: cell.vertices()[global].clone(),
                normal: axis.clone(),
                radius_sq: distance.square(),
                ref_dir: source.plane.normal.clone(),
            }));
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(vertices[local][1]),
                end: Some(vertices[local][0]),
            })
        });
        vertex_arcs[i][0] = Some((caps[0], Sense::Reversed));
        vertex_arcs[j][1] = Some((caps[1], Sense::Same));
        let neighbour = cell
            .polygons
            .iter()
            .enumerate()
            .find_map(|(f, p)| (f != selected && p.contains(&a) && p.contains(&b)).then_some(f))
            .ok_or_else(|| bad("a kept polygon edge has no input neighbour"))?;
        let a_name = source
            .origins
            .first()
            .cloned()
            .ok_or_else(|| bad("a kept polygon requires an origin"))?;
        let b_name = cell.faces[neighbour]
            .origins
            .first()
            .cloned()
            .ok_or_else(|| bad("a kept polygon neighbour requires an origin"))?;
        let edge = EdgeName::new(a_name, b_name).map_err(|e| bad(&e.to_string()))?;
        let cylinder = Cylinder::new(
            cell.vertices()[a].clone(),
            axis,
            distance.clone(),
            source.plane.normal.clone(),
            dec,
        )?;
        let path = vec![
            (lines[i][1], Sense::Reversed),
            (caps[0], Sense::Same),
            (lines[i][0], Sense::Same),
            (caps[1], Sense::Reversed),
        ];
        faces.push(face(
            model,
            Surface::Cylinder(cylinder),
            Sense::Same,
            vec![path],
            None,
            vec![Origin {
                step: step.clone(),
                role: Role::Round { edge },
            }],
        ));
    }
    for (i, uses) in vertex_arcs.into_iter().enumerate() {
        let [Some(outgoing), Some(incoming)] = uses else {
            return Err(bad(
                "a convex polygon vertex requires two semicircular ends",
            ));
        };
        let path = vec![incoming, outgoing];
        let vertex = polygon[i];
        let mut origins: Vec<_> = cell
            .polygons
            .iter()
            .enumerate()
            .filter(|(_, p)| p.contains(&vertex))
            .filter_map(|(f, _)| cell.faces[f].origins.first().cloned())
            .collect();
        origins.sort();
        origins.dedup();
        let vertex_name = VertexName::new(origins).map_err(|e| bad(&e.to_string()))?;
        let sphere = Sphere::new(
            cell.vertices()[vertex].clone(),
            distance.clone(),
            source.plane.normal.clone(),
            source.plane.ref_dir.clone(),
            dec,
        )?;
        faces.push(face(
            model,
            Surface::Sphere(sphere),
            Sense::Same,
            vec![path],
            None,
            vec![Origin {
                step: step.clone(),
                role: Role::RoundCorner {
                    vertex: vertex_name,
                },
            }],
        ));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    Ok(model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    }))
}
