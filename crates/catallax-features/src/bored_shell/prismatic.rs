//! Shared cap wires for strictly separated, constant-section cavity openings.

use std::collections::BTreeMap;

use catallax_numbers::MarkLog;
use catallax_shapes::{Curve, Plane, Point3};
use catallax_topology::{
    Body, BodyKind, Edge, EdgeId, FaceId, Loop, Lump, Sense, Shell, Vertex, VertexId,
};

use super::*;
use crate::draft::bored::Bore;

type Wire = Vec<(EdgeId, Sense)>;
type OffsetWalls = ([Wire; 2], Vec<FaceId>);

pub(super) struct Holes<'a> {
    pub bores: &'a [Bore],
    pub polygons: &'a [ConvexCell],
    pub caps: [usize; 2],
    pub marks: MarkLog,
}

fn value(plane: &Plane, point: &Point3) -> Real {
    point.minus(&plane.origin).dot(&plane.normal)
}

fn outside_margin(value: &Real, squared: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Positive
        && sign(&(&value.square() - squared), dec)? == Sign::Positive)
}

fn clear(
    cell: &ConvexCell,
    holes: &Holes<'_>,
    t: &Real,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let axis = &holes.bores[0].cylinder.axis;
    for (index, wall) in cell.faces.iter().enumerate() {
        if holes.caps.contains(&index) {
            continue;
        }
        let squared = t.square() * wall.plane.normal.norm_squared();
        for opening in holes.polygons {
            for point in opening.vertices() {
                if !outside_margin(&(-value(&wall.plane, point)), &squared, dec)? {
                    return Ok(false);
                }
            }
        }
    }
    for (index, opening) in holes.polygons.iter().enumerate() {
        for original in &opening.faces {
            if sign(&original.plane.normal.cross(axis).norm_squared(), dec)? == Sign::Zero {
                for cap in holes.caps {
                    if sign(&value(&original.plane, &cell.faces[cap].plane.origin), dec)?
                        == Sign::Positive
                    {
                        return Ok(false);
                    }
                }
            }
        }
        for bore in holes.bores {
            let mut separated = false;
            let radius = &bore.cylinder.radius_sq.sqrt() + t;
            for wall in &opening.faces {
                if sign(&wall.plane.normal.dot(axis), dec)? != Sign::Zero {
                    continue;
                }
                if outside_margin(
                    &value(&wall.plane, &bore.cylinder.origin),
                    &(&radius.square() * &wall.plane.normal.norm_squared()),
                    dec,
                )? {
                    separated = true;
                    break;
                }
            }
            if !separated {
                return Ok(false);
            }
        }
        for other in &holes.polygons[..index] {
            let mut separated = false;
            for (a, b) in [(opening, other), (other, opening)] {
                for wall in &a.faces {
                    if sign(&wall.plane.normal.dot(axis), dec)? != Sign::Zero {
                        continue;
                    }
                    let squared =
                        &t.square() * &(&Real::from_i64(4) * &wall.plane.normal.norm_squared());
                    let mut all = true;
                    for point in b.vertices() {
                        all &= outside_margin(&value(&wall.plane, point), &squared, dec)?;
                    }
                    if all {
                        separated = true;
                        break;
                    }
                }
                if separated {
                    break;
                }
            }
            if !separated {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn ordered_caps(
    planes: [Plane; 2],
    axis: &Vector3,
    dec: &mut Decider,
) -> Result<[Plane; 2], Refusal> {
    let mut planes = planes;
    if sign(&planes[0].normal.dot(axis), dec)? == Sign::Positive {
        planes.swap(0, 1);
    }
    Ok(planes)
}

fn row(point: &Point3, caps: &[Plane; 2], dec: &mut Decider) -> Result<usize, Refusal> {
    for (i, cap) in caps.iter().enumerate() {
        if sign(&value(cap, point), dec)? == Sign::Zero {
            return Ok(i);
        }
    }
    Err(bad("a prismatic opening edge is off both cap rows"))
}

fn path(model: &Model, wire: catallax_topology::LoopId) -> Vec<(EdgeId, Sense)> {
    let Loop::Coedges(coedges) = model.loop_(wire) else {
        unreachable!("prismatic cavity boundaries have coedges")
    };
    coedges
        .iter()
        .map(|c| {
            let c = model.coedge(*c);
            (c.edge, c.sense)
        })
        .collect()
}

fn cycle(model: &Model, uses: Vec<(EdgeId, Sense)>) -> Result<Vec<(EdgeId, Sense)>, Refusal> {
    if uses.len() == 1 && model.edge(uses[0].0).start.is_none() {
        return Ok(uses);
    }
    let mut next = BTreeMap::new();
    for (id, sense) in uses {
        let edge = model.edge(id);
        let (Some(a), Some(b)) = (edge.start, edge.end) else {
            return Err(bad("a polygonal cavity wire needs finite endpoints"));
        };
        let (a, b) = if sense == Sense::Same { (a, b) } else { (b, a) };
        if next.insert(a, (id, sense, b)).is_some() {
            return Err(bad("a cavity cap wire branches"));
        }
    }
    let Some(start) = next.keys().next().copied() else {
        return Err(bad("a cavity cap wire is empty"));
    };
    let mut at = start;
    let mut result = Vec::with_capacity(next.len());
    loop {
        let Some((edge, sense, end)) = next.remove(&at) else {
            return Err(bad("a cavity cap wire is open"));
        };
        result.push((edge, sense));
        at = end;
        if at == start {
            break;
        }
    }
    if !next.is_empty() {
        return Err(bad("a convex opening has several cap wires"));
    }
    Ok(result)
}

/// Retains vertical offset faces and moves their two end rows by common axial
/// translations. Circle centres assign rows without testing algebraic endpoints.
fn opening(
    model: &mut Model,
    tool: BodyId,
    old: &[Plane; 2],
    target: &[Plane; 2],
    axis: &Vector3,
    dec: &mut Decider,
) -> Result<OffsetWalls, Refusal> {
    let shifts: [Vector3; 2] = std::array::from_fn(|i| {
        axis.scaled(
            &(target[i]
                .origin
                .minus(&old[i].origin)
                .dot(&target[i].normal)
                / target[i].normal.dot(axis)),
        )
    });
    let mut walls = Vec::new();
    for face in model.faces_of_body(tool) {
        let keep = match model.geom().surface(model.face(face).surface) {
            Surface::Plane(plane) => sign(&plane.normal.dot(axis), dec)? == Sign::Zero,
            Surface::Cylinder(cylinder) => {
                sign(&cylinder.axis.cross(axis).norm_squared(), dec)? == Sign::Zero
            }
            _ => false,
        };
        if keep {
            walls.push(face);
        }
    }
    let mut uses: BTreeMap<EdgeId, Vec<Sense>> = BTreeMap::new();
    for face in &walls {
        for wire in &model.face(*face).loops {
            for (edge, sense) in path(model, *wire) {
                uses.entry(edge).or_default().push(sense);
            }
        }
    }
    let mut rows: BTreeMap<VertexId, usize> = BTreeMap::new();
    let mut edge_rows = BTreeMap::new();
    for edge in uses.keys() {
        let e = model.edge(*edge);
        if let Curve::Circle(circle) = model.geom().curve(e.curve) {
            let r = row(&circle.center, old, dec)?;
            edge_rows.insert(*edge, r);
            for vertex in [e.start, e.end].into_iter().flatten() {
                rows.insert(vertex, r);
            }
        }
    }
    let mut vertices = BTreeMap::new();
    for (vertex, r) in &rows {
        let point = model
            .geom()
            .position(model.vertex(*vertex).point)
            .plus(&shifts[*r]);
        let point = model
            .geom_mut()
            .add_point(catallax_shapes::PointDef::Explicit(point));
        vertices.insert(*vertex, model.add_vertex(Vertex { point }));
    }
    let mut edges = BTreeMap::new();
    let mut boundary: [Vec<(EdgeId, Sense)>; 2] = Default::default();
    for (edge, senses) in &uses {
        let e = model.edge(*edge).clone();
        let curve = match model.geom().curve(e.curve) {
            Curve::Circle(circle) => {
                let mut c = circle.clone();
                c.center = c.center.plus(&shifts[edge_rows[edge]]);
                Curve::Circle(c)
            }
            Curve::Line(line) => {
                let a = rows[&e.start.expect("rounded prism line start")];
                let b = rows[&e.end.expect("rounded prism line end")];
                Curve::Line(catallax_shapes::Line {
                    origin: line.origin.plus(&shifts[a]),
                    dir: line.dir.plus(&shifts[b].minus(&shifts[a])),
                })
            }
            _ => {
                return Err(bad(
                    "a rounded prismatic opening needs lines and circle sectors",
                ));
            }
        };
        let curve = model.geom_mut().add_curve(curve);
        let copied = model.add_edge(Edge {
            curve,
            sense: e.sense,
            start: e.start.map(|v| vertices[&v]),
            end: e.end.map(|v| vertices[&v]),
        });
        edges.insert(*edge, copied);
        match senses.as_slice() {
            [sense] => {
                let r = edge_rows
                    .get(edge)
                    .copied()
                    .unwrap_or_else(|| rows[&e.start.unwrap()]);
                boundary[r].push((copied, *sense));
            }
            [a, b] if a != b => {}
            _ => return Err(bad("rounded prism walls have inconsistent edge uses")),
        }
    }
    let mut result = Vec::new();
    for wall in walls {
        let f = model.face(wall).clone();
        let paths = f
            .loops
            .iter()
            .map(|wire| {
                path(model, *wire)
                    .into_iter()
                    .rev()
                    .map(|(edge, sense)| (edges[&edge], crate::extrude::reverse(sense)))
                    .collect()
            })
            .collect();
        result.push(crate::extrude::face(
            model,
            model.geom().surface(f.surface).clone(),
            crate::extrude::reverse(f.sense),
            paths,
            None,
            model.face_origins(wall).to_vec(),
        ));
    }
    let [a, b] = boundary;
    Ok(([cycle(model, a)?, cycle(model, b)?], result))
}

pub(super) fn build(
    model: &mut Model,
    cell: &ConvexCell,
    holes: Holes<'_>,
    t: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<BodyId>, Refusal> {
    if !clear(cell, &holes, t, dec)? {
        return Ok(None);
    }
    let axis = &holes.bores[0].cylinder.axis;
    let mut caps = holes.caps;
    if sign(&cell.faces[caps[0]].plane.normal.dot(axis), dec)? == Sign::Positive {
        caps.swap(0, 1);
    }
    let target = caps.map(|i| cell.faces[i].plane.clone());
    let mut cap_holes: [Vec<Wire>; 2] = Default::default();
    let mut inner = Vec::new();
    for bore in holes.bores {
        let mut ends = target.each_ref().map(|cap| axis_point(&bore.cylinder, cap));
        if sign(&ends[1].minus(&ends[0]).dot(&bore.cylinder.axis), dec)? == Sign::Negative {
            ends.swap(0, 1);
        }
        let tool = super::super::cylinder(
            model,
            &CylinderSkin {
                cylinder: bore.cylinder.clone(),
                ends,
                origins: [Vec::new(), Vec::new(), bore.origins.clone()],
            },
        );
        let (wires, walls) = opening(model, tool, &target, &target, axis, dec)?;
        for (holes, wire) in cap_holes.iter_mut().zip(wires) {
            holes.push(wire);
        }
        inner.extend(walls);
    }
    for polygon in holes.polygons {
        let old: Vec<_> = polygon
            .faces
            .iter()
            .filter_map(
                |face| match sign(&face.plane.normal.cross(axis).norm_squared(), dec) {
                    Ok(Sign::Zero) => Some(Ok(face.plane.clone())),
                    Ok(_) => None,
                    Err(error) => Some(Err(error)),
                },
            )
            .collect::<Result<_, _>>()?;
        let old = ordered_caps(
            old.try_into()
                .map_err(|_| bad("a polygonal prism needs two caps"))?,
            axis,
            dec,
        )?;
        let tool =
            crate::rounded_polyhedron::build(model, polygon, t, step, holes.marks.clone(), dec)?;
        model.map_face_origins(tool, |origin| match origin.role {
            Role::Round { .. } | Role::RoundCorner { .. } => origin.clone(),
            _ => Origin {
                step: step.clone(),
                role: Role::Inner {
                    of: Box::new(origin.clone()),
                },
            },
        });
        let (wires, walls) = opening(model, tool, &old, &target, axis, dec)?;
        for (holes, wire) in cap_holes.iter_mut().zip(wires) {
            holes.push(wire);
        }
        inner.extend(walls);
    }
    let outer = cell.build(model);
    let mut faces = model.faces_of_body(outer);
    let mut at = 0;
    for (index, cycle) in cell.polygons.iter().enumerate() {
        if cycle.is_empty() {
            continue;
        }
        if let Some(row) = caps.iter().position(|cap| *cap == index) {
            let original = faces[at];
            let f = model.face(original).clone();
            let mut paths: Vec<_> = f.loops.iter().map(|wire| path(model, *wire)).collect();
            paths.append(&mut cap_holes[row]);
            faces[at] = crate::extrude::face(
                model,
                model.geom().surface(f.surface).clone(),
                f.sense,
                paths,
                None,
                model.face_origins(original).to_vec(),
            );
        }
        at += 1;
    }
    faces.extend(inner);
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: holes.marks,
    });
    if !model.check(body).is_empty() {
        return Err(bad("the shared prismatic cavity has invalid incidence"));
    }
    Ok(Some(body))
}
