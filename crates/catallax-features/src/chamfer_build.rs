//! Rail-first elementary chamfer construction over an immutable input skin.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign, StrictOrder};
use catallax_recipe::{EdgeSelection, Origin, Role, StepId};
use catallax_shapes::{Curve, Line, Point3, PointDef, Surface};
use catallax_topology::{BodyId, Edge, EdgeId, FaceId, Loop, Model, Sense, Vertex, VertexId};

use crate::{
    chamfer::{ChamferDefinition, arc::BoundedArc},
    edge_blend::{self, Replacement},
    extrude::reverse,
    profile::sign,
};

type Uses = BTreeMap<EdgeId, Vec<(FaceId, Sense)>>;

fn bad(detail: &str) -> Refusal {
    Refusal::new("chamfer.unsupported-topology", detail)
}
fn vertex(model: &mut Model, p: Point3) -> VertexId {
    let point = model.geom_mut().add_point(PointDef::Explicit(p));
    model.add_vertex(Vertex { point })
}
fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}
fn uses(model: &Model, body: BodyId) -> Uses {
    let mut uses: Uses = BTreeMap::new();
    for face in model.faces_of_body(body) {
        for lp in &model.face(face).loops {
            if let Loop::Coedges(cs) = model.loop_(*lp) {
                for c in cs {
                    let c = model.coedge(*c);
                    uses.entry(c.edge).or_default().push((face, c.sense));
                }
            }
        }
    }
    uses
}
fn matching(model: &Model, face: FaceId, pattern: &catallax_recipe::OriginPattern) -> bool {
    model.face_origins(face).iter().any(|o| pattern.matches(o))
}
fn resolve(
    model: &Model,
    selection: &EdgeSelection,
    uses: &Uses,
    distances: &[Real; 2],
    unequal: bool,
    dec: &mut Decider,
) -> Result<BTreeMap<EdgeId, [Real; 2]>, Refusal> {
    let mut selected = BTreeMap::new();
    let EdgeSelection::Refs(refs) = selection else {
        if unequal {
            return Err(Refusal::new(
                "chamfer.ambiguous-order",
                "two setbacks require ordered edge references",
            ));
        }
        return Ok(uses.keys().map(|e| (*e, distances.clone())).collect());
    };
    for reference in refs {
        if reference.near.is_some() {
            return Err(bad("nearest-edge references need an exact distance query"));
        }
        let mut found = false;
        for (edge, incident) in uses {
            let [(a, _), (b, _)] = incident.as_slice() else {
                continue;
            };
            let forward = matching(model, *a, &reference.between[0])
                && matching(model, *b, &reference.between[1]);
            let backward = matching(model, *b, &reference.between[0])
                && matching(model, *a, &reference.between[1]);
            if !forward && !backward {
                continue;
            }
            let e = model.edge(*edge);
            if !reference.touching.iter().all(|pattern| {
                uses.iter().any(|(other, faces)| {
                    let other = model.edge(*other);
                    [e.start, e.end]
                        .into_iter()
                        .flatten()
                        .any(|v| other.start == Some(v) || other.end == Some(v))
                        && faces.iter().any(|(f, _)| matching(model, *f, pattern))
                })
            }) {
                continue;
            }
            if unequal && forward && backward {
                return Err(Refusal::new(
                    "chamfer.ambiguous-order",
                    "an edge matches its two-distance reference both ways",
                ));
            }
            let ordered = if forward {
                distances.clone()
            } else {
                [distances[1].clone(), distances[0].clone()]
            };
            if let Some(previous) = selected.get(edge) {
                let previous: &[Real; 2] = previous;
                for (a, b) in previous.iter().zip(&ordered) {
                    if sign(&(a - b), dec)? != Sign::Zero {
                        return Err(Refusal::new(
                            "chamfer.conflicting-distances",
                            "overlapping references assign different setbacks to an input face",
                        ));
                    }
                }
            }
            selected.insert(*edge, ordered);
            found = true;
        }
        if !found {
            return Err(Refusal::new(
                "chamfer.lost-edge",
                "an edge reference has no input boundary",
            ));
        }
    }
    Ok(selected)
}

struct Run {
    edge: EdgeId,
    faces: [FaceId; 2],
    senses: [Sense; 2],
    definition: ChamferDefinition,
    ends: BTreeMap<(VertexId, FaceId), VertexId>,
    caps: BTreeMap<VertexId, EdgeId>,
    end_planes: BTreeMap<VertexId, catallax_shapes::Plane>,
}
fn outward(
    model: &Model,
    face: FaceId,
    at: &Point3,
    dec: &mut Decider,
) -> Result<catallax_shapes::Vector3, Refusal> {
    let f = model.face(face);
    let surface = model.geom().surface(f.surface);
    let n = surface.normal_at(at, dec)?;
    if sign(&n.norm_squared(), dec)? == Sign::Zero {
        return Err(bad("singular support has no normal"));
    }
    Ok(n.scaled(&Real::from_i64(if f.sense == Sense::Same { 1 } else { -1 })))
}
fn convex(
    model: &Model,
    edge: EdgeId,
    incident: &[(FaceId, Sense)],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let [(a, sense), (b, _)] = incident else {
        return Err(bad("a chamfer edge must have two incident faces"));
    };
    let e = model.edge(edge);
    let orientation = if *sense == e.sense { 1 } else { -1 };
    match model.geom().curve(e.curve) {
        Curve::Circle(c) => {
            let plane = [*a, *b]
                .into_iter()
                .find(|f| {
                    matches!(
                        model.geom().surface(model.face(*f).surface),
                        Surface::Plane(_)
                    )
                })
                .ok_or_else(|| bad("circular elementary chamfers need one plane"))?;
            let coedge = incident.iter().find(|(f, _)| *f == plane).unwrap().1;
            let orientation = if coedge == e.sense { 1 } else { -1 };
            let n = outward(model, plane, &c.center, dec)?;
            let cylinder = [*a, *b]
                .into_iter()
                .find(|f| {
                    matches!(
                        model.geom().surface(model.face(*f).surface),
                        Surface::Cylinder(_)
                    )
                })
                .ok_or_else(|| bad("circular bevels need a cylindrical incident face"))?;
            let radial_orientation = if model.face(cylinder).sense == Sense::Same {
                1
            } else {
                -1
            };
            Ok(sign(
                &(c.normal.dot(&n) * Real::from_i64(orientation * radial_orientation)),
                dec,
            )? == Sign::Positive)
        }
        Curve::Line(l) => {
            let na = outward(model, *a, &l.origin, dec)?;
            let nb = outward(model, *b, &l.origin, dec)?;
            let value = nb.dot(&na.cross(&l.dir)) * Real::from_i64(orientation);
            match sign(&value, dec)? {
                Sign::Negative => Ok(true),
                Sign::Positive => Ok(false),
                Sign::Zero => Err(bad("tangent faces do not define a regular sharp edge")),
            }
        }
        _ => Err(bad("elementary chamfers need a line or circle")),
    }
}
fn line_plane(
    line: &Line,
    plane: &catallax_shapes::Plane,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let divisor = line.dir.dot(&plane.normal);
    if sign(&divisor, dec)? == Sign::Zero {
        return Err(bad("end face is not transversal to the contact rail"));
    }
    Ok(line.origin.plus(
        &line
            .dir
            .scaled(&(plane.origin.minus(&line.origin).dot(&plane.normal) / divisor)),
    ))
}
fn line_line(a: &Line, b: &Line, dec: &mut Decider) -> Result<Point3, Refusal> {
    let cross = a.dir.cross(&b.dir);
    let divisor = cross.norm_squared();
    if sign(&divisor, dec)? == Sign::Zero {
        return Err(bad("parallel contact rails need chain continuation"));
    }
    let t = b.origin.minus(&a.origin).cross(&b.dir).dot(&cross) / divisor;
    let p = a.origin.plus(&a.dir.scaled(&t));
    if sign(&p.minus(&b.origin).cross(&b.dir).norm_squared(), dec)? != Sign::Zero {
        return Err(bad("contact rails do not meet"));
    }
    Ok(p)
}
fn rail(run: &Run, face: FaceId) -> Result<&Line, Refusal> {
    let index = run
        .faces
        .iter()
        .position(|f| *f == face)
        .ok_or_else(|| bad("rail has no input support"))?;
    let Curve::Line(line) = &run.definition.contacts[index] else {
        return Err(bad("vertex trimming needs straight rails"));
    };
    Ok(line)
}
fn plane(model: &Model, face: FaceId) -> Result<&catallax_shapes::Plane, Refusal> {
    let Surface::Plane(p) = model.geom().surface(model.face(face).surface) else {
        return Err(bad(
            "transversal line-strip ends currently require planar neighbours",
        ));
    };
    Ok(p)
}

fn end_contact(
    model: &Model,
    run: &Run,
    face: FaceId,
    old: VertexId,
    end_plane: &catallax_shapes::Plane,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let index = run.faces.iter().position(|f| *f == face).unwrap();
    match (&run.definition.edge, &run.definition.contacts[index]) {
        (Curve::Line(_), Curve::Line(line)) => line_plane(line, end_plane, dec),
        (Curve::Circle(circle), Curve::Circle(contact)) => {
            let input = model.edge(run.edge);
            let mut ends = [input.start, input.end].map(|v| {
                v.map(|v| model.geom().position(model.vertex(v).point))
                    .ok_or_else(|| bad("a finite circular run needs both endpoints"))
            });
            if input.sense == Sense::Reversed {
                ends.swap(0, 1);
            }
            let [a, b] = ends;
            let arc = BoundedArc::new(circle.clone(), [a?, b?], dec)?;
            let end = usize::from(if input.sense == Sense::Same {
                input.end == Some(old)
            } else {
                input.start == Some(old)
            });
            arc.contact_at(contact, end, end_plane, dec)
        }
        _ => Err(bad("a contact rail has no elementary end trim")),
    }
}

fn conic_connector(
    model: &mut Model,
    first: &Surface,
    second: &Surface,
    from: VertexId,
    to: VertexId,
    dec: &mut Decider,
) -> Result<EdgeId, Refusal> {
    use catallax_intersect::{Component, SurfacePair};
    let a = model.geom_mut().add_surface(first.clone());
    let b = model.geom_mut().add_surface(second.clone());
    let SurfacePair::Meets(components) =
        catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value
    else {
        return Err(bad("a conical end has no exact plane intersection"));
    };
    let points = [from, to].map(|v| model.vertex(v).point);
    for component in components {
        let Component::Curve {
            curve, singular, ..
        } = component
        else {
            continue;
        };
        if !singular.is_empty() || catallax_intersect::curve_is_closed(model.geom(), curve, dec)? {
            continue;
        }
        let mut contains = true;
        for point in points {
            match catallax_intersect::on_curve(model.geom(), curve, point, dec) {
                Outcome::ProvenCoincident(_) => (),
                Outcome::Proven(_) => contains = false,
                Outcome::NotExact(n) => return Err(bad(&n.detail)),
                Outcome::Refused(r) => return Err(r),
            }
        }
        if !contains {
            continue;
        }
        let sense = match catallax_intersect::order_on_curve(
            model.geom(),
            curve,
            points[0],
            points[1],
            dec,
        ) {
            Outcome::Proven(StrictOrder::Less) => Sense::Same,
            Outcome::Proven(StrictOrder::Greater) => Sense::Reversed,
            Outcome::Refused(r) => return Err(r),
            _ => return Err(bad("a conical end boundary has unresolved ordering")),
        };
        return Ok(model.add_edge(Edge {
            curve,
            sense,
            start: Some(from),
            end: Some(to),
        }));
    }
    Err(bad(
        "the two contacts do not share one regular conical end branch",
    ))
}

fn connector(
    model: &mut Model,
    first: &Surface,
    second: &Surface,
    from: VertexId,
    to: VertexId,
    dec: &mut Decider,
) -> Result<EdgeId, Refusal> {
    let p = model.geom().position(model.vertex(from).point);
    let q = model.geom().position(model.vertex(to).point);
    if same(&p, &q, dec)? {
        return Err(bad("a collapsed end needs corner continuation"));
    }
    for surface in [first, second] {
        match surface {
            Surface::Plane(plane) => {
                for point in [&p, &q] {
                    if sign(&point.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero {
                        return Err(bad("connector endpoints do not lie on both supports"));
                    }
                }
            }
            Surface::Cone(cone) => {
                let a = p.minus(&cone.apex);
                let b = q.minus(&cone.apex);
                if sign(&a.cross(&b).norm_squared(), dec)? != Sign::Zero {
                    return conic_connector(model, first, second, from, to, dec);
                }
                if sign(&a.dot(&cone.axis), dec)? != Sign::Positive
                    || sign(&b.dot(&cone.axis), dec)? != Sign::Positive
                {
                    return Err(bad(
                        "a conical strip end needs one regular shared generator",
                    ));
                }
                let polynomial = surface.implicit().expect("cone is algebraic");
                for point in [&p, &q] {
                    if sign(
                        &polynomial.eval_real(&[point.x.clone(), point.y.clone(), point.z.clone()]),
                        dec,
                    )? != Sign::Zero
                    {
                        return Err(bad("generator misses its conical support"));
                    }
                }
            }
            _ => {
                return Err(bad(
                    "elementary connectors need planes or conical generators",
                ));
            }
        }
    }
    let curve = model.geom_mut().add_curve(Curve::Line(Line {
        origin: p.clone(),
        dir: q.minus(&p),
    }));
    Ok(model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: Some(from),
        end: Some(to),
    }))
}
fn region_vertices(model: &Model, run: &Run) -> Vec<Point3> {
    let input = model.edge(run.edge);
    [input.start, input.end]
        .into_iter()
        .flatten()
        .chain(run.ends.values().copied())
        .map(|v| model.geom().position(model.vertex(v).point))
        .collect()
}

fn guard(
    model: &mut Model,
    body: BodyId,
    run: &Run,
    uses: &Uses,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let edge = model.edge(run.edge);
    let ends = [edge.start, edge.end];
    let neighbours: BTreeSet<_> = uses
        .iter()
        .filter(|(e, _)| {
            let e = model.edge(**e);
            ends.into_iter()
                .flatten()
                .any(|v| e.start == Some(v) || e.end == Some(v))
        })
        .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
        .collect();
    for face in model.faces_of_body(body) {
        if run.faces.contains(&face) || neighbours.contains(&face) {
            continue;
        }
        let other = model.geom().surface(model.face(face).surface);
        if !run
            .definition
            .avoids_support(other, &region_vertices(model, run), dec)?
        {
            return Err(Refusal::new(
                "chamfer.unsupported-overflow",
                "the affected region may reach a non-end support; neighbour continuation needs proof",
            ));
        }
    }
    Ok(())
}

pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    selection: &EdgeSelection,
    distances: &[Real; 2],
    unequal: bool,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<ChamferDefinition>), Refusal> {
    for distance in distances {
        if sign(distance, dec)? != Sign::Positive {
            return Err(Refusal::new(
                "chamfer.non-positive-distance",
                "chamfer distances must be positive",
            ));
        }
    }
    let uses = uses(model, body);
    let selected = resolve(model, selection, &uses, distances, unequal, dec)?;
    if selected.is_empty() {
        return Ok((body, Vec::new()));
    }
    let mut runs = Vec::with_capacity(selected.len());
    for (edge, distances) in selected {
        let incident = &uses[&edge];
        let [(a, sa), (b, sb)] = incident.as_slice() else {
            return Err(bad("a selected edge must have two incident faces"));
        };
        let name = model
            .edge_name(body, edge)
            .ok_or_else(|| bad("a selected edge needs an unambiguous origin-pair name"))?;
        let faces = [*a, *b];
        let definition = ChamferDefinition::new(
            Origin {
                step: step.clone(),
                role: Role::Bevel { edge: name },
            },
            faces.map(|f| model.geom().surface(model.face(f).surface).clone()),
            faces.map(|f| model.face(f).sense == Sense::Same),
            model.geom().curve(model.edge(edge).curve).clone(),
            distances,
            convex(model, edge, incident, dec)?,
            dec,
        )?;
        runs.push(Run {
            edge,
            faces,
            senses: [*sa, *sb],
            definition,
            ends: BTreeMap::new(),
            caps: BTreeMap::new(),
            end_planes: BTreeMap::new(),
        });
    }
    let mut edit = Replacement::default();
    let mut end_changes = BTreeMap::new();
    for old in model.vertices_of_body(body) {
        let active: Vec<_> = runs
            .iter()
            .enumerate()
            .filter_map(|(i, run)| {
                let e = model.edge(run.edge);
                (e.start == Some(old) || e.end == Some(old)).then_some(i)
            })
            .collect();
        if active.is_empty() {
            continue;
        }
        let incident: BTreeSet<_> = uses
            .iter()
            .filter(|(edge, _)| {
                let e = model.edge(**edge);
                e.start == Some(old) || e.end == Some(old)
            })
            .flat_map(|(_, incident)| incident.iter().map(|(f, _)| *f))
            .collect();
        if incident.len() != 3 {
            return Err(bad("line-strip end needs a regular three-face vertex"));
        }
        match active.as_slice() {
            [i] => {
                let faces = runs[*i].faces;
                let third = *incident.iter().find(|f| !faces.contains(f)).unwrap();
                let end_plane = plane(model, third)?.clone();
                runs[*i].end_planes.insert(old, end_plane.clone());
                let a = end_contact(model, &runs[*i], faces[0], old, &end_plane, dec)?;
                let b = end_contact(model, &runs[*i], faces[1], old, plane(model, third)?, dec)?;
                let [a, b] = [vertex(model, a), vertex(model, b)];
                for (face, v) in faces.into_iter().zip([a, b]) {
                    runs[*i].ends.insert((old, face), v);
                }
                let surface = runs[*i].definition.surface.clone();
                let third_surface = model.geom().surface(model.face(third).surface).clone();
                let cap = connector(model, &surface, &third_surface, a, b, dec)?;
                runs[*i].caps.insert(old, cap);
                edit.gaps.insert((third, old), cap);
                for (edge, adjacent) in &uses {
                    if *edge == runs[*i].edge {
                        continue;
                    }
                    let e = model.edge(*edge);
                    if e.start != Some(old) && e.end != Some(old) {
                        continue;
                    }
                    let v = if adjacent.iter().any(|(f, _)| *f == faces[0]) {
                        a
                    } else {
                        b
                    };
                    end_changes.insert((*edge, old), v);
                }
            }
            [i, j] => {
                let common = *runs[*i]
                    .faces
                    .iter()
                    .find(|f| runs[*j].faces.contains(f))
                    .ok_or_else(|| bad("selected strips do not share a support at their end"))?;
                let a = *runs[*i].faces.iter().find(|f| **f != common).unwrap();
                let b = *runs[*j].faces.iter().find(|f| **f != common).unwrap();
                let p = line_line(rail(&runs[*i], common)?, rail(&runs[*j], common)?, dec)?;
                let qa = line_plane(rail(&runs[*i], a)?, plane(model, b)?, dec)?;
                let qb = line_plane(rail(&runs[*j], b)?, plane(model, a)?, dec)?;
                if !same(&qa, &qb, dec)? {
                    return Err(bad(
                        "unequal endpoint contacts require a corner continuation patch",
                    ));
                }
                let p = vertex(model, p);
                let q = vertex(model, qa);
                for index in [*i, *j] {
                    for face in runs[index].faces {
                        runs[index]
                            .ends
                            .insert((old, face), if face == common { p } else { q });
                    }
                }
                let cap = connector(
                    model,
                    &runs[*i].definition.surface,
                    &runs[*j].definition.surface,
                    p,
                    q,
                    dec,
                )?;
                for index in [*i, *j] {
                    runs[index].caps.insert(old, cap);
                }
                for edge in uses.keys() {
                    if *edge == runs[*i].edge || *edge == runs[*j].edge {
                        continue;
                    }
                    let e = model.edge(*edge);
                    if e.start == Some(old) || e.end == Some(old) {
                        end_changes.insert((*edge, old), q);
                    }
                }
            }
            _ => return Err(bad("three selected strips at a vertex need a corner patch")),
        }
    }
    for run in &runs {
        let input = model.edge(run.edge);
        if let (Curve::Circle(circle), Some(start), Some(end)) =
            (&run.definition.edge, input.start, input.end)
        {
            let mut ends = [start, end];
            if input.sense == Sense::Reversed {
                ends.swap(0, 1);
            }
            let arc = BoundedArc::new(
                circle.clone(),
                ends.map(|v| model.geom().position(model.vertex(v).point)),
                dec,
            )?;
            let planes = ends.map(|v| {
                run.end_planes
                    .get(&v)
                    .cloned()
                    .ok_or_else(|| bad("a circular corner needs continuation"))
            });
            let [a, b] = planes;
            arc.trim(&run.definition, &[a?, b?], dec)?;
        }
        guard(model, body, run, &uses, dec)?;
    }
    for (i, a) in runs.iter().enumerate() {
        for b in &runs[i + 1..] {
            let ae = model.edge(a.edge);
            let be = model.edge(b.edge);
            if [ae.start, ae.end]
                .into_iter()
                .flatten()
                .any(|v| be.start == Some(v) || be.end == Some(v))
            {
                continue;
            }
            if !a.definition.avoids_region(
                &b.definition,
                &region_vertices(model, a),
                &region_vertices(model, b),
                dec,
            )? {
                return Err(Refusal::new(
                    "chamfer.unsupported-overflow",
                    "nonincident affected regions need a simultaneous overlap proof",
                ));
            }
        }
    }
    for (old, incident) in &uses {
        let e = model.edge(*old).clone();
        let start = e
            .start
            .map(|v| end_changes.get(&(*old, v)).copied().unwrap_or(v));
        let end = e
            .end
            .map(|v| end_changes.get(&(*old, v)).copied().unwrap_or(v));
        if start == e.start && end == e.end {
            continue;
        }
        if let (Some(old_start), Some(old_end), Some(start), Some(end)) =
            (e.start, e.end, start, end)
        {
            let position = |v| model.geom().position(model.vertex(v).point);
            let old_direction = position(old_end).minus(&position(old_start));
            if sign(
                &position(end).minus(&position(start)).dot(&old_direction),
                dec,
            )? != Sign::Positive
            {
                return Err(Refusal::new(
                    "chamfer.unsupported-overflow",
                    "end trimming consumes or reverses a neighbouring edge",
                ));
            }
        }
        let edge = model.add_edge(Edge { start, end, ..e });
        for (face, _) in incident {
            edit.edges.insert((*face, *old), edge);
        }
    }
    for run in &runs {
        let input = model.edge(run.edge).clone();
        let mut rails = Vec::with_capacity(2);
        for (i, face) in run.faces.iter().enumerate() {
            for old in [input.start, input.end].into_iter().flatten() {
                let v = run.ends[&(old, *face)];
                let p = model.geom().position(model.vertex(v).point);
                if !catallax_boolean::face_contains(model, *face, &p, dec)? {
                    return Err(Refusal::new(
                        "chamfer.unsupported-overflow",
                        "a contact foot leaves the input face and needs neighbour continuation",
                    ));
                }
            }
            if let (Some(start), Some(end)) = (input.start, input.end) {
                let position = |v| model.geom().position(model.vertex(v).point);
                let direction = position(end).minus(&position(start));
                let a = run.ends[&(start, *face)];
                let b = run.ends[&(end, *face)];
                if sign(&position(b).minus(&position(a)).dot(&direction), dec)? != Sign::Positive {
                    return Err(Refusal::new(
                        "chamfer.unsupported-overflow",
                        "contact runs collapse or reverse and need continuation",
                    ));
                }
            }
            let curve = model
                .geom_mut()
                .add_curve(run.definition.contacts[i].clone());
            let edge = model.add_edge(Edge {
                curve,
                sense: input.sense,
                start: input.start.map(|v| run.ends[&(v, *face)]),
                end: input.end.map(|v| run.ends[&(v, *face)]),
            });
            edit.edges.insert((*face, run.edge), edge);
            rails.push((edge, reverse(run.senses[i])));
        }
        let paths = if input.start.is_none() && input.end.is_none() {
            vec![vec![rails[0]], vec![rails[1]]]
        } else {
            let mut path = vec![rails[0]];
            let end = edge_blend::ends(model, rails[0].0, rails[0].1).1.unwrap();
            let old = if rails[0].1 == Sense::Same {
                input.end.unwrap()
            } else {
                input.start.unwrap()
            };
            let cap = run.caps[&old];
            let sense = if model.edge(cap).start == Some(end) {
                Sense::Same
            } else {
                Sense::Reversed
            };
            path.push((cap, sense));
            path.push(rails[1]);
            let end = edge_blend::ends(model, rails[1].0, rails[1].1).1.unwrap();
            let old = if rails[1].1 == Sense::Same {
                input.end.unwrap()
            } else {
                input.start.unwrap()
            };
            let cap = run.caps[&old];
            let sense = if model.edge(cap).start == Some(end) {
                Sense::Same
            } else {
                Sense::Reversed
            };
            path.push((cap, sense));
            vec![path]
        };
        let face = crate::extrude::face(
            model,
            run.definition.surface.clone(),
            if run.definition.surface_same_sense {
                Sense::Same
            } else {
                Sense::Reversed
            },
            paths,
            None,
            vec![run.definition.origin.clone()],
        );
        edit.faces.push((run.faces[0], face));
    }
    let result = edge_blend::apply(model, body, &edit)?;
    Ok((result, runs.into_iter().map(|run| run.definition).collect()))
}
