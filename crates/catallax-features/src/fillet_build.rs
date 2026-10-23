//! Rail-first elementary fillet construction over an immutable input skin.

use std::collections::{BTreeMap, BTreeSet};

use catallax_intersect::{Component, SurfacePair};
use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign};
use catallax_recipe::{EdgeSelection, Origin, Role, StepId};
use catallax_shapes::{Curve, Line, Point3, PointDef, Surface};
use catallax_topology::{BodyId, Edge, EdgeId, FaceId, Loop, Model, Sense, Vertex, VertexId};

use crate::{
    edge_blend::{self, Replacement},
    extrude::reverse,
    fillet::{FilletDefinition, separated},
    profile::sign,
};

type Uses = BTreeMap<EdgeId, Vec<(FaceId, Sense)>>;

fn bad(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-topology", detail)
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
fn resolve(model: &Model, selection: &EdgeSelection, uses: &Uses) -> Result<Vec<EdgeId>, Refusal> {
    match selection {
        EdgeSelection::All => Ok(uses.keys().copied().collect()),
        EdgeSelection::Refs(refs) => {
            let mut selected = BTreeSet::new();
            for reference in refs {
                if reference.near.is_some() {
                    return Err(bad(
                        "nearest-edge disambiguation needs an exact trimmed-edge distance query",
                    ));
                }
                let mut found = false;
                for (edge, incident) in uses {
                    let [(a, _), (b, _)] = incident.as_slice() else {
                        continue;
                    };
                    if !((matching(model, *a, &reference.between[0])
                        && matching(model, *b, &reference.between[1]))
                        || (matching(model, *b, &reference.between[0])
                            && matching(model, *a, &reference.between[1])))
                    {
                        continue;
                    }
                    let e = model.edge(*edge);
                    let touches = reference.touching.iter().all(|pattern| {
                        uses.iter().any(|(other, incident)| {
                            let other = model.edge(*other);
                            [e.start, e.end]
                                .into_iter()
                                .flatten()
                                .any(|v| other.start == Some(v) || other.end == Some(v))
                                && incident.iter().any(|(f, _)| matching(model, *f, pattern))
                        })
                    });
                    if touches {
                        selected.insert(*edge);
                        found = true;
                    }
                }
                if !found {
                    return Err(Refusal::new(
                        "fillet.lost-edge",
                        format!("no input edge resolves {reference:?}"),
                    ));
                }
            }
            Ok(selected.into_iter().collect())
        }
    }
}

struct Run {
    edge: EdgeId,
    faces: [FaceId; 2],
    senses: [Sense; 2],
    definition: FilletDefinition,
    ends: BTreeMap<(VertexId, FaceId), VertexId>,
    caps: BTreeMap<VertexId, EdgeId>,
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
        return Err(bad("a fillet edge must have two incident faces"));
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
                .ok_or_else(|| bad("circular elementary fillets need one plane"))?;
            let coedge = incident.iter().find(|(f, _)| *f == plane).unwrap().1;
            let orientation = if coedge == e.sense { 1 } else { -1 };
            let n = outward(model, plane, &c.center, dec)?;
            Ok(sign(&(c.normal.dot(&n) * Real::from_i64(orientation)), dec)? == Sign::Positive)
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
        _ => Err(bad("elementary fillets need a line or circle")),
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

fn connector(
    model: &mut Model,
    a: &Surface,
    b: &Surface,
    from: VertexId,
    to: VertexId,
    sectors: &[&FilletDefinition],
    dec: &mut Decider,
) -> Result<EdgeId, Refusal> {
    let a = model.geom_mut().add_surface(a.clone());
    let b = model.geom_mut().add_surface(b.clone());
    let meeting = catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value;
    let SurfacePair::Meets(components) = meeting else {
        return Err(bad("blend end supports have no regular common curve"));
    };
    let p = model.vertex(from).point;
    let q = model.vertex(to).point;
    for component in components {
        let Component::Curve { curve, .. } = component else {
            continue;
        };
        let on = |model: &Model, p, dec: &mut Decider| -> Result<bool, Refusal> {
            match catallax_intersect::on_curve(model.geom(), curve, p, dec) {
                Outcome::ProvenCoincident(_) => Ok(true),
                Outcome::Proven(_) => Ok(false),
                Outcome::Refused(e) => Err(e),
                Outcome::NotExact(e) => Err(Refusal::new(e.code, e.detail)),
            }
        };
        if !on(model, p, dec)? || !on(model, q, dec)? {
            continue;
        }
        for (sense, start, end) in [(Sense::Same, p, q), (Sense::Reversed, q, p)] {
            let sample =
                catallax_intersect::arc_point(model.geom(), curve, Some(start), Some(end), dec)?
                    .value;
            let mut inside = true;
            for sector in sectors {
                inside &= sector.in_sector(&sample, dec)?;
            }
            if inside {
                return Ok(model.add_edge(Edge {
                    curve,
                    sense,
                    start: Some(from),
                    end: Some(to),
                }));
            }
        }
    }
    Err(bad("no end-curve arc lies in every active contact sector"))
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
        let other = model.geom().surface(model.face(face).surface).clone();
        if separated(&run.definition.surface, &other, dec)? {
            continue;
        }
        let a = model.geom_mut().add_surface(run.definition.surface.clone());
        let b = model.geom_mut().add_surface(other);
        if !matches!(
            catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value,
            SurfacePair::Disjoint
        ) {
            return Err(Refusal::new(
                "fillet.unsupported-overflow",
                "a non-end neighbour intersects the blend carrier; its trimmed contact continuation must be resolved before publishing",
            ));
        }
    }
    Ok(())
}

pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    selection: &EdgeSelection,
    radius: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(BodyId, Vec<FilletDefinition>), Refusal> {
    if sign(radius, dec)? != Sign::Positive {
        return Err(Refusal::new(
            "fillet.non-positive-radius",
            "fillet radius must be positive",
        ));
    }
    let uses = uses(model, body);
    let selected = resolve(model, selection, &uses)?;
    if selected.is_empty() {
        return Ok((body, Vec::new()));
    }
    let mut runs = Vec::with_capacity(selected.len());
    for edge in selected {
        let incident = &uses[&edge];
        let [(a, sa), (b, sb)] = incident.as_slice() else {
            return Err(bad("a selected edge must have two incident faces"));
        };
        let name = model
            .edge_name(body, edge)
            .ok_or_else(|| bad("a selected edge needs an unambiguous origin-pair name"))?;
        let faces = [*a, *b];
        let definition = FilletDefinition::new(
            Origin {
                step: step.clone(),
                role: Role::Blend { edge: name },
            },
            faces.map(|f| model.geom().surface(model.face(f).surface).clone()),
            faces.map(|f| model.face(f).sense == Sense::Same),
            model.geom().curve(model.edge(edge).curve).clone(),
            radius.clone(),
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
        });
    }
    for run in &runs {
        guard(model, body, run, &uses, dec)?;
    }
    for (i, a) in runs.iter().enumerate() {
        for b in &runs[i + 1..] {
            let a_edge = model.edge(a.edge);
            let b_edge = model.edge(b.edge);
            if [a_edge.start, a_edge.end]
                .into_iter()
                .flatten()
                .any(|v| b_edge.start == Some(v) || b_edge.end == Some(v))
            {
                continue;
            }
            if separated(&a.definition.surface, &b.definition.surface, dec)? {
                continue;
            }
            let a = model.geom_mut().add_surface(a.definition.surface.clone());
            let b = model.geom_mut().add_surface(b.definition.surface.clone());
            if !matches!(
                catallax_intersect::intersect_surfaces(model.geom_mut(), a, b, dec)?.value,
                SurfacePair::Disjoint
            ) {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "nonincident blend carriers meet and require simultaneous contact-region trimming",
                ));
            }
        }
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
                let a = line_plane(rail(&runs[*i], faces[0])?, plane(model, third)?, dec)?;
                let b = line_plane(rail(&runs[*i], faces[1])?, plane(model, third)?, dec)?;
                let [a, b] = [vertex(model, a), vertex(model, b)];
                for (face, v) in faces.into_iter().zip([a, b]) {
                    runs[*i].ends.insert((old, face), v);
                }
                let surface = runs[*i].definition.surface.clone();
                let third_surface = model.geom().surface(model.face(third).surface).clone();
                let cap = connector(
                    model,
                    &surface,
                    &third_surface,
                    a,
                    b,
                    &[&runs[*i].definition],
                    dec,
                )?;
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
                    &[&runs[*i].definition, &runs[*j].definition],
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
                        "fillet.unsupported-overflow",
                        "a contact foot leaves the input face and needs neighbour continuation",
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
            if run.definition.convex {
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
