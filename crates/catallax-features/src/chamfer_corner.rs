//! A planar corner sewn from the three incident bevel contact intersections.
use super::*;
use catallax_shapes::Plane;

pub(super) fn three(
    model: &mut Model,
    body: BodyId,
    old: VertexId,
    active: [usize; 3],
    runs: &mut [Run],
    step: &StepId,
    dec: &mut Decider,
) -> Result<FaceId, Refusal> {
    if active.iter().any(|i| !runs[*i].definition.convex) {
        return Err(bad("a mixed concave three-strip corner needs continuation"));
    }
    let faces: BTreeSet<_> = active.iter().flat_map(|i| runs[*i].faces).collect();
    let mut points = Vec::with_capacity(3);
    let mut vertices = BTreeMap::new();
    for face in faces {
        plane(model, face)?;
        let touching: Vec<_> = active
            .iter()
            .copied()
            .filter(|i| runs[*i].faces.contains(&face))
            .collect();
        let [a, b] = touching.as_slice() else {
            return Err(bad("three strips must meet pairwise on three input planes"));
        };
        let point = line_line(rail(&runs[*a], face)?, rail(&runs[*b], face)?, dec)?;
        let vertex = vertex(model, point.clone());
        points.push((vertex, point));
        vertices.insert(face, vertex);
        for i in [*a, *b] {
            runs[i].ends.insert((old, face), vertex);
        }
    }
    let [(a, p), (b, q), (c, r)] = points.as_slice() else {
        return Err(bad("a three-strip corner needs three incident input faces"));
    };
    let reference = q.minus(p);
    let mut normal = reference.cross(&r.minus(p));
    if sign(&normal.norm_squared(), dec)? == Sign::Zero {
        return Err(bad("the three corner contacts collapse to a line"));
    }
    let old_point = model.geom().position(model.vertex(old).point);
    let reverse = match sign(&old_point.minus(p).dot(&normal), dec)? {
        Sign::Positive => false,
        Sign::Negative => true,
        Sign::Zero => return Err(bad("the corner plane has no removed-material side")),
    };
    if reverse {
        normal = normal.scaled(&Real::from_i64(-1));
    }
    let surface = Surface::Plane(Plane {
        origin: p.clone(),
        normal,
        ref_dir: reference,
    });
    let mut edges = Vec::with_capacity(3);
    for i in active {
        let run = &mut runs[i];
        let edge = connector(
            model,
            &run.definition.surface,
            &surface,
            vertices[&run.faces[0]],
            vertices[&run.faces[1]],
            dec,
        )?;
        run.caps.insert(old, edge);
        edges.push(edge);
    }
    let order = if reverse { [*a, *c, *b] } else { [*a, *b, *c] };
    let mut path = Vec::with_capacity(3);
    for i in 0..3 {
        let ends = (Some(order[i]), Some(order[(i + 1) % 3]));
        let edge = *edges
            .iter()
            .find(|edge| {
                edge_blend::ends(model, **edge, Sense::Same) == ends
                    || edge_blend::ends(model, **edge, Sense::Reversed) == ends
            })
            .ok_or_else(|| bad("corner connectors do not form one triangle"))?;
        let sense = if edge_blend::ends(model, edge, Sense::Same) == ends {
            Sense::Same
        } else {
            Sense::Reversed
        };
        path.push((edge, sense));
    }
    let name = model
        .vertex_name(body, old)
        .ok_or_else(|| bad("corner vertex needs a persistent input name"))?;
    Ok(crate::extrude::face(
        model,
        surface,
        Sense::Same,
        vec![path],
        None,
        vec![Origin {
            step: step.clone(),
            role: Role::Corner { vertex: name },
        }],
    ))
}
