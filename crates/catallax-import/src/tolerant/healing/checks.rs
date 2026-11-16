//! Wire, face, shell and solid checks on the unpublished repair candidate.
use super::super::TolerantUse;
use super::*;
use std::collections::VecDeque;

fn endpoints(body: &TolerantBody, coedge: &TolerantUse) -> [usize; 2] {
    let ends = body.edges[coedge.edge].vertices;
    if coedge.forward {
        ends
    } else {
        [ends[1], ends[0]]
    }
}

pub(super) fn wires(
    body: &mut TolerantBody,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    for s in 0..body.shells.len() {
        for f in 0..body.shells[s].faces.len() {
            for b in 0..body.shells[s].faces[f].loops.len() {
                let boundary = &body.shells[s].faces[f].loops[b];
                let name = body.shells[s].faces[f].name.clone();
                if boundary.edges.is_empty() {
                    if boundary.vertex.is_some() {
                        plan.passed(body, &name, "wire.retained-vertex-loop");
                    } else {
                        plan.unresolved(body, &name, "wire.empty-boundary");
                    }
                    continue;
                }
                let original = boundary.edges.clone();
                let mut remaining: Vec<_> = original.iter().cloned().enumerate().collect();
                let mut ordered = vec![remaining.remove(0)];
                let first = endpoints(body, &ordered[0].1)[0];
                let mut failed = false;
                while !remaining.is_empty() {
                    let wanted = endpoints(body, &ordered.last().expect("first use").1)[1];
                    let mut choices = Vec::new();
                    for (i, (_, coedge)) in remaining.iter().enumerate() {
                        dec.charge_work(1)?;
                        if endpoints(body, coedge)[0] == wanted {
                            choices.push(i);
                        }
                    }
                    if choices.len() != 1 {
                        failed = true;
                        break;
                    }
                    ordered.push(remaining.remove(choices[0]));
                }
                if failed || endpoints(body, &ordered.last().expect("first use").1)[1] != first {
                    plan.unresolved(body, &name, "wire.open-or-ambiguous-order");
                    continue;
                }
                let order: Vec<_> = ordered.iter().map(|(i, _)| *i).collect();
                if order.iter().enumerate().any(|(i, j)| i != *j) {
                    body.shells[s].faces[f].loops[b].edges =
                        ordered.into_iter().map(|(_, coedge)| coedge).collect();
                    plan.changes.push(HealingChange::WireReordered {
                        body: body.name.clone(),
                        face: name.clone(),
                        boundary: b,
                        order,
                    });
                }
                plan.passed(body, &name, "wire.closed-shared-vertex-walk");
            }
        }
    }
    Ok(())
}

fn project_line(
    point: &Vec3,
    origin: &Vec3,
    direction: &Vec3,
    dec: &mut Decider,
) -> Result<Vec3, Error> {
    let length = dot(direction, direction, dec)?;
    if length.is_zero() {
        return Err(failure(
            "import.healing-degenerate-curve",
            "line carrier has zero direction",
        ));
    }
    let offset = sub(point, origin, dec)?;
    let parameter = arithmetic(&dot(&offset, direction, dec)?, &length, '/', dec)?;
    let mut result = origin.clone();
    for i in 0..3 {
        let part = arithmetic(&direction[i], &parameter, '*', dec)?;
        result[i] = arithmetic(&origin[i], &part, '+', dec)?;
    }
    Ok(result)
}

fn plane_squared(point: &Vec3, plane: &Placement, dec: &mut Decider) -> Result<Rational, Error> {
    let length = dot(&plane.axis, &plane.axis, dec)?;
    if length.is_zero() {
        return Err(failure(
            "import.healing-degenerate-surface",
            "plane carrier has zero normal",
        ));
    }
    let distance = dot(&sub(point, &plane.origin, dec)?, &plane.axis, dec)?;
    arithmetic(
        &arithmetic(&distance, &distance, '*', dec)?,
        &length,
        '/',
        dec,
    )
}

fn cross(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Vec3, Error> {
    let mut result = std::array::from_fn(|_| Rational::zero());
    for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
        let p = arithmetic(&a[j], &b[k], '*', dec)?;
        let q = arithmetic(&a[k], &b[j], '*', dec)?;
        result[i] = arithmetic(&p, &q, '-', dec)?;
    }
    Ok(result)
}

fn orientation(
    a: &Vec3,
    b: &Vec3,
    c: &Vec3,
    normal: &Vec3,
    dec: &mut Decider,
) -> Result<Rational, Error> {
    let u = sub(b, a, dec)?;
    let v = sub(c, a, dec)?;
    dot(&cross(&u, &v, dec)?, normal, dec)
}

fn opposite(a: &Rational, b: &Rational) -> bool {
    !a.is_zero() && !b.is_zero() && a.sign() != b.sign()
}

// Validate the immutable source too. A carrier removed by sewing must never
// escape its original incidence checks or acquire an unmeasured correspondence.
pub(super) fn source_geometry(
    body_index: usize,
    source: &TolerantBody,
    candidate: &TolerantBody,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    for edge in &source.edges {
        let Curve::Line(origin, direction) = &edge.curve else {
            continue;
        };
        for (end, vertex) in edge.vertices.iter().enumerate() {
            let point = &source.vertices[*vertex].source.point;
            let projection = project_line(point, origin, direction, dec)?;
            let upper = distance_upper(
                &distance_squared(point, &projection, dec)?,
                &source.options.maximum_mm,
                dec,
            )?;
            plan.measured(
                source,
                &format!("original {} endpoint {end}", edge.name),
                "original-source-point-to-original-line-carrier",
                &upper,
            );
            plan.need(
                body_index,
                ToleranceEntity::Vertex(source.vertices[*vertex].source.name.clone()),
                upper,
                "healing.original-source-vertex-incidence",
            );
        }
    }
    for shell in &source.shells {
        for face in &shell.faces {
            let Surface::Plane(plane) = &face.surface else {
                continue;
            };
            for boundary in &face.loops {
                for coedge in &boundary.edges {
                    let edge = &source.edges[coedge.edge];
                    let Curve::Line(origin, direction) = &edge.curve else {
                        continue;
                    };
                    let mut maximum = Rational::zero();
                    for vertex in edge.vertices {
                        let point = project_line(
                            &source.vertices[vertex].source.point,
                            origin,
                            direction,
                            dec,
                        )?;
                        let upper = distance_upper(
                            &plane_squared(&point, plane, dec)?,
                            &source.options.maximum_mm,
                            dec,
                        )?;
                        if compare(&upper, &maximum, dec)?.is_gt() {
                            maximum = upper;
                        }
                    }
                    plan.measured(
                        source,
                        &format!("original {} on {}", edge.name, face.name),
                        "original-finite-line-segment-to-original-plane",
                        &maximum,
                    );
                    plan.need(
                        body_index,
                        ToleranceEntity::Edge(edge.name.clone()),
                        maximum,
                        "healing.original-source-finite-line-incidence",
                    );
                }
            }
        }
    }
    let sewn: Vec<_> = plan
        .changes
        .iter()
        .filter_map(|change| match change {
            HealingChange::EdgesSewn {
                body,
                retained,
                replaced,
                reversed,
                ..
            } if *body == source.name => Some((retained.clone(), replaced.clone(), *reversed)),
            _ => None,
        })
        .collect();
    for (retained, replaced, reversed) in sewn {
        let original = &source.edges[source.edge_indices[&replaced]];
        let result = &candidate.edges[candidate.edge_indices[&retained]];
        let (Curve::Line(old_origin, old_direction), Curve::Line(origin, direction)) =
            (&original.curve, &result.curve)
        else {
            continue;
        };
        let mut maximum = Rational::zero();
        for end in 0..2 {
            let old = &source.vertices[original.vertices[if reversed { 1 - end } else { end }]]
                .source
                .point;
            let current = &candidate.vertices[result.vertices[end]].source.point;
            let a = project_line(old, old_origin, old_direction, dec)?;
            let b = project_line(current, origin, direction, dec)?;
            let upper = distance_upper(
                &distance_squared(&a, &b, dec)?,
                &source.options.maximum_mm,
                dec,
            )?;
            plan.measured(
                source,
                &format!("{replaced} -> {retained} endpoint {end}"),
                "original-finite-line-to-final-retained-line-endpoint-correspondence",
                &upper,
            );
            if compare(&upper, &maximum, dec)?.is_gt() {
                maximum = upper;
            }
        }
        for name in [retained, replaced] {
            plan.need(
                body_index,
                ToleranceEntity::Edge(name),
                maximum.clone(),
                "healing.measured-finite-source-line-correspondence",
            );
        }
    }
    Ok(())
}

pub(super) fn faces(
    body_index: usize,
    body: &mut TolerantBody,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    let uses = usage(body, dec)?;
    for (index, edge) in body
        .edges
        .iter()
        .enumerate()
        .filter(|(i, _)| !uses[*i].is_empty())
    {
        if let Curve::Line(origin, direction) = &edge.curve {
            for (end, vertex) in edge.vertices.iter().enumerate() {
                let point = &body.vertices[*vertex].source.point;
                let projected = project_line(point, origin, direction, dec)?;
                let upper = distance_upper(
                    &distance_squared(point, &projected, dec)?,
                    &body.options.maximum_mm,
                    dec,
                )?;
                plan.measured(
                    body,
                    &format!("{} endpoint {end}", edge.name),
                    "point-to-retained-line-carrier",
                    &upper,
                );
                plan.need(
                    body_index,
                    ToleranceEntity::Vertex(body.vertices[*vertex].source.name.clone()),
                    upper,
                    "healing.measured-vertex-curve-incidence",
                );
            }
        } else {
            plan.unresolved(
                body,
                &body.edges[index].name,
                "edge.native-parameter-incidence-pending",
            );
        }
    }
    for s in 0..body.shells.len() {
        for f in 0..body.shells[s].faces.len() {
            let face = body.shells[s].faces[f].clone();
            let Surface::Plane(plane) = &face.surface else {
                plan.unresolved(body, &face.name, "face.finite-curved-embedding-pending");
                continue;
            };
            if face.loops.len() > 1 {
                plan.unresolved(body, &face.name, "face.inter-boundary-containment-pending");
            }
            for (b, boundary) in face.loops.iter().enumerate() {
                if boundary.edges.len() < 3
                    || !boundary
                        .edges
                        .iter()
                        .all(|u| matches!(body.edges[u.edge].curve, Curve::Line(..)))
                {
                    plan.unresolved(body, &face.name, "face.nonpolygonal-boundary-pending");
                    continue;
                }
                let points: Vec<_> = boundary
                    .edges
                    .iter()
                    .map(|u| body.vertices[endpoints(body, u)[0]].source.point.clone())
                    .collect();
                for coedge in &boundary.edges {
                    let edge = &body.edges[coedge.edge];
                    let Curve::Line(origin, direction) = &edge.curve else {
                        unreachable!()
                    };
                    let mut maximum = Rational::zero();
                    for vertex in edge.vertices {
                        let point = project_line(
                            &body.vertices[vertex].source.point,
                            origin,
                            direction,
                            dec,
                        )?;
                        let upper = distance_upper(
                            &plane_squared(&point, plane, dec)?,
                            &body.options.maximum_mm,
                            dec,
                        )?;
                        if compare(&upper, &maximum, dec)?.is_gt() {
                            maximum = upper;
                        }
                    }
                    plan.measured(
                        body,
                        &format!("{} on {}", edge.name, face.name),
                        "complete-line-segment-between-projected-endpoints-to-plane",
                        &maximum,
                    );
                    plan.need(
                        body_index,
                        ToleranceEntity::Edge(edge.name.clone()),
                        maximum,
                        "healing.measured-finite-line-plane-incidence",
                    );
                }
                let mut crossing = false;
                for i in 0..points.len() {
                    for j in i + 1..points.len() {
                        dec.charge_work(1)?;
                        if j == i + 1 || (i == 0 && j + 1 == points.len()) {
                            continue;
                        }
                        let a = &points[i];
                        let c = &points[j];
                        let d = &points[(j + 1) % points.len()];
                        let p = &points[(i + 1) % points.len()];
                        let ab_c = orientation(a, p, c, &plane.axis, dec)?;
                        let ab_d = orientation(a, p, d, &plane.axis, dec)?;
                        let cd_a = orientation(c, d, a, &plane.axis, dec)?;
                        let cd_b = orientation(c, d, p, &plane.axis, dec)?;
                        // A zero determinant is conservatively unresolved too;
                        // this check never turns a touching polygon into a pass.
                        crossing |= ab_c.is_zero()
                            || ab_d.is_zero()
                            || cd_a.is_zero()
                            || cd_b.is_zero()
                            || (opposite(&ab_c, &ab_d) && opposite(&cd_a, &cd_b));
                    }
                }
                if crossing {
                    plan.unresolved(
                        body,
                        &face.name,
                        "face.crossing-or-collinear-nonadjacent-boundaries",
                    );
                    continue;
                }
                let mut area = Rational::zero();
                for i in 0..points.len() {
                    let value = dot(
                        &cross(&points[i], &points[(i + 1) % points.len()], dec)?,
                        &plane.axis,
                        dec,
                    )?;
                    area = arithmetic(&area, &value, '+', dec)?;
                }
                if area.is_zero() {
                    plan.unresolved(body, &face.name, "face.zero-projected-area");
                    continue;
                }
                let positive =
                    (area.sign() == catallax_numbers::Sign::Positive) == boundary.orientation;
                let desired = boundary.outer == face.same_sense;
                if positive != desired {
                    body.shells[s].faces[f].loops[b].orientation = !boundary.orientation;
                    plan.changes.push(HealingChange::BoundaryReversed {
                        body: body.name.clone(),
                        face: face.name.clone(),
                        boundary: b,
                    });
                }
                plan.passed(
                    body,
                    &face.name,
                    "face.simple-nondegenerate-oriented-planar-line-boundary",
                );
            }
        }
    }
    Ok(())
}

pub(super) fn shells(
    body: &mut TolerantBody,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    let uses = usage(body, dec)?;
    for s in 0..body.shells.len() {
        let count = body.shells[s].faces.len();
        let mut adjacent = vec![Vec::new(); count];
        let mut closed = true;
        for edge in &uses {
            let local: Vec<_> = edge.iter().filter(|edge| edge.shell == s).collect();
            if local.is_empty() {
                continue;
            }
            dec.charge_work(1)?;
            if local.len() != 2 {
                closed = false;
                continue;
            }
            let (a, direction_a) = (local[0].face, local[0].forward);
            let (b, direction_b) = (local[1].face, local[1].forward);
            if a == b {
                if direction_a == direction_b {
                    closed = false;
                    plan.unresolved(body, &body.shells[s].name, "shell.inconsistent-seam-use");
                }
            } else {
                adjacent[a].push((b, direction_a == direction_b));
                adjacent[b].push((a, direction_a == direction_b));
            }
        }
        let mut flipped = vec![None; count];
        let mut components = 0;
        for first in 0..count {
            if flipped[first].is_some() {
                continue;
            }
            components += 1;
            let mut queue = VecDeque::from([first]);
            let mut component = Vec::new();
            let mut consistent = true;
            flipped[first] = Some(false);
            while let Some(face) = queue.pop_front() {
                dec.charge_work(1)?;
                component.push(face);
                for &(other, different) in &adjacent[face] {
                    let expected = flipped[face].expect("queued face") ^ different;
                    match flipped[other] {
                        Some(actual) => consistent &= actual == expected,
                        None => {
                            flipped[other] = Some(expected);
                            queue.push_back(other);
                        }
                    }
                }
            }
            if !consistent {
                plan.unresolved(body, &body.shells[s].name, "shell.non-orientable-edge-uses");
                continue;
            }
            for index in component {
                if flipped[index] == Some(true) {
                    reverse_face(body, s, index, "shell-consistent-edge-uses", plan, dec)?;
                }
            }
        }
        let name = body.shells[s].name.clone();
        if closed {
            plan.passed(body, &name, "shell.two-opposite-uses-per-active-edge");
        } else {
            plan.unresolved(body, &name, "shell.open-or-nonmanifold");
        }
        if components == 1 {
            plan.passed(body, &name, "shell.connected-face-incidence");
        } else {
            plan.unresolved(body, &name, "shell.multiple-face-components");
        }
    }
    Ok(())
}

fn reverse_face(
    body: &mut TolerantBody,
    s: usize,
    f: usize,
    cause: &str,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    dec.charge_work(1 + body.shells[s].faces[f].loops.len() as u64)?;
    let face = &mut body.shells[s].faces[f];
    face.same_sense = !face.same_sense;
    for boundary in &mut face.loops {
        boundary.orientation = !boundary.orientation;
    }
    plan.changes.push(HealingChange::FaceReversed {
        body: body.name.clone(),
        face: face.name.clone(),
        cause: cause.into(),
    });
    Ok(())
}

pub(super) fn solid(
    body: &mut TolerantBody,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    if !body.declared_solid {
        plan.passed(body, &body.name, "solid.source-surface-kind-retained");
        return Ok(());
    }
    if plan.issues.iter().any(|issue| issue.body == body.name) {
        plan.unresolved(body, &body.name, "solid.unresolved-prerequisites");
        return Ok(());
    }
    let mut volume = Rational::zero();
    for shell in &body.shells {
        for face in &shell.faces {
            for boundary in &face.loops {
                let points: Vec<_> = boundary
                    .edges
                    .iter()
                    .map(|u| body.vertices[endpoints(body, u)[0]].source.point.clone())
                    .collect();
                if points.len() < 3 {
                    plan.unresolved(body, &body.name, "solid.volume-domain-pending");
                    return Ok(());
                }
                for i in 1..points.len() - 1 {
                    let part = dot(&points[0], &cross(&points[i], &points[i + 1], dec)?, dec)?;
                    volume = arithmetic(
                        &volume,
                        &part,
                        if boundary.orientation { '+' } else { '-' },
                        dec,
                    )?;
                }
            }
        }
    }
    if volume.is_zero() {
        plan.unresolved(body, &body.name, "solid.zero-oriented-polyhedral-volume");
        return Ok(());
    }
    if volume.sign() == catallax_numbers::Sign::Negative {
        for s in 0..body.shells.len() {
            for f in 0..body.shells[s].faces.len() {
                reverse_face(
                    body,
                    s,
                    f,
                    "solid-positive-oriented-polyhedral-volume",
                    plan,
                    dec,
                )?;
            }
        }
    }
    plan.passed(
        body,
        &body.name,
        "solid.positive-oriented-polyhedral-volume",
    );
    // Cross-face separation and a checked B-rep remain admission obligations.
    plan.unresolved(body, &body.name, "solid.cross-face-validity-pending");
    Ok(())
}
