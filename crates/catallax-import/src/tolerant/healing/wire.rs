//! Ordered-wire connection and reorder, following ShapeFix_Wire's separation
//! between topology order and measured endpoint connection.
use super::*;

struct VertexJoin {
    source: usize,
    retained: usize,
    upper: Rational,
}
struct OrderedUse {
    original: usize,
    coedge: TolerantUse,
    reversed: bool,
}

fn resolve(
    aliases: &BTreeMap<usize, usize>,
    mut index: usize,
    dec: &mut Decider,
) -> Result<usize, Error> {
    let mut steps = 0;
    while let Some(next) = aliases.get(&index) {
        dec.charge_work(1)?;
        steps += 1;
        if steps > aliases.len() {
            return Err(failure(
                "import.healing-topology",
                "cyclic wire endpoint plan",
            ));
        }
        index = *next;
    }
    Ok(index)
}

fn connection_plan(
    body: &TolerantBody,
    uses: &[TolerantUse],
    dec: &mut Decider,
) -> Result<Option<Vec<VertexJoin>>, Error> {
    let mut aliases = BTreeMap::new();
    let mut joins = Vec::new();
    let squared_cap = arithmetic(&body.options.maximum_mm, &body.options.maximum_mm, '*', dec)?;
    for i in 0..uses.len() {
        dec.charge_work(1)?;
        let a = resolve(&aliases, endpoints(body, &uses[i])[1], dec)?;
        let b = resolve(
            &aliases,
            endpoints(body, &uses[(i + 1) % uses.len()])[0],
            dec,
        )?;
        if a == b {
            continue;
        }
        let squared = distance_squared(
            &body.vertices[a].source.point,
            &body.vertices[b].source.point,
            dec,
        )?;
        if compare(&squared, &squared_cap, dec)?.is_gt() {
            return Ok(None);
        }
        // Preserve the first source anchor at closure; preceding connections
        // preserve the already-walked endpoint. No nearest entity is searched.
        let (source, retained) = if i + 1 == uses.len() { (a, b) } else { (b, a) };
        let upper = distance_upper(&squared, &body.options.maximum_mm, dec)?;
        aliases.insert(source, retained);
        joins.push(VertexJoin {
            source,
            retained,
            upper,
        });
    }
    Ok(Some(joins))
}

fn reorder(
    body: &TolerantBody,
    uses: &[TolerantUse],
    dec: &mut Decider,
) -> Result<Option<Vec<OrderedUse>>, Error> {
    let mut remaining: Vec<_> = uses
        .iter()
        .cloned()
        .enumerate()
        .map(|(original, coedge)| OrderedUse {
            original,
            coedge,
            reversed: false,
        })
        .collect();
    let mut ordered = vec![remaining.remove(0)];
    while !remaining.is_empty() {
        let wanted = endpoints(body, &ordered.last().expect("first use").coedge)[1];
        let mut choices = Vec::new();
        for (index, item) in remaining.iter().enumerate() {
            dec.charge_work(1)?;
            if endpoints(body, &item.coedge)[0] == wanted {
                choices.push((index, false));
            }
        }
        if choices.is_empty() {
            for (index, item) in remaining.iter().enumerate() {
                dec.charge_work(1)?;
                if endpoints(body, &item.coedge)[1] == wanted {
                    choices.push((index, true));
                }
            }
        }
        if choices.len() != 1 {
            return Ok(None);
        }
        let (index, reversed) = choices[0];
        let mut item = remaining.remove(index);
        item.coedge.forward ^= reversed;
        item.reversed = reversed;
        ordered.push(item);
    }
    Ok(Some(ordered))
}

fn apply_join(
    body_index: usize,
    body: &mut TolerantBody,
    join: VertexJoin,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    let source = body.vertices[join.source].source.name.clone();
    let retained = body.vertices[join.retained].source.name.clone();
    plan.measured(
        body,
        &format!("wire {source} -> {retained}"),
        "ordered-adjacent-source-endpoint-gap",
        &join.upper,
    );
    for name in [&source, &retained] {
        plan.need(
            body_index,
            ToleranceEntity::Vertex(name.clone()),
            join.upper.clone(),
            "healing.measured-wire-endpoint-connection",
        );
    }
    for edge in &mut body.edges {
        for vertex in &mut edge.vertices {
            dec.charge_work(1)?;
            if *vertex == join.source {
                *vertex = join.retained;
            }
        }
    }
    for shell in &mut body.shells {
        for face in &mut shell.faces {
            for boundary in &mut face.loops {
                dec.charge_work(1)?;
                if boundary.vertex == Some(join.source) {
                    boundary.vertex = Some(join.retained);
                }
            }
        }
    }
    plan.changes.push(HealingChange::VerticesJoined {
        body: body.name.clone(),
        source,
        retained,
        distance_upper_mm: join.upper.to_string(),
    });
    Ok(())
}

pub(super) fn fix(
    body_index: usize,
    source: &TolerantBody,
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
                let mut ordered: Vec<_> = original
                    .iter()
                    .cloned()
                    .enumerate()
                    .map(|(original, coedge)| OrderedUse {
                        original,
                        coedge,
                        reversed: false,
                    })
                    .collect();
                let connected = original.iter().enumerate().all(|(index, usage)| {
                    endpoints(body, usage)[1]
                        == endpoints(body, &original[(index + 1) % original.len()])[0]
                });
                let reordered_uses;
                let uses = if !connected && let Some(reordered) = reorder(body, &original, dec)? {
                    ordered = reordered;
                    reordered_uses = ordered
                        .iter()
                        .map(|item| item.coedge.clone())
                        .collect::<Vec<_>>();
                    &reordered_uses
                } else {
                    &original
                };
                let Some(joins) = connection_plan(body, uses, dec)? else {
                    plan.unresolved(body, &name, "wire.open-or-over-cap-connection");
                    continue;
                };
                // All connections and final closure have been planned before the first
                // write. Whole-body source/cap/topology checks still precede publication.
                for join in joins {
                    apply_join(body_index, body, join, plan, dec)?;
                }
                let order: Vec<_> = ordered.iter().map(|item| item.original).collect();
                let reversed_use_indices: Vec<_> = ordered
                    .iter()
                    .filter(|item| item.reversed)
                    .map(|item| item.original)
                    .collect();
                if order.iter().enumerate().any(|(i, j)| i != *j)
                    || !reversed_use_indices.is_empty()
                {
                    plan.changes.push(HealingChange::WireReordered {
                        body: body.name.clone(),
                        face: name.clone(),
                        boundary: b,
                        order,
                        reversed_use_indices,
                    });
                }
                body.shells[s].faces[f].loops[b].edges =
                    ordered.into_iter().map(|item| item.coedge).collect();
                plan.passed(body, &name, "wire.closed-shared-vertex-walk");
            }
        }
    }
    remove_small_lines(body_index, source, body, plan, dec)?;
    Ok(())
}

fn remove_small_lines(
    body_index: usize,
    source: &TolerantBody,
    body: &mut TolerantBody,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    let uses = usage(body, dec)?;
    let minimum = body.options.initial_mm.clone();
    let minimum_squared = arithmetic(&minimum, &minimum, '*', dec)?;
    for (index, incidence) in uses.iter().enumerate() {
        if incidence.is_empty() {
            continue;
        }
        let name = body.edges[index].name.clone();
        let original = &source.edges[source.edge_indices[&name]];
        let Curve::Line(origin, direction) = &original.curve else {
            continue;
        };
        let a = &source.vertices[original.vertices[0]].source.point;
        let b = &source.vertices[original.vertices[1]].source.point;
        let source_chord = distance_squared(a, b, dec)?;
        if compare(&source_chord, &minimum_squared, dec)?.is_gt() {
            continue;
        }
        // Smallness uses the immutable finite source interval, not a chord
        // accidentally collapsed by an earlier wire or sewing repair.
        let p = checks::project_line(a, origin, direction, dec)?;
        let q = checks::project_line(b, origin, direction, dec)?;
        let length = distance_squared(&p, &q, dec)?;
        if compare(&length, &minimum_squared, dec)?.is_gt() {
            continue;
        }
        let mut safe = true;
        for shell in &body.shells {
            for face in &shell.faces {
                for boundary in &face.loops {
                    dec.charge_work(1 + boundary.edges.len() as u64)?;
                    if !boundary.edges.iter().any(|u| u.edge == index) {
                        continue;
                    }
                    safe &= matches!(face.surface, Surface::Plane(_))
                        && boundary.edges.iter().any(|u| u.edge != index);
                }
            }
        }
        if !safe {
            plan.unresolved(
                body,
                &name,
                "wire.small-line-requires-curved-or-empty-boundary-analysis",
            );
            continue;
        }
        let [retained, removed] = body.edges[index].vertices;
        if retained != removed {
            let squared = distance_squared(
                &body.vertices[removed].source.point,
                &body.vertices[retained].source.point,
                dec,
            )?;
            let upper = distance_upper(&squared, &body.options.maximum_mm, dec)?;
            apply_join(
                body_index,
                body,
                VertexJoin {
                    source: removed,
                    retained,
                    upper,
                },
                plan,
                dec,
            )?;
        }
        for shell in &mut body.shells {
            for face in &mut shell.faces {
                for boundary in &mut face.loops {
                    dec.charge_work(1 + boundary.edges.len() as u64)?;
                    boundary.edges.retain(|u| u.edge != index);
                }
            }
        }
        let upper = distance_upper(&length, &minimum, dec)?;
        plan.measured(
            body,
            &name,
            "complete-original-finite-line-length-before-removal",
            &upper,
        );
        plan.changes.push(HealingChange::SmallLineRemoved {
            body: body.name.clone(),
            edge: name,
            retained_vertex: body.vertices[retained].source.name.clone(),
            source_length_upper_mm: upper.to_string(),
        });
    }
    Ok(())
}
