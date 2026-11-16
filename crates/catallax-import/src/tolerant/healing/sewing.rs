//! Conservative finite-edge sewing: mutually unique free-boundary matches.
use super::*;

#[derive(Clone)]
struct Match {
    other: usize,
    reversed: bool,
    gaps: [Rational; 2],
}

fn parallel(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<bool, Error> {
    for (i, j) in [(0, 1), (1, 2), (2, 0)] {
        let left = arithmetic(&a[i], &b[j], '*', dec)?;
        let right = arithmetic(&a[j], &b[i], '*', dec)?;
        if compare(&left, &right, dec)?.is_ne() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn correspondence(
    body: &TolerantBody,
    a: usize,
    b: usize,
    dec: &mut Decider,
) -> Result<Option<Match>, Error> {
    let left = &body.edges[a];
    let right = &body.edges[b];
    let (Curve::Line(_, left_direction), Curve::Line(_, right_direction)) =
        (&left.curve, &right.curve)
    else {
        // Carrier and endpoint equality do not identify a finite nonlinear
        // domain: opposite senses can select complementary arcs. Retain both
        // edges until a finite-domain correspondence check is implemented.
        return Ok(None);
    };
    if !parallel(left_direction, right_direction, dec)? {
        return Ok(None);
    }
    let cap = &body.options.maximum_mm;
    let cap_squared = arithmetic(cap, cap, '*', dec)?;
    let mut alternatives = Vec::new();
    for reversed in [false, true] {
        let mut distances = Vec::new();
        for i in 0..2 {
            let p = &body.vertices[left.vertices[i]].source.point;
            let q = &body.vertices[right.vertices[if reversed { 1 - i } else { i }]]
                .source
                .point;
            let squared = distance_squared(p, q, dec)?;
            if compare(&squared, &cap_squared, dec)?.is_gt() {
                break;
            }
            distances.push(squared);
        }
        if distances.len() == 2 {
            alternatives.push(Match {
                other: b,
                reversed,
                gaps: [distances[0].clone(), distances[1].clone()],
            });
        }
    }
    // A tiny edge or coincident endpoints may admit both directions. Retain it
    // for a later small/degenerate-edge analysis rather than choosing a branch.
    if alternatives.len() != 1 {
        return Ok(None);
    }
    Ok(alternatives.pop())
}

fn root(aliases: &[usize], mut v: usize, dec: &mut Decider) -> Result<usize, Error> {
    while aliases[v] != v {
        dec.charge_work(1)?;
        v = aliases[v];
    }
    Ok(v)
}

pub(super) fn sew(
    body_index: usize,
    body: &mut TolerantBody,
    plan: &mut Plan,
    dec: &mut Decider,
) -> Result<(), Error> {
    let before = usage(body, dec)?;
    let free: Vec<_> = before
        .iter()
        .enumerate()
        .filter(|(_, uses)| uses.len() == 1)
        .map(|(i, _)| i)
        .collect();
    let mut choices: Vec<Vec<Match>> = vec![Vec::new(); body.edges.len()];
    for &edge in &free {
        dec.charge_work(1)?;
        if !matches!(body.edges[edge].curve, Curve::Line(..)) {
            plan.unresolved(
                body,
                &body.edges[edge].name,
                "sewing.nonlinear-finite-correspondence-pending",
            );
        }
    }
    for (index, &a) in free.iter().enumerate() {
        for &b in &free[index + 1..] {
            dec.charge_work(1)?;
            if before[a][0].shell != before[b][0].shell || before[a][0].face == before[b][0].face {
                continue;
            }
            if let Some(found) = correspondence(body, a, b, dec)? {
                choices[a].push(found.clone());
                choices[b].push(Match {
                    other: a,
                    reversed: found.reversed,
                    gaps: found.gaps,
                });
            }
        }
    }
    let mut aliases: Vec<_> = (0..body.vertices.len()).collect();
    for &a in &free {
        if choices[a].len() != 1 {
            continue;
        }
        let found = &choices[a][0];
        let b = found.other;
        if a >= b || choices[b].len() != 1 || choices[b][0].other != a {
            continue;
        }
        let left = body.edges[a].clone();
        let right = body.edges[b].clone();
        let mut gaps = Vec::new();
        for squared in &found.gaps {
            gaps.push(distance_upper(squared, &body.options.maximum_mm, dec)?);
        }
        for i in 0..2 {
            let keep = root(&aliases, left.vertices[i], dec)?;
            let remove = root(
                &aliases,
                right.vertices[if found.reversed { 1 - i } else { i }],
                dec,
            )?;
            if keep != remove {
                aliases[remove] = keep;
            }
        }
        for shell in &mut body.shells {
            for face in &mut shell.faces {
                for boundary in &mut face.loops {
                    for coedge in &mut boundary.edges {
                        dec.charge_work(1)?;
                        if coedge.edge == b {
                            coedge.edge = a;
                            coedge.forward ^= found.reversed;
                        }
                    }
                }
            }
        }
        let maximum = gaps[0].clone().max(gaps[1].clone());
        for name in [&left.name, &right.name] {
            plan.need(
                body_index,
                ToleranceEntity::Edge(name.clone()),
                maximum.clone(),
                "healing.finite-edge-sewing",
            );
        }
        plan.changes.push(HealingChange::EdgesSewn {
            body: body.name.clone(),
            retained: left.name,
            replaced: right.name,
            reversed: found.reversed,
            endpoint_distance_upper_mm: [gaps[0].to_string(), gaps[1].to_string()],
            checks: vec![
                "parallel-line-carriers-and-finite-endpoint-correspondence".into(),
                "mutually-unique-free-boundary-match".into(),
            ],
        });
    }
    // Validate final, possibly transitive joins against the immutable source
    // points. A chain of small gaps cannot silently exceed the cap in aggregate.
    for index in 0..aliases.len() {
        dec.charge_work(1)?;
        let retained = root(&aliases, index, dec)?;
        if retained == index {
            continue;
        }
        let squared = distance_squared(
            &body.vertices[index].source.point,
            &body.vertices[retained].source.point,
            dec,
        )?;
        let upper = distance_upper(&squared, &body.options.maximum_mm, dec)?;
        for vertex in [index, retained] {
            plan.need(
                body_index,
                ToleranceEntity::Vertex(body.vertices[vertex].source.name.clone()),
                upper.clone(),
                "healing.shared-vertex-join",
            );
        }
        plan.changes.push(HealingChange::VerticesJoined {
            body: body.name.clone(),
            source: body.vertices[index].source.name.clone(),
            retained: body.vertices[retained].source.name.clone(),
            distance_upper_mm: upper.to_string(),
        });
    }
    for edge in &mut body.edges {
        for vertex in &mut edge.vertices {
            dec.charge_work(1)?;
            *vertex = root(&aliases, *vertex, dec)?;
        }
    }
    for shell in &mut body.shells {
        for face in &mut shell.faces {
            for boundary in &mut face.loops {
                if let Some(vertex) = &mut boundary.vertex {
                    dec.charge_work(1)?;
                    *vertex = root(&aliases, *vertex, dec)?;
                }
            }
        }
    }
    let after = usage(body, dec)?;
    if after.iter().any(|uses| uses.len() > 2) && !before.iter().any(|uses| uses.len() > 2) {
        return Err(failure(
            "import.healing-topology",
            "sewing would introduce a non-manifold edge",
        ));
    }
    plan.passed(
        body,
        &body.name,
        "sewing.mutually-unique-finite-boundaries-with-capped-endpoint-correspondence",
    );
    Ok(())
}
