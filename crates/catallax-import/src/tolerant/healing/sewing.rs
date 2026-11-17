//! Conservative finite-edge sewing: mutually unique free-boundary matches.
use super::*;

#[derive(Clone)]
pub(super) struct Match {
    pub(super) other: usize,
    pub(super) reversed: bool,
    pub(super) gaps: [Rational; 2],
    pub(super) check: &'static str,
    pub(super) maximum_sample_mm: Option<Rational>,
}

fn correspondence(
    body: &TolerantBody,
    a: usize,
    b: usize,
    dec: &mut Decider,
) -> Result<Option<Match>, Error> {
    let left = &body.edges[a];
    let right = &body.edges[b];
    let exact_arc = match (&left.curve, &right.curve) {
        (Curve::Line(..), Curve::Line(..)) => false,
        (Curve::Circle(..), Curve::Circle(..)) | (Curve::Ellipse(..), Curve::Ellipse(..)) => {
            charge_curve(&left.curve, dec)?;
            charge_curve(&right.curve, dec)?;
            if left.curve != right.curve {
                return Ok(None);
            }
            true
        }
        _ => return Ok(None),
    };
    let cap = &body.options.maximum_mm;
    let cap_squared = arithmetic(cap, cap, '*', dec)?;
    let mut alternatives = Vec::new();
    for reversed in [false, true] {
        if exact_arc && (left.curve_same_sense ^ reversed) != right.curve_same_sense {
            continue;
        }
        let mut distances = Vec::new();
        for i in 0..2 {
            let p = &body.vertices[left.vertices[i]].source.point;
            let q = &body.vertices[right.vertices[if reversed { 1 - i } else { i }]]
                .source
                .point;
            let squared = distance_squared(p, q, dec)?;
            if (exact_arc && !squared.is_zero()) || compare(&squared, &cap_squared, dec)?.is_gt() {
                break;
            }
            distances.push(squared);
        }
        if distances.len() == 2 {
            alternatives.push(Match {
                other: b,
                reversed,
                gaps: [distances[0].clone(), distances[1].clone()],
                maximum_sample_mm: None,
                check: if exact_arc {
                    "identical-conic-carrier-endpoints-and-finite-traversal"
                } else {
                    "finite-line-segment-endpoint-correspondence"
                },
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
    sew_with(body_index, body, plan, dec, false, false, correspondence)
}

pub(super) fn sew_with(
    body_index: usize,
    body: &mut TolerantBody,
    plan: &mut Plan,
    dec: &mut Decider,
    allow_same_face: bool,
    allow_cross_shells: bool,
    mut correspondence: impl FnMut(&TolerantBody, usize, usize, &mut Decider) -> Result<Option<Match>, Error>,
) -> Result<(), Error> {
    let before = usage(body, dec)?;
    let free: Vec<_> = before
        .iter()
        .enumerate()
        .filter(|(_, uses)| uses.len() == 1)
        .map(|(i, _)| i)
        .collect();
    if free.is_empty() {
        plan.passed(body, &body.name, "sewing.no-free-boundaries");
        return Ok(());
    }
    let mut choices: Vec<Vec<Match>> = vec![Vec::new(); body.edges.len()];
    for (index, &a) in free.iter().enumerate() {
        for &b in &free[index + 1..] {
            dec.charge_work(1)?;
            if (!allow_cross_shells && before[a][0].shell != before[b][0].shell)
                || (!allow_same_face && before[a][0].shell == before[b][0].shell
                    && before[a][0].face == before[b][0].face)
            {
                continue;
            }
            if let Some(found) = correspondence(body, a, b, dec)? {
                choices[a].push(found.clone());
                choices[b].push(Match {
                    other: a,
                    reversed: found.reversed,
                    gaps: found.gaps,
                    check: found.check,
                    maximum_sample_mm: found.maximum_sample_mm,
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
        let mut maximum = gaps[0].clone().max(gaps[1].clone());
        if let Some(sampled) = &found.maximum_sample_mm {
            maximum = maximum.max(sampled.clone());
            plan.measurements.push(HealingMeasurement {
                body: body.name.clone(),
                scope: format!("{}:{}", left.name, right.name),
                check: "17-native-source-curve-sample-pairs-not-a-whole-edge-bound".into(),
                distance_upper_mm: sampled.to_string(),
            });
        }
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
                found.check.into(),
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
    for &edge in &free {
        dec.charge_work(1)?;
        if after[edge].len() == 1 && !matches!(body.edges[edge].curve, Curve::Line(..)) {
            plan.unresolved(
                body,
                &body.edges[edge].name,
                "sewing.nonlinear-finite-correspondence-pending",
            );
        }
    }
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
