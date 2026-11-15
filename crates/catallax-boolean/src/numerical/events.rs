//! Pair native boundary occurrences by ownership and ordered roots.

use super::*;
use catallax_numbers::NumericalRepair;
use catallax_topology::numerical::NumericalVertex;
use std::collections::BTreeMap;

pub(super) struct Events {
    pub vertices: Vec<NumericalVertex>,
    pub uv: Vec<[Rational; 2]>,
    pub ends: BTreeMap<(usize, usize), usize>,
    pub records: Vec<NumericalFiniteEvent>,
    pub repairs: Vec<NumericalRepair>,
    pub checks: Vec<NumericalCheck>,
}

#[derive(Clone)]
struct Occurrence {
    hit: NumericalFiniteHit,
    end: Option<usize>,
    xyz: [Rational; 3],
    domain: [usize; 2],
}

pub(super) fn rational_point(point: Point3) -> Result<[Rational; 3], Refusal> {
    [&point.x, &point.y, &point.z]
        .map(|q| {
            q.as_rational()
                .ok_or_else(|| bad("retained fitted coordinate lost rational identity"))
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| unreachable!())
}
pub(super) fn box_point(point: &[Rational; 3]) -> [Ball; 3] {
    point.each_ref().map(|q| Ball::from_rational(q, BITS))
}
pub(super) fn point_piece(piece: &catallax_intersect::NumericalIntersectionPiece) -> bool {
    let first = piece.curve().control(0);
    piece.curve().controls().all(|p| {
        p.x.same_definition(&first.x)
            && p.y.same_definition(&first.y)
            && p.z.same_definition(&first.z)
    })
}
fn quarter(domain: &[[Rational; 2]; 2]) -> Result<usize, Refusal> {
    (0..4)
        .find(|i| domain[0] == [q(*i as i64), q(*i as i64 + 1)])
        .ok_or_else(|| bad("finite fit lost its native quarter ownership"))
}
fn normalized(mut uv: [Rational; 2]) -> [Rational; 2] {
    if uv[0] == q(4) {
        uv[0] = q(0)
    }
    uv
}
fn signature(hit: &NumericalFiniteHit) -> [Option<Rational>; 4] {
    std::array::from_fn(|i| {
        let source = i / 2;
        let axis = i % 2;
        let value = &hit.uv[source][axis];
        if hit.requested_domains[source][axis].contains(value) {
            Some(if axis == 0 && value == &q(4) {
                q(0)
            } else {
                value.clone()
            })
        } else {
            None
        }
    })
}
fn order(a: &Occurrence, b: &Occurrence) -> std::cmp::Ordering {
    normalized(a.hit.uv[0].clone())
        .cmp(&normalized(b.hit.uv[0].clone()))
        .then_with(|| normalized(a.hit.uv[1].clone()).cmp(&normalized(b.hit.uv[1].clone())))
        .then_with(|| a.hit.piece.cmp(&b.hit.piece))
}

pub(super) fn collect(
    key: &str,
    intersection: &NumericalNativeIntersection,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Events, Refusal> {
    let source = &intersection.windows().sources()[0];
    let mut result = Events {
        vertices: Vec::new(),
        uv: Vec::new(),
        ends: BTreeMap::new(),
        records: Vec::new(),
        repairs: Vec::new(),
        checks: Vec::new(),
    };
    for u in 0..4 {
        for v in 0..2 {
            let uv = [q(u), q(v)];
            let point = source.evaluate_box([&uv[0], &uv[0]], [&uv[1], &uv[1]], BITS, dec)?;
            result.vertices.push(NumericalVertex {
                key: format!("{key}/corner-{u}-{v}"),
                position: point.each_ref().map(Ball::mid),
            });
            result.uv.push(uv);
        }
    }
    let mut groups = BTreeMap::<String, Vec<Occurrence>>::new();
    for (index, piece) in intersection.pieces().iter().enumerate() {
        let point_only = point_piece(piece);
        for end in 0..if point_only { 1 } else { 2 } {
            dec.charge_work(1)?;
            let parameter = piece.curve().domain()[end].clone();
            let at = Real::from(parameter.clone());
            let a = rational_point(piece.uv()[0].point_at(&at, dec)?)?;
            let b = rational_point(piece.uv()[1].point_at(&at, dec)?)?;
            let hit = NumericalFiniteHit {
                piece: index,
                parameter,
                uv: [[a[0].clone(), a[1].clone()], [b[0].clone(), b[1].clone()]],
                requested_domains: piece.source_domains().clone(),
            };
            let boundary = signature(&hit);
            let domain = [
                quarter(&hit.requested_domains[0])?,
                quarter(&hit.requested_domains[1])?,
            ];
            let group = if boundary.iter().all(Option::is_none) {
                format!("interior-piece-{index}")
            } else {
                let angular_owners = [
                    boundary[0].is_none().then_some(domain[0]),
                    boundary[2].is_none().then_some(domain[1]),
                ];
                format!("boundary-{boundary:?}-angular-owners-{angular_owners:?}")
            };
            groups.entry(group).or_default().push(Occurrence {
                hit,
                end: (!point_only).then_some(end),
                xyz: rational_point(piece.curve().point_at(&at, dec)?)?,
                domain,
            });
        }
    }
    for (group, occurrences) in groups {
        let mut arcs = BTreeMap::<[usize; 2], Vec<Occurrence>>::new();
        let mut points = BTreeMap::<[usize; 2], Vec<Occurrence>>::new();
        for occurrence in occurrences {
            if occurrence.end.is_some() {
                arcs.entry(occurrence.domain).or_default().push(occurrence)
            } else {
                points
                    .entry(occurrence.domain)
                    .or_default()
                    .push(occurrence)
            }
        }
        if group.starts_with("interior-piece-") {
            let all: Vec<_> = arcs
                .into_values()
                .chain(points.into_values())
                .flatten()
                .collect();
            add_event(&mut result, key, &group, 0, all, tolerance, dec)?;
            continue;
        }
        let anchor = if arcs.is_empty() { &points } else { &arcs };
        let count = anchor.values().map(Vec::len).max().unwrap_or(0);
        let mut roots = vec![Vec::<Occurrence>::new(); count];
        let mut anchors = Vec::new();
        for bucket in anchor.values() {
            if bucket.len() != count {
                return Err(bad(
                    "adjacent finite boundary root counts disagree; no near-hit merge was applied",
                ));
            }
            let mut sorted = bucket.clone();
            sorted.sort_by(order);
            if anchors.is_empty() {
                anchors = sorted.clone()
            }
            for (root, occurrence) in roots.iter_mut().zip(sorted) {
                root.push(occurrence)
            }
        }
        if !arcs.is_empty() {
            for point in points.into_values().flatten() {
                let mut distances = anchors
                    .iter()
                    .enumerate()
                    .map(|(i, a)| Ok((distance(&box_point(&point.xyz), &box_point(&a.xyz))?, i)))
                    .collect::<Result<Vec<_>, Refusal>>()?;
                distances.sort();
                let Some((gap, index)) = distances.first() else {
                    return Err(bad("contact has no retained finite event"));
                };
                if gap > tolerance {
                    return Err(bad("point report does not match an owned finite arc event"));
                }
                roots[*index].push(point);
            }
        }
        for (ordinal, root) in roots.into_iter().enumerate() {
            add_event(&mut result, key, &group, ordinal, root, tolerance, dec)?
        }
    }
    Ok(result)
}

fn add_event(
    result: &mut Events,
    key: &str,
    group: &str,
    ordinal: usize,
    occurrences: Vec<Occurrence>,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let first = occurrences
        .first()
        .ok_or_else(|| bad("empty finite event"))?;
    let uv = normalized(first.hit.uv[0].clone());
    let point_only = occurrences.iter().all(|o| o.end.is_none());
    let domain = &first.hit.requested_domains[0];
    let on_cone_boundary = (0..2).any(|i| domain[i].contains(&first.hit.uv[0][i]));
    let name = format!("{key}/event/{group}/{ordinal}");
    let vertex = if !point_only || on_cone_boundary {
        if let Some(i) = result.uv.iter().position(|p| p == &uv) {
            Some(i)
        } else {
            let i = result.vertices.len();
            result.vertices.push(NumericalVertex {
                key: name.clone(),
                position: first.xyz.clone(),
            });
            result.uv.push(uv.clone());
            Some(i)
        }
    } else {
        None
    };
    let after = vertex.map_or(&first.xyz, |i| &result.vertices[i].position);
    let mut samples = Vec::new();
    for occurrence in &occurrences {
        dec.charge_work(1)?;
        let gap = distance(&box_point(&occurrence.xyz), &box_point(after))?;
        if gap > tolerance.clone() {
            return Err(bad(
                "owned finite endpoint join exceeds the preset physical tolerance",
            ));
        }
        samples.push(NumericalSample {
            parameters: vec![occurrence.hit.parameter.clone()],
            distance_upper_mm: gap.clone(),
        });
        if occurrence.xyz != *after {
            result.repairs.push(NumericalRepair::VertexMoved {
                vertex: format!(
                    "{name}/piece-{}/end-{:?}",
                    occurrence.hit.piece, occurrence.end
                ),
                before: occurrence.xyz.clone(),
                after: after.clone(),
                distance_upper_mm: gap,
            });
        }
        if let Some(end) = occurrence.end {
            result.ends.insert(
                (occurrence.hit.piece, end),
                vertex.ok_or_else(|| bad("fitted arc has no finite vertex"))?,
            );
        }
    }
    result.checks.push(NumericalCheck::Samples {
        name: "finite.explicit-owned-endpoint-join".into(),
        scope: name.clone(),
        samples,
    });
    result.records.push(NumericalFiniteEvent {
        key: name,
        vertex,
        point_only,
        hits: occurrences.into_iter().map(|o| o.hit).collect(),
    });
    Ok(())
}
