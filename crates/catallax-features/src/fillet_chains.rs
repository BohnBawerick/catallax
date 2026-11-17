//! Topological tangent chains, independent of surface-pair construction recipes.
//! Endpoint direction and support-normal gaps are enclosed in millimetres at r.

use super::*;
use catallax_numbers::{Ball, NumericalCheck, NumericalSample};

type Endpoint = (EdgeId, usize);

pub(crate) struct ChainEdge {
    pub edge: EdgeId,
    pub sense: Sense,
    pub faces: [FaceId; 2],
}

pub(crate) struct Chain {
    pub edges: Vec<ChainEdge>,
    pub closed: bool,
    pub checks: Vec<NumericalCheck>,
}

struct EndGeometry {
    tangent: [Ball; 3],
    normals: [[Ball; 3]; 2],
}

struct Joint {
    ends: [Endpoint; 2],
    swap: bool,
    gap: Rational,
    vertex: VertexId,
}

fn refused(detail: impl Into<String>) -> Refusal {
    Refusal::new("fillet.numerical-chain", detail)
}

fn unit(vector: &Vector3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    let values = [
        dec.enclose_budgeted_at(&vector.x, 128)?,
        dec.enclose_budgeted_at(&vector.y, 128)?,
        dec.enclose_budgeted_at(&vector.z, 128)?,
    ];
    let square = values[0].square() + values[1].square() + values[2].square();
    if square.lo().sign() != Sign::Positive {
        return Err(refused(
            "tangent chain has a singular or unresolved endpoint direction",
        ));
    }
    let inverse = square
        .sqrt()
        .and_then(|length| length.recip())
        .ok_or_else(|| refused("tangent chain direction cannot be normalized"))?;
    Ok(values.map(|value| &value * &inverse))
}

fn geometry(
    model: &Model,
    endpoint: Endpoint,
    vertex: VertexId,
    faces: [FaceId; 2],
    dec: &mut Decider,
) -> Result<EndGeometry, Refusal> {
    let edge = model.edge(endpoint.0);
    let point = model.vertex(vertex).point;
    let positive = (endpoint.1 == 0) == (edge.sense == Sense::Same);
    let tangents = catallax_intersect::tangent_at(model.geom(), edge.curve, point, dec)?;
    let mut tangent = unit(&tangents[usize::from(positive)], dec)?;
    if !positive {
        tangent = tangent.map(|value| -value);
    }
    let point = model.geom().position(point);
    let normal = |face: FaceId, dec: &mut Decider| -> Result<[Ball; 3], Refusal> {
        let definition = model.face(face);
        let normal = model
            .geom()
            .surface(definition.surface)
            .normal_at(&point, dec)?;
        let values = unit(&normal, dec)?;
        Ok(if definition.sense == Sense::Same {
            values
        } else {
            values.map(|value| -value)
        })
    };
    Ok(EndGeometry {
        tangent,
        normals: [normal(faces[0], dec)?, normal(faces[1], dec)?],
    })
}

fn gap(delta: [Ball; 3], radius: &Ball) -> Result<Rational, Refusal> {
    let length = (delta[0].square() + delta[1].square() + delta[2].square())
        .sqrt()
        .ok_or_else(|| refused("endpoint continuity distance is unresolved"))?;
    Ok((&length * radius).hi())
}

fn joining(
    a: &EndGeometry,
    b: &EndGeometry,
    radius: &Ball,
    limit: &Rational,
) -> Result<Option<(bool, Rational)>, Refusal> {
    let tangent = gap(
        std::array::from_fn(|i| &a.tangent[i] + &b.tangent[i]),
        radius,
    )?;
    if &tangent > limit {
        return Ok(None);
    }
    let mut best: Option<(bool, Rational)> = None;
    for swap in [false, true] {
        let mut maximum = tangent.clone();
        for side in 0..2 {
            let other = side ^ usize::from(swap);
            maximum = maximum.max(gap(
                std::array::from_fn(|i| &a.normals[side][i] - &b.normals[other][i]),
                radius,
            )?);
        }
        if &maximum <= limit && best.as_ref().is_none_or(|(_, old)| &maximum < old) {
            best = Some((swap, maximum));
        }
    }
    Ok(best)
}

/// Joins only unique tangent continuations; sharp junctions remain strip ends.
/// The caller supplies positive radius and tolerance before construction.
pub(crate) fn group(
    model: &Model,
    body: BodyId,
    selected: &[EdgeId],
    radius: &Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<Vec<Chain>, Refusal> {
    let incident = uses(model, body);
    let mut remaining: BTreeSet<_> = selected.iter().copied().collect();
    if remaining.is_empty() {
        return Err(refused(
            "a fillet requires at least one selected sharp edge",
        ));
    }
    let mut faces = BTreeMap::new();
    let mut vertices: BTreeMap<VertexId, Vec<Endpoint>> = BTreeMap::new();
    for edge in &remaining {
        dec.charge_work(1)?;
        let Some([(a, _), (b, _)]) = incident.get(edge).map(Vec::as_slice) else {
            return Err(refused(
                "selected edge needs exactly two source-body face uses",
            ));
        };
        if a == b {
            return Err(refused("a same-face seam is not a sharp fillet edge"));
        }
        faces.insert(*edge, [*a, *b]);
        let source = model.edge(*edge);
        match [source.start, source.end] {
            [None, None] => {}
            [Some(a), Some(b)] if a == b => {}
            [Some(a), Some(b)] => {
                vertices.entry(a).or_default().push((*edge, 0));
                vertices.entry(b).or_default().push((*edge, 1));
            }
            _ => return Err(refused("an unbounded half-edge has no finite fillet chain")),
        }
    }
    let scale = Ball::from_rational(radius, 128);
    let limit = tolerance.mul_pow2(-3);
    let mut joints = Vec::new();
    let mut links = BTreeMap::new();
    let mut evaluated = Vec::new();
    for (vertex, ends) in vertices {
        if ends.len() < 2 {
            continue;
        }
        evaluated.clear();
        for endpoint in &ends {
            evaluated.push(geometry(model, *endpoint, vertex, faces[&endpoint.0], dec)?);
        }
        for a in 0..ends.len() {
            for b in a + 1..ends.len() {
                dec.charge_work(1)?;
                if let Some((swap, gap)) = joining(&evaluated[a], &evaluated[b], &scale, &limit)? {
                    if links.contains_key(&ends[a]) || links.contains_key(&ends[b]) {
                        return Err(refused(format!(
                            "tangent continuation branches at {vertex:?}"
                        )));
                    }
                    let index = joints.len();
                    joints.push(Joint {
                        ends: [ends[a], ends[b]],
                        swap,
                        gap,
                        vertex,
                    });
                    links.insert(ends[a], index);
                    links.insert(ends[b], index);
                }
            }
        }
    }
    let mut result = Vec::new();
    while !remaining.is_empty() {
        let start = remaining
            .iter()
            .find_map(|edge| {
                (0..2)
                    .find(|end| !links.contains_key(&(*edge, *end)))
                    .map(|end| (*edge, end))
            })
            .unwrap_or_else(|| (*remaining.first().expect("nonempty selection"), 0));
        let mut current = start;
        let mut swapped = false;
        let mut chain = Chain {
            edges: Vec::new(),
            closed: false,
            checks: Vec::new(),
        };
        loop {
            if !remaining.remove(&current.0) {
                return Err(refused("a tangent chain revisits an already consumed edge"));
            }
            let mut sides = faces[&current.0];
            if swapped {
                sides.reverse();
            }
            chain.edges.push(ChainEdge {
                edge: current.0,
                sense: if current.1 == 0 {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
                faces: sides,
            });
            let source = model.edge(current.0);
            if source.start == source.end {
                chain.closed = true;
                break;
            }
            let exit = (current.0, 1 - current.1);
            let Some(index) = links.get(&exit) else { break };
            let joint = &joints[*index];
            let next = if joint.ends[0] == exit {
                joint.ends[1]
            } else {
                joint.ends[0]
            };
            chain.checks.push(NumericalCheck::Samples {
                name: "fillet.chain.endpoint-continuity-at-radius".into(),
                scope: format!("{:?}/{:?}@{:?}", exit.0, next.0, joint.vertex),
                samples: vec![NumericalSample {
                    parameters: vec![radius.clone()],
                    distance_upper_mm: joint.gap.clone(),
                }],
            });
            swapped ^= joint.swap;
            if next == start {
                if swapped {
                    return Err(refused("a closed tangent chain reverses its support sides"));
                }
                chain.closed = true;
                break;
            }
            current = next;
        }
        result.push(chain);
    }
    Ok(result)
}
