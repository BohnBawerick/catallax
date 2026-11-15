//! Straight through openings beside complete circular bores.

use super::*;
use crate::{
    polyhedron::{ConvexCell, PlaneFace},
    profile::sign,
};
use catallax_numbers::Sign;
use catallax_shapes::Plane;

#[derive(Clone)]
pub(super) struct Opening {
    pub cell: ConvexCell,
    wires: [(FaceId, LoopId); 2],
    walls: Vec<FaceId>,
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-arrangement", detail)
}

fn edges(model: &Model, wire: LoopId) -> Result<Vec<EdgeId>, Refusal> {
    let Loop::Coedges(coedges) = model.loop_(wire) else {
        return Err(bad("a polygonal opening needs finite straight wires"));
    };
    coedges
        .iter()
        .map(|coedge| {
            let edge = model.coedge(*coedge).edge;
            let e = model.edge(edge);
            if !matches!(model.geom().curve(e.curve), Curve::Line(_))
                || e.start.is_none()
                || e.end.is_none()
            {
                return Err(bad("a polygonal opening needs line segments"));
            }
            Ok(edge)
        })
        .collect()
}

fn plane(model: &Model, face: FaceId) -> PlaneFace {
    let f = model.face(face);
    let Surface::Plane(source) = model.geom().surface(f.surface) else {
        unreachable!("the source inventory contains only planes")
    };
    let mut plane = source.clone();
    if f.sense == Sense::Reversed {
        plane.normal = plane.normal.scaled(&Real::from_i64(-1));
    }
    PlaneFace {
        plane,
        origins: model.face_origins(face).to_vec(),
    }
}

/// Removes only matched, convex prismatic voids from the planar envelope. Both
/// physical cap wires and every side wall must match the reconstructed cell.
pub(super) fn extract(
    model: &Model,
    planes: &[FaceId],
    caps: [FaceId; 2],
    holes: &mut BTreeSet<LoopId>,
    axis: &Vector3,
    dec: &mut Decider,
) -> Result<(Vec<FaceId>, Vec<Opening>), Refusal> {
    let cap_planes = caps.map(|cap| plane(model, cap));
    for cap in &cap_planes {
        if sign(&cap.plane.normal.cross(axis).norm_squared(), dec)? != Sign::Zero {
            return Err(bad("polygonal through openings need perpendicular caps"));
        }
    }
    let wires = |cap: FaceId| {
        model
            .face(cap)
            .loops
            .iter()
            .skip(1)
            .filter(|wire| !holes.contains(wire))
            .copied()
            .collect::<Vec<_>>()
    };
    let first = wires(caps[0]);
    let second = wires(caps[1]);
    if first.len() != second.len() {
        return Err(bad("polygonal through openings need matching cap wires"));
    }
    let side_faces: Vec<_> = planes
        .iter()
        .filter(|face| !caps.contains(face))
        .map(|face| {
            let wires = &model.face(*face).loops;
            if wires.len() != 1 {
                return Err(bad("a through opening side needs one finite wire"));
            }
            Ok((*face, edges(model, wires[0])?))
        })
        .collect::<Result<_, Refusal>>()?;
    let owners = |wire: LoopId| -> Result<BTreeSet<FaceId>, Refusal> {
        edges(model, wire)?
            .iter()
            .map(|edge| {
                let mut owners = side_faces.iter().filter(|(_, edges)| edges.contains(edge));
                let Some((face, _)) = owners.next() else {
                    return Err(bad("a polygonal opening edge has no side wall"));
                };
                if owners.next().is_some() {
                    return Err(bad("a polygonal opening edge has several side walls"));
                }
                Ok(*face)
            })
            .collect()
    };
    let mut removed = BTreeSet::new();
    let mut matched = BTreeSet::new();
    let mut openings = Vec::new();
    for a in first {
        let walls = owners(a)?;
        if walls.len() < 3 || walls.len() != edges(model, a)?.len() {
            return Err(bad("a polygonal opening needs distinct side walls"));
        }
        let mut candidates = Vec::new();
        for b in &second {
            if owners(*b)? == walls {
                candidates.push(*b);
            }
        }
        let [b] = candidates.as_slice() else {
            return Err(bad("a polygonal opening has no unique opposite cap wire"));
        };
        if !matched.insert(*b) || walls.iter().any(|wall| !removed.insert(*wall)) {
            return Err(bad("polygonal openings share a physical wall"));
        }
        let mut faces = cap_planes.to_vec();
        for wall in &walls {
            let mut face = plane(model, *wall);
            if sign(&face.plane.normal.dot(axis), dec)? != Sign::Zero
                || side_faces.iter().find(|(f, _)| f == wall).unwrap().1.len() != 4
            {
                return Err(bad(
                    "a polygonal opening needs prismatic quadrilateral sides",
                ));
            }
            face.plane.normal = face.plane.normal.scaled(&Real::from_i64(-1));
            faces.push(face);
        }
        let cell = ConvexCell::new(faces, dec)?;
        let mut physical = BTreeSet::new();
        for wire in [a, *b] {
            for edge in edges(model, wire)? {
                let edge = model.edge(edge);
                physical.extend([edge.start.unwrap(), edge.end.unwrap()]);
            }
        }
        if cell.vertices().len() != physical.len() || cell.is_empty() {
            return Err(bad("a polygonal opening is not one convex prism"));
        }
        for vertex in physical {
            let point = model.geom().position(model.vertex(vertex).point);
            let mut found = false;
            for other in cell.vertices() {
                found |= sign(&point.minus(other).norm_squared(), dec)? == Sign::Zero;
            }
            if !found {
                return Err(bad(
                    "a polygonal opening differs from its physical cap wires",
                ));
            }
        }
        holes.extend([a, *b]);
        openings.push(Opening {
            cell,
            wires: [(caps[0], a), (caps[1], *b)],
            walls: walls.into_iter().collect(),
        });
    }
    Ok((
        planes
            .iter()
            .filter(|face| !removed.contains(face))
            .copied()
            .collect(),
        openings,
    ))
}

/// Strict side clearance keeps each unchanged void internal after drafting.
/// Affine plane values attain their maxima on the finite prism's vertices.
pub(super) fn certify(
    planes: impl Iterator<Item = Plane>,
    axis: &Vector3,
    openings: &[Opening],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if openings.is_empty() {
        return Ok(());
    }
    for plane in planes {
        if sign(&plane.normal.cross(axis).norm_squared(), dec)? == Sign::Zero {
            continue;
        }
        for opening in openings {
            for point in opening.cell.vertices() {
                if sign(&point.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Negative {
                    return Err(bad(
                        "a polygonal opening touches or crosses the drafted envelope",
                    ));
                }
            }
        }
    }
    Ok(())
}

/// The caps stay fixed, so their physical hole edges and the connecting walls
/// can be sewn directly into the drafted envelope without a second Boolean.
pub(super) fn append(
    model: &Model,
    opening: &Opening,
    bore: &BoundBore,
    holes: &mut Vec<HolePath>,
    walls: &mut Vec<WallFace>,
) {
    let path = |wire| {
        let Loop::Coedges(coedges) = model.loop_(wire) else {
            unreachable!("the opening inventory checked finite wires")
        };
        coedges
            .iter()
            .map(|coedge| {
                let coedge = model.coedge(*coedge);
                (coedge.edge, coedge.sense)
            })
            .collect()
    };
    for (cap, wire) in opening.wires {
        let support = bore
            .rings
            .iter()
            .find(|(ring, _)| ring.cap == cap)
            .expect("polygonal and circular openings share both caps")
            .1;
        holes.push((support, path(wire)));
    }
    for face in &opening.walls {
        let f = model.face(*face);
        walls.push((
            model.geom().surface(f.surface).clone(),
            f.sense,
            f.loops.iter().map(|wire| path(*wire)).collect(),
            model.face_origins(*face).to_vec(),
        ));
    }
}

pub(super) fn outside(
    cell: &ConvexCell,
    point: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for face in &cell.faces {
        if sign(
            &point.minus(&face.plane.origin).dot(&face.plane.normal),
            dec,
        )? == Sign::Positive
        {
            return Ok(true);
        }
    }
    Ok(false)
}
