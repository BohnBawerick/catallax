//! Feature recognition: the holes, bosses and slots of a body, read off its
//! faces.
//!
//! The method is the rule-based one of the literature and of open recognisers
//! such as quiddity: features are patterns in the faces of the B-rep, their
//! surface kinds, which side of each surface holds material and which faces
//! meet along which edges (the attributed adjacency graph of Joshi and Chang,
//! 1988). The body already states all of that, so nothing is fitted or
//! sampled. Every geometric condition a feature rests on is decided exactly:
//! round faces share one axis and radius, axes are parallel, a slot's flats
//! are tangent to its ends.
//!
//! The rules read planes and circular cylinders. A round is a cylinder face
//! or several of one cylinder that meet along edges. It is a hole's wall
//! when material lies outside it and a boss's when material lies inside it.
//! A round that closes on itself has ends; an end is closed when the face it
//! meets stays within the round's radius and open when that face reaches
//! beyond it.
//!
//! - A hole is a closed round with material outside. With no closed end it
//!   goes through.
//! - Holes of one radius with parallel axes that open onto the same faces
//!   are a pattern.
//! - A boss is a closed round with material inside, one end on a face that
//!   reaches beyond it and the other end closed.
//! - A slot is two open rounds of one radius with parallel axes, material
//!   outside, joined by two parallel flats tangent to both.
//!
//! A face no rule reads is in no feature. Nothing is guessed from it.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Outcome, Real, Refusal, StrictSign};
use catallax_shapes::{Curve, Point3, Surface, Vector3};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

/// The features of one body.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Recognised {
    /// Round holes, each once.
    pub holes: Vec<Hole>,
    /// Groups of equal holes, as indices into `holes`. A hole is in at most one.
    pub hole_patterns: Vec<Vec<usize>>,
    /// Round bosses.
    pub bosses: Vec<Boss>,
    /// Slots with round ends.
    pub slots: Vec<Slot>,
}

/// A round hole.
#[derive(Clone, Debug, PartialEq)]
pub struct Hole {
    /// Its wall: the faces of one cylinder.
    pub faces: Vec<FaceId>,
    /// The nearest double to its diameter in millimetres.
    pub diameter: f64,
    /// True when no end is closed.
    pub through: bool,
    /// The nearest double to the distance between its two ends, when both are
    /// circles about its axis.
    pub depth: Option<f64>,
    /// The faces it opens onto.
    pub opens_onto: Vec<FaceId>,
    /// The faces that close it.
    pub bottoms: Vec<FaceId>,
}

/// A round boss.
#[derive(Clone, Debug, PartialEq)]
pub struct Boss {
    /// Its wall: the faces of one cylinder.
    pub faces: Vec<FaceId>,
    /// The nearest double to its diameter in millimetres.
    pub diameter: f64,
    /// The nearest double to its height, when both ends are circles about its
    /// axis.
    pub height: Option<f64>,
    /// The face it stands on.
    pub base: FaceId,
    /// The face that closes it.
    pub top: FaceId,
}

/// A slot with round ends.
#[derive(Clone, Debug, PartialEq)]
pub struct Slot {
    /// Its two flat sides.
    pub flats: [FaceId; 2],
    /// Its two round ends, each the faces of one cylinder.
    pub ends: [Vec<FaceId>; 2],
    /// The nearest double to its width, the distance between the flats.
    pub width: f64,
    /// The nearest double to its length over the ends.
    pub length: f64,
    /// True when no face floors it.
    pub through: bool,
}

fn is_zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::Refused(refusal) => Err(refusal),
        Outcome::NotExact(inexact) => Err(Refusal::from(inexact)),
    }
}

fn is_positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(sign) => Ok(sign == StrictSign::Positive),
        Outcome::ProvenCoincident(_) => Ok(false),
        Outcome::Refused(refusal) => Err(refusal),
        Outcome::NotExact(inexact) => Err(Refusal::from(inexact)),
    }
}

fn is_zero_vector(value: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(is_zero(&value.x, dec)? && is_zero(&value.y, dec)? && is_zero(&value.z, dec)?)
}

/// A cylinder as the rules read it.
#[derive(Clone)]
struct Round {
    point: Point3,
    direction: Vector3,
    radius_sq: Real,
    /// Material lies outside the cylinder: a hole's wall.
    hollow: bool,
    faces: Vec<FaceId>,
}

impl Round {
    fn parallel(&self, other: &Round, dec: &mut Decider) -> Result<bool, Refusal> {
        is_zero_vector(&self.direction.cross(&other.direction), dec)
    }

    fn same_radius(&self, other: &Round, dec: &mut Decider) -> Result<bool, Refusal> {
        is_zero(&(&self.radius_sq - &other.radius_sq), dec)
    }

    fn coaxial(&self, other: &Round, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.parallel(other, dec)?
            && is_zero_vector(&other.point.minus(&self.point).cross(&self.direction), dec)?)
    }

    /// The squared distance of a point from the axis, times the squared
    /// length of the axis direction, less the same multiple of the squared
    /// radius: positive beyond the round, zero on it, negative within.
    fn beyond(&self, point: &Point3) -> Real {
        let across = point.minus(&self.point).cross(&self.direction);
        &across.norm_squared() - &(&self.radius_sq * &self.direction.norm_squared())
    }

    fn radius(&self) -> f64 {
        self.radius_sq.to_f64().sqrt()
    }

    /// The squared distance along the axis between two points.
    fn span_sq(&self, from: &Point3, to: &Point3) -> Real {
        to.minus(from).dot(&self.direction).square() / self.direction.norm_squared()
    }
}

/// Which faces use each edge, once per use.
fn edge_uses(model: &Model, faces: &[FaceId]) -> BTreeMap<EdgeId, Vec<FaceId>> {
    let mut uses: BTreeMap<EdgeId, Vec<FaceId>> = BTreeMap::new();
    for face in faces {
        for id in &model.face(*face).loops {
            if let Loop::Coedges(coedges) = model.loop_(*id) {
                for coedge in coedges {
                    uses.entry(model.coedge(*coedge).edge)
                        .or_default()
                        .push(*face);
                }
            }
        }
    }
    uses
}

fn face_edges(model: &Model, face: FaceId) -> Vec<EdgeId> {
    let mut edges = Vec::new();
    for id in &model.face(face).loops {
        if let Loop::Coedges(coedges) = model.loop_(*id) {
            edges.extend(coedges.iter().map(|coedge| model.coedge(*coedge).edge));
        }
    }
    edges
}

/// The cylinders of the body: faces of one axis, radius and material side
/// together.
fn rounds(model: &Model, faces: &[FaceId], dec: &mut Decider) -> Result<Vec<Round>, Refusal> {
    let mut rounds: Vec<Round> = Vec::new();
    for face in faces {
        let entity = model.face(*face);
        let Surface::Cylinder(cylinder) = model.geom().surface(entity.surface) else {
            continue;
        };
        let round = Round {
            point: cylinder.origin.clone(),
            direction: cylinder.axis.clone(),
            radius_sq: cylinder.radius_sq.clone(),
            // A cylinder points away from its axis; a face against it has its
            // material on the far side.
            hollow: entity.sense == Sense::Reversed,
            faces: vec![*face],
        };
        let mut placed = false;
        for known in &mut rounds {
            if known.hollow == round.hollow
                && known.same_radius(&round, dec)?
                && known.coaxial(&round, dec)?
            {
                known.faces.push(*face);
                placed = true;
                break;
            }
        }
        if !placed {
            rounds.push(round);
        }
    }
    // Faces of one cylinder are one wall only where they meet: two holes on
    // one axis with material between them are two rounds.
    Ok(rounds
        .into_iter()
        .flat_map(|round| connected(model, round))
        .collect())
}

/// Splits a round into its edge-connected walls.
fn connected(model: &Model, round: Round) -> Vec<Round> {
    let edges: Vec<BTreeSet<EdgeId>> = round
        .faces
        .iter()
        .map(|face| face_edges(model, *face).into_iter().collect())
        .collect();
    let mut component = vec![usize::MAX; round.faces.len()];
    let mut count = 0;
    for seed in 0..round.faces.len() {
        if component[seed] != usize::MAX {
            continue;
        }
        component[seed] = count;
        let mut stack = vec![seed];
        while let Some(at) = stack.pop() {
            for next in 0..round.faces.len() {
                if component[next] == usize::MAX && !edges[at].is_disjoint(&edges[next]) {
                    component[next] = count;
                    stack.push(next);
                }
            }
        }
        count += 1;
    }
    (0..count)
        .map(|part| Round {
            faces: round
                .faces
                .iter()
                .zip(&component)
                .filter(|(_, of)| **of == part)
                .map(|(face, _)| *face)
                .collect(),
            ..round.clone()
        })
        .collect()
}

/// True when an edge is a straight line along the round's axis.
fn is_ruling(
    model: &Model,
    round: &Round,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    match model.geom().curve(model.edge(edge).curve) {
        Curve::Line(line) => is_zero_vector(&line.dir.cross(&round.direction), dec),
        _ => Ok(false),
    }
}

/// Whether a face reaches beyond the round, judged on its edges other than
/// those it shares with the round. `None` when none of them says.
fn reaches_beyond(
    model: &Model,
    round: &Round,
    face: FaceId,
    shared: &BTreeSet<EdgeId>,
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let mut known = false;
    for edge in face_edges(model, face) {
        if shared.contains(&edge) {
            continue;
        }
        let entity = model.edge(edge);
        let ends: Vec<Point3> = [entity.start, entity.end]
            .into_iter()
            .flatten()
            .map(|vertex| model.geom().position(model.vertex(vertex).point))
            .collect();
        if ends.is_empty() {
            // A whole circle about the round's axis is beyond it by its radius.
            if let Curve::Circle(circle) = model.geom().curve(entity.curve)
                && is_zero_vector(&circle.normal.cross(&round.direction), dec)?
                && is_zero_vector(
                    &circle.center.minus(&round.point).cross(&round.direction),
                    dec,
                )?
            {
                known = true;
                if is_positive(&(&circle.radius_sq - &round.radius_sq), dec)? {
                    return Ok(Some(true));
                }
            }
            continue;
        }
        known = true;
        for point in &ends {
            if is_positive(&round.beyond(point), dec)? {
                return Ok(Some(true));
            }
        }
    }
    // A face bounded by the round alone, such as a hole's flat bottom, stays
    // within it.
    Ok(
        if known
            || face_edges(model, face)
                .iter()
                .all(|edge| shared.contains(edge))
        {
            Some(false)
        } else {
            None
        },
    )
}

/// One end of a closed round: the face it meets and whether that face closes
/// it.
struct End {
    face: FaceId,
    closed: bool,
    /// The centre of the end, when it is a circle about the axis.
    at: Option<Point3>,
}

/// The ends of a round, or `None` when a side of it is open, as a slot's end
/// is.
fn ends(
    model: &Model,
    round: &Round,
    uses: &BTreeMap<EdgeId, Vec<FaceId>>,
    dec: &mut Decider,
) -> Result<Option<Vec<End>>, Refusal> {
    let mut by_face: BTreeMap<FaceId, BTreeSet<EdgeId>> = BTreeMap::new();
    let mut own: BTreeSet<EdgeId> = BTreeSet::new();
    for face in &round.faces {
        own.extend(face_edges(model, *face));
    }
    for edge in &own {
        let users = uses.get(edge).map_or(&[][..], Vec::as_slice);
        let outside: Vec<FaceId> = users
            .iter()
            .copied()
            .filter(|face| !round.faces.contains(face))
            .collect();
        if outside.is_empty() {
            // A seam of the round itself.
            continue;
        }
        if is_ruling(model, round, *edge, dec)? {
            return Ok(None);
        }
        for face in outside {
            by_face.entry(face).or_default().insert(*edge);
        }
    }
    let mut found = Vec::new();
    for (face, shared) in by_face {
        let Some(beyond) = reaches_beyond(model, round, face, &shared, dec)? else {
            return Ok(None);
        };
        let mut at = None;
        for edge in &shared {
            if let Curve::Circle(circle) = model.geom().curve(model.edge(*edge).curve)
                && is_zero_vector(&circle.normal.cross(&round.direction), dec)?
            {
                at = Some(circle.center.clone());
            }
        }
        found.push(End {
            face,
            closed: !beyond,
            at,
        });
    }
    Ok(Some(found))
}

/// The squared distance between the two ends of a round, when it has two and
/// both are circles about its axis.
fn span_sq(round: &Round, ends: &[End]) -> Option<Real> {
    match ends {
        [first, second] => Some(round.span_sq(first.at.as_ref()?, second.at.as_ref()?)),
        _ => None,
    }
}

fn root(value: &Option<Real>) -> Option<f64> {
    value.as_ref().map(|value| value.to_f64().sqrt())
}

/// The plane of a face: a point and a normal.
fn plane(model: &Model, face: FaceId) -> Option<(Point3, Vector3)> {
    match model.geom().surface(model.face(face).surface) {
        Surface::Plane(plane) => Some((plane.origin.clone(), plane.normal.clone())),
        _ => None,
    }
}

/// The flats that join two open rounds into a slot: each meets both rounds
/// along a ruling, they are parallel, and they lie a diameter apart.
fn slot_flats(
    model: &Model,
    first: &Round,
    second: &Round,
    uses: &BTreeMap<EdgeId, Vec<FaceId>>,
    dec: &mut Decider,
) -> Result<Option<[FaceId; 2]>, Refusal> {
    let mut beside = |round: &Round| -> Result<BTreeSet<FaceId>, Refusal> {
        let mut faces = BTreeSet::new();
        for face in &round.faces {
            for edge in face_edges(model, *face) {
                if !is_ruling(model, round, edge, dec)? {
                    continue;
                }
                for user in uses.get(&edge).map_or(&[][..], Vec::as_slice) {
                    if !round.faces.contains(user) {
                        faces.insert(*user);
                    }
                }
            }
        }
        Ok(faces)
    };
    let (near, far) = (beside(first)?, beside(second)?);
    let shared: Vec<FaceId> = near.intersection(&far).copied().collect();
    let [low, high] = shared[..] else {
        return Ok(None);
    };
    let (Some((low_point, low_normal)), Some((high_point, high_normal))) =
        (plane(model, low), plane(model, high))
    else {
        return Ok(None);
    };
    if !is_zero_vector(&low_normal.cross(&high_normal), dec)?
        || !is_zero(&low_normal.dot(&first.direction), dec)?
    {
        return Ok(None);
    }
    // The flats are tangent to the ends: a diameter apart.
    let gap = high_point.minus(&low_point).dot(&low_normal);
    let four = Real::from_i64(4);
    let tangent = &gap.square() - &(&(&four * &first.radius_sq) * &low_normal.norm_squared());
    Ok(is_zero(&tangent, dec)?.then_some([low, high]))
}

/// Recognises the holes, hole patterns, bosses and slots of a body.
///
/// The body is read as it stands; nothing is built or changed. An exact
/// decision that cannot be made within the decider's allowance refuses with
/// the decider's own refusal.
pub fn recognise(model: &Model, body: BodyId, dec: &mut Decider) -> Result<Recognised, Refusal> {
    let faces = model.faces_of_body(body);
    let uses = edge_uses(model, &faces);
    let rounds = rounds(model, &faces, dec)?;
    let mut out = Recognised::default();
    let mut open: Vec<&Round> = Vec::new();
    let mut hole_rounds: Vec<(&Round, Option<Real>)> = Vec::new();
    for round in &rounds {
        let Some(ends) = ends(model, round, &uses, dec)? else {
            if round.hollow {
                open.push(round);
            }
            continue;
        };
        let closed: Vec<&End> = ends.iter().filter(|end| end.closed).collect();
        let reaching: Vec<&End> = ends.iter().filter(|end| !end.closed).collect();
        if round.hollow {
            if reaching.is_empty() {
                // A cavity with no way in is not a hole.
                continue;
            }
            let depth_sq = span_sq(round, &ends);
            out.holes.push(Hole {
                faces: round.faces.clone(),
                diameter: 2.0 * round.radius(),
                through: closed.is_empty(),
                depth: root(&depth_sq),
                opens_onto: reaching.iter().map(|end| end.face).collect(),
                bottoms: closed.iter().map(|end| end.face).collect(),
            });
            hole_rounds.push((round, depth_sq));
        } else if let ([top], [base]) = (&closed[..], &reaching[..]) {
            out.bosses.push(Boss {
                faces: round.faces.clone(),
                diameter: 2.0 * round.radius(),
                height: root(&span_sq(round, &ends)),
                base: base.face,
                top: top.face,
            });
        }
    }

    // Equal holes with parallel axes that open onto the same faces and end
    // the same way.
    let mut grouped = vec![false; out.holes.len()];
    for first in 0..out.holes.len() {
        if grouped[first] {
            continue;
        }
        let mut members = vec![first];
        for other in first + 1..out.holes.len() {
            if grouped[other] {
                continue;
            }
            let (a, b) = (&out.holes[first], &out.holes[other]);
            let ((near, near_depth), (far, far_depth)) = (&hole_rounds[first], &hole_rounds[other]);
            let same_depth = match (near_depth, far_depth) {
                (Some(near), Some(far)) => is_zero(&(near - far), dec)?,
                (None, None) => true,
                _ => false,
            };
            let alike = a.through == b.through
                && a.opens_onto == b.opens_onto
                && same_depth
                && near.same_radius(far, dec)?
                && near.parallel(far, dec)?;
            if alike {
                members.push(other);
            }
        }
        if members.len() > 1 {
            for member in &members {
                grouped[*member] = true;
            }
            out.hole_patterns.push(members);
        }
    }

    // Two open rounds joined by two tangent flats.
    let mut used = vec![false; open.len()];
    for first in 0..open.len() {
        for second in first + 1..open.len() {
            if used[first] || used[second] {
                continue;
            }
            let (near, far) = (open[first], open[second]);
            if !near.same_radius(far, dec)?
                || !near.parallel(far, dec)?
                || near.coaxial(far, dec)?
            {
                continue;
            }
            let Some(flats) = slot_flats(model, near, far, &uses, dec)? else {
                continue;
            };
            used[first] = true;
            used[second] = true;
            let across = far.point.minus(&near.point).cross(&near.direction);
            let apart =
                (across.norm_squared().to_f64() / near.direction.norm_squared().to_f64()).sqrt();
            // A floor meets both ends and both flats and has no edge of its
            // own; a face the slot opens onto has.
            let own: BTreeSet<FaceId> = near
                .faces
                .iter()
                .chain(&far.faces)
                .chain(&flats)
                .copied()
                .collect();
            let mut walls: BTreeSet<EdgeId> = BTreeSet::new();
            let mut around: Option<BTreeSet<FaceId>> = None;
            for face in &own {
                let edges = face_edges(model, *face);
                let beside: BTreeSet<FaceId> = edges
                    .iter()
                    .flat_map(|edge| uses.get(edge).cloned().unwrap_or_default())
                    .collect();
                walls.extend(edges);
                around = Some(match around {
                    None => beside,
                    Some(known) => known.intersection(&beside).copied().collect(),
                });
            }
            let floored = around.unwrap_or_default().into_iter().any(|face| {
                !own.contains(&face)
                    && face_edges(model, face)
                        .iter()
                        .all(|edge| walls.contains(edge))
            });
            out.slots.push(Slot {
                flats,
                ends: [near.faces.clone(), far.faces.clone()],
                width: 2.0 * near.radius(),
                length: apart + 2.0 * near.radius(),
                through: !floored,
            });
        }
    }
    Ok(out)
}
