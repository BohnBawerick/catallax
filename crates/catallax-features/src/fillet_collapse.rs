//! Exact closed-rim consumption at the first finite-cylinder contact event.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{Circle, Curve, Surface};
use catallax_topology::{BodyId, BodyKind, Edge, EdgeId, FaceId, Loop, LoopId, Model, Sense};

use super::{FilletDefinition, Replacement, edge_blend, reverse, sign};

#[derive(Clone, Copy)]
enum End {
    Regular,
    Cap,
    Wall,
    Both,
}

pub(super) struct Limit {
    pub(super) radius: Real,
    pub(super) region: Region,
}

pub(super) struct Region {
    pub(super) beyond: bool,
    cap: FaceId,
    wall: FaceId,
    cap_index: usize,
    selected: EdgeId,
    opposite: EdgeId,
    cap_sense: Sense,
    wall_sense: Sense,
    base_sense: Sense,
    end: End,
}

pub(super) struct Collapse {
    pub(super) body: BodyId,
    pub(super) consumed: [Option<FaceId>; 2],
}

fn ring(model: &Model, loop_: LoopId) -> Option<(EdgeId, Sense)> {
    let Loop::Coedges(coedges) = model.loop_(loop_) else {
        return None;
    };
    let [coedge] = coedges.as_slice() else {
        return None;
    };
    let coedge = model.coedge(*coedge);
    let edge = model.edge(coedge.edge);
    (edge.start.is_none()
        && edge.end.is_none()
        && matches!(model.geom().curve(edge.curve), Curve::Circle(_)))
    .then_some((coedge.edge, coedge.sense))
}

fn disk(model: &Model, face: FaceId) -> Option<(EdgeId, Sense)> {
    let [loop_] = model.face(face).loops.as_slice() else {
        return None;
    };
    ring(model, *loop_)
}

/// Proves the complete input before choosing its first contact event. For
/// rho<R,H both feet are interior. At R the centre circle contracts; at H
/// the wall foot reaches the opposite cap. Neither permits further contact:
/// closest cylinder distance is at most R, and beyond H the centre leaves
/// the opposite edge's inward normal cone. Its flat cap has no offset fold
/// to trigger the wider-ball cliff rule.
/// An ordinary band also uses this finite proof when its untrimmed carrier
/// reaches the base but its exposed quarter remains above it.
pub(super) fn terminal(
    model: &Model,
    body: BodyId,
    selected: EdgeId,
    sources: [FaceId; 2],
    requested: &Real,
    convex: bool,
    dec: &mut Decider,
) -> Result<Option<Limit>, Refusal> {
    let input = model.body(body);
    if !convex || input.kind != BodyKind::Solid {
        return Ok(None);
    }
    // Other lumps have disjoint interiors. Every affected region lies in
    // this complete cylinder, so those components keep their whole sets.
    let Some(faces) = input.lumps.iter().find_map(|lump| {
        let [shell] = model.lump(*lump).shells.as_slice() else {
            return None;
        };
        let faces = &model.shell(*shell).faces;
        (faces.len() == 3 && sources.iter().all(|face| faces.contains(face))).then_some(faces)
    }) else {
        return Ok(None);
    };
    let supports = sources.map(|face| model.geom().surface(model.face(face).surface));
    let (cap_index, plane, cylinder) = match supports {
        [Surface::Plane(plane), Surface::Cylinder(cylinder)] => (0, plane, cylinder),
        [Surface::Cylinder(cylinder), Surface::Plane(plane)] => (1, plane, cylinder),
        _ => return Ok(None),
    };
    let cap = sources[cap_index];
    let wall = sources[1 - cap_index];
    let Some((cap_edge, cap_sense)) = disk(model, cap) else {
        return Ok(None);
    };
    if cap_edge != selected || model.face(wall).sense != Sense::Same {
        return Ok(None);
    }
    let [first, last] = model.face(wall).loops.as_slice() else {
        return Ok(None);
    };
    let (Some(first), Some(last)) = (ring(model, *first), ring(model, *last)) else {
        return Ok(None);
    };
    let (opposite, wall_sense) = if first.0 == selected {
        (last.0, first.1)
    } else if last.0 == selected {
        (first.0, last.1)
    } else {
        return Ok(None);
    };
    let Some(base) = faces.iter().copied().find(|face| !sources.contains(face)) else {
        return Ok(None);
    };
    let Some((base_edge, base_sense)) = disk(model, base) else {
        return Ok(None);
    };
    if base_edge != opposite {
        return Ok(None);
    }
    let (Surface::Plane(base_plane), Curve::Circle(base_circle), Curve::Circle(cap_circle)) = (
        model.geom().surface(model.face(base).surface),
        model.geom().curve(model.edge(opposite).curve),
        model.geom().curve(model.edge(selected).curve),
    ) else {
        return Ok(None);
    };
    let reversed =
        (model.face(cap).sense != Sense::Same).then(|| plane.normal.scaled(&Real::from_i64(-1)));
    let axis = reversed.as_ref().unwrap_or(&plane.normal);
    let height = cap_circle.center.minus(&base_circle.center);
    if sign(&cylinder.radius_sq, dec)? != Sign::Positive
        || sign(&height.dot(axis), dec)? != Sign::Positive
        || sign(&height.cross(axis).norm_squared(), dec)? != Sign::Zero
        || sign(&cylinder.axis.cross(axis).norm_squared(), dec)? != Sign::Zero
        || sign(&base_plane.normal.cross(axis).norm_squared(), dec)? != Sign::Zero
        || sign(&base_circle.normal.cross(axis).norm_squared(), dec)? != Sign::Zero
        || sign(&cap_circle.normal.cross(axis).norm_squared(), dec)? != Sign::Zero
        || sign(
            &base_circle
                .center
                .minus(&base_plane.origin)
                .dot(&base_plane.normal),
            dec,
        )? != Sign::Zero
        || sign(
            &cap_circle.center.minus(&plane.origin).dot(&plane.normal),
            dec,
        )? != Sign::Zero
        || sign(&(&base_circle.radius_sq - &cylinder.radius_sq), dec)? != Sign::Zero
        || sign(&(&cap_circle.radius_sq - &cylinder.radius_sq), dec)? != Sign::Zero
        || sign(&base_plane.normal.dot(axis), dec)?
            != if model.face(base).sense == Sense::Same {
                Sign::Negative
            } else {
                Sign::Positive
            }
    {
        return Ok(None);
    }
    let height_sq = height.norm_squared();
    let end = match sign(&(&cylinder.radius_sq - &height_sq), dec)? {
        Sign::Negative => End::Cap,
        Sign::Positive => End::Wall,
        Sign::Zero => End::Both,
    };
    let limit_sq = if matches!(end, End::Wall) {
        &height_sq
    } else {
        &cylinder.radius_sq
    };
    let requested_sq = requested.square();
    let comparison = sign(&(&requested_sq - limit_sq), dec)?;
    if comparison == Sign::Negative {
        let diameter_sq = &requested_sq * &Real::from_i64(4);
        if sign(&(&diameter_sq - &height_sq), dec)? == Sign::Negative {
            return Ok(None);
        }
    }
    Ok(Some(Limit {
        radius: if comparison == Sign::Negative {
            requested.clone()
        } else if matches!(end, End::Wall) {
            super::fillet_envelope::real_sqrt(&height_sq)
        } else {
            cylinder.radius()
        },
        region: Region {
            beyond: comparison == Sign::Positive,
            cap,
            wall,
            cap_index,
            selected,
            opposite,
            cap_sense,
            wall_sense,
            base_sense,
            end: if comparison == Sign::Negative {
                End::Regular
            } else {
                end
            },
        },
    }))
}

fn same_circle(
    model: &Model,
    circle: &Circle,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let Curve::Circle(other) = model.geom().curve(model.edge(edge).curve) else {
        return Ok(false);
    };
    Ok(
        sign(&circle.center.minus(&other.center).norm_squared(), dec)? == Sign::Zero
            && sign(&(&circle.radius_sq - &other.radius_sq), dec)? == Sign::Zero
            && sign(&circle.normal.cross(&other.normal).norm_squared(), dec)? == Sign::Zero,
    )
}

fn rail(model: &mut Model, source: EdgeId, contact: &Curve) -> EdgeId {
    let sense = model.edge(source).sense;
    let curve = model.geom_mut().add_curve(contact.clone());
    model.add_edge(Edge {
        curve,
        sense,
        start: None,
        end: None,
    })
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-overflow", detail)
}

pub(super) fn build(
    model: &mut Model,
    body: BodyId,
    region: Region,
    definition: &FilletDefinition,
    dec: &mut Decider,
) -> Result<Collapse, Refusal> {
    let (cap_gone, wall_gone) = match region.end {
        End::Regular => (false, false),
        End::Cap => (true, false),
        End::Wall => (false, true),
        End::Both => (true, true),
    };
    let Curve::Circle(wall_contact) = &definition.contacts[1 - region.cap_index] else {
        return Err(bad("the finite cylinder's wall contact is not circular"));
    };
    let mut edit = Replacement::default();
    let consumed = if cap_gone && wall_gone {
        [Some(region.cap), Some(region.wall)]
    } else if cap_gone {
        [Some(region.cap), None]
    } else if wall_gone {
        [Some(region.wall), None]
    } else {
        [None, None]
    };
    edit.consumed.extend(consumed.into_iter().flatten());
    let loops = if cap_gone {
        if definition.point_contact() != Some(region.cap_index) {
            return Err(bad(
                "cap consumption has no matching point contact on its plane",
            ));
        }
        let (rail, sense) = if wall_gone {
            if !same_circle(model, wall_contact, region.opposite, dec)? {
                return Err(bad(
                    "the terminal wall contact does not coincide with the opposite cap",
                ));
            }
            (region.opposite, reverse(region.base_sense))
        } else {
            let rail = rail(
                model,
                region.selected,
                &definition.contacts[1 - region.cap_index],
            );
            edit.edges.insert((region.wall, region.selected), rail);
            (rail, reverse(region.wall_sense))
        };
        vec![vec![(rail, sense)]]
    } else {
        if !matches!(
            definition.surface,
            Surface::Torus(_) | Surface::TorusSheet(_)
        ) {
            return Err(bad(
                "the finite band has no exact torus or outer-sheet carrier",
            ));
        }
        let Curve::Circle(cap_contact) = &definition.contacts[region.cap_index] else {
            return Err(bad("the finite band's plane contact is not circular"));
        };
        if sign(&cap_contact.radius_sq, dec)? != Sign::Positive {
            return Err(bad("a collapsed plane contact must consume the cap"));
        }
        let cap_rail = rail(
            model,
            region.selected,
            &definition.contacts[region.cap_index],
        );
        edit.edges.insert((region.cap, region.selected), cap_rail);
        let (wall_rail, wall_sense) = if wall_gone {
            if !same_circle(model, wall_contact, region.opposite, dec)? {
                return Err(bad(
                    "the terminal wall contact does not coincide with the opposite cap",
                ));
            }
            (region.opposite, reverse(region.base_sense))
        } else {
            let wall_rail = rail(
                model,
                region.selected,
                &definition.contacts[1 - region.cap_index],
            );
            edit.edges.insert((region.wall, region.selected), wall_rail);
            (wall_rail, reverse(region.wall_sense))
        };
        vec![
            vec![(cap_rail, reverse(region.cap_sense))],
            vec![(wall_rail, wall_sense)],
        ]
    };
    let patch = crate::extrude::face(
        model,
        definition.surface.clone(),
        Sense::Same,
        loops,
        None,
        vec![definition.origin.clone()],
    );
    edit.faces.push((region.cap, patch));
    let result = edge_blend::apply(model, body, &edit)?;
    Ok(Collapse {
        body: result,
        consumed,
    })
}
