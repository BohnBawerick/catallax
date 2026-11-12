//! Sew an already-proved convex inset to unchanged complete planar openings.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::Surface;
use catallax_topology::{
    Body, BodyId, BodyKind, EdgeId, FaceId, Loop, Lump, Model, Sense, Shell, VertexId,
};

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.analytic-rim-domain", detail)
}

pub(crate) fn outward(model: &Model, face: FaceId) -> Option<Surface> {
    let face = model.face(face);
    let surface = model.geom().surface(face.surface);
    if face.sense == Sense::Same {
        return Some(surface.clone());
    }
    if let Surface::Plane(plane) = surface {
        let mut plane = plane.clone();
        plane.normal = plane.normal.scaled(&Real::from_i64(-1));
        return Some(Surface::Plane(plane));
    }
    None
}

struct FaceUse {
    surface: Surface,
    sense: Sense,
    paths: Vec<Vec<(EdgeId, Sense)>>,
    apex: Option<VertexId>,
    origins: Vec<Origin>,
}
impl FaceUse {
    fn read(model: &Model, face: FaceId, reverse: bool) -> Result<Self, Refusal> {
        let stored = model.face(face);
        let mut paths = Vec::new();
        let mut apex = None;
        for lp in &stored.loops {
            match model.loop_(*lp) {
                Loop::Coedges(coedges) => {
                    let mut path: Vec<_> = coedges
                        .iter()
                        .map(|coedge| {
                            let coedge = model.coedge(*coedge);
                            (
                                coedge.edge,
                                if reverse {
                                    crate::extrude::reverse(coedge.sense)
                                } else {
                                    coedge.sense
                                },
                            )
                        })
                        .collect();
                    if reverse {
                        path.reverse();
                    }
                    paths.push(path);
                }
                Loop::Vertex(vertex) if apex.is_none() => {
                    apex = Some(*vertex);
                }
                _ => {
                    return Err(bad(
                        "analytic inset sewing needs complete loops and at most one cone tip",
                    ));
                }
            }
        }
        Ok(Self {
            surface: model.geom().surface(stored.surface).clone(),
            sense: if reverse {
                crate::extrude::reverse(stored.sense)
            } else {
                stored.sense
            },
            paths,
            apex,
            origins: model.face_origins(face).to_vec(),
        })
    }
    fn build(self, model: &mut Model) -> FaceId {
        crate::extrude::face(
            model,
            self.surface,
            self.sense,
            self.paths,
            self.apex,
            self.origins,
        )
    }
}

fn matching_planes(
    model: &Model,
    ids: &[FaceId],
    opening: &Surface,
    dec: &mut Decider,
) -> Result<Vec<usize>, Refusal> {
    let mut matches = Vec::new();
    let Surface::Plane(expected) = opening else {
        return Err(bad("an opening must be planar"));
    };
    for (i, face) in ids.iter().enumerate() {
        let Some(Surface::Plane(actual)) = outward(model, *face) else {
            continue;
        };
        if crate::profile::sign(&actual.normal.cross(&expected.normal).norm_squared(), dec)?
            == Sign::Zero
            && crate::profile::sign(
                &actual.origin.minus(&expected.origin).dot(&expected.normal),
                dec,
            )? == Sign::Zero
            && crate::profile::sign(&actual.normal.dot(&expected.normal), dec)? == Sign::Positive
        {
            matches.push(i);
        }
    }
    Ok(matches)
}

// The caller proves the cavity's distance definition and complete convex section
// inventory. Each opened support stays fixed and cuts two nested closed disks.
// Its outer disk minus its inset disk is the exact planar rim. All other cavity
// faces bound material on their opposite side, with the same native meeting edges.
pub(super) fn sew(
    model: &mut Model,
    body: BodyId,
    cavity: BodyId,
    openings: &[Surface],
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let outer_ids = model.faces_of_body(body);
    let inner_ids = model.faces_of_body(cavity);
    let mut outer = outer_ids
        .iter()
        .map(|face| FaceUse::read(model, *face, false))
        .collect::<Result<Vec<_>, _>>()?;
    let mut omitted = vec![false; inner_ids.len()];
    for opening in openings {
        if !matches!(opening, Surface::Plane(_)) {
            return Err(bad(
                "only a proved complete planar opening can be sewn directly",
            ));
        }
        let outer_matches = matching_planes(model, &outer_ids, opening, dec)?;
        let inner_matches = matching_planes(model, &inner_ids, opening, dec)?;
        let ([outside], [inside]) = (outer_matches.as_slice(), inner_matches.as_slice()) else {
            return Err(bad(
                "an unchanged opening must resolve to exactly one outer and inset disk",
            ));
        };
        if omitted[*inside] || outer[*outside].paths.len() != 1 || outer[*outside].apex.is_some() {
            return Err(bad("an analytic rim requires one complete outer disk"));
        }
        let rim = FaceUse::read(model, inner_ids[*inside], true)?;
        if rim.paths.len() != 1 || rim.apex.is_some() {
            return Err(bad("an analytic rim requires one complete inset disk"));
        }
        outer[*outside].paths.extend(rim.paths);
        omitted[*inside] = true;
    }
    let mut inner = Vec::new();
    for (face, omit) in inner_ids.iter().zip(omitted) {
        if !omit {
            inner.push(FaceUse::read(model, *face, true)?.build(model));
        }
    }
    let mut outer: Vec<_> = outer.into_iter().map(|face| face.build(model)).collect();
    let shells = if openings.is_empty() {
        vec![
            model.add_shell(Shell { faces: outer }),
            model.add_shell(Shell { faces: inner }),
        ]
    } else {
        outer.extend(inner);
        vec![model.add_shell(Shell { faces: outer })]
    };
    let lump = model.add_lump(Lump { shells });
    Ok(model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: model.body(body).marks.clone(),
    }))
}
