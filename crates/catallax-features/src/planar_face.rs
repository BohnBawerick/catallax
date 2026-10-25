//! Straight normal sweeps of one selected planar face with a simple, possibly nonconvex, boundary.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, FaceRef, Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Curve, Line, Plane, Point3, PointDef, Surface, Vector3};
use catallax_topology::{
    Body, BodyId, BodyKind, Edge, EdgeId, FaceId, Loop, Lump, Model, Sense, Shell, Vertex, VertexId,
};

use crate::{
    extrude::{face, reverse},
    profile::sign,
};

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.planar-domain", detail)
}

/// One boundary edge in traversal order, with the face across it.
struct Rim {
    from: Point3,
    to: Point3,
    across: FaceId,
}

fn across(model: &Model, coedge: catallax_topology::CoedgeId, own: FaceId) -> Option<FaceId> {
    let edge = model.coedge(coedge).edge;
    model
        .coedges_of_edge(edge)
        .into_iter()
        .filter(|c| *c != coedge)
        .find_map(|c| model.owner_of_coedge(c).map(|(_, f)| f))
        .filter(|f| *f != own)
}

fn neighbours(model: &Model, face: FaceId) -> Vec<FaceId> {
    let mut result = Vec::new();
    for lp in &model.face(face).loops {
        if let Loop::Coedges(coedges) = model.loop_(*lp) {
            for c in coedges {
                result.extend(across(model, *c, face));
            }
        }
    }
    result
}

/// The single face every reference resolves to, or none when the selection is not one face.
fn select(
    model: &Model,
    body: BodyId,
    refs: &[FaceRef],
    dec: &mut Decider,
) -> Result<Option<FaceId>, Refusal> {
    let mut chosen = None;
    for reference in refs {
        let mut found: Vec<_> = model
            .faces_of_body(body)
            .into_iter()
            .filter(|f| {
                model
                    .face_origins(*f)
                    .iter()
                    .any(|o| reference.face.matches(o))
            })
            .filter(|f| {
                reference.touching.iter().all(|pattern| {
                    neighbours(model, *f)
                        .into_iter()
                        .any(|n| model.face_origins(n).iter().any(|o| pattern.matches(o)))
                })
            })
            .collect();
        if found.len() > 1 && reference.near.is_some() {
            use crate::draft::arrangement::{Patch, Skin};
            let mut skin = Skin {
                vertices: Vec::new(),
                patches: Vec::new(),
            };
            for f in &found {
                let Surface::Plane(plane) = model.geom().surface(model.face(*f).surface) else {
                    return Ok(None);
                };
                let loops = boundary(model, *f)?;
                let mut cycles = Vec::new();
                for rims in loops {
                    let first = skin.vertices.len();
                    skin.vertices.extend(rims.into_iter().map(|r| r.from));
                    cycles.push((first..skin.vertices.len()).collect());
                }
                skin.patches.push(Patch {
                    plane: plane.clone(),
                    loops: cycles,
                    origins: model.face_origins(*f).to_vec(),
                });
            }
            let mut pick = reference.clone();
            // Adjacency has already filtered the physical faces above.
            pick.touching.clear();
            let selected = skin.resolve_faces(&[pick], dec).map_err(|e| match e.code {
                "recipe.ambiguous-reference" => Refusal::new("reference.ambiguous-face", e.detail),
                _ => e,
            })?;
            found = selected.into_iter().map(|i| found[i]).collect();
        }
        match found.as_slice() {
            [] => {
                return Err(Refusal::new(
                    "reference.lost-face",
                    format!("no input face resolves {reference:?}"),
                ));
            }
            [one] if chosen.is_none_or(|c| c == *one) => chosen = Some(*one),
            _ => return Ok(None),
        }
    }
    Ok(chosen)
}

/// Boundary loops of a planar face, each in counter-clockwise order about `normal`.
fn boundary(model: &Model, face: FaceId) -> Result<Vec<Vec<Rim>>, Refusal> {
    let mut loops = Vec::new();
    for lp in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Err(bad("a swept planar face needs a polygonal boundary"));
        };
        let mut rims = Vec::with_capacity(coedges.len());
        for c in coedges {
            let coedge = model.coedge(*c);
            let edge = model.edge(coedge.edge);
            if !matches!(model.geom().curve(edge.curve), Curve::Line(_)) {
                return Err(bad("a swept planar face needs straight rims"));
            }
            let (Some(start), Some(end)) = (edge.start, edge.end) else {
                return Err(bad("a polygon edge requires finite endpoints"));
            };
            let (from, to) = match coedge.sense {
                Sense::Same => (start, end),
                Sense::Reversed => (end, start),
            };
            let Some(neighbour) = across(model, *c, face) else {
                return Err(bad(
                    "a swept planar face needs a distinct face across each rim",
                ));
            };
            rims.push(Rim {
                from: model.geom().position(model.vertex(from).point),
                to: model.geom().position(model.vertex(to).point),
                across: neighbour,
            });
        }
        loops.push(rims);
    }
    Ok(loops)
}

fn rim_origins(
    model: &Model,
    face: FaceId,
    across: FaceId,
    step: &StepId,
    complete: bool,
) -> Result<Vec<Origin>, Refusal> {
    let mut names = Vec::new();
    let source = model.face_origins(face);
    let neighbour = model.face_origins(across);
    let source = if complete {
        source
    } else {
        &source[..source.len().min(1)]
    };
    let neighbour = if complete {
        neighbour
    } else {
        &neighbour[..neighbour.len().min(1)]
    };
    for a in source {
        for b in neighbour {
            let edge = EdgeName::new(a.clone(), b.clone()).map_err(|e| bad(&e.to_string()))?;
            names.push(Origin {
                step: step.clone(),
                role: Role::Rim { edge },
            });
        }
    }
    names.sort();
    names.dedup();
    if names.is_empty() {
        return Err(bad("a normal-sweep rim requires named input faces"));
    }
    Ok(names)
}

/// Whether the body is one profile extrusion and `face` is one of its end caps.
/// The field evaluator sweeps exactly these caps from its extrusion definition.
fn plain_prism(model: &Model, body: BodyId, face: FaceId) -> bool {
    let faces = model.faces_of_body(body);
    let Some(step) = faces
        .first()
        .and_then(|f| model.face_origins(*f).first())
        .map(|o| &o.step)
    else {
        return false;
    };
    faces.iter().all(|f| {
        let origins = model.face_origins(*f);
        !origins.is_empty()
            && origins.iter().all(|o| {
                &o.step == step && matches!(o.role, Role::Start | Role::End | Role::Side { .. })
            })
    }) && model
        .face_origins(face)
        .iter()
        .all(|o| matches!(o.role, Role::Start | Role::End))
}

/// Bounded additional family: polygonal faces retained from several elementary
/// planar creators. A plain authored nonconvex side/multi-face skin stays outside
/// this fallback, as do curved neighbours and derived feature envelopes.
fn planar_boolean_skin(model: &Model, body: BodyId) -> bool {
    let mut creators = std::collections::BTreeSet::new();
    model.faces_of_body(body).into_iter().all(|f| {
        let origins = model.face_origins(f);
        origins.iter().for_each(|o| {
            creators.insert(o.step.clone());
        });
        !origins.is_empty()
            && origins
                .iter()
                .all(|o| matches!(o.role, Role::Start | Role::End | Role::Side { .. }))
            && matches!(
                model.geom().surface(model.face(f).surface),
                Surface::Plane(_)
            )
            && boundary(model, f).is_ok()
    }) && creators.len() > 1
}

/// Sweeps the one selected nonconvex planar face of `body` along its normal.
/// A body, selection or boundary outside that family reports `fallback`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn thicken(
    model: &mut Model,
    body: BodyId,
    refs: &[FaceRef],
    thickness: &Real,
    side: ThickenSide,
    step: &StepId,
    fallback: Refusal,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let Some(selected) = select(model, body, refs, dec)? else {
        return Err(fallback);
    };
    let plain = plain_prism(model, body, selected);
    if !plain && !planar_boolean_skin(model, body) {
        return Err(fallback);
    }
    let f = model.face(selected);
    let Surface::Plane(plane) = model.geom().surface(f.surface) else {
        return Err(fallback);
    };
    let mut plane: Plane = plane.clone();
    if f.sense == Sense::Reversed {
        plane.normal = plane.normal.scaled(&Real::from_i64(-1));
    }
    let loops = match boundary(model, selected) {
        Ok(loops) if !loops.is_empty() => loops,
        _ => return Err(fallback),
    };
    // The traversal contract: material lies to the left of every rim, so the
    // loops' signed area about the outward normal is positive.
    let mut area = Vector3::zero();
    for rims in &loops {
        for rim in rims {
            let a = rim.from.minus(&plane.origin);
            let b = rim.to.minus(&plane.origin);
            area = area.plus(&a.cross(&b));
        }
    }
    if sign(&area.dot(&plane.normal), dec)? != Sign::Positive {
        return Err(fallback);
    }
    let names = model.face_origins(selected).to_vec();
    let front = matches!(side, ThickenSide::Front);
    let distance = if front { thickness.clone() } else { -thickness };
    let Some(Surface::Plane(mut far)) =
        Surface::Plane(plane.clone()).parallel_support(&distance, dec)?
    else {
        unreachable!("a plane parallel never collapses")
    };
    let shift = far.origin.minus(&plane.origin);

    // Runs of collinear rims over one neighbour are one side face.
    struct Run {
        ring: usize,
        first: usize,
        last: usize,
        origins: Vec<Origin>,
    }
    let mut runs = Vec::new();
    for (ring, rims) in loops.iter().enumerate() {
        let m = rims.len();
        let direction = |i: usize| rims[i].to.minus(&rims[i].from);
        let mut breaks = Vec::new();
        for i in 0..m {
            let p = (i + m - 1) % m;
            let joined = rims[p].across == rims[i].across
                && sign(&direction(p).cross(&direction(i)).norm_squared(), dec)? == Sign::Zero;
            if !joined {
                breaks.push(i);
            }
        }
        if breaks.is_empty() {
            return Err(fallback);
        }
        for (k, first) in breaks.iter().enumerate() {
            let next = breaks[(k + 1) % breaks.len()];
            let last = (next + m - 1) % m;
            runs.push(Run {
                ring,
                first: *first,
                last,
                origins: rim_origins(model, selected, rims[*first].across, step, !plain)?,
            });
        }
    }

    let vertex = |model: &mut Model, p: &Point3| -> VertexId {
        let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
        model.add_vertex(Vertex { point })
    };
    struct Ring {
        base_edges: Vec<EdgeId>,
        top_edges: Vec<EdgeId>,
        rise: Vec<Option<EdgeId>>,
    }
    let line = |model: &mut Model, origin: &Point3, to: &Point3, a: VertexId, b: VertexId| {
        let curve = model.geom_mut().add_curve(Curve::Line(Line {
            origin: origin.clone(),
            dir: to.minus(origin),
        }));
        model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(a),
            end: Some(b),
        })
    };
    let mut rings = Vec::new();
    for (index, rims) in loops.iter().enumerate() {
        let m = rims.len();
        let tops: Vec<_> = rims.iter().map(|r| r.from.plus(&shift)).collect();
        let base: Vec<_> = rims.iter().map(|r| vertex(model, &r.from)).collect();
        let top: Vec<_> = tops.iter().map(|p| vertex(model, p)).collect();
        let mut base_edges = Vec::with_capacity(m);
        let mut top_edges = Vec::with_capacity(m);
        for i in 0..m {
            let j = (i + 1) % m;
            base_edges.push(line(model, &rims[i].from, &rims[j].from, base[i], base[j]));
            top_edges.push(line(model, &tops[i], &tops[j], top[i], top[j]));
        }
        let mut rise = vec![None; m];
        for run in runs.iter().filter(|r| r.ring == index) {
            for v in [run.first, (run.last + 1) % m] {
                if rise[v].is_none() {
                    rise[v] = Some(line(model, &rims[v].from, &tops[v], base[v], top[v]));
                }
            }
        }
        rings.push(Ring {
            base_edges,
            top_edges,
            rise,
        });
    }

    let mut faces = Vec::new();
    for run in &runs {
        let ring = &rings[run.ring];
        let m = loops[run.ring].len();
        let end = (run.last + 1) % m;
        let mut path = Vec::new();
        let mut i = run.first;
        loop {
            path.push((ring.base_edges[i], Sense::Same));
            if i == run.last {
                break;
            }
            i = (i + 1) % m;
        }
        path.push((ring.rise[end].expect("run end rises"), Sense::Same));
        let mut i = run.last;
        loop {
            path.push((ring.top_edges[i], Sense::Reversed));
            if i == run.first {
                break;
            }
            i = (i + m - 1) % m;
        }
        path.push((
            ring.rise[run.first].expect("run start rises"),
            Sense::Reversed,
        ));
        if !front {
            path = path
                .into_iter()
                .rev()
                .map(|(e, s)| (e, reverse(s)))
                .collect();
        }
        let rim = &loops[run.ring][run.first];
        let surface = Surface::Plane(Plane {
            origin: rim.from.clone(),
            normal: rim.to.minus(&rim.from).cross(&plane.normal),
            ref_dir: plane.normal.clone(),
        });
        faces.push(face(
            model,
            surface,
            Sense::Same,
            vec![path],
            None,
            run.origins.clone(),
        ));
    }

    // Outward normals: the face being thickened faces away from the sweep.
    let ordered = |edges: fn(&Ring) -> &Vec<EdgeId>, flip: bool| -> Vec<Vec<(EdgeId, Sense)>> {
        rings
            .iter()
            .map(|ring| {
                let path: Vec<_> = edges(ring).iter().map(|e| (*e, Sense::Same)).collect();
                if flip {
                    path.into_iter()
                        .rev()
                        .map(|(e, s)| (e, reverse(s)))
                        .collect()
                } else {
                    path
                }
            })
            .collect()
    };
    let negated = |p: &Plane| Plane {
        origin: p.origin.clone(),
        normal: p.normal.scaled(&Real::from_i64(-1)),
        ref_dir: p.ref_dir.clone(),
    };
    far.ref_dir = plane.ref_dir.clone();
    let (start_plane, end_plane) = if front {
        (negated(&plane), far)
    } else {
        (plane.clone(), negated(&far))
    };
    let offset: Vec<_> = names
        .iter()
        .map(|of| Origin {
            step: step.clone(),
            role: Role::Offset {
                of: Box::new(of.clone()),
            },
        })
        .collect();
    faces.push(face(
        model,
        Surface::Plane(start_plane),
        Sense::Same,
        ordered(|r| &r.base_edges, front),
        None,
        names,
    ));
    faces.push(face(
        model,
        Surface::Plane(end_plane),
        Sense::Same,
        ordered(|r| &r.top_edges, !front),
        None,
        offset,
    ));
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    Ok(model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    }))
}
