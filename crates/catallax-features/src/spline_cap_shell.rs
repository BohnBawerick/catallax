//! Exact cap-distance envelopes of complete polynomial profile prisms.

use std::collections::BTreeSet;

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign, ZeroProof};
use catallax_recipe::{FaceRef, Origin, Role, StepId};
use catallax_shapes::{Curve, Plane, Point3, SplineSupport, Surface, Vector3};
use catallax_topology::{BezierExtrusionBand, BodyId, FaceId, Loop, Model, Sense};

use crate::polyhedron::{ConvexCell, PlaneFace};

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.unsupported-envelope", detail)
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {
            Ok(Sign::Zero)
        }
        Outcome::ProvenCoincident(_) => Err(bad("spline cap incidence needs an exact proof")),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

fn outward(model: &Model, face: FaceId, plane: &Plane) -> Plane {
    let mut plane = plane.clone();
    if model.face(face).sense == Sense::Reversed {
        plane.normal = plane.normal.scaled(&Real::from_i64(-1));
    }
    plane
}

fn select(
    model: &Model,
    faces: &[FaceId],
    references: &[FaceRef],
) -> Result<BTreeSet<FaceId>, Refusal> {
    let mut opened = BTreeSet::new();
    for reference in references {
        let found: Vec<_> = faces
            .iter()
            .copied()
            .filter(|f| reference.matches_origins(model.face_origins(*f)))
            .filter(|f| {
                reference.touching.iter().all(|pattern| {
                    model.face(*f).loops.iter().any(|lp| {
                        let Loop::Coedges(cs) = model.loop_(*lp) else {
                            return false;
                        };
                        cs.iter().any(|c| {
                            model
                                .coedges_of_edge(model.coedge(*c).edge)
                                .iter()
                                .any(|other| {
                                    model.owner_of_coedge(*other).is_some_and(|(_, neighbour)| {
                                        neighbour != *f
                                            && faces.contains(&neighbour)
                                            && model
                                                .face_origins(neighbour)
                                                .iter()
                                                .any(|o| pattern.matches(o))
                                    })
                                })
                        })
                    })
                })
            })
            .collect();
        if found.is_empty() {
            return Err(Refusal::new(
                "reference.lost-face",
                "no spline prism face resolves the open reference",
            ));
        }
        if reference.near.is_some() && found.len() > 1 {
            return Err(bad(
                "nearest spline prism selection needs a certified finite-face distance comparison",
            ));
        }
        opened.extend(found);
    }
    Ok(opened)
}

// Complete cap loops and rectangular sewn side bands prove the same finite
// cross-section at every axial level, including every hole in the cap domain.
fn caps(
    model: &Model,
    body: BodyId,
    axis: &Vector3,
    dec: &mut Decider,
) -> Result<Vec<(FaceId, Plane)>, Refusal> {
    let entity = model.body(body);
    let [lump] = entity.lumps.as_slice() else {
        return Err(bad("a spline cap shell requires one material lump"));
    };
    if model.lump(*lump).shells.len() != 1 || !model.check(body).is_empty() {
        return Err(bad(
            "a spline cap shell requires one sewn shell without voids",
        ));
    }
    let faces = model.faces_of_body(body);
    let mut caps = Vec::new();
    for &face in &faces {
        if let Surface::Plane(plane) = model.geom().surface(model.face(face).surface)
            && sign(&plane.normal.cross(axis).norm_squared(), dec)? == Sign::Zero
        {
            caps.push((face, outward(model, face, plane)));
        }
    }
    if caps.len() != 2 {
        return Err(bad(
            "a spline cap shell needs two complete parallel end caps",
        ));
    }
    if sign(&caps[0].1.normal.dot(&caps[1].1.normal), dec)? != Sign::Negative {
        return Err(bad("spline prism cap normals must oppose one another"));
    }
    let normal = &caps[0].1.normal;
    let level = |point: &catallax_shapes::Point3, dec: &mut Decider| -> Result<usize, Refusal> {
        for (i, (_, cap)) in caps.iter().enumerate() {
            if sign(&point.minus(&cap.origin).dot(normal), dec)? == Sign::Zero {
                return Ok(i);
            }
        }
        Err(bad(
            "a spline prism boundary vertex is detached from its caps",
        ))
    };
    if sign(&caps[1].1.origin.minus(&caps[0].1.origin).dot(normal), dec)? != Sign::Negative {
        return Err(bad(
            "spline prism caps do not enclose positive material width",
        ));
    }
    for &face in &faces {
        dec.check_deadline()?;
        if caps.iter().any(|(cap, _)| *cap == face) {
            continue;
        }
        let f = model.face(face);
        let rails = match model.geom().surface(f.surface) {
            Surface::BSpline(SplineSupport::BezierExtrusion(support)) => {
                if sign(&support.direction().cross(axis).norm_squared(), dec)? != Sign::Zero {
                    return Err(bad(
                        "spline prism sides must share one normal extrusion axis",
                    ));
                }
                BezierExtrusionBand::of_face(model, face, dec)?
                    .ok_or_else(|| bad("missing exact polynomial band"))?
                    .rails()
            }
            Surface::Plane(plane) if sign(&plane.normal.dot(axis), dec)? == Sign::Zero => {
                let [lp] = f.loops.as_slice() else {
                    return Err(bad("a planar prism side needs one rectangular loop"));
                };
                let Loop::Coedges(cs) = model.loop_(*lp) else {
                    return Err(bad("a planar prism side needs finite edges"));
                };
                if cs.len() != 4 {
                    return Err(bad("a planar prism side needs four complete edges"));
                }
                let mut rails = Vec::new();
                for c in cs {
                    let edge = model.edge(model.coedge(*c).edge);
                    let Curve::Line(line) = model.geom().curve(edge.curve) else {
                        return Err(bad("a planar prism side needs straight edges"));
                    };
                    let (Some(a), Some(b)) = (edge.start, edge.end) else {
                        return Err(bad("a planar prism edge needs finite endpoints"));
                    };
                    let a = model.geom().position(model.vertex(a).point);
                    let b = model.geom().position(model.vertex(b).point);
                    if level(&a, dec)? == level(&b, dec)? {
                        rails.push(model.coedge(*c).edge);
                    } else if sign(&line.dir.cross(axis).norm_squared(), dec)? != Sign::Zero {
                        return Err(bad(
                            "a planar prism generator must be parallel to the extrusion axis",
                        ));
                    }
                }
                rails
                    .try_into()
                    .map_err(|_| bad("a planar prism side needs two cap rails"))?
            }
            _ => {
                return Err(bad(
                    "spline cap shells accept only complete polynomial and planar prism sides",
                ));
            }
        };
        let mut used = [false; 2];
        for rail in rails {
            let edge = model.edge(rail);
            let vertex = edge
                .start
                .ok_or_else(|| bad("a cap rail needs an endpoint"))?;
            let slot = level(&model.geom().position(model.vertex(vertex).point), dec)?;
            if std::mem::replace(&mut used[slot], true)
                || !model.coedges_of_edge(rail).iter().any(|c| {
                    model
                        .owner_of_coedge(*c)
                        .is_some_and(|(_, f)| f == caps[slot].0)
                })
            {
                return Err(bad(
                    "every prism rail must be sewn to its distinct complete cap",
                ));
            }
        }
    }
    Ok(caps)
}

pub(super) fn shell(
    model: &mut Model,
    body: BodyId,
    open: &[FaceRef],
    thickness: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<Option<BodyId>, Refusal> {
    let faces = model.faces_of_body(body);
    let axis =
        faces.iter().find_map(
            |face| match model.geom().surface(model.face(*face).surface) {
                Surface::BSpline(SplineSupport::BezierExtrusion(support)) => {
                    Some(support.direction().clone())
                }
                _ => None,
            },
        );
    let Some(axis) = axis else {
        return Ok(None);
    };
    let opened = select(model, &faces, open)?;
    if opened.len() == faces.len() {
        return Ok(Some(catallax_topology::empty(model)));
    }
    let caps = caps(model, body, &axis, dec)?;
    let kept: Vec<_> = caps
        .iter()
        .filter(|(face, _)| !opened.contains(face))
        .collect();
    if faces
        .iter()
        .any(|f| !opened.contains(f) && !kept.iter().any(|(cap, _)| cap == f))
    {
        return Err(bad(
            "retained spline prism sides need normal offsets and finite rim envelopes",
        ));
    }
    let width = caps[1]
        .1
        .origin
        .minus(&caps[0].1.origin)
        .dot(&caps[1].1.normal)
        / caps[1].1.normal.norm_squared().sqrt();
    if sign(
        &(&width - &(thickness * &Real::from_i64(kept.len() as i64))),
        dec,
    )? != Sign::Positive
    {
        return Ok(Some(body));
    }
    let mut hull: Vec<_> = model
        .vertices_of_body(body)
        .into_iter()
        .map(|v| model.geom().position(model.vertex(v).point))
        .collect();
    for face in &faces {
        if let Some(band) = BezierExtrusionBand::of_face(model, *face, dec)? {
            for height in band.axial() {
                let shift = band.support().direction().scaled(height);
                hull.extend(band.support().controls().iter().map(|p| p.plus(&shift)));
            }
        }
    }
    // Polynomial convex hulls and planar polygon vertices enclose the whole
    // checked prism. Only an outer box is needed, not solved extrema.
    let mut minimum = [None, None, None];
    let mut maximum = [None, None, None];
    for point in hull {
        for (coordinate, value) in [&point.x, &point.y, &point.z].into_iter().enumerate() {
            let bound = dec.enclose_at(value, 128)?;
            minimum[coordinate] = Some(
                minimum[coordinate]
                    .take()
                    .map_or(bound.lo(), |old: Rational| old.min(bound.lo())),
            );
            maximum[coordinate] = Some(
                maximum[coordinate]
                    .take()
                    .map_or(bound.hi(), |old: Rational| old.max(bound.hi())),
            );
        }
    }
    let mut planes = Vec::new();
    for coordinate in 0..3 {
        for upper in [false, true] {
            let mut origin = [Real::zero(), Real::zero(), Real::zero()];
            origin[coordinate] = Real::from(if upper {
                maximum[coordinate]
                    .as_ref()
                    .ok_or_else(|| bad("missing prism hull"))?
                    + &Rational::one()
            } else {
                minimum[coordinate]
                    .as_ref()
                    .ok_or_else(|| bad("missing prism hull"))?
                    - &Rational::one()
            });
            let mut normal = [Real::zero(), Real::zero(), Real::zero()];
            normal[coordinate] = Real::from_i64(if upper { 1 } else { -1 });
            let mut reference = [Real::zero(), Real::zero(), Real::zero()];
            reference[(coordinate + 1) % 3] = Real::one();
            let [x, y, z] = origin;
            let [nx, ny, nz] = normal;
            let [rx, ry, rz] = reference;
            planes.push(PlaneFace {
                plane: Plane {
                    origin: Point3::new(x, y, z),
                    normal: Vector3::new(nx, ny, nz),
                    ref_dir: Vector3::new(rx, ry, rz),
                },
                origins: Vec::new(),
            });
        }
    }
    // Orthogonal projection stays in the complete finite cap domain. The
    // cavity is exactly the material farther than thickness from every kept cap.
    for (face, cap) in &kept {
        let shift = cap
            .normal
            .scaled(&(-thickness / cap.normal.norm_squared().sqrt()));
        planes.push(PlaneFace {
            plane: Plane {
                origin: cap.origin.plus(&shift),
                normal: cap.normal.clone(),
                ref_dir: cap.ref_dir.clone(),
            },
            origins: model
                .face_origins(*face)
                .iter()
                .map(|of| Origin {
                    step: step.clone(),
                    role: Role::Inner {
                        of: Box::new(of.clone()),
                    },
                })
                .collect(),
        });
    }
    let cavity = ConvexCell::new(planes, dec)?.build(model);
    let answer =
        catallax_boolean::boolean(model, body, cavity, catallax_boolean::BoolOp::Subtract, dec)
            .map_err(|error| error.refusal)?;
    Ok(Some(answer.body))
}
