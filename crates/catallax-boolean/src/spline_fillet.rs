//! Trim routing for patch faces that share a physical rail with a fillet.

use crate::region::Location;
use catallax_numbers::{Decider, Outcome, Refusal, Sign, ZeroProof};
use catallax_shapes::{Curve, FilletCut, Point3, SplineFillet, SplineSupport, Surface};
use catallax_topology::{FaceId, Loop, Model, Sense};

/// Selects general certified winding instead of the complete-four-row test.
/// This recognises an adjacent face, and grants no membership or incidence.
pub(crate) fn borders_fillet(model: &Model, face: FaceId) -> bool {
    model.face(face).loops.iter().any(|lp| {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return false;
        };
        coedges.iter().any(|coedge| {
            model
                .coedges_of_edge(model.coedge(*coedge).edge)
                .into_iter()
                .any(|other| {
                    model.owner_of_coedge(other).is_some_and(|(_, neighbor)| {
                        neighbor != face
                            && matches!(
                                model.geom().surface(model.face(neighbor).surface),
                                Surface::BSpline(SplineSupport::SplineFillet(_))
                            )
                    })
                })
        })
    })
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("boolean.spline-fillet-domain", detail)
}

fn exact_sign(value: &catallax_numbers::Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    use catallax_numbers::StrictSign;
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {
            Ok(Sign::Zero)
        }
        Outcome::ProvenCoincident(_) => Err(bad("fillet cap membership needs an exact plane sign")),
        Outcome::Refused(error) => Err(error),
        Outcome::NotExact(error) => Err(error.into()),
    }
}

fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(&a.minus(b).norm_squared()) {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Ok(true),
        Outcome::ProvenCoincident(_) => {
            Err(bad("fillet endpoint correspondence needs exact identity"))
        }
        Outcome::Refused(error) => Err(error),
        Outcome::NotExact(error) => Err(error.into()),
    }
}

/// A complete physical band between two disjoint transverse cap cuts.
/// This grants no support incidence; the caller must prove that separately.
pub(crate) fn locate(
    model: &Model,
    face: FaceId,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    let f = model.face(face);
    let Surface::BSpline(SplineSupport::SplineFillet(fillet)) = model.geom().surface(f.surface)
    else {
        return Ok(None);
    };
    let [lp] = f.loops.as_slice() else {
        return Err(bad("a fillet band needs one physical loop"));
    };
    let Loop::Coedges(cs) = model.loop_(*lp) else {
        return Err(bad("a fillet band needs physical edges"));
    };
    if cs.len() != 4 {
        return Err(bad(
            "a fillet band needs exactly two rails and two cap cuts",
        ));
    }
    let mut cuts = Vec::new();
    let mut rails = Vec::new();
    let mut cycle = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for &coedge in cs {
        let coedge = model.coedge(coedge);
        if !seen.insert(coedge.edge) {
            return Err(bad("a fillet band repeats an edge"));
        }
        let edge = model.edge(coedge.edge);
        let (Some(a), Some(b)) = (edge.start, edge.end) else {
            return Err(bad("a fillet edge needs two finite vertices"));
        };
        cycle.push(if coedge.sense == Sense::Same {
            [a, b]
        } else {
            [b, a]
        });
        match model.geom().curve(edge.curve) {
            Curve::Intersection(branch) if branch.surfaces.contains(&f.surface) => {
                let plane = branch
                    .surfaces
                    .iter()
                    .find_map(|surface| match model.geom().surface(*surface) {
                        Surface::Plane(plane) => Some(plane),
                        _ => None,
                    })
                    .ok_or_else(|| bad("a fillet cap cut needs a plane"))?;
                let normal =
                    [&plane.normal.x, &plane.normal.y, &plane.normal.z].map(|v| v.as_rational());
                let [Some(x), Some(y), Some(z)] = normal else {
                    return Err(bad("a fillet cap plane needs rational coefficients"));
                };
                let normal = [x, y, z];
                let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
                let offset = normal
                    .iter()
                    .zip(origin)
                    .fold(catallax_numbers::Real::zero(), |sum, (n, x)| {
                        sum + x * &catallax_numbers::Real::from(n.clone())
                    })
                    .as_rational()
                    .ok_or_else(|| bad("a fillet cap plane needs a rational offset"))?;
                let cut = fillet
                    .plane_cut(normal, offset, dec)?
                    .ok_or_else(|| bad("a cap plane misses the fillet"))?;
                let witness = model.geom().position(branch.witness);
                if fillet.side_at(&witness, dec)? != Sign::Zero
                    || exact_sign(&witness.minus(&plane.origin).dot(&plane.normal), dec)?
                        != Sign::Zero
                {
                    return Err(bad("a cap-cut witness is not on both supports"));
                }
                let endpoints = [0, 1].map(|side| {
                    fillet
                        .rail_point(side, &cut.ends[side])
                        .expect("a cut ends at rail nodes")
                });
                let p = model.geom().position(model.vertex(a).point);
                let q = model.geom().position(model.vertex(b).point);
                let forward = if same(&p, &endpoints[0], dec)? && same(&q, &endpoints[1], dec)? {
                    true
                } else if same(&p, &endpoints[1], dec)? && same(&q, &endpoints[0], dec)? {
                    false
                } else {
                    return Err(bad("a cap cut does not join both exact rail endpoints"));
                };
                cuts.push((cut, (coedge.sense == Sense::Same) == forward));
            }
            Curve::BSpline(curve) => rails.push((edge, coedge.sense, curve)),
            _ => return Err(bad("a fillet band has an unproved boundary curve")),
        }
    }
    if cuts.len() != 2 || rails.len() != 2 {
        return Err(bad("a fillet band needs distinct rails and caps"));
    }
    if cycle
        .iter()
        .zip(cycle.iter().cycle().skip(1))
        .any(|(a, b)| a[1] != b[0])
    {
        return Err(bad("the fillet boundary is not sewn in loop order"));
    }
    cuts.sort_by_key(|(cut, _)| cut.spans);
    if cuts[0].0.spans[1] >= cuts[1].0.spans[0] {
        return Err(bad("the cap cuts are not proved disjoint and ordered"));
    }
    let ccw = f.sense == Sense::Same;
    if cuts[0].1 == ccw || cuts[1].1 != ccw {
        return Err(bad("cap directions disagree with the fillet chart"));
    }
    let mut found = [false; 2];
    for (edge, sense, curve) in rails {
        let mut matched = false;
        for (side, found) in found.iter_mut().enumerate() {
            let low = &cuts[0].0.ends[side];
            let high = &cuts[1].0.ends[side];
            let expected = fillet.rail_curve(side, low, high, dec)?;
            if !curve.same_definition(&expected) {
                continue;
            }
            if *found {
                return Err(bad("a fillet band repeats a rail"));
            }
            *found = true;
            matched = true;
            let endpoints = [low, high].map(|tau| fillet.rail_point(side, tau).unwrap());
            let p = model
                .geom()
                .position(model.vertex(edge.start.unwrap()).point);
            let q = model.geom().position(model.vertex(edge.end.unwrap()).point);
            let forward = if same(&p, &endpoints[0], dec)? && same(&q, &endpoints[1], dec)? {
                true
            } else if same(&p, &endpoints[1], dec)? && same(&q, &endpoints[0], dec)? {
                false
            } else {
                return Err(bad("a rail does not end at both cap nodes"));
            };
            if (edge.sense == Sense::Same) != forward
                || ((sense == Sense::Same) == forward) != ((side == 0) == ccw)
            {
                return Err(bad("rail directions disagree with the fillet chart"));
            }
        }
        if !matched {
            return Err(bad(
                "a fillet rail is not the carrier's exact restricted row",
            ));
        }
    }
    if !found.into_iter().all(|v| v) {
        return Err(bad("a fillet band is missing a rail"));
    }
    location(fillet, [&cuts[0].0, &cuts[1].0], point, dec).map(Some)
}

fn location(
    fillet: &SplineFillet,
    cuts: [&FilletCut; 2],
    point: &Point3,
    dec: &mut Decider,
) -> Result<Location, Refusal> {
    let mut boundary = false;
    for (index, cut) in cuts.into_iter().enumerate() {
        let value = [&point.x, &point.y, &point.z]
            .into_iter()
            .zip(&cut.normal)
            .fold(
                -catallax_numbers::Real::from(cut.offset.clone()),
                |sum, (x, n)| sum + x * &catallax_numbers::Real::from(n.clone()),
            );
        let sign = exact_sign(&value, dec)?;
        if sign == Sign::Zero {
            boundary = true;
        } else if (sign == cut.before) == (index == 0) {
            return Ok(Location::Outside);
        }
    }
    Ok(if boundary || fillet.on_boundary(point, dec)? {
        Location::Boundary
    } else {
        Location::Inside
    })
}
