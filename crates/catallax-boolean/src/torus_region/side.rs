//! The face's winding cell from a complete transverse meridian crossing.
//!
//! A certified tube chart omits a parallel disjoint from all boundaries.
//! Its zero cell is the bottom of the cylinder chart. The lowest crossing
//! of a transverse meridian has that cell immediately below it. The coedge
//! orientation tells which adjacent cell is the face's left side. Complete
//! finite curve/plane intersections supply the crossing order, including
//! boundaries with no implicit curve equations.

use crate::geometry::{coincident, natural_ends, on_arc, oriented, product, refusal, sign};
use crate::project::Projection;
use catallax_intersect::{CurveSurface, intersect_curve_surface_framed, tangent_at};
use catallax_numbers::{
    Decider, Fact, FactKind, FactOrigin, Rational, Real, RealPoly, Refusal, Sign,
};
use catallax_shapes::{
    Curve, CurveId, GeomKey, Plane, PointDef, PointId, Surface, SurfaceId, Vector3,
};
use catallax_topology::{FaceId, Loop, Model, Sense};

// The shared finite spline/plane entry point currently admits rational cut
// coefficients. Torus charts also need algebraic coefficients: use the
// existing exact RealPoly root solver, retaining each finite knot domain.
// A knot hit, tangent or contained span rejects this meridian, so every
// accepted crossing has a regular, unique local parameter side.
fn meridian_hits(
    model: &mut Model,
    curve: CurveId,
    plane: SurfaceId,
    dec: &mut Decider,
) -> Result<Option<Vec<(PointId, Vector3)>>, Refusal> {
    if let Curve::BSpline(spline) = model.geom().curve(curve).clone() {
        let Surface::Plane(cut) = model.geom().surface(plane).clone() else {
            unreachable!("meridian plane");
        };
        let normal = [&cut.normal.x, &cut.normal.y, &cut.normal.z];
        let offset = cut.origin.x.clone() * cut.normal.x.clone()
            + cut.origin.y.clone() * cut.normal.y.clone()
            + cut.origin.z.clone() * cut.normal.z.clone();
        let mut hits = Vec::new();
        let before = dec.stats().coincident_separation;
        for (span, [lo, hi]) in spline.spans() {
            let homogeneous = spline.query_span_polynomials(span, dec)?;
            let equation = (0..3).fold(RealPoly::zero(1), |sum, i| {
                sum + homogeneous[i].scale(normal[i])
            }) - homogeneous[3].scale(&offset);
            if equation.is_zero(dec)? {
                return Ok(None);
            }
            for parameter in equation.real_roots(dec)? {
                let lower = sign(
                    &(&parameter - &Real::from(lo.clone())),
                    dec,
                    "boolean.classify-ambiguous",
                )?;
                let upper = sign(
                    &(&Real::from(hi.clone()) - &parameter),
                    dec,
                    "boolean.classify-ambiguous",
                )?;
                if lower == Sign::Negative || upper == Sign::Negative {
                    continue;
                }
                if lower == Sign::Zero
                    || upper == Sign::Zero
                    || sign(
                        &equation
                            .partial(0)
                            .eval_real(std::slice::from_ref(&parameter)),
                        dec,
                        "boolean.classify-ambiguous",
                    )? == Sign::Zero
                {
                    return Ok(None);
                }
                let position = spline.point_at(&parameter, dec)?;
                let tangent = spline.first_at(&parameter, dec)?;
                if dec.stats().coincident_separation != before {
                    return Err(refusal(
                        "numbers.geometric-gap-unproven",
                        "a torus meridian needs exact finite parameter roots",
                    ));
                }
                let point = model.geom_mut().add_point(PointDef::Explicit(position));
                for b in [GeomKey::Curve(curve), GeomKey::Surface(plane)] {
                    model.geom_mut().record_fact(Fact {
                        kind: FactKind::Incident,
                        a: GeomKey::Point(point),
                        b,
                        at: None,
                        origin: FactOrigin::Construction("boolean.torus-meridian"),
                    });
                }
                hits.push((point, tangent));
            }
        }
        return Ok(Some(hits));
    }
    let meeting = intersect_curve_surface_framed(model.geom_mut(), curve, plane, dec)?;
    if !meeting.marks.is_empty() {
        return Err(refusal(
            "numbers.geometric-gap-unproven",
            "a torus side needs exact meridian intersections",
        ));
    }
    match meeting.value {
        CurveSurface::Contained => Ok(None),
        CurveSurface::Meets(hits) => hits
            .into_iter()
            .map(|hit| {
                Ok((
                    hit.point,
                    tangent_at(model.geom(), curve, hit.point, dec)?[1].clone(),
                ))
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some),
        _ => Ok(Some(Vec::new())),
    }
}

pub(crate) fn target(
    model: &Model,
    face: FaceId,
    chart: &Projection,
    dec: &mut Decider,
) -> Result<i64, Refusal> {
    let f = model.face(face);
    let surface = crate::geometry::chart(model.geom().surface(f.surface));
    let Surface::Torus(torus) = surface.as_ref() else {
        unreachable!("torus face");
    };
    // Meridian planes need directions, not unit coordinates. Preserve exact
    // rational frames instead of introducing cancelling square roots.
    let axis = torus.axis.clone();
    let u = torus.ref_dir.clone();
    let v = axis.cross(&u);
    let mut fork = model.clone();
    for choice in 0..crate::geometry::DIRECTION_ATTEMPTS {
        dec.check_deadline()?;
        let k = i64::from(choice);
        let t = Rational::from_ratio(1, k / 4 + 2).expect("positive denominator");
        let over = &Rational::one() + &(&t * &t);
        let c = &(&Rational::one() - &(&t * &t)) / &over;
        let s = &(Rational::from_i64(2) * t) / &over;
        let (c, s) = match k % 4 {
            0 => (c, s),
            1 => (-s, c),
            2 => (-c, -s),
            _ => (s, -c),
        };
        let direction = u.scaled(&Real::from(c)).plus(&v.scaled(&Real::from(s)));
        let mut meridian = Plane {
            origin: torus.center.clone(),
            normal: axis.cross(&direction),
            ref_dir: axis.clone(),
        };
        for coordinate in [
            &mut meridian.origin.x,
            &mut meridian.origin.y,
            &mut meridian.origin.z,
            &mut meridian.normal.x,
            &mut meridian.normal.y,
            &mut meridian.normal.z,
            &mut meridian.ref_dir.x,
            &mut meridian.ref_dir.y,
            &mut meridian.ref_dir.z,
        ] {
            if coordinate.rational_ref().is_none()
                && let Some(value) = dec.rational_value(coordinate)?
            {
                *coordinate = value.into();
            }
        }
        let plane = fork.geom_mut().add_surface(Surface::Plane(meridian));
        let mut lowest: Option<(Real, Sign)> = None;
        let mut transverse = true;
        'boundary: for lp in &f.loops {
            let Loop::Coedges(coedges) = model.loop_(*lp) else {
                continue;
            };
            for coedge in coedges {
                let coedge = model.coedge(*coedge);
                let edge = model.edge(coedge.edge).clone();
                let (from, to) = natural_ends(model, coedge.edge);
                let hits = match meridian_hits(&mut fork, edge.curve, plane, dec)? {
                    None => {
                        transverse = false;
                        break 'boundary;
                    }
                    Some(hits) => hits,
                };
                for (point, tangent) in hits {
                    let position = fork.geom().position(point);
                    if sign(
                        &position.minus(&torus.center).dot(&direction),
                        dec,
                        "boolean.classify-ambiguous",
                    )? != Sign::Positive
                        || !on_arc(
                            &fork,
                            edge.curve,
                            from,
                            to,
                            point,
                            dec,
                            "boolean.classify-ambiguous",
                        )?
                    {
                        continue;
                    }
                    for endpoint in [from, to].into_iter().flatten() {
                        if coincident(
                            fork.geom().same_point(endpoint, point, dec),
                            "boolean.classify-ambiguous",
                        )? {
                            transverse = false;
                            break 'boundary;
                        }
                    }
                    let tangent = oriented(&tangent, product(edge.sense, coedge.sense));
                    let crossing = sign(
                        &axis.dot(&position.minus(&torus.center).cross(&tangent)),
                        dec,
                        "boolean.classify-ambiguous",
                    )?;
                    if crossing == Sign::Zero {
                        transverse = false;
                        break 'boundary;
                    }
                    let height = chart.torus_coordinate(&position);
                    match &lowest {
                        None => lowest = Some((height, crossing)),
                        Some((previous, _)) => {
                            match sign(&(&height - previous), dec, "boolean.classify-ambiguous")? {
                                Sign::Negative => lowest = Some((height, crossing)),
                                Sign::Zero => {
                                    transverse = false;
                                    break 'boundary;
                                }
                                Sign::Positive => {}
                            }
                        }
                    }
                }
            }
        }
        if transverse && let Some((_, crossing)) = lowest {
            // The tube stereographic coordinate reverses orientation.
            // Winding jumps by sign(d azimuth/d edge) on crossing upward.
            return Ok(match (f.sense, crossing) {
                (Sense::Same, Sign::Negative) => -1,
                (Sense::Reversed, Sign::Positive) => 1,
                _ => 0,
            });
        }
    }
    Err(refusal(
        "boolean.torus-chart-required",
        "no complete transverse meridian selected a torus face side",
    ))
}
