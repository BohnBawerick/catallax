//! Exact face points beside a fitted contour generator.
//!
//! A fitted generator lies within its deviation of the surface, on either
//! side, so a ray from one of its points may start just inside the solid.
//! The ray instead starts from an exact point of the face whose outward
//! normal strictly faces the viewer, found as the face point nearest to the
//! fitted point moved along the view, toward the viewer or away from it:
//! the facing side of a generator lies ahead of it where the face bends
//! away from the view and behind it where the face bends toward the view.
//! From there the ray leaves the surface where the true generator's ray
//! grazes it, and meets the solid again wherever the generator's would (a
//! concave generator meets its own face at once, as it should).
//!
//! An exact generator on a spline face has the same trouble the other way:
//! its ray is tangent to the patch, which the patch intersector refuses, so
//! it starts beside the face too.
//!
//! A fitted edge has the same trouble with both its faces. Its ray starts
//! from an exact point of one of them moved across the edge, within the
//! face's tangent plane, to the face's own side (`across`), where the face's
//! facing is read exactly. Exact face membership finds the side once per
//! edge and face; later pieces of the edge reuse it.
//!
//! The offset starts at the fit's rule T4 tolerance, or rule T1's floor for
//! an exact generator, and doubles until the exact checks pass. Quadric and
//! torus faces take the Boolean engine's rational witness
//! (`surface_witness`); rational patch faces invert their patch by
//! Gauss-Newton steps (Piegl and Tiller, "The NURBS Book", 2nd ed., section
//! 6.1) from the nearest knot-span corner or middle, then evaluate it
//! exactly at the rational parameters reached; an extruded profile takes its
//! height along the extrusion and the nearest profile point.

use catallax_boolean::{face_contains, surface_witness};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{Point3, RationalPatch, SplineSupport, Surface, Vector3};
use catallax_topology::{FaceId, Model, RationalPatchFace, Sense};

use super::cross::SOLVE_STEPS;
use super::geo::{Box3, V3, box_distance, distance};
use crate::queries::sign;
use crate::views::float::{above, dot, rational};

fn rational_point(p: V3) -> Result<Point3, Refusal> {
    Ok(Point3::new(
        Real::from(rational(p[0])?),
        Real::from(rational(p[1])?),
        Real::from(rational(p[2])?),
    ))
}

fn to64(p: &Point3, floor: &Rational, dec: &mut Decider) -> Result<V3, Refusal> {
    Ok([
        crate::views::float::double(&p.x, floor, dec)?,
        crate::views::float::double(&p.y, floor, dec)?,
        crate::views::float::double(&p.z, floor, dec)?,
    ])
}

/// The patch point nearest `target`, exactly on the patch.
fn invert(
    patch: &RationalPatch,
    target: V3,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<(f64, Point3), Refusal> {
    let carrier = patch.carrier();
    let domain = carrier.domain();
    let bounds = domain.map(|[a, b]| [a.to_f64_nearest(), b.to_f64_nearest()]);
    let starts = |axis: usize| -> Vec<f64> {
        let mut out = Vec::new();
        for [a, b] in carrier.spans(axis) {
            let (a, b) = (a.to_f64_nearest(), b.to_f64_nearest());
            out.extend([a, (a + b) / 2.0]);
        }
        out.push(bounds[axis][1]);
        out
    };
    // Bounds rounded to binary64 can reach past the exact domain.
    let at = |u: f64, v: f64| -> Result<[Real; 2], Refusal> {
        let clamp = |x: f64, [lo, hi]: [&Rational; 2]| -> Result<Real, Refusal> {
            Ok(Real::from(rational(x)?.max(lo.clone()).min(hi.clone())))
        };
        Ok([clamp(u, domain[0])?, clamp(v, domain[1])?])
    };
    let mut best = (f64::INFINITY, bounds[0][0], bounds[1][0]);
    for &u in &starts(0) {
        for &v in &starts(1) {
            dec.check_deadline()?;
            let p = to64(&patch.point_at(&at(u, v)?, dec)?, floor, dec)?;
            let gap = distance(p, target);
            if gap < best.0 {
                best = (gap, u, v);
            }
        }
    }
    let (_, mut u, mut v) = best;
    for _ in 0..SOLVE_STEPS {
        dec.check_deadline()?;
        let parameters = at(u, v)?;
        let p = to64(&patch.point_at(&parameters, dec)?, floor, dec)?;
        let [du, dv] = patch.derivatives_at(&parameters, dec)?;
        let (du, dv) = (to64v(&du, floor, dec)?, to64v(&dv, floor, dec)?);
        let r = [target[0] - p[0], target[1] - p[1], target[2] - p[2]];
        // Normal equations of the linearised least-squares step.
        let (a, b, c) = (dot(du, du), dot(du, dv), dot(dv, dv));
        let (e, f) = (dot(du, r), dot(dv, r));
        let det = a * c - b * b;
        if det == 0.0 || !det.is_finite() {
            break;
        }
        let next = (
            (u + (e * c - b * f) / det).clamp(bounds[0][0], bounds[0][1]),
            (v + (a * f - b * e) / det).clamp(bounds[1][0], bounds[1][1]),
        );
        let resolution = |x: [f64; 2]| 2.0 * f64::EPSILON * x[0].abs().max(x[1].abs()).max(1.0);
        let done = (next.0 - u).abs() <= resolution(bounds[0])
            && (next.1 - v).abs() <= resolution(bounds[1]);
        (u, v) = next;
        if done {
            break;
        }
    }
    let q = patch.point_at(&at(u, v)?, dec)?;
    Ok((distance(to64(&q, floor, dec)?, target), q))
}

fn to64v(v: &Vector3, floor: &Rational, dec: &mut Decider) -> Result<V3, Refusal> {
    Ok([
        crate::views::float::double(&v.x, floor, dec)?,
        crate::views::float::double(&v.y, floor, dec)?,
        crate::views::float::double(&v.z, floor, dec)?,
    ])
}

/// An exact point of `face` near `target`, or `None` where its surface has
/// no supported construction there.
fn near(
    scratch: &Model,
    face: FaceId,
    target: V3,
    width: f64,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<Option<Point3>, Refusal> {
    let surface = scratch.geom().surface(scratch.face(face).surface);
    match surface {
        Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
            Ok(Some(invert(patch, target, floor, dec)?.1))
        }
        Surface::BSpline(SplineSupport::RationalAtlas(_)) => {
            let cells = scratch
                .face_chart_cover(face)
                .map_or_else(|| vec![face], <[FaceId]>::to_vec);
            let mut best: Option<(f64, Point3)> = None;
            for cell in cells {
                let Some(chart) = RationalPatchFace::of_chart_face(scratch, cell, dec)? else {
                    return Ok(None);
                };
                let found = invert(chart.support(), target, floor, dec)?;
                if best.as_ref().is_none_or(|b| found.0 < b.0) {
                    best = Some(found);
                }
            }
            Ok(best.map(|b| b.1))
        }
        Surface::BSpline(SplineSupport::BezierExtrusion(extrusion)) => {
            // S(u, v) = B(u) + v w, with w normal to the profile's plane: the
            // nearest point has v from the target's height along w and u
            // nearest on the profile to the target moved down to it.
            let profile = extrusion.lower_edge();
            let w = extrusion.direction();
            let w64 = to64v(w, floor, dec)?;
            let base = to64(profile.control(0), floor, dec)?;
            let rise = std::array::from_fn(|k| target[k] - base[k]);
            let v = (dot(rise, w64) / dot(w64, w64)).clamp(0.0, 1.0);
            let knots: Vec<f64> = profile
                .knots()
                .iter()
                .map(Rational::to_f64_nearest)
                .collect();
            let mut ctrl = Vec::with_capacity(profile.control_count());
            for i in 0..profile.control_count() {
                let c = to64(profile.control(i), floor, dec)?;
                let weight = crate::views::float::double(profile.weight(i), floor, dec)?;
                ctrl.push([weight * c[0], weight * c[1], weight * c[2], weight]);
            }
            let curve = super::geo::Geo::Spline {
                spans: super::geo::decompose(profile.degree(), &knots, &ctrl),
            };
            let u = curve.param(std::array::from_fn(|k| target[k] - v * w64[k]));
            let on_profile = profile.point_at(&Real::from(rational(u)?), dec)?;
            Ok(Some(on_profile.plus(&w.scaled(&Real::from(rational(v)?)))))
        }
        Surface::BSpline(_) => Ok(None),
        other => {
            let width = rational(width)?;
            Ok(surface_witness(other, &rational_point(target)?, &width, dec)?.map(|p| p.to_real()))
        }
    }
}

/// An exact point of `face` just beside the fitted edge point `p64`, with
/// the sign of its outward normal along the view direction `d`, and the
/// side it lies on. The point is the face point nearest to `p64` moved
/// across the edge, perpendicular to the edge's binary64 `tangent` within
/// the face's tangent plane: along `normal x tangent` for side `1`, against
/// it for `-1`. The offset starts at `start`, which must clear the fit's
/// tolerance, and doubles while it stays within `limit`. No attempt goes
/// past `limit`: a `start` beyond it refuses before any.
///
/// A face lies on one side of its boundary edge all along it, so `known`,
/// the side an earlier piece of the edge confirmed, is tried first. Near a
/// point of a boundary curve and away from the face's other boundaries, a
/// face is the region on its own side of the curve: where `others`, the
/// boxes of the face's other boundary edges, all lie farther than twice the
/// offset from `p64`, and the surface point lies within half the offset of
/// its target, the point on the known side is in the face. Elsewhere, with
/// no side known yet, or with `others` `None` for a face with a point
/// boundary, exact face membership decides.
#[allow(clippy::too_many_arguments)]
pub(super) fn across(
    scratch: &mut Model,
    face: FaceId,
    p64: V3,
    tangent: V3,
    known: Option<f64>,
    others: Option<&[Box3]>,
    d: &Vector3,
    start: f64,
    limit: f64,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<(Point3, i8, f64), Refusal> {
    let refused = || {
        Refusal::new(
            "app.view-sample",
            "no exact face point beside a fitted edge lies in its face",
        )
    };
    if start > limit {
        return Err(refused());
    }
    let outward: i8 = if scratch.face(face).sense == Sense::Same {
        1
    } else {
        -1
    };
    let surface = scratch.face(face).surface;
    // The face's tangent plane at the edge, from its normal at a surface
    // point next to the edge.
    let foot = near(scratch, face, p64, start / 4.0, floor, dec)?.ok_or_else(refused)?;
    let normal = to64v(
        &scratch.geom().surface(surface).normal_at(&foot, dec)?,
        floor,
        dec,
    )?;
    let side = crate::views::float::cross(normal, tangent);
    let width = dot(side, side).sqrt();
    if !above(width, 0.0) || !width.is_finite() {
        return Err(refused());
    }
    let side = side.map(|x| x / width);
    let first = known.unwrap_or(1.0);
    let mut offset = start;
    while offset <= limit {
        for (k, way) in [first, -first].into_iter().enumerate() {
            dec.check_deadline()?;
            let target = std::array::from_fn(|k| p64[k] + way * offset * side[k]);
            let Some(q) = near(scratch, face, target, offset / 4.0, floor, dec)? else {
                continue;
            };
            let clear = known.is_some()
                && k == 0
                && others
                    .is_some_and(|boxes| boxes.iter().all(|b| box_distance(b, p64) > 2.0 * offset))
                && distance(to64(&q, floor, dec)?, target) <= offset / 2.0;
            if clear || face_contains(scratch, face, &q, dec)? {
                let n = scratch.geom().surface(surface).normal_at(&q, dec)?;
                let facing = sign(&n.dot(d), dec)? * outward;
                return Ok((q, facing, way));
            }
        }
        offset *= 2.0;
    }
    Err(refused())
}

/// An exact point of `face` beside the fitted point `p64`, whose outward
/// normal strictly faces the view direction `d`. The offset starts at
/// `start` and doubles while it stays within `limit`, the body's size.
#[allow(clippy::too_many_arguments)]
pub(super) fn beside(
    scratch: &mut Model,
    face: FaceId,
    p64: V3,
    d64: V3,
    d: &Vector3,
    start: f64,
    limit: f64,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let outward = if scratch.face(face).sense == Sense::Same {
        1
    } else {
        -1
    };
    let mut offset = start;
    while offset <= limit {
        // Toward the viewer the nearest face point faces it where the face
        // bends away from the view (a visible generator); where it bends
        // toward it, the facing side lies behind the generator.
        for side in [1.0, -1.0] {
            dec.check_deadline()?;
            let target = std::array::from_fn(|k| p64[k] + side * offset * d64[k]);
            if let Some(q) = near(scratch, face, target, offset / 4.0, floor, dec)?
                && face_contains(scratch, face, &q, dec)?
            {
                let surface = scratch.face(face).surface;
                let n = scratch.geom().surface(surface).normal_at(&q, dec)?;
                if sign(&n.dot(d), dec)? * outward > 0 {
                    return Ok(q);
                }
            }
        }
        offset *= 2.0;
    }
    Err(Refusal::new(
        "app.view-sample",
        "no exact face point beside a fitted contour generator faces the viewer",
    ))
}
