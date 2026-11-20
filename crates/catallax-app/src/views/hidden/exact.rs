//! Exact cuts where one stroke's projection meets another's, for lines and
//! conics. Each reduces to a curve against a plane.
//!
//! A point of A lies over B on the sheet when the line through it along the
//! view meets B, that is when it lies on B's extrusion along the view.
//! OpenCASCADE's `HLRBRep_Intersector` intersects the projected curves in
//! the sheet instead (`IntCurve_IntConicConic` for lines and conics); the
//! extrusion keeps the points exact in model space.
//!
//! - B a line: its extrusion is the plane through it holding the view.
//! - B a conic seen edge on (its plane holds the view): that plane.
//! - A a line, B a conic seen at an angle: B against A's extrusion plane,
//!   each point carried along the view onto A.
//! - A and B conics in parallel planes, seen at an angle: B slid along the
//!   view into A's plane; two coplanar circles meet on their radical plane
//!   (Salmon, "A Treatise on Conic Sections", chapter 6, the radical axis).
//!
//! When A lies in B's extrusion plane, their projections overlap, and A is
//! cut where it lies over B's ends and folds instead. Other pairs go to the
//! binary64 search.

use catallax_intersect::{
    CurveSurface, intersect_curve_surface_framed, intersect_curve_surface_with_line_parameter,
};
use catallax_numbers::{Decider, MarkLog, Real, Refusal};
use catallax_shapes::{CurveId, Plane, Point3, PointDef, PointId, Surface, Vector3};
use catallax_topology::Model;

use super::{Form, Stroke};
use crate::queries::{perpendicular, sign};

/// Points of `curve` in the plane through `origin` normal to `normal`;
/// `None` when the curve lies in the plane.
pub(super) fn in_plane(
    scratch: &mut Model,
    curve: CurveId,
    line: bool,
    origin: Point3,
    normal: Vector3,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Option<Vec<PointId>>, Refusal> {
    let ref_dir = perpendicular(&normal, dec)?;
    let plane = scratch.geom_mut().add_surface(Surface::Plane(Plane {
        origin,
        normal,
        ref_dir,
    }));
    let hits = if line {
        intersect_curve_surface_with_line_parameter(scratch.geom_mut(), curve, plane, dec)?
    } else {
        intersect_curve_surface_framed(scratch.geom_mut(), curve, plane, dec)?
    };
    marks.extend(&hits.marks);
    Ok(match hits.value {
        CurveSurface::Disjoint => Some(Vec::new()),
        CurveSurface::Contained => None,
        CurveSurface::Meets(hits) => Some(hits.into_iter().map(|h| h.point).collect()),
    })
}

/// The points of `a` over the point `q` on the sheet, when `a` lies in a
/// plane holding the view line through `q`: `a` against another plane
/// through that line, with normal `m` perpendicular to the view.
fn over(
    scratch: &mut Model,
    a: &Stroke,
    q: &Point3,
    m: Vector3,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Vec<PointId>, Refusal> {
    Ok(in_plane(scratch, a.curve, a.is_line(), q.clone(), m, marks, dec)?.unwrap_or_default())
}

/// The fold points of a conic seen edge on, where its tangent runs along
/// the view: its section by the plane through its centre normal to the view.
pub(super) fn folds(
    scratch: &mut Model,
    a: &Stroke,
    d: &Vector3,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Vec<PointId>, Refusal> {
    let Form::Conic {
        center,
        edge_on: true,
        ..
    } = &a.form
    else {
        return Ok(Vec::new());
    };
    Ok(in_plane(
        scratch,
        a.curve,
        false,
        center.clone(),
        d.clone(),
        marks,
        dec,
    )?
    .unwrap_or_default())
}

/// Where `b`'s projection ends: its end points and, seen edge on, its folds.
fn outline_ends(
    scratch: &mut Model,
    b: &Stroke,
    d: &Vector3,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Vec<Point3>, Refusal> {
    let mut points: Vec<Point3> = b
        .ends
        .into_iter()
        .flatten()
        .map(|p| scratch.geom().position(p))
        .collect();
    for fold in folds(scratch, b, d, marks, dec)? {
        points.push(scratch.geom().position(fold));
    }
    Ok(points)
}

/// Exact cut points on `a` where its projection meets `b`'s, or `None` when
/// the pair needs the binary64 search. Points may lie outside either span;
/// the caller filters them.
pub(super) fn cuts(
    scratch: &mut Model,
    a: &Stroke,
    b: &Stroke,
    d: &Vector3,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Option<Vec<PointId>>, Refusal> {
    if matches!(a.form, Form::Spline) {
        return Ok(None);
    }
    match &b.form {
        Form::Line { origin, dir } => {
            let normal = dir.cross(d);
            match in_plane(
                scratch,
                a.curve,
                a.is_line(),
                origin.clone(),
                normal.clone(),
                marks,
                dec,
            )? {
                Some(points) => Ok(Some(points)),
                None => {
                    let m = normal.cross(d);
                    let mut points = Vec::new();
                    for q in outline_ends(scratch, b, d, marks, dec)? {
                        points.extend(over(scratch, a, &q, m.clone(), marks, dec)?);
                    }
                    Ok(Some(points))
                }
            }
        }
        Form::Conic {
            center,
            normal,
            edge_on: true,
            ..
        } => match in_plane(
            scratch,
            a.curve,
            a.is_line(),
            center.clone(),
            normal.clone(),
            marks,
            dec,
        )? {
            Some(points) => Ok(Some(points)),
            None => {
                let m = normal.cross(d);
                let mut points = Vec::new();
                for q in outline_ends(scratch, b, d, marks, dec)? {
                    points.extend(over(scratch, a, &q, m.clone(), marks, dec)?);
                }
                Ok(Some(points))
            }
        },
        Form::Conic {
            center: cb,
            normal: nb,
            radius_sq: rb,
            edge_on: false,
        } => match &a.form {
            Form::Line { .. } | Form::Conic { edge_on: true, .. } => {
                // B against A's extrusion plane, then along the view onto A.
                let (origin, normal) = match &a.form {
                    Form::Line { origin, dir } => (origin.clone(), dir.cross(d)),
                    Form::Conic { center, normal, .. } => (center.clone(), normal.clone()),
                    Form::Spline => unreachable!("a spline stroke returned above"),
                };
                let Some(on_b) =
                    in_plane(scratch, b.curve, false, origin, normal.clone(), marks, dec)?
                else {
                    return Ok(None);
                };
                let m = normal.cross(d);
                let mut points = Vec::new();
                for q in on_b {
                    let q = scratch.geom().position(q);
                    points.extend(over(scratch, a, &q, m.clone(), marks, dec)?);
                }
                Ok(Some(points))
            }
            Form::Conic {
                center: ca,
                normal: na,
                radius_sq: ra,
                edge_on: false,
            } => {
                let (Some(ra), Some(rb)) = (ra, rb) else {
                    return Ok(None);
                };
                if sign(&na.cross(nb).norm_squared(), dec)? != 0 {
                    return Ok(None);
                }
                // Slide B along the view into A's plane.
                let slide = &ca.minus(cb).dot(na) / &d.dot(na);
                let moved = cb.plus(&d.scaled(&slide));
                let w = moved.minus(ca);
                if sign(&w.norm_squared(), dec)? == 0 {
                    if sign(&(ra - rb), dec)? != 0 {
                        return Ok(Some(Vec::new()));
                    }
                    // The same circle on the sheet: cut at B's ends, slid
                    // along the view onto A, where they lie exactly.
                    let mut points = Vec::new();
                    for q in b.ends.into_iter().flatten() {
                        let q = scratch.geom().position(q).plus(&d.scaled(&slide));
                        points.push(scratch.geom_mut().add_point(PointDef::Explicit(q)));
                    }
                    return Ok(Some(points));
                }
                // |x - ca|^2 - ra = |x - moved|^2 - rb, linear in x:
                // 2 x . w = |moved|^2 - |ca|^2 + ra - rb.
                let origin_vector = |p: &Point3| p.minus(&Point3::origin());
                let (ma, mb) = (origin_vector(ca), origin_vector(&moved));
                let level = &(&(&mb.norm_squared() - &ma.norm_squared()) + ra) - rb;
                let along = &(&level - &(&Real::from_i64(2) * &ma.dot(&w)))
                    / &(&Real::from_i64(2) * &w.norm_squared());
                let origin = ca.plus(&w.scaled(&along));
                Ok(Some(
                    in_plane(scratch, a.curve, false, origin, w, marks, dec)?.unwrap_or_default(),
                ))
            }
            _ => Ok(None),
        },
        Form::Spline => Ok(None),
    }
}
