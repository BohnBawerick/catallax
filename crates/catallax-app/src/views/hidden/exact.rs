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

use catallax_intersect::intersect_curve_plane_points;
use catallax_numbers::{Decider, MarkLog, Real, Refusal};
use catallax_shapes::{Meet, MeetOf, Plane, Point3, PointDef, PointId, Position, Surface, Vector3};
use catallax_topology::Model;

use super::{Form, Stroke};
use crate::queries::{perpendicular, sign};

/// Intersections with the cutting plane. A finite exact line only needs
/// interior cuts because its endpoints already split the view.
fn in_plane(
    scratch: &mut Model,
    stroke: &Stroke,
    origin: Point3,
    normal: Vector3,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Option<Vec<PointId>>, Refusal> {
    let exact_line = match (&stroke.form, stroke.ends) {
        (
            Form::Line {
                origin: line_origin,
                dir,
            },
            Some(ends),
        ) if stroke.exact_line_trim
            && (line_origin.as_rational().is_none()
                || dir.as_rational().is_none()
                || origin.as_rational().is_none()
                || normal.as_rational().is_none()) =>
        {
            Some((dir, ends))
        }
        _ => None,
    };
    let (origin, normal) = if exact_line.is_some() {
        let [nx, ny, nz] = reduced([&normal.x, &normal.y, &normal.z], dec)?;
        (reduced_point(&origin, dec)?, Vector3::new(nx, ny, nz))
    } else {
        (origin, normal)
    };
    let ref_dir = perpendicular(&normal, dec)?;
    let plane = Plane {
        origin,
        normal,
        ref_dir,
    };
    if let Some((dir, [from, to])) = exact_line {
        let from = reduced_point(&scratch.geom().position(from), dec)?;
        let to = reduced_point(&scratch.geom().position(to), dec)?;
        let [a, b] = [&from, &to].map(|p| p.minus(&plane.origin).dot(&plane.normal));
        let [a_side, b_side] = [sign(&a, dec)?, sign(&b, dec)?];
        if a_side == 0 && b_side == 0 {
            return Ok(None);
        }
        if a_side == b_side || a_side == 0 || b_side == 0 {
            return Ok(Some(Vec::new()));
        }
        let [dx, dy, dz] = reduced([&dir.x, &dir.y, &dir.z], dec)?;
        let point = axial_cut(&from, &Vector3::new(dx, dy, dz), &plane)
            .unwrap_or_else(|| from.plus(&to.minus(&from).scaled(&(-&a / (&b - &a)))));
        let geom = scratch.geom_mut();
        let surface = geom.add_surface(Surface::Plane(plane));
        let cut = geom.add_point(PointDef::Meet(Meet {
            of: MeetOf::CurveSurface(stroke.curve, surface),
            position: Position::Closed(point),
        }));
        return Ok(Some(vec![cut]));
    }
    if let Form::OpenConic(open) = &stroke.form {
        let points = open.plane_points(&plane.origin, &plane.normal, dec)?;
        return Ok(points.map(|points| {
            let geom = scratch.geom_mut();
            let surface = geom.add_surface(Surface::Plane(plane));
            points
                .into_iter()
                .map(|point| {
                    geom.add_point(PointDef::Meet(Meet {
                        of: MeetOf::CurveSurface(stroke.curve, surface),
                        position: Position::Closed(point),
                    }))
                })
                .collect()
        }));
    }
    let hits = intersect_curve_plane_points(scratch.geom_mut(), stroke.curve, plane, dec)?;
    marks.extend(&hits.marks);
    Ok(hits.value)
}

fn reduced(values: [&Real; 3], dec: &mut Decider) -> Result<[Real; 3], Refusal> {
    let mut result = values.map(Clone::clone);
    for value in &mut result {
        if value.rational_ref().is_none()
            && let Some(rational) = dec.rational_value(value)?
        {
            *value = rational.into();
        }
    }
    Ok(result)
}

fn reduced_point(p: &Point3, dec: &mut Decider) -> Result<Point3, Refusal> {
    let [x, y, z] = reduced([&p.x, &p.y, &p.z], dec)?;
    Ok(Point3::new(x, y, z))
}

/// Solve for the one free coordinate of an axial line. Its origin along
/// that axis and the scale of its direction cannot affect the solution.
fn axial_cut(from: &Point3, dir: &Vector3, plane: &Plane) -> Option<Point3> {
    let direction = [&dir.x, &dir.y, &dir.z];
    let axis = (0..3).find(|&axis| {
        direction
            .iter()
            .enumerate()
            .all(|(i, value)| i == axis || value.rational_ref().is_some_and(|q| q.is_zero()))
    })?;
    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
    let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
    let mut coords = [from.x.clone(), from.y.clone(), from.z.clone()];
    let offset = (0..3).filter(|&i| i != axis).fold(Real::zero(), |sum, i| {
        sum + &(origin[i] - &coords[i]) * normal[i]
    });
    coords[axis] = origin[axis] + &(&offset / normal[axis]);
    let [x, y, z] = coords;
    Some(Point3::new(x, y, z))
}

/// Factor out a line endpoint whose projection lies on the circle. Solving
/// for it again hides its equality to the endpoint behind a square root.
fn circle_at_line_end(
    scratch: &mut Model,
    line: &Stroke,
    center: &Point3,
    normal: &Vector3,
    radius_sq: &Real,
    view: &Vector3,
    dec: &mut Decider,
) -> Result<Option<Vec<PointId>>, Refusal> {
    let Some([from, to]) = line.ends.filter(|_| line.exact_line_trim) else {
        return Ok(None);
    };
    let from = reduced_point(&scratch.geom().position(from), dec)?;
    let to = reduced_point(&scratch.geom().position(to), dec)?;
    let rate = view.dot(normal);
    let project = |p: &Point3| {
        let offset = p.minus(center);
        offset.minus(&view.scaled(&(offset.dot(normal) / rate.clone())))
    };
    let [a, b] = [&from, &to].map(project);
    let ends = [
        sign(&(&a.norm_squared() - radius_sq), dec)?,
        sign(&(&b.norm_squared() - radius_sq), dec)?,
    ];
    let (start, end, at, other) = match ends {
        [0, 0] => return Ok(Some(Vec::new())),
        [0, _] => (from, to, a, b),
        [_, 0] => (to, from, b, a),
        _ => return Ok(None),
    };
    let step = other.minus(&at);
    let quadratic = step.norm_squared();
    let advance = -Real::from_i64(2) * at.dot(&step);
    if sign(&advance, dec)? <= 0 || sign(&(&quadratic - &advance), dec)? <= 0 {
        return Ok(Some(Vec::new()));
    }
    let point = start.plus(&end.minus(&start).scaled(&(&advance / &quadratic)));
    Ok(Some(vec![
        scratch.geom_mut().add_point(PointDef::Explicit(point)),
    ]))
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
    Ok(in_plane(scratch, a, q.clone(), m, marks, dec)?.unwrap_or_default())
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
    if let Form::OpenConic(open) = &a.form {
        return if open.edge_on {
            Ok(open
                .folds(d, dec)?
                .into_iter()
                .map(|p| scratch.geom_mut().add_point(PointDef::Explicit(p)))
                .collect())
        } else {
            Ok(Vec::new())
        };
    }
    let Form::Conic {
        center,
        edge_on: true,
        ..
    } = &a.form
    else {
        return Ok(Vec::new());
    };
    Ok(in_plane(scratch, a, center.clone(), d.clone(), marks, dec)?.unwrap_or_default())
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

/// A plane extruding a projected line along the view.
fn sheet_plane(form: &Form, d: &Vector3) -> Option<(Point3, Vector3)> {
    match form {
        Form::Line { origin, dir } => Some((origin.clone(), dir.cross(d))),
        Form::Conic {
            center,
            normal,
            edge_on: true,
            ..
        } => Some((center.clone(), normal.clone())),
        Form::OpenConic(open) if open.edge_on => Some((open.origin.clone(), open.normal.clone())),
        _ => None,
    }
}

/// B against A's extrusion plane, then along the view onto A.
fn conic_cuts(
    scratch: &mut Model,
    a: &Stroke,
    b: &Stroke,
    d: &Vector3,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Option<Vec<PointId>>, Refusal> {
    let Some((origin, normal)) = sheet_plane(&a.form, d) else {
        return Ok(None);
    };
    let Some(on_b) = in_plane(scratch, b, origin, normal.clone(), marks, dec)? else {
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
    if let Some((origin, normal)) = sheet_plane(&b.form, d) {
        return match in_plane(scratch, a, origin, normal.clone(), marks, dec)? {
            Some(points) => Ok(Some(points)),
            None => {
                let m = normal.cross(d);
                let mut points = Vec::new();
                for q in outline_ends(scratch, b, d, marks, dec)? {
                    points.extend(over(scratch, a, &q, m.clone(), marks, dec)?);
                }
                Ok(Some(points))
            }
        };
    }
    match &b.form {
        Form::Conic {
            center: cb,
            normal: nb,
            radius_sq: rb,
            edge_on: false,
        } => match &a.form {
            Form::Line { .. } | Form::Conic { edge_on: true, .. } | Form::OpenConic(_)
                if sheet_plane(&a.form, d).is_some() =>
            {
                if let Some(radius_sq) = rb
                    && let Some(points) = circle_at_line_end(scratch, a, cb, nb, radius_sq, d, dec)?
                {
                    return Ok(Some(points));
                }
                conic_cuts(scratch, a, b, d, marks, dec)
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
                    in_plane(scratch, a, origin, w, marks, dec)?.unwrap_or_default(),
                ))
            }
            _ => Ok(None),
        },
        Form::OpenConic(_) => conic_cuts(scratch, a, b, d, marks, dec),
        _ => Ok(None),
    }
}
