//! Contour generators of a ring torus seen obliquely: neither along its axis
//! nor across it.
//!
//! On the torus `c + (R + r cos v) e(u) + r sin v a`, with `a` the unit axis
//! and `e(u)` the unit radial direction at angle `u`, the normal is
//! `cos v e(u) + sin v a` (do Carmo, "Differential Geometry of Curves and
//! Surfaces", section 2-2, the torus). The generator `N . d = 0` therefore has
//! a closed form at every angle:
//!
//! ```text
//! (cos v, sin v) = s (a.d, -(e(u).d)) / sqrt((a.d)^2 + (e(u).d)^2),  s = +1 or -1
//! ```
//!
//! With `a.d` nonzero the root never vanishes, so each sign `s` gives one
//! smooth closed branch over a full turn of `u`. The branches are degree-8
//! curves, not conics, so they are fitted (`fit.rs`) and marked approximate
//! (decision Q14). OpenCASCADE has no closed form for tori either: its
//! `Contap_ContAna` covers spheres, cylinders and cones only, a torus goes
//! to the marching path of `Contap_Contour`, and `HLRTopoBRep_DSFiller`
//! approximates the marched points with `BRepApprox_Approx`. Here the closed
//! form replaces the marching; the fit and its measured deviation follow.
//!
//! Trimming stays exact. A boundary circle of a torus face is a parallel
//! (centred on the axis, normal to it) or a meridian (in a plane through the
//! axis). On a parallel of radius `rho` at height `h`, the generator points
//! are where the circle meets the plane `(p - c).d = -h R (a.d) / (rho - R)`;
//! on a meridian with tube centre `m`, where it meets the plane
//! `(p - m).d = 0`. Both follow from `N . d = 0` with `N` proportional to `p`
//! minus its nearest core-circle point. The exact Boolean engine's
//! curve-plane intersection gives those points. Each lies on the branch
//! `s = sign(rho - R) sign(a.d)`. Cuts are ordered by exact angle around the
//! axis; exact face membership of one exact generator point inside each
//! piece keeps or drops it, as for the quadrics. Other boundary curves refuse
//! with `app.contour-trim`.

use catallax_boolean::face_contains;
use catallax_intersect::{CurveForm, CurveSurface, curve_form, intersect_curve_surface};
use catallax_numbers::{Decider, MarkLog, Rational, Real, Refusal};
use catallax_shapes::{
    BSplineCurve, Circle, Curve, Plane, Point3, PointId, Surface, Torus, Vector3,
};
use catallax_topology::numerical::{gap, tolerance};
use catallax_topology::{FaceId, Loop, Model};

use super::fit::{self, Jet};
use super::float::{cross, dot, double, rational, unit, unit_direction};
use super::{ContourArc, ContourFit, same_point};
use crate::queries::{perpendicular, sign};

fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == 0)
}

/// The torus and view in binary64, relative to the torus centre.
struct Frame {
    /// Unit seam direction, `e(0)`.
    reference: [f64; 3],
    /// Unit `a x reference`, `e(pi / 2)`.
    across: [f64; 3],
    axis: [f64; 3],
    view: [f64; 3],
    /// `a . d` for the unit vectors.
    alpha: f64,
    major: f64,
    minor: f64,
}

impl Frame {
    /// The point of branch `s` at angle `u` and its derivative in `u`,
    /// relative to the centre. With `b = e.d`, `n^2 = alpha^2 + b^2`:
    /// `(cos v)' = -s alpha b b' / n^3` and `(sin v)' = -s alpha^2 b' / n^3`.
    fn jet(&self, s: f64, u: f64) -> Jet {
        let (sin_u, cos_u) = (libm::sin(u), libm::cos(u));
        let e: [f64; 3] =
            std::array::from_fn(|k| cos_u * self.reference[k] + sin_u * self.across[k]);
        let de: [f64; 3] =
            std::array::from_fn(|k| -sin_u * self.reference[k] + cos_u * self.across[k]);
        let b = dot(e, self.view);
        let db = dot(de, self.view);
        let n2 = self.alpha * self.alpha + b * b;
        let n = n2.sqrt();
        let n3 = n2 * n;
        let cos_v = s * self.alpha / n;
        let sin_v = -s * b / n;
        let dcos_v = -s * self.alpha * b * db / n3;
        let dsin_v = -s * self.alpha * self.alpha * db / n3;
        let rho = self.major + self.minor * cos_v;
        (
            std::array::from_fn(|k| rho * e[k] + self.minor * sin_v * self.axis[k]),
            std::array::from_fn(|k| {
                self.minor * dcos_v * e[k] + rho * de[k] + self.minor * dsin_v * self.axis[k]
            }),
        )
    }
}

/// One place where the face boundary crosses a branch.
struct Cut {
    point: PointId,
    /// Exact position around the axis, scaled per axis by positive
    /// constants: `x` along the seam direction, `y` along `a x seam`.
    x: Real,
    y: Real,
    /// Angle from the seam direction in `[0, 2 pi)`.
    angle: f64,
}

/// Exact data of the torus and view.
struct Exact<'a> {
    torus: &'a Torus,
    view: &'a Vector3,
    /// `a x seam`.
    across: Vector3,
    axis_len: Real,
}

impl Exact<'_> {
    /// The part of `p - c` perpendicular to the axis.
    fn radial(&self, p: &Point3) -> Vector3 {
        let t = self.torus;
        let offset = p.minus(&t.center);
        offset.minus(
            &t.axis
                .scaled(&(offset.dot(&t.axis) / t.axis.norm_squared())),
        )
    }

    /// Position around the axis of a radial vector, as `(x, y)` with `y`
    /// rescaled so that both share the factor `|seam|`.
    fn around(&self, radial: &Vector3) -> (Real, Real) {
        (
            radial.dot(&self.torus.ref_dir),
            radial.dot(&self.across) / self.axis_len.clone(),
        )
    }

    /// The exact generator point of branch `s` in the radial direction
    /// `radial`.
    fn point(&self, s: i64, radial: &Vector3) -> Point3 {
        let t = self.torus;
        let e = radial.scaled(&(Real::one() / radial.norm_squared().sqrt()));
        let alpha = t.axis.dot(self.view) / self.axis_len.clone();
        let beta = e.dot(self.view);
        let n = (alpha.square() + beta.square()).sqrt();
        let s = Real::from_i64(s);
        let cos_v = &(&s * &alpha) / &n;
        let sin_v = -(&(&s * &beta) / &n);
        let axis = t.axis.scaled(&(Real::one() / self.axis_len.clone()));
        t.center
            .plus(&e.scaled(&(&t.major + &(&t.minor * &cos_v))))
            .plus(&axis.scaled(&(&t.minor * &sin_v)))
    }
}

/// Whether direction `p` comes strictly before `q` turning anticlockwise from
/// `from`, with `from` itself first. Exact.
fn before(
    from: (&Real, &Real),
    p: (&Real, &Real),
    q: (&Real, &Real),
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let cross = |a: (&Real, &Real), b: (&Real, &Real)| &(a.0 * b.1) - &(a.1 * b.0);
    let dot = |a: (&Real, &Real), b: (&Real, &Real)| &(a.0 * b.0) + &(a.1 * b.1);
    // Half 0 holds angles in [0, pi) from `from`, half 1 the rest.
    let half = |v: (&Real, &Real), dec: &mut Decider| -> Result<u8, Refusal> {
        let turn = sign(&cross(from, v), dec)?;
        Ok(u8::from(
            !(turn > 0 || (turn == 0 && sign(&dot(from, v), dec)? > 0)),
        ))
    };
    let (hp, hq) = (half(p, dec)?, half(q, dec)?);
    if hp != hq {
        return Ok(hp < hq);
    }
    Ok(sign(&cross(p, q), dec)? > 0)
}

/// The plane whose meeting with a boundary circle gives the circle's
/// generator points, or `None` for a circle the generator never meets.
fn cut_plane(
    exact: &Exact<'_>,
    (center, normal, radius_sq): (&Point3, &Vector3, &Real),
    dec: &mut Decider,
) -> Result<Option<(Point3, Vector3)>, Refusal> {
    let t = exact.torus;
    let d = exact.view;
    let offset = center.minus(&t.center);
    let parallel = zero(&normal.cross(&t.axis).norm_squared(), dec)?
        && zero(&offset.cross(&t.axis).norm_squared(), dec)?;
    if parallel {
        // The circle lies on the torus, so (rho - R)^2 + H^2 = r^2 with H the
        // signed height; rho - R takes the sign of rho^2 - R^2. With
        // h = H |a|, the plane is (p - c).d = -h (a.d) R / ((rho - R) |a|^2).
        let side = sign(&(radius_sq - &t.major.square()), dec)?;
        if side == 0 {
            // The top and bottom circles, where N = +-a and a.d is nonzero.
            return Ok(None);
        }
        let h = offset.dot(&t.axis);
        let root = (&t.minor.square() - &(&h.square() / &t.axis.norm_squared())).sqrt();
        let rho_gap = &Real::from_i64(i64::from(side)) * &root;
        let level = -(&(&(&h * &t.axis.dot(d)) * &t.major) / &(&rho_gap * &t.axis.norm_squared()));
        let origin = t.center.plus(&d.scaled(&(&level / &d.norm_squared())));
        return Ok(Some((origin, d.clone())));
    }
    let meridian = zero(&normal.dot(&t.axis), dec)? && zero(&offset.dot(normal), dec)?;
    if meridian {
        return Ok(Some((center.clone(), d.clone())));
    }
    Err(Refusal::new(
        "app.contour-trim",
        "a torus boundary circle that is neither a parallel nor a meridian is not supported",
    ))
}

/// Every boundary crossing of the face with each branch, by branch sign
/// `+1` then `-1`.
fn cuts(
    scratch: &mut Model,
    face: FaceId,
    exact: &Exact<'_>,
    marks: &mut MarkLog,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<[Vec<Cut>; 2], Refusal> {
    let t = exact.torus;
    let mut points: Vec<PointId> = Vec::new();
    for lp in scratch.face(face).loops.clone() {
        match scratch.loop_(lp).clone() {
            Loop::Vertex(vertex) => {
                let point = scratch.vertex(vertex).point;
                // On the generator iff ((p - c).d) |q| = R (q.d), q radial.
                let p = scratch.geom().position(point);
                let q = exact.radial(&p);
                let along = &p.minus(&t.center).dot(exact.view) * &q.norm_squared().sqrt();
                if zero(&(&along - &(&t.major * &q.dot(exact.view))), dec)? {
                    points.push(point);
                }
            }
            Loop::Coedges(coedges) => {
                for coedge in coedges {
                    let curve = scratch.edge(scratch.coedge(coedge).edge).curve;
                    let CurveForm::Circle {
                        center,
                        normal,
                        radius_sq,
                        ref_dir,
                    } = curve_form(scratch.geom(), curve, dec)?
                    else {
                        return Err(Refusal::new(
                            "app.contour-trim",
                            "a torus face bounded by a curve other than a circle is not supported",
                        ));
                    };
                    let Some((origin, plane_normal)) =
                        cut_plane(exact, (&center, &normal, &radius_sq), dec)?
                    else {
                        continue;
                    };
                    let circle = Circle {
                        center,
                        normal,
                        radius_sq,
                        ref_dir,
                    };
                    let normal = plane_normal;
                    let ref_dir = perpendicular(&normal, dec)?;
                    let plane = scratch.geom_mut().add_surface(Surface::Plane(Plane {
                        origin,
                        normal,
                        ref_dir,
                    }));
                    // The boundary circle by its exact form: an edge stored as
                    // an intersection component would take the general path.
                    let circle = scratch.geom_mut().add_curve(Curve::Circle(circle));
                    let hits = intersect_curve_surface(scratch.geom_mut(), circle, plane, dec)?;
                    marks.extend(&hits.marks);
                    match hits.value {
                        CurveSurface::Disjoint => {}
                        CurveSurface::Meets(hits) => {
                            points.extend(hits.into_iter().map(|h| h.point))
                        }
                        CurveSurface::Contained => {
                            return Err(Refusal::new(
                                "app.contour-degenerate",
                                "a torus boundary circle lies in its generator plane",
                            ));
                        }
                    }
                }
            }
        }
    }
    let ad = sign(&t.axis.dot(exact.view), dec)?;
    let mut branches: [Vec<Cut>; 2] = [Vec::new(), Vec::new()];
    for point in points {
        dec.check_deadline()?;
        let p = scratch.geom().position(point);
        let q = exact.radial(&p);
        let outer = sign(&(&q.norm_squared() - &t.major.square()), dec)?;
        let branch = &mut branches[usize::from(outer * ad < 0)];
        let mut seen = false;
        for kept in branch.iter() {
            if same_point(scratch, point, kept.point, dec)? {
                seen = true;
                break;
            }
        }
        if seen {
            continue;
        }
        let (x, y) = exact.around(&q);
        let angle = libm::atan2(double(&y, floor, dec)?, double(&x, floor, dec)?)
            .rem_euclid(std::f64::consts::TAU);
        branch.push(Cut { point, x, y, angle });
    }
    // Exact anticlockwise order from the seam direction.
    let seam = (Real::one(), Real::zero());
    for branch in &mut branches {
        for i in 1..branch.len() {
            let mut j = i;
            while j > 0
                && before(
                    (&seam.0, &seam.1),
                    (&branch[j].x, &branch[j].y),
                    (&branch[j - 1].x, &branch[j - 1].y),
                    dec,
                )?
            {
                branch.swap(j, j - 1);
                j -= 1;
            }
        }
    }
    Ok(branches)
}

/// Fit settings shared by every arc of one torus face.
struct Settings {
    frame: Frame,
    target_mm: f64,
    floor_mm: f64,
}

/// Rule T1's floor and rule T2's target for a torus fitted in a frame
/// centred on it: the coordinate size is `R + r`, and the tube radius is the
/// smallest feature the fit runs along.
fn settings(exact: &Exact<'_>, dec: &mut Decider) -> Result<(Settings, Rational), Refusal> {
    let t = exact.torus;
    // An upper bound within rule T1's fixed floor of an irrational size.
    let upper = |value: &Real, dec: &mut Decider| -> Result<Rational, Refusal> {
        Ok(match value.rational_ref() {
            Some(r) => r.clone(),
            None => dec.enclose(value, &tolerance::floor_minimum_mm())?.hi(),
        })
    };
    let size = upper(&(&t.major + &t.minor), dec)?;
    let tube = upper(&t.minor, dec)?;
    let floor = tolerance::floor_mm(&size);
    let target = tolerance::own_fitted(&floor, Some(&tube)).target_mm;
    let axis = unit_direction(&t.axis, &floor, dec)?;
    let reference = unit_direction(&t.ref_dir, &floor, dec)?;
    let view = unit_direction(exact.view, &floor, dec)?;
    let frame = Frame {
        across: unit(cross(axis, reference)),
        reference,
        axis,
        alpha: dot(axis, view),
        view,
        major: double(&t.major, &floor, dec)?,
        minor: double(&t.minor, &floor, dec)?,
    };
    Ok((
        Settings {
            frame,
            target_mm: target.to_f64_nearest(),
            floor_mm: floor.to_f64_nearest().next_up(),
        },
        floor,
    ))
}

/// Fits branch `s` over `[from, to]` and records the arc.
#[allow(clippy::too_many_arguments)]
fn push_fitted(
    scratch: &mut Model,
    face: FaceId,
    exact: &Exact<'_>,
    settings: &Settings,
    s: i64,
    (from, to): (f64, f64),
    ends: (Option<PointId>, Option<PointId>),
    arcs: &mut Vec<ContourArc>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let frame = &settings.frame;
    let sf = s as f64;
    let first = frame.jet(sf, from);
    // A whole branch closes on its own first point.
    let last = if ends == (None, None) {
        first
    } else {
        frame.jet(sf, to)
    };
    let fitted = fit::fit(
        |u| frame.jet(sf, u),
        (from, first),
        (to, last),
        settings.target_mm,
        settings.floor_mm,
        gap::EDGE_WORK,
        dec,
    )?;
    let mut knots = Vec::with_capacity(2 * fitted.breaks.len() + 4);
    for (i, &u) in fitted.breaks.iter().enumerate() {
        let copies = if i == 0 || i + 1 == fitted.breaks.len() {
            4
        } else {
            2
        };
        let u = rational(u)?;
        knots.extend(std::iter::repeat_n(u, copies));
    }
    let centre = &exact.torus.center;
    let mut controls = Vec::with_capacity(fitted.controls.len());
    for local in &fitted.controls {
        controls.push(Point3::new(
            &centre.x + &Real::from(rational(local[0])?),
            &centre.y + &Real::from(rational(local[1])?),
            &centre.z + &Real::from(rational(local[2])?),
        ));
    }
    let weights = vec![Real::one(); controls.len()];
    let curve = BSplineCurve::new(3, knots, controls, weights, dec)?;
    let curve = scratch.geom_mut().add_curve(Curve::BSpline(curve));
    arcs.push(ContourArc {
        face,
        curve,
        start: ends.0,
        end: ends.1,
        on_faces: vec![face],
        fit: Some(ContourFit {
            deviation_mm: fitted.deviation_mm,
            tolerance_mm: (gap::MARGIN * fitted.deviation_mm).next_up(),
            target_mm: settings.target_mm,
            spans: fitted.breaks.len() - 1,
        }),
    });
    Ok(())
}

/// The fitted, exactly trimmed generators of one ring-torus face for a view
/// direction with `a.d` and `a x d` both nonzero.
pub(super) fn oblique_generators(
    scratch: &mut Model,
    face: FaceId,
    torus: &Torus,
    view: &Vector3,
    arcs: &mut Vec<ContourArc>,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let exact = Exact {
        torus,
        view,
        across: torus.axis.cross(&torus.ref_dir),
        axis_len: torus.axis.norm_squared().sqrt(),
    };
    let (settings, floor) = settings(&exact, dec)?;
    let branches = cuts(scratch, face, &exact, marks, &floor, dec)?;
    let tau = std::f64::consts::TAU;
    for (cuts, s) in branches.into_iter().zip([1_i64, -1]) {
        dec.check_deadline()?;
        let m = cuts.len();
        // Pieces between consecutive cuts, each with its angular extent; a
        // branch with no cut is one piece from the seam all the way round.
        let extents: Vec<f64> = if m < 2 {
            vec![tau]
        } else {
            (0..m)
                .map(|i| (cuts[(i + 1) % m].angle - cuts[i].angle).rem_euclid(tau))
                .collect()
        };
        let mut inside = Vec::with_capacity(extents.len());
        for (i, &extent) in extents.iter().enumerate() {
            let start = if m == 0 { 0.0 } else { cuts[i].angle };
            let middle = start + extent / 2.0;
            let (sin_m, cos_m) = (libm::sin(middle), libm::cos(middle));
            // An exact radial direction at that angle: x along the seam and
            // y along a x seam, sharing the factor |seam|.
            let radial = torus.ref_dir.scaled(&Real::from(rational(cos_m)?)).plus(
                &exact
                    .across
                    .scaled(&(&Real::from(rational(sin_m)?) / &exact.axis_len)),
            );
            if m > 0 {
                let (x, y) = exact.around(&radial);
                let next = &cuts[(i + 1) % m];
                let from = (&cuts[i].x, &cuts[i].y);
                let strictly_after = before(from, from, (&x, &y), dec)?;
                let strictly_before = m == 1 || before(from, (&x, &y), (&next.x, &next.y), dec)?;
                if !(strictly_after && strictly_before) {
                    return Err(Refusal::new(
                        "app.contour-cuts",
                        "torus generator cuts lie closer than binary64 angles separate",
                    ));
                }
            }
            let probe = exact.point(s, &radial);
            inside.push(face_contains(scratch, face, &probe, dec)?);
        }
        let n = extents.len();
        if inside.iter().all(|&kept| kept) {
            let from = if m == 0 { 0.0 } else { cuts[0].angle };
            push_fitted(
                scratch,
                face,
                &exact,
                &settings,
                s,
                (from, from + tau),
                (None, None),
                arcs,
                dec,
            )?;
            continue;
        }
        // Runs of kept pieces, starting after a dropped one.
        let first = inside
            .iter()
            .position(|&kept| !kept)
            .map_or(0, |i| (i + 1) % n);
        let mut run: Option<(usize, f64, f64)> = None;
        for step in 0..n {
            let i = (first + step) % n;
            if inside[i] {
                run = Some(match run {
                    None => (i, cuts[i].angle, cuts[i].angle + extents[i]),
                    Some((s0, a, b)) => (s0, a, b + extents[i]),
                });
            }
            if (!inside[i] || step + 1 == n)
                && let Some((s0, a, b)) = run.take()
            {
                // The run ends at the cut after its last kept piece.
                let last = cuts[(i + usize::from(inside[i])) % m].point;
                push_fitted(
                    scratch,
                    face,
                    &exact,
                    &settings,
                    s,
                    (a, b),
                    (Some(cuts[s0].point), Some(last)),
                    arcs,
                    dec,
                )?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::RVector3;

    fn v(x: i64, y: i64, z: i64) -> Vector3 {
        RVector3::new(x.into(), y.into(), z.into()).to_real()
    }

    fn circle(center: Vector3, normal: Vector3, radius_sq: i64) -> (Point3, Vector3, Real) {
        (
            Point3::origin().plus(&center),
            normal,
            Real::from_i64(radius_sq),
        )
    }

    fn plane(
        exact: &Exact<'_>,
        (center, normal, radius_sq): (Point3, Vector3, Real),
        dec: &mut Decider,
    ) -> Result<Option<(Point3, Vector3)>, Refusal> {
        cut_plane(exact, (&center, &normal, &radius_sq), dec)
    }

    #[test]
    fn boundary_circles_are_classified_as_parallels_meridians_or_refused() {
        let torus = Torus {
            center: Point3::origin(),
            axis: v(0, 0, 1),
            major: Real::from_i64(5),
            minor: Real::from_i64(1),
            ref_dir: v(1, 0, 0),
        };
        let view = v(0, 1, 1);
        let exact = Exact {
            torus: &torus,
            across: torus.axis.cross(&torus.ref_dir),
            axis_len: Real::one(),
            view: &view,
        };
        let mut dec = Decider::default();
        // The outer equator: h = 0, so the plane passes through the centre.
        let (origin, normal) = plane(&exact, circle(v(0, 0, 0), v(0, 0, 2), 36), &mut dec)
            .unwrap()
            .unwrap();
        assert_eq!(
            sign(&origin.minus(&Point3::origin()).dot(&view), &mut dec).unwrap(),
            0
        );
        assert_eq!(
            sign(&normal.cross(&view).norm_squared(), &mut dec).unwrap(),
            0
        );
        // The top circle, where the normal is the axis: no generator point.
        assert!(
            plane(&exact, circle(v(0, 0, 1), v(0, 0, -1), 25), &mut dec)
                .unwrap()
                .is_none()
        );
        // A meridian: the plane through its tube centre normal to the view.
        let (origin, _) = plane(&exact, circle(v(0, 5, 0), v(1, 0, 0), 1), &mut dec)
            .unwrap()
            .unwrap();
        assert_eq!(
            sign(&(&origin.y - &Real::from_i64(5)), &mut dec).unwrap(),
            0
        );
        // A tilted circle is neither.
        let refusal = plane(&exact, circle(v(0, 5, 0), v(1, 0, 1), 1), &mut dec).unwrap_err();
        assert_eq!(refusal.code, "app.contour-trim");
    }

    #[test]
    fn the_closed_form_generator_has_a_normal_perpendicular_to_the_view() {
        let frame = Frame {
            reference: [1.0, 0.0, 0.0],
            across: [0.0, 1.0, 0.0],
            axis: [0.0, 0.0, 1.0],
            view: unit([1.0, 2.0, 2.0]),
            alpha: 2.0 / 3.0,
            major: 5.0,
            minor: 1.0,
        };
        for s in [1.0, -1.0] {
            for k in 0..64 {
                let u = k as f64 * std::f64::consts::TAU / 64.0;
                let (p, dp) = frame.jet(s, u);
                // The tube centre is R e(u); the normal points from it to p.
                let core = [5.0 * libm::cos(u), 5.0 * libm::sin(u), 0.0];
                let normal = [p[0] - core[0], p[1] - core[1], p[2] - core[2]];
                assert!(dot(normal, frame.view).abs() < 1e-14);
                assert!((dot(normal, normal) - 1.0).abs() < 1e-14);
                // The derivative matches a central difference.
                let h = 1e-6;
                let (a, _) = frame.jet(s, u + h);
                let (b, _) = frame.jet(s, u - h);
                for axis in 0..3 {
                    assert!(((a[axis] - b[axis]) / (2.0 * h) - dp[axis]).abs() < 1e-8);
                }
            }
        }
    }
}
