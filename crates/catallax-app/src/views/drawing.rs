//! Drawing curves: the binary64 lines, arcs, circles and NURBS a 2D drawing
//! takes, from the exact curves of edges, face boundaries, hidden-line
//! strokes and sections.
//!
//! Every curve keeps its exact shape where a drawing format can carry it,
//! and its numbers are binary64 values of exact ones: lines and circular
//! arcs as themselves; ellipse, parabola and hyperbola arcs in the rational
//! quadratic form that represents a conic exactly (Piegl and Tiller, "The
//! NURBS Book", 2nd ed., section 7.5 and algorithm A7.1, as OpenCASCADE's
//! `Convert_EllipseToBSplineCurve` with its `TgtThetaOver2` parametrisation
//! writes them); B-spline arcs as their carrier's Bezier spans between the
//! ends (algorithm A5.6, the cut OpenCASCADE's `Geom_BSplineCurve::Segment`
//! makes by knot insertion). An arc's direction is decided exactly: from an
//! exact point inside it, never from a rounded angle.
//!
//! Intersection branches, patch sections and patch pairs have no closed form.
//! They are traced (`catallax_intersect::trace`, the certified trace the
//! STEP writer also fits) and fitted by cubic Hermite spans as the contour
//! generators are (`fit.rs`, SMLib's `SmHermiteCurve` form; OpenCASCADE
//! approximates its walked intersection lines into B-splines the same way
//! in `ApproxInt`). Such a curve carries a [`ContourFit`] and the answer a
//! [`DRAWING_FIT_MARK`] mark: decision Q14 marks it approximate.

use std::cell::RefCell;

use catallax_intersect::{CurveForm, TracePiece, arc_point, curve_form, curve_is_closed, trace};
use catallax_numbers::{
    Answer, Answered, Ball, BoundKind, Decider, Interval, Mark, MarkLog, Rational, Real, Refusal,
};
use catallax_shapes::{BSplineCurve, CurveId, Point3, PointId, Vector3};
use catallax_topology::numerical::{gap, tolerance};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

use super::ContourFit;
use super::fit::{self, Jet};
use super::float::{above, cross, dot, unit, unit_direction};
use super::hidden::geo::{Geo, Span, decompose, restrict};
use super::hidden::{HiddenLineView, ViewStroke};
use crate::ApplicationBody;
use crate::queries::{QueryRefusal, query, sign};

/// Code of the mark an answer carries when a drawing curve was fitted. Its
/// bound is the largest [`ContourFit::tolerance_mm`]: a measured value, not
/// a proof.
pub const DRAWING_FIT_MARK: &str = "app.drawing-fitted";

/// A curve for a 2D drawing, in model coordinates in millimetres.
#[derive(Clone, Debug, PartialEq)]
pub enum DrawingCurve {
    /// A straight segment from `a` to `b`.
    Line {
        /// Start point.
        a: [f64; 3],
        /// End point.
        b: [f64; 3],
    },
    /// A circular arc, counter-clockwise about `normal` from `start_deg` to
    /// `end_deg`: its point at angle `t` is
    /// `center + radius (cos t x_axis + sin t (normal x x_axis))`.
    Arc {
        /// Centre.
        center: [f64; 3],
        /// Unit axis; the arc runs counter-clockwise seen from its tip.
        normal: [f64; 3],
        /// Unit direction of angle zero, perpendicular to `normal`: the
        /// direction of the start point.
        x_axis: [f64; 3],
        /// Radius.
        radius: f64,
        /// Start angle in degrees: zero.
        start_deg: f64,
        /// End angle in degrees, above `start_deg` and at most 360 above it.
        end_deg: f64,
    },
    /// A full circle, counter-clockwise about `normal`.
    Circle {
        /// Centre.
        center: [f64; 3],
        /// Unit axis.
        normal: [f64; 3],
        /// Radius.
        radius: f64,
    },
    /// A clamped non-uniform rational B-spline: `knots` has
    /// `points.len() + degree + 1` entries; `weights` is empty when every
    /// weight is 1, otherwise one per point.
    Nurbs {
        /// Polynomial degree.
        degree: usize,
        /// Full knot vector, non-decreasing.
        knots: Vec<f64>,
        /// Control points.
        points: Vec<[f64; 3]>,
        /// Control point weights, or empty.
        weights: Vec<f64>,
    },
}

impl DrawingCurve {
    /// The same curve run the other way.
    pub fn reversed(&self) -> DrawingCurve {
        let neg = |v: &[f64; 3]| [-v[0], -v[1], -v[2]];
        match self {
            DrawingCurve::Line { a, b } => DrawingCurve::Line { a: *b, b: *a },
            // Counter-clockwise about the opposite axis from the end's
            // mirrored angle: angle `t` about `-n` is angle `-t` about `n`.
            DrawingCurve::Arc {
                center,
                normal,
                x_axis,
                radius,
                start_deg,
                end_deg,
            } => {
                let sweep = end_deg - start_deg;
                let first = (-end_deg).rem_euclid(360.0);
                DrawingCurve::Arc {
                    center: *center,
                    normal: neg(normal),
                    x_axis: *x_axis,
                    radius: *radius,
                    start_deg: first,
                    end_deg: first + sweep,
                }
            }
            DrawingCurve::Circle {
                center,
                normal,
                radius,
            } => DrawingCurve::Circle {
                center: *center,
                normal: neg(normal),
                radius: *radius,
            },
            DrawingCurve::Nurbs {
                degree,
                knots,
                points,
                weights,
            } => {
                let (lo, hi) = (knots[0], knots[knots.len() - 1]);
                DrawingCurve::Nurbs {
                    degree: *degree,
                    knots: knots.iter().rev().map(|k| lo + hi - k).collect(),
                    points: points.iter().rev().copied().collect(),
                    weights: weights.iter().rev().copied().collect(),
                }
            }
        }
    }
}

/// A drawing curve and, for a fitted one, how far it may lie from the true
/// curve.
#[derive(Clone, Debug)]
pub struct DrawnCurve {
    /// The curve.
    pub curve: DrawingCurve,
    /// `None` when the curve is the exact one in binary64.
    pub fit: Option<ContourFit>,
}

impl DrawnCurve {
    fn reversed(self) -> DrawnCurve {
        DrawnCurve {
            curve: self.curve.reversed(),
            fit: self.fit,
        }
    }
}

/// One edge of a face boundary, drawn in the direction the loop travels.
#[derive(Clone, Debug)]
pub struct DrawnEdge {
    /// The edge, in the queried model.
    pub edge: EdgeId,
    /// Its curve in travel order.
    pub drawn: DrawnCurve,
}

fn refused(detail: &str) -> Refusal {
    Refusal::new("app.drawing-curve", detail.to_string())
}

/// Source: IEEE 754 binary64 carries 53 significant bits. An enclosure 2^-60
/// of a value's size wide pins its nearest binary64 value to within one
/// unit in the last place.
const OUTPUT_BITS: i64 = 60;

/// The binary64 value of an exact number: the nearest one of a rational,
/// otherwise the middle of an enclosure narrower than its last place.
fn out(value: &Real, dec: &mut Decider) -> Result<f64, Refusal> {
    if let Some(rational) = value.rational_ref() {
        return Ok(rational.to_f64_nearest());
    }
    let approx = value.approx();
    let size = approx.lo().abs().max(approx.hi().abs());
    let width = if above(size, 0.0) && size.is_finite() {
        Rational::from_f64_exact(size)
            .ok_or_else(|| refused("a drawing coordinate is not finite"))?
            .mul_pow2(-OUTPUT_BITS)
    } else {
        Rational::one().mul_pow2(-OUTPUT_BITS)
    };
    Ok(dec.enclose(value, &width)?.mid().to_f64_nearest())
}

fn point(p: &Point3, dec: &mut Decider) -> Result<[f64; 3], Refusal> {
    Ok([out(&p.x, dec)?, out(&p.y, dec)?, out(&p.z, dec)?])
}

/// The unit binary64 direction of an exact nonzero vector.
fn direction(v: &Vector3, dec: &mut Decider) -> Result<[f64; 3], Refusal> {
    if [&v.x, &v.y, &v.z]
        .iter()
        .all(|c| c.rational_ref().is_some())
    {
        return unit_direction(v, &tolerance::floor_minimum_mm(), dec);
    }
    Ok(unit([out(&v.x, dec)?, out(&v.y, dec)?, out(&v.z, dec)?]))
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// The sweep of an arc from `a` to `b` about its centre, counter-clockwise
/// about the axis whose exact signs of `sin` and `cos` of the sweep are
/// `s` and `c`, in `(0, 2 pi)`. The binary64 value `estimate` is clamped to
/// the quadrant those signs prove.
fn sweep(s: i8, c: i8, estimate: f64) -> Result<f64, Refusal> {
    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    let [lo, hi] = match (s, c) {
        (1, 1) => [0.0, FRAC_PI_2],
        (1, 0) => [FRAC_PI_2, FRAC_PI_2],
        (1, _) => [FRAC_PI_2, PI],
        (0, -1) => [PI, PI],
        (-1, -1) => [PI, PI + FRAC_PI_2],
        (-1, 0) => [PI + FRAC_PI_2, PI + FRAC_PI_2],
        (-1, _) => [PI + FRAC_PI_2, TAU],
        _ => return Err(refused("an arc's two ends are one point")),
    };
    // The representative of the estimate, modulo a turn, nearest the
    // proven range: a sweep just short of a whole turn can round to just
    // past zero.
    let est = estimate.rem_euclid(TAU);
    let off = |x: f64| (lo - x).max(x - hi).max(0.0);
    let best = [est - TAU, est, est + TAU]
        .into_iter()
        .min_by(|a, b| off(*a).total_cmp(&off(*b)))
        .expect("three candidates");
    Ok(best.clamp(lo, hi))
}

/// How a conic arc turns: the exact axis it runs counter-clockwise about,
/// and the exact signs of the sine and cosine of its sweep in the conic's
/// own angle.
struct Turn {
    normal: Vector3,
    sin: i8,
    cos: i8,
}

/// Decides how an arc of a circle or ellipse runs from `s` to `e` through
/// the interior point `m`: counter-clockwise about `normal` exactly when
/// `m - s` turns to `e - s` that way. For an ellipse, `axes` holds its major
/// direction and squared semi-axes, and the cosine's sign is read in the
/// coordinates that make it a circle, `(u . x) / major` and
/// `(v . x) / minor` with `u` and `v` the unit axes; `None` for a circle.
fn turn(
    center: &Point3,
    normal: &Vector3,
    axes: Option<(&Vector3, &Real, &Real)>,
    [s, m, e]: [&Point3; 3],
    dec: &mut Decider,
) -> Result<Turn, Refusal> {
    let way = sign(&normal.dot(&m.minus(s).cross(&e.minus(s))), dec)?;
    let normal = match way {
        1 => normal.clone(),
        -1 => normal.scaled(&Real::from_i64(-1)),
        _ => return Err(refused("an arc's interior point is not between its ends")),
    };
    let (a, b) = (s.minus(center), e.minus(center));
    let sin = sign(&normal.dot(&a.cross(&b)), dec)?;
    let cos = match axes {
        None => sign(&a.dot(&b), dec)?,
        Some((major_dir, major_sq, minor_sq)) => {
            let minor_dir = normal.cross(major_dir);
            let along =
                a.dot(major_dir) * b.dot(major_dir) / (major_dir.norm_squared() * major_sq.clone());
            let across = a.dot(&minor_dir) * b.dot(&minor_dir)
                / (minor_dir.norm_squared() * minor_sq.clone());
            sign(&(along + across), dec)?
        }
    };
    Ok(Turn { normal, sin, cos })
}

/// A circle, or with `ends` its arc between two points through an exact
/// interior point.
fn circle(
    center: &Point3,
    normal: &Vector3,
    radius_sq: &Real,
    ends: Option<([&Point3; 2], Point3)>,
    dec: &mut Decider,
) -> Result<DrawingCurve, Refusal> {
    let c = point(center, dec)?;
    let radius = out(&radius_sq.sqrt(), dec)?;
    let Some(([s, e], m)) = ends else {
        return Ok(DrawingCurve::Circle {
            center: c,
            normal: direction(normal, dec)?,
            radius,
        });
    };
    let t = turn(center, normal, None, [s, &m, e], dec)?;
    let n = direction(&t.normal, dec)?;
    let (a, b) = (sub(point(s, dec)?, c), sub(point(e, dec)?, c));
    let estimate = dot(n, cross(a, b)).atan2(dot(a, b));
    let angle = sweep(t.sin, t.cos, estimate)?;
    Ok(DrawingCurve::Arc {
        center: c,
        normal: n,
        x_axis: unit(a),
        radius,
        start_deg: 0.0,
        end_deg: angle.to_degrees(),
    })
}

/// The rational quadratic arcs of the ellipse `c + cos t u + sin t v` from
/// angle `from` through `sweep`, at most a quarter turn each (Piegl and
/// Tiller, algorithm A7.1, which an affine map carries from the circle to
/// the ellipse): the middle control of each arc lies on its middle angle at
/// `1 / cos(half)` and weighs `cos(half)`.
fn ellipse_nurbs(
    c: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
    from: f64,
    sweep: f64,
    ends: Option<[[f64; 3]; 2]>,
) -> DrawingCurve {
    use std::f64::consts::FRAC_PI_2;
    let at = |t: f64, r: f64| add(c, add(scale(u, r * t.cos()), scale(v, r * t.sin())));
    let arcs = ((sweep / FRAC_PI_2).ceil() as usize).clamp(1, 4);
    let step = sweep / arcs as f64;
    let w = (step / 2.0).cos();
    let mut points = vec![at(from, 1.0)];
    let mut weights = vec![1.0];
    let mut knots = vec![0.0; 3];
    for i in 1..=arcs {
        let angle = from + step * i as f64;
        points.push(at(angle - step / 2.0, 1.0 / w));
        points.push(at(angle, 1.0));
        weights.extend([w, 1.0]);
        let k = i as f64 / arcs as f64;
        knots.extend(std::iter::repeat_n(k, if i == arcs { 3 } else { 2 }));
    }
    if let Some([s, e]) = ends {
        points[0] = s;
        let last = points.len() - 1;
        points[last] = e;
    }
    DrawingCurve::Nurbs {
        degree: 2,
        knots,
        points,
        weights,
    }
}

/// An ellipse, or with `ends` its arc between two points through an exact
/// interior point; a whole ellipse starts at `anchor` when it has one.
#[allow(clippy::too_many_arguments)]
fn ellipse(
    center: &Point3,
    normal: &Vector3,
    major_dir: &Vector3,
    major_sq: &Real,
    minor_sq: &Real,
    ends: Option<([&Point3; 2], Point3)>,
    anchor: Option<&Point3>,
    dec: &mut Decider,
) -> Result<DrawingCurve, Refusal> {
    let c = point(center, dec)?;
    let (a, b) = (out(&major_sq.sqrt(), dec)?, out(&minor_sq.sqrt(), dec)?);
    let u = direction(major_dir, dec)?;
    // The eccentric angle of a point in the frame `u`, `n x u`.
    let angle = |p: [f64; 3], n: [f64; 3]| {
        let d = sub(p, c);
        (dot(d, cross(n, u)) / b).atan2(dot(d, u) / a)
    };
    let Some(([s, e], m)) = ends else {
        let n = direction(normal, dec)?;
        let v = cross(n, u);
        let (from, at) = match anchor {
            Some(p) => {
                let p = point(p, dec)?;
                (angle(p, n), Some([p, p]))
            }
            None => (0.0, None),
        };
        return Ok(ellipse_nurbs(
            c,
            scale(u, a),
            scale(v, b),
            from,
            std::f64::consts::TAU,
            at,
        ));
    };
    let t = turn(
        center,
        normal,
        Some((major_dir, major_sq, minor_sq)),
        [s, &m, e],
        dec,
    )?;
    let n = direction(&t.normal, dec)?;
    let (ps, pe) = (point(s, dec)?, point(e, dec)?);
    let (from, to) = (angle(ps, n), angle(pe, n));
    let swept = sweep(t.sin, t.cos, to - from)?;
    Ok(ellipse_nurbs(
        c,
        scale(u, a),
        scale(cross(n, u), b),
        from,
        swept,
        Some([ps, pe]),
    ))
}

/// One quadratic Bezier arc from `p0` to `p2` with middle control `p1`,
/// rational when the middle control has a weight.
fn quadratic(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], w: Option<f64>) -> DrawingCurve {
    DrawingCurve::Nurbs {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        points: vec![p0, p1, p2],
        weights: w.map_or_else(Vec::new, |w| vec![1.0, w, 1.0]),
    }
}

/// The plane frame of a parabola or hyperbola: the unit axis and the unit
/// transverse direction `normal x axis`.
fn open_frame(
    normal: &Vector3,
    axis_dir: &Vector3,
    dec: &mut Decider,
) -> Result<([f64; 3], [f64; 3]), Refusal> {
    Ok((
        direction(axis_dir, dec)?,
        direction(&normal.cross(axis_dir), dec)?,
    ))
}

/// The arc of `curve`, whose closed form is `form`, from `ends[0]` to
/// `ends[1]` along its natural direction; `None` for the whole closed
/// curve, and two equal ends for the whole closed curve from that point.
fn natural(
    model: &Model,
    curve: CurveId,
    form: CurveForm,
    ends: Option<[PointId; 2]>,
    dec: &mut Decider,
) -> Answered<DrawnCurve> {
    let geom = model.geom();
    let exact = |curve| Answer::exact(DrawnCurve { curve, fit: None });
    let positions = ends.map(|[a, b]| [geom.position(a), geom.position(b)]);
    let anchored = match &positions {
        Some([a, b]) => sign(&a.minus(b).norm_squared(), dec)? == 0,
        None => false,
    };
    // An exact point strictly inside a proper arc decides its direction.
    let mut marks = MarkLog::new();
    let mut interior = |dec: &mut Decider| -> Result<Option<Point3>, Refusal> {
        let Some([a, b]) = ends else { return Ok(None) };
        if anchored {
            return Ok(None);
        }
        let inside = arc_point(geom, curve, Some(a), Some(b), dec)?;
        marks.extend(&inside.marks);
        Ok(Some(inside.value))
    };
    let drawn = match form {
        CurveForm::Line { .. } => {
            let Some([a, b]) = &positions else {
                return Err(refused("a whole line has no drawing"));
            };
            exact(DrawingCurve::Line {
                a: point(a, dec)?,
                b: point(b, dec)?,
            })
        }
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } => {
            let m = interior(dec)?;
            let arc = match (&positions, m) {
                (Some([a, b]), Some(m)) => Some(([a, b], m)),
                _ => None,
            };
            exact(circle(&center, &normal, &radius_sq, arc, dec)?)
        }
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        } => {
            let m = interior(dec)?;
            let arc = match (&positions, m) {
                (Some([a, b]), Some(m)) => Some(([a, b], m)),
                _ => None,
            };
            let anchor = positions.as_ref().filter(|_| anchored).map(|p| &p[0]);
            exact(ellipse(
                &center, &normal, &major_dir, &major_sq, &minor_sq, arc, anchor, dec,
            )?)
        }
        CurveForm::Parabola {
            vertex,
            normal,
            axis_dir,
            focal,
        } => {
            let Some([a, b]) = &positions else {
                return Err(refused("a whole parabola has no drawing"));
            };
            // On `y^2 = 4 f x` the tangents at `y0` and `y1` meet at
            // `(y0 y1 / (4 f), (y0 + y1) / 2)`; a parabola arc is a
            // polynomial quadratic Bezier with that middle control (Farin,
            // "Curves and Surfaces for CAGD", 5th ed., section 12.2).
            let (x, y) = open_frame(&normal, &axis_dir, dec)?;
            let v = point(&vertex, dec)?;
            let f = out(&focal, dec)?;
            let (p0, p2) = (point(a, dec)?, point(b, dec)?);
            let (y0, y1) = (dot(sub(p0, v), y), dot(sub(p2, v), y));
            let p1 = add(
                v,
                add(scale(x, y0 * y1 / (4.0 * f)), scale(y, (y0 + y1) / 2.0)),
            );
            exact(quadratic(p0, p1, p2, None))
        }
        CurveForm::Hyperbola {
            center,
            normal,
            axis_dir,
            transverse_sq,
            conjugate_sq,
            branch,
        } => {
            let Some([a, b]) = &positions else {
                return Err(refused("a whole hyperbola has no drawing"));
            };
            // On `(s a cosh p, b sinh p)` the arc from `p0` to `p1` has its
            // middle control at `(s a cosh m, b sinh m) / cosh h` with weight
            // `cosh h`, `m` the middle and `h` the half span: a hyperbolic
            // rotation carries it to the symmetric arc, whose tangents meet
            // on the axis at `a / cosh h` and whose vertex fixes the weight
            // (Piegl and Tiller, section 7.3, the conic shoulder point).
            let (x, y) = open_frame(&normal, &axis_dir, dec)?;
            let c = point(&center, dec)?;
            let (ta, tb) = (
                out(&transverse_sq.sqrt(), dec)?,
                out(&conjugate_sq.sqrt(), dec)?,
            );
            let s = match branch {
                catallax_numbers::StrictSign::Positive => 1.0,
                catallax_numbers::StrictSign::Negative => -1.0,
            };
            let (p0, p2) = (point(a, dec)?, point(b, dec)?);
            let (h0, h1) = (
                (dot(sub(p0, c), y) / tb).asinh(),
                (dot(sub(p2, c), y) / tb).asinh(),
            );
            let (mid, half) = ((h0 + h1) / 2.0, (h1 - h0).abs() / 2.0);
            let w = half.cosh();
            let p1 = add(
                c,
                add(
                    scale(x, s * ta * mid.cosh() / w),
                    scale(y, tb * mid.sinh() / w),
                ),
            );
            exact(quadratic(p0, p1, p2, Some(w)))
        }
        CurveForm::BSpline(spline) => {
            let closed = curve_is_closed(geom, curve, dec)?;
            exact(spline_arc(
                &spline,
                positions.as_ref(),
                anchored,
                closed,
                dec,
            )?)
        }
        CurveForm::Branch | CurveForm::PatchSection(_) | CurveForm::PatchPair(_) => {
            let ends = ends.map(|[a, b]| if anchored { [a, a] } else { [a, b] });
            let drawn = traced(model, curve, ends, dec)?;
            marks.extend(&drawn.marks);
            if let Some(fit) = &drawn.value.fit {
                marks.push(fit_mark(fit)?);
            }
            Answer::exact(drawn.value)
        }
        _ => return Err(refused("this curve has no drawing conversion")),
    };
    Ok(Answer::with_marks(drawn.value, marks))
}

pub(in crate::views) fn fit_mark(fit: &ContourFit) -> Result<Mark, Refusal> {
    let bound = Rational::from_f64_exact(fit.tolerance_mm)
        .ok_or_else(|| refused("a fitted drawing curve has a nonfinite tolerance"))?;
    Ok(Mark::new(
        DRAWING_FIT_MARK,
        BoundKind::Distance,
        bound,
        "fitted drawing curve: rule T4 tolerance of the largest measured deviation, not a proof",
    ))
}

/// The arc of a B-spline carrier between `ends`, as its Bezier spans there.
fn spline_arc(
    spline: &BSplineCurve,
    ends: Option<&[Point3; 2]>,
    anchored: bool,
    closed: bool,
    dec: &mut Decider,
) -> Result<DrawingCurve, Refusal> {
    let p = spline.degree();
    let knots: Vec<f64> = spline
        .knots()
        .iter()
        .map(Rational::to_f64_nearest)
        .collect();
    let mut ctrl = Vec::with_capacity(spline.control_count());
    for i in 0..spline.control_count() {
        let x = point(spline.control(i), dec)?;
        let w = out(spline.weight(i), dec)?;
        ctrl.push([w * x[0], w * x[1], w * x[2], w]);
    }
    let spans = decompose(p, &knots, &ctrl);
    if spans.is_empty() {
        return Err(refused("a B-spline carrier has no span"));
    }
    let (lo, hi) = (spans[0].lo, spans[spans.len() - 1].hi);
    let geo = Geo::Spline {
        spans: spans.clone(),
    };
    // Parameter ranges in order; a closed carrier may wrap past its end.
    let mut exact_ends = None;
    let ranges: Vec<[f64; 2]> = match ends {
        None => vec![[lo, hi]],
        Some([a, b]) => {
            let (pa, pb) = (point(a, dec)?, point(b, dec)?);
            exact_ends = Some([pa, pb]);
            let (ta, tb) = (geo.param(pa), geo.param(pb));
            let near =
                |t: f64, k: f64| (t - k).abs() <= 4.0 * f64::EPSILON * (hi - lo).abs().max(1.0);
            if anchored {
                if near(ta, lo) || near(ta, hi) {
                    vec![[lo, hi]]
                } else {
                    vec![[ta, hi], [lo, ta]]
                }
            } else {
                let ta = if closed && near(ta, hi) { lo } else { ta };
                let tb = if closed && near(tb, lo) { hi } else { tb };
                if tb > ta {
                    vec![[ta, tb]]
                } else if closed {
                    vec![[ta, hi], [lo, tb]]
                } else {
                    return Err(refused("a B-spline arc's ends are out of order"));
                }
            }
        }
    };
    let mut pieces: Vec<Span> = Vec::new();
    let mut offset = 0.0;
    for [a, b] in ranges {
        for span in &spans {
            let (from, to) = (a.max(span.lo), b.min(span.hi));
            if to <= from {
                continue;
            }
            let width = span.hi - span.lo;
            pieces.push(Span {
                lo: from + offset,
                hi: to + offset,
                ctrl: restrict(&span.ctrl, (from - span.lo) / width, (to - span.lo) / width),
            });
        }
        offset += hi - lo;
    }
    if pieces.is_empty() {
        return Err(refused("a B-spline arc is empty"));
    }
    // A clamped curve of Bezier pieces: every joint knot `p` times.
    let mut knots = vec![pieces[0].lo; p + 1];
    let mut homogeneous = vec![pieces[0].ctrl[0]];
    for (i, piece) in pieces.iter().enumerate() {
        homogeneous.extend_from_slice(&piece.ctrl[1..]);
        let copies = if i + 1 == pieces.len() { p + 1 } else { p };
        knots.extend(std::iter::repeat_n(piece.hi, copies));
    }
    let mut points: Vec<[f64; 3]> = homogeneous
        .iter()
        .map(|h| [h[0] / h[3], h[1] / h[3], h[2] / h[3]])
        .collect();
    let weights: Vec<f64> = homogeneous.iter().map(|h| h[3]).collect();
    if let Some([a, b]) = exact_ends {
        let last = points.len() - 1;
        points[0] = a;
        points[last] = b;
    }
    let rational = weights.iter().any(|&w| above((w - weights[0]).abs(), 0.0));
    Ok(DrawingCurve::Nurbs {
        degree: p,
        knots,
        points,
        weights: if rational { weights } else { Vec::new() },
    })
}

/// Source: binary128-like working precision, the cap
/// `TracePiece::taylor_interval` uses for patch sections.
const TRACE_BITS: u32 = 128;

/// A trace piece's point and derivative at `t` in binary64: the binary64
/// interval Taylor filter, or the certified ball when the filter cannot.
fn piece_jet(piece: &TracePiece, t: f64, dec: &Decider) -> Result<Jet, Refusal> {
    if let Some(at) = Interval::point(t)
        && let Some(jets) = piece.taylor_interval(at, 1, dec)?
    {
        let (p, d) = (jets[0].map(|x| x.mid()), jets[1].map(|x| x.mid()));
        if p.iter().chain(&d).all(|x| x.is_finite()) {
            return Ok((p, d));
        }
    }
    let t =
        Rational::from_f64_exact(t).ok_or_else(|| refused("a trace parameter is not finite"))?;
    let jets = piece.taylor(&Ball::from_rational(&t, TRACE_BITS), 1, dec)?;
    let mid = |b: &Ball| b.mid().to_f64_nearest();
    Ok((jets[0].each_ref().map(mid), jets[1].each_ref().map(mid)))
}

/// The arc of the traced curve `curve` from `ends[0]` to `ends[1]` along its
/// natural direction, fitted as [`fitted`] fits it: both ends the same point
/// for the whole closed curve from there, `None` for the whole closed curve
/// from the trace's own start. The fitted curve runs through the binary64
/// images of the ends, and a whole closed curve closes on its first point.
/// The answer carries the trace's marks, not the fit's.
pub(in crate::views) fn traced(
    model: &Model,
    curve: CurveId,
    ends: Option<[PointId; 2]>,
    dec: &mut Decider,
) -> Answered<DrawnCurve> {
    let geom = model.geom();
    let traced = trace(geom, curve, ends.map(|e| e[0]), ends.map(|e| e[1]), dec)?;
    let ends = match ends {
        Some([a, b]) => Some([
            point(&geom.position(a), dec)?,
            point(&geom.position(b), dec)?,
        ]),
        None => None,
    };
    let drawn = fitted(traced.value.pieces(), ends, dec)?;
    Ok(Answer::with_marks(drawn, traced.marks))
}

/// Fits each trace piece by cubic Hermite spans in its own parameter and
/// joins the pieces, end to end, into one cubic B-spline: a double knot at
/// every span joint (C1 within a piece) and a triple knot where pieces meet
/// (C0, since the pieces' parameters run at different speeds). `ends`
/// replace the first and last points by the exact ends; without them the
/// trace is a whole closed curve, and the last point becomes the first. A
/// moved end control moves the curve by at most its own move, since the
/// cubic Bernstein weights are at most one, so the deviation grows by it.
fn fitted(
    pieces: &[TracePiece],
    ends: Option<[[f64; 3]; 2]>,
    dec: &mut Decider,
) -> Result<DrawnCurve, Refusal> {
    if pieces.is_empty() {
        return Err(refused("a traced curve has no piece"));
    }
    let mut probe = Decider::new(*dec.policy());
    probe.set_deadline(dec.deadline());
    // Rule T1's floor for the curve's coordinate size, read at the piece
    // ends; the fit runs in a frame at the first point, so that binary64
    // keeps its digits for the curve's own extent.
    let mut size = 0.0_f64;
    let mut origin = None;
    for piece in pieces {
        for t in [0.0, 1.0] {
            let (p, _) = piece_jet(piece, t, &probe)?;
            size = p.iter().fold(size, |m, x| m.max(x.abs()));
            origin.get_or_insert(p);
        }
    }
    let origin = origin.expect("nonempty");
    let floor = tolerance::floor_mm(
        &Rational::from_f64_exact(size.next_up())
            .ok_or_else(|| refused("a traced curve is not finite"))?,
    );
    let target_mm = tolerance::own_fitted(&floor, None)
        .target_mm
        .to_f64_nearest();
    let floor_mm = floor.to_f64_nearest().next_up();
    let mut knots = vec![0.0; 4];
    let mut points: Vec<[f64; 3]> = Vec::new();
    let mut deviation = 0.0_f64;
    let mut spans = 0;
    for (k, piece) in pieces.iter().enumerate() {
        let failed: RefCell<Option<Refusal>> = RefCell::new(None);
        let jet = |t: f64| -> Jet {
            match piece_jet(piece, t, &probe) {
                Ok((p, d)) => (sub(p, origin), d),
                Err(e) => {
                    failed.borrow_mut().get_or_insert(e);
                    ([f64::NAN; 3], [f64::NAN; 3])
                }
            }
        };
        let fit = fit::fit(
            jet,
            (0.0, jet(0.0)),
            (1.0, jet(1.0)),
            target_mm,
            floor_mm,
            gap::EDGE_WORK,
            dec,
        );
        if let Some(e) = failed.into_inner() {
            return Err(e);
        }
        let fit = fit?;
        deviation = deviation.max(fit.deviation_mm);
        spans += fit.breaks.len() - 1;
        let shift = k as f64;
        for &b in &fit.breaks[1..fit.breaks.len() - 1] {
            knots.extend([b + shift; 2]);
        }
        let last = k + 1 == pieces.len();
        knots.extend(std::iter::repeat_n(1.0 + shift, if last { 4 } else { 3 }));
        let skip = usize::from(k > 0);
        points.extend(fit.controls[skip..].iter().map(|&c| add(c, origin)));
    }
    let last = points.len() - 1;
    let [a, b] = ends.unwrap_or([points[0], points[0]]);
    let moved = |from: [f64; 3], to: [f64; 3]| {
        let gap = sub(to, from);
        dot(gap, gap).sqrt()
    };
    deviation += moved(points[0], a).max(moved(points[last], b));
    points[0] = a;
    points[last] = b;
    Ok(DrawnCurve {
        curve: DrawingCurve::Nurbs {
            degree: 3,
            knots,
            points,
            weights: Vec::new(),
        },
        fit: Some(ContourFit {
            deviation_mm: deviation,
            tolerance_mm: (gap::MARGIN * deviation).next_up(),
            target_mm,
            spans,
        }),
    })
}

/// The drawing of `curve` in `model` from `start` to `end` along the
/// curve's natural direction, as a [`super::ViewStroke`] or a
/// [`super::ContourArc`] holds it: both `None` for the whole closed curve,
/// and both the same point for the whole closed curve from there.
///
/// Lines, circles, ellipses, parabolas, hyperbolas and B-splines are exact
/// to binary64. Intersection branches, patch sections and patch pairs are
/// fitted, carry their fit, and the answer carries a [`DRAWING_FIT_MARK`]
/// mark. A whole line, parabola or hyperbola refuses with
/// `app.drawing-curve`.
pub fn drawing_curve(
    model: &Model,
    curve: CurveId,
    start: Option<PointId>,
    end: Option<PointId>,
    dec: &mut Decider,
) -> Answered<DrawnCurve> {
    dec.check_deadline()?;
    let ends = match (start, end) {
        (Some(a), Some(b)) => Some([a, b]),
        (None, None) => None,
        _ => return Err(refused("an arc needs both ends or neither")),
    };
    let form = curve_form(model.geom(), curve, dec)?;
    natural(model, curve, form, ends, dec)
}

/// An edge drawn from its start vertex to its end vertex.
fn edge_drawn(model: &Model, edge: EdgeId, dec: &mut Decider) -> Answered<DrawnCurve> {
    let stored = model.edge(edge);
    let (start, end) = match (stored.start, stored.end) {
        (Some(a), Some(b)) => {
            let (a, b) = (model.vertex(a).point, model.vertex(b).point);
            match stored.sense {
                Sense::Same => (Some(a), Some(b)),
                Sense::Reversed => (Some(b), Some(a)),
            }
        }
        _ => (None, None),
    };
    let drawn = drawing_curve(model, stored.curve, start, end, dec)?;
    Ok(match stored.sense {
        Sense::Same => drawn,
        Sense::Reversed => Answer::with_marks(drawn.value.reversed(), drawn.marks),
    })
}

/// The drawing of one edge of `body`, from its start vertex to its end
/// vertex, as [`drawing_curve`] draws its carrier. An edge outside the body
/// refuses with `app.curve-edge`.
pub fn edge_drawing(
    model: &Model,
    body: BodyId,
    edge: EdgeId,
    dec: &mut Decider,
) -> Result<Answer<DrawnCurve>, QueryRefusal> {
    query(model, &[body], || {
        if !model.edges_of_body(body).contains(&edge) {
            return Err(Refusal::new(
                "app.curve-edge",
                "edge is not in the queried body",
            ));
        }
        edge_drawn(model, edge, dec)
    })
}

/// The boundary loops of one face of `body`, in stored order, each a chain
/// of edge drawings in the direction the loop travels, so that each curve
/// ends where the next starts. Point boundaries, such as a cone apex, draw
/// nothing and are left out. Loops are not sorted outer first. A face
/// outside the body refuses with `app.loop-face`.
pub fn face_drawing_loops(
    model: &Model,
    body: BodyId,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Answer<Vec<Vec<DrawnEdge>>>, QueryRefusal> {
    query(model, &[body], || {
        if !model.faces_of_body(body).contains(&face) {
            return Err(Refusal::new(
                "app.loop-face",
                "face is not in the queried body",
            ));
        }
        let mut marks = MarkLog::new();
        let mut loops = Vec::new();
        for &lp in &model.face(face).loops {
            let Loop::Coedges(coedges) = model.loop_(lp) else {
                continue;
            };
            let mut chain = Vec::with_capacity(coedges.len());
            for &coedge in coedges {
                dec.check_deadline()?;
                let use_ = model.coedge(coedge);
                let drawn = edge_drawn(model, use_.edge, dec)?;
                marks.extend(&drawn.marks);
                chain.push(DrawnEdge {
                    edge: use_.edge,
                    drawn: match use_.sense {
                        Sense::Same => drawn.value,
                        Sense::Reversed => drawn.value.reversed(),
                    },
                });
            }
            loops.push(chain);
        }
        Ok(Answer::with_marks(loops, marks))
    })
}

/// The drawing of one stroke of a hidden-line view, from its start to its
/// end. A fitted contour generator keeps its [`ContourFit`]: its carrier is
/// drawn exactly, and the fit says how far that carrier lies from the true
/// generator.
pub fn stroke_drawing(
    view: &HiddenLineView,
    stroke: &ViewStroke,
    dec: &mut Decider,
) -> Answered<DrawnCurve> {
    let mut drawn = drawing_curve(&view.model, stroke.curve, stroke.start, stroke.end, dec)?;
    if drawn.value.fit.is_none() {
        drawn.value.fit.clone_from(&stroke.fit);
    }
    Ok(drawn)
}

impl ApplicationBody {
    /// The drawing of one edge of this result, from its start to its end.
    pub fn edge_drawing(
        &self,
        edge: EdgeId,
        dec: &mut Decider,
    ) -> Result<Answer<DrawnCurve>, QueryRefusal> {
        edge_drawing(&self.built.model, self.built.result, edge, dec)
    }

    /// The boundary loops of one face of this result as drawing curves.
    pub fn face_drawing_loops(
        &self,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Answer<Vec<Vec<DrawnEdge>>>, QueryRefusal> {
        face_drawing_loops(&self.built.model, self.built.result, face, dec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_reversed_arc_runs_the_same_points_the_other_way() {
        let arc = DrawingCurve::Arc {
            center: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            x_axis: [1.0, 0.0, 0.0],
            radius: 2.0,
            start_deg: 30.0,
            end_deg: 120.0,
        };
        let DrawingCurve::Arc {
            normal,
            start_deg,
            end_deg,
            ..
        } = arc.reversed()
        else {
            panic!("an arc");
        };
        // About -z, angle 240 is angle 120 about +z: the old end.
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12;
        assert!(close(normal[2], -1.0) && close(normal[0], 0.0) && close(normal[1], 0.0));
        assert!(close(start_deg, 240.0) && close(end_deg, 330.0));
    }

    #[test]
    fn the_sweep_stays_in_its_proven_quadrant() {
        use std::f64::consts::{FRAC_PI_2, PI};
        // A sweep just under a half turn whose binary64 estimate rounded
        // past it stays below the half turn.
        let close = |a: f64, b: f64| (a - b).abs() < 1e-12;
        assert!(close(sweep(1, -1, -PI).unwrap(), PI));
        assert!(close(sweep(1, 0, 1.0).unwrap(), FRAC_PI_2));
        // Just short of a whole turn, rounded past zero.
        assert!(close(sweep(-1, 1, 0.1).unwrap(), 2.0 * PI));
        assert!(close(sweep(-1, 1, 6.0).unwrap(), 6.0));
        assert!(sweep(0, 1, 0.0).is_err());
    }
}
