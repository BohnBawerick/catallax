//! Hidden-line views: every edge and contour generator of a solid, cut where
//! visibility can change and marked visible or hidden.
//!
//! The method is OpenCASCADE's `HLRBRep_Algo`, which SMLib's `SmHCR` shares:
//! add the contour generators to the edges, cut every curve where its
//! projection crosses another's, and classify each piece by one interior
//! point. Visibility along a curve changes only where its projection passes
//! a visual contour of the solid (Appel, "The notion of quantitative
//! invisibility and the machine rendering of solids", ACM National
//! Conference 1967): a contour generator, or an edge whose two faces do not
//! both face the viewer or both face away. An edge between two planes that
//! both strictly face the same way is not one, and cuts nothing.
//!
//! Cuts are exact where the pair reduces to a curve against a plane: lines,
//! conics seen edge on, a line against a conic, and conics in parallel planes
//! (`exact.rs`). Other pairs, and every fitted arc, are cut by a binary64
//! search (`cross.rs`) at an exact point of the stroke's own carrier near the
//! crossing; the answer then carries a [`VIEW_SPLIT_MARK`] mark, decision
//! Q14. Each piece is classified by an exact ray toward the viewer from an
//! exact point inside it (`ray.rs`), or for a fitted generator from an exact
//! point of its face beside it (`surface.rs`). Rational patch faces meet the
//! rays through SMLib's bounded binary64 solver, checked on the exact patch,
//! and the answer then carries a [`VIEW_RAY_MARK`] mark, decision Q12's
//! filtered decision for a face no exact intersection serves cheaply.
//!
//! An edge keeps its id along every piece, including where a contour
//! generator runs along it: that generator is not drawn again, and the
//! edge's pieces under it are outline, not smooth.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use catallax_intersect::{
    CurveForm, Cyclic, arc_point, curve_form, curve_is_closed, cyclic_order, order_on_curve,
};
use catallax_numbers::{
    Answer, BoundKind, Decider, Integer, Mark, MarkLog, Outcome, Rational, Real, Refusal, Sign,
    StrictOrder,
};
use catallax_shapes::{
    Circle, Curve, CurveId, Ellipse, Line, Point3, PointId, RVector3, Surface, Vector3,
};
use catallax_topology::numerical::tolerance;
use catallax_topology::{BodyId, BodyKind, EdgeId, FaceId, Model, Sense};

use super::float::{double, rational, unit_direction};
use super::{ContourFit, Contours, contour_generators, same_point};
use crate::ApplicationBody;
use crate::queries::{QueryRefusal, query, sign};

mod cross;
mod exact;
mod geo;
mod ray;
mod surface;

use geo::{Box3, Geo, Sheet, V3};

/// Code of the mark a [`HiddenLineView`] carries when a cut was placed by
/// the binary64 search. Its bound is the largest sheet distance left between
/// the two curves at such a cut: a measured value, not a proof.
pub const VIEW_SPLIT_MARK: &str = "app.view-split";

/// Code of the mark a [`HiddenLineView`] carries when a piece's visibility
/// was decided by a binary64 ray against a rational patch face. Its bound is
/// the tolerance those meetings were checked to: a measured value, not a
/// proof.
pub const VIEW_RAY_MARK: &str = "app.view-ray";

/// What a stroke draws.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StrokeSource {
    /// All or part of an edge of the body, under its own id.
    Edge(EdgeId),
    /// All or part of a contour generator of a curved face that runs along
    /// no edge.
    Contour(FaceId),
}

/// One piece of a hidden-line view: a stretch of one edge or contour
/// generator with one visibility.
#[derive(Clone, Debug)]
pub struct ViewStroke {
    /// The edge or face it comes from.
    pub source: StrokeSource,
    /// Carrier in [`HiddenLineView::model`]: a plain line, circle, ellipse
    /// or B-spline holding the stroke.
    pub curve: CurveId,
    /// Start in the carrier's natural direction; `None` with `end` for the
    /// whole closed carrier.
    pub start: Option<PointId>,
    /// End in the carrier's natural direction.
    pub end: Option<PointId>,
    /// False where the solid hides it from the viewer.
    pub visible: bool,
    /// An edge between faces tangent along it, away from the outline:
    /// drawings show it thin or leave it out.
    pub smooth: bool,
    /// An edge with the same face on both sides.
    pub seam: bool,
    /// The fit of a fitted contour generator, as in
    /// [`super::ContourArc::fit`].
    pub fit: Option<ContourFit>,
}

/// A hidden-line view of one solid for one direction.
#[derive(Clone, Debug)]
pub struct HiddenLineView {
    /// Owns every curve and point handle below; edge and face ids are the
    /// input body's.
    pub model: Model,
    /// Strokes in edge order, then contour order, each run along its
    /// carrier.
    pub strokes: Vec<ViewStroke>,
}

/// Exact shape of a stroke's carrier, as far as the exact cuts read it.
pub(super) enum Form {
    Line {
        origin: Point3,
        dir: Vector3,
    },
    Conic {
        center: Point3,
        normal: Vector3,
        /// For a circle; an ellipse has none.
        radius_sq: Option<Real>,
        /// The conic's plane holds the view direction.
        edge_on: bool,
    },
    Spline,
}

struct Carrier {
    curve: CurveId,
    form: Form,
    /// For a conic, `u` and `w` of `c + cos t u + sin t w`, exactly.
    axes: Option<[Vector3; 2]>,
    geo: Geo,
    closed: bool,
}

pub(super) struct Stroke {
    source: StrokeSource,
    drawn: bool,
    occluder: bool,
    faces: Vec<FaceId>,
    seam: bool,
    fit: Option<ContourFit>,
    pub(super) curve: CurveId,
    pub(super) form: Form,
    axes: Option<[Vector3; 2]>,
    geo: Geo,
    /// Ends in the carrier's natural direction; `None` for a whole closed
    /// carrier.
    pub(super) ends: Option<[PointId; 2]>,
    /// The two ends are one point: the whole closed carrier from there.
    anchored: bool,
    closed: bool,
    /// Binary64 parameters of the span, unwrapped.
    range: [f64; 2],
    /// Sheet bounds of the span.
    sheet: [f64; 4],
}

impl Stroke {
    pub(super) fn is_line(&self) -> bool {
        matches!(self.form, Form::Line { .. })
    }
}

/// The square root of an exact square rational stays rational; anything
/// else becomes a lazy root.
fn root(value: &Real) -> Real {
    fn integer_sqrt(n: &Integer) -> Option<Integer> {
        if n.sign() == Sign::Negative {
            return None;
        }
        if n.is_zero() {
            return Some(Integer::zero());
        }
        let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
        loop {
            let y = (&x + &n.div_rem_floor(&x)?.0)
                .div_rem_floor(&Integer::from_i64(2))?
                .0;
            if y >= x {
                return (&x * &x == *n).then_some(x);
            }
            x = y;
        }
    }
    value
        .rational_ref()
        .and_then(|q| Rational::new(integer_sqrt(q.numerator())?, integer_sqrt(q.denominator())?))
        .map_or_else(|| value.sqrt(), Real::from)
}

fn point64(p: &Point3, floor: &Rational, dec: &mut Decider) -> Result<V3, Refusal> {
    Ok([
        double(&p.x, floor, dec)?,
        double(&p.y, floor, dec)?,
        double(&p.z, floor, dec)?,
    ])
}

fn vector64(v: &Vector3, floor: &Rational, dec: &mut Decider) -> Result<V3, Refusal> {
    Ok([
        double(&v.x, floor, dec)?,
        double(&v.y, floor, dec)?,
        double(&v.z, floor, dec)?,
    ])
}

/// A plain carrier for `curve`: its closed form as a line, circle, ellipse
/// or B-spline. Intersection components are read through their closed form,
/// so later curve-plane meetings take the direct conic paths.
fn carrier(
    scratch: &mut Model,
    curve: CurveId,
    d: &Vector3,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<Carrier, Refusal> {
    match curve_form(scratch.geom(), curve, dec)? {
        CurveForm::Line { origin, dir } => {
            let geo = Geo::Line {
                o: point64(&origin, floor, dec)?,
                v: vector64(&dir, floor, dec)?,
            };
            let plain = scratch.geom_mut().add_curve(Curve::Line(Line {
                origin: origin.clone(),
                dir: dir.clone(),
            }));
            Ok(Carrier {
                curve: plain,
                form: Form::Line { origin, dir },
                axes: None,
                geo,
                closed: false,
            })
        }
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        } => {
            let side = normal.cross(&ref_dir);
            let u = ref_dir.scaled(&root(&(&radius_sq / &ref_dir.norm_squared())));
            let w = side.scaled(&root(&(&radius_sq / &side.norm_squared())));
            let geo = Geo::Conic {
                c: point64(&center, floor, dec)?,
                u: vector64(&u, floor, dec)?,
                w: vector64(&w, floor, dec)?,
            };
            let edge_on = sign(&normal.dot(d), dec)? == 0;
            let plain = scratch.geom_mut().add_curve(Curve::Circle(Circle {
                center: center.clone(),
                normal: normal.clone(),
                radius_sq: radius_sq.clone(),
                ref_dir,
            }));
            Ok(Carrier {
                curve: plain,
                form: Form::Conic {
                    center,
                    normal,
                    radius_sq: Some(radius_sq),
                    edge_on,
                },
                axes: Some([u, w]),
                geo,
                closed: true,
            })
        }
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        } => {
            let side = normal.cross(&major_dir);
            let u = major_dir.scaled(&root(&(&major_sq / &major_dir.norm_squared())));
            let w = side.scaled(&root(&(&minor_sq / &side.norm_squared())));
            let geo = Geo::Conic {
                c: point64(&center, floor, dec)?,
                u: vector64(&u, floor, dec)?,
                w: vector64(&w, floor, dec)?,
            };
            let edge_on = sign(&normal.dot(d), dec)? == 0;
            let plain = scratch.geom_mut().add_curve(Curve::Ellipse(Ellipse {
                center: center.clone(),
                normal: normal.clone(),
                major_dir,
                major_sq,
                minor_sq,
            }));
            Ok(Carrier {
                curve: plain,
                form: Form::Conic {
                    center,
                    normal,
                    radius_sq: None,
                    edge_on,
                },
                axes: Some([u, w]),
                geo,
                closed: true,
            })
        }
        CurveForm::BSpline(spline) => {
            let plain = if matches!(scratch.geom().curve(curve), Curve::BSpline(_)) {
                curve
            } else {
                scratch.geom_mut().add_curve(Curve::BSpline(spline.clone()))
            };
            let knots: Vec<f64> = spline
                .knots()
                .iter()
                .map(Rational::to_f64_nearest)
                .collect();
            let mut ctrl = Vec::with_capacity(spline.control_count());
            for i in 0..spline.control_count() {
                let p = point64(spline.control(i), floor, dec)?;
                let w = double(spline.weight(i), floor, dec)?;
                ctrl.push([w * p[0], w * p[1], w * p[2], w]);
            }
            let geo = Geo::Spline {
                spans: geo::decompose(spline.degree(), &knots, &ctrl),
            };
            let closed = curve_is_closed(scratch.geom(), plain, dec)?;
            Ok(Carrier {
                curve: plain,
                form: Form::Spline,
                axes: None,
                geo,
                closed,
            })
        }
        _ => Err(Refusal::new(
            "app.view-curve",
            "hidden-line views of edges on intersection branches, patch sections and traced \
             curves need their drawing conversion; unfinished",
        )),
    }
}

/// Builds a stroke on `curve` between `ends`, or `None` when it projects to
/// a point: a line along the view.
#[allow(clippy::too_many_arguments)]
fn stroke(
    scratch: &mut Model,
    source: StrokeSource,
    curve: CurveId,
    ends: Option<[PointId; 2]>,
    sheet: &Sheet,
    d: &Vector3,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<(Stroke, Box3, bool), Refusal> {
    let c = carrier(scratch, curve, d, floor, dec)?;
    let on_sheet = match &c.form {
        Form::Line { dir, .. } => sign(&dir.cross(d).norm_squared(), dec)? != 0,
        _ => true,
    };
    let anchored = match ends {
        Some([a, b]) => c.closed && same_point(scratch, a, b, dec)?,
        None => false,
    };
    let range = match ends {
        None => c.geo.whole(),
        Some(_) if anchored && matches!(c.geo, Geo::Spline { .. }) => c.geo.whole(),
        Some([a, b]) => {
            let a = point64(&scratch.geom().position(a), floor, dec)?;
            let start = c.geo.param(a);
            if anchored {
                [start, start + TAU]
            } else {
                let mut end = c
                    .geo
                    .param(point64(&scratch.geom().position(b), floor, dec)?);
                if matches!(c.geo, Geo::Conic { .. }) {
                    while end <= start {
                        end += TAU;
                    }
                }
                [start, end]
            }
        }
    };
    let bounds = sheet.bounds(&c.geo, range[0], range[1]);
    let box3 = geo::box3(&c.geo, range[0], range[1]);
    Ok((
        Stroke {
            source,
            drawn: on_sheet,
            occluder: on_sheet,
            faces: Vec::new(),
            seam: false,
            fit: None,
            curve: c.curve,
            form: c.form,
            axes: c.axes,
            geo: c.geo,
            ends,
            anchored,
            closed: c.closed,
            range,
            sheet: bounds,
        },
        box3,
        on_sheet,
    ))
}

/// The outward normal's facing of a planar face: the sign of `n . d`, or
/// `None` for a curved face.
fn plane_facing(
    model: &Model,
    face: FaceId,
    d: &Vector3,
    dec: &mut Decider,
) -> Result<Option<i8>, Refusal> {
    let stored = model.face(face);
    let Surface::Plane(plane) = model.geom().surface(stored.surface) else {
        return Ok(None);
    };
    let s = sign(&plane.normal.dot(d), dec)?;
    Ok(Some(if stored.sense == Sense::Same { s } else { -s }))
}

/// A conservative binary64 box of a face, from its boundary and its
/// surface: a plane, cylinder or cone face lies on segments between its
/// boundary points (and a cone's apex); spheres and tori use their whole
/// box. `None` where no cheap bound exists.
fn face_box(
    model: &Model,
    face: FaceId,
    edge_boxes: &std::collections::BTreeMap<EdgeId, Box3>,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<Option<Box3>, Refusal> {
    fn grow(b: Box3, bounds: &mut Option<Box3>) {
        *bounds = Some(bounds.map_or(b, |a| geo::union(a, b)));
    }
    let stored = model.face(face);
    let mut bounds: Option<Box3> = None;
    let whole = |c: V3, r: f64| -> Box3 {
        [
            [c[0] - r, c[1] - r, c[2] - r],
            [c[0] + r, c[1] + r, c[2] + r],
        ]
    };
    match model.geom().surface(stored.surface) {
        Surface::Plane(_) | Surface::Cylinder(_) | Surface::EllipticCylinder(_) => {}
        Surface::Cone(cone) => grow(
            geo::point_box(point64(&cone.apex, floor, dec)?),
            &mut bounds,
        ),
        Surface::Sphere(s) => {
            let c = point64(&s.center, floor, dec)?;
            return Ok(Some(whole(c, double(&s.radius, floor, dec)?.abs())));
        }
        Surface::Torus(t) => {
            let c = point64(&t.center, floor, dec)?;
            let r = double(&t.major, floor, dec)?.abs() + double(&t.minor, floor, dec)?.abs();
            return Ok(Some(whole(c, r)));
        }
        _ => return Ok(None),
    }
    for &lp in &stored.loops {
        match model.loop_(lp) {
            catallax_topology::Loop::Vertex(v) => {
                let p = point64(&model.geom().position(model.vertex(*v).point), floor, dec)?;
                grow(geo::point_box(p), &mut bounds);
            }
            catallax_topology::Loop::Coedges(coedges) => {
                for &coedge in coedges {
                    let Some(&b) = edge_boxes.get(&model.coedge(coedge).edge) else {
                        return Ok(None);
                    };
                    grow(b, &mut bounds);
                }
            }
        }
    }
    Ok(bounds)
}

/// An exact point of the stroke's carrier at binary64 parameter `t`.
fn exact_point(scratch: &Model, s: &Stroke, t: f64, dec: &mut Decider) -> Result<Point3, Refusal> {
    match &s.form {
        Form::Line { origin, dir } => Ok(origin.plus(&dir.scaled(&Real::from(rational(t)?)))),
        Form::Conic { center, .. } => {
            let [u, w] = s.axes.as_ref().expect("a conic stroke keeps its axes");
            // The half-angle form keeps the point exactly on the conic:
            // cos t = (1 - s^2) / (1 + s^2), sin t = 2 s / (1 + s^2).
            let mut angle = t.rem_euclid(TAU);
            let mut flip = Real::one();
            if angle > PI {
                angle -= TAU;
            }
            if angle.abs() > FRAC_PI_2 {
                angle -= PI.copysign(angle);
                flip = -Real::one();
            }
            let s = Real::from(rational((angle / 2.0).tan())?);
            let one = Real::one();
            let den = &one + &s.square();
            let cos = &(&(&one - &s.square()) / &den) * &flip;
            let sin = &(&(&Real::from_i64(2) * &s) / &den) * &flip;
            Ok(center.plus(&u.scaled(&cos).plus(&w.scaled(&sin))))
        }
        Form::Spline => {
            let Curve::BSpline(curve) = scratch.geom().curve(s.curve) else {
                unreachable!("a spline stroke's carrier is a B-spline");
            };
            // Knots rounded to binary64 can reach past the exact domain.
            let [lo, hi] = curve.domain();
            let at = rational(t)?.max(lo.clone()).min(hi.clone());
            curve.point_at(&Real::from(at), dec)
        }
    }
}

/// The parameter margin at `t`: rule T1's floor over the carrier's speed.
fn slack(s: &Stroke, t: f64, margin: f64) -> f64 {
    let (_, dp) = s.geo.jet(t);
    let speed = super::float::dot(dp, dp).sqrt();
    if speed > 0.0 {
        margin / speed
    } else {
        f64::INFINITY
    }
}

/// Whether `a` comes before `b` along the stroke: on a closed carrier,
/// going round from `anchor`.
fn ordered(
    scratch: &Model,
    s: &Stroke,
    anchor: Option<PointId>,
    a: PointId,
    b: PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let outcome = if s.closed {
        let Some(anchor) = anchor else {
            return Ok(true);
        };
        if a == anchor {
            return Ok(true);
        }
        match cyclic_order(scratch.geom(), s.curve, [anchor, a, b], dec) {
            Outcome::Proven(order) => Outcome::Proven(order == Cyclic::Forward),
            Outcome::ProvenCoincident(c) => Outcome::ProvenCoincident(c),
            Outcome::NotExact(n) => Outcome::NotExact(n),
            Outcome::Refused(r) => Outcome::Refused(r),
        }
    } else {
        match order_on_curve(scratch.geom(), s.curve, a, b, dec) {
            Outcome::Proven(order) => Outcome::Proven(order == StrictOrder::Less),
            Outcome::ProvenCoincident(c) => Outcome::ProvenCoincident(c),
            Outcome::NotExact(n) => Outcome::NotExact(n),
            Outcome::Refused(r) => Outcome::Refused(r),
        }
    };
    match outcome {
        Outcome::Proven(before) => Ok(before),
        Outcome::ProvenCoincident(_) => Ok(true),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(r) => Err(r),
    }
}

/// A cut candidate: its point, its binary64 parameter, and the residual of
/// a binary64 cut.
type Cut = (PointId, f64, Option<f64>);

/// The cuts strictly inside a stroke's span, distinct and in order.
///
/// Binary64 parameters decide wherever they are farther apart than twice
/// the parameter margin (rule T1's floor over the carrier's speed): every
/// binary64 image here is within the floor of its exact point. Closer than
/// that, exact predicates decide: membership near the span's ends, identity
/// and order of near neighbours.
fn settle(
    scratch: &mut Model,
    s: &Stroke,
    candidates: Vec<(PointId, Option<f64>)>,
    margin: f64,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<Vec<Cut>, Refusal> {
    let [lo, hi] = s.range;
    let periodic = s.closed && !matches!(s.geo, Geo::Spline { .. });
    // Spline strokes are cut by the binary64 search only, and decide by
    // their binary64 parameters.
    let exact = !matches!(s.form, Form::Spline);
    // On an anchored loop both ends are one point; a point near it is
    // placed by its order against the loop's middle.
    let mut middle: Option<PointId> = None;
    let mut kept: Vec<Cut> = Vec::new();
    for (point, residual) in candidates {
        dec.check_deadline()?;
        let p64 = point64(&scratch.geom().position(point), floor, dec)?;
        let mut t = s.geo.param(p64);
        if periodic {
            t = lo + (t - lo).rem_euclid(TAU);
        }
        let Some([a, b]) = s.ends else {
            kept.push((point, t, residual));
            continue;
        };
        let mu = slack(s, t, margin);
        let near_start = t - lo <= 2.0 * mu || (periodic && lo + TAU - t <= 2.0 * mu);
        let near_end = (hi - t).abs() <= 2.0 * mu || (periodic && t > hi);
        if !near_start && !near_end {
            if t > lo && t < hi {
                kept.push((point, t, residual));
            }
            continue;
        }
        if !exact {
            // A spline stroke's ends need not lie on its fitted carrier;
            // a cut this near one is that end.
            continue;
        }
        if same_point(scratch, point, a, dec)? || same_point(scratch, point, b, dec)? {
            continue;
        }
        let inside = if s.anchored {
            true
        } else if s.closed {
            match cyclic_order(scratch.geom(), s.curve, [a, point, b], dec) {
                Outcome::Proven(order) => order == Cyclic::Forward,
                Outcome::ProvenCoincident(_) => false,
                Outcome::NotExact(n) => return Err(Refusal::from(n)),
                Outcome::Refused(r) => return Err(r),
            }
        } else {
            ordered(scratch, s, None, a, point, dec)? && ordered(scratch, s, None, point, b, dec)?
        };
        if !inside {
            continue;
        }
        let t = if s.anchored {
            let m = match middle {
                Some(m) => m,
                None => {
                    let p = exact_point(scratch, s, lo + PI, dec)?;
                    let m = scratch
                        .geom_mut()
                        .add_point(catallax_shapes::PointDef::Explicit(p));
                    middle = Some(m);
                    m
                }
            };
            if ordered(scratch, s, Some(a), point, m, dec)? {
                lo
            } else {
                hi
            }
        } else if periodic && t > hi {
            if t - hi < lo + TAU - t { hi } else { lo }
        } else {
            t.clamp(lo, hi)
        };
        kept.push((point, t, residual));
    }
    kept.sort_by(|x, y| x.1.total_cmp(&y.1));
    let mut out: Vec<Cut> = Vec::with_capacity(kept.len());
    'next: for cut in kept {
        let mu = slack(s, cut.1, margin);
        for previous in out.iter_mut().rev() {
            if cut.1 - previous.1 > 2.0 * mu.max(slack(s, previous.1, margin)) {
                break;
            }
            if !exact || same_point(scratch, cut.0, previous.0, dec)? {
                if previous.2.is_some() && cut.2.is_none() {
                    *previous = cut;
                }
                continue 'next;
            }
        }
        let anchor = s.ends.map(|e| e[0]).or(out.first().map(|c| c.0));
        let mut at = out.len();
        while exact
            && at > 0
            && cut.1 - out[at - 1].1 <= 2.0 * mu.max(slack(s, out[at - 1].1, margin))
            && !ordered(scratch, s, anchor, out[at - 1].0, cut.0, dec)?
        {
            at -= 1;
        }
        out.insert(at, cut);
    }
    Ok(out)
}

/// The pieces of a stroke between its ends and cuts, with their binary64
/// parameter intervals.
fn pieces(s: &Stroke, cuts: &[Cut]) -> Vec<(Option<PointId>, Option<PointId>, f64, f64)> {
    let [lo, hi] = s.range;
    match s.ends {
        Some([a, b]) if !s.anchored => {
            let mut points = vec![(a, lo)];
            points.extend(cuts.iter().map(|c| (c.0, c.1)));
            points.push((b, hi));
            points
                .windows(2)
                .map(|w| (Some(w[0].0), Some(w[1].0), w[0].1, w[1].1))
                .collect()
        }
        Some([a, _]) => {
            let mut points = vec![(a, lo)];
            points.extend(cuts.iter().map(|c| (c.0, c.1)));
            let n = points.len();
            (0..n)
                .map(|i| {
                    let (p, t0) = points[i];
                    let (q, t1) = if i + 1 < n { points[i + 1] } else { (a, hi) };
                    (Some(p), Some(q), t0, t1)
                })
                .collect()
        }
        None if cuts.is_empty() => vec![(None, None, lo, hi)],
        None => {
            let n = cuts.len();
            let period = hi - lo;
            (0..n)
                .map(|i| {
                    let (p, t0) = (cuts[i].0, cuts[i].1);
                    let (q, t1) = if i + 1 < n {
                        (cuts[i + 1].0, cuts[i + 1].1)
                    } else {
                        (cuts[0].0, cuts[0].1 + period)
                    };
                    (Some(p), Some(q), t0, t1)
                })
                .collect()
        }
    }
}

/// An exact point strictly inside a piece: the carrier's point at the
/// middle binary64 parameter, when the margin or exact order proves it
/// inside; otherwise the intersection crate's arc point.
fn sample(
    scratch: &mut Model,
    s: &Stroke,
    piece: (Option<PointId>, Option<PointId>, f64, f64),
    margin: f64,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let (from, to, t0, t1) = piece;
    let mid = (t0 + t1) / 2.0;
    let mu = slack(s, mid, margin)
        .max(slack(s, t0, margin))
        .max(slack(s, t1, margin));
    let whole = from.is_none() || from == to;
    let spline = matches!(s.form, Form::Spline);
    if whole || spline || (mid - t0 > 2.0 * mu && t1 - mid > 2.0 * mu) {
        return exact_point(scratch, s, mid, dec);
    }
    // Too short for binary64 to place the middle: the carrier's point there
    // still keeps short numbers when exact order confirms it.
    let candidate = exact_point(scratch, s, mid, dec)?;
    let (Some(a), Some(b)) = (from, to) else {
        unreachable!("a whole piece returned above");
    };
    let c = scratch
        .geom_mut()
        .add_point(catallax_shapes::PointDef::Explicit(candidate.clone()));
    let inside = if s.closed {
        match cyclic_order(scratch.geom(), s.curve, [a, c, b], dec) {
            Outcome::Proven(order) => order == Cyclic::Forward,
            Outcome::ProvenCoincident(_) => false,
            Outcome::NotExact(n) => return Err(Refusal::from(n)),
            Outcome::Refused(r) => return Err(r),
        }
    } else {
        !same_point(scratch, a, c, dec)?
            && !same_point(scratch, c, b, dec)?
            && ordered(scratch, s, None, a, c, dec)?
            && ordered(scratch, s, None, c, b, dec)?
    };
    if inside {
        return Ok(candidate);
    }
    let point = arc_point(scratch.geom(), s.curve, from, to, dec)?;
    Ok(point.value)
}

/// The visibility and smoothness of one piece at its sample point.
#[allow(clippy::too_many_arguments)]
fn classify(
    scratch: &mut Model,
    s: &Stroke,
    p: &Point3,
    scene: &ray::Scene,
    d: &Vector3,
    d64: V3,
    floor: &Rational,
    marks: &mut MarkLog,
    numeric: &mut bool,
    dec: &mut Decider,
) -> Result<(bool, bool), Refusal> {
    // A fitted generator's ray, and an exact one's on a spline face, starts
    // from the face beside it.
    let start = match s.source {
        StrokeSource::Contour(face)
            if s.fit.is_some()
                || matches!(
                    scratch.geom().surface(scratch.face(face).surface),
                    Surface::BSpline(_)
                ) =>
        {
            let p64 = point64(p, floor, dec)?;
            let offset = s.fit.as_ref().map_or(scene.margin, |f| f.tolerance_mm);
            surface::beside(scratch, face, p64, d64, d, offset, scene.size, floor, dec)?
        }
        _ => p.clone(),
    };
    let p64 = point64(&start, floor, dec)?;
    let hidden = ray::hidden(scratch, scene, &start, p64, d, d64, marks, numeric, dec)?;
    let smooth = match (s.source, s.faces.as_slice()) {
        (StrokeSource::Edge(_), &[f, g]) if f != g => {
            let surfaces = [scratch.face(f).surface, scratch.face(g).surface];
            let n1 = scratch.geom().surface(surfaces[0]).normal_at(p, dec)?;
            let n2 = scratch.geom().surface(surfaces[1]).normal_at(p, dec)?;
            sign(&n1.cross(&n2).norm_squared(), dec)? == 0 && sign(&n1.dot(d), dec)? != 0
        }
        _ => false,
    };
    Ok((!hidden, smooth))
}

/// Hidden-line view of `body` seen along `direction`, which points toward
/// the viewer.
///
/// Every edge and every contour generator ([`contour_generators`]) comes
/// out in pieces, each wholly visible or wholly hidden. Pieces of edges and
/// exact generators lie on exact carriers; where both curves of a crossing
/// are lines or conics that reduce to a plane, the cut is exact too, and
/// elsewhere it is a binary64 cut on the exact carrier and the answer
/// carries a [`VIEW_SPLIT_MARK`] mark. Fitted generators keep their
/// [`ContourFit`] and the contour answer's mark. A piece is hidden when the
/// open ray from it toward the viewer meets the solid's inside; a ray that
/// only grazes the boundary, as along a face seen edge on, leaves it
/// visible. Rays meet rational patch faces through SMLib's bounded binary64
/// solver, checked on the exact patch, and the answer then carries a
/// [`VIEW_RAY_MARK`] mark. Lines along the view draw nothing and are left
/// out.
///
/// Edges keep their ids. A contour generator that runs along an edge is
/// drawn as that edge, whose pieces there are not smooth. Seams are marked,
/// not dropped. Edges on intersection branches without a closed form, patch
/// sections and traced curves refuse with `app.view-curve`, as do the
/// faces [`contour_generators`] refuses; a body that is not one solid
/// refuses with `app.view-solid`.
pub fn hidden_line_view(
    model: &Model,
    body: BodyId,
    direction: &RVector3,
    dec: &mut Decider,
) -> Result<Answer<HiddenLineView>, QueryRefusal> {
    query(model, &[body], || {
        dec.check_deadline()?;
        if direction.is_zero() {
            return Err(Refusal::new("app.view-direction", "view direction is zero"));
        }
        if model.body(body).kind != BodyKind::Solid {
            return Err(Refusal::new(
                "app.view-solid",
                "a hidden-line view needs one solid body",
            ));
        }
        let contours = contour_generators(model, body, direction, dec).map_err(|q| q.refusal)?;
        let mut marks = contours.marks;
        let Contours {
            model: mut scratch,
            arcs,
        } = contours.value;
        let d = direction.to_real();
        let fine = tolerance::floor_minimum_mm();
        let d64 = unit_direction(&d, &fine, dec)?;
        let sheet = Sheet::new(d64);

        let mut strokes: Vec<Stroke> = Vec::new();
        let mut edge_boxes = std::collections::BTreeMap::new();
        for edge in model.edges_of_body(body) {
            dec.check_deadline()?;
            let stored = scratch.edge(edge).clone();
            let mut faces: Vec<FaceId> = Vec::new();
            let uses = scratch.coedges_of_edge(edge);
            for &coedge in &uses {
                if let Some((_, face)) = scratch.owner_of_coedge(coedge)
                    && !faces.contains(&face)
                {
                    faces.push(face);
                }
            }
            let ends = match (stored.start, stored.end) {
                (Some(a), Some(b)) => {
                    let (a, b) = (scratch.vertex(a).point, scratch.vertex(b).point);
                    Some(match stored.sense {
                        Sense::Same => [a, b],
                        Sense::Reversed => [b, a],
                    })
                }
                _ => None,
            };
            let (mut s, box3, on_sheet) = stroke(
                &mut scratch,
                StrokeSource::Edge(edge),
                stored.curve,
                ends,
                &sheet,
                &d,
                &fine,
                dec,
            )?;
            edge_boxes.insert(edge, box3);
            if !on_sheet {
                continue;
            }
            s.seam = faces.len() == 1 && uses.len() >= 2;
            if let [f, g] = faces.as_slice() {
                let (a, b) = (
                    plane_facing(&scratch, *f, &d, dec)?,
                    plane_facing(&scratch, *g, &d, dec)?,
                );
                if let (Some(a), Some(b)) = (a, b)
                    && a == b
                    && a != 0
                {
                    s.occluder = false;
                }
            }
            s.faces = faces;
            strokes.push(s);
        }
        for arc in &arcs {
            dec.check_deadline()?;
            let ends = match (arc.start, arc.end) {
                (Some(a), Some(b)) => Some([a, b]),
                _ => None,
            };
            let (mut s, _, on_sheet) = stroke(
                &mut scratch,
                StrokeSource::Contour(arc.face),
                arc.curve,
                ends,
                &sheet,
                &d,
                &fine,
                dec,
            )?;
            if !on_sheet {
                continue;
            }
            // A generator along an edge is drawn as that edge; it still
            // bounds what it hides.
            s.drawn = arc.on_faces.len() == 1;
            s.faces = arc.on_faces.clone();
            s.fit = arc.fit.clone();
            strokes.push(s);
        }

        // Rule T1's floor for the body's coordinate size bounds every
        // binary64 image above.
        let mut size = 0.0_f64;
        for b in edge_boxes.values() {
            for corner in b {
                size = corner.iter().fold(size, |m, x| m.max(x.abs()));
            }
        }
        let mut faces = Vec::new();
        for face in model.faces_of_body(body) {
            let b = face_box(&scratch, face, &edge_boxes, &fine, dec)?;
            if let Some(b) = b {
                for corner in &b {
                    size = corner.iter().fold(size, |m, x| m.max(x.abs()));
                }
            }
            faces.push((face, b));
        }
        let floor = tolerance::floor_mm(&rational(size.next_up())?);
        let margin = floor.to_f64_nearest().next_up();
        let mut vertices = Vec::new();
        for vertex in model.vertices_of_body(body) {
            let point = scratch.vertex(vertex).point;
            vertices.push((point, point64(&scratch.geom().position(point), &fine, dec)?));
        }
        let mut patches = std::collections::BTreeMap::new();
        for face in faces.iter().map(|f| f.0).collect::<Vec<_>>() {
            let cells = match scratch.geom().surface(scratch.face(face).surface) {
                Surface::BSpline(catallax_shapes::SplineSupport::RationalPatch(_)) => vec![face],
                Surface::BSpline(catallax_shapes::SplineSupport::RationalAtlas(_)) => scratch
                    .face_chart_cover(face)
                    .map_or_else(|| vec![face], <[FaceId]>::to_vec),
                _ => continue,
            };
            let mut charts = Vec::with_capacity(cells.len());
            for cell in cells {
                let chart =
                    catallax_topology::RationalPatchFace::of_chart_face(&scratch, cell, dec)?
                        .ok_or_else(|| {
                            Refusal::new(
                                "app.view-trim",
                                "a hidden-line view of a rational spline face needs complete patch \
                             charts",
                            )
                        })?;
                let support = chart.support().clone();
                // A positive-weight patch lies in its control hull (Piegl
                // and Tiller, section 4.4), so its control box bounds it.
                let mut bounds: Option<Box3> = None;
                for c in support.carrier().controls() {
                    let b = geo::point_box(point64(c, &fine, dec)?);
                    bounds = Some(bounds.map_or(b, |a| geo::union(a, b)));
                }
                charts.push((support, bounds.expect("a patch has controls")));
            }
            patches.insert(face, charts);
        }
        // A spline face is bounded by its charts' boxes.
        for (face, bounds) in &mut faces {
            if let Some(charts) = patches.get(face) {
                *bounds = charts.iter().map(|c| c.1).reduce(geo::union);
            }
        }
        let scene = ray::Scene {
            faces,
            patches,
            vertices,
            margin,
            size,
        };

        // Cuts: every drawn stroke against every occluder whose sheet
        // bounds meet its own.
        let mut work = 0_u64;
        let mut candidates: Vec<Vec<(PointId, Option<f64>)>> = vec![Vec::new(); strokes.len()];
        for i in 0..strokes.len() {
            if !strokes[i].drawn {
                continue;
            }
            dec.check_deadline()?;
            let a = &strokes[i];
            for fold in exact::folds(&mut scratch, a, &d, &mut marks, dec)? {
                candidates[i].push((fold, None));
            }
            if matches!(a.form, Form::Spline) {
                for t in cross::slow_points(
                    &sheet,
                    &a.geo,
                    a.range,
                    a.ends.is_none() || a.anchored,
                    &mut work,
                ) {
                    let p = exact_point(&scratch, a, t, dec)?;
                    let p = scratch
                        .geom_mut()
                        .add_point(catallax_shapes::PointDef::Explicit(p));
                    candidates[i].push((p, Some(0.0)));
                }
                for c in cross::crossings(
                    &sheet, &a.geo, a.range, &a.geo, a.range, true, margin, &mut work,
                ) {
                    for t in [c.t, c.u] {
                        let p = exact_point(&scratch, a, t, dec)?;
                        let p = scratch
                            .geom_mut()
                            .add_point(catallax_shapes::PointDef::Explicit(p));
                        candidates[i].push((p, Some(c.residual)));
                    }
                }
            }
            for j in 0..strokes.len() {
                let (a, b) = (&strokes[i], &strokes[j]);
                if i == j || !b.occluder {
                    continue;
                }
                let (x, y) = (a.sheet, b.sheet);
                if x[0] > y[2] + margin
                    || y[0] > x[2] + margin
                    || x[1] > y[3] + margin
                    || y[1] > x[3] + margin
                {
                    continue;
                }
                dec.check_deadline()?;
                match exact::cuts(&mut scratch, a, b, &d, &mut marks, dec)? {
                    Some(points) => {
                        for p in points {
                            let q =
                                sheet.project(point64(&scratch.geom().position(p), &fine, dec)?);
                            if q[0] >= y[0] - margin
                                && q[0] <= y[2] + margin
                                && q[1] >= y[1] - margin
                                && q[1] <= y[3] + margin
                            {
                                candidates[i].push((p, None));
                            }
                        }
                    }
                    None => {
                        // A fitted arc lies within its tolerance of the
                        // generator it stands for, whose ends lie on edges.
                        let reach = margin
                            + a.fit.as_ref().map_or(0.0, |f| f.tolerance_mm)
                            + b.fit.as_ref().map_or(0.0, |f| f.tolerance_mm);
                        let found = cross::crossings(
                            &sheet, &a.geo, a.range, &b.geo, b.range, false, reach, &mut work,
                        );
                        // A fitted generator's exact ends that lie on A stand
                        // for the binary64 cuts near them.
                        let mut ends = Vec::new();
                        if b.fit.is_some() && !matches!(a.form, Form::Spline) {
                            for e in b.ends.into_iter().flatten() {
                                if super::lies_on(&scratch, a.curve, e, dec)? {
                                    ends.push((
                                        e,
                                        point64(&scratch.geom().position(e), &fine, dec)?,
                                    ));
                                }
                            }
                        }
                        for c in found {
                            let at = a.geo.point(c.t);
                            if let Some(&(e, _)) = ends
                                .iter()
                                .find(|(_, e64)| geo::distance(at, *e64) <= reach)
                            {
                                candidates[i].push((e, None));
                                continue;
                            }
                            let p = exact_point(&scratch, a, c.t, dec)?;
                            let p = scratch
                                .geom_mut()
                                .add_point(catallax_shapes::PointDef::Explicit(p));
                            candidates[i].push((p, Some(c.residual)));
                        }
                    }
                }
            }
            dec.charge_work(std::mem::take(&mut work))?;
        }

        // Pieces, classified and merged.
        let mut numeric = false;
        let mut out = Vec::new();
        let mut worst: Option<f64> = None;
        for (i, s) in strokes.iter().enumerate() {
            if !s.drawn {
                continue;
            }
            let cuts = settle(
                &mut scratch,
                s,
                std::mem::take(&mut candidates[i]),
                margin,
                &fine,
                dec,
            )?;
            for cut in &cuts {
                if let Some(r) = cut.2 {
                    worst = Some(worst.map_or(r, |w: f64| w.max(r)));
                }
            }
            let parts = pieces(s, &cuts);
            let mut labels = Vec::with_capacity(parts.len());
            for &part in &parts {
                dec.check_deadline()?;
                let p = sample(&mut scratch, s, part, margin, dec)?;
                labels.push(classify(
                    &mut scratch,
                    s,
                    &p,
                    &scene,
                    &d,
                    d64,
                    &fine,
                    &mut marks,
                    &mut numeric,
                    dec,
                )?);
            }
            let n = parts.len();
            let emit =
                |from: Option<PointId>, to: Option<PointId>, label: (bool, bool)| ViewStroke {
                    source: s.source,
                    curve: s.curve,
                    start: from,
                    end: to,
                    visible: label.0,
                    smooth: label.1,
                    seam: s.seam,
                    fit: s.fit.clone(),
                };
            let closed_loop = s.ends.is_none() || s.anchored;
            if closed_loop && labels.iter().all(|l| *l == labels[0]) {
                let (from, to) = if s.anchored {
                    (parts[0].0, parts[0].0)
                } else {
                    (None, None)
                };
                out.push(emit(from, to, labels[0]));
                continue;
            }
            // Runs of equal labels; a closed loop starts where one run ends.
            let first = if closed_loop {
                (0..n)
                    .find(|&k| labels[k] != labels[(k + n - 1) % n])
                    .unwrap_or(0)
            } else {
                0
            };
            let mut run: Option<(usize, usize)> = None;
            for step in 0..n {
                let k = (first + step) % n;
                if let Some((start, end)) = run
                    && labels[k] != labels[start]
                {
                    out.push(emit(parts[start].0, parts[end].1, labels[start]));
                    run = None;
                }
                run = Some(run.map_or((k, k), |(start, _)| (start, k)));
            }
            if let Some((start, end)) = run {
                out.push(emit(parts[start].0, parts[end].1, labels[start]));
            }
        }
        if numeric {
            marks.push(Mark::new(
                VIEW_RAY_MARK,
                BoundKind::Distance,
                rational(margin)?,
                "a visibility decided by a binary64 ray against a rational patch face: the \
                 tolerance its meetings were checked to, not a proof",
            ));
        }
        if let Some(worst) = worst {
            marks.push(Mark::new(
                VIEW_SPLIT_MARK,
                BoundKind::Distance,
                rational(worst.next_up())?,
                "a visibility change placed by a binary64 crossing search: the sheet distance \
                 left between the two curves there, not a proof",
            ));
        }
        Ok(Answer::with_marks(
            HiddenLineView {
                model: scratch,
                strokes: out,
            },
            marks,
        ))
    })
}

impl ApplicationBody {
    /// Hidden-line view of this result for a direction toward the viewer.
    pub fn hidden_line_view(
        &self,
        direction: &RVector3,
        dec: &mut Decider,
    ) -> Result<Answer<HiddenLineView>, QueryRefusal> {
        hidden_line_view(&self.built.model, self.built.result, direction, dec)
    }
}
