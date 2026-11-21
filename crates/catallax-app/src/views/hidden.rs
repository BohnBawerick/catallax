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
//!
//! Edges on intersection branches, patch sections and patch pairs have no
//! closed form. Each is traced and fitted over its span as its drawing is
//! (`drawing.rs`), as OpenCASCADE's HLR works on an intersection edge's
//! approximated 3D curve, and the stroke keeps the fit. Such a piece is
//! hidden by its own faces when both face away beside it (SMLib's
//! `SmHCR::IsBackFacingEdge`); otherwise its ray starts from an exact point
//! of a facing face just beside it (`surface.rs`), and the answer carries a
//! [`super::DRAWING_FIT_MARK`] mark. That point lies no farther from the
//! edge than the piece's own ends. A piece whose ends lie nearer its sample
//! than twice the fit's tolerance has no such point: it takes the visibility
//! of the piece before it, and the answer carries a
//! [`VIEW_SHORT_PIECE_MARK`] mark. OpenCASCADE's `HLRBRep_Hider::Hide` leaves
//! visible a hidden interval no longer than 1e-7 in parameter, and SMLib's
//! `SmHCR` places a projected vertex within its zone tolerance of an interval
//! end on that end (`FindInterval`); the tolerance here is the fit's own.

use std::f64::consts::{PI, TAU};

use catallax_intersect::{
    CurveForm, Cyclic, arc_point, curve_form, curve_is_closed, cyclic_order, order_on_curve,
};
use catallax_numbers::{
    Answer, BoundKind, Decider, Mark, MarkLog, Outcome, Rational, Real, Refusal, StrictOrder,
};
use catallax_shapes::{
    BSplineCurve, Circle, Curve, CurveId, Ellipse, Line, Point3, PointId, RVector3, Surface,
    Vector3,
};
use catallax_topology::numerical::tolerance;
use catallax_topology::{BodyId, BodyKind, EdgeId, FaceId, Model, Sense};

use super::float::{
    conic_point, double, native_parameter, rational, root, shifted_knots, unit_direction,
};
use super::{ContourFit, Contours, contour_generators, same_point};
use crate::ApplicationBody;
use crate::queries::{QueryRefusal, query, sign};

mod cross;
mod exact;
pub(super) mod geo;
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

/// Code of the mark a [`HiddenLineView`] carries when a piece of a fitted
/// edge was too short to classify by itself: no point beside it clears the
/// fit's tolerance within the piece's own ends, so it took the visibility of
/// the piece before it. Its bound is the longest such piece: a measured
/// value, not a proof.
pub const VIEW_SHORT_PIECE_MARK: &str = "app.view-short-piece";

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
    /// [`super::ContourArc::fit`], or of a fitted edge's carrier: an edge on
    /// an intersection branch, patch section or patch pair.
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
    /// The fit of a traced curve's fitted carrier.
    fit: Option<ContourFit>,
}

/// The carrier of an exact B-spline already in `scratch`.
fn spline_carrier(
    scratch: &Model,
    curve: CurveId,
    spline: &BSplineCurve,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<Carrier, Refusal> {
    let knots = shifted_knots(spline).ok_or_else(|| {
        Refusal::new(
            "app.view-curve",
            "two B-spline knots are closer than binary64 separates",
        )
    })?;
    let mut ctrl = Vec::with_capacity(spline.control_count());
    for i in 0..spline.control_count() {
        let p = point64(spline.control(i), floor, dec)?;
        let w = double(spline.weight(i), floor, dec)?;
        ctrl.push([w * p[0], w * p[1], w * p[2], w]);
    }
    let geo = Geo::Spline {
        spans: geo::decompose(spline.degree(), &knots, &ctrl),
    };
    let closed = curve_is_closed(scratch.geom(), curve, dec)?;
    Ok(Carrier {
        curve,
        form: Form::Spline,
        axes: None,
        geo,
        closed,
        fit: None,
    })
}

/// An intersection branch, patch section or patch pair between `ends`,
/// traced and fitted as its drawing is (`drawing::traced`), held as an
/// exact cubic B-spline whose controls are the fit's binary64 values. The
/// fit's knots and controls are binary64 numbers, so the carrier is exactly
/// the fitted curve.
fn fitted_carrier(
    scratch: &mut Model,
    curve: CurveId,
    ends: Option<[PointId; 2]>,
    floor: &Rational,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<Carrier, Refusal> {
    // Two ends at one point ask for the whole closed curve from there.
    let ends = match ends {
        Some([a, b]) if same_point(scratch, a, b, dec)? => Some([a, a]),
        other => other,
    };
    let drawn = super::drawing::traced(scratch, curve, ends, dec)?;
    marks.extend(&drawn.marks);
    let super::DrawingCurve::Nurbs {
        degree,
        knots,
        points,
        ..
    } = drawn.value.curve
    else {
        unreachable!("a traced curve is fitted as a NURBS");
    };
    let knots = knots
        .iter()
        .map(|&k| rational(k))
        .collect::<Result<Vec<_>, _>>()?;
    let mut controls = Vec::with_capacity(points.len());
    for p in &points {
        controls.push(Point3::new(
            Real::from(rational(p[0])?),
            Real::from(rational(p[1])?),
            Real::from(rational(p[2])?),
        ));
    }
    let weights = vec![Real::one(); controls.len()];
    let spline = BSplineCurve::new(degree, knots, controls, weights, dec)?;
    let plain = scratch.geom_mut().add_curve(Curve::BSpline(spline.clone()));
    let mut c = spline_carrier(scratch, plain, &spline, floor, dec)?;
    c.fit = drawn.value.fit;
    Ok(c)
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

/// A plain carrier for `curve` between `ends`: its closed form as a line,
/// circle, ellipse or B-spline. Intersection components are read through
/// their closed form, so later curve-plane meetings take the direct conic
/// paths. Intersection branches, patch sections and patch pairs have none,
/// and take a fitted B-spline over exactly the span (`fitted_carrier`).
#[allow(clippy::too_many_arguments)]
fn carrier(
    scratch: &mut Model,
    curve: CurveId,
    ends: Option<[PointId; 2]>,
    d: &Vector3,
    floor: &Rational,
    marks: &mut MarkLog,
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
                fit: None,
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
                fit: None,
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
                fit: None,
            })
        }
        CurveForm::BSpline(spline) => {
            let plain = if matches!(scratch.geom().curve(curve), Curve::BSpline(_)) {
                curve
            } else {
                scratch.geom_mut().add_curve(Curve::BSpline(spline.clone()))
            };
            spline_carrier(scratch, plain, &spline, floor, dec)
        }
        CurveForm::Branch | CurveForm::PatchSection(_) | CurveForm::PatchPair(_) => {
            fitted_carrier(scratch, curve, ends, floor, marks, dec)
        }
        _ => Err(Refusal::new(
            "app.view-curve",
            "hidden-line views of edges on parabolas and hyperbolas are unfinished",
        )),
    }
}

/// Builds a stroke on `curve` between `ends`, with its binary64 box and
/// whether it draws anything: a line along the view projects to a point.
#[allow(clippy::too_many_arguments)]
fn stroke(
    scratch: &mut Model,
    source: StrokeSource,
    curve: CurveId,
    ends: Option<[PointId; 2]>,
    sheet: &Sheet,
    d: &Vector3,
    floor: &Rational,
    marks: &mut MarkLog,
    dec: &mut Decider,
) -> Result<(Stroke, Box3, bool), Refusal> {
    let c = carrier(scratch, curve, ends, d, floor, marks, dec)?;
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
    let mut box3 = geo::box3(&c.geo, range[0], range[1]);
    // The true curve lies within the fit's tolerance of a fitted carrier.
    if let Some(fit) = &c.fit {
        for k in 0..3 {
            box3[0][k] -= fit.tolerance_mm;
            box3[1][k] += fit.tolerance_mm;
        }
    }
    Ok((
        Stroke {
            source,
            drawn: on_sheet,
            occluder: on_sheet,
            faces: Vec::new(),
            seam: false,
            fit: c.fit,
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
            conic_point(center, [u, w], t)
        }
        Form::Spline => {
            let Curve::BSpline(curve) = scratch.geom().curve(s.curve) else {
                unreachable!("a spline stroke's carrier is a B-spline");
            };
            curve.point_at(&Real::from(native_parameter(curve, t)?), dec)
        }
    }
}

/// `t` moved into a closed spline carrier's domain: the last piece of a
/// closed stroke cut into pieces runs past the domain's end.
fn within(s: &Stroke, t: f64) -> f64 {
    let [lo, hi] = s.geo.whole();
    if s.closed && matches!(s.geo, Geo::Spline { .. }) && t > hi {
        t - (hi - lo)
    } else {
        t
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
        return exact_point(scratch, s, within(s, mid), dec);
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

/// Where the ray of a fitted edge's piece starts, from the binary64 image
/// `p64` of its sample point at the middle of the binary64 interval
/// `(t0, t1)`: an exact point of one of the edge's faces just beside it,
/// moved across the edge by an offset from `start` up to `limit` (see
/// [`classify`]), on a face that faces the viewer there, or failing that one
/// seen edge on; `None` when every face of the edge faces away there, so
/// that they hide the piece themselves. `sides` keeps the side of the edge
/// each face lies on, as `surface::across` confirms it, for the edge's later
/// pieces.
///
/// SMLib's `SmHCR` decides an edge the same way. `IsBackFacingEdge` hides an
/// edge whose faces both face away, and for a concave edge whose faces face
/// opposite ways it tests whether the view runs into the solid's sector; an
/// edge it keeps is classified by a ray that stops short of the edge's own
/// faces by an offset. OpenCASCADE's `HLRBRep_Data` skips the edge's own
/// face instead. Here the ray starts on the face, exactly: from a facing
/// face it leaves the solid at once, then meets the other face just beside
/// the edge exactly when the view runs into a concave edge's sector, so the
/// sector test needs no separate code. A face seen edge on keeps the ray
/// along it, where the ray's own grazing rules decide.
#[allow(clippy::too_many_arguments)]
fn edge_start(
    scratch: &mut Model,
    s: &Stroke,
    (t0, t1): (f64, f64),
    p64: V3,
    (start, limit): (f64, f64),
    scene: &ray::Scene,
    sides: &mut std::collections::BTreeMap<FaceId, f64>,
    d: &Vector3,
    floor: &Rational,
    dec: &mut Decider,
) -> Result<Option<Point3>, Refusal> {
    let (_, tangent) = s.geo.jet(within(s, (t0 + t1) / 2.0));
    let StrokeSource::Edge(edge) = s.source else {
        unreachable!("only edge strokes start beside their faces");
    };
    let mut edge_on = None;
    for &face in &s.faces {
        // The face's other boundary edges, and whether it has a point
        // boundary: where none comes near, the face beside the edge is the
        // side the topology names.
        let mut others = Vec::new();
        let mut point_boundary = false;
        for &lp in &scratch.face(face).loops {
            match scratch.loop_(lp) {
                catallax_topology::Loop::Vertex(_) => point_boundary = true,
                catallax_topology::Loop::Coedges(coedges) => {
                    for &c in coedges {
                        let e = scratch.coedge(c).edge;
                        if e != edge
                            && let Some(&b) = scene.edges.get(&e)
                        {
                            others.push(b);
                        }
                    }
                }
            }
        }
        let others = (!point_boundary).then_some(others.as_slice());
        let known = sides.get(&face).copied();
        let (q, facing, side) = surface::across(
            scratch, face, p64, tangent, known, others, d, start, limit, floor, dec,
        )?;
        sides.insert(face, side);
        if facing > 0 {
            return Ok(Some(q));
        }
        if facing == 0 && edge_on.is_none() {
            edge_on = Some(q);
        }
    }
    Ok(edge_on)
}

/// What [`classify`] finds for one piece.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Label {
    /// Whether the piece is visible, and whether it is smooth.
    Seen(bool, bool),
    /// A fitted edge's piece too short to start its ray beside it within
    /// its own ends, with its length: it takes the label of the piece
    /// before it ([`carry_labels`]).
    Short(f64),
}

/// The label of one piece, over the binary64 interval `part`, at its sample
/// point `p`; `sides` as for [`edge_start`].
#[allow(clippy::too_many_arguments)]
fn classify(
    scratch: &mut Model,
    s: &Stroke,
    part: (f64, f64),
    p: &Point3,
    scene: &ray::Scene,
    sides: &mut std::collections::BTreeMap<FaceId, f64>,
    d: &Vector3,
    d64: V3,
    floor: &Rational,
    marks: &mut MarkLog,
    numeric: &mut bool,
    dec: &mut Decider,
) -> Result<Label, Refusal> {
    // A fitted edge is never smooth: smoothness compares its faces' normals
    // at an exact point of the edge, and a fitted edge has none at hand.
    if let (StrokeSource::Edge(_), Some(fit)) = (s.source, &s.fit) {
        let p64 = point64(p, floor, dec)?;
        let ends = [part.0, part.1].map(|t| geo::distance(p64, s.geo.point(within(s, t))));
        // The carrier lies within the fit's tolerance of the edge (rule
        // T4), and its binary64 images within the floor (rule T1): a point
        // twice as far across the edge lies on the face's own side. The
        // offset grows no farther from the edge than the piece's own ends,
        // where an outline may change what the piece sees; a piece whose
        // ends lie nearer than the first offset has no such point.
        let start = 2.0 * fit.tolerance_mm.max(scene.margin);
        let limit = ends[0].min(ends[1]);
        if limit < start {
            return Ok(Label::Short(ends[0] + ends[1]));
        }
        let Some(q) = edge_start(
            scratch,
            s,
            part,
            p64,
            (start, limit),
            scene,
            sides,
            d,
            floor,
            dec,
        )?
        else {
            return Ok(Label::Seen(false, false));
        };
        let q64 = point64(&q, floor, dec)?;
        let hidden = ray::hidden(scratch, scene, &q, q64, d, d64, marks, numeric, dec)?;
        return Ok(Label::Seen(!hidden, false));
    }
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
    Ok(Label::Seen(!hidden, smooth))
}

/// The visibility and smoothness of a stroke's pieces from their labels. A
/// short piece takes the label of the piece before it; on an open stroke,
/// pieces before the first decided one take its label, and on a closed one,
/// the last decided one's. `longest` keeps the longest short piece. Refuses
/// when no piece is decided.
fn carry_labels(
    labels: &[Label],
    closed: bool,
    longest: &mut Option<f64>,
) -> Result<Vec<(bool, bool)>, Refusal> {
    let seen = |label: &Label| match *label {
        Label::Seen(visible, smooth) => Some((visible, smooth)),
        Label::Short(_) => None,
    };
    let carried = if closed {
        labels.iter().rev().find_map(seen)
    } else {
        labels.iter().find_map(seen)
    };
    let Some(mut carry) = carried else {
        return Err(Refusal::new(
            "app.view-sample",
            "every piece of a fitted edge is too short to start a ray beside it",
        ));
    };
    let mut out = Vec::with_capacity(labels.len());
    for label in labels {
        match *label {
            Label::Seen(visible, smooth) => carry = (visible, smooth),
            Label::Short(length) => *longest = Some(longest.map_or(length, |l| l.max(length))),
        }
        out.push(carry);
    }
    Ok(out)
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
/// sections and patch pairs lie on fitted carriers, keep their fit, are
/// never smooth, and the answer carries a [`super::DRAWING_FIT_MARK`] mark
/// with the largest tolerance. A piece of such an edge too short to classify
/// within the fit's tolerance takes the visibility of the piece before it,
/// and the answer carries a [`VIEW_SHORT_PIECE_MARK`] mark; an edge with no
/// longer piece refuses with `app.view-sample`. Edges on parabolas and
/// hyperbolas refuse with
/// `app.view-curve`, as do the faces [`contour_generators`] refuses; a body
/// that is not one solid refuses with `app.view-solid`.
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
        // The fit with the largest tolerance among fitted edges.
        let mut fitted_edges: Option<ContourFit> = None;
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
                &mut marks,
                dec,
            )?;
            edge_boxes.insert(edge, box3);
            if let Some(fit) = &s.fit
                && fitted_edges
                    .as_ref()
                    .is_none_or(|f: &ContourFit| f.tolerance_mm < fit.tolerance_mm)
            {
                fitted_edges = Some(fit.clone());
            }
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
                &mut marks,
                dec,
            )?;
            if !on_sheet {
                continue;
            }
            // A generator along an edge is drawn as that edge; it still
            // bounds what it hides.
            s.drawn = arc.on_faces.len() == 1;
            s.faces = arc.on_faces.clone();
            if arc.fit.is_some() {
                s.fit = arc.fit.clone();
            }
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
            edges: edge_boxes,
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
                // Two stretches of one fitted carrier each lie within the
                // fit's tolerance of the curve it stands for.
                let reach = margin + 2.0 * a.fit.as_ref().map_or(0.0, |f| f.tolerance_mm);
                for c in cross::crossings(
                    &sheet, &a.geo, a.range, &a.geo, a.range, true, reach, &mut work,
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
                // A fitted carrier lies within its tolerance of the curve
                // it stands for: an edge, or a generator whose ends lie on
                // edges.
                let reach = margin
                    + a.fit.as_ref().map_or(0.0, |f| f.tolerance_mm)
                    + b.fit.as_ref().map_or(0.0, |f| f.tolerance_mm);
                let (x, y) = (a.sheet, b.sheet);
                if x[0] > y[2] + reach
                    || y[0] > x[2] + reach
                    || x[1] > y[3] + reach
                    || y[1] > x[3] + reach
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
                        let found = cross::crossings(
                            &sheet, &a.geo, a.range, &b.geo, b.range, false, reach, &mut work,
                        );
                        // An end of another stroke on a fitted edge, within
                        // reach in space, is a T-junction: an outline ending
                        // on the edge, whose projection the edge's touches
                        // there. The leaf search can pass such a tangent
                        // meeting by, so the end cuts the carrier itself.
                        if a.fit.is_some() && matches!(a.source, StrokeSource::Edge(_)) {
                            for e in b.ends.into_iter().flatten() {
                                let e64 = point64(&scratch.geom().position(e), &fine, dec)?;
                                let t = a.geo.param(e64);
                                let gap = geo::distance(a.geo.point(t), e64);
                                if gap <= reach {
                                    let p = exact_point(&scratch, a, t, dec)?;
                                    let p = scratch
                                        .geom_mut()
                                        .add_point(catallax_shapes::PointDef::Explicit(p));
                                    candidates[i].push((p, Some(gap)));
                                }
                            }
                        }
                        // A fitted carrier's exact ends that lie on A, a
                        // fitted edge's vertices or a fitted generator's
                        // boundary points, stand for the binary64 cuts near
                        // them.
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
        let mut short: Option<f64> = None;
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
            let mut sides = std::collections::BTreeMap::new();
            for &part in &parts {
                dec.check_deadline()?;
                let p = sample(&mut scratch, s, part, margin, dec)?;
                labels.push(classify(
                    &mut scratch,
                    s,
                    (part.2, part.3),
                    &p,
                    &scene,
                    &mut sides,
                    &d,
                    d64,
                    &fine,
                    &mut marks,
                    &mut numeric,
                    dec,
                )?);
            }
            let closed_loop = s.ends.is_none() || s.anchored;
            let labels = carry_labels(&labels, closed_loop, &mut short)?;
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
        if let Some(length) = short {
            marks.push(Mark::new(
                VIEW_SHORT_PIECE_MARK,
                BoundKind::Distance,
                rational(length.next_up())?,
                "a piece of a fitted edge too short to start a ray beside it took the visibility \
                 of the piece before it: the longest such piece, measured in binary64, not a \
                 proof",
            ));
        }
        if let Some(fit) = &fitted_edges {
            marks.push(super::drawing::fit_mark(fit)?);
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
