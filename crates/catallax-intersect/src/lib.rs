//! Piece 4 of catallax: exact intersections.
//!
//! Plane/cylinder/sphere intersections and closed-form cone/torus sections,
//! stored by definition; curve-surface meetings, branch-aware membership,
//! ordering and tangents; certified finite curve traces.
//!
//! The contract is `docs/pieces/04-intersections.md`.

use catallax_numbers::{
    Answered, Ball, Decider, Interval, Outcome, Rational, Refusal, StrictOrder,
};
use catallax_shapes::{
    Apart, CurveId, GeomStore, Point3, PointId, SplineSupport, Surface, SurfaceId, SurfaceKind,
    Vector3,
};
use std::{
    collections::VecDeque,
    sync::{Arc, Mutex},
};

mod algebraic_solve;
mod arc;
mod bezier_curve;
mod chart;
mod fillet_cut;
mod hits;
mod math;
mod pair;
mod query;
mod solve;
mod span_hull;
mod spline_arithmetic;
mod spline_germ;
mod spline_line;
pub use spline_germ::native_spline_normal_acceleration;
mod spline_plane;
pub use spline_line::native_spline_line_points;

/// How two objects meet along a component or at a point.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Contact {
    /// They cross.
    Crossing,
    /// They touch without crossing.
    Tangent,
}

/// One connected piece of where two surfaces meet.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Component {
    /// A curve.
    Curve {
        /// The curve, closed form or intersection branch.
        curve: CurveId,
        /// Crossing or tangent along the curve.
        contact: Contact,
        /// Singular points of the curve, as vertices.
        singular: Vec<PointId>,
    },
    /// An isolated touching point.
    Point {
        /// The point.
        point: PointId,
    },
}

/// Where two surfaces meet.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum SurfacePair {
    /// Nowhere.
    Disjoint,
    /// They are the same point set.
    Coincident {
        /// Whether their normals agree.
        same_sense: bool,
    },
    /// In these components, in the order of the piece description.
    Meets(Vec<Component>),
}

/// One meeting point along a curve.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CurveHit {
    /// The point.
    pub point: PointId,
    /// Crossing or touching.
    pub contact: Contact,
}

/// Where a curve meets a surface.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum CurveSurface {
    /// Nowhere.
    Disjoint,
    /// The curve lies in the surface.
    Contained,
    /// At these points, along the curve's natural direction.
    Meets(Vec<CurveHit>),
}

/// Cyclic order of three points on a closed curve.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Cyclic {
    /// From the first, the second comes before the third in the natural direction.
    Forward,
    /// From the first, the third comes before the second.
    Backward,
}

/// The exact shape of a curve's carrier.
#[non_exhaustive]
#[derive(Clone, Debug)]
pub enum CurveForm {
    /// A straight line through `origin` along `dir`.
    Line {
        /// A point of the line.
        origin: Point3,
        /// Its natural direction.
        dir: Vector3,
    },
    /// A circle, counter-clockwise seen from the tip of `normal`.
    Circle {
        /// Centre.
        center: Point3,
        /// Normal of its plane.
        normal: Vector3,
        /// Positive squared radius.
        radius_sq: catallax_numbers::Real,
        /// Direction of the point at `t = 0`.
        ref_dir: Vector3,
    },
    /// An ellipse, counter-clockwise seen from the tip of `normal`.
    Ellipse {
        /// Centre.
        center: Point3,
        /// Normal of its plane.
        normal: Vector3,
        /// Direction of the major axis and of the point at `t = 0`.
        major_dir: Vector3,
        /// Positive squared semi-major axis.
        major_sq: catallax_numbers::Real,
        /// Positive squared semi-minor axis.
        minor_sq: catallax_numbers::Real,
    },
    /// A parabola `y² = 4 focal x` in its oriented plane.
    Parabola {
        /// Vertex.
        vertex: Point3,
        /// Plane normal; fixes the natural direction at the vertex.
        normal: Vector3,
        /// Direction from the vertex towards the focus.
        axis_dir: Vector3,
        /// Positive distance from the vertex to the focus.
        focal: catallax_numbers::Real,
    },
    /// One selected branch of `x²/transverse_sq - y²/conjugate_sq = 1`.
    Hyperbola {
        /// Centre of the two-branch carrier.
        center: Point3,
        /// Plane normal; fixes the natural direction at the vertex.
        normal: Vector3,
        /// Direction of the transverse axis.
        axis_dir: Vector3,
        /// Positive squared transverse semi-axis.
        transverse_sq: catallax_numbers::Real,
        /// Positive squared conjugate semi-axis.
        conjugate_sq: catallax_numbers::Real,
        /// Sign of the transverse coordinate on this component.
        branch: catallax_numbers::StrictSign,
    },
    /// A finite positive-weight spline, with its actual knot domain and controls.
    /// The tag does not assert that the curve is embedded or bounds a valid face.
    BSpline(catallax_shapes::BSplineCurve),
    /// Exact patch-parameter branch with a complete regular cell cover.
    PatchSection(catallax_shapes::PatchSectionBranch),
    /// An intersection branch with no conic closed form.
    Branch,
    /// Exact meeting selected by a complete two-source parameter cover.
    PatchPair(catallax_shapes::PatchPairCurve),
}

type JetKey = (Rational, Rational, u32, usize);
type JetCache = VecDeque<(JetKey, Vec<[Ball; 3]>)>;

#[derive(Clone, Debug)]
enum TraceCoordinates {
    Explicit(ExplicitTrace),
    Section(Arc<catallax_shapes::SectionTracePiece>),
    PatchPair(Arc<catallax_shapes::PairTracePiece>),
    FilletCut(Arc<catallax_shapes::FilletCutTrace>, bool),
}

#[derive(Clone, Debug)]
struct ExplicitTrace {
    coords: [chart::E; 3],
    jets: Arc<Mutex<JetCache>>,
    program: Arc<std::sync::OnceLock<chart::JetProgram>>,
    span: Option<Arc<span_hull::SpanHull>>,
    // Work-span observation only: one generation per cached jet, in the
    // cache's order.
    #[cfg(feature = "work-span")]
    jet_generations: Arc<Mutex<VecDeque<u64>>>,
}

/// One parametrised piece of a traced arc.
#[derive(Clone, Debug)]
pub struct TracePiece {
    coords: TraceCoordinates,
}

impl TracePiece {
    fn new(coords: [chart::E; 3]) -> Self {
        Self {
            coords: TraceCoordinates::Explicit(ExplicitTrace {
                coords,
                jets: Arc::new(Mutex::new(VecDeque::new())),
                program: Arc::new(std::sync::OnceLock::new()),
                span: None,
                #[cfg(feature = "work-span")]
                jet_generations: Arc::new(Mutex::new(VecDeque::new())),
            }),
        }
    }

    fn section(piece: catallax_shapes::SectionTracePiece) -> Self {
        Self {
            coords: TraceCoordinates::Section(Arc::new(piece)),
        }
    }

    /// Use a producer-owned four-parameter implicit graph, never a fitted curve.
    pub fn from_patch_pair(piece: catallax_shapes::PairTracePiece) -> Self {
        Self {
            coords: TraceCoordinates::PatchPair(Arc::new(piece)),
        }
    }

    fn spline(
        curve: &catallax_shapes::BSplineCurve,
        span: usize,
        from: catallax_numbers::Real,
        to: catallax_numbers::Real,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let coords = chart::spline_span(curve, span, from.clone(), to.clone(), dec)?;
        let mut piece = Self::new(coords);
        if let TraceCoordinates::Explicit(definition) = &mut piece.coords {
            definition.span = Some(Arc::new(span_hull::SpanHull::new(curve, span, [from, to])));
        }
        Ok(piece)
    }

    fn enclose_jets(
        &self,
        t: &Ball,
        order: usize,
        dec: &Decider,
    ) -> Result<Vec<[Ball; 3]>, Refusal> {
        match &self.coords {
            TraceCoordinates::Explicit(definition) => definition.enclose_jets(t, order, dec),
            TraceCoordinates::Section(piece) => piece.taylor(t, order, dec),
            TraceCoordinates::PatchPair(piece) => piece.taylor(t, order, dec),
            TraceCoordinates::FilletCut(cut, forward) => {
                dec.check_deadline()?;
                if order != 0 {
                    return Err(Refusal::new(
                        "intersect.fillet-cut-derivative",
                        "a cap-cut trace encloses positions; exact derivative enclosures are not implemented",
                    ));
                }
                let point = if *forward {
                    cut.enclose(t, dec)?
                } else {
                    let parameter = &Ball::from_rational(&Rational::one(), t.precision()) - t;
                    cut.enclose(&parameter, dec)?
                };
                Ok(vec![point])
            }
        }
    }

    /// Parameter range of the piece.
    pub fn range(&self) -> (Rational, Rational) {
        (Rational::zero(), Rational::one())
    }

    /// Proven enclosure of the point at parameter `t`, or a deadline refusal.
    pub fn evaluate(&self, t: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        Ok(self.enclose_jets(t, 0, dec)?.remove(0))
    }

    /// A proven box round the piece over all of `t`. Where [`Self::evaluate`]
    /// reads a spline span's polynomial, which can overstate a wide `t` by orders
    /// of magnitude, this cuts down its controls, so the box follows the curve.
    pub fn enclose(&self, t: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        if let TraceCoordinates::Explicit(ExplicitTrace {
            span: Some(span), ..
        }) = &self.coords
            && t.lo() < t.hi()
            && let Some(found) = span.enclose(t, dec)?
        {
            return Ok(found);
        }
        self.evaluate(t, dec)
    }

    /// Exact algebraic point at a parameter in `[0, 1]`.
    /// Retains incidence when display rounding collapses distinct samples.
    pub fn evaluate_exact(&self, t: &Rational, dec: &mut Decider) -> Result<Point3, Refusal> {
        match &self.coords {
            TraceCoordinates::Explicit(definition) => {
                chart::exact_position(&definition.coords, t, dec)
            }
            TraceCoordinates::Section(piece) => piece.point_at(t, dec),
            TraceCoordinates::PatchPair(piece) => piece.point_at(t, dec),
            TraceCoordinates::FilletCut(cut, forward) => {
                if t == &Rational::zero() || t == &Rational::one() {
                    return Ok(cut.endpoint(usize::from((t == &Rational::one()) == *forward)));
                }
                Err(Refusal::new(
                    "intersect.fillet-cut-exact-point",
                    "an interior cap-cut point needs an exact root representation",
                ))
            }
        }
    }

    /// Proven enclosure of the derivative at parameter `t`, or a deadline refusal.
    pub fn derivative(&self, t: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        Ok(self.enclose_jets(t, 1, dec)?.remove(1))
    }

    /// Proven enclosures of the Taylor coefficients `C^(k)(s) / k!` for
    /// `k = 0..=order`, valid for every `s` in `t`, or a deadline refusal.
    pub fn taylor(&self, t: &Ball, order: usize, dec: &Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        self.enclose_jets(t, order, dec)
    }

    /// Outward-rounded binary64 Taylor enclosures, using the same trace
    /// definition and recurrence as [`Self::taylor`]. `None` means the filter
    /// could not enclose a denominator or finite result; use `taylor` then.
    pub fn taylor_interval(
        &self,
        t: Interval,
        order: usize,
        dec: &Decider,
    ) -> Result<Option<Vec<[Interval; 3]>>, Refusal> {
        dec.check_deadline()?;
        let explicit = match &self.coords {
            TraceCoordinates::Explicit(explicit) => explicit,
            TraceCoordinates::FilletCut(..) => return Ok(None),
            TraceCoordinates::Section(_) | TraceCoordinates::PatchPair(_) => {
                let (Some(lo), Some(hi)) = (
                    Rational::from_f64_exact(t.lo()),
                    Rational::from_f64_exact(t.hi()),
                ) else {
                    return Ok(None);
                };
                let bits = dec.policy().max_precision_bits.min(128);
                let parameter =
                    Ball::from_rational(&lo, bits).hull(&Ball::from_rational(&hi, bits));
                return Ok(Some(
                    self.taylor(&parameter, order, dec)?
                        .into_iter()
                        .map(|xyz| xyz.each_ref().map(Ball::to_interval))
                        .collect(),
                ));
            }
        };
        let program = explicit
            .program
            .get_or_init(|| chart::JetProgram::new(&explicit.coords));
        chart::interval_jets(program, t, order, dec)
    }
}

impl ExplicitTrace {
    fn enclose_jets(
        &self,
        t: &Ball,
        order: usize,
        dec: &Decider,
    ) -> Result<Vec<[Ball; 3]>, Refusal> {
        dec.check_deadline()?;
        if !catallax_numbers::shared_caches() {
            // A worker leaves the shared panel cache unfilled for the writer.
            let program = self
                .program
                .get_or_init(|| chart::JetProgram::new(&self.coords));
            let jets = chart::ball_jets(program, t, order, dec)?;
            dec.check_deadline()?;
            return Ok(jets);
        }
        let key = (t.mid(), t.rad(), t.precision(), order);
        let mut cache = self.jets.lock().expect("trace Taylor enclosure lock");
        dec.check_deadline()?;
        if let Some(position) = cache.iter().position(|(input, _)| input == &key) {
            #[cfg(feature = "work-span")]
            if let Some(&generation) = self.jet_generations.lock().expect("lock").get(position) {
                catallax_numbers::work_span::read_at("shared cache: trace jets", generation, 0);
            }
            let result = cache[position].1.clone();
            dec.check_deadline()?;
            return Ok(result);
        }
        let program = self
            .program
            .get_or_init(|| chart::JetProgram::new(&self.coords));
        let jets = chart::ball_jets(program, t, order, dec)?;
        if cache.len() == 128 {
            cache.pop_front();
            #[cfg(feature = "work-span")]
            self.jet_generations.lock().expect("lock").pop_front();
        }
        cache.push_back((key, jets.clone()));
        #[cfg(feature = "work-span")]
        {
            let generation = catallax_numbers::work_span::Instance::new().id();
            self.jet_generations
                .lock()
                .expect("lock")
                .push_back(generation);
            catallax_numbers::work_span::write_at("shared cache: trace jets", generation, 0);
        }
        dec.check_deadline()?;
        Ok(jets)
    }
}

/// A certified trace of a curve arc, in pieces, in the natural direction.
#[derive(Clone, Debug)]
pub struct CurveTrace {
    pieces: Vec<TracePiece>,
}

impl CurveTrace {
    /// The pieces in order along the arc.
    pub fn pieces(&self) -> &[TracePiece] {
        &self.pieces
    }
}

/// Whether these retained carriers have an intersection path; individual poses can still refuse.
pub fn supports(a: &Surface, b: &Surface) -> bool {
    match (a, b) {
        (
            analytic,
            Surface::BSpline(SplineSupport::RationalPatch(_) | SplineSupport::RationalAtlas(_)),
        )
        | (
            Surface::BSpline(SplineSupport::RationalPatch(_) | SplineSupport::RationalAtlas(_)),
            analytic,
        ) => {
            return catallax_shapes::SectionSurface::of(analytic).is_some();
        }
        _ => (),
    }
    let (a, b) = (a.kind(), b.kind());
    (matches!(
        a,
        SurfaceKind::Plane
            | SurfaceKind::Cylinder
            | SurfaceKind::Cone
            | SurfaceKind::Sphere
            | SurfaceKind::Torus
    ) && matches!(
        b,
        SurfaceKind::Plane
            | SurfaceKind::Cylinder
            | SurfaceKind::Cone
            | SurfaceKind::Sphere
            | SurfaceKind::Torus
    )) || matches!(
        (a, b),
        (SurfaceKind::Plane, SurfaceKind::EllipticCylinder)
            | (SurfaceKind::EllipticCylinder, SurfaceKind::Plane)
    )
}

/// Where two surfaces meet; new curves and points are added to `geom`.
pub fn intersect_surfaces(
    geom: &mut GeomStore,
    a: SurfaceId,
    b: SurfaceId,
    dec: &mut Decider,
) -> Answered<SurfacePair> {
    pair::construct(geom, a, b, dec)
}

/// Where a curve meets a surface; new points are added to `geom`.
pub fn intersect_curve_surface(
    geom: &mut GeomStore,
    c: CurveId,
    s: SurfaceId,
    dec: &mut Decider,
) -> Answered<CurveSurface> {
    hits::intersect(geom, c, s, false, dec)
}

/// Prefer direct conic-section constructions for curve-surface hits.
/// The exact point definitions can differ from the default solver's definitions.
pub fn intersect_curve_surface_framed(
    geom: &mut GeomStore,
    c: CurveId,
    s: SurfaceId,
    dec: &mut Decider,
) -> Answered<CurveSurface> {
    hits::intersect(geom, c, s, true, dec)
}

/// Uses a checked line parameter instead of eliminating world coordinates.
/// Other curve families keep their canonical constructions. A line's exact
/// point definitions can differ from the canonical solver's definitions.
pub fn intersect_curve_surface_with_line_parameter(
    geom: &mut GeomStore,
    c: CurveId,
    s: SurfaceId,
    dec: &mut Decider,
) -> Answered<CurveSurface> {
    hits::with_line_parameter(geom, c, s, dec)
}

/// Whether exact controls identify a complete section of a checked Bezier extrusion.
/// Incidence narrows the search but does not replace the control and domain proof.
pub fn is_bezier_extrusion_section(
    geom: &GeomStore,
    curve: CurveId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    Ok(bezier_curve::section(geom, curve, dec)?.is_some())
}

/// All finite meetings of a checked Bezier extrusion section with a whole line.
/// Returns `None` outside this curve family. Collinear overlap is refused by name.
pub fn bezier_section_line_points(
    geom: &GeomStore,
    section: CurveId,
    line: CurveId,
    dec: &mut Decider,
) -> Result<Option<Vec<Point3>>, Refusal> {
    bezier_curve::line_points(geom, section, line, dec)
}

/// Complete meetings of two checked finite Bezier sections.
#[derive(Clone, Debug)]
pub enum BezierSectionPair {
    /// The complete finite point sets coincide.
    Coincident,
    /// All shared endpoints after exact control-hull separation of the interiors.
    Points(Vec<Point3>),
}

/// Proves coincidence or separation except at shared endpoints. Interior
/// intersections and exhausted subdivision proofs are refused by name.
/// Returns `None` unless both curves have checked extrusion-section ownership.
pub fn bezier_section_meetings(
    geom: &GeomStore,
    a: CurveId,
    b: CurveId,
    dec: &mut Decider,
) -> Result<Option<BezierSectionPair>, Refusal> {
    bezier_curve::meetings(geom, a, b, dec)
}

/// Normal acceleration of a checked Bezier section for the supplied tangent.
/// The tangent must be nonzero and parallel to the curve derivative. Its scale
/// fixes the acceleration scale. Returns `None` outside this curve family.
pub fn bezier_section_normal_acceleration(
    geom: &GeomStore,
    section: CurveId,
    position: &Point3,
    tangent: &Vector3,
    dec: &mut Decider,
) -> Result<Option<Vector3>, Refusal> {
    bezier_curve::normal_acceleration(geom, section, position, tangent, dec)
}

/// The closed form of a curve's carrier.
pub fn curve_form(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<CurveForm, Refusal> {
    query::form(geom, c, dec)
}

/// Whether a curve closes on itself, rather than requiring ordered finite endpoints.
///
/// Lines, parabolas and hyperbola branches are open. Circles, ellipses and
/// the supported cylinder/sphere, cylinder/cylinder and parallel-axis
/// cylinder/torus and cone/torus branches are closed.
/// Unsupported curve definitions return the same refusal as [`curve_form`].
pub fn curve_is_closed(geom: &GeomStore, c: CurveId, dec: &mut Decider) -> Result<bool, Refusal> {
    if fillet_cut::of_curve(geom, c, dec)?.is_some() {
        return Ok(false);
    }
    if let Some(closed) = geom.curve(c).is_closed() {
        return Ok(closed);
    }
    match query::form(geom, c, dec)? {
        CurveForm::Line { .. } | CurveForm::Parabola { .. } | CurveForm::Hyperbola { .. } => {
            Ok(false)
        }
        CurveForm::Circle { .. } | CurveForm::Ellipse { .. } => Ok(true),
        CurveForm::BSpline(curve) => {
            let [lo, hi] = curve.domain();
            let first = curve.point_at(&lo.clone().into(), dec)?;
            let last = curve.point_at(&hi.clone().into(), dec)?;
            Ok(
                math::sign(&last.minus(&first).norm_squared(), dec)?
                    == catallax_numbers::Sign::Zero,
            )
        }
        CurveForm::PatchSection(branch) => Ok(branch.is_closed()),
        CurveForm::PatchPair(branch) => Ok(branch.is_closed()),
        CurveForm::Branch => {
            let plan = query::validate(geom, c, dec)?.expect("intersection branch");
            match &plan.frame {
                Some(
                    pair::Frame::CylinderSphere { .. }
                    | pair::Frame::Cylinders { .. }
                    | pair::Frame::CylinderTorus { .. }
                    | pair::Frame::ConeTorus(_),
                ) => Ok(true),
                None => Err(math::refused_curve()),
            }
        }
    }
}

/// Whether a point lies on the closure of the curve itself (for an
/// intersection curve, on its branch).
pub fn on_curve(geom: &GeomStore, c: CurveId, p: PointId, dec: &mut Decider) -> Outcome<Apart> {
    query::on(geom, c, p, dec)
}

/// Whether the closures of two curves are the same point set.
pub fn same_curve(geom: &GeomStore, a: CurveId, b: CurveId, dec: &mut Decider) -> Outcome<Apart> {
    query::same(geom, a, b, dec)
}

/// Order of two points along an open curve.
pub fn order_on_curve(
    geom: &GeomStore,
    c: CurveId,
    p: PointId,
    q: PointId,
    dec: &mut Decider,
) -> Outcome<StrictOrder> {
    query::order(geom, c, p, q, dec)
}

/// Cyclic order of three points on a closed curve.
pub fn cyclic_order(
    geom: &GeomStore,
    c: CurveId,
    points: [PointId; 3],
    dec: &mut Decider,
) -> Outcome<Cyclic> {
    query::cyclic(geom, c, points, dec)
}

/// Exact tangents in the natural direction at a point of the curve:
/// `[arriving, leaving]`, equal at regular points. Deciding whether the point
/// is singular can end not exact, which refuses.
pub fn tangent_at(
    geom: &GeomStore,
    c: CurveId,
    p: PointId,
    dec: &mut Decider,
) -> Result<[Vector3; 2], Refusal> {
    query::tangent(geom, c, p, dec)
}

/// An exact point strictly inside the arc from `from` to `to` (both `None`:
/// the whole closed curve).
pub fn arc_point(
    geom: &GeomStore,
    c: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Answered<Point3> {
    arc::point(geom, c, from, to, dec)
}

/// Certified trace of the arc from `from` to `to` (both `None`: the whole closed curve).
pub fn trace(
    geom: &GeomStore,
    c: CurveId,
    from: Option<PointId>,
    to: Option<PointId>,
    dec: &mut Decider,
) -> Answered<CurveTrace> {
    arc::trace(geom, c, from, to, dec)
}

/// Exact normal-offset profile quarter at an axial coordinate of its original
/// elliptic cylinder, before any spline approximation. Quarter `k` runs from
/// angle `k*pi/2` to `(k+1)*pi/2`, including both ends, for `k` in `0..4`.
/// Inward profiles at or beyond the normal reach refuse.
///
/// This traces the definition, not the represented boundary of a fitted body.
pub fn trace_elliptic_offset_quarter(
    definition: &catallax_shapes::EllipticOffsetDefinition,
    axial: &catallax_numbers::Real,
    quarter: usize,
    dec: &mut Decider,
) -> Answered<TracePiece> {
    let coords = chart::offset_profile_quarter(definition, axial, quarter, dec)?;
    Ok(catallax_numbers::Answer::exact(TracePiece::new(coords)))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod spline_curve_tests;

#[cfg(test)]
mod section_trace_tests;
