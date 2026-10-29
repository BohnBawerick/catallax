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
    Apart, CurveId, GeomStore, Point3, PointId, SurfaceId, SurfaceKind, Vector3,
};
use std::{
    collections::VecDeque,
    sync::{Arc, Mutex},
};

mod algebraic_solve;
mod arc;
mod chart;
mod hits;
mod math;
mod pair;
mod query;
mod solve;

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
    /// An intersection branch with no conic closed form.
    Branch,
}

type JetKey = (Rational, Rational, u32, usize);
type JetCache = VecDeque<(JetKey, Vec<[Ball; 3]>)>;

/// One parametrised piece of a traced arc.
#[derive(Clone, Debug)]
pub struct TracePiece {
    coords: [chart::E; 3],
    jets: Arc<Mutex<JetCache>>,
}

impl TracePiece {
    fn new(coords: [chart::E; 3]) -> Self {
        Self {
            coords,
            jets: Arc::new(Mutex::new(VecDeque::new())),
        }
    }

    fn enclose_jets(
        &self,
        t: &Ball,
        order: usize,
        dec: &Decider,
    ) -> Result<Vec<[Ball; 3]>, Refusal> {
        dec.check_deadline()?;
        let key = (t.mid(), t.rad(), t.precision(), order);
        let mut cache = self.jets.lock().expect("trace Taylor enclosure lock");
        dec.check_deadline()?;
        if let Some((_, jets)) = cache.iter().find(|(input, _)| input == &key) {
            let result = jets.clone();
            dec.check_deadline()?;
            return Ok(result);
        }
        let jets = chart::ball_jets(&self.coords, t, order, dec)?;
        // Retain a bounded recent working set without flushing all the panels
        // when the cache fills. Eviction changes cost, never an enclosure.
        if cache.len() == 128 {
            cache.pop_front();
        }
        cache.push_back((key, jets.clone()));
        dec.check_deadline()?;
        Ok(jets)
    }

    /// Parameter range of the piece.
    pub fn range(&self) -> (Rational, Rational) {
        (Rational::zero(), Rational::one())
    }

    /// Proven enclosure of the point at parameter `t`, or a deadline refusal.
    pub fn evaluate(&self, t: &Ball, dec: &Decider) -> Result<[Ball; 3], Refusal> {
        Ok(self.enclose_jets(t, 0, dec)?.remove(0))
    }

    /// Exact algebraic point at a parameter in `[0, 1]`.
    /// Retains incidence when display rounding collapses distinct samples.
    pub fn evaluate_exact(&self, t: &Rational, dec: &Decider) -> Result<Point3, Refusal> {
        chart::exact_position(&self.coords, t, dec)
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
        chart::interval_jets(&self.coords, t, order, dec)
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

/// Whether these kinds have an intersection path; individual poses can still refuse.
pub fn supports(a: SurfaceKind, b: SurfaceKind) -> bool {
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
    hits::intersect(geom, c, s, dec)
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
