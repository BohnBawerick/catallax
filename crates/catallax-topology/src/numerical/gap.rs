//! The edge-gap checker of roadmap item 3.0A (decision Q6, setup A+).
//!
//! An edge of a numerical body has two copies: a 3D curve, and a face copy
//! lifted through the face's surface at the same parameter. This module finds
//! the largest distance between them, cheaply and without a pass or fail
//! limit of its own. The caller stores the result times a margin (rule T4).
//!
//! The method, with its sources:
//!
//! 1. Round both copies to binary64 in a local frame reached by an exact
//!    shift (rule T1; `lift.rs`).
//! 2. A straight residual has its maximum at an end. This is the closed form
//!    for a line on a plane or along a generator. A circle or an ellipse
//!    against a lift that is one too, at the same angular rate, has a conic
//!    residual: its stationary points are the roots of a quartic. Both are
//!    the pairs of OpenCASCADE's `ElementaryMaximum`.
//! 3. Otherwise cut the edge at every break: knots of the 3D curve, knots of
//!    the face copy, every crossing of the face copy with a surface knot
//!    line, seam or pole line, the surface's domain border, and any joint the
//!    builder reports. OpenCASCADE's `GeomLib_CheckCurveOnSurface` cuts at
//!    the same places since its October 2026 fix.
//! 4. On each piece take `max(17, 3 x composed degree)` samples at Chebyshev
//!    points of the first kind, plus both ends. The count depends only on
//!    degrees and is made even. No point of an even grid is a rational share
//!    of its piece, so none coincides with a builder that checks at eighths.
//! 5. Read the Chebyshev coefficient tail of the residual on each piece. It
//!    estimates what the samples leave unresolved. While that estimate does
//!    not fit inside the storage margin, the points of the piece are
//!    doubled. Every sample already taken stays in the search. The tail can
//!    only add points; it never ends the search early.
//! 6. Polish every local maximum of the samples with Brent's method on the
//!    true gap.
//! 7. Report the largest value, where it is, and the measurement error of
//!    rule T4: the binary64 rounding estimate plus what the last tail says
//!    the samples left unresolved. Nothing is escalated to higher precision
//!    here.
//!
//! Step 2 of section 4.2, imports whose copies run at different speeds, is
//! outside this search: `geometric.rs` reads the geometric gap of an edge use
//! both ways, `refit.rs` gives a face copy the parameter of its 3D curve, and
//! `numerical::admit_import` decides between them (decision Q9). The search
//! above then measures the copy that is kept. All of that draws on the one
//! work allowance of the edge use ([`Limits`]).
//!
//! A carrier with no kernel form, such as an imported surface a reader holds
//! in the file's own terms, is read through its own binary64 evaluator
//! ([`SampledCurve`], [`SampledSurface`]). The search is the same: the
//! carrier states its breaks and its degree, and where it has no degree the
//! least count of step 4 is taken and step 5 adds points.
//!
//! All arithmetic is IEEE basic operations in one thread, so the result is
//! the same bits on every machine and at any thread count.

use std::sync::Arc;

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_shapes::{BSplineCurve, BSplineSurface, Curve, Surface};

mod geometric;
mod lift;
mod numeric;
/// New face copies of an edge use: a refit from its 3D curve (rules T3, T7
/// and T8) and a change of parameter of the copy it has (decision Q9).
pub mod refit;

pub use geometric::{GeometricMeasure, measure_geometric};

use lift::{Curve64, Face64, Lifted, Revolved, Scratch, Surface64};

pub(in crate::numerical) use lift::defining_points;
pub(in crate::numerical) use numeric::{brent_maximum, chebyshev_shares, try_brent_maximum};

/// A face's surface in binary64, in a frame a whole body shares (rule T1),
/// with its parameters unshifted, for the distance readings of the
/// tolerance rules.
pub(in crate::numerical) struct FramedSurface {
    surface: Surface64,
    scratch: Scratch,
    /// The period of u on a surface of revolution.
    pub period_u: Option<f64>,
}

impl FramedSurface {
    /// Rounds `surface` in the frame shifted to `origin`.
    pub(in crate::numerical) fn new(
        surface: &Surface,
        origin: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let surface = lift::framed_surface(surface, origin, dec)?;
        let period_u =
            matches!(surface, Surface64::Revolved { .. }).then_some(std::f64::consts::TAU);
        Ok(Self {
            surface,
            scratch: Scratch::default(),
            period_u,
        })
    }

    /// The point at `uv`.
    pub(in crate::numerical) fn point(&mut self, uv: [f64; 2]) -> [f64; 3] {
        self.surface.point(uv, &mut self.scratch).value
    }

    /// The foot of `target` from a seed, and the evaluations spent.
    pub(in crate::numerical) fn foot(
        &mut self,
        target: &[f64; 3],
        seed: [f64; 2],
    ) -> ([f64; 2], u64) {
        self.surface.foot(target, seed, &mut self.scratch)
    }

    /// A box that holds the surface over the u/v rectangle `rect`, or `None`
    /// for a plane, whose faces lie inside the box of their edges. A surface
    /// of revolution is held by the box of its local frame's extent: a ring
    /// of its largest radius over the v range, a sphere or a torus whole. A
    /// spline is held by the controls whose basis functions are nonzero on
    /// the rectangle (Piegl and Tiller, "The NURBS Book", section 4.4, and the
    /// local support of section 2.2; positive weights).
    pub(in crate::numerical) fn bounding_box(&self, rect: [[f64; 2]; 2]) -> Option<[[f64; 2]; 3]> {
        let v_range = rect[1];
        let mut low = [f64::INFINITY; 3];
        let mut high = [f64::NEG_INFINITY; 3];
        let mut take = |p: [f64; 3]| {
            for axis in 0..3 {
                low[axis] = low[axis].min(p[axis]);
                high[axis] = high[axis].max(p[axis]);
            }
        };
        match &self.surface {
            Surface64::Plane { .. } => return None,
            Surface64::Sampled(_) => {
                unreachable!("a framed surface is rounded from a kernel surface")
            }
            Surface64::Spline(net) => {
                for control in net.controls_over(rect) {
                    let w = control[3];
                    take([control[0] / w, control[1] / w, control[2] / w]);
                }
            }
            Surface64::Revolved {
                kind,
                origin,
                frame,
                sizes,
            } => {
                let [v0, v1] = v_range;
                let local: [[f64; 2]; 3] = match kind {
                    Revolved::Cylinder => [[-sizes[0], sizes[0]], [-sizes[0], sizes[0]], [v0, v1]],
                    Revolved::EllipticCylinder => {
                        [[-sizes[0], sizes[0]], [-sizes[1], sizes[1]], [v0, v1]]
                    }
                    Revolved::Cone => {
                        let r = v0.abs().max(v1.abs()) * sizes[0].abs();
                        [[-r, r], [-r, r], [v0, v1]]
                    }
                    Revolved::Sphere => [[-sizes[0], sizes[0]]; 3],
                    Revolved::Torus => {
                        let r = sizes[0].abs() + sizes[1].abs();
                        [[-r, r], [-r, r], [-sizes[1].abs(), sizes[1].abs()]]
                    }
                };
                for a in local[0] {
                    for b in local[1] {
                        for c in local[2] {
                            take(std::array::from_fn(|k| {
                                origin[k] + a * frame[0][k] + b * frame[1][k] + c * frame[2][k]
                            }));
                        }
                    }
                }
            }
        }
        Some(std::array::from_fn(|axis| {
            let slack = BOX_ROUNDINGS
                * numeric::UNIT_ROUNDOFF
                * low[axis].abs().max(high[axis].abs()).max(f64::MIN_POSITIVE);
            [low[axis] - slack, high[axis] + slack]
        }))
    }
}

/// A retained u/v curve in binary64, for the tolerance rules' reading of a
/// face's boundary in its parameters.
pub(in crate::numerical) struct FramedUvCurve {
    curve: lift::Spline<3>,
    scratch: Scratch,
}

impl FramedUvCurve {
    /// Rounds `source`.
    pub(in crate::numerical) fn new(
        source: &BSplineCurve,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Ok(Self {
            curve: lift::uv_curve(source, dec)?,
            scratch: Scratch::default(),
        })
    }

    /// The controls in u/v: their hull holds the curve (Piegl and Tiller,
    /// "The NURBS Book", section 4.2; positive weights).
    pub(in crate::numerical) fn controls(&self) -> Vec<[f64; 2]> {
        self.curve
            .controls
            .iter()
            .map(|c| [c[0] / c[2], c[1] / c[2]])
            .collect()
    }

    /// The u/v at `t`.
    pub(in crate::numerical) fn point(&mut self, t: f64) -> [f64; 2] {
        let value = self.curve.point(t, &mut self.scratch).value;
        [value[0], value[1]]
    }
}

/// Source: the allowance `lift.rs` gives a domain border, sixteen roundings
/// of the largest coordinate, here for the corners of a box.
const BOX_ROUNDINGS: f64 = 16.0;

/// An edge's 3D curve in binary64, in a frame a whole body shares (rule T1),
/// for the distance readings of the tolerance rules.
pub(in crate::numerical) struct FramedCurve {
    curve: Curve64,
    scratch: Scratch,
    /// The parameter range, low end first.
    pub range: [f64; 2],
    /// The low end, every distinct knot strictly inside the range, and the
    /// high end.
    pub breaks: Vec<f64>,
    /// The polynomial degree of the pieces between breaks.
    pub degree: usize,
}

impl FramedCurve {
    /// Rounds `curve` between `parameters` in the frame shifted to `origin`.
    pub(in crate::numerical) fn new(
        curve: &Curve,
        parameters: [&Rational; 2],
        origin: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let curve = lift::framed_curve(curve, origin, dec)?;
        let ends = parameters.map(Rational::to_f64_nearest);
        let range = [ends[0].min(ends[1]), ends[0].max(ends[1])];
        let mut breaks = vec![range[0]];
        let degree = match &curve {
            Curve64::Line { .. } => 1,
            Curve64::Conic { .. } => 2,
            Curve64::Sampled(_) => unreachable!("a framed curve is rounded from a kernel curve"),
            Curve64::Spline(spline) => {
                for knot in &spline.knots {
                    if *knot > range[0] && *knot < range[1] && breaks.last() != Some(knot) {
                        breaks.push(*knot);
                    }
                }
                spline.degree
            }
        };
        breaks.push(range[1]);
        Ok(Self {
            curve,
            scratch: Scratch::default(),
            range,
            breaks,
            degree,
        })
    }

    /// The point at `t`.
    pub(in crate::numerical) fn point(&mut self, t: f64) -> [f64; 3] {
        self.curve.point(t, &mut self.scratch).value
    }

    /// A box that holds the curve over its range: the two ends of a line,
    /// a conic's centre widened by its two semi-axes, or the controls of a
    /// spline (Piegl and Tiller, "The NURBS Book", section 4.2; positive weights). It is
    /// widened by the binary64 rounding of its own corners.
    pub(in crate::numerical) fn bounding_box(&mut self) -> [[f64; 2]; 3] {
        let mut low = [f64::INFINITY; 3];
        let mut high = [f64::NEG_INFINITY; 3];
        let mut take = |point: [f64; 3], reach: [f64; 3]| {
            for axis in 0..3 {
                low[axis] = low[axis].min(point[axis] - reach[axis]);
                high[axis] = high[axis].max(point[axis] + reach[axis]);
            }
        };
        match &self.curve {
            Curve64::Line { .. } => {
                let ends = self.range;
                for t in ends {
                    let p = self.curve.point(t, &mut self.scratch).value;
                    take(p, [0.0; 3]);
                }
            }
            Curve64::Conic {
                center,
                major,
                minor,
            } => take(
                *center,
                std::array::from_fn(|axis| major[axis].abs() + minor[axis].abs()),
            ),
            Curve64::Spline(spline) => {
                for control in &spline.controls {
                    let w = control[3];
                    take([control[0] / w, control[1] / w, control[2] / w], [0.0; 3]);
                }
            }
            Curve64::Sampled(_) => unreachable!("a framed curve is rounded from a kernel curve"),
        }
        std::array::from_fn(|axis| {
            let slack = BOX_ROUNDINGS
                * numeric::UNIT_ROUNDOFF
                * low[axis].abs().max(high[axis].abs()).max(f64::MIN_POSITIVE);
            [low[axis] - slack, high[axis] + slack]
        })
    }
}

/// The binary64 evaluator of a 3D curve with no kernel form: the point at a
/// parameter, and a bound on the rounding error of that point.
pub type CurveEvaluator = Arc<dyn Fn(f64) -> ([f64; 3], f64) + Send + Sync>;

/// The binary64 evaluator of a surface with no kernel form: the point at a
/// u/v, and a bound on the rounding error of that point. It answers at every
/// u/v it is asked, continuing past its own borders as it sees fit.
pub type SurfaceEvaluator = Arc<dyn Fn([f64; 2]) -> ([f64; 3], f64) + Send + Sync>;

/// A 3D curve read through its own evaluator.
///
/// A reader of a foreign file holds carriers the kernel has no form for: an
/// offset curve, a curve on a swept surface, a curve in the file's own
/// parameter units. The checker reads such a curve at the points it chooses,
/// as it reads a kernel curve. What it needs of the shape is stated here, so
/// the point count still depends on the shape alone (decision Q22, guard 2).
///
/// The evaluator gives points in a frame of its owner's choice. The surface
/// of the same edge use is then a [`SampledSurface`] in the same frame: the
/// checker cannot shift an evaluator to a local frame of its own (rule T1).
#[derive(Clone)]
pub struct SampledCurve {
    /// The point at a parameter and its rounding bound.
    pub point: CurveEvaluator,
    /// Parameters at which the curve changes piece: its knots, the borders
    /// of its charts, the joints of a traced curve. In any order; those
    /// outside the edge are not used.
    pub breaks: Vec<f64>,
    /// The polynomial degree of its pieces. Zero for a curve that has none,
    /// as for a circle.
    pub degree: usize,
    /// The work one evaluation counts against the edge's allowance, on top
    /// of the checker's own unit. The checker takes it from the allowance
    /// before each evaluation and charges it to the operation, so the hard
    /// work cap of the edge covers it (section 4.2 step 8). Zero for an
    /// evaluator that costs what a kernel curve does.
    pub work: u64,
}

impl std::fmt::Debug for SampledCurve {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        out.debug_struct("SampledCurve")
            .field("breaks", &self.breaks)
            .field("degree", &self.degree)
            .field("work", &self.work)
            .finish_non_exhaustive()
    }
}

/// A surface read through its own evaluator. See [`SampledCurve`].
#[derive(Clone)]
pub struct SampledSurface {
    /// The point at a u/v and its rounding bound.
    pub point: SurfaceEvaluator,
    /// The values of u, then of v, at which the surface changes piece: its
    /// knot lines, the borders of its charts and its pole lines. A face copy
    /// is cut where it crosses one (section 4.2 step 3).
    pub lines: [Vec<f64>; 2],
    /// The period of u, then of v, where the surface closes on itself. Every
    /// whole number of periods is a seam, and a face copy is cut there too.
    pub periods: [Option<f64>; 2],
    /// The sum of its two polynomial degrees. Zero for a surface that has
    /// none, as for a cylinder.
    pub degree: usize,
    /// The work one evaluation counts against the edge's allowance, on top
    /// of the checker's own unit. See [`SampledCurve::work`].
    pub work: u64,
}

impl std::fmt::Debug for SampledSurface {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        out.debug_struct("SampledSurface")
            .field("lines", &self.lines)
            .field("periods", &self.periods)
            .field("degree", &self.degree)
            .field("work", &self.work)
            .finish_non_exhaustive()
    }
}

/// The 3D copy of an edge.
#[derive(Clone, Copy, Debug)]
pub enum GapCurve<'a> {
    /// A kernel curve: a line, a circle, an ellipse or a B-spline.
    Kernel(&'a Curve),
    /// A curve read through its own evaluator.
    Sampled(&'a SampledCurve),
}

/// The surface a face copy is lifted through.
#[derive(Clone, Copy, Debug)]
pub enum GapSurface<'a> {
    /// A kernel surface. A spline support is evaluated through its carrier.
    Kernel(&'a Surface),
    /// A spline carrier on its own.
    Carrier(&'a BSplineSurface),
    /// A surface read through its own evaluator.
    Sampled(&'a SampledSurface),
}

/// The copy of an edge in a face's parameters.
#[derive(Clone, Copy, Debug)]
pub enum FaceCopyRef<'a> {
    /// A retained u/v curve in the 3D curve's parameter.
    Native(&'a BSplineCurve),
    /// Straight in u/v between two end values, in traversal order.
    Linear(&'a [[Rational; 2]; 2]),
    /// The foot of each 3D point on a plane. The gap is then the distance of
    /// the 3D curve from the plane.
    Projected,
}

/// The two copies of one edge use.
#[derive(Clone, Copy, Debug)]
pub struct EdgeCopies<'a> {
    /// The 3D copy.
    pub curve: GapCurve<'a>,
    /// Edge parameters in traversal order.
    pub parameters: [&'a Rational; 2],
    /// The surface of the face.
    pub surface: GapSurface<'a>,
    /// The face copy.
    pub face_copy: FaceCopyRef<'a>,
}

/// What a builder knows about its own edge.
#[derive(Clone, Debug, Default)]
pub struct Hints {
    /// Parameters the builder wants looked at: joints of a traced curve that
    /// are not knots, and the place of its own largest residual if it has
    /// one. Each becomes a break, so it is sampled and bracketed.
    pub parameters: Vec<f64>,
}

/// Limits fixed before a measurement starts.
#[derive(Clone, Copy, Debug)]
pub struct Limits {
    /// Most work one edge may spend before a named refusal. A unit is one
    /// break or one halving of a crossing search while planning, one residual
    /// evaluation, one sample and mode pair of a tail transform, or one
    /// domain probe: everything the measurement does that grows with the
    /// edge.
    ///
    /// An edge use has one allowance. Where a use is read more than once,
    /// as the gap rule for imports reads it, each reading is given what the
    /// readings before it left, so the cap holds for all of them together.
    pub work: u64,
}

/// Source: judgement; no source gives a per-edge cap. One million units is
/// five percent of the 20M operation target (ROADMAP, rules that stay fixed
/// between checkpoints).
pub const EDGE_WORK: u64 = crate::numerical::DEFAULT_WORK / 20;

impl Default for Limits {
    fn default() -> Self {
        Self { work: EDGE_WORK }
    }
}

/// How a result was reached.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Method {
    /// A straight residual, read at its two ends.
    StraightResidual,
    /// A conic residual, read at its ends and its stationary points.
    ConicResidual,
    /// The two ends alone, by request.
    EndsOnly,
    /// Cut, sampled, refined and polished.
    Sampled,
}

/// Why a piece was sampled again.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Reason {
    /// Its Chebyshev coefficient tail was large against the largest gap seen.
    LoudTail,
}

/// One doubling of the points of a piece after its first sampling.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Refinement {
    /// The piece whose points were doubled.
    pub piece: [f64; 2],
    /// Why.
    pub reason: Reason,
    /// The interior sample points of its new grid.
    pub points: usize,
}

/// One piece of the first sample plan.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PlannedPiece {
    /// The piece in edge-parameter order.
    pub piece: [f64; 2],
    /// Interior sample points. Both ends are sampled as well.
    pub points: usize,
}

/// The largest gap found on one edge use.
#[derive(Clone, Debug)]
pub struct GapMeasure {
    /// The largest same-parameter distance found, in millimetres.
    pub raw_mm: f64,
    /// The measurement error of rule T4: an estimate of the binary64
    /// rounding error of that value, plus an estimate of how far the true
    /// maximum may lie above it because of what the samples did not resolve.
    /// The second part is zero for a straight residual.
    pub error_mm: f64,
    /// The noise floor of step 7: the binary64 rounding part of
    /// `error_mm`. Only this part can be removed by evaluating at a higher
    /// precision.
    pub noise_mm: f64,
    /// The edge parameter where it was found.
    pub parameter: f64,
    /// The first parameter in traversal order at which the face copy is
    /// outside a spline surface's domain. The gap there is measured against
    /// the end span's polynomial continuation.
    pub domain_exit: Option<f64>,
    /// How the result was reached.
    pub method: Method,
    /// The first sample plan, which depends only on the shape of the edge.
    pub plan: Vec<PlannedPiece>,
    /// Every later doubling of a piece's points, with its reason.
    pub refinements: Vec<Refinement>,
    /// Residual evaluations spent.
    pub evaluations: u64,
    /// Work counted against the edge's cap: the planning, the evaluations,
    /// the tail transforms and the domain probes.
    pub work: u64,
}

/// Source: rule T4 of the tolerance rules (Checkpoint 2 decisions, section
/// 4.3): the stored tolerance is 1.1 times the measured gap plus the
/// measurement error.
pub const MARGIN: f64 = 1.1;

impl GapMeasure {
    /// The reported gap: the raw value plus its measurement error.
    pub fn reported_mm(&self) -> f64 {
        (self.raw_mm + self.error_mm).next_up()
    }
    /// The tolerance rule T4 stores for a fresh measurement, rounded up.
    pub fn stored_mm(&self) -> f64 {
        (MARGIN * self.reported_mm()).next_up()
    }
}

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("numerical.gap-unsupported", detail)
}

/// The checker's refusal when an edge would spend more than its allowance.
/// Source: decision Q6, section 4.2 step 8: a hard per-edge work cap refuses
/// by name.
pub(super) const WORK_CAP: &str = "numerical.gap-work-cap";

fn work_cap(detail: &str) -> Refusal {
    Refusal::new(WORK_CAP, detail)
}

/// Source: decision Q6, section 4.2 step 4, after Chebfun's seventeen-point
/// first grid (Aurentz and Trefethen, "Chopping a Chebyshev series", ACM
/// TOMS 43, 2017).
const LEAST_POINTS: usize = 17;
/// Source: OpenCASCADE `GeomLib_CheckCurveOnSurface::FillSubIntervals`:
/// three samples per unit of composed degree, because a polynomial of degree
/// d has at most d - 1 interior extrema.
const POINTS_PER_DEGREE: usize = 3;

/// The degrees that decide a piece's sample count. Nothing else does.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct PieceShape {
    /// Degree of the 3D copy; zero for a circle or an ellipse.
    pub curve: usize,
    /// Degree of the face copy.
    pub face: usize,
    /// Sum of the surface's two degrees; one for a plane, zero for a curved
    /// analytic surface.
    pub surface: usize,
}

impl PieceShape {
    /// The composed degree of the residual: the 3D copy's degree, or the face
    /// copy's degree times the surface's degree sum (the OpenCASCADE rule).
    pub fn composed_degree(&self) -> usize {
        self.curve.max(self.face * self.surface)
    }
    /// Interior sample points for a piece of this shape:
    /// `max(17, 3 x composed degree)`, made even so that no point is a
    /// rational share of the piece (see `numeric::chebyshev_shares`). An odd
    /// grid would have the middle of the piece, which is one of the eighths
    /// a builder checks at.
    pub fn points(&self) -> usize {
        let count = LEAST_POINTS.max(POINTS_PER_DEGREE * self.composed_degree());
        count + count % 2
    }
}

fn shape(lifted: &Lifted) -> PieceShape {
    let curve = match &lifted.curve {
        Curve64::Line { .. } => 1,
        Curve64::Conic { .. } => 0,
        Curve64::Spline(spline) => spline.degree,
        Curve64::Sampled(sampled) => sampled.degree,
    };
    PieceShape {
        curve,
        face: match &lifted.face {
            Face64::Native(spline) => spline.degree,
            Face64::Linear { .. } => 1,
            // The foot of the 3D point has the 3D copy's degree.
            Face64::Projected => curve,
        },
        surface: match &lifted.surface {
            Surface64::Plane { .. } => 1,
            Surface64::Revolved { .. } => 0,
            Surface64::Spline(net) => net.degrees[0] + net.degrees[1],
            Surface64::Sampled(sampled) => sampled.degree,
        },
    }
}

/// Source: IEEE 754 binary64 has 53 significant bits. Two breaks closer than
/// 2^-45 of the edge are within 256 units in the last place of each other at
/// the edge's own scale and leave no room for an interior sample.
const NARROWEST_PIECE: f64 = 1.0 / ((1_u64 << 45) as f64);

/// The seams of a periodic parameter within `range`: every whole number of
/// `turn`. One spare turn covers the rounding of the two quotients.
///
/// The number of seams is fixed before any is listed. A face copy that winds
/// around its surface more than `most` times is refused, not walked.
fn seams(turn: f64, range: [f64; 2], most: u64, out: &mut Vec<f64>) -> Result<(), Refusal> {
    let first = (range[0] / turn).ceil();
    let turns = (range[1] / turn).floor() - first + 2.0;
    if !turns.is_finite() || turns > most as f64 {
        return Err(work_cap(
            "a face copy winds around its surface more often than its edge may be cut",
        ));
    }
    for index in 0..turns.max(0.0) as u64 {
        let seam = (first + index as f64) * turn;
        if seam <= range[1] {
            out.push(seam);
        }
    }
    Ok(())
}

/// The values of one surface parameter, within `range`, at which the lifted
/// curve is cut: knot lines and domain borders of a spline surface, the seam
/// of a periodic chart, and the pole lines of a sphere or a cone. A sampled
/// surface states its own lines and periods.
fn cut_lines(
    surface: &Surface64,
    axis: usize,
    range: [f64; 2],
    most: u64,
    out: &mut Vec<f64>,
) -> Result<(), Refusal> {
    let within = |value: f64| value >= range[0] && value <= range[1];
    match surface {
        Surface64::Plane { .. } => (),
        Surface64::Spline(net) => {
            out.extend(net.knots[axis].iter().copied().filter(|knot| within(*knot)));
        }
        Surface64::Revolved { kind, .. } => {
            let turn = 2.0 * std::f64::consts::PI;
            if axis == 0 {
                // The seam, u = 0 modulo one turn.
                seams(turn, range, most, out)?;
            } else {
                let poles: &[f64] = match kind {
                    Revolved::Sphere => {
                        &[-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2]
                    }
                    Revolved::Cone => &[0.0],
                    _ => &[],
                };
                out.extend(poles.iter().copied().filter(|pole| within(*pole)));
            }
        }
        Surface64::Sampled(sampled) => {
            out.extend(
                sampled.lines[axis]
                    .iter()
                    .copied()
                    .filter(|line| within(*line)),
            );
            if let Some(turn) = sampled.periods[axis] {
                if !(turn.is_finite() && turn > 0.0) {
                    return Err(unsupported("a sampled surface has no finite period"));
                }
                seams(turn, range, most, out)?;
            }
        }
    }
    out.dedup();
    Ok(())
}

/// Parameters in `(from, to)` where the face copy crosses a cut line of the
/// surface. The face copy is one polynomial or rational piece on the interval.
///
/// Each cut line and each halving of its root search is one unit of the
/// edge's allowance, taken before that work is done. The search stops with
/// the named refusal where the allowance runs out.
fn crossings(
    search: &mut Search<'_>,
    from: f64,
    to: f64,
    out: &mut Vec<f64>,
) -> Result<(), Refusal> {
    let lifted = search.lifted;
    // Homogeneous Bezier coefficients (u w, v w, w) of the face copy.
    let piece: Vec<[f64; 3]> = match &lifted.face {
        Face64::Projected => return Ok(()),
        Face64::Linear { ends, .. } => {
            let [a, b] = lifted.parameters;
            [from, to]
                .map(|t| {
                    let share = (t - a) / (b - a);
                    [
                        ends[0][0] + share * (ends[1][0] - ends[0][0]),
                        ends[0][1] + share * (ends[1][1] - ends[0][1]),
                        1.0,
                    ]
                })
                .into()
        }
        Face64::Native(spline) => {
            let middle = from + (to - from) / 2.0;
            let span =
                numeric::find_span(&spline.knots, spline.degree, spline.controls.len(), middle);
            numeric::bezier_piece(
                &spline.knots,
                spline.degree,
                &spline.controls,
                span,
                from,
                to,
            )
        }
    };
    let mut lines = Vec::new();
    let mut shares = Vec::new();
    for axis in 0..2 {
        // The convex hull of the Cartesian coefficients bounds the piece.
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for control in &piece {
            lo = lo.min(control[axis] / control[2]);
            hi = hi.max(control[axis] / control[2]);
        }
        if !(lo.is_finite() && hi.is_finite()) {
            return Err(unsupported(
                "a face copy has no finite u/v range on one of its pieces",
            ));
        }
        lines.clear();
        cut_lines(&lifted.surface, axis, [lo, hi], search.left(), &mut lines)?;
        for line in &lines {
            // The visit to this line, then its halvings one by one.
            search.reserve(1)?;
            let coefficients: Vec<f64> = piece.iter().map(|c| c[axis] - line * c[2]).collect();
            let size = piece
                .iter()
                .fold(line.abs(), |sum, c| sum.max((c[axis] / c[2]).abs()));
            shares.clear();
            let before = search.left();
            let mut left = before;
            let finished = numeric::sign_changes(
                &coefficients,
                numeric::UNIT_ROUNDOFF * numeric::roundings(piece.len() - 1) * size,
                &mut shares,
                &mut left,
            );
            // The halvings are recorded whether or not the search finished.
            search.reserve(before - left)?;
            if finished.is_none() {
                return Err(work_cap(
                    "an edge spent its work allowance before its breaks were found",
                ));
            }
            out.extend(shares.iter().map(|share| from + share * (to - from)));
        }
    }
    Ok(())
}

/// Every break of the edge in ascending order, both ends included.
///
/// Planning draws on the same allowance as sampling. Each break is one unit,
/// and so is each cut line a piece is tested against and each halving of
/// that test. Every unit is taken from the edge's allowance before its work
/// is done, so planning stops where the allowance runs out and the search
/// that follows starts with what planning left.
fn breaks(search: &mut Search<'_>, hints: &Hints, dec: &mut Decider) -> Result<Vec<f64>, Refusal> {
    let lifted = search.lifted;
    let [lo, hi] = lifted.range();
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return Err(unsupported("an edge has no finite parameter range"));
    }
    let inside = |t: &f64| *t > lo && *t < hi;
    let mut cuts = vec![lo, hi];
    match &lifted.curve {
        Curve64::Spline(spline) => cuts.extend(spline.knots.iter().copied().filter(inside)),
        Curve64::Sampled(sampled) => cuts.extend(sampled.breaks.iter().copied().filter(inside)),
        Curve64::Line { .. } | Curve64::Conic { .. } => (),
    }
    if let Face64::Native(spline) = &lifted.face {
        cuts.extend(spline.knots.iter().copied().filter(inside));
    }
    cuts.extend(hints.parameters.iter().copied().filter(inside));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    search.reserve(cuts.len() as u64)?;
    search.settle(dec)?;
    // Between these cuts the face copy is one polynomial or rational piece.
    let mut crossed = Vec::new();
    for pair in cuts.windows(2) {
        crossings(search, pair[0], pair[1], &mut crossed)?;
        search.settle(dec)?;
    }
    cuts.extend(crossed.into_iter().filter(inside));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    Ok(cuts)
}

/// The first sample plan of an edge. It is a function of the edge's shape:
/// its degrees, its knots and the places where the face copy crosses the
/// surface's knot lines. A scaled or shifted copy of the edge gets the same
/// plan. A straight residual is read at its two ends and has no plan.
pub fn plan(
    copies: &EdgeCopies<'_>,
    hints: &Hints,
    limits: Limits,
    dec: &mut Decider,
) -> Result<Vec<PlannedPiece>, Refusal> {
    let lifted = lift::lift(copies, dec)?;
    if straight_residual(&lifted) || conic_residual(&lifted) {
        return Ok(Vec::new());
    }
    let mut search = Search {
        lifted: &lifted,
        scratch: Scratch::default(),
        limit: limits.work,
        evaluations: 0,
        other: 0,
        charged: 0,
    };
    let breaks = breaks(&mut search, hints, dec);
    if breaks.is_err() {
        search.settle(dec).ok();
    }
    Ok(planned(&lifted, &breaks?))
}

fn planned(lifted: &Lifted, breaks: &[f64]) -> Vec<PlannedPiece> {
    // Every piece gets the full count, however narrow it is.
    let points = shape(lifted).points();
    breaks
        .windows(2)
        .map(|pair| PlannedPiece {
            piece: [pair[0], pair[1]],
            points,
        })
        .collect()
}

/// One sample: its parameter, the residual, the squared gap and the size of
/// the residual's rounding error.
#[derive(Clone, Copy)]
struct Sample {
    t: f64,
    residual: [f64; 3],
    squared: f64,
    error: f64,
}

struct Search<'a> {
    lifted: &'a Lifted,
    scratch: Scratch,
    limit: u64,
    evaluations: u64,
    /// Work other than residual evaluations: transforms and domain probes.
    other: u64,
    charged: u64,
}

impl Search<'_> {
    fn spent(&self) -> u64 {
        self.evaluations + self.other
    }

    /// What is left of the edge's allowance.
    fn left(&self) -> u64 {
        self.limit.saturating_sub(self.spent())
    }

    fn sample(&mut self, t: f64) -> Result<Sample, Refusal> {
        // The evaluation's own unit and what a sampled carrier states for
        // it, taken before the evaluation is made.
        let stated = self.lifted.sample_work() - 1;
        if self.spent() >= self.limit || stated > self.limit - self.spent() - 1 {
            return Err(work_cap(
                "an edge spent its work allowance before its gap was resolved",
            ));
        }
        self.evaluations += 1;
        self.other += stated;
        let lift::Bounded { value, error } = self.lifted.residual(t, &mut self.scratch);
        Ok(Sample {
            t,
            residual: value,
            squared: value[0] * value[0] + value[1] * value[1] + value[2] * value[2],
            error,
        })
    }

    /// Counts work other than a residual evaluation against the cap, before
    /// that work is done.
    fn reserve(&mut self, units: u64) -> Result<(), Refusal> {
        if units > self.left() {
            return Err(work_cap(
                "an edge spent its work allowance before its gap was resolved",
            ));
        }
        self.other += units;
        Ok(())
    }

    /// Charges one work unit for every unit spent since the last call, as
    /// the kernel's other binary64 samplers charge. It is called once per
    /// piece, doubling and polish, so the operation allowance is checked that
    /// often.
    fn settle(&mut self, dec: &mut Decider) -> Result<(), Refusal> {
        let due = self.spent() - self.charged;
        if due == 0 {
            return Ok(());
        }
        self.charged = self.spent();
        dec.charge_work(due)
    }
}

/// Source: Trefethen, "Approximation Theory and Approximation Practice",
/// 2013, chapter 4: a Chebyshev interpolant is off by at most twice the sum
/// of the magnitudes of the coefficients it does not carry.
const INTERPOLATION_FACTOR: f64 = 2.0;

/// Source: the triangle inequality. If the interpolant is within `e` of the
/// residual, the residual's true maximum is at most `e` above the
/// interpolant's, and the polished maximum is at most `e` below it, so the
/// search can be short by two interpolation errors.
const SEARCH_FACTOR: f64 = 2.0;

/// What one piece's samples leave unresolved.
struct Unresolved {
    /// An estimate of how far the true maximum on the piece may lie above
    /// the largest value the search can reach there.
    search: f64,
    /// The part of that estimate that rounding noise alone can produce. Each
    /// coefficient is at most twice the largest sample, so the samples'
    /// rounding errors can add twice their size per tail coefficient.
    noise: f64,
}

/// Reads the tail of one grid. Its work is counted against the cap first.
fn unresolved(search: &mut Search<'_>, grid: &[Sample]) -> Result<Unresolved, Refusal> {
    search.reserve(numeric::chebyshev_tail_pairs(grid.len()))?;
    let residuals: Vec<[f64; 3]> = grid.iter().map(|sample| sample.residual).collect();
    let error = grid
        .iter()
        .fold(0.0_f64, |sum, sample| sum.max(sample.error));
    let scale = SEARCH_FACTOR * INTERPOLATION_FACTOR;
    Ok(Unresolved {
        search: scale * numeric::chebyshev_tail(&residuals),
        noise: scale * grid.len() as f64 * error,
    })
}

impl Unresolved {
    /// Whether the piece leaves more unresolved than the storage margin
    /// covers: the search uncertainty against the margin's share of the
    /// largest gap seen so far, with rounding noise set aside.
    fn loud(&self, largest: f64) -> bool {
        self.search > (MARGIN - 1.0) * largest + self.noise
    }
}

/// Measures the largest gap between the two copies of one edge use.
///
/// Refuses with `numerical.gap-unsupported` for a curve or surface without a
/// parameter evaluator, and with `numerical.gap-work-cap` when the edge would
/// spend more than its work allowance. One work unit is charged per
/// evaluation, as the kernel's other binary64 samplers charge, and one per
/// sample and mode pair of a tail transform and per domain probe, in one
/// charge per piece.
pub fn measure(
    copies: &EdgeCopies<'_>,
    hints: &Hints,
    limits: Limits,
    dec: &mut Decider,
) -> Result<GapMeasure, Refusal> {
    under_allowance(copies, limits, dec, |search, dec| {
        search_edge(search, hints, dec)
    })
}

/// Runs one search of an edge under its allowance.
fn under_allowance<T>(
    copies: &EdgeCopies<'_>,
    limits: Limits,
    dec: &mut Decider,
    run: impl FnOnce(&mut Search<'_>, &mut Decider) -> Result<T, Refusal>,
) -> Result<T, Refusal> {
    let lifted = lift::lift(copies, dec)?;
    // The closing probes of the rounding allowance count inside the same cap,
    // so they are set aside before the first sample.
    let Some(limit) = limits.work.checked_sub(lifted.closing_work()) else {
        return Err(work_cap(
            "the work allowance of an edge is smaller than its closing probes",
        ));
    };
    let mut search = Search {
        lifted: &lifted,
        scratch: Scratch::default(),
        limit,
        evaluations: 0,
        other: 0,
        charged: 0,
    };
    let found = run(&mut search, dec);
    if let Err(refusal) = &found {
        // Work done before a refusal is still work: it is charged to the
        // operation before the refusal is returned. If the operation's own
        // allowance is gone as well, the edge's refusal is still the answer.
        let (spent, limit) = (search.spent(), search.limit);
        search.settle(dec).ok();
        if refusal.code == WORK_CAP {
            return Err(work_cap(&format!(
                "{}; {spent} of {limit} units were spent",
                refusal.detail
            )));
        }
    }
    found
}

/// The search of one lifted edge under its allowance.
fn search_edge(
    search: &mut Search<'_>,
    hints: &Hints,
    dec: &mut Decider,
) -> Result<GapMeasure, Refusal> {
    let lifted = search.lifted;
    let [first, last] = lifted.parameters;
    if straight_residual(lifted) {
        let ends = [search.sample(first)?, search.sample(last)?];
        let best = if ends[1].squared > ends[0].squared {
            ends[1]
        } else {
            ends[0]
        };
        // Convexity makes the end values the maximum: nothing is unresolved.
        return finish(
            search,
            best,
            (0.0, None, Method::StraightResidual),
            (Vec::new(), Vec::new()),
            dec,
        );
    }
    if conic_residual(lifted) {
        let best = conic_maximum(search, dec)?;
        // Every stationary point was read: nothing is unresolved.
        return finish(
            search,
            best,
            (0.0, None, Method::ConicResidual),
            (Vec::new(), Vec::new()),
            dec,
        );
    }
    let breaks = breaks(search, hints, dec)?;
    let plan = planned(lifted, &breaks);
    // Refuse before the first sample when the plan alone is more than
    // planning left of the allowance: its samples, the first tail transform
    // of each piece and, on a spline surface, one domain probe per piece.
    let probes = u64::from(matches!(lifted.surface, Surface64::Spline(_)));
    // Each sample counts its own unit and what a sampled carrier states.
    let each = lifted.sample_work();
    let planned_work = plan.iter().fold(each, |sum, piece| {
        sum.saturating_add((piece.points as u64 + 1).saturating_mul(each))
            .saturating_add(numeric::chebyshev_tail_pairs(piece.points) + probes)
    });
    if planned_work > search.left() {
        return Err(work_cap(
            "the first sample plan of an edge is larger than its work allowance",
        ));
    }

    // First pass: every piece at its planned points, then its far end.
    // Neighbours share an end.
    let mut pieces: Vec<[Vec<Sample>; 2]> = Vec::with_capacity(plan.len());
    let mut largest = 0.0_f64;
    let mut ends = vec![search.sample(breaks[0])?];
    for piece in &plan {
        let [from, to] = piece.piece;
        let mut grid = Vec::with_capacity(piece.points);
        for share in numeric::chebyshev_shares(piece.points) {
            grid.push(search.sample(from + share * (to - from))?);
        }
        ends.push(search.sample(to)?);
        search.settle(dec)?;
        for sample in grid.iter().chain(ends.last()) {
            largest = largest.max(sample.squared.sqrt());
        }
        // The grid whose tail is read, and the samples of earlier grids.
        pieces.push([grid, Vec::new()]);
    }
    largest = largest.max(ends[0].squared.sqrt());

    // Second pass: double the points of every piece whose tail is loud.
    let mut refinements = Vec::new();
    let mut search_error = 0.0_f64;
    let mut samples: Vec<Sample> = Vec::new();
    for (index, [mut grid, mut earlier]) in pieces.into_iter().enumerate() {
        let (from, to) = (ends[index].t, ends[index + 1].t);
        loop {
            let left = unresolved(search, &grid)?;
            let fresh: Vec<f64> = numeric::chebyshev_shares(2 * grid.len())
                .into_iter()
                .map(|share| from + share * (to - from))
                .collect();
            // Once the points of a finer grid no longer fall strictly in
            // order between the ends, binary64 has no room for them.
            let room = fresh[0] > from
                && fresh[fresh.len() - 1] < to
                && fresh.windows(2).all(|pair| pair[0] < pair[1]);
            if !left.loud(largest) || !room {
                search_error = search_error.max(left.search);
                break;
            }
            refinements.push(Refinement {
                piece: [from, to],
                reason: Reason::LoudTail,
                points: fresh.len(),
            });
            earlier.append(&mut grid);
            for t in fresh {
                let sample = search.sample(t)?;
                largest = largest.max(sample.squared.sqrt());
                grid.push(sample);
            }
            search.settle(dec)?;
        }
        search.settle(dec)?;
        if index == 0 {
            samples.push(ends[0]);
        }
        let start = samples.len();
        samples.append(&mut earlier);
        samples.append(&mut grid);
        samples[start..].sort_by(|a, b| a.t.total_cmp(&b.t));
        samples.push(ends[index + 1]);
    }

    // Third pass: polish every local maximum of the samples.
    let best = polish(search, &samples, |search, t, _| search.sample(t), dec)?;

    // The first parameter in traversal order outside the surface's domain.
    let mut exit = None;
    if matches!(lifted.surface, Surface64::Spline(_)) {
        let forward = first <= last;
        let mut order: Vec<&[f64]> = breaks.windows(2).collect();
        if !forward {
            order.reverse();
        }
        for pair in order {
            let middle = pair[0] + (pair[1] - pair[0]) / 2.0;
            search.reserve(1)?;
            if lifted.outside(middle, &mut search.scratch) {
                exit = Some(if forward { pair[0] } else { pair[1] });
                break;
            }
        }
    }
    finish(
        search,
        best,
        (search_error, exit, Method::Sampled),
        (plan, refinements),
        dec,
    )
}

/// A sampled value that step 6 can bracket and polish.
trait Peak: Copy {
    /// The parameter it was read at.
    fn at(&self) -> f64;
    /// The squared gap there.
    fn squared(&self) -> f64;
    /// The size of its rounding error.
    fn error(&self) -> f64;
}

impl Peak for Sample {
    fn at(&self) -> f64 {
        self.t
    }
    fn squared(&self) -> f64 {
        self.squared
    }
    fn error(&self) -> f64 {
        self.error
    }
}

/// Step 6 of section 4.2: polishes every local maximum of `samples`, which
/// ascend in parameter, with Brent's method on the true gap, and returns the
/// largest value read. `read` evaluates the gap at a parameter inside the
/// bracket of the sample it is given.
fn polish<S: Peak>(
    search: &mut Search<'_>,
    samples: &[S],
    read: impl Fn(&mut Search<'_>, f64, &S) -> Result<S, Refusal>,
    dec: &mut Decider,
) -> Result<S, Refusal> {
    let [lo, hi] = search.lifted.range();
    let mut best = samples[0];
    for sample in samples {
        if sample.squared() > best.squared() {
            best = *sample;
        }
    }
    let mut i = 0;
    while i < samples.len() {
        // A run of equal values is one candidate.
        let mut j = i;
        while j + 1 < samples.len()
            && samples[j + 1].squared().to_bits() == samples[i].squared().to_bits()
        {
            j += 1;
        }
        let left = i.checked_sub(1);
        let right = (j + 1 < samples.len()).then_some(j + 1);
        let peak = left.is_none_or(|l| samples[l].squared() < samples[i].squared())
            && right.is_none_or(|r| samples[r].squared() < samples[i].squared());
        // A value inside its own rounding error is noise, not a maximum.
        let above_noise = samples[i].squared().sqrt() > samples[i].error();
        let (from, to) = (
            samples[left.unwrap_or(i)].at(),
            samples[right.unwrap_or(j)].at(),
        );
        let resolvable = to - from > 2.0 * NARROWEST_PIECE * (hi - lo);
        if peak && above_noise && resolvable {
            let mut failure = None;
            let mut polished = samples[i];
            let (_, _, _) = numeric::brent_maximum(&mut |share| match read(
                search,
                from + share * (to - from),
                &samples[i],
            ) {
                Ok(sample) => {
                    if sample.squared() > polished.squared() {
                        polished = sample;
                    }
                    sample.squared()
                }
                Err(refusal) => {
                    failure.get_or_insert(refusal);
                    f64::NEG_INFINITY
                }
            });
            if let Some(refusal) = failure {
                return Err(refusal);
            }
            search.settle(dec)?;
            if polished.squared() > best.squared() {
                best = polished;
            }
        }
        i = j + 1;
    }
    Ok(best)
}

/// The gap at the two ends of an edge alone, for a face copy that is
/// declared only there.
///
/// `limits` is what the edge use has left of its one allowance when another
/// measurement of the same use came first. The two end samples and the
/// closing probes draw on it, and it refuses by name like [`measure`].
pub fn measure_ends(
    copies: &EdgeCopies<'_>,
    limits: Limits,
    dec: &mut Decider,
) -> Result<GapMeasure, Refusal> {
    under_allowance(copies, limits, dec, |search, dec| {
        let [first, last] = search.lifted.parameters;
        let ends = [search.sample(first)?, search.sample(last)?];
        let best = if ends[1].squared > ends[0].squared {
            ends[1]
        } else {
            ends[0]
        };
        finish(
            search,
            best,
            (0.0, None, Method::EndsOnly),
            (Vec::new(), Vec::new()),
            dec,
        )
    })
}

fn finish(
    search: &mut Search<'_>,
    best: Sample,
    (search_error, domain_exit, method): (f64, Option<f64>, Method),
    (plan, refinements): (Vec<PlannedPiece>, Vec<Refinement>),
    dec: &mut Decider,
) -> Result<GapMeasure, Refusal> {
    // Half the narrowest piece next to the reported point, or half the edge.
    let [lo, hi] = search.lifted.range();
    let reach = plan
        .iter()
        .filter(|piece| piece.piece[0] <= best.t && best.t <= piece.piece[1])
        .map(|piece| piece.piece[1] - piece.piece[0])
        .fold(hi - lo, f64::min)
        / 2.0;
    // These probes were set aside when the search started, so the total
    // stays inside the edge's work allowance.
    let (allowance, evaluations) = search.lifted.allowance(best.t, reach, &mut search.scratch);
    debug_assert!(evaluations <= lift::ALLOWANCE_EVALUATIONS);
    search.evaluations += evaluations;
    // What the sampled carriers state for the points those probes read.
    search.other += search.lifted.closing_work() - lift::ALLOWANCE_EVALUATIONS;
    search.settle(dec)?;
    Ok(GapMeasure {
        raw_mm: best.squared.sqrt(),
        error_mm: allowance + search_error,
        noise_mm: allowance,
        parameter: best.t,
        domain_exit,
        method,
        plan,
        refinements,
        evaluations: search.evaluations,
        work: search.spent(),
    })
}

/// Whether the residual is a conic in the parameter,
/// `O + A cos t + B sin t`: a circle or an ellipse against a lift that is one
/// too, at the same angular rate. That holds for the foot of the curve on a
/// plane, and on a surface of revolution for a straight copy that keeps its
/// v and moves its u at exactly the rate of the edge parameter. OpenCASCADE's
/// `ElementaryMaximum` takes the same pairs: both lifted curves are circles
/// or ellipses.
fn conic_residual(lifted: &Lifted) -> bool {
    if !matches!(lifted.curve, Curve64::Conic { .. }) {
        return false;
    }
    match (&lifted.surface, &lifted.face) {
        (Surface64::Plane { .. }, Face64::Projected) => true,
        (
            Surface64::Revolved { .. },
            Face64::Linear {
                turns_with_edge, ..
            },
        ) => *turns_with_edge,
        _ => false,
    }
}

/// The largest sample of a conic residual `r = O + A cos t + B sin t`.
///
/// `O`, `A` and `B` come from three evaluations a quarter turn apart. The
/// squared length is stationary where `r . r' = 0`. On a quarter turn from
/// angle `a`, with `z = tan((t - a) / 2)` in `[0, 1]`, that is the quartic
///
/// ```text
/// (p + m) + 2 (k - q) z - 6 m z^2 - 2 (q + k) z^3 + (m - p) z^4 = 0
/// ```
///
/// with `p = O . B`, `q = O . A`, `m = A . B` and `k = |B|^2 - |A|^2` for the
/// `A` and `B` of that quarter. Its roots are isolated in Bernstein form by
/// halving. The residual is then evaluated at each root, at each quarter
/// border and at the two ends, and the largest value is returned. One turn
/// holds every stationary value of a longer edge.
/// Source: OpenCASCADE `ElementaryMaximum` (the same stationary equation,
/// solved there by `math_TrigonometricFunctionRoots`).
fn conic_maximum(search: &mut Search<'_>, dec: &mut Decider) -> Result<Sample, Refusal> {
    let [lo, hi] = search.lifted.range();
    let quarter = std::f64::consts::FRAC_PI_2;
    let basis = [
        search.sample(lo)?,
        search.sample(lo + quarter)?,
        search.sample(lo + 2.0 * quarter)?,
    ];
    let offset: [f64; 3] =
        std::array::from_fn(|k| (basis[0].residual[k] + basis[2].residual[k]) / 2.0);
    let cosine: [f64; 3] =
        std::array::from_fn(|k| (basis[0].residual[k] - basis[2].residual[k]) / 2.0);
    let sine: [f64; 3] = std::array::from_fn(|k| basis[1].residual[k] - offset[k]);
    let dot = |a: &[f64; 3], b: &[f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let negated = |a: &[f64; 3]| a.map(|x| -x);
    let mut best = basis[0];
    let mut read = |search: &mut Search<'_>, t: f64| -> Result<(), Refusal> {
        if t > lo && t <= hi {
            let sample = search.sample(t)?;
            if sample.squared > best.squared {
                best = sample;
            }
        }
        Ok(())
    };
    read(search, hi)?;
    let quarters = ((hi - lo) / quarter).ceil().min(numeric::QUADRANTS) as usize;
    let mut shares = Vec::new();
    for index in 0..quarters {
        let from = lo + index as f64 * quarter;
        // The cosine and sine vectors, turned by whole quarters.
        let (a, b) = match index {
            0 => (cosine, sine),
            1 => (sine, negated(&cosine)),
            2 => (negated(&cosine), negated(&sine)),
            _ => (negated(&sine), cosine),
        };
        let (p, q, m) = (dot(&offset, &b), dot(&offset, &a), dot(&a, &b));
        let k = dot(&b, &b) - dot(&a, &a);
        let power = [
            p + m,
            2.0 * (k - q),
            -2.0 * (m + m + m),
            -2.0 * (q + k),
            m - p,
        ];
        let size = power.iter().fold(0.0_f64, |sum, c| sum.max(c.abs()));
        shares.clear();
        // The quartic of this quarter, then its halvings one by one, each
        // taken from the edge's allowance before it is done.
        search.reserve(1)?;
        let before = search.left();
        let mut left = before;
        let finished = numeric::sign_changes(
            &numeric::bernstein_from_power(&power),
            numeric::UNIT_ROUNDOFF * numeric::roundings(power.len() - 1) * size,
            &mut shares,
            &mut left,
        );
        // The halvings are recorded whether or not the search finished.
        search.reserve(before - left)?;
        if finished.is_none() {
            return Err(work_cap(
                "an edge spent its work allowance before its gap was resolved",
            ));
        }
        read(search, from)?;
        for z in &shares {
            read(search, from + 2.0 * numeric::arctangent(*z))?;
        }
    }
    search.settle(dec)?;
    Ok(best)
}

/// Whether the residual is a straight function of the parameter: a line
/// against a face copy whose lift is also straight. That holds on a plane for
/// a straight or projected copy, and along a generator of a cylinder or cone,
/// where the copy keeps its angle. The norm of a straight function is convex,
/// so its maximum is at an end.
fn straight_residual(lifted: &Lifted) -> bool {
    if !matches!(lifted.curve, Curve64::Line { .. }) {
        return false;
    }
    match (&lifted.surface, &lifted.face) {
        (Surface64::Plane { .. }, Face64::Linear { .. } | Face64::Projected) => true,
        (
            Surface64::Revolved {
                kind: Revolved::Cylinder | Revolved::EllipticCylinder | Revolved::Cone,
                ..
            },
            Face64::Linear { ends, .. },
        ) => (ends[1][0] - ends[0][0]).abs() <= 0.0,
        _ => false,
    }
}
