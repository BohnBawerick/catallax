//! Numerical finite trace/patch meetings; never exact incidence or an empty-set proof.
//!
//! OCCT HInter accepts polygonal curve proposals. Regular retained traces with
//! second derivatives are adapted to the existing SMLib GlobalCurveIntersect
//! port using derivative-bounded polylines. Its finite subdivision rejects boxes
//! and uses guarded local solves. Reported hits are checked on both originals.

use crate::{CurveTrace, TracePiece};
use catallax_numbers::{
    Ball, Decider, NumericalCheck, NumericalEvidence, NumericalSample, Rational, Real, Refusal,
    Sign,
};
use catallax_shapes::RationalPatch;
use catallax_trace::{
    curve_surface_global::global_curve_intersect_with_work,
    extent::{Extent1, Extent2},
    solution::SolutionType,
    spline::{BSplineCurve, BSplineSurface, ControlPoint},
    surface_work::WorkReport,
    vector::Vector2,
};

const BITS: u32 = 128;

/// One numerical meeting, retaining the original finite source parameters.
#[derive(Clone, Debug)]
pub struct NumericalCurveSurfaceHit {
    /// Index in the original trace's ordered pieces.
    pub piece: usize,
    /// Parameter in that piece's closed unit interval.
    pub parameter: Rational,
    /// Parameters of the original rational patch.
    pub uv: [Rational; 2],
    /// Original trace position enclosure, not an exact meeting definition.
    pub position: [Ball; 3],
    /// Strict source-piece tangent/target-normal enclosure used for admission.
    pub normal_dot: Ball,
}

/// Checked numerical reports, with no exact incidence or complete hit-count claim.
#[derive(Clone, Debug)]
pub struct NumericalCurveSurfaceMeetings {
    /// Every retained finite report, including shared-piece endpoint reports.
    pub hits: Vec<NumericalCurveSurfaceHit>,
    /// Fixed tolerance and residuals evaluated on both original definitions.
    pub evidence: NumericalEvidence,
    /// Charged SMLib search work for all source pieces.
    pub search_work: WorkReport,
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("intersect.numerical-curve-surface", detail)
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value).ok_or_else(|| bad("non-finite numerical source value"))
}

fn approximate(value: &Real, dec: &mut Decider) -> Result<f64, Refusal> {
    dec.charge_work(1)?;
    let value = dec.enclose_budgeted_at(value, BITS)?.mid().to_f64_nearest();
    exact(value)?;
    Ok(value)
}

fn upper(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}

fn source_surface(patch: &RationalPatch, dec: &mut Decider) -> Result<BSplineSurface, Refusal> {
    let source = patch.carrier();
    let mut controls = Vec::with_capacity(source.controls().len());
    for (point, weight) in source.controls().iter().zip(source.weights()) {
        dec.charge_work(1)?;
        let weight = approximate(weight, dec)?;
        if weight <= 0.0 {
            return Err(bad("source adaptation lost a positive patch weight"));
        }
        let coordinates = [
            approximate(&point.x, dec)? * weight,
            approximate(&point.y, dec)? * weight,
            approximate(&point.z, dec)? * weight,
            weight,
        ];
        if coordinates.iter().any(|value| !value.is_finite()) {
            return Err(bad("a homogeneous patch control overflowed binary64"));
        }
        controls.push(ControlPoint { coordinates });
    }
    let mut knots = Vec::with_capacity(2);
    for row in source.knots() {
        let mut values = Vec::with_capacity(row.len());
        for knot in row {
            dec.charge_work(1)?;
            let value = knot.to_f64_nearest();
            exact(value)?;
            values.push(value);
        }
        knots.push(values);
    }
    BSplineSurface::new(
        source.degrees(),
        knots.try_into().unwrap(),
        source.control_count(),
        controls,
    )
    .map_err(|error| bad(format!("finite patch adaptation failed: {error:?}")))
}

#[derive(Clone)]
struct RoundedPoint {
    coordinates: [f64; 4],
    error: Rational,
}

fn rounded_point(piece: &TracePiece, t: f64, dec: &mut Decider) -> Result<RoundedPoint, Refusal> {
    dec.charge_work(1)?;
    let parameter = Ball::from_rational(&exact(t)?, BITS);
    let position = piece.evaluate(&parameter, dec)?;
    let mut coordinates = [0.0; 4];
    coordinates[3] = 1.0;
    let mut error = Rational::zero();
    for axis in 0..3 {
        coordinates[axis] = position[axis].mid().to_f64_nearest();
        let written = Ball::from_rational(&exact(coordinates[axis])?, BITS);
        error = error + upper(&(&position[axis] - &written));
    }
    Ok(RoundedPoint { coordinates, error })
}

struct Segment {
    lo: f64,
    hi: f64,
    start: RoundedPoint,
    end: RoundedPoint,
}

fn source_curve(
    piece: &TracePiece,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<(BSplineCurve, [Ball; 3]), Refusal> {
    let start = rounded_point(piece, 0.0, dec)?;
    let end = rounded_point(piece, 1.0, dec)?;
    let mut controls = vec![ControlPoint {
        coordinates: start.coordinates,
    }];
    let mut knots = vec![0.0, 0.0];
    let mut pending = vec![Segment {
        lo: 0.0,
        hi: 1.0,
        start,
        end,
    }];
    let mut derivative = None;
    while let Some(segment) = pending.pop() {
        dec.charge_work(1)?;
        let lo = exact(segment.lo)?;
        let hi = exact(segment.hi)?;
        let interval = Ball::from_rational(&lo, BITS).hull(&Ball::from_rational(&hi, BITS));
        let mut jets = piece.taylor(&interval, 2, dec)?;
        let curvature = jets[2]
            .iter()
            .fold(Rational::zero(), |sum, value| sum + upper(value));
        if derivative.is_none() {
            derivative = Some(jets.swap_remove(1));
        }
        let width = &hi - &lo;
        // Taylor coefficient two is C''/2. Linear interpolation error is
        // bounded by sup|C''| h²/8, plus the rounded endpoints' displacement.
        let error = curvature * (&width * &width).mul_pow2(-2)
            + segment.start.error.clone().max(segment.end.error.clone());
        if error <= *tolerance {
            controls.push(ControlPoint {
                coordinates: segment.end.coordinates,
            });
            knots.push(segment.hi);
            continue;
        }
        let middle = segment.lo + (segment.hi - segment.lo) * 0.5;
        if middle <= segment.lo || middle >= segment.hi {
            return Err(bad(
                "the curve proposal cannot meet the preset tolerance at binary64 resolution",
            ));
        }
        let point = rounded_point(piece, middle, dec)?;
        pending.push(Segment {
            lo: middle,
            hi: segment.hi,
            start: point.clone(),
            end: segment.end,
        });
        pending.push(Segment {
            lo: segment.lo,
            hi: middle,
            start: segment.start,
            end: point,
        });
    }
    knots.push(1.0);
    let curve = BSplineCurve::new(1, knots, controls)
        .map_err(|error| bad(format!("finite trace adaptation failed: {error:?}")))?;
    Ok((
        curve,
        derivative.expect("the initial trace interval was enclosed"),
    ))
}

fn pay(units: u64, deadline: Option<std::time::Instant>) -> bool {
    let mut dec = Decider::default();
    dec.set_deadline(deadline);
    dec.charge_work(units).is_ok()
}

fn add_report(total: &mut WorkReport, report: WorkReport) {
    total.attempted += report.attempted;
    total.accepted += report.accepted;
    total.search_branches += report.search_branches;
    total.local_iterations += report.local_iterations;
    total.surface_evaluations += report.surface_evaluations;
    total.curve_evaluations += report.curve_evaluations;
    total.source_items += report.source_items;
    total.trace_visits += report.trace_visits;
    total.cancelled |= report.cancelled;
}

/// Searches finite pieces using SMLib's global subdivision and guarded local solver.
/// Curve proposals have derivative-bounded chord error below tolerance/8; every
/// reported incidence is checked on the original definitions at 128 bits.
/// Each report requires a strict source-piece tangent/target-normal enclosure.
/// Unresolved contacts, overlap, incomplete search or excessive residual refuse.
/// Requires one caller-owned work allowance; no tolerance or cap is escalated.
/// Neither an empty result nor its reports create exact geometry or incidence.
pub fn patch_meetings(
    trace: &CurveTrace,
    patch: &RationalPatch,
    tolerance_mm: &Rational,
    dec: &mut Decider,
) -> Result<NumericalCurveSurfaceMeetings, Refusal> {
    if !catallax_numbers::operation_budget_active() {
        return Err(bad(
            "finite numerical intersections require a caller-owned work allowance",
        ));
    }
    if tolerance_mm.sign() != Sign::Positive {
        return Err(bad("the intersection tolerance must be positive"));
    }
    let tolerance = tolerance_mm.to_f64_nearest();
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(bad("the intersection tolerance is outside finite binary64"));
    }
    let surface = source_surface(patch, dec)?;
    let domain = patch.domain();
    let float_domain = Extent2::new(
        Vector2::new(domain[0][0].to_f64_nearest(), domain[1][0].to_f64_nearest()),
        Vector2::new(domain[0][1].to_f64_nearest(), domain[1][1].to_f64_nearest()),
    );
    let mut hits = Vec::new();
    let mut samples = Vec::new();
    let mut search_work = WorkReport::default();
    for (index, piece) in trace.pieces().iter().enumerate() {
        let (curve, derivative) = source_curve(piece, &tolerance_mm.mul_pow2(-3), dec)?;
        let (solutions, work) = global_curve_intersect_with_work(
            &surface,
            float_domain,
            &curve,
            Extent1::new(0.0, 1.0),
            tolerance / 4.0,
            dec.deadline(),
            pay,
        );
        add_report(&mut search_work, work);
        if search_work.cancelled {
            dec.charge_work(0)?;
            return Err(bad(
                "finite intersection cancelled without a retained result",
            ));
        }
        let solutions = solutions.map_err(|error| {
            Refusal::new(
                "intersect.numerical-search-incomplete",
                format!("finite curve/surface search failed: {error:?}"),
            )
        })?;
        for solution in solutions.as_slice() {
            dec.charge_work(1)?;
            if !matches!(
                solution.kind,
                SolutionType::SingleValue | SolutionType::PolySingleValue
            ) {
                return Err(Refusal::new(
                    "intersect.numerical-overlap",
                    "a finite curve/surface overlap needs interval ownership",
                ));
            }
            let [t, u, v] = solution.start.parameters[..3].try_into().unwrap();
            let parameter = exact(t)?;
            let uv = [exact(u)?, exact(v)?];
            if parameter < Rational::zero()
                || parameter > Rational::one()
                || (0..2).any(|axis| uv[axis] < domain[axis][0] || uv[axis] > domain[axis][1])
            {
                return Err(bad("a numerical meeting leaves its finite source domains"));
            }
            let position = piece.evaluate(&Ball::from_rational(&parameter, BITS), dec)?;
            let at = uv.clone().map(Real::from);
            let point = patch.point_at(&at, dec)?;
            let target = [
                dec.enclose_budgeted_at(&point.x, BITS)?,
                dec.enclose_budgeted_at(&point.y, BITS)?,
                dec.enclose_budgeted_at(&point.z, BITS)?,
            ];
            let distance_upper_mm = (0..3).fold(Rational::zero(), |sum, axis| {
                sum + upper(&(&position[axis] - &target[axis]))
            });
            if distance_upper_mm > *tolerance_mm {
                return Err(bad(
                    "a reported incidence exceeds the preset original-source tolerance",
                ));
            }
            let [du, dv] = patch.derivatives_at(&at, dec)?;
            let normal = du.cross(&dv);
            let normal = [
                dec.enclose_budgeted_at(&normal.x, BITS)?,
                dec.enclose_budgeted_at(&normal.y, BITS)?,
                dec.enclose_budgeted_at(&normal.z, BITS)?,
            ];
            let normal_dot = (&normal[0] * &derivative[0])
                + (&normal[1] * &derivative[1])
                + (&normal[2] * &derivative[2]);
            if normal_dot.contains_zero() {
                return Err(Refusal::new(
                    "intersect.numerical-nontransverse",
                    "the original source-piece crossing is unresolved, tangent or overlapping",
                ));
            }
            samples.push(NumericalSample {
                parameters: vec![parameter.clone(), uv[0].clone(), uv[1].clone()],
                distance_upper_mm,
            });
            hits.push(NumericalCurveSurfaceHit {
                piece: index,
                parameter,
                uv,
                position,
                normal_dot,
            });
        }
    }
    let mut checks = vec![NumericalCheck::Passed {
        name: "curve-surface.completed-finite-numerical-search".into(),
        scope: "retained-trace-pieces-no-exact-hit-count-proof".into(),
    }];
    if !samples.is_empty() {
        checks.push(NumericalCheck::Passed {
            name: "curve-surface.strict-oriented-crossing-enclosure".into(),
            scope: "every-retained-source-piece-at-the-reported-target-normal".into(),
        });
        checks.push(NumericalCheck::Samples {
            name: "curve-surface.original-source-incidence".into(),
            scope: "finite-original-trace-and-rational-patch".into(),
            samples,
        });
    }
    let evidence = NumericalEvidence::new(tolerance_mm.clone(), checks, Vec::new())?;
    dec.check_deadline()?;
    Ok(NumericalCurveSurfaceMeetings {
        hits,
        evidence,
        search_work,
    })
}
