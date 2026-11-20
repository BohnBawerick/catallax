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
//!    for a line on a plane or along a generator, as in OpenCASCADE's
//!    `ElementaryMaximum`.
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
//! All arithmetic is IEEE basic operations in one thread, so the result is
//! the same bits on every machine and at any thread count.

use catallax_numbers::{Decider, Rational, Refusal};
use catallax_shapes::{BSplineCurve, BSplineSurface, Curve, Surface};

mod lift;
mod numeric;

use lift::{Curve64, Face64, Lifted, Revolved, Scratch, Surface64};

/// The surface a face copy is lifted through.
#[derive(Clone, Copy, Debug)]
pub enum GapSurface<'a> {
    /// A kernel surface. A spline support is evaluated through its carrier.
    Kernel(&'a Surface),
    /// A spline carrier on its own.
    Carrier(&'a BSplineSurface),
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
    /// The 3D copy: a line, a circle, an ellipse or a B-spline.
    pub curve: &'a Curve,
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
    /// residual evaluation, one sample and mode pair of a tail transform, or
    /// one domain probe: everything the measurement does that grows with the
    /// edge.
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
    /// Work counted against the edge's cap: the evaluations, the tail
    /// transforms and the domain probes.
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

fn work_cap(detail: &str) -> Refusal {
    Refusal::new("numerical.gap-work-cap", detail)
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
        },
    }
}

/// Source: IEEE 754 binary64 has 53 significant bits. Two breaks closer than
/// 2^-45 of the edge are within 256 units in the last place of each other at
/// the edge's own scale and leave no room for an interior sample.
const NARROWEST_PIECE: f64 = 1.0 / ((1_u64 << 45) as f64);

/// The values of one surface parameter, within `range`, at which the lifted
/// curve is cut: knot lines and domain borders of a spline surface, the seam
/// of a periodic chart, and the pole lines of a sphere or a cone.
///
/// The number of seams is fixed before any is listed. A face copy that winds
/// around its surface more than `most` times is refused, not walked.
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
                // The seam, u = 0 modulo one turn. One spare turn covers
                // the rounding of the two quotients.
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
    }
    out.dedup();
    Ok(())
}

/// Parameters in `(from, to)` where the face copy crosses a cut line of the
/// surface. The face copy is one polynomial or rational piece on the interval.
/// Returns the work spent.
fn crossings(
    lifted: &Lifted,
    from: f64,
    to: f64,
    most: u64,
    out: &mut Vec<f64>,
) -> Result<u64, Refusal> {
    let mut work = 0;
    // Homogeneous Bezier coefficients (u w, v w, w) of the face copy.
    let piece: Vec<[f64; 3]> = match &lifted.face {
        Face64::Projected => return Ok(work),
        Face64::Linear { ends } => {
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
        cut_lines(&lifted.surface, axis, [lo, hi], most, &mut lines)?;
        for line in &lines {
            let coefficients: Vec<f64> = piece.iter().map(|c| c[axis] - line * c[2]).collect();
            let size = piece
                .iter()
                .fold(line.abs(), |sum, c| sum.max((c[axis] / c[2]).abs()));
            shares.clear();
            let halvings = numeric::sign_changes(
                &coefficients,
                numeric::UNIT_ROUNDOFF * numeric::roundings(piece.len() - 1) * size,
                &mut shares,
            );
            work += 1 + halvings;
            out.extend(shares.iter().map(|share| from + share * (to - from)));
        }
    }
    Ok(work)
}

/// Every break of the edge in ascending order, both ends included. Work is
/// charged window by window, and an edge with more than `most` breaks is
/// refused: each piece takes at least one evaluation, so it could not be
/// sampled within the allowance.
fn breaks(
    lifted: &Lifted,
    hints: &Hints,
    most: u64,
    dec: &mut Decider,
) -> Result<Vec<f64>, Refusal> {
    let [lo, hi] = lifted.range();
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return Err(unsupported("an edge has no finite parameter range"));
    }
    let inside = |t: &f64| *t > lo && *t < hi;
    let mut cuts = vec![lo, hi];
    if let Curve64::Spline(spline) = &lifted.curve {
        cuts.extend(spline.knots.iter().copied().filter(inside));
    }
    if let Face64::Native(spline) = &lifted.face {
        cuts.extend(spline.knots.iter().copied().filter(inside));
    }
    cuts.extend(hints.parameters.iter().copied().filter(inside));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let too_many = || work_cap("an edge has more breaks than its evaluation allowance can sample");
    if cuts.len() as u64 > most {
        return Err(too_many());
    }
    dec.charge_work(cuts.len() as u64)?;
    // Between these cuts the face copy is one polynomial or rational piece.
    let mut crossed = Vec::new();
    for pair in cuts.windows(2) {
        let work = crossings(lifted, pair[0], pair[1], most, &mut crossed)?;
        if work > 0 {
            dec.charge_work(work)?;
        }
        if (cuts.len() + crossed.len()) as u64 > most {
            return Err(too_many());
        }
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
    if straight_residual(&lifted) {
        return Ok(Vec::new());
    }
    Ok(planned(&lifted, &breaks(&lifted, hints, limits.work, dec)?))
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

    fn sample(&mut self, t: f64) -> Result<Sample, Refusal> {
        if self.spent() >= self.limit {
            return Err(work_cap(
                "an edge spent its work allowance before its gap was resolved",
            ));
        }
        self.evaluations += 1;
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
        if units > self.limit.saturating_sub(self.spent()) {
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
    let lifted = lift::lift(copies, dec)?;
    // The closing probes of the rounding allowance count inside the same cap,
    // so they are set aside before the first sample.
    let Some(limit) = limits.work.checked_sub(lift::ALLOWANCE_EVALUATIONS) else {
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
    let [first, last] = lifted.parameters;
    if straight_residual(&lifted) {
        let ends = [search.sample(first)?, search.sample(last)?];
        let best = if ends[1].squared > ends[0].squared {
            ends[1]
        } else {
            ends[0]
        };
        // Convexity makes the end values the maximum: nothing is unresolved.
        return finish(
            &mut search,
            best,
            (0.0, None, Method::StraightResidual),
            (Vec::new(), Vec::new()),
            dec,
        );
    }
    let breaks = breaks(&lifted, hints, limit, dec)?;
    let plan = planned(&lifted, &breaks);
    let [lo, hi] = lifted.range();
    // Refuse before the first sample when the plan alone is over the cap:
    // its samples, the first tail transform of each piece and, on a spline
    // surface, one domain probe per piece.
    let probes = u64::from(matches!(lifted.surface, Surface64::Spline(_)));
    let planned_work = plan.iter().fold(1_u64, |sum, piece| {
        sum + piece.points as u64 + 1 + numeric::chebyshev_tail_pairs(piece.points) + probes
    });
    if planned_work > limit {
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
            let left = unresolved(&mut search, &grid)?;
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
    let mut best = samples[0];
    for sample in &samples {
        if sample.squared > best.squared {
            best = *sample;
        }
    }
    let mut i = 0;
    while i < samples.len() {
        // A run of equal values is one candidate.
        let mut j = i;
        while j + 1 < samples.len()
            && samples[j + 1].squared.to_bits() == samples[i].squared.to_bits()
        {
            j += 1;
        }
        let left = i.checked_sub(1);
        let right = (j + 1 < samples.len()).then_some(j + 1);
        let peak = left.is_none_or(|l| samples[l].squared < samples[i].squared)
            && right.is_none_or(|r| samples[r].squared < samples[i].squared);
        // A value inside its own rounding error is noise, not a maximum.
        let above_noise = samples[i].squared.sqrt() > samples[i].error;
        let (from, to) = (samples[left.unwrap_or(i)].t, samples[right.unwrap_or(j)].t);
        let resolvable = to - from > 2.0 * NARROWEST_PIECE * (hi - lo);
        if peak && above_noise && resolvable {
            let mut failure = None;
            let mut polished = samples[i];
            let (_, _, _) = numeric::brent_maximum(&mut |share| match search
                .sample(from + share * (to - from))
            {
                Ok(sample) => {
                    if sample.squared > polished.squared {
                        polished = sample;
                    }
                    sample.squared
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
            if polished.squared > best.squared {
                best = polished;
            }
        }
        i = j + 1;
    }

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
        &mut search,
        best,
        (search_error, exit, Method::Sampled),
        (plan, refinements),
        dec,
    )
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
    search.settle(dec)?;
    Ok(GapMeasure {
        raw_mm: best.squared.sqrt(),
        error_mm: allowance + search_error,
        parameter: best.t,
        domain_exit,
        method,
        plan,
        refinements,
        evaluations: search.evaluations,
        work: search.spent(),
    })
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
            Face64::Linear { ends },
        ) => (ends[1][0] - ends[0][0]).abs() <= 0.0,
        _ => false,
    }
}
