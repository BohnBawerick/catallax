//! The geometric gap of one edge use (decision Q9, and section 4.2 step 2 of
//! the Checkpoint 2 decisions): the largest distance from a point of one copy
//! to the other copy as a curve in space, taken both ways.
//!
//! [`measure`](super::measure) pairs the two copies at equal parameters. An
//! imported edge promises no such pairing, so its copies can lie on each
//! other and still be far apart at equal parameters. This reading pairs each
//! point with its foot on the other copy instead.
//!
//! The method, with its sources:
//!
//! 1. The two ends are paired with each other, first with first and last
//!    with last, as OpenCASCADE's `BRepLib_ValidateEdge::processApprox` pairs
//!    them. A copy that is slid along the other one shows here.
//! 2. Between the ends, each point of one copy is paired with its foot on the
//!    other. The foot is found by Gauss-Newton from the point at the same
//!    parameter, stays inside the edge's range, and takes no step that does
//!    not come nearer. So no pair is further apart than the same-parameter
//!    pair it started from, and the geometric gap is never read larger than
//!    the same-parameter gap. OpenCASCADE does the same with
//!    `Extrema_LocateExtPC` seeded at the proportional parameter, at 22 even
//!    points. Here the points are those of the same-parameter plan (steps 3
//!    and 4 of the checker), so their count depends on the shape of the edge
//!    alone.
//! 3. Every local maximum of the sampled foot distances is polished with
//!    Brent's method (step 6 of the checker).
//! 4. Both directions are read, 3D curve to face copy and face copy to 3D
//!    curve. The larger one is the geometric gap.
//!
//! The result decides whether an import is reparametrised. It is never
//! stored: what is stored is the same-parameter gap of the copy admission
//! keeps, measured by the full checker. A foot that settles on a wrong branch
//! is still a point of the other copy, so it can only read a distance too
//! large, which keeps a copy as it was imported.

use catallax_numbers::{Decider, Refusal};

use super::lift::{Bounded, Lifted, Scratch};
use super::{
    EdgeCopies, FaceCopyRef, Hints, Limits, Peak, PlannedPiece, Search, breaks, numeric, planned,
    polish, under_allowance, work_cap,
};

/// The geometric gap of one edge use.
#[derive(Clone, Copy, Debug)]
pub struct GeometricMeasure {
    /// The largest distance found from a point of the 3D curve to the lifted
    /// face copy, in millimetres.
    pub forward_mm: f64,
    /// The largest distance found from a point of the lifted face copy to
    /// the 3D curve, in millimetres.
    pub backward_mm: f64,
    /// An estimate of the binary64 rounding error of the larger of the two.
    pub error_mm: f64,
    /// The edge parameter of the point the larger one was read at, on that
    /// point's own copy.
    pub parameter: f64,
    /// Evaluations of either copy, each with the work a sampled carrier
    /// states for it.
    pub evaluations: u64,
    /// Work counted against the edge's cap: the planning and the evaluations.
    pub work: u64,
}

impl GeometricMeasure {
    /// The larger of the two directions.
    pub fn raw_mm(&self) -> f64 {
        self.forward_mm.max(self.backward_mm)
    }
    /// The reported geometric gap: the raw value plus its rounding error.
    pub fn reported_mm(&self) -> f64 {
        (self.raw_mm() + self.error_mm).next_up()
    }
}

/// Which copy the points are taken from.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Way {
    /// Points of the 3D curve against the lifted face copy.
    Forward,
    /// Points of the lifted face copy against the 3D curve.
    Backward,
}

/// One point of a copy and the point of the other copy it is paired with.
#[derive(Clone, Copy)]
pub(super) struct Foot {
    /// The parameter of the point on its own copy.
    pub at: f64,
    /// The parameter of its partner on the other copy.
    pub foot: f64,
    /// Their squared distance.
    pub squared: f64,
    /// The size of the rounding error of the two evaluations.
    pub error: f64,
}

impl Peak for Foot {
    fn at(&self) -> f64 {
        self.at
    }
    fn squared(&self) -> f64 {
        self.squared
    }
    fn error(&self) -> f64 {
        self.error
    }
}

/// Source: Dennis and Schnabel 1983, section 10.1, as for the foot on a
/// surface in `lift.rs`: Gauss-Newton converges fast from a seed near a foot
/// whose residual is small; eight steps settle such a seed, and more cannot
/// help a stalled one.
const FOOT_STEPS: usize = 8;

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn squared_distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let d: [f64; 3] = std::array::from_fn(|k| a[k] - b[k]);
    dot(&d, &d)
}

/// The point of one copy at `x`: the 3D curve, or the face copy lifted
/// through the surface.
fn copy_point(lifted: &Lifted, on_face: bool, x: f64, scratch: &mut Scratch) -> Bounded<[f64; 3]> {
    if !on_face {
        return lifted.curve.point(x, scratch);
    }
    let uv = lifted
        .face_at(x, scratch)
        .expect("a geometric reading is taken of a copy with its own u/v");
    lifted.surface.point(uv.value, scratch)
}

/// Takes one evaluation of a copy from what is left of an allowance, before
/// it is made: of the 3D curve, or of the face copy lifted through the
/// surface. A sampled carrier's stated work is taken with it.
fn spend(lifted: &Lifted, on_face: bool, left: &mut u64) -> Option<()> {
    if on_face {
        lifted.surface.spend(left)
    } else {
        lifted.curve.spend(left)
    }
}

/// The two copies at the same parameter. `None` when `left` runs out.
pub(super) fn pair(
    lifted: &Lifted,
    scratch: &mut Scratch,
    way: Way,
    x: f64,
    left: &mut u64,
) -> Option<Foot> {
    let on_face = way == Way::Backward;
    spend(lifted, on_face, left)?;
    let own = copy_point(lifted, on_face, x, scratch);
    spend(lifted, !on_face, left)?;
    let other = copy_point(lifted, !on_face, x, scratch);
    Some(Foot {
        at: x,
        foot: x,
        squared: squared_distance(&own.value, &other.value),
        error: own.error + other.error,
    })
}

/// The foot on the other copy of the point of one copy at `x`, inside the
/// edge's range, from the other copy's point at the same parameter:
/// Gauss-Newton on the squared distance with a forward-difference rate,
/// whose step is the square root of the machine epsilon times the size of
/// the parameter (Dennis and Schnabel 1983, sections 5.4 and 10.1).
/// OpenCASCADE's `ProjectPointOnCurve` in `Approx_SameParameter.cxx` is the
/// same iteration with the range clamp. A step that does not come nearer is
/// not taken, so the pair is never further apart than the same-parameter
/// pair. Each evaluation is taken from `left` before it is made; `None` when
/// `left` runs out.
pub(super) fn foot(
    lifted: &Lifted,
    scratch: &mut Scratch,
    way: Way,
    x: f64,
    left: &mut u64,
) -> Option<Foot> {
    let [lo, hi] = lifted.range();
    let on_face = way == Way::Backward;
    spend(lifted, on_face, left)?;
    let target = copy_point(lifted, on_face, x, scratch);
    let mut s = x.clamp(lo, hi);
    spend(lifted, !on_face, left)?;
    let mut base = copy_point(lifted, !on_face, s, scratch);
    let mut squared = squared_distance(&target.value, &base.value);
    for _ in 0..FOOT_STEPS {
        let scale = s.abs().max(hi - lo);
        let mut step = f64::EPSILON.sqrt() * scale;
        if s + step > hi {
            step = -step;
        }
        spend(lifted, !on_face, left)?;
        let shifted = copy_point(lifted, !on_face, s + step, scratch);
        let rate: [f64; 3] = std::array::from_fn(|k| (shifted.value[k] - base.value[k]) / step);
        let residual: [f64; 3] = std::array::from_fn(|k| target.value[k] - base.value[k]);
        let (along, speed) = (dot(&rate, &residual), dot(&rate, &rate));
        if !(speed > 0.0 && along.is_finite() && speed.is_finite()) {
            break;
        }
        let next = (s + along / speed).clamp(lo, hi);
        if (next - s).abs() <= 0.0 {
            break;
        }
        spend(lifted, !on_face, left)?;
        let there = copy_point(lifted, !on_face, next, scratch);
        let moved = squared_distance(&target.value, &there.value);
        if moved >= squared {
            break;
        }
        let settled = (next - s).abs() <= f64::EPSILON * scale;
        (s, base, squared) = (next, there, moved);
        if settled {
            break;
        }
    }
    Some(Foot {
        at: x,
        foot: s,
        squared,
        error: target.error + base.error,
    })
}

/// Records what a reading took of the edge's allowance, whether or not it
/// finished, and refuses by name where the allowance ran out.
fn recorded(
    search: &mut Search<'_>,
    before: u64,
    left: u64,
    found: Option<Foot>,
) -> Result<Foot, Refusal> {
    search.reserve(before - left)?;
    found.ok_or_else(|| {
        work_cap("an edge spent its work allowance before its geometric gap was read")
    })
}

fn read_pair(search: &mut Search<'_>, way: Way, x: f64) -> Result<Foot, Refusal> {
    let before = search.left();
    let mut left = before;
    let found = pair(search.lifted, &mut search.scratch, way, x, &mut left);
    recorded(search, before, left, found)
}

fn read_foot(search: &mut Search<'_>, way: Way, x: f64) -> Result<Foot, Refusal> {
    let before = search.left();
    let mut left = before;
    let found = foot(search.lifted, &mut search.scratch, way, x, &mut left);
    recorded(search, before, left, found)
}

/// The largest distance from a point of one copy to the other copy.
fn sweep(
    search: &mut Search<'_>,
    plan: &[PlannedPiece],
    way: Way,
    dec: &mut Decider,
) -> Result<Foot, Refusal> {
    let [lo, hi] = search.lifted.range();
    let mut samples = vec![read_pair(search, way, lo)?];
    for piece in plan {
        let [from, to] = piece.piece;
        for share in numeric::chebyshev_shares(piece.points) {
            samples.push(read_foot(search, way, from + share * (to - from))?);
        }
        // The far end of the piece: a break inside the edge, or its last end.
        samples.push(if to >= hi {
            read_pair(search, way, to)?
        } else {
            read_foot(search, way, to)?
        });
        search.settle(dec)?;
    }
    polish(
        search,
        &samples,
        |search, x, _| read_foot(search, way, x),
        dec,
    )
}

/// Measures the geometric gap of one edge use, both ways.
///
/// `None` for a curve projected on a plane: that copy has no u/v of its own,
/// and its same-parameter gap is already the distance of the 3D curve from
/// the plane.
///
/// `limits` is what the use's edge has left of its one allowance when another
/// measurement of the same edge came first. One unit is counted per break and
/// per halving while planning, as [`measure`](super::measure) counts them,
/// and one per evaluation of either copy. It refuses with
/// `numerical.gap-work-cap` where the allowance runs out, and the same units
/// are charged to the operation's work.
pub fn measure_geometric(
    copies: &EdgeCopies<'_>,
    hints: &Hints,
    limits: Limits,
    dec: &mut Decider,
) -> Result<Option<GeometricMeasure>, Refusal> {
    if matches!(copies.face_copy, FaceCopyRef::Projected) {
        return Ok(None);
    }
    under_allowance(copies, limits, dec, |search, dec| {
        let breaks = breaks(search, hints, dec)?;
        let plan = planned(search.lifted, &breaks);
        let planning = search.spent();
        let forward = sweep(search, &plan, Way::Forward, dec)?;
        let backward = sweep(search, &plan, Way::Backward, dec)?;
        search.settle(dec)?;
        let larger = if backward.squared > forward.squared {
            backward
        } else {
            forward
        };
        Ok(GeometricMeasure {
            forward_mm: forward.squared.sqrt(),
            backward_mm: backward.squared.sqrt(),
            error_mm: larger.error,
            parameter: larger.at,
            evaluations: search.spent() - planning,
            work: search.spent(),
        })
    })
    .map(Some)
}
