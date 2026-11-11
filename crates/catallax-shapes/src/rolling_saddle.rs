//! A ball rolling between nonparallel cylinders.
//!
//! Take the wide axis through the origin along `axis`, the common perpendicular
//! along `across` and `along` square to both, so the narrow axis runs through
//! `(gap, 0, 0)` in the direction `(0, S, C)`. Offsetting each cylinder by the
//! ball radius gives radii `P` and `p`. A ball centre on the narrow offset at
//! angle `t` has `x = gap + p cos t`, and the wide offset fixes
//!
//! ```text
//! y = branch sqrt(P^2 - x^2)
//! h = (y - p C sin t) / S
//! z = (C y - p sin t) / S
//! ```
//!
//! with `h` the station along the narrow axis. When `|gap| + p < P` the root
//! never closes, so `t` charts one whole ring. Two rational angles and that one
//! root are the whole law, and the root makes the ring an elliptic curve, so no
//! rational surface matches it and the carrier is a fit.
//!
//! The carrier has three rows. Its two rails are rational curves exactly on
//! their cylinders: the narrow rail keeps the exact angle `t` with a fitted
//! height, and the wide rail takes a fitted half-angle tangent and height. The
//! middle row is the apex of each section's conic, held in the narrow cylinder's
//! tangent plane with a fitted middle weight, so every section is a true conic.
//! The two-sided bound is recomputed from the retained node data and interval
//! Taylor jets of the law, as for [`RollingEnvelope`](crate::RollingEnvelope).
//!
//! `u` runs over `[-1, 3]`. On `[-1, 1]` the half-angle tangent of `t` is `u`, and
//! on `[1, 3]` the same half turn repeats from the opposite side, so the rail's
//! denominator is two at every chart end and the closed carrier's weights close.
//!
//! When the two offsets' shadows cross instead, so `|P - p| < |gap| < P + p`, one
//! closed loop runs out to where `y = 0` and back round to where `sin t = 0`, and
//! both `t` and the wide angle turn round on it. With `e` the sign of the gap, take
//! `X = e x`, which runs between `P` and `|gap| - p`, and write
//!
//! ```text
//! X = c + d cos 2 w    c = (P + |gap| - p) / 2    d = (P - |gap| + p) / 2
//! y = sin w sqrt(2 d (P + X))
//! p sin t = cos w sqrt(2 d (|gap| + p - X))
//! ```
//!
//! Both roots stay positive, so `w` charts the loop with no turning point, as `t`
//! charts a through ring. The rails take fitted half-angle tangents measured from
//! the far side of each cylinder, `2 d sin w / sqrt(2 d (P + X))` on the wide one
//! and `-2 d cos w / sqrt(2 d (|gap| + p - X))` on the narrow one, which stay
//! finite because the loop never reaches those far sides. Each rail's denominator
//! is then of degree six, so the middle row reaches degree twelve.
//!
//! Equal radii on crossing axes instead give four pole-to-pole bands. Each
//! has rational contact rays and an exact degree (6,2) carrier, with two Bezier
//! spans and collapsed parameter ends. Its whole-polynomial certificate and
//! strict tube-reach proof replace the regular ring's fit and width estimate;
//! see `collapsed`.

use crate::rational_sheet::{LocalFrame, Patch, along, cover, split_net};
use crate::rolling_envelope::{Jet, combine, cubic, elevate, hermite, product};
use crate::{
    BSplineCurve, BSplineSurface, Cylinder, EnvelopeGrid, Isometry, Plane, Point3, Pose,
    SplineBoundary, Vector3, math,
};
use catallax_numbers::{Ball, Decider, EvidenceBound, EvidenceUnit, Rational, Real, Refusal, Sign};
use std::sync::{Arc, Mutex};

#[path = "rolling_saddle/collapsed.rs"]
mod collapsed;
#[path = "rolling_saddle/grid.rs"]
mod grid;
#[path = "rolling_saddle/line.rs"]
mod line;
#[path = "rolling_saddle/table.rs"]
mod table;

/// Kernel-made approximation target, as for a rolling envelope.
pub const ROLLING_SADDLE_TARGET_NM: i64 = 100;
const HARD_CAP_NM: i64 = 1_000;
// Spans per half turn: the fit doubles from the first to the last.
const FIRST_SPANS: usize = 8;
const LAST_SPANS: usize = 128;
// Halvings of one of a loop's first spans. A loop's section angle can come far
// nearer a half turn, where the rows bend hardest, so its spans follow the bound.
const MAX_SPLITS: u32 = 7;
const DEGREE: usize = 9;
const OVERLAP_DEGREE: usize = 12;
const BITS: u32 = 128;
// Node data keep this many fractional bits, far below the target.
const NODE_BITS: i64 = 32;
// Cells of one span in the regularity proof, by doubling.
const MAX_CELLS: usize = 64;

const ROWS: usize = 7;
const HEIGHT: usize = 0;
const LIFT: usize = 1;
const HALF: usize = 2;
const TANGENT: usize = 3;
const AXIAL: usize = 4;
const WEIGHT: usize = 5;
// The narrow rail's half-angle tangent, on an overlap's loop only.
const TURNED: usize = 6;

// Halvings of the open patches before a range gives up.
const MAX_RANGE_SPLITS: usize = 40;

fn invalid(detail: impl Into<String>) -> Refusal {
    Refusal::new("shapes.rolling-saddle", detail)
}
fn unresolved(detail: impl Into<String>) -> Refusal {
    Refusal::new("numbers.geometric-gap-unproven", detail)
}
fn b(value: &Rational, bits: u32) -> Ball {
    Ball::from_rational(value, bits)
}
fn q(value: i64) -> Rational {
    Rational::from_i64(value)
}
fn magnitude(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}
fn snap(value: &Rational) -> Rational {
    Rational::from(value.mul_pow2(NODE_BITS).floor()).mul_pow2(-NODE_BITS)
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &math::exact_sqrt(&v.norm_squared())))
}
fn real(value: &Rational) -> Real {
    Real::from(value.clone())
}

/// Which contact component the law charts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SaddleRing {
    /// The narrow offset passes through the wide one: `t` charts one of two rings.
    Through,
    /// The offsets' shadows cross, so one loop turns round in `t` and in the wide
    /// angle.
    Overlap,
    /// One pole-to-pole band of equal cylinders on crossing axes. The two
    /// parameter ends collapse to shared contacts; `branch` selects the wide
    /// side and this flag selects the narrow side.
    Crossing {
        /// Whether the selected half lies on the positive narrow-side axis.
        narrow_positive: bool,
    },
}

/// The side of a cylinder the ball sits on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BallSide {
    /// Outside the cylinder, so the offset radius grows.
    Outside,
    /// Inside the cylinder, so the offset radius shrinks.
    Inside,
}
impl BallSide {
    fn sign(self) -> i64 {
        match self {
            Self::Outside => 1,
            Self::Inside => -1,
        }
    }
}

/// The numbers of a [`RollingSaddleDefinition`], for code that integrates over the
/// law. `across`, `along` and `axis` are the frame of the module note, and may be
/// left-handed on a mirrored copy.
#[derive(Clone, Debug)]
pub struct SaddleLaw {
    /// The narrow axis's closest point to the wide axis.
    pub foot: Point3,
    /// Unit, from the wide axis to the narrow axis.
    pub across: Vector3,
    /// Unit, the narrow axis's part square to the wide axis.
    pub along: Vector3,
    /// The unit wide axis.
    pub axis: Vector3,
    /// The wide cylinder's radius.
    pub wide: Real,
    /// The narrow cylinder's radius.
    pub narrow: Real,
    /// The ball's radius.
    pub radius: Real,
    /// The signed distance between the axes, along `across`.
    pub gap: Real,
    /// Cosine of the angle between the axes.
    pub cosine: Real,
    /// Sine of the angle between the axes, positive.
    pub sine: Real,
    /// The ball's side of the wide and narrow cylinders: 1 outside, -1 inside.
    pub sides: [i64; 2],
    /// The sign of the ring's `along` coordinate. An overlap's loop has no side.
    pub branch: i64,
    /// Which ring this is.
    pub ring: SaddleRing,
    /// The gap's sign, `1` when the narrow axis lies toward `across`.
    pub gap_sign: i64,
}

impl SaddleLaw {
    /// The narrow axis.
    pub fn narrow_axis(&self) -> Vector3 {
        self.axis
            .scaled(&self.cosine)
            .plus(&self.along.scaled(&self.sine))
    }

    /// The direction square to the narrow axis and to `across`.
    pub fn side_axis(&self) -> Vector3 {
        self.along
            .scaled(&self.cosine)
            .minus(&self.axis.scaled(&self.sine))
    }

    /// For a rail's table, the direction from its cylinder's axis at angle zero,
    /// and the way a growing angle turns it. `0` is the wide rail.
    pub fn rail_zero(&self, rail: usize) -> [Vector3; 2] {
        let facing = Real::from_i64(self.gap_sign);
        match (self.ring, rail) {
            (SaddleRing::Through, 0) => [
                self.along.scaled(&Real::from_i64(self.branch)),
                self.across.clone(),
            ],
            (SaddleRing::Through, _) => [self.across.clone(), self.side_axis()],
            (SaddleRing::Overlap, 0) => [self.across.scaled(&facing), self.along.clone()],
            (SaddleRing::Overlap, _) => [
                self.across.scaled(&-&facing),
                self.side_axis().scaled(&Real::from_i64(-1)),
            ],
            (SaddleRing::Crossing { .. }, 0) => [
                self.along.scaled(&Real::from_i64(self.branch)),
                self.across.clone(),
            ],
            (SaddleRing::Crossing { narrow_positive }, _) => [
                self.side_axis()
                    .scaled(&Real::from_i64(if narrow_positive { 1 } else { -1 })),
                self.across.clone(),
            ],
        }
    }
}

/// The exact law of a constant-radius ball resting on two nonparallel
/// cylinders: a regular closed ring or an equal-radius pole-to-pole band.
#[derive(Clone, Debug)]
pub struct RollingSaddleDefinition {
    origin: Point3,
    foot: Point3,
    // The cylinders' own axes, which keep a rational line solver in reach.
    axes: [Vector3; 2],
    across: Vector3,
    along: Vector3,
    axis: Vector3,
    wide: Real,
    narrow: Real,
    radius: Real,
    gap: Real,
    cosine: Real,
    sine: Real,
    sides: [i64; 2],
    branch: i64,
    ring: SaddleRing,
    // The sign of `gap`, which a loop's chart reads.
    gap_positive: bool,
}

impl RollingSaddleDefinition {
    /// Checks the regular case: the offsets are positive, the axes are not
    /// parallel, and the narrow offset and cylinder both lie strictly inside the
    /// wide one when seen along the narrow axis. `leaving` picks the ring the
    /// narrow axis direction points out of, and `false` picks the other.
    pub fn new(
        wide: &Cylinder,
        narrow: &Cylinder,
        radius: &Real,
        sides: [BallSide; 2],
        leaving: bool,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        if math::sign(radius, dec)? != Sign::Positive {
            return Err(invalid("ball radius must be positive"));
        }
        let frame = Self::frame_of(wide, narrow, radius, sides, dec)?;
        let reach = math::exact_sqrt(&frame.gap.square());
        for (big, small) in [
            (frame.offset(0), frame.offset(1)),
            (frame.wide.clone(), frame.narrow.clone()),
        ] {
            if math::sign(&(&(&big - &small) - &reach), dec)? != Sign::Positive {
                return Err(invalid(
                    "the narrow cylinder does not pass through the wide one",
                ));
            }
        }
        Ok(Self {
            branch: if leaving { 1 } else { -1 },
            ..frame
        })
    }

    /// The one loop of two cylinders whose shadows cross along the narrow axis.
    /// Checks that the offsets' shadows and the cylinders' own both cross, so
    /// both circles reach past each other and neither holds the other.
    pub fn overlapping(
        wide: &Cylinder,
        narrow: &Cylinder,
        radius: &Real,
        sides: [BallSide; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        if math::sign(radius, dec)? != Sign::Positive {
            return Err(invalid("ball radius must be positive"));
        }
        let frame = Self::frame_of(wide, narrow, radius, sides, dec)?;
        let reach = math::exact_sqrt(&frame.gap.square());
        for (big, small) in [
            (frame.offset(0), frame.offset(1)),
            (frame.wide.clone(), frame.narrow.clone()),
        ] {
            let apart = math::exact_sqrt(&(&big - &small).square());
            if math::sign(&(&reach - &apart), dec)? != Sign::Positive
                || math::sign(&(&(&big + &small) - &reach), dec)? != Sign::Positive
            {
                return Err(invalid("the cylinders do not cross in one loop"));
            }
        }
        Ok(Self {
            ring: SaddleRing::Overlap,
            ..frame
        })
    }

    // The degree of a span of this ring along `u`.
    fn degree(&self) -> usize {
        match self.ring {
            SaddleRing::Through => DEGREE,
            SaddleRing::Overlap => OVERLAP_DEGREE,
            SaddleRing::Crossing { .. } => 6,
        }
    }

    // The frame of the common perpendicular, with both offsets checked positive.
    fn frame_of(
        wide: &Cylinder,
        narrow: &Cylinder,
        radius: &Real,
        sides: [BallSide; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let a = unit(&wide.axis);
        let n = unit(&narrow.axis);
        let cosine = a.dot(&n);
        let square = &Real::one() - &cosine.square();
        if math::sign(&square, dec)? != Sign::Positive {
            return Err(invalid("parallel cylinders have no saddle"));
        }
        let sine = math::exact_sqrt(&square);
        let along = n.minus(&a.scaled(&cosine)).scaled(&(&Real::one() / &sine));
        let across = along.cross(&a);
        // The closest points of the two axes, as rationals when the input is.
        let between = narrow.origin.minus(&wide.origin);
        let (long, short) = (&wide.axis, &narrow.axis);
        let (aa, cc, ac) = (long.dot(long), short.dot(short), long.dot(short));
        let (ad, cd) = (long.dot(&between), short.dot(&between));
        let det = &(&aa * &cc) - &ac.square();
        let on_wide = &(&(&cc * &ad) - &(&ac * &cd)) / &det;
        let on_narrow = &(&(&ac * &ad) - &(&aa * &cd)) / &det;
        let origin = wide.origin.plus(&long.scaled(&on_wide));
        let foot = narrow.origin.plus(&short.scaled(&on_narrow));
        let gap = foot.minus(&origin).dot(&across);
        let gap_positive = math::sign(&gap, dec)? != Sign::Negative;
        let radii = [&wide.radius_sq, &narrow.radius_sq].map(math::exact_sqrt);
        let [wide_radius, narrow_radius] = radii;
        let offset =
            |radius_of: &Real, side: BallSide| radius_of + &(radius * &Real::from_i64(side.sign()));
        let offsets = [
            offset(&wide_radius, sides[0]),
            offset(&narrow_radius, sides[1]),
        ];
        for value in &offsets {
            if math::sign(value, dec)? != Sign::Positive {
                return Err(invalid("the ball does not fit inside a cylinder"));
            }
        }
        Ok(Self {
            origin,
            foot,
            axes: [wide.axis.clone(), narrow.axis.clone()],
            across,
            along,
            axis: a,
            wide: wide_radius,
            narrow: narrow_radius,
            radius: radius.clone(),
            gap,
            cosine,
            sine,
            sides: sides.map(BallSide::sign),
            branch: 1,
            ring: SaddleRing::Through,
            gap_positive,
        })
    }

    /// The ring through a point of both cylinders, under the checks of
    /// [`Self::new`]. A regular ring never reaches the plane through the wide
    /// axis and the common perpendicular, so the point's side of it names the ring.
    pub fn through(
        wide: &Cylinder,
        narrow: &Cylinder,
        radius: &Real,
        sides: [BallSide; 2],
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let definition = Self::new(wide, narrow, radius, sides, true, dec)?;
        let [_, y, _] = definition.local(point);
        match math::sign(&y, dec)? {
            Sign::Positive => Ok(definition),
            Sign::Negative => Ok(Self {
                branch: -1,
                ..definition
            }),
            Sign::Zero => Err(invalid("the point is not on a regular ring")),
        }
    }

    /// The contact component through a point of both cylinders: a through ring,
    /// an overlap loop, or an equal-radius band between shared poles. The last
    /// case requires equal signed offsets and a strict centre-ellipse reach proof.
    pub fn crossing(
        wide: &Cylinder,
        narrow: &Cylinder,
        radius: &Real,
        sides: [BallSide; 2],
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let frame = Self::frame_of(wide, narrow, radius, sides, dec)?;
        let same_sides = frame.sides[0] == frame.sides[1];
        let collapsed =
            match collapsed::all_zero(&[&frame.wide - &frame.narrow, frame.gap.clone()], dec)? {
                collapsed::Zero::Yes => same_sides,
                collapsed::Zero::No => false,
                // Different sides rule the equal-crossing band out by themselves.
                collapsed::Zero::Undecided(_) if !same_sides => false,
                collapsed::Zero::Undecided(refusal) => return Err(refusal),
            };
        if collapsed {
            return collapsed::through(frame, point, dec);
        }
        let reach = math::exact_sqrt(&frame.gap.square());
        let clear = &(&frame.wide - &frame.narrow) - &reach;
        if math::sign(&clear, dec)? == Sign::Positive {
            Self::through(wide, narrow, radius, sides, point, dec)
        } else {
            Self::overlapping(wide, narrow, radius, sides, dec)
        }
    }

    /// Whether a tangent at a point of the sharp edge the chart rounds runs the
    /// way a growing `u` does.
    ///
    /// A ring through the narrow cylinder turns counter-clockwise about its axis
    /// throughout, so one tangent decides the direction of a whole ring.
    ///
    /// A loop has no such turn, because it does not wind round either cylinder.
    /// With `s` a point's height on the side axis and `y` its height on `along`,
    /// the law gives `(s, y) = (B cos w, A sin w)`, with
    /// `A^2 = 2 d (P + X)` and `B^2 = 2 d (g + p - X)`. With `S = sin^2 w`,
    ///
    /// ```text
    /// A B (s y' - y s') = 16 d^2 (P (p - d) - 2 d (p - d) S + d g S^2)
    /// ```
    ///
    /// where `0 < d < min(P, p)` and `p - d < g`. Its discriminant is
    /// `4 d (p - d) (d (p - d) - g P) < 0`, so the shadow on `(s, y)` turns one
    /// way about the origin throughout, the way `w` grows. The sharp edge has the
    /// same law with the cylinders' own radii, so the same test holds on it.
    pub fn runs_with_chart(
        &self,
        point: &Point3,
        tangent: &Vector3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let turn = match self.ring {
            SaddleRing::Through => self
                .narrow_axis()
                .cross(&self.radial(1, point))
                .dot(tangent),
            SaddleRing::Overlap => {
                let [_, y, _] = self.local(point);
                let side = self.law().side_axis();
                let s = point.minus(&self.foot).dot(&side);
                &(&s * &tangent.dot(&self.along)) - &(&y * &tangent.dot(&side))
            }
            SaddleRing::Crossing { .. } => {
                let [x, y, _] = self.local(point);
                (&y * &tangent.dot(&self.across) - &x * &tangent.dot(&self.along))
                    * Real::from_i64(self.branch)
            }
        };
        match math::sign(&turn, dec)? {
            Sign::Positive => Ok(true),
            Sign::Negative => Ok(false),
            Sign::Zero => Err(invalid("the tangent does not turn round the narrow axis")),
        }
    }

    fn mapped(&self, pose: &Isometry) -> Self {
        Self {
            origin: pose.apply_point(&self.origin),
            foot: pose.apply_point(&self.foot),
            axes: [
                pose.apply_vector(&self.axes[0]),
                pose.apply_vector(&self.axes[1]),
            ],
            across: pose.apply_vector(&self.across),
            along: pose.apply_vector(&self.along),
            axis: pose.apply_vector(&self.axis),
            ..self.clone()
        }
    }

    /// The ball radius.
    pub fn radius(&self) -> &Real {
        &self.radius
    }

    /// The law's origin, on the wide axis.
    pub fn origin(&self) -> &Point3 {
        &self.origin
    }

    /// The law's numbers.
    pub fn law(&self) -> SaddleLaw {
        SaddleLaw {
            foot: self.foot.clone(),
            across: self.across.clone(),
            along: self.along.clone(),
            axis: self.axis.clone(),
            wide: self.wide.clone(),
            narrow: self.narrow.clone(),
            radius: self.radius.clone(),
            gap: self.gap.clone(),
            cosine: self.cosine.clone(),
            sine: self.sine.clone(),
            sides: self.sides,
            branch: self.branch,
            ring: self.ring,
            gap_sign: if self.gap_positive { 1 } else { -1 },
        }
    }

    fn offset(&self, which: usize) -> Real {
        let radius = if which == 0 { &self.wide } else { &self.narrow };
        radius + &(&self.radius * &Real::from_i64(self.sides[which]))
    }
    fn frame(&self, bits: u32, dec: &Decider) -> Result<LocalFrame, Refusal> {
        LocalFrame::new(
            &self.origin,
            [&self.across, &self.along, &self.axis],
            bits,
            dec,
        )
    }
    fn budgeted_frame(&self, bits: u32, dec: &mut Decider) -> Result<LocalFrame, Refusal> {
        LocalFrame::budgeted(
            &self.origin,
            [&self.across, &self.along, &self.axis],
            bits,
            dec,
        )
    }
    fn local(&self, point: &Point3) -> [Real; 3] {
        let d = point.minus(&self.origin);
        [d.dot(&self.across), d.dot(&self.along), d.dot(&self.axis)]
    }
    // A point's offset from a cylinder's axis, as a vector.
    fn radial(&self, which: usize, point: &Point3) -> Vector3 {
        let (at, axis) = if which == 0 {
            (&self.origin, self.axis.clone())
        } else {
            (&self.foot, self.narrow_axis())
        };
        let d = point.minus(at);
        d.minus(&axis.scaled(&d.dot(&axis)))
    }
    fn axis_gap(&self, which: usize, point: &Point3) -> Real {
        self.radial(which, point).norm_squared()
    }
    fn side_axis(&self) -> Vector3 {
        self.along
            .scaled(&self.cosine)
            .plus(&self.axis.scaled(&-&self.sine))
    }
    fn narrow_axis(&self) -> Vector3 {
        self.axis
            .scaled(&self.cosine)
            .plus(&self.along.scaled(&self.sine))
    }
    fn place(&self, x: &Real, y: &Real, z: &Real) -> Point3 {
        self.origin
            .plus(&self.across.scaled(x))
            .plus(&self.along.scaled(y))
            .plus(&self.axis.scaled(z))
    }

    /// The law at a rational chart parameter of `[-1, 3]`, exactly: the ball
    /// centre and its contacts on the wide and narrow cylinders.
    pub fn section(&self, u: &Rational) -> Result<[Point3; 3], Refusal> {
        if matches!(self.ring, SaddleRing::Crossing { .. }) {
            return collapsed::section(self, u);
        }
        let chart = chart_of(u)?;
        let s = u - &q(2 * chart);
        let wide = Rational::one() + &s * &s;
        let turn = q(1 - 2 * chart);
        let cos = real(&((Rational::one() - &s * &s) * turn.clone() / wide.clone()));
        let sin = real(&(&s * &q(2) * turn / wide));
        let [big, small] = [0, 1].map(|which| self.offset(which));
        let (x, y, cos, sin) = match self.ring {
            SaddleRing::Through => {
                let x = &self.gap + &(&small * &cos);
                let root = math::exact_sqrt(&(&big.square() - &x.square()));
                (x, &root * &Real::from_i64(self.branch), cos, sin)
            }
            SaddleRing::Overlap => {
                let reach = math::exact_sqrt(&self.gap.square());
                let middle = &(&(&big + &reach) - &small) / &Real::from_i64(2);
                let half = &(&(&big - &reach) + &small) / &Real::from_i64(2);
                let far = &middle + &(&half * &(&cos.square() - &sin.square()));
                let twice = &half * &Real::from_i64(2);
                let up = math::exact_sqrt(&(&twice * &(&big + &far)));
                let round = math::exact_sqrt(&(&twice * &(&(&reach + &small) - &far)));
                let x = &(&far * &self.gap) / &reach;
                let across = &(&x - &self.gap) / &small;
                (x, &sin * &up, across, &(&cos * &round) / &small)
            }
            SaddleRing::Crossing { .. } => unreachable!("handled before the closed-ring charts"),
        };
        let z = &(&(&self.cosine * &y) - &(&small * &sin)) / &self.sine;
        let scale = &self.wide / &big;
        let step = &self.narrow - &small;
        Ok([
            self.place(&x, &y, &z),
            self.place(&(&x * &scale), &(&y * &scale), &z),
            self.place(
                &(&x + &(&step * &cos)),
                &(&y + &(&(&step * &self.cosine) * &sin)),
                &(&z - &(&(&step * &self.sine) * &sin)),
            ),
        ])
    }
}

fn chart_of(u: &Rational) -> Result<i64, Refusal> {
    if *u < q(-1) || *u > q(3) {
        return Err(invalid("chart parameter is outside [-1, 3]"));
    }
    Ok(i64::from(*u > Rational::one()))
}

// Enclosed law constants.
struct Constants {
    gap: Ball,
    cosine: Ball,
    sine: Ball,
    over_sine: Ball,
    wide: Ball,
    narrow: Ball,
    big: Ball,
    small: Ball,
    radius: Ball,
    sides: [Ball; 2],
    branch: Ball,
    overlap: Option<Overlap>,
    one: Ball,
    two: Ball,
}

// An overlap's loop: `X = middle + half cos 2 w`, `facing` the gap's sign and `end`
// the value `|gap| + p` that `X` stays below.
struct Overlap {
    facing: Ball,
    middle: Ball,
    half: Ball,
    end: Ball,
}

impl Constants {
    fn new(def: &RollingSaddleDefinition, bits: u32, dec: &Decider) -> Result<Self, Refusal> {
        let e = |value: &Real| dec.enclose_at(value, bits);
        let one = b(&Rational::one(), bits);
        let sine = e(&def.sine)?;
        let over_sine = one
            .checked_div(&sine)
            .ok_or_else(|| unresolved("saddle axes are not separated from parallel"))?;
        let (big, small, gap) = (e(&def.offset(0))?, e(&def.offset(1))?, e(&def.gap)?);
        let overlap = (def.ring == SaddleRing::Overlap).then(|| {
            let facing = b(&q(if def.gap_positive { 1 } else { -1 }), bits);
            let reach = &gap * &facing;
            let half = b(&Rational::from_ratio(1, 2).expect("nonzero"), bits);
            Overlap {
                middle: &(&(&big + &reach) - &small) * &half,
                half: &(&(&big - &reach) + &small) * &half,
                end: &reach + &small,
                facing,
            }
        });
        Ok(Self {
            gap,
            cosine: e(&def.cosine)?,
            sine,
            over_sine,
            wide: e(&def.wide)?,
            narrow: e(&def.narrow)?,
            big,
            small,
            overlap,
            radius: e(&def.radius)?,
            sides: def.sides.map(|side| b(&q(side), bits)),
            branch: b(&q(def.branch), bits),
            two: b(&q(2), bits),
            one,
        })
    }
    fn c(&self, value: &Ball) -> Jet {
        Jet::constant(value.clone())
    }
}

// The law at one station, as jets in the chart parameter.
struct Rows {
    height: Jet,
    lift: Jet,
    half: Jet,
    turned: Option<Jet>,
    tangent: Jet,
    axial: Jet,
    weight: Option<Jet>,
    cosine: Jet,
    centre: [Jet; 3],
    wide: [Jet; 3],
    narrow: [Jet; 3],
}

// The ring's plan coordinates at a chart station, with `cos t`, `sin t` and the
// half-angle tangents of the rails.
struct Plan {
    x: Jet,
    y: Jet,
    cos: Jet,
    sin: Jet,
    half: Jet,
    turned: Option<Jet>,
}

fn plan(k: &Constants, turn_cos: Jet, turn_sin: Jet) -> Result<Plan, Refusal> {
    let Some(overlap) = &k.overlap else {
        let (cos, sin) = (turn_cos, turn_sin);
        let x = k.c(&k.gap).add(&cos.scale(&k.small));
        let root = k.c(&k.big.square()).sub(&x.square()).sqrt()?;
        let half = x.div(&k.c(&k.big).add(&root))?;
        return Ok(Plan {
            y: root.scale(&k.branch),
            x,
            cos,
            sin,
            half,
            turned: None,
        });
    };
    let wave = turn_cos.square().sub(&turn_sin.square());
    let far = k.c(&overlap.middle).add(&wave.scale(&overlap.half));
    let twice = &overlap.half * &k.two;
    let up = far.add(&k.c(&k.big)).scale(&twice).sqrt()?;
    let round = k.c(&overlap.end).sub(&far).scale(&twice).sqrt()?;
    let over_small = k
        .one
        .checked_div(&k.small)
        .ok_or_else(|| unresolved("narrow offset is not separated from zero"))?;
    let x = far.scale(&overlap.facing);
    Ok(Plan {
        y: turn_sin.mul(&up),
        cos: x.sub(&k.c(&k.gap)).scale(&over_small),
        sin: turn_cos.mul(&round).scale(&over_small),
        half: turn_sin.scale(&twice).div(&up)?,
        turned: Some(turn_cos.scale(&twice).div(&round)?.neg()),
        x,
    })
}

// `fitted` holds the rails' fitted half-angle tangents, the narrow one on a loop.
fn rows(
    k: &Constants,
    chart: i64,
    s: &Jet,
    fitted: Option<(&Jet, Option<&Jet>)>,
) -> Result<Rows, Refusal> {
    let bits = k.one.precision();
    let one = k.c(&k.one);
    let turn = b(&q(1 - 2 * chart), bits);
    let square = s.square();
    let den = one.add(&square);
    let turn_cos = one.sub(&square).scale(&turn).div(&den)?;
    let turn_sin = s.scale(&k.two).scale(&turn).div(&den)?;
    let Plan {
        x,
        y,
        cos,
        sin,
        half,
        turned,
    } = plan(k, turn_cos, turn_sin)?;
    let height = y
        .sub(&sin.scale(&(&k.small * &k.cosine)))
        .scale(&k.over_sine);
    let lift = y
        .scale(&k.cosine)
        .sub(&sin.scale(&k.small))
        .scale(&k.over_sine);
    let over_big = k
        .one
        .checked_div(&k.big)
        .ok_or_else(|| unresolved("wide offset is not separated from zero"))?;
    let (qx, qy) = (x.scale(&over_big), y.scale(&over_big));
    // Unit rays from the centre: -side (qx, qy, 0) to the wide contact and
    // -side (cos, C sin, -S sin) to the narrow one.
    let both = &k.sides[0] * &k.sides[1];
    let meet = qx
        .mul(&cos)
        .add(&qy.mul(&sin).scale(&k.cosine))
        .scale(&both);
    let one_d = one.add(&meet);
    // Apex M = C + r (a + c) / (1 + d). Against the narrow contact it moves
    // r (a - d c) / (1 + d), which has no part along c.
    let toward_tangent = qx
        .mul(&sin)
        .sub(&qy.mul(&cos).scale(&k.cosine))
        .scale(&k.sides[0]);
    let toward_axis = qy.scale(&(&k.sine * &k.sides[0])).neg();
    let tangent = toward_tangent.scale(&k.radius).div(&one_d)?;
    let axial = toward_axis.scale(&k.radius).div(&one_d)?;
    let half_turn = k
        .one
        .checked_div(&k.two)
        .expect("two is separated from zero");
    let cosine = one_d.scale(&half_turn).sqrt()?;
    // The middle weight cos(half) sqrt(w0 w2) over the narrow rail's denominator.
    let weight = match fitted {
        None => None,
        Some((w, narrow)) => {
            let rail = narrow.map_or_else(|| Jet(den.0.clone()), |n| one.add(&n.square()));
            Some(
                one_d
                    .mul(&one.add(&w.square()))
                    .div(&rail.scale(&k.two))?
                    .sqrt()?,
            )
        }
    };
    let step = &k.narrow - &k.small;
    let scale = k
        .wide
        .checked_div(&k.big)
        .ok_or_else(|| unresolved("wide offset is not separated from zero"))?;
    Ok(Rows {
        centre: [x.clone(), y.clone(), lift.clone()],
        wide: [x.scale(&scale), y.scale(&scale), lift.clone()],
        narrow: [
            x.add(&cos.scale(&step)),
            y.add(&sin.scale(&(&step * &k.cosine))),
            lift.sub(&sin.scale(&(&step * &k.sine))),
        ],
        height,
        lift,
        half,
        turned,
        tangent,
        axial,
        weight,
        cosine,
    })
}

/// Hermite data of the fitted rows at one chart node.
#[derive(Clone, Debug)]
struct Node {
    u: Rational,
    values: [Rational; ROWS],
    slopes: [Rational; ROWS],
    // The weight's slope on the span ending here. The rail denominator turns
    // round at a chart end, so the weight has a corner there and only there.
    arriving: Rational,
}

// The fitted rows' data at one chart station. Its weight has no corner, which
// holds at every node but a through ring's chart ends.
fn node_at(k: &Constants, chart: i64, s: &Rational, bits: u32) -> Result<Node, Refusal> {
    let point = Jet::linear(b(s, bits), k.one.clone());
    let exact = rows(k, chart, &point, None)?;
    let mut values: [Rational; ROWS] = std::array::from_fn(|_| Rational::zero());
    let mut slopes = values.clone();
    let plain = [
        &exact.height,
        &exact.lift,
        &exact.half,
        &exact.tangent,
        &exact.axial,
    ];
    for (j, jet) in plain.into_iter().enumerate() {
        let [value, slope] = jet.data();
        values[j] = snap(&value.mid());
        slopes[j] = snap(&slope.mid());
    }
    if let Some(jet) = &exact.turned {
        let [value, slope] = jet.data();
        values[TURNED] = snap(&value.mid());
        slopes[TURNED] = snap(&slope.mid());
    }
    let fitted = Jet::linear(b(&values[HALF], bits), b(&slopes[HALF], bits));
    let turned = Jet::linear(b(&values[TURNED], bits), b(&slopes[TURNED], bits));
    let narrow = k.overlap.is_some().then_some(&turned);
    let weight = rows(k, chart, &point, Some((&fitted, narrow)))?.weight;
    let [value, slope] = weight.expect("weight requested").data();
    values[WEIGHT] = snap(&value.mid());
    slopes[WEIGHT] = snap(&slope.mid());
    Ok(Node {
        u: s + &q(2 * chart),
        arriving: slopes[WEIGHT].clone(),
        values,
        slopes,
    })
}

fn produce(
    def: &RollingSaddleDefinition,
    spans: usize,
    bits: u32,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    let k = Constants::new(def, bits, dec)?;
    let mut nodes = Vec::with_capacity(2 * spans + 1);
    for chart in 0..2 {
        for i in 0..spans {
            dec.check_deadline()?;
            let s =
                Rational::from_ratio(2 * i as i64 - spans as i64, spans as i64).expect("nonzero");
            let mut node = node_at(&k, chart, &s, bits)?;
            // A through ring's rail denominator turns round at a chart end.
            if i == 0 && chart == 1 && k.overlap.is_none() {
                let fitted = Jet::linear(b(&node.values[HALF], bits), b(&node.slopes[HALF], bits));
                let end = Jet::linear(k.one.clone(), k.one.clone());
                let weight = rows(&k, 0, &end, Some((&fitted, None)))?.weight;
                node.arriving = snap(&weight.expect("weight requested").data()[1].mid());
            }
            nodes.push(node);
        }
    }
    let first = nodes[0].clone();
    let arriving = if k.overlap.is_some() {
        first.slopes[WEIGHT].clone()
    } else {
        let fitted = Jet::linear(b(&first.values[HALF], bits), b(&first.slopes[HALF], bits));
        let end = Jet::linear(k.one.clone(), k.one.clone());
        let weight = rows(&k, 1, &end, Some((&fitted, None)))?.weight;
        snap(&weight.expect("weight requested").data()[1].mid())
    };
    nodes.push(Node {
        u: q(3),
        arriving,
        ..first
    });
    Ok(nodes)
}

// The chart of a span and its ends in that chart's parameter.
fn span_chart(pair: &[Node]) -> (i64, [Rational; 2]) {
    let chart = i64::from(pair[0].u >= Rational::one());
    let shift = q(2 * chart);
    (chart, [&pair[0].u - &shift, &pair[1].u - &shift])
}

type RowData = [[Ball; 2]; ROWS];

// The same bound as a rolling envelope's. Both surfaces are rational quadratics
// in v with shared rail weights, so with `eps` the relative middle-weight error
//   |S - E| <= max_i |s_i - e_i| + (eps / (2 (1 - eps))) r / cos(half)
// and each row error is the cubic Hermite remainder plus the node data's miss.
fn recompute_bound(
    def: &RollingSaddleDefinition,
    nodes: &[Node],
    bits: u32,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let k = Constants::new(def, bits, dec)?;
    let mut worst = Rational::zero();
    for pair in nodes.windows(2) {
        dec.check_deadline()?;
        worst = worst.max(span_bound(&k, pair, bits)?);
    }
    Ok(worst)
}

fn span_bound(k: &Constants, pair: &[Node], bits: u32) -> Result<Rational, Refusal> {
    let enclosed = |node: &Node, ending: bool| -> RowData {
        std::array::from_fn(|j| {
            let slope = if ending && j == WEIGHT {
                &node.arriving
            } else {
                &node.slopes[j]
            };
            [b(&node.values[j], bits), b(slope, bits)]
        })
    };
    let exact_at = |chart: i64, s: &Rational, data: &RowData| -> Result<RowData, Refusal> {
        let point = Jet::linear(b(s, bits), k.one.clone());
        let fitted = Jet::linear(data[HALF][0].clone(), data[HALF][1].clone());
        let turned = Jet::linear(data[TURNED][0].clone(), data[TURNED][1].clone());
        let narrow = k.overlap.is_some().then_some(&turned);
        let exact = rows(k, chart, &point, Some((&fitted, narrow)))?;
        let zero = [b(&Rational::zero(), bits), b(&Rational::zero(), bits)];
        Ok([
            exact.height.data(),
            exact.lift.data(),
            exact.half.data(),
            exact.tangent.data(),
            exact.axial.data(),
            exact.weight.expect("weight requested").data(),
            exact.turned.as_ref().map_or(zero, Jet::data),
        ])
    };
    let (chart, ends) = span_chart(pair);
    let span = [&ends[0], &ends[1]];
    let left = enclosed(&pair[0], false);
    let right = enclosed(&pair[1], true);
    let left_exact = exact_at(chart, &ends[0], &left)?;
    let right_exact = exact_at(chart, &ends[1], &right)?;
    let h = &ends[1] - &ends[0];
    let s = Jet::linear(b(&ends[0], bits).hull(&b(&ends[1], bits)), k.one.clone());
    let fitted = hermite(
        &s,
        span,
        [
            [&left[HALF][0], &left[HALF][1]],
            [&right[HALF][0], &right[HALF][1]],
        ],
        bits,
    )?;
    let turned = hermite(
        &s,
        span,
        [
            [&left[TURNED][0], &left[TURNED][1]],
            [&right[TURNED][0], &right[TURNED][1]],
        ],
        bits,
    )?;
    let narrow_half = k.overlap.is_some().then_some(&turned);
    let whole = rows(k, chart, &s, Some((&fitted, narrow_half)))?;
    let weight_jet = whole.weight.as_ref().expect("weight requested");
    let nothing = b(&Rational::zero(), bits);
    let fourth = [
        &whole.height.0[4],
        &whole.lift.0[4],
        &whole.half.0[4],
        &whole.tangent.0[4],
        &whole.axial.0[4],
        &weight_jet.0[4],
        whole.turned.as_ref().map_or(&nothing, |jet| &jet.0[4]),
    ];
    let h4 = h.pow(4).expect("nonzero span");
    let error = |j: usize| -> Rational {
        let remainder = magnitude(fourth[j]) * h4.clone() / q(16);
        let values = magnitude(&(&left[j][0] - &left_exact[j][0]))
            .max(magnitude(&(&right[j][0] - &right_exact[j][0])));
        let slopes = magnitude(&(&left[j][1] - &left_exact[j][1]))
            + magnitude(&(&right[j][1] - &right_exact[j][1]));
        remainder + values + q(4) * h.clone() * slopes / q(27)
    };
    // An angle moves at most 2 per unit of its half-angle tangent.
    let wide = q(2) * k.wide.hi() * error(HALF) + error(LIFT);
    // On a loop the narrow rail's own angle is fitted too, and it swings the
    // apex about the axis by its distance from there.
    let swing = if k.overlap.is_some() {
        q(2) * error(TURNED)
    } else {
        Rational::zero()
    };
    let narrow = error(HEIGHT) + k.narrow.hi() * swing.clone();
    let apex = error(TANGENT)
        + error(AXIAL)
        + error(HEIGHT)
        + (k.narrow.hi() + magnitude(&whole.tangent.0[0])) * swing;
    let weight = weight_jet.0[0].lo();
    if weight <= Rational::zero() {
        return Err(unresolved("exact middle weight is not separated from zero"));
    }
    let eps = error(WEIGHT) / weight;
    if eps >= Rational::one() {
        return Err(invalid("middle weight error is not below its exact value"));
    }
    let cosine = whole.cosine.0[0].lo();
    if cosine <= Rational::zero() {
        return Err(unresolved(
            "contact rays are not separated from a half turn",
        ));
    }
    let spread = eps.clone() / (q(2) * (Rational::one() - eps)) * k.radius.hi() / cosine;
    Ok(wide.max(narrow).max(apex) + spread)
}

// Fits a loop: each span over the target halves, so a hard stretch costs only its
// own spans. A span the bound cannot resolve halves too, until it runs out of
// halvings. The nodes come out in chart order, after `pair`'s first.
fn settle(
    k: &Constants,
    pair: [Node; 2],
    target: &Rational,
    depth: u32,
    bits: u32,
    dec: &mut Decider,
    nodes: &mut Vec<Node>,
) -> Result<Rational, Refusal> {
    dec.check_deadline()?;
    match span_bound(k, &pair, bits) {
        Ok(upper) if upper <= *target || depth >= MAX_SPLITS => {
            let [_, end] = pair;
            nodes.push(end);
            return Ok(upper);
        }
        Err(error) if depth >= MAX_SPLITS => return Err(error),
        _ => {}
    }
    let (chart, ends) = span_chart(&pair);
    let middle = node_at(k, chart, &((&ends[0] + &ends[1]) / q(2)), bits)?;
    let [start, end] = pair;
    let left = settle(
        k,
        [start, middle.clone()],
        target,
        depth + 1,
        bits,
        dec,
        nodes,
    )?;
    let right = settle(k, [middle, end], target, depth + 1, bits, dec, nodes)?;
    Ok(left.max(right))
}

fn refine(
    def: &RollingSaddleDefinition,
    target: &Rational,
    bits: u32,
    dec: &mut Decider,
) -> Result<(Vec<Node>, Rational), Refusal> {
    let k = Constants::new(def, bits, dec)?;
    let first = produce(def, FIRST_SPANS, bits, dec)?;
    let mut nodes = vec![first[0].clone()];
    let mut worst = Rational::zero();
    for pair in first.windows(2) {
        let pair = [pair[0].clone(), pair[1].clone()];
        worst = worst.max(settle(&k, pair, target, 0, bits, dec, &mut nodes)?);
    }
    Ok((nodes, worst))
}

// Proves the exact surface regular. With rail changes `Ai = r ei' . c'`, both
// rays square to `c'`, and `c` the cosine of half the section angle, the ray `h`
// from the section's middle gives
//   S_u . c' = |c'|^2 + (A0 + A1) cos h / (2 c) + (A1 - A0) sin h / (2 s)
// where `cos h` lies in `[c, 1]` and `|sin h| <= s`. `S_v` is square to `c'`, so a
// positive result makes the two derivatives independent.
fn check_regular(
    def: &RollingSaddleDefinition,
    nodes: &[Node],
    bits: u32,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let k = Constants::new(def, bits, dec)?;
    let zero = Rational::zero();
    let first = |jets: &[Jet; 3]| -> [Ball; 3] { std::array::from_fn(|i| jets[i].0[1].clone()) };
    let dot =
        |a: &[Ball; 3], c: &[Ball; 3]| &(&(&a[0] * &c[0]) + &(&a[1] * &c[1])) + &(&a[2] * &c[2]);
    for pair in nodes.windows(2) {
        let (chart, ends) = span_chart(pair);
        let mut cells = 1;
        loop {
            dec.check_deadline()?;
            let width = (&ends[1] - &ends[0]) / q(cells as i64);
            let proved = (0..cells).try_fold(true, |proved, cell| -> Result<bool, Refusal> {
                if !proved {
                    return Ok(false);
                }
                let low = &ends[0] + &(&width * &q(cell as i64));
                let high = &low + &width;
                let s = Jet::linear(b(&low, bits).hull(&b(&high, bits)), k.one.clone());
                let at = rows(&k, chart, &s, None)?;
                let path = first(&at.centre);
                let speed = dot(&path, &path);
                let cosine = &at.cosine.0[0];
                if cosine.hi() >= Rational::one() || cosine.lo() <= zero || speed.lo() <= zero {
                    return Ok(false);
                }
                let [wide, narrow] =
                    [&at.wide, &at.narrow].map(|rail| &dot(&first(rail), &path) - &speed);
                let half = b(&Rational::from_ratio(1, 2).expect("nonzero"), bits);
                let middle = &(&wide + &narrow) * &half;
                let over = middle
                    .checked_div(cosine)
                    .ok_or_else(|| unresolved("section angle is not resolved"))?;
                let spread = magnitude(&(&(&narrow - &wide) * &half));
                Ok(speed.lo() + middle.lo().min(over.lo()) - spread > zero)
            })?;
            if proved {
                break;
            }
            cells *= 2;
            if cells > MAX_CELLS {
                return Err(unresolved(
                    "the ball's surface is not proved regular on a span",
                ));
            }
        }
    }
    Ok(())
}

fn scaled(poly: &[Real], factor: &Real) -> Vec<Real> {
    poly.iter().map(|value| value * factor).collect()
}

// Local homogeneous Bernstein rows of one span, [wide, middle, narrow], the wide
// rail in the wide frame and the others in the narrow one.
fn span_rows(def: &RollingSaddleDefinition, pair: &[Node]) -> [[Vec<Real>; 4]; 3] {
    let (chart, ends) = span_chart(pair);
    let span = [&ends[0], &ends[1]];
    let one = Real::one();
    let minus = -&one;
    let f = |j: usize| {
        let ending = if j == WEIGHT {
            &pair[1].arriving
        } else {
            &pair[1].slopes[j]
        };
        cubic(
            span,
            [&real(&pair[0].values[j]), &real(&pair[1].values[j])],
            [&pair[0].slopes[j], ending],
        )
    };
    let (cos, sin, den) = if def.ring == SaddleRing::Overlap {
        let turned = f(TURNED);
        let square = product(&turned, &turned);
        let ones = vec![one.clone(); 7];
        let facing = Real::from_i64(if def.gap_positive { 1 } else { -1 });
        (
            combine(&ones, &square, &-&facing, &facing),
            scaled(&elevate(&turned, 6), &Real::from_i64(-2)),
            combine(&ones, &square, &one, &one),
        )
    } else {
        let s = vec![real(&ends[0]), real(&ends[1])];
        let square = product(&s, &s);
        let turn = Real::from_i64(1 - 2 * chart);
        let ones = vec![one.clone(); 3];
        (
            combine(&ones, &square, &turn, &-&turn),
            scaled(&elevate(&s, 2), &(&turn * &Real::from_i64(2))),
            combine(&ones, &square, &one, &one),
        )
    };
    let rail_x = scaled(&cos, &def.narrow);
    let rail_e = scaled(&sin, &def.narrow);
    let narrow = [
        rail_x.clone(),
        rail_e.clone(),
        product(&f(HEIGHT), &den),
        den.clone(),
    ];
    let tangent = f(TANGENT);
    let x = combine(&rail_x, &product(&tangent, &sin), &one, &minus);
    let e = combine(&rail_e, &product(&tangent, &cos), &one, &one);
    let height = combine(&f(HEIGHT), &f(AXIAL), &one, &one);
    let weight = f(WEIGHT);
    let middle = [
        product(&weight, &x),
        product(&weight, &e),
        product(&weight, &product(&height, &den)),
        product(&weight, &den),
    ];
    let half = f(HALF);
    let half_square = product(&half, &half);
    let ones = vec![one.clone(); 7];
    let wide_den = combine(&ones, &half_square, &one, &one);
    let twice = scaled(&elevate(&half, 6), &(&def.wide * &Real::from_i64(2)));
    let flat = combine(&ones, &half_square, &one, &minus);
    let lifted = product(&f(LIFT), &wide_den);
    // A through ring's angle runs from `along`, and a loop's from the gap's side.
    let wide = if def.ring == SaddleRing::Overlap {
        let facing = Real::from_i64(if def.gap_positive { 1 } else { -1 });
        [
            scaled(&flat, &(&def.wide * &facing)),
            twice,
            lifted,
            wide_den,
        ]
    } else {
        let branch = Real::from_i64(def.branch);
        [
            twice,
            scaled(&flat, &(&def.wide * &branch)),
            lifted,
            wide_den,
        ]
    };
    let degree = def.degree();
    [wide, middle, narrow].map(|row| row.map(|c| elevate(&c, degree)))
}

fn build(
    def: &RollingSaddleDefinition,
    nodes: &[Node],
    dec: &mut Decider,
) -> Result<(BSplineSurface, Vec<[Real; 4]>), Refusal> {
    let foot = &def.foot;
    let frame = [def.across.clone(), def.side_axis(), def.narrow_axis()];
    let place = |row: usize, c: [&Real; 4]| -> Point3 {
        let [x, y, z] = [c[0] / c[3], c[1] / c[3], c[2] / c[3]];
        if row == 0 {
            def.place(&x, &y, &z)
        } else {
            foot.plus(&frame[0].scaled(&x))
                .plus(&frame[1].scaled(&y))
                .plus(&frame[2].scaled(&z))
        }
    };
    let spans = nodes.len() - 1;
    let degree = def.degree();
    let count = spans * degree + 1;
    let mut controls = Vec::with_capacity(count * 3);
    let mut weights = Vec::with_capacity(count * 3);
    let mut homogeneous = Vec::with_capacity(count * 3);
    let mut knots = vec![q(-1); degree + 1];
    for (index, pair) in nodes.windows(2).enumerate() {
        dec.check_deadline()?;
        let rows = span_rows(def, pair);
        let last = index + 1 == spans;
        for i in 0..degree + usize::from(last) {
            for (row, data) in rows.iter().enumerate() {
                controls.push(place(
                    row,
                    [&data[0][i], &data[1][i], &data[2][i], &data[3][i]],
                ));
                weights.push(data[3][i].clone());
                let local = std::array::from_fn(|coordinate| data[coordinate][i].clone());
                // The other rows use the narrow frame, whose origin is
                // origin + gap across. Rotate their homogeneous coordinates
                // into the wide frame without dividing by the weight.
                homogeneous.push(if row == 0 {
                    local
                } else {
                    let [x, y, z, w] = local;
                    [
                        &x + &(&def.gap * &w),
                        &(&def.cosine * &y) + &(&def.sine * &z),
                        &(&def.cosine * &z) - &(&def.sine * &y),
                        w,
                    ]
                });
            }
        }
        let repeat = degree + usize::from(last);
        knots.extend(std::iter::repeat_n(pair[1].u.clone(), repeat));
    }
    let v = [0, 0, 0, 1, 1, 1].map(q).to_vec();
    let carrier = BSplineSurface::new([degree, 2], [knots, v], [count, 3], controls, weights, dec)?;
    Ok((carrier, homogeneous))
}

// Enclosed homogeneous controls in the wide frame, u-major. Rigid copies share
// one net, so no consumer pays to enclose a moved carrier again.
#[derive(Debug)]
struct Net {
    bits: u32,
    controls: Vec<[Ball; 4]>,
    grid: Mutex<Option<(Rational, Arc<grid::LocalGrid>)>>,
    turns: Mutex<Option<Vec<[Rational; 3]>>>,
}

impl Net {
    fn from_homogeneous(controls: &[[Real; 4]], dec: &Decider) -> Result<Self, Refusal> {
        let bits = BITS.min(dec.policy().max_precision_bits);
        let values: Vec<_> = controls.iter().flat_map(|point| point.iter()).collect();
        let enclosed = dec.enclose_all_at(&values, bits)?;
        let controls = enclosed
            .as_chunks::<4>()
            .0
            .iter()
            .map(|point| std::array::from_fn(|coordinate| point[coordinate].clone()))
            .collect();
        Ok(Self {
            bits,
            controls,
            grid: Mutex::new(None),
            turns: Mutex::new(None),
        })
    }

    fn new(
        def: &RollingSaddleDefinition,
        carrier: &BSplineSurface,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let bits = BITS.min(dec.policy().max_precision_bits);
        let frame = def.frame(bits, dec)?;
        let controls = carrier
            .controls()
            .iter()
            .zip(carrier.weights())
            .map(|(point, weight)| {
                let w = dec.enclose_at(weight, bits)?;
                let [x, y, z] = frame.local(&math::point(point, bits, dec)?);
                Ok([&x * &w, &y * &w, &z * &w, w])
            })
            .collect::<Result<Vec<_>, Refusal>>()?;
        Ok(Self {
            bits,
            controls,
            grid: Mutex::new(None),
            turns: Mutex::new(None),
        })
    }
}

// Proves the carrier itself regular, patch by patch. With P the homogeneous
// point less `m W` for a near constant `m`, the derivatives over a patch are
// `D = P_x W - P W_x` up to the positive `W^2`, and a positive `(D_u x D_v) . n`
// for a fixed `n` keeps them independent.
fn check_carrier(net: &Net, carrier: &BSplineSurface, dec: &Decider) -> Result<(), Refusal> {
    let bits = net.bits;
    let degree = carrier.degrees()[0];
    let zero = b(&Rational::zero(), bits);
    let hull = |values: &mut dyn Iterator<Item = Ball>| {
        values.fold(None::<Ball>, |all, value| {
            Some(all.map_or(value.clone(), |all| all.hull(&value)))
        })
    };
    let mut stack: Vec<(Patch, u32)> = Patch::all(carrier, &net.controls, degree)
        .into_iter()
        .map(|patch| (patch, 0))
        .collect();
    while let Some((patch, depth)) = stack.pop() {
        dec.check_deadline()?;
        let [du, dv] = [0, 1].map(|axis| &patch.cell[axis][1] - &patch.cell[axis][0]);
        let first = &patch.net[0];
        let middle = [0, 1, 2].map(|k| {
            let at = first[k].checked_div(&first[3]);
            at.map(|value| b(&value.mid(), bits))
        });
        let [Some(mx), Some(my), Some(mz)] = middle else {
            return Err(unresolved("subdivided saddle weight is not positive"));
        };
        let middle = [mx, my, mz];
        let shifted: Vec<[Ball; 4]> = patch
            .net
            .iter()
            .map(|h| {
                [
                    &h[0] - &(&h[3] * &middle[0]),
                    &h[1] - &(&h[3] * &middle[1]),
                    &h[2] - &(&h[3] * &middle[2]),
                    h[3].clone(),
                ]
            })
            .collect();
        let cells = |step: [usize; 2], scale: &Rational, k: usize| {
            let mut all = Vec::new();
            for i in 0..=degree - step[0] {
                for j in 0..=2 - step[1] {
                    let near = &shifted[i * 3 + j];
                    let far = &shifted[(i + step[0]) * 3 + j + step[1]];
                    all.push(&(&far[k] - &near[k]) * &b(scale, bits));
                }
            }
            hull(&mut all.into_iter()).expect("a patch has controls")
        };
        let value = |k: usize| {
            hull(&mut shifted.iter().map(|h| h[k].clone())).expect("a patch has controls")
        };
        let scales = [q(degree as i64) / du.clone(), q(2) / dv.clone()];
        let [weight, along, across] = [
            value(3),
            cells([1, 0], &scales[0], 3),
            cells([0, 1], &scales[1], 3),
        ];
        let derivative = |step: [usize; 2], scale: &Rational, change: &Ball| -> [Ball; 3] {
            [0, 1, 2].map(|k| &(&cells(step, scale, k) * &weight) - &(&value(k) * change))
        };
        let u = derivative([1, 0], &scales[0], &along);
        let v = derivative([0, 1], &scales[1], &across);
        let cross = |i: usize, j: usize| &(&u[i] * &v[j]) - &(&u[j] * &v[i]);
        let normal = [cross(1, 2), cross(2, 0), cross(0, 1)];
        let direction = normal.each_ref().map(|value| b(&value.mid(), bits));
        let along_normal = normal
            .iter()
            .zip(&direction)
            .fold(zero.clone(), |sum, (value, unit)| sum + value * unit);
        if along_normal.lo() > Rational::zero() {
            continue;
        }
        if depth >= 12 {
            return Err(unresolved("the fitted carrier is not proved regular"));
        }
        // Compare in the same metric as a patch split, across the longer side.
        let axis = usize::from(dv / q(8) > du);
        let half = (&patch.cell[axis][0] + &patch.cell[axis][1]) / q(2);
        let [left, right] = split_net(&patch.net, [degree, 2], axis);
        let (mut low, mut high) = (patch.cell.clone(), patch.cell.clone());
        low[axis][1] = half.clone();
        high[axis][0] = half;
        for (cell, net) in [(low, left), (high, right)] {
            stack.push((Patch { degree, cell, net }, depth + 1));
        }
    }
    Ok(())
}

/// A rolling-ball band over `[-1, 3] x [0, 1]`, from the wide rail to the
/// narrow one. Regular rings close in `u`; an equal-cylinder crossing band
/// instead collapses each `u` end to one contact.
#[derive(Clone, Debug)]
pub struct RollingSaddle {
    data: Arc<Data>,
}

#[derive(Debug)]
struct Data {
    definition: RollingSaddleDefinition,
    construction: Construction,
    carrier: BSplineSurface,
    cached_rails: Option<Arc<[BSplineCurve; 2]>>,
    net: Arc<Net>,
    bound: EvidenceBound,
    reversing: bool,
}

#[derive(Clone, Debug)]
enum Construction {
    Fitted(Vec<Node>),
    Crossing(Arc<collapsed::Exact>),
}

impl RollingSaddle {
    /// Fits spans whose recomputed bound meets the 100 nm target: the fewest
    /// uniform spans, by doubling, on a through ring, and spans halved where the
    /// bound needs them on a loop. A fit that meets only the 1 um cap is accepted
    /// and reported by its bound, and anything larger refuses.
    /// Equal crossing bands reduce exactly to two degree-six spans instead,
    /// with a zero bound proved over their complete collapsed chart.
    pub fn fit(definition: RollingSaddleDefinition, dec: &mut Decider) -> Result<Self, Refusal> {
        if matches!(definition.ring, SaddleRing::Crossing { .. }) {
            return collapsed::build(definition, dec);
        }
        let bits = BITS.min(dec.policy().max_precision_bits);
        let target = Rational::from_ratio(ROLLING_SADDLE_TARGET_NM, 1_000_000).expect("nonzero");
        let cap = Rational::from_ratio(HARD_CAP_NM, 1_000_000).expect("nonzero");
        let (nodes, upper) = if definition.ring == SaddleRing::Overlap {
            refine(&definition, &target, bits, dec)?
        } else {
            let mut spans = FIRST_SPANS;
            loop {
                let nodes = produce(&definition, spans, bits, dec)?;
                match recompute_bound(&definition, &nodes, bits, dec) {
                    Ok(upper) if upper <= target || spans >= LAST_SPANS => break (nodes, upper),
                    Err(error) if spans >= LAST_SPANS => return Err(error),
                    _ => spans *= 2,
                }
            }
        };
        if upper > cap {
            return Err(invalid("the saddle fit misses the 1 um approximation cap"));
        }
        check_regular(&definition, &nodes, bits, dec)?;
        let (carrier, homogeneous) = build(&definition, &nodes, dec)?;
        let net = Arc::new(Net::from_homogeneous(&homogeneous, dec)?);
        check_carrier(&net, &carrier, dec)?;
        Ok(Self {
            data: Arc::new(Data {
                definition,
                construction: Construction::Fitted(nodes),
                carrier,
                cached_rails: None,
                net,
                bound: EvidenceBound {
                    upper,
                    unit: EvidenceUnit::Millimetres,
                },
                reversing: false,
            }),
        })
    }

    /// The represented geometry.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }

    /// The exact law.
    pub fn definition(&self) -> &RollingSaddleDefinition {
        &self.data.definition
    }

    /// Two-sided distance bound between the carrier and the exact surface.
    pub fn bound(&self) -> &EvidenceBound {
        &self.data.bound
    }

    /// Whether this exact crossing band rounds the given whole source ellipse.
    /// Other saddle families have non-elliptic contact curves and return false.
    pub fn is_edge(&self, edge: &crate::Ellipse, dec: &mut Decider) -> Result<bool, Refusal> {
        if !matches!(self.definition().ring, SaddleRing::Crossing { .. }) {
            return Ok(false);
        }
        collapsed::same_ellipse(self.definition(), &self.definition().wide, edge, dec)
    }

    /// Whether this exact crossing band's centres lie on the given ellipse.
    /// Its strict reach proof makes the rolling radius the distance to that
    /// whole ellipse, not merely to one corresponding centre.
    pub fn centres_on(&self, path: &crate::Ellipse, dec: &mut Decider) -> Result<bool, Refusal> {
        if !matches!(self.definition().ring, SaddleRing::Crossing { .. }) {
            return Ok(false);
        }
        collapsed::same_ellipse(self.definition(), &self.definition().offset(0), path, dec)
    }

    /// Number of Bezier spans along the band.
    pub fn span_count(&self) -> usize {
        match &self.data.construction {
            Construction::Fitted(nodes) => nodes.len() - 1,
            Construction::Crossing(_) => 2,
        }
    }

    fn knot(&self, index: usize) -> &Rational {
        match &self.data.construction {
            Construction::Fitted(nodes) => &nodes[index].u,
            Construction::Crossing(exact) => &exact.knots[index],
        }
    }

    /// The widest and narrowest spans, in the chart's `u`.
    pub fn span_widths(&self) -> [Rational; 2] {
        let mut widths = (0..self.span_count()).map(|i| self.knot(i + 1) - self.knot(i));
        let first = widths.next().expect("a ring has spans");
        widths.fold([first.clone(), first], |[wide, narrow], width| {
            [wide.max(width.clone()), narrow.min(width)]
        })
    }

    /// The rails on the wide and narrow cylinders, the carrier's own boundary
    /// rows, so a neighbour can share them exactly.
    pub fn rails(&self) -> [BSplineCurve; 2] {
        if let Some(rails) = &self.data.cached_rails {
            return [rails[0].clone(), rails[1].clone()];
        }
        [SplineBoundary::VMin, SplineBoundary::VMax]
            .map(|side| BSplineCurve::boundary(&self.data.carrier, side))
    }

    /// Positive identity of the shared law and carrier, for cache reuse.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data.net, &other.data.net)
    }

    fn query_bits(&self, dec: &Decider) -> Result<u32, Refusal> {
        match self.data.net.bits.min(dec.policy().max_precision_bits) {
            0 => Err(Refusal::new(
                "numbers.precision-limit",
                "no saddle query precision available",
            )),
            bits => Ok(bits),
        }
    }

    fn patches(&self, bits: u32) -> Vec<Patch> {
        Patch::at(
            &self.data.carrier,
            &self.data.net.controls,
            self.data.net.bits,
            self.data.carrier.degrees()[0],
            bits,
        )
    }

    /// Enclosed world controls and weights, from the shared net.
    pub(crate) fn query_controls(
        &self,
        dec: &mut Decider,
    ) -> Result<Vec<([Ball; 3], Ball)>, Refusal> {
        let bits = self.query_bits(dec)?;
        let frame = self.data.definition.budgeted_frame(bits, dec)?;
        self.data
            .net
            .controls
            .iter()
            .map(|h| {
                let local = [0, 1, 2].map(|k| {
                    h[k].with_precision(bits)
                        .checked_div(&h[3].with_precision(bits))
                        .ok_or_else(|| unresolved("saddle weight is not positive"))
                });
                let [x, y, z] = local;
                Ok((frame.world(&[x?, y?, z?]), h[3].clone()))
            })
            .collect()
    }

    /// Enclosed world controls and weights, for a writer that rounds them.
    pub fn enclosed_controls(&self, dec: &mut Decider) -> Result<Vec<([Ball; 3], Ball)>, Refusal> {
        self.query_controls(dec)
    }

    /// Encloses every chart preimage of a spatial box.
    pub fn parameters_of_box(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        dec.check_deadline()?;
        let bits = self.query_bits(dec)?;
        let ceiling = dec.policy().max_precision_bits;
        let point = point
            .each_ref()
            .map(|value| value.with_precision(value.precision().min(ceiling)));
        let point = self.data.definition.frame(bits, dec)?.local(&point);
        let found = cover(self.patches(bits), &point, dec)?;
        Ok(found.map(|h| h.map(|[lo, hi]| b(&lo, bits).hull(&b(&hi, bits)))))
    }

    /// Whether a point is exactly on a rail. `Ok(false)` says only that.
    pub fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        math::validate_point(point, dec)?;
        Ok(self.rail_of(point, dec)?.is_some())
    }

    fn rail_of(&self, point: &Point3, dec: &mut Decider) -> Result<Option<usize>, Refusal> {
        let def = &self.data.definition;
        let rails = self.rails();
        for (which, rail) in rails.iter().enumerate() {
            let gap = def.axis_gap(which, point);
            let radius = if which == 0 { &def.wide } else { &def.narrow };
            if math::sign(&(&gap - &radius.square()), dec)? == Sign::Zero
                && rail.parameter_of(point, dec)?.is_some()
            {
                return Ok(Some(which));
            }
        }
        Ok(None)
    }

    /// Exact normal at a rail point, pointing to the ball's centre. Interior
    /// points have no exact chart parameters and refuse.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        let def = &self.data.definition;
        let Some(which) = self.rail_of(point, dec)? else {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "a saddle normal needs exact chart parameters, and only its rails have them",
            ));
        };
        let radial = def.radial(which, point);
        let normal = radial.scaled(&Real::from_i64(def.sides[which]));
        Ok(if self.data.reversing {
            normal.scaled(&Real::from_i64(-1))
        } else {
            normal
        })
    }

    /// Whether a point is on the carrier. A nonzero sign proves it off, and
    /// names no global side of the finite support. Exact zero is decided on the
    /// rails only, so an interior point within the enclosure's width refuses.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        if self.on_boundary(point, dec)? {
            return Ok(Sign::Zero);
        }
        let bits = self.query_bits(dec)?;
        let at = math::budgeted_point(point, bits, dec)?;
        if self.parameters_of_box(&at, dec)?.is_none() {
            return Ok(Sign::Positive);
        }
        Err(unresolved(
            "a point is within the carrier's enclosure width",
        ))
    }

    /// Whether every control of the carrier lies strictly on one side of a
    /// plane, which keeps the whole band off that plane.
    pub fn clear_of(&self, plane: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        let mut side = None;
        for control in self.data.carrier.controls() {
            let value = math::sign(&control.minus(&plane.origin).dot(&plane.normal), dec)?;
            if value == Sign::Zero || side.replace(value).is_some_and(|last| last != value) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Bounds of the carrier along a world coordinate, from its control hull.
    pub fn extent(&self, axis: usize, dec: &mut Decider) -> Result<[Rational; 2], Refusal> {
        let mut range: Option<[Rational; 2]> = None;
        for (control, _) in self.query_controls(dec)? {
            let [lo, hi] = [control[axis].lo(), control[axis].hi()];
            range = Some(match range {
                None => [lo, hi],
                Some([low, high]) => [low.min(lo), high.max(hi)],
            });
        }
        range.ok_or_else(|| unresolved("a saddle has controls"))
    }

    /// The smallest and largest of a world coordinate on the carrier, each as a
    /// pair of bounds no more than `width` apart. Positive weights keep each
    /// patch in its control hull, and its corners are points of the carrier.
    pub fn range(
        &self,
        axis: usize,
        width: &Rational,
        dec: &Decider,
    ) -> Result<[[Rational; 2]; 2], Refusal> {
        let bits = self.query_bits(dec)?;
        let frame = self.data.definition.frame(bits, dec)?;
        let mut answer = Vec::with_capacity(2);
        for lowest in [true, false] {
            // Search for the largest of `s * coordinate`, with `s` set by the end.
            // Each patch gives its hull's top, and the best proved of its corners.
            let measure = |patch: &Patch| -> Result<(Rational, Rational), Refusal> {
                let (mut top, mut best) = (None::<Rational>, None::<Rational>);
                let last = patch.net.len() - 1;
                for (index, h) in patch.net.iter().enumerate() {
                    let local = [0, 1, 2].map(|k| h[k].checked_div(&h[3]));
                    let [Some(x), Some(y), Some(z)] = local else {
                        return Err(unresolved("subdivided saddle weight is not positive"));
                    };
                    let world = frame.world(&[x, y, z]);
                    let (lo, hi) = if lowest {
                        (-world[axis].hi(), -world[axis].lo())
                    } else {
                        (world[axis].lo(), world[axis].hi())
                    };
                    top = Some(top.map_or(hi.clone(), |value| value.max(hi)));
                    if index == 0 || index == 2 || index == last - 2 || index == last {
                        best = Some(best.map_or(lo.clone(), |value| value.max(lo)));
                    }
                }
                top.zip(best)
                    .ok_or_else(|| unresolved("a patch has controls"))
            };
            let mut open = Vec::new();
            let mut floor = None::<Rational>;
            for patch in self.patches(bits) {
                let (top, best) = measure(&patch)?;
                floor = Some(floor.map_or(best.clone(), |value| value.max(best)));
                open.push((top, patch));
            }
            let mut floor = floor.ok_or_else(|| unresolved("a saddle has patches"))?;
            for _ in 0..MAX_RANGE_SPLITS {
                dec.check_deadline()?;
                open.retain(|(top, _)| *top > floor);
                if open.iter().all(|(top, _)| top - &floor <= *width) {
                    break;
                }
                let mut next = Vec::with_capacity(open.len() * 2);
                for (top, patch) in std::mem::take(&mut open) {
                    if &top - &floor <= *width {
                        next.push((top, patch));
                        continue;
                    }
                    for piece in along(along(vec![patch], 0), 1) {
                        let (top, best) = measure(&piece)?;
                        floor = floor.max(best);
                        next.push((top, piece));
                    }
                }
                open = next;
            }
            let top = open
                .iter()
                .map(|(top, _)| top)
                .max()
                .cloned()
                .unwrap_or_else(|| floor.clone())
                .max(floor.clone());
            if &top - &floor > *width {
                return Err(unresolved("a saddle's range did not reach its width"));
            }
            answer.push(if lowest { [-top, -floor] } else { [floor, top] });
        }
        let [low, high]: [[Rational; 2]; 2] = answer.try_into().expect("one answer for each end");
        Ok([low, high])
    }

    // The grid for a tolerance, cached on the net.
    fn cached_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<Arc<grid::LocalGrid>, Refusal> {
        if !catallax_numbers::shared_caches() {
            // A worker leaves the shared cache unfilled for the writer.
            return Ok(Arc::new(self.local_grid(tolerance, max_points, dec)?));
        }
        let mut cache = self.data.net.grid.lock().expect("no panic holds this lock");
        #[cfg(feature = "work-span")]
        let span_owner = std::ptr::from_ref(&self.data.net.grid) as usize as u64;
        match cache.as_ref() {
            Some((key, grid)) if key == tolerance => {
                // A hit skips the grid's charge, so it waits for the fill.
                #[cfg(feature = "work-span")]
                catallax_numbers::work_span::read_at("shared cache: saddle grid", span_owner, 0);
                Ok(grid.clone())
            }
            _ => {
                let grid = Arc::new(self.local_grid(tolerance, max_points, dec)?);
                *cache = Some((tolerance.clone(), grid.clone()));
                #[cfg(feature = "work-span")]
                catallax_numbers::work_span::write_at("shared cache: saddle grid", span_owner, 0);
                Ok(grid)
            }
        }
    }

    /// Grid with `chord` at most `tolerance`, cached per tolerance on the net.
    pub fn display_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<EnvelopeGrid, Refusal> {
        let local = self.cached_grid(tolerance, max_points, dec)?;
        let frame = self.data.definition.frame(self.data.net.bits, dec)?;
        Ok(EnvelopeGrid {
            columns: local.grid.columns,
            rows: local.grid.rows,
            points: local.grid.points.iter().map(|p| frame.world(p)).collect(),
            chord: local.grid.chord.clone(),
        })
    }

    /// The chart parameter of each column of [`Self::display_grid`], which need
    /// not be even.
    pub fn display_columns(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<Vec<Rational>, Refusal> {
        Ok(self
            .cached_grid(tolerance, max_points, dec)?
            .parameters
            .clone())
    }

    /// Exact rigid transport of the law and carrier.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }
    pub(crate) fn mapped(&self, pose: &Isometry) -> Self {
        let carrier = self.data.carrier.mapped(pose);
        let cached_rails = self.data.cached_rails.as_ref().map(|rails| {
            let sides = [SplineBoundary::VMin, SplineBoundary::VMax];
            Arc::new(std::array::from_fn(|index| {
                rails[index].mapped_boundary(&carrier, sides[index], pose)
            }))
        });
        Self {
            data: Arc::new(Data {
                definition: self.data.definition.mapped(pose),
                construction: self.data.construction.clone(),
                carrier,
                cached_rails,
                net: self.data.net.clone(),
                bound: self.data.bound.clone(),
                reversing: self.data.reversing != pose.reversing(),
            }),
        }
    }
}

#[cfg(test)]
mod homogeneous_tests {
    use super::*;
    use catallax_numbers::{Outcome, reuse_operation_enclosures, with_operation_budget};

    #[test]
    fn retained_rows_match_the_world_carrier_in_both_axis_directions() {
        let vector =
            |x, y, z| Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z));
        let point = |x, y, z| Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z));
        for (gap, cosine, leaving) in [(1, 4, true), (-1, -4, false)] {
            let wide = Cylinder {
                origin: point(0, 0, 0),
                axis: vector(0, 0, 1),
                radius_sq: Real::from_i64(25),
                ref_dir: vector(1, 0, 0),
            };
            let narrow = Cylinder {
                origin: point(gap, 0, 0),
                axis: vector(0, 3, cosine),
                radius_sq: Real::from_i64(4),
                ref_dir: vector(1, 0, 0),
            };
            let mut dec = Decider::default();
            let def = RollingSaddleDefinition::new(
                &wide,
                &narrow,
                &Real::from(Rational::from_ratio(1, 2).unwrap()),
                [BallSide::Outside; 2],
                leaving,
                &mut dec,
            )
            .unwrap();
            let nodes = produce(&def, FIRST_SPANS, BITS, &mut dec).unwrap();
            let (carrier, homogeneous) = build(&def, &nodes, &mut dec).unwrap();
            let (direct, direct_work) = with_operation_budget(1_000_000, || {
                reuse_operation_enclosures();
                Net::from_homogeneous(&homogeneous, &dec).unwrap()
            });
            let (round_trip, round_trip_work) = with_operation_budget(1_000_000, || {
                reuse_operation_enclosures();
                Net::new(&def, &carrier, &dec).unwrap()
            });
            assert!(!direct_work.exhausted && !round_trip_work.exhausted);
            assert!(direct_work.used < round_trip_work.used);
            for (a, b) in direct.controls.iter().zip(&round_trip.controls) {
                for (a, b) in a.iter().zip(b) {
                    assert!(a.lo() <= b.hi() && b.lo() <= a.hi());
                    assert!(a.hi() - a.lo() <= b.hi() - b.lo());
                }
            }
            for index in [0, 1, 2, homogeneous.len() / 2, homogeneous.len() - 1] {
                let local = def.local(&carrier.controls()[index]);
                for (axis, value) in local.iter().enumerate() {
                    let residual = &(value * &carrier.weights()[index]) - &homogeneous[index][axis];
                    assert!(matches!(dec.sign(&residual), Outcome::ProvenCoincident(_)));
                }
            }
            dec.set_deadline(Some(
                std::time::Instant::now() - std::time::Duration::from_secs(1),
            ));
            assert_eq!(
                Net::from_homogeneous(&homogeneous, &dec).unwrap_err().code,
                "numbers.deadline"
            );
        }
    }
}
