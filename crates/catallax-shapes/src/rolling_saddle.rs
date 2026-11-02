//! A ball rolling in the closed valley where a narrow cylinder passes through a
//! wide one at an angle.
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

use crate::rolling_envelope::{Jet, combine, cubic, elevate, hermite, product};
use crate::{BSplineCurve, BSplineSurface, Cylinder, Isometry, Point3, Pose, SplineBoundary};
use crate::{Vector3, math};
use catallax_numbers::{Ball, Decider, EvidenceBound, EvidenceUnit, Rational, Real, Refusal, Sign};
use std::sync::Arc;

/// Kernel-made approximation target, as for a rolling envelope.
pub const ROLLING_SADDLE_TARGET_NM: i64 = 100;
const HARD_CAP_NM: i64 = 1_000;
// Spans per half turn: the fit doubles from the first to the last.
const FIRST_SPANS: usize = 8;
const LAST_SPANS: usize = 128;
const DEGREE: usize = 9;
const BITS: u32 = 128;
// Node data keep this many fractional bits, far below the target.
const NODE_BITS: i64 = 48;
// Cells of one span in the regularity proof, by doubling.
const MAX_CELLS: usize = 64;

const ROWS: usize = 6;
const HEIGHT: usize = 0;
const LIFT: usize = 1;
const HALF: usize = 2;
const TANGENT: usize = 3;
const AXIAL: usize = 4;
const WEIGHT: usize = 5;

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

/// The exact law of a ball of radius `radius` resting on two nonparallel
/// cylinders, on one of the two rings a narrow cylinder cuts through a wide one.
#[derive(Clone, Debug)]
pub struct RollingSaddleDefinition {
    origin: Point3,
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
        let between = narrow.origin.minus(&wide.origin);
        let gap = between.dot(&across);
        let station = -&(&between.dot(&along) / &sine);
        let origin = wide
            .origin
            .plus(&a.scaled(&(&between.dot(&a) + &(&station * &cosine))));
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
        let reach = math::exact_sqrt(&gap.square());
        for (big, small) in [(&offsets[0], &offsets[1]), (&wide_radius, &narrow_radius)] {
            if math::sign(&(&(big - small) - &reach), dec)? != Sign::Positive {
                return Err(invalid(
                    "the narrow cylinder does not pass through the wide one",
                ));
            }
        }
        Ok(Self {
            origin,
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
            branch: if leaving { 1 } else { -1 },
        })
    }

    fn mapped(&self, pose: &Isometry) -> Self {
        Self {
            origin: pose.apply_point(&self.origin),
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

    fn offset(&self, which: usize) -> Real {
        let radius = if which == 0 { &self.wide } else { &self.narrow };
        radius + &(&self.radius * &Real::from_i64(self.sides[which]))
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
    fn foot(&self) -> Point3 {
        self.origin.plus(&self.across.scaled(&self.gap))
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
        let chart = chart_of(u)?;
        let s = u - &q(2 * chart);
        let wide = Rational::one() + &s * &s;
        let turn = q(1 - 2 * chart);
        let cos = real(&((Rational::one() - &s * &s) * turn.clone() / wide.clone()));
        let sin = real(&(&s * &q(2) * turn / wide));
        let [big, small] = [0, 1].map(|which| self.offset(which));
        let x = &self.gap + &(&small * &cos);
        let root = math::exact_sqrt(&(&big.square() - &x.square()));
        let y = &root * &Real::from_i64(self.branch);
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
    one: Ball,
    two: Ball,
}

impl Constants {
    fn new(def: &RollingSaddleDefinition, bits: u32, dec: &Decider) -> Result<Self, Refusal> {
        let e = |value: &Real| dec.enclose_at(value, bits);
        let one = b(&Rational::one(), bits);
        let sine = e(&def.sine)?;
        let over_sine = one
            .checked_div(&sine)
            .ok_or_else(|| unresolved("saddle axes are not separated from parallel"))?;
        Ok(Self {
            gap: e(&def.gap)?,
            cosine: e(&def.cosine)?,
            sine,
            over_sine,
            wide: e(&def.wide)?,
            narrow: e(&def.narrow)?,
            big: e(&def.offset(0))?,
            small: e(&def.offset(1))?,
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
    tangent: Jet,
    axial: Jet,
    weight: Option<Jet>,
    cosine: Jet,
    centre: [Jet; 3],
    wide: [Jet; 3],
    narrow: [Jet; 3],
}

fn rows(k: &Constants, chart: i64, s: &Jet, fitted: Option<&Jet>) -> Result<Rows, Refusal> {
    let bits = k.one.precision();
    let one = k.c(&k.one);
    let turn = b(&q(1 - 2 * chart), bits);
    let square = s.square();
    let den = one.add(&square);
    let cos = one.sub(&square).scale(&turn).div(&den)?;
    let sin = s.scale(&k.two).scale(&turn).div(&den)?;
    let x = k.c(&k.gap).add(&cos.scale(&k.small));
    let root = k.c(&k.big.square()).sub(&x.square()).sqrt()?;
    let y = root.scale(&k.branch);
    let height = y
        .sub(&sin.scale(&(&k.small * &k.cosine)))
        .scale(&k.over_sine);
    let lift = y
        .scale(&k.cosine)
        .sub(&sin.scale(&k.small))
        .scale(&k.over_sine);
    let half = x.div(&k.c(&k.big).add(&root))?;
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
        Some(w) => Some(
            one_d
                .mul(&one.add(&w.square()))
                .div(&den.scale(&k.two))?
                .sqrt()?,
        ),
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
            let point = Jet::linear(b(&s, bits), k.one.clone());
            let exact = rows(&k, chart, &point, None)?;
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
            let fitted = Jet::linear(b(&values[HALF], bits), b(&slopes[HALF], bits));
            let weight = rows(&k, chart, &point, Some(&fitted))?.weight;
            let [value, slope] = weight.expect("weight requested").data();
            values[WEIGHT] = snap(&value.mid());
            slopes[WEIGHT] = snap(&slope.mid());
            let arriving = if i == 0 && chart == 1 {
                let end = Jet::linear(k.one.clone(), k.one.clone());
                let weight = rows(&k, 0, &end, Some(&fitted))?.weight;
                snap(&weight.expect("weight requested").data()[1].mid())
            } else {
                slopes[WEIGHT].clone()
            };
            nodes.push(Node {
                u: &s + &q(2 * chart),
                values,
                slopes,
                arriving,
            });
        }
    }
    let first = nodes[0].clone();
    let fitted = Jet::linear(b(&first.values[HALF], bits), b(&first.slopes[HALF], bits));
    let end = Jet::linear(k.one.clone(), k.one.clone());
    let weight = rows(&k, 1, &end, Some(&fitted))?.weight;
    let arriving = snap(&weight.expect("weight requested").data()[1].mid());
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
        let exact = rows(&k, chart, &point, Some(&fitted))?;
        Ok([
            exact.height.data(),
            exact.lift.data(),
            exact.half.data(),
            exact.tangent.data(),
            exact.axial.data(),
            exact.weight.expect("weight requested").data(),
        ])
    };
    let mut worst = Rational::zero();
    for pair in nodes.windows(2) {
        dec.check_deadline()?;
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
        let whole = rows(&k, chart, &s, Some(&fitted))?;
        let weight_jet = whole.weight.as_ref().expect("weight requested");
        let fourth = [
            &whole.height.0[4],
            &whole.lift.0[4],
            &whole.half.0[4],
            &whole.tangent.0[4],
            &whole.axial.0[4],
            &weight_jet.0[4],
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
        let narrow = error(HEIGHT);
        let apex = error(TANGENT) + error(AXIAL) + error(HEIGHT);
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
        worst = worst.max(wide.max(narrow).max(apex) + spread);
    }
    Ok(worst)
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
    let s = vec![real(&ends[0]), real(&ends[1])];
    let square = product(&s, &s);
    let turn = Real::from_i64(1 - 2 * chart);
    let ones = vec![one.clone(); 3];
    let den = combine(&ones, &square, &one, &one);
    let cos = combine(&ones, &square, &turn, &-&turn);
    let sin = scaled(&elevate(&s, 2), &(&turn * &Real::from_i64(2)));
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
    let branch = Real::from_i64(def.branch);
    let wide_den = combine(&ones, &half_square, &one, &one);
    let wide = [
        scaled(&elevate(&half, 6), &(&def.wide * &Real::from_i64(2))),
        scaled(
            &combine(&ones, &half_square, &one, &minus),
            &(&def.wide * &branch),
        ),
        product(&f(LIFT), &wide_den),
        wide_den,
    ];
    [wide, middle, narrow].map(|row| row.map(|c| elevate(&c, DEGREE)))
}

fn build(
    def: &RollingSaddleDefinition,
    nodes: &[Node],
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let foot = def.foot();
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
    let count = spans * DEGREE + 1;
    let mut controls = Vec::with_capacity(count * 3);
    let mut weights = Vec::with_capacity(count * 3);
    let mut knots = vec![q(-1); DEGREE + 1];
    for (index, pair) in nodes.windows(2).enumerate() {
        dec.check_deadline()?;
        let rows = span_rows(def, pair);
        let last = index + 1 == spans;
        for i in 0..DEGREE + usize::from(last) {
            for (row, data) in rows.iter().enumerate() {
                controls.push(place(
                    row,
                    [&data[0][i], &data[1][i], &data[2][i], &data[3][i]],
                ));
                weights.push(data[3][i].clone());
            }
        }
        let repeat = DEGREE + usize::from(last);
        knots.extend(std::iter::repeat_n(pair[1].u.clone(), repeat));
    }
    let v = [0, 0, 0, 1, 1, 1].map(q).to_vec();
    BSplineSurface::new([DEGREE, 2], [knots, v], [count, 3], controls, weights, dec)
}

/// A closed ring of the ball's surface. The carrier is the represented geometry
/// over `[-1, 3] x [0, 1]`: `u` runs once round the ring and `v` from the wide
/// rail to the narrow one.
#[derive(Clone, Debug)]
pub struct RollingSaddle {
    data: Arc<Data>,
}

#[derive(Debug)]
struct Data {
    definition: RollingSaddleDefinition,
    nodes: Vec<Node>,
    carrier: BSplineSurface,
    bound: EvidenceBound,
}

impl RollingSaddle {
    /// Fits the fewest uniform spans, by doubling, whose recomputed bound meets
    /// the 100 nm target. A fit that meets only the 1 um cap is accepted and
    /// reported by its bound, and anything larger refuses.
    pub fn fit(definition: RollingSaddleDefinition, dec: &mut Decider) -> Result<Self, Refusal> {
        let bits = BITS.min(dec.policy().max_precision_bits);
        let target = Rational::from_ratio(ROLLING_SADDLE_TARGET_NM, 1_000_000).expect("nonzero");
        let cap = Rational::from_ratio(HARD_CAP_NM, 1_000_000).expect("nonzero");
        let mut spans = FIRST_SPANS;
        let (nodes, upper) = loop {
            let nodes = produce(&definition, spans, bits, dec)?;
            match recompute_bound(&definition, &nodes, bits, dec) {
                Ok(upper) if upper <= target || spans >= LAST_SPANS => break (nodes, upper),
                Err(error) if spans >= LAST_SPANS => return Err(error),
                _ => spans *= 2,
            }
        };
        if upper > cap {
            return Err(invalid("the saddle fit misses the 1 um approximation cap"));
        }
        check_regular(&definition, &nodes, bits, dec)?;
        let carrier = build(&definition, &nodes, dec)?;
        Ok(Self {
            data: Arc::new(Data {
                definition,
                nodes,
                carrier,
                bound: EvidenceBound {
                    upper,
                    unit: EvidenceUnit::Millimetres,
                },
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

    /// Fitted spans round the ring.
    pub fn span_count(&self) -> usize {
        self.data.nodes.len() - 1
    }

    /// The rails on the wide and narrow cylinders, the carrier's own boundary
    /// rows, so a neighbour can share them exactly.
    pub fn rails(&self) -> [BSplineCurve; 2] {
        [SplineBoundary::VMin, SplineBoundary::VMax]
            .map(|side| BSplineCurve::boundary(&self.data.carrier, side))
    }

    /// Exact rigid transport of the law and carrier.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }
    pub(crate) fn mapped(&self, pose: &Isometry) -> Self {
        Self {
            data: Arc::new(Data {
                definition: self.data.definition.mapped(pose),
                nodes: self.data.nodes.clone(),
                carrier: self.data.carrier.mapped(pose),
                bound: self.data.bound.clone(),
            }),
        }
    }
}
