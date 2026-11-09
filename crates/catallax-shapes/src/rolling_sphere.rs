//! A ball rolling between a sphere and a cylinder whose axis misses the
//! sphere's centre.
//!
//! Take the sphere's centre as the origin, the unit cylinder axis as `axis`,
//! `across` from the centre square to the axis, and `along = axis x across`,
//! so the axis runs through `(gap, 0, z)`. Offsetting the sphere and the
//! cylinder by the ball radius gives radii `P` and `p`. A ball centre on the
//! cylinder's offset at angle `t` about the axis has
//!
//! ```text
//! x = gap + p cos t,  y = p sin t,  z = branch sqrt(P^2 - x^2 - y^2)
//! ```
//!
//! When `gap + p < P` the root never closes, so `t` charts one of two whole
//! rings, one above and one below the centre. The root makes each ring an
//! elliptic curve, so no rational surface matches it and the carrier is a fit.
//!
//! The carrier has three rows. The cylinder rail keeps the exact angle `t` with
//! a fitted height, so it lies exactly on the cylinder. The sphere rail takes
//! fitted stereographic coordinates of its contact ray from the pole on the
//! far side of its ring, so it lies exactly on the sphere. The middle row is
//! the apex of each section's circle, held in the cylinder's tangent plane at
//! the rail with a fitted middle weight, so every section is a true conic. The
//! two-sided bound is recomputed from the retained node data and interval
//! Taylor jets of the law, as for [`RollingSaddle`](crate::RollingSaddle).
//!
//! `u` runs over `[-1, 3]`. On `[-1, 1]` the half-angle tangent of `t` is `u`,
//! and on `[1, 3]` the same half turn repeats from the opposite side, so the
//! cylinder rail's denominator is two at every chart end and the closed
//! carrier's weights close.

use crate::rational_sheet::{LocalFrame, Patch, along, cover, split_net};
use crate::rolling_envelope::{Jet, combine, cubic, elevate, hermite, product};
use crate::{
    BSplineCurve, BSplineSurface, BallSide, Cylinder, Isometry, ParameterTable, Plane, Point3,
    Pose, Sphere, SplineBoundary, Vector3, math,
};
use catallax_numbers::{Ball, Decider, EvidenceBound, EvidenceUnit, Rational, Real, Refusal, Sign};
use std::{collections::BTreeMap, sync::Arc};

mod line;
#[cfg(test)]
mod tests;

/// Kernel-made approximation target, as for a rolling saddle.
pub const ROLLING_SPHERE_TARGET_NM: i64 = 100;
const HARD_CAP_NM: i64 = 1_000;
// Spans per half turn: the fit doubles from the first to the last.
const FIRST_SPANS: usize = 8;
const LAST_SPANS: usize = 128;
const DEGREE: usize = 8;
const BITS: u32 = 128;
// Node data keep this many fractional bits, far below the target.
const NODE_BITS: i64 = 32;
// Cells of one span in the regularity proof, by doubling.
const MAX_CELLS: usize = 64;
// Halvings of the open patches before a range gives up.
const MAX_RANGE_SPLITS: usize = 40;

const ROWS: usize = 6;
const HEIGHT: usize = 0;
const SIGMA: usize = 1;
const TAU: usize = 2;
const TANGENT: usize = 3;
const AXIAL: usize = 4;
const WEIGHT: usize = 5;

fn invalid(detail: impl Into<String>) -> Refusal {
    Refusal::new("shapes.rolling-sphere", detail)
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
fn real(value: &Rational) -> Real {
    Real::from(value.clone())
}
fn side_sign(side: BallSide) -> i64 {
    match side {
        BallSide::Outside => 1,
        BallSide::Inside => -1,
    }
}

/// The exact law of a constant-radius ball resting on a sphere and on a
/// cylinder whose axis misses the sphere's centre, along one regular ring.
#[derive(Clone, Debug)]
pub struct RollingSphereDefinition {
    centre: Point3,
    across: Vector3,
    along: Vector3,
    axis: Vector3,
    sphere: Real,
    cylinder: Real,
    radius: Real,
    gap: Real,
    // The ball's side of the sphere and of the cylinder: 1 outside, -1 inside.
    sides: [i64; 2],
    branch: i64,
}

impl RollingSphereDefinition {
    /// The ring through a point of both supports. Checks that the offsets are
    /// positive, that the axis misses the centre, and that the cylinder and its
    /// offset both pass through the sphere and its offset without touching
    /// their far side, so each meets it in two separate rings. The point's side
    /// of the plane through the centre square to the axis names the ring.
    pub fn new(
        sphere: &Sphere,
        cylinder: &Cylinder,
        radius: &Real,
        sides: [BallSide; 2],
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        if math::sign(radius, dec)? != Sign::Positive {
            return Err(invalid("ball radius must be positive"));
        }
        let length = math::exact_sqrt(&cylinder.axis.norm_squared());
        if math::sign(&length, dec)? != Sign::Positive {
            return Err(invalid("the cylinder axis is zero"));
        }
        let axis = cylinder.axis.scaled(&(&Real::one() / &length));
        let offset = cylinder.origin.minus(&sphere.center);
        let square = offset.minus(&axis.scaled(&offset.dot(&axis)));
        let gap = math::exact_sqrt(&square.norm_squared());
        if math::sign(&gap, dec)? != Sign::Positive {
            return Err(invalid(
                "a cylinder through the sphere's centre gives a coaxial ring",
            ));
        }
        let across = square.scaled(&(&Real::one() / &gap));
        let along = axis.cross(&across);
        let sides = sides.map(side_sign);
        let cylinder_radius = math::exact_sqrt(&cylinder.radius_sq);
        let definition = Self {
            centre: sphere.center.clone(),
            across,
            along,
            axis,
            sphere: sphere.radius.clone(),
            cylinder: cylinder_radius,
            radius: radius.clone(),
            gap,
            sides,
            branch: 1,
        };
        for value in [definition.offset(0), definition.offset(1)] {
            if math::sign(&value, dec)? != Sign::Positive {
                return Err(invalid("the ball does not fit inside a support"));
            }
        }
        for (big, small) in [
            (definition.offset(0), definition.offset(1)),
            (definition.sphere.clone(), definition.cylinder.clone()),
        ] {
            if math::sign(&(&(&big - &small) - &definition.gap), dec)? != Sign::Positive {
                return Err(invalid(
                    "the cylinder does not pass cleanly through the sphere",
                ));
            }
        }
        let [_, _, z] = definition.local(point);
        match math::sign(&z, dec)? {
            Sign::Positive => Ok(definition),
            Sign::Negative => Ok(Self {
                branch: -1,
                ..definition
            }),
            Sign::Zero => Err(invalid("the point is not on a regular ring")),
        }
    }

    /// Whether a tangent at a point of the sharp edge the chart rounds runs the
    /// way a growing `u` does. A ring turns from `across` toward `along` about
    /// the cylinder axis throughout, so one tangent decides the whole ring. The
    /// test reads the chart's own frame, so a mirrored copy keeps its answer.
    pub fn runs_with_chart(
        &self,
        point: &Point3,
        tangent: &Vector3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let radial = self.radial(point);
        let turn = &(&radial.dot(&self.across) * &tangent.dot(&self.along))
            - &(&radial.dot(&self.along) * &tangent.dot(&self.across));
        match math::sign(&turn, dec)? {
            Sign::Positive => Ok(true),
            Sign::Negative => Ok(false),
            Sign::Zero => Err(invalid("the tangent does not turn round the cylinder axis")),
        }
    }

    // A mirror leaves the frame left-handed. Every row is placed through the
    // same mapped frame, so the charts keep their meaning.
    fn mapped(&self, pose: &Isometry) -> Self {
        Self {
            centre: pose.apply_point(&self.centre),
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

    /// The sphere's centre, the law's origin.
    pub fn centre(&self) -> &Point3 {
        &self.centre
    }

    /// The sphere's radius.
    pub fn sphere_radius(&self) -> &Real {
        &self.sphere
    }

    /// The cylinder's radius.
    pub fn cylinder_radius(&self) -> &Real {
        &self.cylinder
    }

    /// The unit cylinder axis.
    pub fn axis(&self) -> &Vector3 {
        &self.axis
    }

    /// The ball's side of the sphere and of the cylinder: 1 outside, -1 inside.
    pub fn sides(&self) -> [i64; 2] {
        self.sides
    }

    /// The sign of the ring's coordinate along the axis.
    pub fn branch(&self) -> i64 {
        self.branch
    }

    /// The distance from the sphere's centre to the cylinder's axis.
    pub fn gap(&self) -> &Real {
        &self.gap
    }

    /// The unit direction from the sphere's centre to the cylinder's axis,
    /// square to the axis. The chart angle is zero along it.
    pub fn across(&self) -> &Vector3 {
        &self.across
    }

    /// `axis x across`, where the chart angle is a quarter turn. A mirrored
    /// copy keeps the mapped vector, so there it is `across x axis`.
    pub fn along(&self) -> &Vector3 {
        &self.along
    }

    /// `across x along`, the unit direction the chart angle rises about: a
    /// growing `u` turns the ring counterclockwise seen from its tip. It is the
    /// axis in the law's own frame and the reversed axis on a mirrored copy.
    pub fn turn_axis(&self) -> Vector3 {
        self.across.cross(&self.along)
    }

    /// The radius of the ball centre's offset of the sphere (`0`) or the
    /// cylinder (`1`).
    pub fn offset_radius(&self, which: usize) -> Real {
        self.offset(which)
    }

    // 0 is the sphere's offset, 1 the cylinder's.
    fn offset(&self, which: usize) -> Real {
        let radius = if which == 0 {
            &self.sphere
        } else {
            &self.cylinder
        };
        radius + &(&self.radius * &Real::from_i64(self.sides[which]))
    }
    fn frame(&self, bits: u32, dec: &Decider) -> Result<LocalFrame, Refusal> {
        LocalFrame::new(
            &self.centre,
            [&self.across, &self.along, &self.axis],
            bits,
            dec,
        )
    }
    fn budgeted_frame(&self, bits: u32, dec: &mut Decider) -> Result<LocalFrame, Refusal> {
        LocalFrame::budgeted(
            &self.centre,
            [&self.across, &self.along, &self.axis],
            bits,
            dec,
        )
    }
    fn local(&self, point: &Point3) -> [Real; 3] {
        let d = point.minus(&self.centre);
        [d.dot(&self.across), d.dot(&self.along), d.dot(&self.axis)]
    }
    // A point's offset from the cylinder axis.
    fn radial(&self, point: &Point3) -> Vector3 {
        let d = point
            .minus(&self.centre)
            .minus(&self.across.scaled(&self.gap));
        d.minus(&self.axis.scaled(&d.dot(&self.axis)))
    }
    fn place(&self, x: &Real, y: &Real, z: &Real) -> Point3 {
        self.centre
            .plus(&self.across.scaled(x))
            .plus(&self.along.scaled(y))
            .plus(&self.axis.scaled(z))
    }

    /// The law at a rational chart parameter of `[-1, 3]`, exactly: the ball
    /// centre and its contacts on the sphere and the cylinder.
    pub fn section(&self, u: &Rational) -> Result<[Point3; 3], Refusal> {
        let chart = chart_of(u)?;
        let s = u - &q(2 * chart);
        let den = Rational::one() + &s * &s;
        let turn = q(1 - 2 * chart);
        let cos = real(&((Rational::one() - &s * &s) * turn.clone() / den.clone()));
        let sin = real(&(&s * &q(2) * turn / den));
        let [big, small] = [0, 1].map(|which| self.offset(which));
        let x = &self.gap + &(&small * &cos);
        let y = &small * &sin;
        let root = math::exact_sqrt(&(&(&big.square() - &x.square()) - &y.square()));
        let z = &root * &Real::from_i64(self.branch);
        let scale = &self.sphere / &big;
        Ok([
            self.place(&x, &y, &z),
            self.place(&(&x * &scale), &(&y * &scale), &(&z * &scale)),
            self.place(
                &(&self.gap + &(&self.cylinder * &cos)),
                &(&self.cylinder * &sin),
                &z,
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
    big: Ball,
    small: Ball,
    sphere: Ball,
    cylinder: Ball,
    radius: Ball,
    sides: [Ball; 2],
    branch: Ball,
    one: Ball,
    two: Ball,
}

impl Constants {
    fn new(def: &RollingSphereDefinition, bits: u32, dec: &Decider) -> Result<Self, Refusal> {
        let e = |value: &Real| dec.enclose_at(value, bits);
        Ok(Self {
            gap: e(&def.gap)?,
            big: e(&def.offset(0))?,
            small: e(&def.offset(1))?,
            sphere: e(&def.sphere)?,
            cylinder: e(&def.cylinder)?,
            radius: e(&def.radius)?,
            sides: def.sides.map(|side| b(&q(side), bits)),
            branch: b(&q(def.branch), bits),
            one: b(&Rational::one(), bits),
            two: b(&q(2), bits),
        })
    }
    fn c(&self, value: &Ball) -> Jet {
        Jet::constant(value.clone())
    }
}

// The law at one station, as jets in the chart parameter.
struct Rows {
    height: Jet,
    sigma: Jet,
    tau: Jet,
    tangent: Jet,
    axial: Jet,
    weight: Option<Jet>,
    cosine: Jet,
    centre: [Jet; 3],
    sphere: [Jet; 3],
    cylinder: [Jet; 3],
}

// `fitted` holds the sphere rail's fitted stereographic coordinates.
fn rows(k: &Constants, chart: i64, s: &Jet, fitted: Option<[&Jet; 2]>) -> Result<Rows, Refusal> {
    let bits = k.one.precision();
    let one = k.c(&k.one);
    let turn = b(&q(1 - 2 * chart), bits);
    let square = s.square();
    let den = one.add(&square);
    let cos = one.sub(&square).scale(&turn).div(&den)?;
    let sin = s.scale(&k.two).scale(&turn).div(&den)?;
    let x = k.c(&k.gap).add(&cos.scale(&k.small));
    let y = sin.scale(&k.small);
    let root = k
        .c(&k.big.square())
        .sub(&x.square())
        .sub(&y.square())
        .sqrt()?;
    let height = root.scale(&k.branch);
    let over_big = k
        .one
        .checked_div(&k.big)
        .ok_or_else(|| unresolved("sphere offset is not separated from zero"))?;
    // Stereographic coordinates of the unit ray (x, y, z) / P from the pole on
    // the far side of the ring, where the ray's lift `1 + |z| / P` exceeds one.
    let lift = one.add(&root.scale(&over_big));
    let sigma = x.scale(&over_big).div(&lift)?;
    let tau = y.scale(&over_big).div(&lift)?;
    // Unit rays from the centre: -side_c (cos, sin, 0) to the cylinder contact
    // and -side_s (x, y, z) / P to the sphere contact.
    let both = &k.sides[0] * &k.sides[1];
    let meet = k
        .c(&k.gap)
        .mul(&cos)
        .add(&k.c(&k.small))
        .scale(&(&both * &over_big));
    let one_d = one.add(&meet);
    // Apex M = C + r (a_c + a_s) / (1 + d). Against the cylinder contact it
    // moves r (a_s - d a_c) / (1 + d), which has no part along a_c.
    let reach = &(&k.radius * &k.sides[0]) * &over_big;
    let tangent = sin.scale(&(&reach * &k.gap)).div(&one_d)?;
    let axial = height.scale(&reach).neg().div(&one_d)?;
    let half_turn = k
        .one
        .checked_div(&k.two)
        .expect("two is separated from zero");
    let cosine = one_d.scale(&half_turn).sqrt()?;
    // The middle weight cos(half) sqrt(w0 w2) over the cylinder rail's
    // denominator, with the sphere rail's denominator from its fitted rows.
    let weight = match fitted {
        None => None,
        Some([fs, ft]) => Some(
            one_d
                .mul(&one.add(&fs.square()).add(&ft.square()))
                .div(&den.scale(&k.two))?
                .sqrt()?,
        ),
    };
    let scale = k
        .sphere
        .checked_div(&k.big)
        .ok_or_else(|| unresolved("sphere offset is not separated from zero"))?;
    Ok(Rows {
        centre: [x.clone(), y.clone(), height.clone()],
        sphere: [x.scale(&scale), y.scale(&scale), height.scale(&scale)],
        cylinder: [
            k.c(&k.gap).add(&cos.scale(&k.cylinder)),
            sin.scale(&k.cylinder),
            height.clone(),
        ],
        height,
        sigma,
        tau,
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
    // The weight's slope on the span ending here. The cylinder rail's
    // denominator turns round at a chart end, so the weight has a corner there
    // and only there.
    arriving: Rational,
}

fn fitted_pair(data: &[[Ball; 2]; ROWS]) -> [Jet; 2] {
    [SIGMA, TAU].map(|j| Jet::linear(data[j][0].clone(), data[j][1].clone()))
}

// The weight's slope where a chart ends, at `s = 1` of `chart`, given the
// fitted rows' data there.
fn ending_weight(k: &Constants, chart: i64, node: &Node, bits: u32) -> Result<Rational, Refusal> {
    let data: [[Ball; 2]; ROWS] =
        std::array::from_fn(|j| [b(&node.values[j], bits), b(&node.slopes[j], bits)]);
    let [fs, ft] = fitted_pair(&data);
    let end = Jet::linear(k.one.clone(), k.one.clone());
    let weight = rows(k, chart, &end, Some([&fs, &ft]))?.weight;
    Ok(snap(&weight.expect("weight requested").data()[1].mid()))
}

fn node_at(k: &Constants, chart: i64, s: &Rational, bits: u32) -> Result<Node, Refusal> {
    let point = Jet::linear(b(s, bits), k.one.clone());
    let exact = rows(k, chart, &point, None)?;
    let mut values: [Rational; ROWS] = std::array::from_fn(|_| Rational::zero());
    let mut slopes = values.clone();
    let plain = [
        (HEIGHT, &exact.height),
        (SIGMA, &exact.sigma),
        (TAU, &exact.tau),
        (TANGENT, &exact.tangent),
        (AXIAL, &exact.axial),
    ];
    for (j, jet) in plain {
        let [value, slope] = jet.data();
        values[j] = snap(&value.mid());
        slopes[j] = snap(&slope.mid());
    }
    let [fs, ft] = [SIGMA, TAU].map(|j| Jet::linear(b(&values[j], bits), b(&slopes[j], bits)));
    let weight = rows(k, chart, &point, Some([&fs, &ft]))?.weight;
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

// One definition and precision own these samples. Endpoint rows depend on
// the snapped sigma/tau data, which are unchanged by span refinement.
#[derive(Default)]
struct FitSamples {
    nodes: BTreeMap<(i64, Rational), Node>,
    rows: BTreeMap<(i64, Rational), [[Ball; 2]; ROWS]>,
    closing_weight: Option<Rational>,
}

fn produce(
    k: &Constants,
    spans: usize,
    bits: u32,
    samples: &mut FitSamples,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    let mut nodes = Vec::with_capacity(2 * spans + 1);
    for chart in 0..2 {
        for i in 0..spans {
            dec.check_deadline()?;
            let s =
                Rational::from_ratio(2 * i as i64 - spans as i64, spans as i64).expect("nonzero");
            let key = (chart, s.clone());
            let node = match samples.nodes.entry(key) {
                std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
                std::collections::btree_map::Entry::Vacant(entry) => {
                    let mut node = node_at(k, chart, &s, bits)?;
                    // The rail denominator turns round where the first chart ends.
                    if i == 0 && chart == 1 {
                        node.arriving = ending_weight(k, 0, &node, bits)?;
                    }
                    entry.insert(node)
                }
            };
            nodes.push(node.clone());
        }
    }
    let first = nodes[0].clone();
    let arriving = match &samples.closing_weight {
        Some(value) => value.clone(),
        None => {
            let value = ending_weight(k, 1, &first, bits)?;
            samples.closing_weight = Some(value.clone());
            value
        }
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

fn recompute_bound(
    k: &Constants,
    nodes: &[Node],
    bits: u32,
    samples: &mut FitSamples,
    stop_above: Option<&Rational>,
    dec: &mut Decider,
) -> Result<Option<Rational>, Refusal> {
    let mut worst = Rational::zero();
    for pair in nodes.windows(2) {
        dec.check_deadline()?;
        worst = worst.max(span_bound(k, pair, bits, samples)?);
        if stop_above.is_some_and(|target| worst > *target) {
            return Ok(None);
        }
    }
    Ok(Some(worst))
}

// The same bound as a rolling saddle's. Both surfaces are rational quadratics
// in v with shared rail weights, so with `eps` the relative middle-weight error
//   |S - E| <= max_i |s_i - e_i| + (eps / (2 (1 - eps))) r / cos(half)
// and each row error is the cubic Hermite remainder plus the node data's miss.
// The inverse stereographic map scales lengths by at most two, so the sphere
// rail moves at most `2 R (|d sigma| + |d tau|)`.
fn span_bound(
    k: &Constants,
    pair: &[Node],
    bits: u32,
    samples: &mut FitSamples,
) -> Result<Rational, Refusal> {
    let enclosed = |node: &Node, ending: bool| -> [[Ball; 2]; ROWS] {
        std::array::from_fn(|j| {
            let slope = if ending && j == WEIGHT {
                &node.arriving
            } else {
                &node.slopes[j]
            };
            [b(&node.values[j], bits), b(slope, bits)]
        })
    };
    let mut exact_at = |chart: i64,
                        s: &Rational,
                        data: &[[Ball; 2]; ROWS]|
     -> Result<[[Ball; 2]; ROWS], Refusal> {
        let key = (chart, s.clone());
        if let Some(exact) = samples.rows.get(&key) {
            return Ok(exact.clone());
        }
        let point = Jet::linear(b(s, bits), k.one.clone());
        let [fs, ft] = fitted_pair(data);
        let exact = rows(k, chart, &point, Some([&fs, &ft]))?;
        let data = [
            exact.height.data(),
            exact.sigma.data(),
            exact.tau.data(),
            exact.tangent.data(),
            exact.axial.data(),
            exact.weight.expect("weight requested").data(),
        ];
        samples.rows.insert(key, data.clone());
        Ok(data)
    };
    let (chart, ends) = span_chart(pair);
    let span = [&ends[0], &ends[1]];
    let left = enclosed(&pair[0], false);
    let right = enclosed(&pair[1], true);
    let left_exact = exact_at(chart, &ends[0], &left)?;
    let right_exact = exact_at(chart, &ends[1], &right)?;
    let h = &ends[1] - &ends[0];
    let s = Jet::linear(b(&ends[0], bits).hull(&b(&ends[1], bits)), k.one.clone());
    let [fs, ft] = [SIGMA, TAU].map(|j| {
        hermite(
            &s,
            span,
            [[&left[j][0], &left[j][1]], [&right[j][0], &right[j][1]]],
            bits,
        )
    });
    let (fs, ft) = (fs?, ft?);
    let whole = rows(k, chart, &s, Some([&fs, &ft]))?;
    let weight_jet = whole.weight.as_ref().expect("weight requested");
    let fourth = [
        &whole.height.0[4],
        &whole.sigma.0[4],
        &whole.tau.0[4],
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
    let sphere = q(2) * k.sphere.hi() * (error(SIGMA) + error(TAU));
    let cylinder = error(HEIGHT);
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
    Ok(sphere.max(cylinder).max(apex) + spread)
}

// Proves the law regular on every span: the centre path moves, the contact
// rays are neither equal nor opposite, and every point of each section moves
// along the path, so the swept surface has independent derivatives.
fn check_regular(
    def: &RollingSphereDefinition,
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
                let [sphere, cylinder] =
                    [&at.sphere, &at.cylinder].map(|rail| &dot(&first(rail), &path) - &speed);
                let half = b(&Rational::from_ratio(1, 2).expect("nonzero"), bits);
                let middle = &(&sphere + &cylinder) * &half;
                let over = middle
                    .checked_div(cosine)
                    .ok_or_else(|| unresolved("section angle is not resolved"))?;
                let spread = magnitude(&(&(&cylinder - &sphere) * &half));
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

// Local homogeneous Bernstein rows of one span, [sphere, middle, cylinder],
// all in the sphere's frame.
fn span_rows(def: &RollingSphereDefinition, pair: &[Node]) -> [[Vec<Real>; 4]; 3] {
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
    let cos = combine(&ones, &square, &turn, &-&turn);
    let sin = scaled(&elevate(&s, 2), &(&turn * &Real::from_i64(2)));
    let den = combine(&ones, &square, &one, &one);
    let rail_x = combine(
        &scaled(&den, &def.gap),
        &scaled(&cos, &def.cylinder),
        &one,
        &one,
    );
    let rail_y = scaled(&sin, &def.cylinder);
    let cylinder = [
        rail_x.clone(),
        rail_y.clone(),
        product(&f(HEIGHT), &den),
        den.clone(),
    ];
    let tangent = f(TANGENT);
    let x = combine(&rail_x, &product(&tangent, &sin), &one, &minus);
    let y = combine(&rail_y, &product(&tangent, &cos), &one, &one);
    let height = combine(&f(HEIGHT), &f(AXIAL), &one, &one);
    let weight = f(WEIGHT);
    let middle = [
        product(&weight, &x),
        product(&weight, &y),
        product(&weight, &product(&height, &den)),
        product(&weight, &den),
    ];
    let (sigma, tau) = (f(SIGMA), f(TAU));
    let both = combine(&product(&sigma, &sigma), &product(&tau, &tau), &one, &one);
    let ones = vec![one.clone(); 7];
    let twice = &def.sphere * &Real::from_i64(2);
    let sphere = [
        scaled(&elevate(&sigma, 6), &twice),
        scaled(&elevate(&tau, 6), &twice),
        scaled(
            &combine(&ones, &both, &one, &minus),
            &(&def.sphere * &Real::from_i64(def.branch)),
        ),
        combine(&ones, &both, &one, &one),
    ];
    [sphere, middle, cylinder].map(|row| row.map(|c| elevate(&c, DEGREE)))
}

fn build(
    def: &RollingSphereDefinition,
    nodes: &[Node],
    dec: &mut Decider,
) -> Result<(BSplineSurface, Vec<[Real; 4]>), Refusal> {
    let spans = nodes.len() - 1;
    let count = spans * DEGREE + 1;
    let mut controls = Vec::with_capacity(count * 3);
    let mut weights = Vec::with_capacity(count * 3);
    let mut homogeneous = Vec::with_capacity(count * 3);
    let mut knots = vec![q(-1); DEGREE + 1];
    for (index, pair) in nodes.windows(2).enumerate() {
        dec.check_deadline()?;
        let rows = span_rows(def, pair);
        let last = index + 1 == spans;
        for i in 0..DEGREE + usize::from(last) {
            for data in &rows {
                let w = &data[3][i];
                controls.push(def.place(
                    &(&data[0][i] / w),
                    &(&data[1][i] / w),
                    &(&data[2][i] / w),
                ));
                weights.push(w.clone());
                homogeneous.push(std::array::from_fn(|coordinate| {
                    data[coordinate][i].clone()
                }));
            }
        }
        let repeat = DEGREE + usize::from(last);
        knots.extend(std::iter::repeat_n(pair[1].u.clone(), repeat));
    }
    let v = [0, 0, 0, 1, 1, 1].map(q).to_vec();
    let carrier = BSplineSurface::new([DEGREE, 2], [knots, v], [count, 3], controls, weights, dec)?;
    Ok((carrier, homogeneous))
}

// Enclosed homogeneous controls in the sphere's frame, u-major. Rigid copies
// share one net, so no consumer pays to enclose a moved carrier again.
#[derive(Debug)]
struct Net {
    bits: u32,
    controls: Vec<[Ball; 4]>,
}

impl Net {
    fn new(controls: &[[Real; 4]], dec: &Decider) -> Result<Self, Refusal> {
        let bits = BITS.min(dec.policy().max_precision_bits);
        // These are the exact coefficients build() places in the world frame.
        // Enclose them before placement instead of undoing that placement with
        // interval arithmetic and multiplying the Euclidean controls by W again.
        let values: Vec<_> = controls.iter().flat_map(|point| point.iter()).collect();
        let enclosed = dec.enclose_all_by_operations_at(&values, bits)?;
        let controls = enclosed
            .chunks_exact(4)
            .map(|point| std::array::from_fn(|coordinate| point[coordinate].clone()))
            .collect();
        Ok(Self { bits, controls })
    }
}

// Proves the carrier itself regular, patch by patch, as for a saddle. With P
// the homogeneous point less `m W` for a near constant `m`, the derivatives
// over a patch are `D = P_x W - P W_x` up to the positive `W^2`, and a positive
// `(D_u x D_v) . n` for a fixed `n` keeps them independent.
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
            return Err(unresolved("subdivided band weight is not positive"));
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
        let [weight, along_u, across_v] = [
            value(3),
            cells([1, 0], &scales[0], 3),
            cells([0, 1], &scales[1], 3),
        ];
        let derivative = |step: [usize; 2], scale: &Rational, change: &Ball| -> [Ball; 3] {
            [0, 1, 2].map(|k| &(&cells(step, scale, k) * &weight) - &(&value(k) * change))
        };
        let u = derivative([1, 0], &scales[0], &along_u);
        let v = derivative([0, 1], &scales[1], &across_v);
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

/// A rolling-ball band over `[-1, 3] x [0, 1]`, from the sphere rail to the
/// cylinder rail, closed in `u`.
#[derive(Clone, Debug)]
pub struct RollingSphere {
    data: Arc<Data>,
}

#[derive(Debug)]
struct Data {
    definition: RollingSphereDefinition,
    nodes: Vec<Node>,
    carrier: BSplineSurface,
    net: Arc<Net>,
    bound: EvidenceBound,
    reversing: bool,
}

impl RollingSphere {
    /// Fits the fewest uniform spans, by doubling, whose recomputed bound meets
    /// the 100 nm target. A fit that meets only the 1 um cap is accepted and
    /// reported by its bound, and anything larger refuses.
    pub fn fit(definition: RollingSphereDefinition, dec: &mut Decider) -> Result<Self, Refusal> {
        let bits = BITS.min(dec.policy().max_precision_bits);
        let target = Rational::from_ratio(ROLLING_SPHERE_TARGET_NM, 1_000_000).expect("nonzero");
        let cap = Rational::from_ratio(HARD_CAP_NM, 1_000_000).expect("nonzero");
        let k = Constants::new(&definition, bits, dec)?;
        let mut samples = FitSamples::default();
        let mut spans = FIRST_SPANS;
        let (nodes, upper) = loop {
            let nodes = produce(&k, spans, bits, &mut samples, dec)?;
            let stop_above = (spans < LAST_SPANS).then_some(&target);
            match recompute_bound(&k, &nodes, bits, &mut samples, stop_above, dec) {
                Ok(Some(upper)) => break (nodes, upper),
                Err(error) if spans >= LAST_SPANS => return Err(error),
                _ => spans *= 2,
            }
        };
        if upper > cap {
            return Err(invalid("the band fit misses the 1 um approximation cap"));
        }
        check_regular(&definition, &nodes, bits, dec)?;
        let (carrier, homogeneous) = build(&definition, &nodes, dec)?;
        let net = Arc::new(Net::new(&homogeneous, dec)?);
        check_carrier(&net, &carrier, dec)?;
        Ok(Self {
            data: Arc::new(Data {
                definition,
                nodes,
                carrier,
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
    pub fn definition(&self) -> &RollingSphereDefinition {
        &self.data.definition
    }

    /// Two-sided distance bound between the carrier and the exact surface.
    pub fn bound(&self) -> &EvidenceBound {
        &self.data.bound
    }

    /// Number of Bezier spans along the band.
    pub fn span_count(&self) -> usize {
        self.data.nodes.len() - 1
    }

    /// The span ends along `u`.
    pub fn knots(&self) -> Vec<Rational> {
        self.data.nodes.iter().map(|node| node.u.clone()).collect()
    }

    /// The rails on the sphere and the cylinder, the carrier's own boundary
    /// rows, so a neighbour can share them exactly.
    pub fn rails(&self) -> [BSplineCurve; 2] {
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
                "no band query precision available",
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
                        .ok_or_else(|| unresolved("band weight is not positive"))
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

    // The rail a point lies on, with the knot span and exact parameter there.
    fn rail_of(
        &self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<(usize, usize, Real)>, Refusal> {
        let def = &self.data.definition;
        let rails = self.rails();
        for (which, rail) in rails.iter().enumerate() {
            let (gap, radius) = if which == 0 {
                (point.minus(&def.centre).norm_squared(), &def.sphere)
            } else {
                (def.radial(point).norm_squared(), &def.cylinder)
            };
            if math::sign(&(&gap - &radius.square()), dec)? != Sign::Zero {
                continue;
            }
            if let Some(at) = rail.parameter_of(point, dec)? {
                return Ok(Some((which, at.span(), at.value().clone())));
            }
        }
        Ok(None)
    }

    /// Exact normal of the carrier at a rail point, pointing to the ball's
    /// centre. Interior points have no exact chart parameters and refuse.
    ///
    /// The cylinder rail lies on the cylinder and the middle row lies in the
    /// cylinder's tangent plane at it, so both carrier derivatives there are
    /// square to the cylinder's radial, and that radial is the normal. The
    /// sphere rail and the middle row come from separately rounded data, so
    /// the carrier is only within its bound of tangent to the sphere there,
    /// and the normal is the rail's tangent crossed with the step to the
    /// middle row, which is the carrier's `v` derivative up to a positive
    /// factor.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        let def = &self.data.definition;
        let Some((which, span, at)) = self.rail_of(point, dec)? else {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "a band normal needs exact chart parameters, and only its rails have them",
            ));
        };
        let normal = if which == 0 {
            self.sphere_rail_normal(point, span, &at, dec)?
        } else {
            def.radial(point).scaled(&Real::from_i64(def.sides[1]))
        };
        Ok(if self.data.reversing {
            normal.scaled(&Real::from_i64(-1))
        } else {
            normal
        })
    }

    // The carrier's normal at a sphere rail point, from the exact span the
    // rail parameter lies in. Every node is a knot of full multiplicity, so a
    // span's controls are its own Bezier controls, and the rail and middle rows
    // of a span are each a Bezier curve. The sphere's radial is within the
    // bound of the normal, so it settles the orientation.
    fn sphere_rail_normal(
        &self,
        point: &Point3,
        span: usize,
        at: &Real,
        dec: &mut Decider,
    ) -> Result<Vector3, Refusal> {
        let def = &self.data.definition;
        let carrier = &self.data.carrier;
        let [_, rows] = carrier.control_count();
        let knots = carrier.knots()[0];
        let ends = [&knots[span], &knots[span + 1]];
        let span_knots: Vec<Rational> = ends
            .iter()
            .flat_map(|end| std::iter::repeat_n((*end).clone(), DEGREE + 1))
            .collect();
        let mut row = |which: usize| {
            let indices = (span - DEGREE..=span).map(|i| i * rows + which);
            BSplineCurve::new(
                DEGREE,
                span_knots.clone(),
                indices
                    .clone()
                    .map(|i| carrier.controls()[i].clone())
                    .collect(),
                indices.map(|i| carrier.weights()[i].clone()).collect(),
                dec,
            )
        };
        let rail = row(0)?;
        let middle = row(1)?;
        let tangent = rail.first_at(at, dec)?;
        let normal = tangent.cross(&middle.point_at(at, dec)?.minus(point));
        let radial = point
            .minus(&def.centre)
            .scaled(&Real::from_i64(def.sides[0]));
        match math::sign(&normal.dot(&radial), dec)? {
            Sign::Positive => Ok(normal),
            Sign::Negative => Ok(normal.scaled(&Real::from_i64(-1))),
            Sign::Zero => Err(invalid("the band is square to the sphere at its rail")),
        }
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
        range.ok_or_else(|| unresolved("a band has controls"))
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
            let measure = |patch: &Patch| -> Result<(Rational, Rational), Refusal> {
                let (mut top, mut best) = (None::<Rational>, None::<Rational>);
                let last = patch.net.len() - 1;
                for (index, h) in patch.net.iter().enumerate() {
                    let local = [0, 1, 2].map(|k| h[k].checked_div(&h[3]));
                    let [Some(x), Some(y), Some(z)] = local else {
                        return Err(unresolved("subdivided band weight is not positive"));
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
            let mut floor = floor.ok_or_else(|| unresolved("a band has patches"))?;
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
                return Err(unresolved("a band's range did not reach its width"));
            }
            answer.push(if lowest { [-top, -floor] } else { [floor, top] });
        }
        let [low, high]: [[Rational; 2]; 2] = answer.try_into().expect("one answer for each end");
        Ok([low, high])
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
                net: self.data.net.clone(),
                bound: self.data.bound.clone(),
                reversing: self.data.reversing != pose.reversing(),
            }),
        }
    }
}

// d/ds 2 atan(f) = 2 f' / (1 + f^2), valid through the jet's third coefficient.
fn turning(f: &Jet, bits: u32) -> Result<Jet, Refusal> {
    derivative(f, bits)
        .scale(&b(&q(2), bits))
        .div(&Jet::constant(b(&q(1), bits)).add(&f.square()))
}

// Row `row` of a span's Hermite fit at the chart parameter `s`.
fn fitted_row(s: &Jet, pair: &[Node], row: usize, bits: u32) -> Result<Jet, Refusal> {
    let (_, ends) = span_chart(pair);
    let data = [0, 1].map(|end| {
        [
            b(&pair[end].values[row], bits),
            b(&pair[end].slopes[row], bits),
        ]
    });
    hermite(
        s,
        [&ends[0], &ends[1]],
        [[&data[0][0], &data[0][1]], [&data[1][0], &data[1][1]]],
        bits,
    )
}

// A rail's chart on its support at `s` of a span: the angle's offset from the
// law's `t` with its slope jet, and the second coordinate with its slope jet.
// On the sphere rail the angle's offset is `atan` of the fitted contact's turn
// from the direction `t`, which stays inside a quarter turn because the axis
// passes inside the cylinder; the second coordinate is the latitude
// `branch (pi/2 - 2 atan rho)` with `rho^2 = sigma^2 + tau^2`. On the cylinder
// rail the angle is `t` itself and the second coordinate is the fitted height.
struct RailJets {
    offset: Option<Jet>,
    slope: Jet,
    second: Jet,
    rise: Jet,
}

fn rail_jets(
    k: &Constants,
    rail: usize,
    s: &Jet,
    pair: &[Node],
    bits: u32,
) -> Result<RailJets, Refusal> {
    let (chart, _) = span_chart(pair);
    let turn = b(&q(1 - 2 * chart), bits);
    let along = turning(s, bits)?;
    if rail == 1 {
        let height = fitted_row(s, pair, HEIGHT, bits)?;
        let rise = derivative(&height, bits);
        return Ok(RailJets {
            offset: None,
            slope: along,
            second: height,
            rise,
        });
    }
    let sigma = fitted_row(s, pair, SIGMA, bits)?;
    let tau = fitted_row(s, pair, TAU, bits)?;
    let one = k.c(&k.one);
    let square = s.square();
    let flat = one.sub(&square);
    let double = s.scale(&k.two);
    // The contact's components along and across the direction `t`, both over
    // the same positive factor `(1 + s^2) / turn^2`.
    let ahead = sigma.mul(&flat).add(&tau.mul(&double)).scale(&turn);
    let aside = tau.mul(&flat).sub(&sigma.mul(&double)).scale(&turn);
    if ahead.0[0].lo() <= Rational::zero() {
        return Err(unresolved(
            "the sphere rail does not stay within a quarter turn of its cylinder angle",
        ));
    }
    let lean = aside.div(&ahead)?;
    let half = k
        .one
        .checked_div(&k.two)
        .expect("two is separated from zero");
    let slope = along.add(&turning(&lean, bits)?.scale(&half));
    let rho = sigma.square().add(&tau.square()).sqrt()?;
    let rise = turning(&rho, bits)?.scale(&k.branch).neg();
    Ok(RailJets {
        offset: Some(lean),
        slope,
        second: rho,
        rise,
    })
}

// The jet of f' through its third coefficient.
fn derivative(f: &Jet, bits: u32) -> Jet {
    let zero = b(&Rational::zero(), bits);
    Jet(std::array::from_fn(|i| {
        if i == 4 {
            zero.clone()
        } else {
            &f.0[i + 1] * &b(&q(i as i64 + 1), bits)
        }
    }))
}

impl RollingSphere {
    /// Fits rail `rail`'s chart on its support within `tolerance` radians:
    /// `0` is the sphere rail and `1` the cylinder rail. The angle turns about
    /// the cylinder's axis from `across` toward `along` and runs from about
    /// `-pi/2` to `3 pi/2`; on the cylinder rail it is the law's `t`. The
    /// height is the cylinder rail's station along the axis from the sphere's
    /// centre, and on the sphere rail the latitude from the plane through the
    /// centre square to the axis, positive along the axis. A sphere rail that
    /// does not wind once round the axis, or a rail the jet cannot resolve
    /// within 2^12 pieces per span, refuses.
    pub fn parameter_table(
        &self,
        rail: usize,
        tolerance: &Rational,
        dec: &Decider,
    ) -> Result<ParameterTable, Refusal> {
        let nodes = &self.data.nodes;
        let bits = self.data.net.bits;
        let k = Constants::new(&self.data.definition, bits, dec)?;
        let pi = Ball::pi_checked(bits, dec)?;
        let quarter = &pi * &b(&Rational::from_ratio(1, 2).expect("nonzero"), bits);
        let mut pieces = 1i64;
        loop {
            dec.check_deadline()?;
            let mut worst = Rational::zero();
            let mut height = Rational::zero();
            for pair in nodes.windows(2) {
                let (_, ends) = span_chart(pair);
                let h = (&ends[1] - &ends[0]) / q(pieces);
                let h4 = h.pow(4).expect("nonzero span");
                for piece in 0..pieces {
                    dec.check_deadline()?;
                    let lo = &ends[0] + &(&h * &q(piece));
                    let hi = &lo + &h;
                    let s = Jet::linear(b(&lo, bits).hull(&b(&hi, bits)), k.one.clone());
                    let jets = rail_jets(&k, rail, &s, pair, bits)?;
                    // g'''' = 6 c3 and the cubic remainder is h^4 |g''''| / 384.
                    worst = worst.max(magnitude(&jets.slope.0[3]) * h4.clone() / q(64));
                    if rail == 0 {
                        height = height.max(magnitude(&jets.rise.0[3]) * h4.clone() / q(64));
                    }
                }
            }
            if (worst > *tolerance || height > *tolerance) && pieces < 4096 {
                pieces *= 2;
                continue;
            }
            if worst > *tolerance || height > *tolerance {
                return Err(unresolved("rail parameter table misses its tolerance"));
            }
            let mut table = ParameterTable {
                knots: Vec::new(),
                angle: Vec::new(),
                height: Vec::new(),
                bound: worst,
                height_bound: height,
            };
            let spans = nodes.len() - 1;
            for (index, pair) in nodes.windows(2).enumerate() {
                let (chart, ends) = span_chart(pair);
                let h = (&ends[1] - &ends[0]) / q(pieces);
                let turn = &pi * &b(&q(chart), bits);
                let last = index + 1 == spans;
                for piece in 0..pieces + i64::from(last) {
                    dec.check_deadline()?;
                    let x = &ends[0] + &(&h * &q(piece));
                    let s = Jet::linear(b(&x, bits), k.one.clone());
                    let jets = rail_jets(&k, rail, &s, pair, bits)?;
                    let mut value = &(&k.two * &s.0[0].atan_checked(dec)?) + &turn;
                    if let Some(lean) = &jets.offset {
                        value = &value + &lean.0[0].atan_checked(dec)?;
                    }
                    let second = if rail == 0 {
                        let polar = &k.two * &jets.second.0[0].atan_checked(dec)?;
                        &(&quarter - &polar) * &k.branch
                    } else {
                        jets.second.0[0].clone()
                    };
                    let [slope, rise] = [&jets.slope, &jets.rise].map(|jet| jet.0[0].clone());
                    table.bound = &table.bound + &(value.rad() + slope.rad());
                    table.height_bound = &table.height_bound + &(second.rad() + rise.rad());
                    table.angle.push([value.mid(), slope.mid()]);
                    table.height.push([second.mid(), rise.mid()]);
                    table.knots.push(&x + &q(2 * chart));
                }
            }
            return Ok(table);
        }
    }
}
