//! Checked rolling-ball envelope while a constant-radius ball rests on a
//! cylinder and on a parallel circular mouth edge.
//!
//! The ball centre follows the intersection of the cylinder's offset and the
//! tube around the mouth circle. Its envelope between the two contacts is not
//! rational, so the represented face is a genuine rational B-spline whose
//! boundary rows are constructed exactly on their supports: the cylinder rail,
//! the mouth circle and the two meridian arcs where the ball leaves the edge for
//! the plane. Only the middle row and two scalar reparametrisations are fitted.
//! The two-sided position bound is recomputed from the retained node data and
//! interval Taylor jets of the exact definition, never from producer residuals.

use crate::{BSplineSurface, Isometry, Point3, Pose, Rotation, Vector3, math};
use catallax_numbers::{
    Ball, Decider, EvidenceBound, EvidenceUnit, Outcome, Rational, Real, Refusal, Sign, StrictSign,
};
use std::sync::Arc;

fn refusal(detail: impl Into<String>) -> Refusal {
    Refusal::new("shapes.rolling-envelope", detail)
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

// Coefficients are Taylor derivatives divided by factorial, over the entire
// supplied parameter ball. Order four is what the Hermite remainder needs.
#[derive(Clone)]
struct Jet([Ball; 5]);
impl Jet {
    fn constant(value: Ball) -> Self {
        let bits = value.precision();
        Self(std::array::from_fn(|i| {
            if i == 0 {
                value.clone()
            } else {
                b(&Rational::zero(), bits)
            }
        }))
    }
    fn linear(value: Ball, slope: Ball) -> Self {
        let mut result = Self::constant(value);
        result.0[1] = slope;
        result
    }
    fn add(&self, other: &Self) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] + &other.0[i]))
    }
    fn neg(&self) -> Self {
        Self(std::array::from_fn(|i| -&self.0[i]))
    }
    fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }
    fn scale(&self, factor: &Ball) -> Self {
        Self(std::array::from_fn(|i| &self.0[i] * factor))
    }
    fn mul(&self, other: &Self) -> Self {
        let bits = self.0[0].precision();
        Self(std::array::from_fn(|n| {
            (0..=n).fold(b(&Rational::zero(), bits), |sum, i| {
                sum + &self.0[i] * &other.0[n - i]
            })
        }))
    }
    fn square(&self) -> Self {
        let mut result = self.mul(self);
        result.0[0] = self.0[0].square();
        result
    }
    fn reciprocal(&self) -> Result<Self, Refusal> {
        let bits = self.0[0].precision();
        let zero = b(&Rational::zero(), bits);
        let mut result: [Ball; 5] = std::array::from_fn(|_| zero.clone());
        result[0] = b(&Rational::one(), bits)
            .checked_div(&self.0[0])
            .ok_or_else(|| unresolved("envelope chart divisor is not separated from zero"))?;
        for n in 1..5 {
            let sum = (1..=n).fold(zero.clone(), |sum, i| sum + &self.0[i] * &result[n - i]);
            result[n] = (-sum)
                .checked_div(&self.0[0])
                .ok_or_else(|| unresolved("envelope chart derivative divisor is unresolved"))?;
        }
        Ok(Self(result))
    }
    fn div(&self, other: &Self) -> Result<Self, Refusal> {
        Ok(self.mul(&other.reciprocal()?))
    }
    fn sqrt(&self) -> Result<Self, Refusal> {
        let bits = self.0[0].precision();
        if self.0[0].lo() <= Rational::zero() {
            return Err(unresolved(
                "envelope chart requires a positive square-root lower bound",
            ));
        }
        let zero = b(&Rational::zero(), bits);
        let mut result: [Ball; 5] = std::array::from_fn(|_| zero.clone());
        result[0] = self.0[0]
            .sqrt()
            .ok_or_else(|| unresolved("undefined envelope chart square root"))?;
        let denominator = &b(&q(2), bits) * &result[0];
        for n in 1..5 {
            let sum = (1..n).fold(zero.clone(), |sum, i| sum + &result[i] * &result[n - i]);
            result[n] = (&self.0[n] - &sum)
                .checked_div(&denominator)
                .ok_or_else(|| unresolved("envelope square-root derivative is unresolved"))?;
        }
        Ok(Self(result))
    }
    /// Value and first derivative, the Hermite data at a point ball.
    fn data(&self) -> [Ball; 2] {
        [self.0[0].clone(), self.0[1].clone()]
    }
}

/// Exact definition of the envelope of a ball of radius `r` that touches a
/// cylinder of radius `body` from outside and the edge of a mouth circle of
/// radius `mouth` whose centre lies `offset` from the cylinder axis, in the
/// plane through `origin` perpendicular to `axis`.
///
/// `axis` points to the ball side of that plane and `toward` from the cylinder
/// axis to the mouth centre; `across` completes the frame. The mouth edge's
/// second face has outward normal `-edge[0] * radial + edge[1] * axis` in the
/// mouth's meridian half-plane, so the ball stays on the edge, not that face.
#[derive(Clone, Debug)]
pub struct RollingEnvelopeDefinition {
    origin: Point3,
    axis: Vector3,
    toward: Vector3,
    across: Vector3,
    body: Real,
    mouth: Real,
    offset: Real,
    radius: Real,
    edge: [Real; 2],
}

// Exact transition constants: the ball leaves the mouth edge for the plane at
// cylinder angle +-theta0, with t0 = tan(theta0 / 2). `square` is t0^2 and
// `mouth` the mouth reparametrisation over t0 at the negative transition, both
// free of the square root in `z`.
#[derive(Clone, Debug)]
struct Transition {
    centre: Real,
    cos: Real,
    x: Real,
    z: Real,
    t: Real,
    square: Real,
    mouth: Real,
}

impl RollingEnvelopeDefinition {
    /// Validates an orthonormal frame and the exact rim-contact regime: a
    /// mouth outside the cylinder that the offset cylinder crosses on its near
    /// side, with the ball strictly inside the edge's normal cone between the
    /// two transitions.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        origin: Point3,
        axis: Vector3,
        toward: Vector3,
        body: Real,
        mouth: Real,
        offset: Real,
        radius: Real,
        edge: [Real; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        math::validate_point(&origin, dec)?;
        math::validate_vector(&axis, dec)?;
        math::validate_vector(&toward, dec)?;
        for value in [&body, &mouth, &offset, &radius, &edge[0], &edge[1]] {
            math::validate(value, dec)?;
        }
        let one = Real::one();
        for (value, what) in [
            (
                axis.norm_squared() - one.clone(),
                "envelope axis must be a unit vector",
            ),
            (
                toward.norm_squared() - one.clone(),
                "mouth direction must be a unit vector",
            ),
            (
                axis.dot(&toward),
                "mouth direction must be perpendicular to the axis",
            ),
        ] {
            if math::sign(&value, dec)? != Sign::Zero {
                return Err(refusal(what));
            }
        }
        let positive = |value: &Real, what: &str, dec: &mut Decider| -> Result<(), Refusal> {
            if math::sign(value, dec)? != Sign::Positive {
                return Err(refusal(what.to_owned()));
            }
            Ok(())
        };
        positive(&body, "cylinder radius must be positive", dec)?;
        positive(&mouth, "mouth radius must be positive", dec)?;
        positive(&radius, "ball radius must be positive", dec)?;
        positive(
            &edge[0],
            "edge face normal must lean toward the mouth axis",
            dec,
        )?;
        positive(
            &edge[1],
            "edge face normal must lean toward the ball side",
            dec,
        )?;
        let centre = &body + &radius;
        positive(
            &(&(&offset - &mouth) - &body),
            "the mouth must lie wholly outside the cylinder",
            dec,
        )?;
        positive(
            &(&centre - &(&offset - &mouth)),
            "the ball centre path must reach the mouth edge",
            dec,
        )?;
        positive(
            &(&offset - &centre),
            "the rim contact must stay on the near side of the mouth axis",
            dec,
        )?;
        // The edge-face margin p*sqrt(r^2-g^2)-q*g decreases in the gap g=m-rho,
        // which is largest at the station nearest the cylinder.
        let gap = &mouth - &(&offset - &centre);
        let height = (&radius.square() - &gap.square()).sqrt();
        positive(
            &(&(&edge[0] * &height) - &(&edge[1] * &gap)),
            "the ball would leave the mouth edge for its second face",
            dec,
        )?;
        let across = toward.cross(&axis);
        Ok(Self {
            origin,
            axis,
            toward,
            across,
            body,
            mouth,
            offset,
            radius,
            edge,
        })
    }
    /// Point on the cylinder axis in the mouth plane.
    pub fn origin(&self) -> &Point3 {
        &self.origin
    }
    /// Unit cylinder axis, pointing to the ball side of the mouth plane.
    pub fn axis(&self) -> &Vector3 {
        &self.axis
    }
    /// Unit direction from the cylinder axis to the mouth centre.
    pub fn toward(&self) -> &Vector3 {
        &self.toward
    }
    /// Third frame vector; increasing chart parameter `u` follows it.
    pub fn across(&self) -> &Vector3 {
        &self.across
    }
    /// Cylinder radius.
    pub fn body(&self) -> &Real {
        &self.body
    }
    /// Mouth circle radius.
    pub fn mouth(&self) -> &Real {
        &self.mouth
    }
    /// Distance from the cylinder axis to the mouth centre.
    pub fn offset(&self) -> &Real {
        &self.offset
    }
    /// Ball radius.
    pub fn radius(&self) -> &Real {
        &self.radius
    }
    /// Radial and axial parts of the mouth edge's second outward face normal.
    pub fn edge(&self) -> &[Real; 2] {
        &self.edge
    }
    pub(crate) fn mapped(&self, pose: &Isometry) -> Self {
        Self {
            origin: pose.apply_point(&self.origin),
            axis: pose.apply_vector(&self.axis),
            toward: pose.apply_vector(&self.toward),
            across: pose.apply_vector(&self.across),
            ..self.clone()
        }
    }
    fn transition(&self) -> Transition {
        let centre = &self.body + &self.radius;
        let cos = &(&(&centre.square() + &self.offset.square()) - &self.mouth.square())
            / &(&(&Real::from_i64(2) * &self.offset) * &centre);
        let x = &centre * &cos;
        let z = (&centre.square() - &x.square()).sqrt();
        let t = &z / &(&centre + &x);
        let square = &(&centre - &x) / &(&centre + &x);
        let mouth = &(&centre + &x) / &(&(&self.mouth + &self.offset) - &x);
        Transition {
            centre,
            cos,
            x,
            z,
            t,
            square,
            mouth,
        }
    }
    fn world(&self, local: [&Real; 3]) -> Point3 {
        self.origin.plus(
            &self
                .toward
                .scaled(local[0])
                .plus(&self.axis.scaled(local[1]))
                .plus(&self.across.scaled(local[2])),
        )
    }
    fn local(&self, point: &Point3) -> [Real; 3] {
        let d = point.minus(&self.origin);
        [d.dot(&self.toward), d.dot(&self.axis), d.dot(&self.across)]
    }
    // The carrier control and weight of one homogeneous row entry. Origin and
    // frame apply to each part, so rational parts stay single leaves.
    fn control(&self, row: &[Split; 4], i: usize, t: &Real) -> (Point3, Real) {
        let o = &self.origin;
        let frame = [
            [&o.x, &self.toward.x, &self.axis.x, &self.across.x],
            [&o.y, &self.toward.y, &self.axis.y, &self.across.y],
            [&o.z, &self.toward.z, &self.axis.z, &self.across.z],
        ];
        let weight = row[3].at(i, t);
        let [x, y, z] = frame.map(|[origin, toward, axis, across]| {
            let part = |c: [&Real; 4]| {
                [(origin, c[3]), (toward, c[0]), (axis, c[1]), (across, c[2])]
                    .iter()
                    .fold(Real::zero(), |sum, (f, v)| plus(&sum, &times(f, v)))
            };
            let even = part([0, 1, 2, 3].map(|c| &row[c].even[i]));
            let odd = part([0, 1, 2, 3].map(|c| &row[c].odd[i]));
            &join(&even, &odd, t) / &weight
        });
        (Point3::new(x, y, z), weight)
    }
    /// Exact transition corners: cylinder foot and mouth/plane foot at `sign`.
    pub fn corners(&self, sign: Sign) -> [Point3; 2] {
        let t = self.transition();
        let z = if sign == Sign::Negative {
            -&t.z
        } else {
            t.z.clone()
        };
        let scale = &self.body / &t.centre;
        [
            self.world([&(&t.x * &scale), &self.radius, &(&z * &scale)]),
            self.world([&t.x, &Real::zero(), &z]),
        ]
    }
}

// Ball enclosures of the definition constants at one working precision.
struct Constants {
    centre: Ball,
    mouth: Ball,
    offset: Ball,
    radius: Ball,
    t: Ball,
    one: Ball,
    two: Ball,
}

impl Constants {
    fn new(
        definition: &RollingEnvelopeDefinition,
        bits: u32,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let t = definition.transition();
        let e = |value: &Real| dec.enclose_at(value, bits);
        Ok(Self {
            centre: e(&t.centre)?,
            mouth: e(&definition.mouth)?,
            offset: e(&definition.offset)?,
            radius: e(&definition.radius)?,
            t: e(&t.t)?,
            one: b(&Rational::one(), bits),
            two: b(&q(2), bits),
        })
    }
    fn c(&self, value: &Ball) -> Jet {
        Jet::constant(value.clone())
    }
}

// Exact rows of the envelope at chart parameter s, as jets in s. `fitted` is
// the represented mouth reparametrisation, which fixes the mouth row weight
// and therefore the exact middle weight. The apex is kept in the cylinder's
// tangent plane at the rail point: `tangent` along the rail circle's tangent
// and `axial` along the axis, so the represented middle row leaves the rail
// tangentially and the carrier stays outside the cylinder.
struct Rows {
    height: Jet,
    mouth: Jet,
    tangent: Jet,
    axial: Jet,
    weight: Option<Jet>,
    cosine: Jet,
}

fn rows(k: &Constants, s: &Jet, fitted: Option<&Jet>) -> Result<Rows, Refusal> {
    let one = k.c(&k.one);
    let u = s.scale(&k.t);
    let u2 = u.square();
    let den = one.add(&u2);
    let cos = one.sub(&u2).div(&den)?;
    let sin = u.scale(&k.two).div(&den)?;
    let cx = cos.scale(&k.centre);
    let cz = sin.scale(&k.centre);
    let d = &k.offset;
    let rho = k
        .c(&(&k.centre.square() + &d.square()))
        .sub(&cos.scale(&(&(&k.two * d) * &k.centre)))
        .sqrt()?;
    let gap = k.c(&k.mouth).sub(&rho);
    let height = k.c(&k.radius.square()).sub(&gap.square()).sqrt()?;
    let rx = k.c(d).add(&cx.sub(&k.c(d)).scale(&k.mouth).div(&rho)?);
    let rz = cz.scale(&k.mouth).div(&rho)?;
    let mouth = cz.neg().div(&rho.add(&k.c(d)).sub(&cx))?;
    let r_inv = b(&Rational::one(), k.one.precision())
        .checked_div(&k.radius)
        .ok_or_else(|| unresolved("ball radius is not separated from zero"))?;
    // a = -(cos, 0, sin) and c = (R - C)/r are the unit contact directions.
    let dot = cos
        .mul(&rx.sub(&cx))
        .add(&sin.mul(&rz.sub(&cz)))
        .scale(&r_inv)
        .neg();
    let one_d = one.add(&dot);
    // Apex M = C + (r a + R - C) / (1 + d), relative to the rail point A.
    let tangent = sin
        .mul(&rx.sub(&cx))
        .neg()
        .add(&cos.mul(&rz.sub(&cz)))
        .div(&one_d)?;
    let axial = height.div(&one_d)?.neg();
    // Middle weight cos(half) sqrt(w0 w2), divided by the rail weight 1 + u^2
    // that the tangent-plane form of the middle row multiplies back.
    let weight = match fitted {
        None => None,
        Some(w) => Some(
            one_d
                .mul(&one.add(&w.square()))
                .div(&den.scale(&k.two))?
                .sqrt()?,
        ),
    };
    let cosine = one_d
        .scale(
            &b(&Rational::one(), k.one.precision())
                .checked_div(&k.two)
                .expect("two is separated from zero"),
        )
        .sqrt()?;
    Ok(Rows {
        height,
        mouth,
        tangent,
        axial,
        weight,
        cosine,
    })
}

/// Hermite data of the fitted scalar rows at one chart node.
#[derive(Clone, Debug)]
struct Node {
    s: Rational,
    // Value then s-derivative: rail height, mouth reparametrisation, apex
    // tangent-plane offsets and middle weight, all in the definition frame.
    values: [Real; ROWS],
    slopes: [Rational; ROWS],
}

const ROWS: usize = 5;
// Enclosed value and slope of each fitted row at one node.
type RowData = [[Ball; 2]; ROWS];

const HEIGHT: usize = 0;
const MOUTH: usize = 1;
const TANGENT: usize = 2;
const AXIAL: usize = 3;
const WEIGHT: usize = 4;

// The Hermite cubic on one span, in exact Taylor form about its midpoint. The
// cardinal-basis form would lose cancellation over an interval and inflate
// every higher Taylor coefficient by inverse powers of the span width.
fn hermite(
    s: &Jet,
    span: [&Rational; 2],
    data: [[&Ball; 2]; 2],
    bits: u32,
) -> Result<Jet, Refusal> {
    let h = span[1] - span[0];
    let inverse = h.recip().ok_or_else(|| refusal("empty envelope span"))?;
    let c = |v: Rational| b(&v, bits);
    let [[v0, d0], [v1, d1]] = data;
    let hb = c(h.clone());
    let a0 = &(&(v0 + v1) * &c(Rational::from_ratio(1, 2).expect("nonzero")))
        + &(&(&hb * &(d0 - d1)) * &c(Rational::from_ratio(1, 8).expect("nonzero")));
    let a1 = &(&(v1 - v0) * &c(q(3) * inverse.clone() / q(2)))
        - &(&(d0 + d1) * &c(Rational::from_ratio(1, 4).expect("nonzero")));
    let a2 = &(d1 - d0) * &c(inverse.clone() / q(2));
    let a3 = &(&(&(v0 - v1) * &c(q(2))) + &(&hb * &(d0 + d1)))
        * &c(inverse.pow(3).expect("nonzero span"));
    let sigma = s.sub(&Jet::constant(c((span[0] + span[1]) / q(2))));
    // Horner form: ((a3 sigma + a2) sigma + a1) sigma + a0.
    let mut value = Jet::constant(a3);
    for coefficient in [a2, a1, a0] {
        value = value.mul(&sigma).add(&Jet::constant(coefficient));
    }
    Ok(value)
}

// Bernstein-form polynomial arithmetic on exact coefficients.
fn binomial(n: usize, k: usize) -> Rational {
    let mut value = Rational::one();
    for i in 0..k {
        value = value * q((n - i) as i64) / q((i + 1) as i64);
    }
    value
}
fn is_zero(value: &Real) -> bool {
    value.rational_ref().is_some_and(Rational::is_zero)
}
fn plus(a: &Real, c: &Real) -> Real {
    if is_zero(a) {
        c.clone()
    } else if is_zero(c) {
        a.clone()
    } else {
        a + c
    }
}
fn times(a: &Real, c: &Real) -> Real {
    if is_zero(a) || is_zero(c) {
        Real::zero()
    } else {
        a * c
    }
}
fn join(even: &Real, odd: &Real, t: &Real) -> Real {
    plus(even, &times(t, odd))
}
fn product(a: &[Real], c: &[Real]) -> Vec<Real> {
    let (m, n) = (a.len() - 1, c.len() - 1);
    (0..=m + n)
        .map(|k| {
            let mut sum = Real::zero();
            for i in k.saturating_sub(n)..=k.min(m) {
                let factor = binomial(m, i) * binomial(n, k - i) / binomial(m + n, k);
                sum = plus(
                    &sum,
                    &times(&times(&a[i], &c[k - i]), &Real::from_rational(factor)),
                );
            }
            sum
        })
        .collect()
}
fn elevate(a: &[Real], degree: usize) -> Vec<Real> {
    let mut a = a.to_vec();
    while a.len() - 1 < degree {
        let n = a.len();
        let mut next = Vec::with_capacity(n + 1);
        next.push(a[0].clone());
        for i in 1..n {
            let left = Real::from_rational(q(i as i64) / q(n as i64));
            let right = Real::from_rational(Rational::one() - q(i as i64) / q(n as i64));
            next.push(plus(&times(&a[i - 1], &left), &times(&a[i], &right)));
        }
        next.push(a[n - 1].clone());
        a = next;
    }
    a
}
fn combine(a: &[Real], c: &[Real], sa: &Real, sc: &Real) -> Vec<Real> {
    let degree = (a.len().max(c.len())) - 1;
    let (a, c) = (elevate(a, degree), elevate(c, degree));
    a.iter()
        .zip(&c)
        .map(|(x, y)| plus(&times(x, sa), &times(y, sc)))
        .collect()
}

// One span coordinate in Bernstein form, as even + t * odd for the transition
// tangent t. With rational node data t is the only irrational in a span and t^2
// is rational, so keeping the parts apart lets Real's exact cache collapse each
// to one leaf instead of retaining a deep expression over t in every control.
#[derive(Clone)]
struct Split {
    even: Vec<Real>,
    odd: Vec<Real>,
}
impl Split {
    fn even(even: Vec<Real>) -> Self {
        let odd = vec![Real::zero(); even.len()];
        Self { even, odd }
    }
    fn odd(odd: Vec<Real>) -> Self {
        let even = vec![Real::zero(); odd.len()];
        Self { even, odd }
    }
    fn scaled(&self, factor: &Real) -> Self {
        let scale = |p: &[Real]| p.iter().map(|x| times(x, factor)).collect();
        Self {
            even: scale(&self.even),
            odd: scale(&self.odd),
        }
    }
    fn combine(&self, other: &Self, sa: &Real, sc: &Real) -> Self {
        Self {
            even: combine(&self.even, &other.even, sa, sc),
            odd: combine(&self.odd, &other.odd, sa, sc),
        }
    }
    fn product(&self, other: &Self, square: &Real) -> Self {
        let one = Real::one();
        Self {
            even: combine(
                &product(&self.even, &other.even),
                &product(&self.odd, &other.odd),
                &one,
                square,
            ),
            odd: combine(
                &product(&self.even, &other.odd),
                &product(&self.odd, &other.even),
                &one,
                &one,
            ),
        }
    }
    fn elevated(&self, degree: usize) -> Self {
        Self {
            even: elevate(&self.even, degree),
            odd: elevate(&self.odd, degree),
        }
    }
    fn at(&self, i: usize, t: &Real) -> Real {
        join(&self.even[i], &self.odd[i], t)
    }
    fn values(&self, t: &Real) -> Vec<Real> {
        (0..self.even.len()).map(|i| self.at(i, t)).collect()
    }
}

// Station parametrisation of one span, [1 + u^2, 1 - u^2, 2u] for u = t s:
// the exact cosine and sine are the last two over the first.
fn station(t: &Transition, span: [&Rational; 2]) -> [Split; 3] {
    let one = Real::one();
    let u = Split::odd(span.map(|s| Real::from_rational(s.clone())).to_vec());
    let u2 = u.product(&u, &t.square);
    let ones = Split::even(vec![one.clone(); 3]);
    [
        ones.combine(&u2, &one, &one),
        ones.combine(&u2, &one, &(-&one)),
        u.elevated(2).scaled(&Real::from_i64(2)),
    ]
}
fn cubic(span: [&Rational; 2], values: [&Real; 2], slopes: [&Rational; 2]) -> Vec<Real> {
    let third = Real::from_rational((span[1] - span[0]) / q(3));
    vec![
        values[0].clone(),
        values[0] + &(&third * &Real::from_rational(slopes[0].clone())),
        values[1] - &(&third * &Real::from_rational(slopes[1].clone())),
        values[1].clone(),
    ]
}

const DEGREE: usize = 8;

/// A finite checked rolling-ball envelope face. The carrier is the represented
/// geometry over `[-1,1] x [0,1]`: `u` runs along the ball path from the
/// negative to the positive transition, `v` from the cylinder rail to the mouth.
#[derive(Clone, Debug)]
pub struct RollingEnvelope {
    data: Arc<Data>,
}

#[derive(Debug)]
struct Data {
    definition: RollingEnvelopeDefinition,
    nodes: Vec<Node>,
    carrier: BSplineSurface,
    bound: EvidenceBound,
    reversing: bool,
}

/// Kernel-made approximation target from ARCHITECTURE's composable-error policy.
pub const ROLLING_ENVELOPE_TARGET_NM: i64 = 100;
const HARD_CAP_NM: i64 = 1_000;
const MAX_SPANS: usize = 128;

impl RollingEnvelope {
    /// Fits the fewest uniform spans, by doubling from eight, whose recomputed
    /// bound meets the 100 nm target. A fit that meets only the 1 um hard cap is
    /// accepted and reported by its bound; anything larger refuses.
    pub fn fit(definition: RollingEnvelopeDefinition, dec: &mut Decider) -> Result<Self, Refusal> {
        let target = Rational::from_ratio(ROLLING_ENVELOPE_TARGET_NM, 1_000_000).expect("nonzero");
        let cap = Rational::from_ratio(HARD_CAP_NM, 1_000_000).expect("nonzero");
        let mut spans = 8;
        let (nodes, upper) = loop {
            let nodes = produce(&definition, spans, dec)?;
            // A coarse cover may leave a span's interval jet unresolved; only
            // the finest attempted cover's refusal is final.
            match recompute_bound(&definition, &nodes, dec) {
                Ok(upper) if upper <= target || spans >= MAX_SPANS => break (nodes, upper),
                Err(error) if spans >= MAX_SPANS => return Err(error),
                _ => spans *= 2,
            }
        };
        if upper > cap {
            return Err(refusal(
                "the envelope fit misses the 1 um kernel approximation cap",
            ));
        }
        let carrier = build(&definition, &nodes, dec)?;
        check_transitions(&definition, &carrier, dec)?;
        check_regular(&definition, &carrier, dec)?;
        check_clearance(&definition, &nodes, dec)?;
        Ok(Self {
            data: Arc::new(Data {
                definition,
                nodes,
                carrier,
                bound: EvidenceBound {
                    upper,
                    unit: EvidenceUnit::Millimetres,
                },
                reversing: false,
            }),
        })
    }
    /// The same certified face for a definition that differs only by a
    /// rotation about the shared cylinder axis, such as a patterned hole. The
    /// copy is rigid transport, so its bound and proofs carry over unchanged.
    pub fn reframed(
        &self,
        definition: RollingEnvelopeDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let own = &self.data.definition;
        let same = |a: &Real, c: &Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(math::sign(&(a - c), dec)? == Sign::Zero)
        };
        let same_vector = |a: &Vector3, c: &Vector3, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(same(&a.x, &c.x, dec)? && same(&a.y, &c.y, dec)? && same(&a.z, &c.z, dec)?)
        };
        let scalars = [
            (&own.body, &definition.body),
            (&own.mouth, &definition.mouth),
            (&own.offset, &definition.offset),
            (&own.radius, &definition.radius),
            (&own.edge[0], &definition.edge[0]),
            (&own.edge[1], &definition.edge[1]),
        ];
        for (a, c) in scalars {
            if !same(a, c, dec)? {
                return Err(refusal("a reframed envelope needs the same contact law"));
            }
        }
        if !same_vector(&own.axis, &definition.axis, dec)?
            || !same_vector(&own.origin.minus(&definition.origin), &Vector3::zero(), dec)?
        {
            return Err(refusal("a reframed envelope must share the cylinder axis"));
        }
        // The half-angle quaternion (1 + t0.t1, t0 x t1) turns t0 onto t1 about
        // the shared axis; opposite directions use a half turn about the axis.
        let mut w = &Real::one() + &own.toward.dot(&definition.toward);
        let mut v = own.toward.cross(&definition.toward);
        if math::sign(&w, dec)? == Sign::Zero {
            w = Real::zero();
            v = own.axis.clone();
        }
        let rotation = Rotation::from_real_quaternion(w, v.x, v.y, v.z, dec)?;
        let centre = own.origin.minus(&Point3::origin());
        let translation = centre.minus(&rotation.apply_real(&centre));
        let copy = self.transformed(&Pose::new(rotation, translation, dec)?);
        let moved = &copy.data.definition;
        if !same_vector(&moved.toward, &definition.toward, dec)?
            || !same_vector(&moved.across, &definition.across, dec)?
            || !same_vector(
                &moved.origin.minus(&definition.origin),
                &Vector3::zero(),
                dec,
            )?
        {
            return Err(refusal("reframing did not reproduce the requested frame"));
        }
        Ok(copy)
    }
    /// Represented geometry over `[-1,1] x [0,1]`.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }
    /// Retained exact envelope definition.
    pub fn definition(&self) -> &RollingEnvelopeDefinition {
        &self.data.definition
    }
    /// Proved two-way position bound between carrier and exact envelope, in millimetres.
    pub fn bound(&self) -> &EvidenceBound {
        &self.data.bound
    }
    /// Number of chart spans along the ball path.
    pub fn span_count(&self) -> usize {
        self.data.nodes.len() - 1
    }
    /// Positive identity of the immutable checked definition.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
    /// Rigid transport of the carrier and frame; the bound is unchanged.
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
                reversing: self.data.reversing != pose.reversing(),
            }),
        }
    }

    /// Exact incidence on a boundary built exactly on its support: the
    /// cylinder rail, the mouth arc or a transition meridian. `Ok(false)` only
    /// says the point is on none of them.
    pub fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        dec.check_deadline()?;
        math::validate_point(point, dec)?;
        let def = &self.data.definition;
        let [x, y, z] = def.local(point);
        let t = def.transition();
        let zero = |value: Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(math::sign(&value, dec)? == Sign::Zero)
        };
        let within = |value: &Real, limit: &Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(math::sign(&(&value.square() - &limit.square()), dec)? != Sign::Positive)
        };
        // Mouth arc between the transitions, on the side nearer the cylinder.
        if zero(y.clone(), dec)?
            && zero(
                &(&(&x - &def.offset).square() + &z.square()) - &def.mouth.square(),
                dec,
            )?
            && within(&z, &t.z, dec)?
            && math::sign(&(&x - &t.x), dec)? != Sign::Positive
        {
            return Ok(true);
        }
        // Cylinder rail: the represented height at the exact chart parameter.
        let rail_z = &(&def.body * &t.z) / &t.centre;
        if zero(&(&x.square() + &z.square()) - &def.body.square(), dec)?
            && math::sign(&x, dec)? == Sign::Positive
            && within(&z, &rail_z, dec)?
        {
            let s = &(&z / &(&def.body + &x)) / &t.t;
            if let Some(height) = self.rail_height(&s, dec)? {
                return zero(&y - &height, dec);
            }
        }
        Ok(self.meridian(&[x, y, z], dec)?.is_some())
    }

    // The transition meridian through a local point, as the signed z of its
    // ball centre: a quarter arc from the cylinder foot to the plane foot.
    fn meridian(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Real>, Refusal> {
        let def = &self.data.definition;
        let t = def.transition();
        let [x, y, z] = local;
        for z_c in [t.z.clone(), -&t.z] {
            let sin = &z_c / &t.centre;
            let radial = &(x * &t.cos) + &(z * &sin);
            let tangential = &(z * &t.cos) - &(x * &sin);
            if math::sign(&tangential, dec)? == Sign::Zero
                && math::sign(
                    &(&(&(&radial - &t.centre).square() + &(y - &def.radius).square())
                        - &def.radius.square()),
                    dec,
                )? == Sign::Zero
                && math::sign(&(&radial - &t.centre), dec)? != Sign::Positive
                && math::sign(&(y - &def.radius), dec)? != Sign::Positive
            {
                return Ok(Some(z_c));
            }
        }
        Ok(None)
    }

    fn rail_height(&self, s: &Real, dec: &mut Decider) -> Result<Option<Real>, Refusal> {
        for pair in self.data.nodes.windows(2) {
            let lo = Real::from_rational(pair[0].s.clone());
            let hi = Real::from_rational(pair[1].s.clone());
            if math::sign(&(s - &lo), dec)? == Sign::Negative
                || math::sign(&(s - &hi), dec)? == Sign::Positive
            {
                continue;
            }
            let span = [&pair[0].s, &pair[1].s];
            let control = cubic(
                span,
                [&pair[0].values[HEIGHT], &pair[1].values[HEIGHT]],
                [&pair[0].slopes[HEIGHT], &pair[1].slopes[HEIGHT]],
            );
            let h = Real::from_rational(span[1] - span[0]);
            let tau = &(s - &lo) / &h;
            let rest = &Real::one() - &tau;
            let mut value = Real::zero();
            for (i, c) in control.iter().enumerate() {
                let weight = &(&Real::from_rational(binomial(3, i)) * &tau.pow_u32(i as u32))
                    * &rest.pow_u32((3 - i) as u32);
                value = &value + &(c * &weight);
            }
            return Ok(Some(value));
        }
        Ok(None)
    }

    /// Exact normal at a transition-meridian point, pointing to the ball side.
    /// Interior points have no exact chart parameters and refuse.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        math::validate_point(point, dec)?;
        let def = &self.data.definition;
        let Some(z_c) = self.meridian(&def.local(point), dec)? else {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "envelope normal needs exact chart parameters; only transition meridians have them",
            ));
        };
        let t = def.transition();
        let normal = def.world([&t.x, &def.radius, &z_c]).minus(point);
        Ok(if self.data.reversing {
            normal.scaled(&Real::from_i64(-1))
        } else {
            normal
        })
    }

    /// Encloses every chart preimage of a spatial box by pruned Bezier
    /// subdivision. Existence still requires a separate incidence proof.
    pub fn parameters_of_box(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        dec.check_deadline()?;
        let bits = point
            .iter()
            .map(Ball::precision)
            .max()
            .expect("three coordinates")
            .max(64)
            .min(dec.policy().max_precision_bits);
        let mut patches = Patch::all(&self.data.carrier, bits, dec)?;
        for _ in 0..24 {
            dec.check_deadline()?;
            let mut next = Vec::with_capacity(patches.len() * 2);
            for patch in patches {
                let hull = patch.hull()?;
                if (0..3).all(|i| meets(&hull[i], &point[i])) {
                    next.extend(patch.split());
                }
            }
            if next.is_empty() {
                return Ok(None);
            }
            if next.len() > 4096 {
                return Err(unresolved("envelope preimage cover does not localise"));
            }
            patches = next;
        }
        let mut hull: Option<[[Rational; 2]; 2]> = None;
        for patch in patches {
            let cell = patch.cell;
            hull = Some(match hull {
                None => cell,
                Some(h) => [
                    [
                        h[0][0].clone().min(cell[0][0].clone()),
                        h[0][1].clone().max(cell[0][1].clone()),
                    ],
                    [
                        h[1][0].clone().min(cell[1][0].clone()),
                        h[1][1].clone().max(cell[1][1].clone()),
                    ],
                ],
            });
        }
        Ok(hull.map(|h| h.map(|[lo, hi]| b(&lo, bits).hull(&b(&hi, bits)))))
    }

    /// Oriented side, positive toward the ball centres unless transported by a
    /// reflection. Only points whose projection along the certified chart
    /// direction lands on the finite face have a side certificate.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        if self.on_boundary(point, dec)? {
            return Ok(Sign::Zero);
        }
        let bits = 128.min(dec.policy().max_precision_bits);
        let def = &self.data.definition;
        let p = math::point(point, bits, dec)?;
        let frame = Frame::new(def, bits, dec)?;
        let target = frame.project(&p);
        let mut patches = Patch::all(&self.data.carrier, bits, dec)?;
        for _ in 0..64 {
            dec.check_deadline()?;
            let mut next = Vec::new();
            let mut signs = Vec::new();
            for patch in patches {
                let image = patch.projected(&frame)?;
                if !meets(&image[0], &target[0]) || !meets(&image[1], &target[1]) {
                    continue;
                }
                match (&target[2] - &image[2]).sign() {
                    Some(sign) => signs.push(sign),
                    None => next.extend(patch.split()),
                }
                if next.len() > 512 {
                    return Err(unresolved("envelope side search does not localise"));
                }
            }
            if next.is_empty() {
                let Some(first_sign) = signs.first().copied() else {
                    return Err(Refusal::new(
                        "shapes.bspline-domain",
                        "side point leaves the finite envelope chart",
                    ));
                };
                if signs.iter().any(|s| *s != first_sign) {
                    return Err(unresolved("envelope side search returned mixed signs"));
                }
                let toward_ball = first_sign == StrictSign::Positive;
                return Ok(if toward_ball != self.data.reversing {
                    Sign::Positive
                } else {
                    Sign::Negative
                });
            }
            patches = next;
        }
        Err(unresolved(
            "envelope side is within the chart enclosure width",
        ))
    }
}

fn meets(a: &Ball, c: &Ball) -> bool {
    a.hi() >= c.lo() && a.lo() <= c.hi()
}
fn dot(a: &[Ball; 3], c: &[Ball; 3]) -> Ball {
    &(&(&a[0] * &c[0]) + &(&a[1] * &c[1])) + &(&a[2] * &c[2])
}

// Enclosed chart directions: F = (across . S, (toward - axis) . S) is the
// certified injective projection; `lift`, along toward + axis, measures the
// side. In a right-handed frame across x (toward - axis) = toward + axis.
struct Frame {
    first: [Ball; 3],
    second: [Ball; 3],
    lift: [Ball; 3],
}
impl Frame {
    fn new(def: &RollingEnvelopeDefinition, bits: u32, dec: &Decider) -> Result<Self, Refusal> {
        let enclose = |v: &Vector3| -> Result<[Ball; 3], Refusal> {
            Ok([
                dec.enclose_at(&v.x, bits)?,
                dec.enclose_at(&v.y, bits)?,
                dec.enclose_at(&v.z, bits)?,
            ])
        };
        Ok(Self {
            first: enclose(&def.across)?,
            second: enclose(&def.toward.minus(&def.axis))?,
            lift: enclose(&def.toward.plus(&def.axis))?,
        })
    }
    fn project(&self, p: &[Ball; 3]) -> [Ball; 3] {
        [
            dot(p, &self.first),
            dot(p, &self.second),
            dot(p, &self.lift),
        ]
    }
}

// One rational Bezier patch of the carrier, degree (6, 2), as enclosed
// homogeneous controls in u-major order. Subdivision is exact up to rounding,
// and positive weights give the convex-hull enclosure.
struct Patch {
    cell: [[Rational; 2]; 2],
    net: Vec<[Ball; 4]>,
}
impl Patch {
    fn all(carrier: &BSplineSurface, bits: u32, dec: &Decider) -> Result<Vec<Self>, Refusal> {
        let mut patches = Vec::new();
        for (index, span) in carrier.spans(0).enumerate() {
            let mut net = Vec::with_capacity((DEGREE + 1) * 3);
            for i in index * DEGREE..=index * DEGREE + DEGREE {
                for j in 0..3 {
                    let p = &carrier.controls()[i * 3 + j];
                    let w = dec.enclose_at(&carrier.weights()[i * 3 + j], bits)?;
                    let c = math::point(p, bits, dec)?;
                    net.push([&c[0] * &w, &c[1] * &w, &c[2] * &w, w]);
                }
            }
            patches.push(Self {
                cell: [
                    [span[0].clone(), span[1].clone()],
                    [Rational::zero(), Rational::one()],
                ],
                net,
            });
        }
        Ok(patches)
    }
    fn split(self) -> [Self; 2] {
        // A span's u extent is about a tenth of v's; compare in a metric that
        // roughly equalises the two chart speeds.
        let du = &self.cell[0][1] - &self.cell[0][0];
        let dv = (&self.cell[1][1] - &self.cell[1][0]) / q(8);
        let axis = usize::from(dv > du);
        let half = (&self.cell[axis][0] + &self.cell[axis][1]) / q(2);
        let mut low = self.cell.clone();
        let mut high = self.cell;
        low[axis][1] = half.clone();
        high[axis][0] = half;
        let [left, right] = split_net(&self.net, [DEGREE, 2], axis);
        [
            Self {
                cell: low,
                net: left,
            },
            Self {
                cell: high,
                net: right,
            },
        ]
    }
    fn cartesian(&self) -> Result<Vec<[Ball; 3]>, Refusal> {
        self.net
            .iter()
            .map(|h| {
                let w = &h[3];
                if w.lo() <= Rational::zero() {
                    return Err(unresolved("subdivided envelope weight is not positive"));
                }
                Ok([0, 1, 2].map(|i| h[i].checked_div(w).expect("positive weight")))
            })
            .collect()
    }
    fn hull(&self) -> Result<[Ball; 3], Refusal> {
        let points = self.cartesian()?;
        Ok(std::array::from_fn(|i| {
            points[1..]
                .iter()
                .fold(points[0][i].clone(), |hull, p| hull.hull(&p[i]))
        }))
    }
    fn projected(&self, frame: &Frame) -> Result<[Ball; 3], Refusal> {
        let points: Vec<_> = self.cartesian()?.iter().map(|p| frame.project(p)).collect();
        Ok(std::array::from_fn(|i| {
            points[1..]
                .iter()
                .fold(points[0][i].clone(), |hull, p| hull.hull(&p[i]))
        }))
    }
}

// De Casteljau halving of a tensor Bernstein net along one axis.
fn split_net<T: Clone + Halve>(net: &[T], degrees: [usize; 2], axis: usize) -> [Vec<T>; 2] {
    let [m, n] = degrees;
    let index = |i: usize, j: usize| i * (n + 1) + j;
    let mut left = net.to_vec();
    let mut right = net.to_vec();
    let (lines, length) = if axis == 0 {
        (n + 1, m + 1)
    } else {
        (m + 1, n + 1)
    };
    for line in 0..lines {
        let at = |k: usize| {
            if axis == 0 {
                index(k, line)
            } else {
                index(line, k)
            }
        };
        let mut row: Vec<T> = (0..length).map(|k| net[at(k)].clone()).collect();
        let mut low = Vec::with_capacity(length);
        let mut high = vec![row[length - 1].clone(); length];
        for level in 0..length {
            low.push(row[0].clone());
            high[length - 1 - level] = row[length - 1 - level].clone();
            row = (0..row.len() - 1)
                .map(|k| row[k].halve(&row[k + 1]))
                .collect();
            if row.is_empty() {
                break;
            }
        }
        for k in 0..length {
            left[at(k)] = low[k].clone();
            right[at(k)] = high[k].clone();
        }
    }
    [left, right]
}

trait Halve {
    fn halve(&self, other: &Self) -> Self;
}
impl Halve for Ball {
    fn halve(&self, other: &Self) -> Self {
        &(self + other)
            * &b(
                &Rational::from_ratio(1, 2).expect("nonzero"),
                self.precision(),
            )
    }
}
impl Halve for [Ball; 4] {
    fn halve(&self, other: &Self) -> Self {
        std::array::from_fn(|i| self[i].halve(&other[i]))
    }
}

// Producer: Hermite node data from point enclosures of the exact rows. Interior
// values are rational midpoints; transition values are exact so that the end
// sections are the meridian arcs shared with the neighbouring torus pieces.
fn produce(
    def: &RollingEnvelopeDefinition,
    spans: usize,
    dec: &mut Decider,
) -> Result<Vec<Node>, Refusal> {
    let bits = 128.min(dec.policy().max_precision_bits);
    let k = Constants::new(def, bits, dec)?;
    let t = def.transition();
    let mut nodes = Vec::with_capacity(spans + 1);
    for i in 0..=spans {
        dec.check_deadline()?;
        let s = Rational::from_ratio(2 * i as i64 - spans as i64, spans as i64).expect("nonzero");
        let point = Jet::linear(b(&s, bits), k.one.clone());
        let exact = rows(&k, &point, None)?;
        let mid = |ball: &Ball| ball.mid();
        let mut slopes: [Rational; ROWS] = std::array::from_fn(|_| Rational::zero());
        let mut values: [Real; ROWS] = std::array::from_fn(|_| Real::zero());
        let scalar = [&exact.height, &exact.mouth, &exact.tangent, &exact.axial];
        for (j, jet) in scalar.into_iter().enumerate() {
            let [value, slope] = jet.data();
            values[j] = Real::from_rational(mid(&value));
            slopes[j] = mid(&slope);
        }
        let end = i == 0 || i == spans;
        if end {
            values[HEIGHT] = def.radius.clone();
            let z = if i == 0 { -&t.z } else { t.z.clone() };
            values[MOUTH] = -(&z / &(&(&def.mouth + &def.offset) - &t.x));
            // The transition apex is the plane/cylinder corner below the rail.
            values[TANGENT] = Real::zero();
            values[AXIAL] = -&def.radius;
        }
        let fitted = Jet::linear(
            dec.enclose_at(&values[MOUTH], bits)?,
            b(&slopes[MOUTH], bits),
        );
        let with_weight = rows(&k, &point, Some(&fitted))?;
        let [value, slope] = with_weight.weight.expect("weight requested").data();
        values[WEIGHT] = if end {
            // d = 0 at a transition, so the weight is sqrt((1+w^2)/(2(1+u^2))).
            (&(&Real::one() + &values[MOUTH].square())
                / &(&(&Real::one() + &t.t.square()) * &Real::from_i64(2)))
                .sqrt()
        } else {
            Real::from_rational(mid(&value))
        };
        slopes[WEIGHT] = mid(&slope);
        nodes.push(Node { s, values, slopes });
    }
    Ok(nodes)
}

// Recomputed two-sided bound. With rows r_i and weights w_i, both carrier and
// exact section are rational quadratics in v with the same Bernstein basis.
// Rail and mouth rows share their weights exactly, so
//   |S - E| <= max_i |s_i - e_i| + (|eps| / (2 (1 - |eps|))) * r / cos(half)
// where eps is the relative middle-weight error and r / cos(half) bounds every
// exact control's distance from the ball centre. Each row error is the cubic
// Hermite remainder plus the effect of inexact node data.
fn recompute_bound(
    def: &RollingEnvelopeDefinition,
    nodes: &[Node],
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let bits = 128.min(dec.policy().max_precision_bits);
    let k = Constants::new(def, bits, dec)?;
    let mut worst = Rational::zero();
    let enclosed = |node: &Node| -> Result<RowData, Refusal> {
        let mut out: RowData = std::array::from_fn(|_| [k.one.clone(), k.one.clone()]);
        for ((slot, value), slope) in out.iter_mut().zip(&node.values).zip(&node.slopes) {
            *slot = [dec.enclose_at(value, bits)?, b(slope, bits)];
        }
        Ok(out)
    };
    // Exact rows at a node, with the mouth row fixed by its represented data.
    let exact_at = |node: &Node, data: &RowData| -> Result<RowData, Refusal> {
        let point = Jet::linear(b(&node.s, bits), k.one.clone());
        let fitted = Jet::linear(data[MOUTH][0].clone(), data[MOUTH][1].clone());
        let exact = rows(&k, &point, Some(&fitted))?;
        Ok([
            exact.height.data(),
            exact.mouth.data(),
            exact.tangent.data(),
            exact.axial.data(),
            exact.weight.expect("weight requested").data(),
        ])
    };
    let mut previous: Option<(RowData, RowData)> = None;
    for pair in nodes.windows(2) {
        dec.check_deadline()?;
        let left = match previous.take() {
            Some(left) => left,
            None => {
                let data = enclosed(&pair[0])?;
                let exact = exact_at(&pair[0], &data)?;
                (data, exact)
            }
        };
        let data = enclosed(&pair[1])?;
        let exact = exact_at(&pair[1], &data)?;
        let right = (data, exact);
        let span = [&pair[0].s, &pair[1].s];
        let h = span[1] - span[0];
        let s = Jet::linear(b(span[0], bits).hull(&b(span[1], bits)), k.one.clone());
        let fitted = hermite(
            &s,
            span,
            [
                [&left.0[MOUTH][0], &left.0[MOUTH][1]],
                [&right.0[MOUTH][0], &right.0[MOUTH][1]],
            ],
            bits,
        )?;
        let whole = rows(&k, &s, Some(&fitted))?;
        let fourth = [
            &whole.height.0[4],
            &whole.mouth.0[4],
            &whole.tangent.0[4],
            &whole.axial.0[4],
            &whole.weight.as_ref().expect("weight requested").0[4],
        ];
        let h4 = h.pow(4).expect("nonzero span");
        let error = |j: usize| -> Rational {
            let remainder = magnitude(fourth[j]) * h4.clone() / q(16);
            let values = magnitude(&(&left.0[j][0] - &left.1[j][0]))
                .max(magnitude(&(&right.0[j][0] - &right.1[j][0])));
            let slopes = magnitude(&(&left.0[j][1] - &left.1[j][1]))
                + magnitude(&(&right.0[j][1] - &right.1[j][1]));
            remainder + values + q(4) * h.clone() * slopes / q(27)
        };
        let rail = error(HEIGHT);
        let mouth = q(2) * k.mouth.hi() * error(MOUTH);
        // The represented apex is A~ + tangent~ T + axial~ axis with an exact
        // unit tangent T, and A~ moves only by the rail height error.
        let apex = error(TANGENT) + error(AXIAL) + error(HEIGHT);
        let weight = whole.weight.as_ref().expect("weight requested").0[0].lo();
        if weight <= Rational::zero() {
            return Err(unresolved("exact middle weight is not separated from zero"));
        }
        let eps = error(WEIGHT) / weight;
        if eps >= Rational::one() {
            return Err(refusal("middle weight error is not below its exact value"));
        }
        let cosine = whole.cosine.0[0].lo();
        if cosine <= Rational::zero() {
            return Err(unresolved(
                "contact angle is not separated from a half turn",
            ));
        }
        let spread = eps.clone() / (q(2) * (Rational::one() - eps)) * k.radius.hi() / cosine;
        worst = worst.max(rail.max(mouth).max(apex) + spread);
        previous = Some(right);
    }
    Ok(worst)
}

// Local homogeneous Bernstein rows of one span, [rail, middle, mouth], each
// as [x, y, z, w] in the definition frame, elevated to the carrier degree.
fn span_rows(def: &RollingEnvelopeDefinition, pair: &[Node]) -> [[Split; 4]; 3] {
    let t = def.transition();
    let one = Real::one();
    let square = &t.square;
    let span = [&pair[0].s, &pair[1].s];
    let f = |j: usize| {
        cubic(
            span,
            [&pair[0].values[j], &pair[1].values[j]],
            [&pair[0].slopes[j], &pair[1].slopes[j]],
        )
    };
    let [den, cos, sin] = station(&t, span);
    let rail_x = cos.scaled(&def.body);
    let rail_z = sin.scaled(&def.body);
    // Rail row: (body (1-u^2), height (1+u^2), 2 body u) / (1+u^2).
    let rail = [
        rail_x.clone(),
        Split::even(f(HEIGHT)).product(&den, square),
        rail_z.clone(),
        den.clone(),
    ];
    // Middle row: the rail point moved within the cylinder's tangent plane by
    // the fitted tangential and axial offsets, times the fitted middle weight.
    let tangent = Split::even(f(TANGENT));
    let x = rail_x.combine(&tangent.product(&sin, square), &one, &(-&one));
    let z = rail_z.combine(&tangent.product(&cos, square), &one, &one);
    let height = Split::even(combine(&f(HEIGHT), &f(AXIAL), &one, &one));
    let weight = Split::even(f(WEIGHT));
    let middle = [
        weight.product(&x, square),
        weight.product(&height.product(&den, square), square),
        weight.product(&z, square),
        weight.product(&den, square),
    ];
    // Mouth row: the mouth circle at angle pi + 2 atan(w). Its parameter is an
    // exact odd multiple of t at a transition and rational at interior nodes.
    let parts = pair.iter().map(|node| {
        if node.s.abs() == Rational::one() {
            let odd = -&(&t.mouth * &Real::from_rational(node.s.clone()));
            (Real::zero(), odd)
        } else {
            (node.values[MOUTH].clone(), Real::zero())
        }
    });
    let (even, odd): (Vec<Real>, Vec<Real>) = parts.unzip();
    let flat = Rational::zero();
    let w = Split {
        even: cubic(
            span,
            [&even[0], &even[1]],
            [&pair[0].slopes[MOUTH], &pair[1].slopes[MOUTH]],
        ),
        odd: cubic(span, [&odd[0], &odd[1]], [&flat, &flat]),
    };
    let w2 = w.product(&w, square);
    let ones6 = Split::even(vec![one.clone(); 7]);
    let mouth_den = ones6.combine(&w2, &one, &one);
    let mouth = [
        mouth_den.combine(
            &ones6.combine(&w2, &one, &(-&one)),
            &def.offset,
            &(-&def.mouth),
        ),
        Split::even(vec![Real::zero(); 7]),
        w.scaled(&(&def.mouth * &Real::from_i64(-2))),
        mouth_den,
    ];
    [rail, middle, mouth].map(|row| row.map(|c| c.elevated(DEGREE)))
}

// Builds the carrier from node data alone: degree (6, 2), C0 knots at nodes.
fn build(
    def: &RollingEnvelopeDefinition,
    nodes: &[Node],
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let t = def.transition();
    let mut controls = Vec::with_capacity((nodes.len() - 1) * DEGREE * 3 + 3);
    let mut weights = Vec::with_capacity(controls.capacity());
    let mut knots = vec![Rational::from_i64(-1); DEGREE + 1];
    let count = (nodes.len() - 1) * DEGREE + 1;
    for (index, pair) in nodes.windows(2).enumerate() {
        dec.check_deadline()?;
        let rows = span_rows(def, pair);
        let last = index + 2 == nodes.len();
        let take = if last { DEGREE + 1 } else { DEGREE };
        for i in 0..take {
            for row in &rows {
                let (control, weight) = def.control(row, i, &t.t);
                controls.push(control);
                weights.push(weight);
            }
        }
        let knot = pair[1].s.clone();
        if last {
            knots.extend(std::iter::repeat_n(knot, DEGREE + 1));
        } else {
            knots.extend(std::iter::repeat_n(knot, DEGREE));
        }
    }
    let v = vec![
        Rational::zero(),
        Rational::zero(),
        Rational::zero(),
        Rational::one(),
        Rational::one(),
        Rational::one(),
    ];
    BSplineSurface::new([DEGREE, 2], [knots, v], [count, 3], controls, weights, dec)
}

// Proves the end sections are the exact transition meridian arcs: their
// corners are the exact transition feet, the apex is the corner of the
// plane/cylinder section, and the conic weight ratio is that of a quarter circle.
fn check_transitions(
    def: &RollingEnvelopeDefinition,
    carrier: &BSplineSurface,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let same = |a: &Real, c: &Real, dec: &mut Decider| -> Result<(), Refusal> {
        match dec.compare(a, c) {
            Outcome::ProvenCoincident(_) => Ok(()),
            Outcome::Proven(_) => Err(refusal(
                "envelope end section is not the transition meridian",
            )),
            Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
            Outcome::Refused(error) => Err(error),
        }
    };
    let count = carrier.control_count()[0];
    for (sign, column) in [(Sign::Negative, 0), (Sign::Positive, count - 1)] {
        let corners = def.corners(sign);
        let controls = &carrier.controls()[column * 3..column * 3 + 3];
        let weights = &carrier.weights()[column * 3..column * 3 + 3];
        for (control, corner) in [(&controls[0], &corners[0]), (&controls[2], &corners[1])] {
            for (a, c) in [
                (&control.x, &corner.x),
                (&control.y, &corner.y),
                (&control.z, &corner.z),
            ] {
                same(a, c, dec)?;
            }
        }
        let apex = &controls[1];
        let t = def.transition();
        let z = if sign == Sign::Negative {
            -&t.z
        } else {
            t.z.clone()
        };
        let expected = def.world([
            &(&(&def.body / &t.centre) * &t.x),
            &Real::zero(),
            &(&(&def.body / &t.centre) * &z),
        ]);
        for (a, c) in [
            (&apex.x, &expected.x),
            (&apex.y, &expected.y),
            (&apex.z, &expected.z),
        ] {
            same(a, c, dec)?;
        }
        same(
            &(&weights[1].square() * &Real::from_i64(2)),
            &(&weights[0] * &weights[2]),
            dec,
        )?;
    }
    Ok(())
}

// Exact clearances that make the face meet its neighbours only along its
// boundary, proved on Bernstein coefficients of local homogeneous rows:
// - Outside the cylinder: with c, s the exact station cosine and sine,
//   x c + z s - body w is identically zero on the rail and middle rows (the
//   middle row lies in the rail's tangent plane), so it suffices that the mouth
//   row's polynomial is positive; then S . radial > body for every v > 0.
// - Above the mouth plane: the rail height and the middle row's height are
//   nonnegative, the rail strictly, and the mouth row lies in the plane.
// - Inside the hole's sector: every control column but the transition ones lies
//   strictly inside both meridian half-spaces; the transition columns are the
//   meridian arcs proved by `check_transitions`.
fn check_clearance(
    def: &RollingEnvelopeDefinition,
    nodes: &[Node],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let bits = 128.min(dec.policy().max_precision_bits);
    let t = def.transition();
    let e = |x: &Real| dec.enclose_at(x, bits);
    let (x0, z0) = (e(&t.x)?, e(&t.z)?);
    let one = Real::one();
    let spans = nodes.len() - 1;
    for (index, pair) in nodes.windows(2).enumerate() {
        dec.check_deadline()?;
        let rows = span_rows(def, pair);
        let span = [&pair[0].s, &pair[1].s];
        let [den, cos, sin] = station(&t, span);
        let [x, _, z, w] = &rows[2];
        let radial = x
            .product(&cos, &t.square)
            .combine(&z.product(&sin, &t.square), &one, &one);
        let gap = radial.combine(&w.product(&den, &t.square), &one, &(-&def.body));
        if !positive_bernstein(&gap.values(&t.t), bits, dec)? {
            return Err(refusal("envelope mouth row is not clear of the cylinder"));
        }
        let cubic_of = |j: usize| {
            cubic(
                span,
                [&pair[0].values[j], &pair[1].values[j]],
                [&pair[0].slopes[j], &pair[1].slopes[j]],
            )
        };
        let height = cubic_of(HEIGHT);
        let middle = combine(&height, &cubic_of(AXIAL), &one, &one);
        for (row, strict) in [(&height, true), (&middle, false), (&cubic_of(WEIGHT), true)] {
            for c in row {
                let c = e(c)?;
                if c.lo() < Rational::zero() || (strict && c.lo() == Rational::zero()) {
                    return Err(refusal("envelope rows are not clear of the mouth plane"));
                }
            }
        }
        for column in 0..=DEGREE {
            let transition =
                (index == 0 && column == 0) || (index + 1 == spans && column == DEGREE);
            if transition {
                continue;
            }
            for row in &rows {
                let (xh, zh) = (e(&row[0].at(column, &t.t))?, e(&row[2].at(column, &t.t))?);
                let plus = &(&x0 * &zh) - &(&z0 * &xh);
                let minus = -&(&(&x0 * &zh) + &(&z0 * &xh));
                if plus.hi() >= Rational::zero() || minus.hi() >= Rational::zero() {
                    return Err(refusal("envelope leaves the sector between its meridians"));
                }
            }
        }
    }
    Ok(())
}

fn positive_bernstein(c: &[Real], bits: u32, dec: &Decider) -> Result<bool, Refusal> {
    let mut stack = vec![(
        c.iter()
            .map(|x| dec.enclose_at(x, bits))
            .collect::<Result<Vec<_>, _>>()?,
        0,
    )];
    while let Some((poly, depth)) = stack.pop() {
        if poly.iter().all(|x| x.lo() > Rational::zero()) {
            continue;
        }
        if depth >= 8 || poly.iter().any(|x| x.hi() <= Rational::zero()) {
            return Ok(false);
        }
        let [left, right] = split_net(&poly, [poly.len() - 1, 0], 0);
        stack.push((left, depth + 1));
        stack.push((right, depth + 1));
    }
    Ok(true)
}

// Gale-Nikaido on the closed chart rectangle: with F = (across . S,
// (toward - axis) . S), a Jacobian that is a P-matrix everywhere makes F, hence
// the carrier, injective; its positive determinant proves regularity and puts
// the chart normal on the ball side, toward + axis. On each Bezier patch the
// signs are those of polynomials in the homogeneous planar map H = (w, G1, G2):
//   dG1/du: w G1_u - w_u G1,  dG2/dv: w G2_v - w_v G2,  det J: det(H, H_u, H_v),
// proved by positive Bernstein coefficients after bounded subdivision.
fn check_regular(
    def: &RollingEnvelopeDefinition,
    carrier: &BSplineSurface,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let bits = 128.min(dec.policy().max_precision_bits);
    let frame = Frame::new(def, bits, dec)?;
    let cross = |a: &[Ball; 3], c: &[Ball; 3]| -> [Ball; 3] {
        [
            &(&a[1] * &c[2]) - &(&a[2] * &c[1]),
            &(&a[2] * &c[0]) - &(&a[0] * &c[2]),
            &(&a[0] * &c[1]) - &(&a[1] * &c[0]),
        ]
    };
    if dot(&cross(&frame.first, &frame.second), &frame.lift).lo() <= Rational::zero() {
        return Err(refusal("envelope frame is not right-handed"));
    }
    let mut stack = Vec::new();
    for patch in Patch::all(carrier, bits, dec)? {
        let planar: Vec<[Ball; 3]> = patch
            .net
            .iter()
            .map(|h| {
                let x = [h[0].clone(), h[1].clone(), h[2].clone()];
                [h[3].clone(), dot(&x, &frame.first), dot(&x, &frame.second)]
            })
            .collect();
        stack.push((planar, patch.cell, 0));
    }
    while let Some((net, cell, depth)) = stack.pop() {
        dec.check_deadline()?;
        // Translating the planar map by an exact constant leaves its Jacobian
        // unchanged and keeps the products free of large cancelling terms.
        let shift = [1, 2].map(|k| {
            let h = &net[0];
            h[k].checked_div(&h[0]).map(|c| b(&c.mid(), bits))
        });
        let [Some(shift1), Some(shift2)] = shift else {
            return Err(unresolved("subdivided envelope weight is not positive"));
        };
        let poly =
            |f: &dyn Fn(&[Ball; 3]) -> Ball| Poly::new([DEGREE, 2], net.iter().map(f).collect());
        let w = poly(&|h| h[0].clone());
        let g1 = poly(&|h| &h[1] - &(&h[0] * &shift1));
        let g2 = poly(&|h| &h[2] - &(&h[0] * &shift2));
        let (wu, wv) = (w.derivative(0).scaled(), w.derivative(1).scaled());
        let (g1u, g1v) = (g1.derivative(0).scaled(), g1.derivative(1).scaled());
        let (g2u, g2v) = (g2.derivative(0).scaled(), g2.derivative(1).scaled());
        let (w, g1, g2) = (w.scaled(), g1.scaled(), g2.scaled());
        let first = w.mul(&g1u).sub(&wu.mul(&g1));
        let second = w.mul(&g2v).sub(&wv.mul(&g2));
        let det = w
            .mul(&g1u.mul(&g2v).sub(&g2u.mul(&g1v)))
            .sub(&g1.mul(&wu.mul(&g2v).sub(&g2u.mul(&wv))))
            .add(&g2.mul(&wu.mul(&g1v).sub(&g1u.mul(&wv))));
        if first.positive() && second.positive() && det.positive() {
            continue;
        }
        if depth >= 10 {
            return Err(unresolved(
                "envelope chart Jacobian is not certified as a P-matrix",
            ));
        }
        let du = &cell[0][1] - &cell[0][0];
        let dv = (&cell[1][1] - &cell[1][0]) / q(8);
        let axis = usize::from(dv > du);
        let half = (&cell[axis][0] + &cell[axis][1]) / q(2);
        let mut low = cell.clone();
        let mut high = cell;
        low[axis][1] = half.clone();
        high[axis][0] = half;
        let [left, right] = split_net(&net, [DEGREE, 2], axis);
        stack.push((left, low, depth + 1));
        stack.push((right, high, depth + 1));
    }
    Ok(())
}

impl Halve for [Ball; 3] {
    fn halve(&self, other: &Self) -> Self {
        std::array::from_fn(|i| self[i].halve(&other[i]))
    }
}

// Tensor Bernstein polynomial with enclosed coefficients, u-major. After
// `scaled`, coefficient (i, j) carries the factor C(m, i) C(n, j): products are
// then plain convolutions, and positive factors leave every sign unchanged.
struct Poly {
    degrees: [usize; 2],
    c: Vec<Ball>,
}
impl Poly {
    fn new(degrees: [usize; 2], c: Vec<Ball>) -> Self {
        Self { degrees, c }
    }
    fn at(&self, i: usize, j: usize) -> &Ball {
        &self.c[i * (self.degrees[1] + 1) + j]
    }
    // Positive multiples of the partial derivatives; every sign test below
    // multiplies the same derivative factors in each of its terms.
    fn derivative(&self, axis: usize) -> Self {
        let [m, n] = self.degrees;
        let degrees = if axis == 0 { [m - 1, n] } else { [m, n - 1] };
        let mut c = Vec::with_capacity((degrees[0] + 1) * (degrees[1] + 1));
        for i in 0..=degrees[0] {
            for j in 0..=degrees[1] {
                let next = if axis == 0 {
                    self.at(i + 1, j)
                } else {
                    self.at(i, j + 1)
                };
                c.push(next - self.at(i, j));
            }
        }
        Self { degrees, c }
    }
    fn scaled(self) -> Self {
        let [m, n] = self.degrees;
        let bits = self.c[0].precision();
        let rows: Vec<Ball> = (0..=m).map(|i| b(&binomial(m, i), bits)).collect();
        let columns: Vec<Ball> = (0..=n).map(|j| b(&binomial(n, j), bits)).collect();
        let c = self
            .c
            .iter()
            .enumerate()
            .map(|(k, c)| &(c * &rows[k / (n + 1)]) * &columns[k % (n + 1)])
            .collect();
        Self {
            degrees: self.degrees,
            c,
        }
    }
    fn mul(&self, other: &Self) -> Self {
        let [m1, n1] = self.degrees;
        let [m2, n2] = other.degrees;
        let degrees = [m1 + m2, n1 + n2];
        let width = degrees[1] + 1;
        let bits = self.c[0].precision();
        let mut out = vec![b(&Rational::zero(), bits); (degrees[0] + 1) * width];
        for i in 0..=m1 {
            for j in 0..=n1 {
                let a = self.at(i, j);
                for k in 0..=m2 {
                    for l in 0..=n2 {
                        let slot = &mut out[(i + k) * width + j + l];
                        *slot = &*slot + &(a * other.at(k, l));
                    }
                }
            }
        }
        Self { degrees, c: out }
    }
    fn add(&self, other: &Self) -> Self {
        Self {
            degrees: self.degrees,
            c: self.c.iter().zip(&other.c).map(|(a, c)| a + c).collect(),
        }
    }
    fn sub(&self, other: &Self) -> Self {
        Self {
            degrees: self.degrees,
            c: self.c.iter().zip(&other.c).map(|(a, c)| a - c).collect(),
        }
    }
    fn positive(&self) -> bool {
        self.c.iter().all(|c| c.lo() > Rational::zero())
    }
}
