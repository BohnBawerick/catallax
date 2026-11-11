use crate::{
    CurveForm,
    math::*,
    pair::{Frame, Plan},
};
use catallax_numbers::{Ball, Decider, Interval, Rational, Real, Refusal, Sign};
use catallax_shapes::{Curve, CurveId, GeomStore, Point3, Vector3};
use std::{
    collections::{BTreeMap, HashMap},
    ops::{Add, Div, Mul, Neg, Sub},
    sync::{Arc, Mutex, OnceLock},
};

mod cone_torus;
mod cylinder_torus;

const INTERVAL_PRECISION: u32 = 128;

#[derive(Debug, Default)]
struct Enclosures {
    balls: Mutex<BTreeMap<u32, Ball>>,
    interval: OnceLock<Interval>,
}

#[derive(Clone, Debug)]
pub(super) struct E(Arc<Node>);
#[derive(Debug)]
enum Node {
    T,
    Constant(Real, Enclosures),
    Add(E, E),
    Mul(E, E),
    Div(E, E),
    Neg(E),
    Sqrt(E),
}
impl E {
    fn t() -> Self {
        Self(Arc::new(Node::T))
    }
    pub fn constant(x: Real) -> Self {
        Self(Arc::new(Node::Constant(x, Enclosures::default())))
    }
    fn integer(x: i64) -> Self {
        Self::constant(Real::from_i64(x))
    }
    fn sqrt(self) -> Self {
        Self(Arc::new(Node::Sqrt(self)))
    }
    fn square(self) -> Self {
        self.clone() * self
    }
    fn compose(&self, t: &E, cache: &mut HashMap<usize, E>) -> E {
        let key = Arc::as_ptr(&self.0) as usize;
        if let Some(e) = cache.get(&key) {
            return e.clone();
        }
        let e = match self.0.as_ref() {
            Node::T => t.clone(),
            Node::Constant(..) => self.clone(),
            Node::Add(a, b) => a.compose(t, cache) + b.compose(t, cache),
            Node::Mul(a, b) => a.compose(t, cache) * b.compose(t, cache),
            Node::Div(a, b) => a.compose(t, cache) / b.compose(t, cache),
            Node::Neg(a) => -a.compose(t, cache),
            Node::Sqrt(a) => a.compose(t, cache).sqrt(),
        };
        cache.insert(key, e.clone());
        e
    }
}
impl Add for E {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self(Arc::new(Node::Add(self, b)))
    }
}
impl Mul for E {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self(Arc::new(Node::Mul(self, b)))
    }
}
impl Div for E {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        Self(Arc::new(Node::Div(self, b)))
    }
}
impl Neg for E {
    type Output = Self;
    fn neg(self) -> Self {
        Self(Arc::new(Node::Neg(self)))
    }
}
impl Sub for E {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}

// The same formal-series recurrence is used for exact derivatives and interval jets.
trait Scalar: Clone {
    fn integer(x: i64, p: u32) -> Self;
    fn plus(&self, b: &Self) -> Self;
    fn minus(&self, b: &Self) -> Self;
    fn times(&self, b: &Self) -> Self;
    fn divide(&self, b: &Self) -> Option<Self>;
    fn root(&self) -> Option<Self>;
    fn is_zero(&self) -> bool {
        false
    }
    fn constant(
        x: &Real,
        cache: &Enclosures,
        p: u32,
        dec: Option<&Decider>,
    ) -> Result<Self, JetError>;
}
impl Scalar for Real {
    fn integer(x: i64, _: u32) -> Self {
        Real::from_i64(x)
    }
    fn plus(&self, b: &Self) -> Self {
        self + b
    }
    fn minus(&self, b: &Self) -> Self {
        self - b
    }
    fn times(&self, b: &Self) -> Self {
        self * b
    }
    fn divide(&self, b: &Self) -> Option<Self> {
        Some(self / b)
    }
    fn root(&self) -> Option<Self> {
        Some(self.sqrt())
    }
    fn constant(x: &Real, _: &Enclosures, _: u32, _: Option<&Decider>) -> Result<Self, JetError> {
        Ok(x.clone())
    }
}
impl Scalar for Ball {
    fn integer(x: i64, p: u32) -> Self {
        Ball::from_rational(&r(x), p)
    }
    fn plus(&self, b: &Self) -> Self {
        if self.is_exact_zero() {
            b.clone()
        } else if b.is_exact_zero() {
            self.clone()
        } else {
            self + b
        }
    }
    fn minus(&self, b: &Self) -> Self {
        if b.is_exact_zero() {
            self.clone()
        } else {
            self - b
        }
    }
    fn times(&self, b: &Self) -> Self {
        self * b
    }
    fn divide(&self, b: &Self) -> Option<Self> {
        self.checked_div(b)
    }
    fn root(&self) -> Option<Self> {
        self.sqrt()
    }
    fn is_zero(&self) -> bool {
        self.is_exact_zero()
    }
    fn constant(
        x: &Real,
        cache: &Enclosures,
        p: u32,
        dec: Option<&Decider>,
    ) -> Result<Self, JetError> {
        if let Some(r) = x.as_rational() {
            return Ok(Ball::from_rational(&r, p));
        }
        if !catallax_numbers::shared_caches() {
            // A worker leaves the shared constant unfilled for the writer.
            return Ok(dec
                .expect("interval jets have a decider")
                .enclose_at(x, p)?);
        }
        let mut cached = cache.balls.lock().expect("trace constant enclosure lock");
        if let Some(ball) = cached.get(&p) {
            return Ok(ball.clone());
        }
        let ball = dec
            .expect("interval jets have a decider")
            .enclose_at(x, p)?;
        cached.insert(p, ball.clone());
        Ok(ball)
    }
}
impl Scalar for Interval {
    fn integer(x: i64, _: u32) -> Self {
        if (-(1_i64 << 53)..=(1_i64 << 53)).contains(&x) {
            Interval::point(x as f64).unwrap()
        } else {
            Interval::from_rational(&r(x))
        }
    }
    fn plus(&self, b: &Self) -> Self {
        *self + *b
    }
    fn minus(&self, b: &Self) -> Self {
        *self - *b
    }
    fn times(&self, b: &Self) -> Self {
        *self * *b
    }
    fn divide(&self, b: &Self) -> Option<Self> {
        if b.contains_zero() {
            return None;
        }
        let value = *self / *b;
        value.is_finite().then_some(value)
    }
    fn root(&self) -> Option<Self> {
        self.sqrt()
    }
    #[allow(clippy::float_cmp)] // Singleton identity, not approximate coincidence.
    fn is_zero(&self) -> bool {
        self.lo() == 0.0 && self.hi() == 0.0
    }
    fn constant(
        x: &Real,
        cache: &Enclosures,
        p: u32,
        dec: Option<&Decider>,
    ) -> Result<Self, JetError> {
        if let Some(value) = x.rational_ref() {
            return Ok(*cache
                .interval
                .get_or_init(|| Interval::from_rational(value)));
        }
        // A worker leaves the shared constant unread and unfilled for the writer.
        let shared = catallax_numbers::shared_caches();
        if shared
            && p == INTERVAL_PRECISION
            && let Some(value) = cache.interval.get()
        {
            return Ok(*value);
        }
        let value = <Ball as Scalar>::constant(x, cache, p, dec)?.to_interval();
        if shared && p == INTERVAL_PRECISION {
            let _ = cache.interval.set(value);
        }
        Ok(value)
    }
}

#[derive(Debug)]
enum JetError {
    Split,
    Deadline(Refusal),
}
impl From<Refusal> for JetError {
    fn from(refusal: Refusal) -> Self {
        Self::Deadline(refusal)
    }
}
fn check_deadline(dec: Option<&Decider>) -> Result<(), JetError> {
    if let Some(dec) = dec {
        dec.check_deadline()?;
    }
    Ok(())
}
#[derive(Debug)]
struct JetNode {
    expression: E,
    children: [usize; 2],
}

/// The pointer-identical expression DAG, indexed once for repeated jet queries.
#[derive(Debug)]
pub(super) struct JetProgram {
    nodes: Vec<JetNode>,
    roots: [usize; 3],
}

impl JetProgram {
    pub fn new(coords: &[E; 3]) -> Self {
        fn insert(e: &E, nodes: &mut Vec<JetNode>, seen: &mut HashMap<usize, usize>) -> usize {
            let key = Arc::as_ptr(&e.0) as usize;
            if let Some(index) = seen.get(&key) {
                return *index;
            }
            let children = match e.0.as_ref() {
                Node::Add(a, b) | Node::Mul(a, b) | Node::Div(a, b) => {
                    [insert(a, nodes, seen), insert(b, nodes, seen)]
                }
                Node::Neg(a) | Node::Sqrt(a) => [insert(a, nodes, seen), 0],
                Node::T | Node::Constant(..) => [0, 0],
            };
            let index = nodes.len();
            nodes.push(JetNode {
                expression: e.clone(),
                children,
            });
            seen.insert(key, index);
            index
        }
        let mut nodes = Vec::new();
        let mut seen = HashMap::new();
        let roots = coords.each_ref().map(|e| insert(e, &mut nodes, &mut seen));
        Self { nodes, roots }
    }
}

struct JetValues<'a, S> {
    program: &'a JetProgram,
    starts: Vec<Option<usize>>,
    coefficients: Vec<S>,
    products: Vec<Option<S>>,
}

impl<'a, S: Scalar> JetValues<'a, S> {
    fn new(program: &'a JetProgram, n: usize) -> Self {
        Self {
            program,
            starts: vec![None; program.nodes.len()],
            coefficients: Vec::with_capacity(program.nodes.len() * (n + 1)),
            products: Vec::new(),
        }
    }

    fn jets(
        &mut self,
        index: usize,
        t: &S,
        n: usize,
        p: u32,
        dec: Option<&Decider>,
    ) -> Result<usize, JetError> {
        check_deadline(dec)?;
        if let Some(start) = self.starts[index] {
            return Ok(start);
        }
        let start = self.coefficients.len();
        self.coefficients
            .extend(std::iter::repeat_n(S::integer(0, p), n + 1));
        let program = self.program;
        let node = &program.nodes[index];
        let [left, right] = node.children;
        match node.expression.0.as_ref() {
            Node::T => {
                self.coefficients[start] = t.clone();
                if n > 0 {
                    self.coefficients[start + 1] = S::integer(1, p);
                }
            }
            Node::Constant(x, b) => self.coefficients[start] = S::constant(x, b, p, dec)?,
            Node::Add(..) => {
                let a = self.jets(left, t, n, p, dec)?;
                let b = self.jets(right, t, n, p, dec)?;
                for k in 0..=n {
                    check_deadline(dec)?;
                    self.coefficients[start + k] =
                        self.coefficients[a + k].plus(&self.coefficients[b + k]);
                }
            }
            Node::Neg(_) => {
                let a = self.jets(left, t, n, p, dec)?;
                for k in 0..=n {
                    check_deadline(dec)?;
                    self.coefficients[start + k] =
                        S::integer(0, p).minus(&self.coefficients[a + k]);
                }
            }
            Node::Mul(a, b) => {
                let square = Arc::ptr_eq(&a.0, &b.0);
                let a = self.jets(left, t, n, p, dec)?;
                let b = self.jets(right, t, n, p, dec)?;
                let a_order = self.coefficients[a..=a + n]
                    .iter()
                    .rposition(|v| !v.is_zero())
                    .unwrap_or(0);
                let b_order = self.coefficients[b..=b + n]
                    .iter()
                    .rposition(|v| !v.is_zero())
                    .unwrap_or(0);
                self.products.clear();
                if square {
                    self.products.resize_with(a_order + 1, || None);
                }
                for k in 0..=a_order.saturating_add(b_order).min(n) {
                    for i in k.saturating_sub(b_order)..=k.min(a_order) {
                        check_deadline(dec)?;
                        if square {
                            let index = i.min(k - i);
                            if i == index {
                                self.products[index] = (!self.coefficients[a + i].is_zero()
                                    && !self.coefficients[a + k - i].is_zero())
                                .then(|| {
                                    self.coefficients[a + i].times(&self.coefficients[a + k - i])
                                });
                            }
                            if let Some(product) = &self.products[i.min(k - i)] {
                                self.coefficients[start + k] =
                                    self.coefficients[start + k].plus(product);
                            }
                        } else if !self.coefficients[a + i].is_zero()
                            && !self.coefficients[b + k - i].is_zero()
                        {
                            let product =
                                self.coefficients[a + i].times(&self.coefficients[b + k - i]);
                            self.coefficients[start + k] =
                                self.coefficients[start + k].plus(&product);
                        }
                    }
                }
            }
            Node::Div(..) => {
                let a = self.jets(left, t, n, p, dec)?;
                let b = self.jets(right, t, n, p, dec)?;
                let inverse = S::integer(1, p)
                    .divide(&self.coefficients[b])
                    .ok_or(JetError::Split)?;
                let order = self.coefficients[b..=b + n]
                    .iter()
                    .rposition(|v| !v.is_zero())
                    .unwrap_or(0);
                for k in 0..=n {
                    let mut v = self.coefficients[a + k].clone();
                    for i in 1..=k.min(order) {
                        check_deadline(dec)?;
                        if !self.coefficients[b + i].is_zero()
                            && !self.coefficients[start + k - i].is_zero()
                        {
                            v = v.minus(
                                &self.coefficients[b + i].times(&self.coefficients[start + k - i]),
                            );
                        }
                    }
                    self.coefficients[start + k] = if v.is_zero() { v } else { v.times(&inverse) };
                }
            }
            Node::Sqrt(_) => {
                let a = self.jets(left, t, n, p, dec)?;
                self.coefficients[start] = self.coefficients[a].root().ok_or(JetError::Split)?;
                let inverse = if n == 0 {
                    None
                } else {
                    Some(
                        S::integer(1, p)
                            .divide(&S::integer(2, p).times(&self.coefficients[start]))
                            .ok_or(JetError::Split)?,
                    )
                };
                for k in 1..=n {
                    let mut v = self.coefficients[a + k].clone();
                    self.products.clear();
                    for i in 1..k {
                        check_deadline(dec)?;
                        if i <= k / 2 {
                            self.products.push(
                                (!self.coefficients[start + i].is_zero()
                                    && !self.coefficients[start + k - i].is_zero())
                                .then(|| {
                                    self.coefficients[start + i]
                                        .times(&self.coefficients[start + k - i])
                                }),
                            );
                        }
                        // Retain each symmetric subtraction's original order and precision.
                        if let Some(product) = &self.products[i.min(k - i) - 1] {
                            v = v.minus(product);
                        }
                    }
                    self.coefficients[start + k] = if v.is_zero() {
                        v
                    } else {
                        v.times(inverse.as_ref().unwrap())
                    };
                }
            }
        }
        self.starts[index] = Some(start);
        Ok(start)
    }
}
#[derive(Clone, Debug)]
pub(super) struct Linear {
    pub origin: Point3,
    pub vector: Vector3,
}
impl Linear {
    fn at(&self, p: &Point3) -> Real {
        p.minus(&self.origin).dot(&self.vector)
    }
}
#[derive(Clone, Debug)]
enum Inverse {
    Circle {
        cos: Linear,
        sin: Linear,
    },
    Squared {
        value: Linear,
        lo: Real,
        span: Real,
    },
    Affine {
        value: Linear,
        lo: Real,
        span: Real,
    },
    TorusArc {
        coordinates: Arc<crate::pair::cone_torus::Coordinates>,
        quarter: usize,
        lo: Real,
        span: Real,
        orders: [usize; 2],
    },
}
#[derive(Clone, Debug)]
pub(super) struct Chart {
    pub coords: [E; 3],
    inverse: Inverse,
    checks: Vec<(Linear, i64)>,
    reverse: bool,
}
impl Chart {
    pub fn position(&self, t: &Real) -> Point3 {
        let j = self.exact_jets(t, 0);
        Point3::new(j[0][0].clone(), j[0][1].clone(), j[0][2].clone())
    }
    pub fn exact_jets(&self, t: &Real, n: usize) -> Vec<[Real; 3]> {
        let program = JetProgram::new(&self.coords);
        let mut values = JetValues::new(&program, n);
        let j = program
            .roots
            .map(|root| values.jets(root, t, n, 0, None).unwrap());
        (0..=n)
            .map(|k| std::array::from_fn(|axis| values.coefficients[j[axis] + k].clone()))
            .collect()
    }
    pub fn endpoint_tangent(&self, end: bool) -> Vector3 {
        let order = match &self.inverse {
            Inverse::TorusArc { orders, .. } if orders[usize::from(end ^ self.reverse)] == 3 => {
                // A triple radicand root at the cone apex uses a quadratic
                // parameter. Its first nonzero position jet has order two.
                2
            }
            _ => 1,
        };
        let jets = self.exact_jets(&Real::from_i64(i64::from(end)), order);
        let tangent = Vector3::new(
            jets[order][0].clone(),
            jets[order][1].clone(),
            jets[order][2].clone(),
        );
        if end && order % 2 == 0 {
            tangent.scaled(&Real::from_i64(-1))
        } else {
            tangent
        }
    }
    pub fn locate(&self, p: &Point3, dec: &mut Decider) -> Result<Option<Real>, Failure> {
        for (check, s) in &self.checks {
            let v = check.at(p) * Real::from_i64(*s);
            if sign(&v, dec)? == Sign::Negative {
                return Ok(None);
            }
        }
        let t = match &self.inverse {
            Inverse::Circle { cos, sin } => sin.at(p) / (Real::one() + cos.at(p)),
            Inverse::Squared { value, lo, span } => {
                let u = (&value.at(p) - lo) / span.clone();
                u.sqrt() / (Real::one() + (Real::one() - u).sqrt())
            }
            Inverse::Affine { value, lo, span } => (&value.at(p) - lo) / span.clone(),
            Inverse::TorusArc {
                coordinates,
                quarter,
                lo,
                span,
                orders,
            } => {
                let angle = coordinates.minor_angle(p);
                let Some(value) =
                    crate::pair::cone_torus::quarter_parameter(&angle, *quarter, dec)?
                else {
                    return Ok(None);
                };
                let u = (&value - lo) / span.clone();
                if sign(&u, dec)? == Sign::Negative
                    || sign(&(&u - &Real::one()), dec)? == Sign::Positive
                {
                    return Ok(None);
                }
                match orders.map(|order| order % 2 == 1) {
                    [false, false] => u,
                    [true, false] => u.sqrt(),
                    [false, true] => Real::one() - (Real::one() - u).sqrt(),
                    [true, true] => u.sqrt() / (Real::one() + (Real::one() - u).sqrt()),
                }
            }
        };
        Ok(Some(if self.reverse { Real::one() - t } else { t }))
    }
    fn reverse(&mut self) {
        let t = E::integer(1) - E::t();
        let mut cache = HashMap::new();
        self.coords = self.coords.each_ref().map(|e| e.compose(&t, &mut cache));
        self.reverse = !self.reverse;
    }
    pub fn restricted(&self, from: Real, to: Real, dec: &mut Decider) -> Result<[E; 3], Refusal> {
        let t = E::constant(from.clone()) + E::constant(to - from) * E::t();
        let mut cache = HashMap::new();
        let coords = self.coords.each_ref().map(|e| e.compose(&t, &mut cache));
        dec.check_deadline()?;
        Ok(coords)
    }
}

/// The exact rational polynomial on one already-selected finite knot span.
/// Reparameterising this expression never extrapolates it to another span.
pub(super) fn spline_span(
    curve: &catallax_shapes::BSplineCurve,
    span: usize,
    from: Real,
    to: Real,
    dec: &Decider,
) -> Result<[E; 3], Refusal> {
    // Re-expanded about the span's start in its own parameter `t`, so the
    // coefficients stay near the controls' size. In the global parameter they
    // cancel heavily on a span far from zero, and Horner's rule over an
    // interval of `t` then encloses far more than the curve.
    let width = &to - &from;
    let polynomial = |p: &catallax_numbers::RealPoly| {
        let terms = p.terms();
        let Some(degree) = terms
            .iter()
            .map(|(monomial, _)| monomial.exponents()[0] as usize)
            .max()
        else {
            return E::integer(0);
        };
        let mut global = vec![Real::zero(); degree + 1];
        for (monomial, coefficient) in &terms {
            global[monomial.exponents()[0] as usize] = coefficient.clone();
        }
        // Taylor shift by `from`: repeated synthetic division.
        for k in 0..degree {
            for j in (k..degree).rev() {
                global[j] = &global[j] + &(&global[j + 1] * &from);
            }
        }
        let mut scale = Real::one();
        let local: Vec<Real> = global
            .iter()
            .map(|value| {
                let scaled = value * &scale;
                scale = &scale * &width;
                scaled
            })
            .collect();
        let mut expression = E::constant(local[degree].clone());
        for value in local[..degree].iter().rev() {
            expression = expression * E::t() + E::constant(value.clone());
        }
        expression
    };
    let [x, y, z, w] = curve
        .query_span_polynomials(span, dec)?
        .each_ref()
        .map(polynomial);
    if matches!(w.0.as_ref(), Node::Constant(value, _)
        if value.rational_ref().is_some_and(|q| q == &Rational::one()))
    {
        return Ok([x, y, z]);
    }
    Ok([x / w.clone(), y / w.clone(), z / w])
}
fn enclose_jets(
    program: &JetProgram,
    t: &Ball,
    n: usize,
    dec: &Decider,
) -> Result<Vec<[Ball; 3]>, JetError> {
    let mut values = JetValues::new(program, n);
    let [a, b, c] = program.roots;
    let a = values.jets(a, t, n, t.precision(), Some(dec))?;
    let b = values.jets(b, t, n, t.precision(), Some(dec))?;
    let c = values.jets(c, t, n, t.precision(), Some(dec))?;
    Ok((0..=n)
        .map(|k| [a, b, c].map(|start| values.coefficients[start + k].clone()))
        .collect())
}
pub(super) fn interval_jets(
    program: &JetProgram,
    t: Interval,
    n: usize,
    dec: &Decider,
) -> Result<Option<Vec<[Interval; 3]>>, Refusal> {
    let mut cache = JetValues::new(program, n);
    let result = program
        .roots
        .map(|root| cache.jets(root, &t, n, INTERVAL_PRECISION, Some(dec)));
    let mut values = Vec::with_capacity(3);
    for value in result {
        match value {
            Ok(value) => values.push(value),
            Err(JetError::Split) => return Ok(None),
            Err(JetError::Deadline(refusal)) if refusal.code == "numbers.deadline" => {
                return Err(refusal);
            }
            Err(JetError::Deadline(_)) => return Ok(None),
        }
    }
    let result: Vec<_> = (0..=n)
        .map(|k| std::array::from_fn(|axis| cache.coefficients[values[axis] + k]))
        .collect();
    Ok(result
        .iter()
        .flatten()
        .all(Interval::is_finite)
        .then_some(result))
}

// Point evaluation needs no Taylor coefficient vectors. Retain the order-zero
// recurrence's exact expression definitions, including its zero/reciprocal nodes.
fn exact_value(
    e: &E,
    parameter: &Real,
    zero: &Real,
    one: &Real,
    dec: &Decider,
    cache: &mut HashMap<usize, Real>,
) -> Result<Real, Refusal> {
    dec.check_deadline()?;
    let key = Arc::as_ptr(&e.0) as usize;
    if let Some(value) = cache.get(&key) {
        return Ok(value.clone());
    }
    let mut child = |e| exact_value(e, parameter, zero, one, dec, cache);
    let value = match e.0.as_ref() {
        Node::T => parameter.clone(),
        Node::Constant(value, _) => value.clone(),
        Node::Add(a, b) => child(a)? + child(b)?,
        Node::Neg(a) => zero - &child(a)?,
        Node::Mul(a, b) => zero + &(child(a)? * child(b)?),
        Node::Div(a, b) => child(a)? * (one / &child(b)?),
        Node::Sqrt(a) => child(a)?.sqrt(),
    };
    cache.insert(key, value.clone());
    Ok(value)
}

pub(super) fn exact_position(
    coords: &[E; 3],
    t: &Rational,
    dec: &Decider,
) -> Result<Point3, Refusal> {
    let parameter = Real::from(t.clone());
    let zero = Real::zero();
    let one = Real::one();
    let mut cache = HashMap::new();
    let mut coordinate = |e| exact_value(e, &parameter, &zero, &one, dec, &mut cache);
    Ok(Point3::new(
        coordinate(&coords[0])?,
        coordinate(&coords[1])?,
        coordinate(&coords[2])?,
    ))
}

pub(super) fn ball_jets(
    program: &JetProgram,
    t: &Ball,
    n: usize,
    dec: &Decider,
) -> Result<Vec<[Ball; 3]>, Refusal> {
    dec.check_deadline()?;
    let p = t.precision();
    let wide = t.with_precision(p.saturating_add(64));
    match enclose_jets(program, &wide, n, dec) {
        Ok(j) => {
            return Ok(j
                .into_iter()
                .map(|v| v.map(|b| b.with_precision(p)))
                .collect());
        }
        Err(JetError::Deadline(refusal)) => return Err(refusal),
        Err(JetError::Split) => {}
    }
    let mid = t.mid();
    let left = Ball::from_rational(&t.lo(), p.saturating_add(64))
        .hull(&Ball::from_rational(&mid, p.saturating_add(64)));
    let right = Ball::from_rational(&mid, p.saturating_add(64))
        .hull(&Ball::from_rational(&t.hi(), p.saturating_add(64)));
    let a = ball_jets(program, &left, n, dec)?;
    let b = ball_jets(program, &right, n, dec)?;
    Ok(a.into_iter()
        .zip(b)
        .map(|(a, b)| std::array::from_fn(|i| a[i].hull(&b[i]).with_precision(p)))
        .collect())
}
fn trig() -> (E, E) {
    let t = E::t();
    let denom = E::integer(1) + t.clone().square();
    (
        (E::integer(1) - t.clone().square()) / denom.clone(),
        E::integer(2) * t / denom,
    )
}
fn point(origin: &Point3, e1: &Vector3, e2: &Vector3, e3: &Vector3, x: E, y: E, z: E) -> [E; 3] {
    let o = [&origin.x, &origin.y, &origin.z];
    let a = [&e1.x, &e1.y, &e1.z];
    let b = [&e2.x, &e2.y, &e2.z];
    let c = [&e3.x, &e3.y, &e3.z];
    std::array::from_fn(|i| {
        E::constant(o[i].clone())
            + E::constant(a[i].clone()) * x.clone()
            + E::constant(b[i].clone()) * y.clone()
            + E::constant(c[i].clone()) * z.clone()
    })
}
fn lin(origin: &Point3, v: &Vector3) -> Linear {
    Linear {
        origin: origin.clone(),
        vector: v.clone(),
    }
}
fn analytic(center: Point3, u: Vector3, v: Vector3) -> Vec<Chart> {
    let zero = Vector3::new(Real::zero(), Real::zero(), Real::zero());
    (0..4)
        .map(|q| {
            let (a, b) = match q {
                0 => (u.clone(), v.clone()),
                1 => (v.clone(), u.scaled(&Real::from_i64(-1))),
                2 => (u.scaled(&Real::from_i64(-1)), v.scaled(&Real::from_i64(-1))),
                _ => (v.scaled(&Real::from_i64(-1)), u.clone()),
            };
            let cos = lin(&center, &a.scaled(&(Real::one() / a.norm_squared())));
            let sin = lin(&center, &b.scaled(&(Real::one() / b.norm_squared())));
            let (c, s) = trig();
            Chart {
                coords: point(&center, &a, &b, &zero, c, s, E::integer(0)),
                inverse: Inverse::Circle {
                    cos: cos.clone(),
                    sin: sin.clone(),
                },
                checks: vec![(cos, 1), (sin, 1)],
                reverse: false,
            }
        })
        .collect()
}
pub(super) fn form_charts(form: &CurveForm) -> Vec<Chart> {
    match form {
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        } => {
            let u = ref_dir
                .clone()
                .scaled(&root_real(&(radius_sq / &real(&ref_dir.norm_squared()))));
            let v = normal.cross(ref_dir).scaled(&root_real(
                &(radius_sq / &real(&normal.cross(ref_dir).norm_squared())),
            ));
            analytic(center.clone(), u, v)
        }
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        } => {
            let u = major_dir
                .clone()
                .scaled(&root_real(&(major_sq / &real(&major_dir.norm_squared()))));
            let cross = normal.cross(major_dir);
            let v = cross
                .clone()
                .scaled(&root_real(&(minor_sq / &real(&cross.norm_squared()))));
            analytic(center.clone(), u, v)
        }
        _ => unreachable!("closed analytic form"),
    }
}
pub(super) fn offset_profile_quarter(
    definition: &catallax_shapes::EllipticOffsetDefinition,
    axial: &Real,
    quarter: usize,
    dec: &mut Decider,
) -> Result<[E; 3], Refusal> {
    use catallax_shapes::{OffsetSide, Surface};
    if quarter >= 4 {
        return Err(Refusal::new(
            "intersect.offset-quarter",
            "profile quarter must be in 0..4",
        ));
    }
    dec.check_deadline()?;
    let Surface::EllipticCylinder(original) = definition.original() else {
        unreachable!("validated elliptic-offset definition");
    };
    if definition.side() == OffsetSide::Inward {
        let reach = &original.minor_sq / &original.major_sq.sqrt();
        if sign(&(&reach - definition.distance()), dec).map_err(Refusal::from)? != Sign::Positive {
            return Err(Refusal::new(
                "intersect.offset-nonregular-branch",
                "inward profile must remain strictly below the original normal reach",
            ));
        }
    }
    if axial.rational_ref().is_none() {
        dec.enclose_budgeted_at(axial, 64u32.min(dec.policy().max_precision_bits))?;
    }
    let unit = |v: &Vector3| v.scaled(&(&Real::one() / &v.norm_squared().sqrt()));
    let major = unit(&original.major_dir);
    let axis = unit(&original.axis);
    let minor = axis.cross(&major);
    let center = original.origin.plus(&axis.scaled(axial));
    let a = E::constant(original.major_sq.sqrt());
    let b = E::constant(original.minor_sq.sqrt());
    let t = E::t();
    let square = t.clone().square();
    let denominator = E::integer(1) + square.clone();
    let (c, s) = rotate_quarter(E::integer(1) - square, E::integer(2) * t, quarter);
    let mut x = a.clone() * c.clone() / denominator.clone();
    let mut y = b.clone() * s.clone() / denominator;
    if !definition
        .distance()
        .rational_ref()
        .is_some_and(Rational::is_zero)
    {
        let d = E::constant(match definition.side() {
            OffsetSide::Outward => definition.distance().clone(),
            OffsetSide::Inward => -definition.distance(),
        });
        let bc = b * c;
        let as_ = a * s;
        let rho = (bc.clone().square() + as_.clone().square()).sqrt();
        x = x + d.clone() * bc / rho.clone();
        y = y + d * as_ / rho;
    }
    let origin = [&center.x, &center.y, &center.z];
    let major = [&major.x, &major.y, &major.z];
    let minor = [&minor.x, &minor.y, &minor.z];
    Ok(std::array::from_fn(|i| {
        let mut value = if origin[i].rational_ref().is_some_and(Rational::is_zero) {
            None
        } else {
            Some(E::constant(origin[i].clone()))
        };
        for (coefficient, expression) in [(major[i], &x), (minor[i], &y)] {
            if coefficient.rational_ref().is_some_and(Rational::is_zero) {
                continue;
            }
            let term = if coefficient
                .rational_ref()
                .is_some_and(|r| r == &Rational::one())
            {
                expression.clone()
            } else {
                E::constant(coefficient.clone()) * expression.clone()
            };
            value = Some(match value {
                Some(value) => value + term,
                None => term,
            });
        }
        value.unwrap_or_else(|| E::integer(0))
    }))
}

fn rotate_quarter(c: E, s: E, q: usize) -> (E, E) {
    match q {
        0 => (c, s),
        1 => (-s, c),
        2 => (-c, -s),
        _ => (s, -c),
    }
}

fn rotated(c: E, s: E, q: usize) -> (E, E) {
    let (c, s) = rotate_quarter(c, s, q);
    // Rational parameter-plane rotation keeps the table's witnesses away
    // from chart seams without introducing angles into any predicate.
    let a = E::constant(real(&Rational::from_ratio(3, 5).unwrap()));
    let b = E::constant(real(&Rational::from_ratio(4, 5).unwrap()));
    (a.clone() * c.clone() - b.clone() * s.clone(), b * c + a * s)
}
pub(super) fn branch_charts(
    geom: &GeomStore,
    c: CurveId,
    plan: &Plan,
    dec: &mut Decider,
) -> Result<Vec<Chart>, Failure> {
    let Curve::Intersection(curve) = geom.curve(c) else {
        unreachable!()
    };
    let witness = geom.position(curve.witness);
    let branch = if let Some(d) = &plan.divider {
        if sign(&d.value(&witness), dec)? == Sign::Positive {
            1
        } else {
            -1
        }
    } else {
        1
    };
    let mut charts = vec![];
    let forward;
    match plan.frame.as_ref().expect("quartic frame") {
        Frame::CylinderSphere {
            origin,
            e1,
            e2,
            e3,
            radius,
            d,
            k,
            mode,
            sense,
        } => {
            forward = match *mode {
                6 => *sense,
                7 => *sense * branch,
                8 => -*sense * branch,
                _ => unreachable!("CS branch configuration"),
            };
            let origin = origin.clone();
            let xcoord = lin(&origin, e1);
            let ycoord = lin(&origin, e2);
            let zcoord = lin(&origin, e3);
            if *mode == 8 {
                for q in 0..4 {
                    let (c, s) = trig();
                    let (c, s) = rotated(c, s, q);
                    let x = E::constant(real(radius)) * c;
                    let y = E::constant(real(radius)) * s;
                    let z = E::integer(branch)
                        * (E::constant(k.clone())
                            + E::constant(&Real::from_i64(2) * d) * x.clone())
                        .sqrt();
                    let (cos, sin, checks) = circle_inverse(&origin, e1, e2, &real(radius), q);
                    charts.push(Chart {
                        coords: point(&origin, e1, e2, e3, x, y, z),
                        inverse: Inverse::Circle { cos, sin },
                        checks,
                        reverse: false,
                    });
                }
            } else {
                let lo = if *mode == 7 {
                    -real(radius)
                } else {
                    -k / (&Real::from_i64(2) * d)
                };
                let span = &real(radius) - &lo;
                let signs = if *mode == 7 {
                    vec![(1, branch, false), (-1, branch, true)]
                } else {
                    vec![(1, 1, false), (-1, 1, true), (-1, -1, false), (1, -1, true)]
                };
                for (ys, zs, reverse) in signs {
                    let (c, s) = trig();
                    let x =
                        E::constant(lo.clone()) + E::constant(span.clone()) * s.clone().square();
                    let y = if *mode == 7 {
                        E::integer(ys) * E::constant(span.clone()) * s.clone() * c
                    } else {
                        E::integer(ys)
                            * c
                            * (E::constant(span.clone()) * (E::constant(real(radius)) + x.clone()))
                                .sqrt()
                    };
                    let z = E::integer(zs)
                        * s
                        * E::constant((&Real::from_i64(2) * d) * span.clone()).sqrt();
                    let mut chart = Chart {
                        coords: point(&origin, e1, e2, e3, x, y, z),
                        inverse: Inverse::Squared {
                            value: xcoord.clone(),
                            lo: lo.clone(),
                            span: span.clone(),
                        },
                        checks: vec![(ycoord.clone(), ys), (zcoord.clone(), zs)],
                        reverse: false,
                    };
                    if reverse {
                        chart.reverse();
                    }
                    charts.push(chart);
                }
            }
        }
        Frame::Cylinders {
            origin,
            e1,
            e2,
            e3,
            big,
            small,
            delta,
            alpha,
            gamma,
            mode,
            sense,
        } => {
            forward = *sense * if *mode == 3 { 1 } else { branch };
            let origin = origin.clone();
            let xcoord = lin(&origin, e1);
            let ycoord = lin(&origin, e2);
            let residual = Linear {
                origin: origin.clone(),
                vector: e3.scaled(alpha).minus(&e1.scaled(gamma)),
            };
            if *mode == 5 {
                for q in 0..4 {
                    let (c, s) = trig();
                    let (c, s) = rotated(c, s, q);
                    let y = E::constant(delta.clone()) + E::constant(real(small)) * s;
                    let x = E::integer(branch)
                        * (E::constant(big.square()) - y.clone().square()).sqrt();
                    let z = (E::constant(gamma.clone()) * x.clone() + E::constant(real(small)) * c)
                        / E::constant(alpha.clone());
                    let shifted = origin.plus(&e2.scaled(delta));
                    let (cos, sin, checks) =
                        circle_inverse(&shifted, &residual.vector, e2, &real(small), q);
                    charts.push(Chart {
                        coords: point(&origin, e1, e2, e3, x, y, z),
                        inverse: Inverse::Circle { cos, sin },
                        checks,
                        reverse: false,
                    });
                }
            } else {
                let lo = delta - &real(small);
                let span = &real(big) - &lo;
                let signs = if *mode == 4 {
                    vec![(branch, 1, false), (branch, -1, true)]
                } else {
                    vec![(1, 1, false), (-1, 1, true), (-1, -1, false), (1, -1, true)]
                };
                for (xs, zs, reverse) in signs {
                    let (c, s) = trig();
                    let y =
                        E::constant(lo.clone()) + E::constant(span.clone()) * s.clone().square();
                    let x = E::integer(xs)
                        * c.clone()
                        * (E::constant(span.clone()) * (E::constant(real(big)) + y.clone())).sqrt();
                    let res = if *mode == 4 {
                        E::integer(zs) * E::constant(span.clone()) * s * c
                    } else {
                        E::integer(zs)
                            * s
                            * (E::constant(span.clone())
                                * (E::constant(delta + &real(small)) - y.clone()))
                            .sqrt()
                    };
                    let z =
                        (E::constant(gamma.clone()) * x.clone() + res) / E::constant(alpha.clone());
                    let mut chart = Chart {
                        coords: point(&origin, e1, e2, e3, x, y, z),
                        inverse: Inverse::Squared {
                            value: ycoord.clone(),
                            lo: lo.clone(),
                            span: span.clone(),
                        },
                        checks: vec![(xcoord.clone(), xs), (residual.clone(), zs)],
                        reverse: false,
                    };
                    if reverse {
                        chart.reverse();
                    }
                    charts.push(chart);
                }
            }
        }
        frame @ Frame::CylinderTorus { .. } => {
            let (pieces, direction) = cylinder_torus::charts(frame, branch);
            charts = pieces;
            forward = direction;
        }
        Frame::ConeTorus(frame) => {
            let component = frame.component(&witness, dec)?;
            return Ok(cone_torus::charts(frame, component));
        }
    }
    // In the cylinder/sphere frame grad(cylinder) x grad(sphere)
    // is -4*z times increasing cylinder angle. For two cylinders its
    // y component is -4*alpha*x*(alpha*z-gamma*x). Positive implicit
    // scale factors and our right-handed frames preserve these signs.
    if forward < 0 {
        charts.reverse();
        for ch in &mut charts {
            ch.reverse();
        }
    }
    Ok(charts)
}
fn circle_inverse(
    origin: &Point3,
    e1: &Vector3,
    e2: &Vector3,
    radius: &Real,
    q: usize,
) -> (Linear, Linear, Vec<(Linear, i64)>) {
    let cosine = real(&Rational::from_ratio(3, 5).unwrap());
    let sine = real(&Rational::from_ratio(4, 5).unwrap());
    let u = e1.scaled(&cosine).plus(&e2.scaled(&sine));
    let v = e2.scaled(&cosine).minus(&e1.scaled(&sine));
    let (e1, e2) = (&u, &v);
    let (a, b) = match q {
        0 => (e1.clone(), e2.clone()),
        1 => (e2.clone(), e1.scaled(&Real::from_i64(-1))),
        2 => (
            e1.scaled(&Real::from_i64(-1)),
            e2.scaled(&Real::from_i64(-1)),
        ),
        _ => (e2.scaled(&Real::from_i64(-1)), e1.clone()),
    };
    let cos = lin(origin, &a.scaled(&(&Real::one() / radius)));
    let sin = lin(origin, &b.scaled(&(&Real::one() / radius)));
    let checks = vec![(cos.clone(), 1), (sin.clone(), 1)];
    (cos, sin, checks)
}
pub(super) fn locate(
    charts: &[Chart],
    p: &Point3,
    dec: &mut Decider,
) -> Result<(usize, Real), Failure> {
    for (i, ch) in charts.iter().enumerate() {
        if let Some(t) = ch.locate(p, dec)? {
            return Ok((i, t));
        }
    }
    panic!("point must lie on curve")
}
pub(super) fn segment(from: &Point3, to: &Point3, dec: &mut Decider) -> Result<[E; 3], Refusal> {
    let a = [&from.x, &from.y, &from.z];
    let b = [&to.x, &to.y, &to.z];
    let coords: [E; 3] =
        std::array::from_fn(|i| E::constant(a[i].clone()) + E::constant(b[i] - a[i]) * E::t());
    dec.check_deadline()?;
    Ok(coords)
}
pub(super) fn conic_arc(
    form: &CurveForm,
    from: &Point3,
    to: &Point3,
    dec: &mut Decider,
) -> Result<Chart, Refusal> {
    let (origin, normal, axis) = match form {
        CurveForm::Parabola {
            vertex,
            normal,
            axis_dir,
            ..
        } => (vertex, normal, axis_dir),
        CurveForm::Hyperbola {
            center,
            normal,
            axis_dir,
            ..
        } => (center, normal, axis_dir),
        _ => unreachable!("open conic"),
    };
    let u = unit(axis);
    let v = unit(&normal.cross(axis));
    let lo = from.minus(origin).dot(&v);
    let span = to.minus(from).dot(&v);
    let y = E::constant(lo.clone()) + E::constant(span.clone()) * E::t();
    let x = match form {
        CurveForm::Parabola { focal, .. } => {
            y.clone().square() / E::constant(&Real::from_i64(4) * focal)
        }
        CurveForm::Hyperbola {
            transverse_sq,
            conjugate_sq,
            branch,
            ..
        } => {
            let x = (E::constant(transverse_sq.clone())
                * (E::integer(1) + y.clone().square() / E::constant(conjugate_sq.clone())))
            .sqrt();
            if *branch == catallax_numbers::StrictSign::Positive {
                x
            } else {
                -x
            }
        }
        _ => unreachable!("open conic"),
    };
    let o = [&origin.x, &origin.y, &origin.z];
    let a = [&u.x, &u.y, &u.z];
    let b = [&v.x, &v.y, &v.z];
    let coords = std::array::from_fn(|i| {
        E::constant(o[i].clone())
            + E::constant(a[i].clone()) * x.clone()
            + E::constant(b[i].clone()) * y.clone()
    });
    dec.check_deadline()?;
    Ok(Chart {
        coords,
        inverse: Inverse::Affine {
            value: lin(origin, &v),
            lo,
            span,
        },
        checks: vec![],
        reverse: false,
    })
}
pub(super) fn conic_local(
    form: &CurveForm,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Chart, Refusal> {
    let (normal, axis) = match form {
        CurveForm::Parabola {
            normal, axis_dir, ..
        }
        | CurveForm::Hyperbola {
            normal, axis_dir, ..
        } => (normal, axis_dir),
        _ => unreachable!("open conic"),
    };
    let v = unit(&normal.cross(axis));
    conic_arc(
        form,
        &point.plus(&v.scaled(&Real::from_i64(-1))),
        &point.plus(&v),
        dec,
    )
}
pub(super) fn line_local(
    origin: &Point3,
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<Chart, Refusal> {
    let v = direction;
    let from = origin.plus(&v.scaled(&Real::from_i64(-1)));
    let to = origin.plus(v);
    let inverse = lin(origin, &v.scaled(&(Real::one() / direction.norm_squared())));
    Ok(Chart {
        coords: segment(&from, &to, dec)?,
        inverse: Inverse::Affine {
            value: inverse,
            lo: Real::from_i64(-1),
            span: Real::from_i64(2),
        },
        checks: vec![],
        reverse: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TracePiece;

    #[test]
    fn algebraic_segment_taylor_matches_its_exact_line() {
        let mut dec = Decider::default();
        let root = Real::from_i64(2).sqrt();
        let start = Point3::new(root.clone(), Real::zero(), Real::zero());
        let end = Point3::new(root, Real::one(), Real::zero());
        let piece = TracePiece::new(segment(&start, &end, &mut dec).unwrap());
        let filtered = piece
            .taylor_interval(Interval::point(0.5).unwrap(), 1, &dec)
            .unwrap()
            .unwrap();
        let lo = Rational::from_f64_exact(filtered[0][0].lo()).unwrap();
        let hi = Rational::from_f64_exact(filtered[0][0].hi()).unwrap();
        assert!(lo >= r(0) && &lo * &lo <= r(2) && &hi * &hi >= r(2));

        for p in [64, 256] {
            let t = Ball::from_rational(&Rational::from_ratio(1, 2).unwrap(), p);
            let jets = piece.taylor(&t, 1, &dec).unwrap();
            assert!(jets[0][0].square().contains(&r(2)));
            assert!(jets[0][1].contains(&Rational::from_ratio(1, 2).unwrap()));
            assert!(jets[1][1].contains(&r(1)));
        }
    }

    #[test]
    fn shared_nested_jets_enclose_analytic_rational_derivatives() {
        let t = E::t();
        let root = (t.clone() + E::integer(1)).square().sqrt();
        let x = root.clone() / (t + E::integer(2));
        let piece = TracePiece::new([x.clone(), x.square(), root]);
        let dec = Decider::default();
        for parameter in [r(0), Rational::from_ratio(1, 2).unwrap(), r(1)] {
            let inverse = (&parameter + &r(2)).recip().unwrap();
            for order in [0, 1, 2, 6] {
                for precision in [64, 192] {
                    let jets = piece
                        .taylor(&Ball::from_rational(&parameter, precision), order, &dec)
                        .unwrap();
                    let x = &r(1) - &inverse;
                    let mut power = &inverse * &inverse;
                    for (k, values) in jets.iter().enumerate() {
                        let expected = if k == 0 {
                            [x.clone(), &x * &x, &parameter + &r(1)]
                        } else {
                            let sign = r(if k % 2 == 1 { 1 } else { -1 });
                            [
                                &sign * &power,
                                &sign * &(&r(2) * &power - &r(k as i64 + 1) * &(&power * &inverse)),
                                r(i64::from(k == 1)),
                            ]
                        };
                        for axis in 0..3 {
                            assert!(
                                values[axis].contains(&expected[axis]),
                                "t={parameter}, order={order}, precision={precision}, k={k}, axis={axis}"
                            );
                        }
                        if k > 0 {
                            power = &power * &inverse;
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn trace_enclosures_and_failures_are_independent_of_cache_warmth() {
        fn piece() -> TracePiece {
            let root = Real::from_i64(2).sqrt();
            let epsilon = Real::from_rational(Rational::parse("1e-100").unwrap());
            let divisor = &(&root + &epsilon) - &root;
            TracePiece::new([E::constant(Real::one() / divisor), E::t(), E::integer(0)])
        }
        let dec = Decider::default();
        let warm = piece();
        let value = Rational::parse("1e100").unwrap();
        let high = Ball::from_rational(&r(0), 512);
        assert!(warm.evaluate(&high, &dec).unwrap()[0].contains(&value));
        assert!(
            warm.taylor_interval(Interval::point(0.0).unwrap(), 0, &dec)
                .unwrap()
                .is_none()
        );
        for precision in [64, 256, 512] {
            let point = high.with_precision(precision);
            match (piece().evaluate(&point, &dec), warm.evaluate(&point, &dec)) {
                (Ok(cold), Ok(warm)) => {
                    for (a, b) in cold.iter().zip(&warm) {
                        assert_eq!((a.lo(), a.hi()), (b.lo(), b.hi()));
                    }
                    assert!(warm[0].contains(&value));
                    assert!(warm[1].contains(&r(0)));
                    assert!(warm[2].contains(&r(0)));
                }
                (Err(cold), Err(warm)) => {
                    assert_eq!(cold.code, "numbers.divisor-not-separated");
                    assert_eq!(warm.code, cold.code);
                }
                outcome => panic!("cache warmth changed the same trace query: {outcome:?}"),
            }
        }
        let mut expired = Decider::default();
        expired.set_deadline(Some(std::time::Instant::now()));
        assert_eq!(
            warm.evaluate(&high, &expired).unwrap_err().code,
            "numbers.deadline"
        );
    }

    #[test]
    fn binary64_taylor_filter_encloses_conic_derivatives_and_defers_singular_panels() {
        let dec = Decider::default();
        let (cosine, sine) = trig();
        let coords = [cosine, sine, E::t()];
        let program = JetProgram::new(&coords);
        for (lo, hi) in [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0)] {
            let filtered = interval_jets(&program, Interval::new(lo, hi).unwrap(), 2, &dec)
                .unwrap()
                .unwrap();
            for t in [lo, (lo + hi) * 0.5, hi] {
                let mut cache = JetValues::new(&program, 2);
                for (axis, root) in program.roots.iter().enumerate() {
                    let exact = cache
                        .jets(
                            *root,
                            &Real::from(Rational::from_f64_exact(t).unwrap()),
                            2,
                            0,
                            None,
                        )
                        .unwrap();
                    for (k, bounds) in filtered.iter().enumerate() {
                        let value = cache.coefficients[exact + k].rational_ref().unwrap();
                        assert!(Rational::from_f64_exact(bounds[axis].lo()).unwrap() <= *value);
                        assert!(Rational::from_f64_exact(bounds[axis].hi()).unwrap() >= *value);
                    }
                }
            }
        }
        let singular = [
            E::integer(1) / (E::t() - E::constant(Real::from(Rational::from_ratio(1, 2).unwrap()))),
            E::t(),
            E::integer(0),
        ];
        assert!(
            interval_jets(
                &JetProgram::new(&singular),
                Interval::new(0.25, 0.75).unwrap(),
                2,
                &dec,
            )
            .unwrap()
            .is_none()
        );
        let root = Real::from_i64(2).sqrt();
        let tiny = Real::from(Rational::parse("1e-50").unwrap());
        let ill_conditioned = [
            E::constant(Real::one() / (&(&root + &tiny) - &root)),
            E::t(),
            E::integer(0),
        ];
        assert!(
            interval_jets(
                &JetProgram::new(&ill_conditioned),
                Interval::point(0.0).unwrap(),
                0,
                &dec,
            )
            .unwrap()
            .is_none()
        );
        let mut expired = Decider::default();
        expired.set_deadline(Some(std::time::Instant::now()));
        assert_eq!(
            interval_jets(&program, Interval::point(0.0).unwrap(), 0, &expired)
                .unwrap_err()
                .code,
            "numbers.deadline"
        );
    }
}
