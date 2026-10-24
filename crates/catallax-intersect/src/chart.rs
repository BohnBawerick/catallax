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
    rc::Rc,
    sync::{Arc, Mutex, OnceLock},
};

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
        if p == INTERVAL_PRECISION
            && let Some(value) = cache.interval.get()
        {
            return Ok(*value);
        }
        let value = <Ball as Scalar>::constant(x, cache, p, dec)?.to_interval();
        if p == INTERVAL_PRECISION {
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
fn jets<S: Scalar>(
    e: &E,
    t: &S,
    n: usize,
    p: u32,
    dec: Option<&Decider>,
    cache: &mut HashMap<usize, Rc<[S]>>,
) -> Result<Rc<[S]>, JetError> {
    check_deadline(dec)?;
    let key = Arc::as_ptr(&e.0) as usize;
    if let Some(j) = cache.get(&key) {
        return Ok(j.clone());
    }
    let mut out = vec![S::integer(0, p); n + 1];
    match e.0.as_ref() {
        Node::T => {
            out[0] = t.clone();
            if n > 0 {
                out[1] = S::integer(1, p);
            }
        }
        Node::Constant(x, b) => out[0] = S::constant(x, b, p, dec)?,
        Node::Add(a, b) => {
            let a = jets(a, t, n, p, dec, cache)?;
            let b = jets(b, t, n, p, dec, cache)?;
            for k in 0..=n {
                check_deadline(dec)?;
                out[k] = a[k].plus(&b[k]);
            }
        }
        Node::Neg(a) => {
            let a = jets(a, t, n, p, dec, cache)?;
            for k in 0..=n {
                check_deadline(dec)?;
                out[k] = S::integer(0, p).minus(&a[k]);
            }
        }
        Node::Mul(a, b) => {
            let square = Arc::ptr_eq(&a.0, &b.0);
            let a = jets(a, t, n, p, dec, cache)?;
            let b = jets(b, t, n, p, dec, cache)?;
            let a_order = a.iter().rposition(|v| !v.is_zero()).unwrap_or(0);
            let b_order = b.iter().rposition(|v| !v.is_zero()).unwrap_or(0);
            let mut products = if square {
                vec![None; a_order + 1]
            } else {
                Vec::new()
            };
            for k in 0..=a_order.saturating_add(b_order).min(n) {
                for i in k.saturating_sub(b_order)..=k.min(a_order) {
                    check_deadline(dec)?;
                    if square {
                        let index = i.min(k - i);
                        if i == index {
                            products[index] = (!a[i].is_zero() && !a[k - i].is_zero())
                                .then(|| a[i].times(&a[k - i]));
                        }
                        if let Some(product) = &products[i.min(k - i)] {
                            out[k] = out[k].plus(product);
                        }
                    } else if !a[i].is_zero() && !b[k - i].is_zero() {
                        out[k] = out[k].plus(&a[i].times(&b[k - i]));
                    }
                }
            }
        }
        Node::Div(a, b) => {
            let a = jets(a, t, n, p, dec, cache)?;
            let b = jets(b, t, n, p, dec, cache)?;
            let inverse = S::integer(1, p).divide(&b[0]).ok_or(JetError::Split)?;
            let order = b.iter().rposition(|v| !v.is_zero()).unwrap_or(0);
            for k in 0..=n {
                let mut v = a[k].clone();
                for i in 1..=k.min(order) {
                    check_deadline(dec)?;
                    if !b[i].is_zero() && !out[k - i].is_zero() {
                        v = v.minus(&b[i].times(&out[k - i]));
                    }
                }
                out[k] = if v.is_zero() { v } else { v.times(&inverse) };
            }
        }
        Node::Sqrt(a) => {
            let a = jets(a, t, n, p, dec, cache)?;
            out[0] = a[0].root().ok_or(JetError::Split)?;
            let inverse = if n == 0 {
                None
            } else {
                Some(
                    S::integer(1, p)
                        .divide(&S::integer(2, p).times(&out[0]))
                        .ok_or(JetError::Split)?,
                )
            };
            let mut products = Vec::with_capacity(n / 2);
            for k in 1..=n {
                let mut v = a[k].clone();
                products.clear();
                for i in 1..k {
                    check_deadline(dec)?;
                    if i <= k / 2 {
                        products.push(
                            (!out[i].is_zero() && !out[k - i].is_zero())
                                .then(|| out[i].times(&out[k - i])),
                        );
                    }
                    // Reuse each symmetric product, but retain every
                    // subtraction in its original order and precision.
                    if let Some(product) = &products[i.min(k - i) - 1] {
                        v = v.minus(product);
                    }
                }
                out[k] = if v.is_zero() {
                    v
                } else {
                    v.times(inverse.as_ref().unwrap())
                };
            }
        }
    }
    let out: Rc<[S]> = out.into();
    cache.insert(key, out.clone());
    Ok(out)
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
    Circle { cos: Linear, sin: Linear },
    Squared { value: Linear, lo: Real, span: Real },
    Affine { value: Linear, lo: Real, span: Real },
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
        let mut cache = HashMap::new();
        let j = self
            .coords
            .each_ref()
            .map(|e| jets(e, t, n, 0, None, &mut cache).unwrap());
        (0..=n)
            .map(|k| [j[0][k].clone(), j[1][k].clone(), j[2][k].clone()])
            .collect()
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
    let parameter = E::constant(from.clone()) + E::constant(to - from) * E::t();
    let polynomial = |p: &catallax_numbers::RealPoly| {
        let terms = p.terms();
        let Some((monomial, coefficient)) = terms.last() else {
            return E::integer(0);
        };
        let mut power = monomial.exponents()[0];
        let mut expression = E::constant(coefficient.clone());
        for (monomial, coefficient) in terms[..terms.len() - 1].iter().rev() {
            let next = monomial.exponents()[0];
            for _ in next..power {
                expression = expression * parameter.clone();
            }
            expression = expression + E::constant(coefficient.clone());
            power = next;
        }
        for _ in 0..power {
            expression = expression * parameter.clone();
        }
        expression
    };
    let [x, y, z, w] = curve
        .span_polynomials(span, dec)?
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
    coords: &[E; 3],
    t: &Ball,
    n: usize,
    dec: &Decider,
) -> Result<Vec<[Ball; 3]>, JetError> {
    let mut cache = HashMap::new();
    let a = jets(&coords[0], t, n, t.precision(), Some(dec), &mut cache)?;
    let b = jets(&coords[1], t, n, t.precision(), Some(dec), &mut cache)?;
    let c = jets(&coords[2], t, n, t.precision(), Some(dec), &mut cache)?;
    Ok((0..=n)
        .map(|k| [a[k].clone(), b[k].clone(), c[k].clone()])
        .collect())
}
pub(super) fn interval_jets(
    coords: &[E; 3],
    t: Interval,
    n: usize,
    dec: &Decider,
) -> Result<Option<Vec<[Interval; 3]>>, Refusal> {
    let mut cache = HashMap::new();
    let result = coords
        .each_ref()
        .map(|e| jets(e, &t, n, INTERVAL_PRECISION, Some(dec), &mut cache));
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
        .map(|k| [values[0][k], values[1][k], values[2][k]])
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
    coords: &[E; 3],
    t: &Ball,
    n: usize,
    dec: &Decider,
) -> Result<Vec<[Ball; 3]>, Refusal> {
    dec.check_deadline()?;
    let p = t.precision();
    let wide = t.with_precision(p.saturating_add(64));
    match enclose_jets(coords, &wide, n, dec) {
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
    let a = ball_jets(coords, &left, n, dec)?;
    let b = ball_jets(coords, &right, n, dec)?;
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
        for (lo, hi) in [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0)] {
            let filtered = interval_jets(&coords, Interval::new(lo, hi).unwrap(), 2, &dec)
                .unwrap()
                .unwrap();
            for t in [lo, (lo + hi) * 0.5, hi] {
                let mut cache = HashMap::new();
                for (axis, expression) in coords.iter().enumerate() {
                    let exact = jets(
                        expression,
                        &Real::from(Rational::from_f64_exact(t).unwrap()),
                        2,
                        0,
                        None,
                        &mut cache,
                    )
                    .unwrap();
                    for k in 0..=2 {
                        let value = exact[k].rational_ref().unwrap();
                        assert!(
                            Rational::from_f64_exact(filtered[k][axis].lo()).unwrap() <= *value
                        );
                        assert!(
                            Rational::from_f64_exact(filtered[k][axis].hi()).unwrap() >= *value
                        );
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
            interval_jets(&singular, Interval::new(0.25, 0.75).unwrap(), 2, &dec)
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
            interval_jets(&ill_conditioned, Interval::point(0.0).unwrap(), 0, &dec)
                .unwrap()
                .is_none()
        );
        let mut expired = Decider::default();
        expired.set_deadline(Some(std::time::Instant::now()));
        assert_eq!(
            interval_jets(&coords, Interval::point(0.0).unwrap(), 0, &expired)
                .unwrap_err()
                .code,
            "numbers.deadline"
        );
    }
}
