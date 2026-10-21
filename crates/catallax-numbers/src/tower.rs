//! Exact signs in real quadratic towers. Every generator is positive and its
//! radicand is positive and nonsquare in the preceding field. Consequently
//! a + b*sqrt(r) is zero exactly when both coefficients are zero.
use crate::{Rational, Real, Sign, decider::Work, real::Expr};
use std::{collections::HashMap, sync::Arc};

const MAX_LEVELS: usize = 8;
const MAX_STEPS: usize = 100_000;
const MAX_WORK: u64 = 1_000_000;
const LIMIT: &str = "quadratic-tower-limit";

// Lower-level elements need no padding when a later generator is added. Shared
// coefficient pairs make expression-cache hits and field embeddings cheap.
#[derive(Clone)]
enum Element {
    Rational(Rational),
    Pair {
        level: usize,
        coefficients: Arc<(Element, Element)>,
    },
}
impl Element {
    fn level(&self) -> usize {
        match self {
            Self::Rational(_) => 0,
            Self::Pair { level, .. } => *level,
        }
    }
    fn is_zero(&self) -> bool {
        matches!(self, Self::Rational(r) if r.is_zero())
    }
    fn pair(level: usize, a: Self, b: Self) -> Self {
        if b.is_zero() {
            a
        } else {
            Self::Pair {
                level,
                coefficients: Arc::new((a, b)),
            }
        }
    }
    fn split<'a>(&'a self, level: usize, zero: &'a Self) -> (&'a Self, &'a Self) {
        match self {
            Self::Pair {
                level: k,
                coefficients,
            } if *k == level => (&coefficients.0, &coefficients.1),
            _ => (self, zero),
        }
    }
}

struct Cost<'a> {
    work: &'a mut Work,
    steps: usize,
}
impl Cost<'_> {
    fn check(&self) -> Result<(), &'static str> {
        if self.steps >= MAX_STEPS {
            Err(LIMIT)
        } else {
            Ok(())
        }
    }
    fn step(&mut self) -> Result<(), &'static str> {
        self.check()?;
        self.steps += 1;
        self.work.charge(1)
    }
    fn exact(&mut self, a: &Rational, b: &Rational, mul: bool) -> Result<(), &'static str> {
        self.step()?;
        self.work.exact(a, b, mul)?;
        self.check()
    }
}

struct Tower {
    radicands: Vec<Element>,
    zero: Element,
}
impl Tower {
    fn add(&self, x: &Element, y: &Element, w: &mut Cost<'_>) -> Result<Element, &'static str> {
        w.step()?;
        if x.is_zero() {
            return Ok(y.clone());
        }
        if y.is_zero() {
            return Ok(x.clone());
        }
        if let (Element::Rational(a), Element::Rational(b)) = (x, y) {
            w.exact(a, b, false)?;
            return Ok(Element::Rational(a + b));
        }
        let level = x.level().max(y.level());
        let (a, b) = x.split(level, &self.zero);
        let (c, d) = y.split(level, &self.zero);
        Ok(Element::pair(level, self.add(a, c, w)?, self.add(b, d, w)?))
    }
    fn neg(x: &Element, w: &mut Cost<'_>) -> Result<Element, &'static str> {
        w.step()?;
        Ok(match x {
            Element::Rational(r) => {
                w.exact(r, r, false)?;
                Element::Rational(-r)
            }
            Element::Pair {
                level,
                coefficients,
            } => Element::pair(
                *level,
                Self::neg(&coefficients.0, w)?,
                Self::neg(&coefficients.1, w)?,
            ),
        })
    }
    fn sub(&self, x: &Element, y: &Element, w: &mut Cost<'_>) -> Result<Element, &'static str> {
        self.add(x, &Self::neg(y, w)?, w)
    }
    fn mul(&self, x: &Element, y: &Element, w: &mut Cost<'_>) -> Result<Element, &'static str> {
        w.step()?;
        if x.is_zero() || y.is_zero() {
            return Ok(self.zero.clone());
        }
        if let (Element::Rational(a), Element::Rational(b)) = (x, y) {
            w.exact(a, b, true)?;
            return Ok(Element::Rational(a * b));
        }
        let level = x.level().max(y.level());
        let (a, b) = x.split(level, &self.zero);
        let (c, d) = y.split(level, &self.zero);
        let ac = self.mul(a, c, w)?;
        let bd = self.mul(b, d, w)?;
        let bdr = self.mul(&bd, &self.radicands[level - 1], w)?;
        let ad = self.mul(a, d, w)?;
        let bc = self.mul(b, c, w)?;
        Ok(Element::pair(
            level,
            self.add(&ac, &bdr, w)?,
            self.add(&ad, &bc, w)?,
        ))
    }
    fn norm(
        &self,
        level: usize,
        a: &Element,
        b: &Element,
        w: &mut Cost<'_>,
    ) -> Result<Element, &'static str> {
        let aa = self.mul(a, a, w)?;
        let bb = self.mul(b, b, w)?;
        let bbr = self.mul(&bb, &self.radicands[level - 1], w)?;
        self.sub(&aa, &bbr, w)
    }
    fn inv(&self, x: &Element, w: &mut Cost<'_>) -> Result<Option<Element>, &'static str> {
        w.step()?;
        match x {
            Element::Rational(r) => {
                w.exact(r, r, true)?;
                Ok(r.recip().map(Element::Rational))
            }
            Element::Pair {
                level,
                coefficients,
            } => {
                let (a, b) = &**coefficients;
                let norm = self.norm(*level, a, b, w)?;
                let Some(inv) = self.inv(&norm, w)? else {
                    return Ok(None);
                };
                Ok(Some(Element::pair(
                    *level,
                    self.mul(a, &inv, w)?,
                    Self::neg(&self.mul(b, &inv, w)?, w)?,
                )))
            }
        }
    }
    fn sign(&self, x: &Element, w: &mut Cost<'_>) -> Result<Sign, &'static str> {
        w.step()?;
        match x {
            Element::Rational(r) => Ok(r.sign()),
            Element::Pair {
                level,
                coefficients,
            } => {
                let (a, b) = &**coefficients;
                let sa = self.sign(a, w)?;
                let sb = self.sign(b, w)?;
                if sa == Sign::Zero || sa == sb {
                    return Ok(sb);
                }
                let norm = self.norm(*level, a, b, w)?;
                Ok(match self.sign(&norm, w)? {
                    Sign::Positive => sa,
                    Sign::Negative => sb,
                    Sign::Zero => Sign::Zero,
                })
            }
        }
    }
    // Solve (p+q√r)^2=a+b√r in the preceding field. When b=0,
    // either q=0 or p=0. Otherwise (p²-q²r)²=a²-b²r and
    // p²=(a ± sqrt(a²-b²r))/2. Both choices must be tested.
    fn sqrt_in(
        &self,
        x: &Element,
        level: usize,
        w: &mut Cost<'_>,
    ) -> Result<Option<Element>, &'static str> {
        w.step()?;
        if level == 0 {
            let Element::Rational(r) = x else {
                unreachable!()
            };
            let root = crate::quadratic::rational_sqrt(r, w.work)?;
            w.check()?;
            return Ok(root.map(Element::Rational));
        }
        let (a, b) = x.split(level, &self.zero);
        let r = &self.radicands[level - 1];
        if b.is_zero() {
            if let Some(root) = self.sqrt_in(a, level - 1, w)? {
                return Ok(Some(root));
            }
            let inv = self.inv(r, w)?.expect("positive radicand");
            let quotient = self.mul(a, &inv, w)?;
            return Ok(self
                .sqrt_in(&quotient, level - 1, w)?
                .map(|q| Element::pair(level, self.zero.clone(), q)));
        }
        let norm = self.norm(level, a, b, w)?;
        let Some(n) = self.sqrt_in(&norm, level - 1, w)? else {
            return Ok(None);
        };
        for positive in [true, false] {
            let sum = if positive {
                self.add(a, &n, w)?
            } else {
                self.sub(a, &n, w)?
            };
            let half = Element::Rational(Rational::one().mul_pow2(-1));
            let p2 = self.mul(&sum, &half, w)?;
            let Some(p) = self.sqrt_in(&p2, level - 1, w)? else {
                continue;
            };
            if p.is_zero() {
                continue;
            }
            let two_p = self.add(&p, &p, w)?;
            let inv = self.inv(&two_p, w)?.expect("nonzero square root");
            let q = self.mul(b, &inv, w)?;
            let candidate = Element::pair(level, p, q);
            // Check the reconstructed square before choosing the real embedding.
            let square = self.mul(&candidate, &candidate, w)?;
            if !self.sub(&square, x, w)?.is_zero() {
                continue;
            }
            return Ok(Some(if self.sign(&candidate, w)? == Sign::Negative {
                Self::neg(&candidate, w)?
            } else {
                candidate
            }));
        }
        Ok(None)
    }
    fn sqrt(&mut self, x: Element, w: &mut Cost<'_>) -> Result<Option<Element>, &'static str> {
        match self.sign(&x, w)? {
            Sign::Negative => return Ok(None),
            Sign::Zero => return Ok(Some(self.zero.clone())),
            Sign::Positive => {}
        }
        let level = self.radicands.len();
        if let Some(root) = self.sqrt_in(&x, level, w)? {
            return Ok(Some(root));
        }
        if level == MAX_LEVELS {
            return Ok(None);
        }
        self.radicands.push(x);
        Ok(Some(Element::pair(
            level + 1,
            self.zero.clone(),
            Element::Rational(Rational::one()),
        )))
    }
    // Constructible quartics split into quadratics over a real quadratic
    // extension. The even depressed form is immediate; otherwise a positive
    // rational root of the cubic resolvent supplies that extension. Select
    // the real embedding through the certified open isolating interval.
    fn quartic(
        &mut self,
        root: &crate::Algebraic,
        polynomial: &crate::UPoly,
        w: &mut Cost<'_>,
    ) -> Result<Option<Element>, &'static str> {
        let coefficients = polynomial.coefficients();
        let limbs = coefficients
            .iter()
            .map(|r| {
                r.numerator()
                    .bit_length()
                    .max(r.denominator().bit_length())
                    .div_ceil(64)
                    .max(1)
            })
            .max()
            .unwrap();
        // The fixed-degree depression uses fewer than 64 rational operations;
        // intermediate coefficient lengths are at most eight input lengths.
        w.work
            .charge(8192u64.saturating_mul(limbs).saturating_mul(limbs))?;
        let leading = &coefficients[4];
        let a = &coefficients[3] / leading;
        let b = &coefficients[2] / leading;
        let c = &coefficients[1] / leading;
        let d = &coefficients[0] / leading;
        let a2 = &a * &a;
        let a3 = &a2 * &a;
        let linear = &(&c - &(&a * &b).mul_pow2(-1)) + &a3.mul_pow2(-3);
        let p = &b - &(&a2 * &Rational::from_i64(3)).mul_pow2(-3);
        let constant = &(&(&d - &(&a * &c).mul_pow2(-2)) + &(&a2 * &b).mul_pow2(-4))
            - &(&(&a2 * &a2) * &Rational::from_i64(3)).mul_pow2(-8);
        let shift = Element::Rational(a.mul_pow2(-2));
        let half = Element::Rational(Rational::one().mul_pow2(-1));
        let discriminant = &(&p * &p) - &constant.mul_pow2(2);
        if !linear.is_zero() {
            // (y²+s*y+t)(y²-s*y+u), z=s²:
            // z³+2p*z²+(p²-4r)*z-q²=0, t,u=(p+z ∓ q/s)/2.
            let resolvent = crate::UPoly::new(vec![
                -(&linear * &linear),
                discriminant,
                p.mul_pow2(1),
                Rational::one(),
            ]);
            if !crate::algebraic::rational_root_possible(&resolvent, w.work)? {
                return Ok(None);
            }
            for candidate in crate::algebraic::isolate(&resolvent, w.work)? {
                let Some(z) = candidate
                    .value
                    .as_rational()
                    .filter(|z| z.sign() == Sign::Positive)
                else {
                    continue;
                };
                let Some(s) = self.sqrt(Element::Rational(z.clone()), w)? else {
                    continue;
                };
                for reverse in [false, true] {
                    let s = if reverse {
                        Self::neg(&s, w)?
                    } else {
                        s.clone()
                    };
                    let inverse = self.inv(&s, w)?.expect("positive resolvent root");
                    let quotient = self.mul(&Element::Rational(linear.clone()), &inverse, w)?;
                    let sum = self.sub(&Element::Rational(&p + z), &quotient, w)?;
                    let t = self.mul(&sum, &half, w)?;
                    let four_t = self.mul(&t, &Element::Rational(Rational::from_i64(4)), w)?;
                    let disc = self.sub(&Element::Rational(z.clone()), &four_t, w)?;
                    let Some(delta) = self.sqrt(disc, w)? else {
                        continue;
                    };
                    for negative in [true, false] {
                        let delta = if negative {
                            Self::neg(&delta, w)?
                        } else {
                            delta.clone()
                        };
                        let numerator = self.sub(&delta, &s, w)?;
                        let value = self.mul(&numerator, &half, w)?;
                        let value = self.sub(&value, &shift, w)?;
                        if self.matches_root(&value, root, polynomial, w)? {
                            return Ok(Some(value));
                        }
                    }
                }
            }
            return Ok(None);
        }
        let Some(delta) = self.sqrt(Element::Rational(discriminant), w)? else {
            return Ok(None);
        };
        for lower in [true, false] {
            let negative_p = Element::Rational(-&p);
            let numerator = if lower {
                self.sub(&negative_p, &delta, w)?
            } else {
                self.add(&negative_p, &delta, w)?
            };
            let square = self.mul(&numerator, &half, w)?;
            let Some(magnitude) = self.sqrt(square, w)? else {
                continue;
            };
            for negative in [true, false] {
                let value = if negative {
                    Self::neg(&magnitude, w)?
                } else {
                    magnitude.clone()
                };
                let value = self.sub(&value, &shift, w)?;
                if self.matches_root(&value, root, polynomial, w)? {
                    return Ok(Some(value));
                }
            }
        }
        Ok(None)
    }
    fn matches_root(
        &self,
        value: &Element,
        root: &crate::Algebraic,
        polynomial: &crate::UPoly,
        w: &mut Cost<'_>,
    ) -> Result<bool, &'static str> {
        let (lo, hi) = root.isolating_interval();
        let above = self.sub(value, &Element::Rational(lo), w)?;
        let below = self.sub(&Element::Rational(hi), value, w)?;
        if self.sign(&above, w)? != Sign::Positive || self.sign(&below, w)? != Sign::Positive {
            return Ok(false);
        }
        let mut residual = self.zero.clone();
        for coefficient in polynomial.coefficients().iter().rev() {
            residual = self.mul(&residual, value, w)?;
            residual = self.add(&residual, &Element::Rational(coefficient.clone()), w)?;
        }
        Ok(residual.is_zero())
    }
    fn visit(
        &mut self,
        x: &Real,
        w: &mut Cost<'_>,
        cache: &mut HashMap<*const crate::real::Node, Option<Element>>,
    ) -> Result<Option<Element>, &'static str> {
        w.step()?;
        let key = Arc::as_ptr(&x.node);
        if let Some(v) = cache.get(&key) {
            return Ok(v.clone());
        }
        // Only a defined rational cache may hide an expression subtree.
        let result = if x.node.fast_defined
            && let Some(r) = &x.node.rational
        {
            w.exact(r, r, false)?;
            Some(Element::Rational(r.clone()))
        } else {
            match &x.node.expr {
                Expr::Rational(r) => Some(Element::Rational(r.clone())),
                Expr::Algebraic(a) => {
                    let p = a.defining_polynomial();
                    if p.degree() == Some(4) {
                        let result = self.quartic(a, &p, w)?;
                        cache.insert(key, result.clone());
                        return Ok(result);
                    }
                    if p.degree() != Some(2) {
                        return Ok(None);
                    }
                    w.work.poly(&p, &p)?;
                    w.check()?;
                    let c = p.coefficients();
                    let disc = &(&c[1] * &c[1]) - &(&(&c[0] * &c[2]) * &Rational::from_i64(4));
                    let (lo, _) = a.isolating_interval();
                    // At the lower endpoint the polynomial has the opposite
                    // sign to its derivative at the isolated simple root.
                    let lower_sign = w.work.eval_poly(&p, &lo)?.sign();
                    w.check()?;
                    let denominator = c[2].mul_pow2(1);
                    let factor = if lower_sign == Sign::Negative {
                        Rational::one()
                    } else {
                        -Rational::one()
                    };
                    let Some(root) = self.sqrt(Element::Rational(disc), w)? else {
                        return Ok(None);
                    };
                    let scaled = self.mul(&root, &Element::Rational(&factor / &denominator), w)?;
                    Some(self.add(&scaled, &Element::Rational(-&c[1] / denominator), w)?)
                }
                Expr::Coordinate(..) => None,
                Expr::Neg(a) => match self.visit(a, w, cache)? {
                    Some(a) => Some(Self::neg(&a, w)?),
                    None => None,
                },
                Expr::Sqrt(a) => match self.visit(a, w, cache)? {
                    Some(a) => self.sqrt(a, w)?,
                    None => None,
                },
                Expr::Square(a) | Expr::Pow(a, _) => {
                    let Some(mut a) = self.visit(a, w, cache)? else {
                        return Ok(None);
                    };
                    let mut n = match &x.node.expr {
                        Expr::Pow(_, n) => *n,
                        _ => 2,
                    };
                    let mut result = Element::Rational(Rational::one());
                    while n != 0 {
                        if n & 1 != 0 {
                            result = self.mul(&result, &a, w)?;
                        }
                        n >>= 1;
                        if n != 0 {
                            a = self.mul(&a, &a, w)?;
                        }
                    }
                    Some(result)
                }
                Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                    let Some(a) = self.visit(a, w, cache)? else {
                        return Ok(None);
                    };
                    let Some(b) = self.visit(b, w, cache)? else {
                        return Ok(None);
                    };
                    match &x.node.expr {
                        Expr::Add(..) => Some(self.add(&a, &b, w)?),
                        Expr::Sub(..) => Some(self.sub(&a, &b, w)?),
                        Expr::Mul(..) => Some(self.mul(&a, &b, w)?),
                        _ => match self.inv(&b, w)? {
                            Some(inv) => Some(self.mul(&a, &inv, w)?),
                            None => None,
                        },
                    }
                }
            }
        };
        cache.insert(key, result.clone());
        Ok(result)
    }
}

fn has_quartic(x: &Real, w: &mut Cost<'_>) -> Result<bool, &'static str> {
    let mut pending = vec![x];
    let mut seen = std::collections::HashSet::new();
    while let Some(x) = pending.pop() {
        w.step()?;
        if !seen.insert(Arc::as_ptr(&x.node)) || (x.node.fast_defined && x.node.rational.is_some())
        {
            continue;
        }
        match &x.node.expr {
            Expr::Algebraic(a) if matches!(&a.repr, crate::algebraic::AlgebraicRepr::Root(root) if root.polynomial.degree() == Some(4)) =>
            {
                return Ok(true);
            }
            Expr::Neg(a) | Expr::Sqrt(a) | Expr::Square(a) | Expr::Pow(a, _) => pending.push(a),
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                pending.extend([a, b])
            }
            _ => {}
        }
    }
    Ok(false)
}

pub(crate) fn sign(x: &Real, work: &mut Work) -> Result<Option<Sign>, &'static str> {
    // Temporarily narrow the work limit, so every reservation is checked before
    // arithmetic. Restore it even when this attempt declines or hits a deadline.
    let decision_limit = work.policy.max_work;
    let allowance = MAX_WORK.min(decision_limit.saturating_sub(work.used) / 8);
    work.policy.max_work = work.used + allowance;
    let mut cost = Cost { work, steps: 0 };
    let mut tower = Tower {
        radicands: Vec::new(),
        zero: Element::Rational(Rational::zero()),
    };
    let attempt = (|| {
        // Polynomial reduction in one quartic generator avoids reconstructing
        // that generator just to evaluate its own defining identities. Pure
        // quadratic-tower inputs still go directly to tower arithmetic.
        if has_quartic(x, &mut cost)?
            && let Some(value) = crate::generator::rational(x, cost.work)?
        {
            return Ok(Some(value.sign()));
        }
        let value = tower.visit(x, &mut cost, &mut HashMap::new())?;
        value.map(|value| tower.sign(&value, &mut cost)).transpose()
    })();
    work.policy.max_work = decision_limit;
    match attempt {
        Err(LIMIT | "numbers.budget-exhausted") => Ok(None),
        other => other,
    }
}
