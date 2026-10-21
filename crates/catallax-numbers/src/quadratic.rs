//! Exact identities in a single real quadratic field, with rational fallback.
use crate::{Rational, Real, Sign, decider::Work, real::Expr};
use std::{collections::HashMap, sync::Arc};

#[derive(Clone)]
struct Value {
    a: Rational,
    b: Rational,
    d: Rational,
}
impl Value {
    fn rational(a: Rational) -> Self {
        Self {
            a,
            b: Rational::zero(),
            d: Rational::zero(),
        }
    }
    fn normalize(mut self) -> Self {
        if self.b.is_zero() {
            self.d = Rational::zero();
        }
        self
    }
    fn reserve(&self, other: &Self, work: &mut Work) -> Result<(), &'static str> {
        let limbs = [&self.a, &self.b, &self.d, &other.a, &other.b, &other.d]
            .into_iter()
            .map(|r| {
                r.numerator()
                    .bit_length()
                    .max(r.denominator().bit_length())
                    .div_ceil(64)
                    .max(1)
            })
            .max()
            .unwrap();
        // Composite field operations reserve 512 L² logical work units, with L
        // the largest input numerator/denominator limb count. Rational-square
        // tests additionally charge their actual operands through Work::exact.
        // This fixed price is deterministic and precedes the arithmetic.
        work.charge(512u64.saturating_mul(limbs).saturating_mul(limbs))
    }
    fn combine(&self, rhs: &Self, op: u8, work: &mut Work) -> Result<Option<Self>, &'static str> {
        self.reserve(rhs, work)?;
        let d = if self.b.is_zero() { &rhs.d } else { &self.d };
        let scaled = if !self.b.is_zero() && !rhs.b.is_zero() && self.d != rhs.d {
            // Q(sqrt(d)) = Q(sqrt(e)) exactly when e/d is a rational square.
            let Some(scale) = rational_sqrt(&(&rhs.d / &self.d), work)? else {
                return Ok(None);
            };
            Some(&rhs.b * &scale)
        } else {
            None
        };
        let a = &self.a;
        let b = &self.b;
        let c = &rhs.a;
        let e = scaled.as_ref().unwrap_or(&rhs.b);
        let (a, b) = match op {
            0 => (a + c, b + e),
            1 => (a - c, b - e),
            2 => (&(a * c) + &(&(b * e) * d), a * e + b * c),
            _ => {
                let norm = &(c * c) - &(&(e * e) * d);
                if norm.is_zero() {
                    return Ok(None);
                }
                (
                    &(&(a * c) - &(&(b * e) * d)) / &norm,
                    (b * c - a * e) / norm,
                )
            }
        };
        Ok(Some(Self { a, b, d: d.clone() }.normalize()))
    }
    fn sign(&self) -> Sign {
        if self.b.is_zero() {
            return self.a.sign();
        }
        if self.a.sign() == self.b.sign() || self.a.is_zero() {
            return self.b.sign();
        }
        match (&self.a * &self.a).cmp(&(&(&self.b * &self.b) * &self.d)) {
            std::cmp::Ordering::Equal => Sign::Zero,
            std::cmp::Ordering::Greater => self.a.sign(),
            std::cmp::Ordering::Less => self.b.sign(),
        }
    }
    fn sqrt(&self, work: &mut Work) -> Result<Option<Self>, &'static str> {
        self.reserve(self, work)?;
        if self.b.is_zero() {
            if self.a.sign() == Sign::Negative {
                return Ok(None);
            }
            return Ok(Some(match rational_sqrt(&self.a, work)? {
                Some(r) => Self::rational(r),
                None => Self {
                    a: Rational::zero(),
                    b: Rational::one(),
                    d: self.a.clone(),
                },
            }));
        }
        // If sqrt(a+b√d)=p+q√d, its norm has rational square
        // a²-b²d=(p²-q²d)² and p²=(a±sqrt(norm))/2.
        // Test both choices exactly and select the nonnegative square root.
        let norm = &(&self.a * &self.a) - &(&(&self.b * &self.b) * &self.d);
        let Some(norm) = rational_sqrt(&norm, work)? else {
            return Ok(None);
        };
        for n in [norm.clone(), -norm] {
            let Some(p) = rational_sqrt(&(&self.a + &n).mul_pow2(-1), work)? else {
                continue;
            };
            if p.is_zero() {
                continue;
            }
            let q = &self.b / &p.mul_pow2(1);
            let mut result = Self {
                a: p,
                b: q,
                d: self.d.clone(),
            };
            if result.sign() == Sign::Negative {
                result.a = -result.a;
                result.b = -result.b;
            }
            return Ok(Some(result));
        }
        Ok(None)
    }
}
pub(crate) fn rational_sqrt(
    r: &Rational,
    work: &mut Work,
) -> Result<Option<Rational>, &'static str> {
    work.exact(r, r, true)?;
    if r.sign() == Sign::Negative {
        return Ok(None);
    }
    let n = r.numerator().0.sqrt();
    let d = r.denominator().0.sqrt();
    if &n * &n != r.numerator().0 || &d * &d != r.denominator().0 {
        return Ok(None);
    }
    Ok(Some(
        Rational::new(crate::Integer(n), crate::Integer(d)).unwrap(),
    ))
}

// Every operation below is an identity in Q(√d), d>0 nonsquare. In particular
// zero means BOTH coefficients vanish, never that a numeric residual is small.
// An unsupported extension declines and leaves the general separation proof
// available. The caller must establish expression definedness first, so these
// identities cannot erase a hidden zero divisor or negative square root.
// This cache is owned by one attempt inside one decision, not by Real nodes.
// Every visit, including a memoized visit, costs one work unit in addition to
// the composite-operation reservations above.
pub(crate) fn sign(x: &Real, work: &mut Work) -> Result<Option<Sign>, &'static str> {
    let Some(value) = visit(x, work, &mut HashMap::new())? else {
        return Ok(None);
    };
    value.reserve(&value, work)?;
    Ok(Some(value.sign()))
}
pub(crate) fn zero(x: &Real, work: &mut Work) -> Result<bool, &'static str> {
    Ok(sign(x, work)? == Some(Sign::Zero))
}
fn visit(
    x: &Real,
    work: &mut Work,
    cache: &mut HashMap<*const crate::real::Node, Option<Value>>,
) -> Result<Option<Value>, &'static str> {
    work.charge(1)?;
    let key = Arc::as_ptr(&x.node);
    if let Some(v) = cache.get(&key) {
        return Ok(v.clone());
    }
    let result = if let Some(r) = &x.node.rational {
        work.exact(r, r, false)?;
        Some(Value::rational(r.clone()))
    } else {
        match &x.node.expr {
            Expr::Rational(_) => unreachable!(),
            Expr::Algebraic(a) => {
                let p = a.defining_polynomial();
                if p.degree() == Some(2) {
                    let c = p.coefficients();
                    work.poly(&p, &p)?;
                    let disc = &(&c[1] * &c[1]) - &(&(&c[0] * &c[2]) * &Rational::from_i64(4));
                    let denominator = c[2].mul_pow2(1);
                    let (lo, _) = a.isolating_interval();
                    // On an isolating interval the lower endpoint sign is
                    // opposite to the derivative sign at its simple root.
                    // For ax²+bx+c, 2ax+b is ±sqrt(b²-4ac), so this
                    // selects the exact root without an approximate threshold.
                    let lower_sign = work.eval_poly(&p, &lo)?.sign();
                    let factor = if lower_sign == Sign::Negative {
                        Rational::one()
                    } else {
                        -Rational::one()
                    };
                    match Value::rational(disc).sqrt(work)? {
                        Some(root) => {
                            let scaled = root
                                .combine(&Value::rational(&factor / &denominator), 2, work)?
                                .unwrap();
                            scaled.combine(&Value::rational(-&c[1] / denominator), 0, work)?
                        }
                        None => None,
                    }
                } else {
                    None
                }
            }
            Expr::Coordinate(..) => None,
            Expr::Neg(a) => visit(a, work, cache)?.map(|v| Value {
                a: -v.a,
                b: -v.b,
                d: v.d,
            }),
            Expr::Sqrt(a) => match visit(a, work, cache)? {
                Some(a) => a.sqrt(work)?,
                None => None,
            },
            Expr::Square(a) | Expr::Pow(a, _) => {
                if let Some(mut a) = visit(a, work, cache)? {
                    let mut n = match &x.node.expr {
                        Expr::Pow(_, n) => *n,
                        _ => 2,
                    };
                    let mut result = Value::rational(Rational::one());
                    while n != 0 {
                        if n & 1 != 0 {
                            result = result.combine(&a, 2, work)?.unwrap();
                        }
                        n >>= 1;
                        if n != 0 {
                            a = a.combine(&a, 2, work)?.unwrap();
                        }
                    }
                    Some(result)
                } else {
                    None
                }
            }
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                let a = visit(a, work, cache)?;
                let b = visit(b, work, cache)?;
                match a.zip(b) {
                    Some((a, b)) => a.combine(
                        &b,
                        match &x.node.expr {
                            Expr::Add(..) => 0,
                            Expr::Sub(..) => 1,
                            Expr::Mul(..) => 2,
                            _ => 3,
                        },
                        work,
                    )?,
                    None => None,
                }
            }
        }
    };
    cache.insert(key, result.clone());
    Ok(result)
}
