//! Immutable shared real expressions with conservative construction-time caches.

use crate::{Algebraic, Ball, Interval, IsolatedRoot, Rational};
use std::sync::Arc;

/// Exact real expression. Cloning shares its immutable definition and caches.
#[derive(Clone, Debug)]
pub struct Real {
    pub(crate) node: Arc<Node>,
}
#[derive(Clone, Debug)]
pub(crate) enum Expr {
    Rational(Rational),
    Algebraic(Algebraic),
    Coordinate(Arc<IsolatedRoot>, usize),
    Add(Real, Real),
    Sub(Real, Real),
    Mul(Real, Real),
    Div(Real, Real),
    Neg(Real),
    Sqrt(Real),
    Square(Real),
    Pow(Real, u32),
}
#[derive(Debug)]
pub(crate) struct Node {
    pub expr: Expr,
    pub approx: Interval,
    pub coarse: Option<Ball>,
    pub rational: Option<Rational>,
    pub fast_defined: bool,
}
impl Real {
    /// Zero.
    pub fn zero() -> Self {
        Self::from_i64(0)
    }
    /// One.
    pub fn one() -> Self {
        Self::from_i64(1)
    }
    /// Exact rational leaf.
    pub fn from_rational(r: Rational) -> Self {
        let approx = Interval::from_rational(&r);
        Self::make(Expr::Rational(r.clone()), approx, None, Some(r), true)
    }
    /// Exact integer leaf.
    pub fn from_i64(v: i64) -> Self {
        Self::from_rational(Rational::from_i64(v))
    }
    /// Exact algebraic leaf.
    pub fn from_algebraic(a: Algebraic) -> Self {
        let (lo, hi) = a.isolating_interval();
        let approx = Interval::from_rational(&lo).hull(&Interval::from_rational(&hi));
        let coarse = Some(Ball::bounds(lo, hi, 1));
        let rational = a.as_rational().cloned();
        Self::make(Expr::Algebraic(a), approx, coarse, rational, true)
    }
    /// Coordinate of a certified system root; panics for an invalid axis.
    pub fn from_root_coordinate(root: &IsolatedRoot, axis: usize) -> Self {
        assert!(axis < root.root_box().dimension());
        let lo = root.root_box().lo()[axis].clone();
        let hi = root.root_box().hi()[axis].clone();
        let approx = Interval::from_rational(&lo).hull(&Interval::from_rational(&hi));
        Self::make(
            Expr::Coordinate(Arc::new(root.clone()), axis),
            approx,
            Some(Ball::bounds(lo, hi, 1)),
            None,
            true,
        )
    }
    /// Lazy square root. Domain validation belongs to the decider.
    pub fn sqrt(&self) -> Self {
        let approx = self.node.approx.sqrt().unwrap_or(Interval::entire());
        Self::make(
            Expr::Sqrt(self.clone()),
            approx,
            self.node.coarse.as_ref().and_then(Ball::sqrt),
            None,
            self.node.fast_defined && self.node.approx.lo() >= 0.0,
        )
    }
    /// Lazy square, retaining the entire operand definition.
    pub fn square(&self) -> Self {
        Self::make(
            Expr::Square(self.clone()),
            self.node.approx.square(),
            self.node.coarse.as_ref().map(Ball::square),
            self.node.rational.as_ref().map(|r| r * r),
            self.node.fast_defined,
        )
    }
    /// Lazy non-negative integer power, retaining definedness even for exponent zero.
    pub fn pow_u32(&self, exp: u32) -> Self {
        let mut a = self.node.approx;
        let mut r = Interval::point(1.0).unwrap();
        let mut n = exp;
        while n != 0 {
            if n & 1 != 0 {
                r = r * a;
            }
            n >>= 1;
            if n != 0 {
                a = a.square();
            }
        }
        let rational = self.node.rational.as_ref().map(|r| {
            let mut a = r.clone();
            let mut v = Rational::one();
            let mut n = exp;
            while n != 0 {
                if n & 1 != 0 {
                    v = &v * &a;
                }
                n >>= 1;
                if n != 0 {
                    a = &a * &a;
                }
            }
            v
        });
        Self::make(
            Expr::Pow(self.clone(), exp),
            r,
            self.node.coarse.as_ref().map(|b| b.pow_u32(exp)),
            rational,
            self.node.fast_defined,
        )
    }
    /// Exact value of a structurally rational expression; does not refine.
    pub fn as_rational(&self) -> Option<Rational> {
        self.node.rational.clone()
    }
    /// Cached, conservative double enclosure.
    pub fn approx(&self) -> Interval {
        self.node.approx
    }
    /// Double approximation for display only.
    pub fn to_f64(&self) -> f64 {
        self.node
            .rational
            .as_ref()
            .map_or_else(|| self.approx().mid(), Rational::to_f64_nearest)
    }
    fn make(
        expr: Expr,
        approx: Interval,
        coarse: Option<Ball>,
        rational: Option<Rational>,
        fast_defined: bool,
    ) -> Self {
        let coarse = rational
            .as_ref()
            .map(|r| Ball::from_rational(r, 1))
            .or(coarse)
            .or_else(|| {
                if fast_defined && approx.is_finite() {
                    Some(Ball::bounds(
                        Rational::from_f64_exact(approx.lo()).unwrap(),
                        Rational::from_f64_exact(approx.hi()).unwrap(),
                        1,
                    ))
                } else {
                    None
                }
            });
        Self {
            node: Arc::new(Node {
                expr,
                approx,
                coarse,
                rational,
                fast_defined,
            }),
        }
    }
    fn add(&self, rhs: &Self) -> Self {
        let coarse = self
            .node
            .coarse
            .as_ref()
            .zip(rhs.node.coarse.as_ref())
            .map(|(a, b)| a + b);
        let rational = self
            .node
            .rational
            .as_ref()
            .zip(rhs.node.rational.as_ref())
            .map(|(a, b)| a + b);
        Self::make(
            Expr::Add(self.clone(), rhs.clone()),
            self.approx() + rhs.approx(),
            coarse,
            rational,
            self.node.fast_defined && rhs.node.fast_defined,
        )
    }
    fn sub(&self, rhs: &Self) -> Self {
        let coarse = self
            .node
            .coarse
            .as_ref()
            .zip(rhs.node.coarse.as_ref())
            .map(|(a, b)| a - b);
        let rational = self
            .node
            .rational
            .as_ref()
            .zip(rhs.node.rational.as_ref())
            .map(|(a, b)| a - b);
        Self::make(
            Expr::Sub(self.clone(), rhs.clone()),
            self.approx() - rhs.approx(),
            coarse,
            rational,
            self.node.fast_defined && rhs.node.fast_defined,
        )
    }
    fn mul(&self, rhs: &Self) -> Self {
        let coarse = self
            .node
            .coarse
            .as_ref()
            .zip(rhs.node.coarse.as_ref())
            .map(|(a, b)| a * b);
        let rational = self
            .node
            .rational
            .as_ref()
            .zip(rhs.node.rational.as_ref())
            .map(|(a, b)| a * b);
        Self::make(
            Expr::Mul(self.clone(), rhs.clone()),
            self.approx() * rhs.approx(),
            coarse,
            rational,
            self.node.fast_defined && rhs.node.fast_defined,
        )
    }
    fn div(&self, rhs: &Self) -> Self {
        let coarse = self
            .node
            .coarse
            .as_ref()
            .zip(rhs.node.coarse.as_ref())
            .and_then(|(a, b)| a.checked_div(b));
        let rational = self
            .node
            .rational
            .as_ref()
            .zip(rhs.node.rational.as_ref())
            .and_then(|(a, b)| a.checked_div(b));
        Self::make(
            Expr::Div(self.clone(), rhs.clone()),
            self.approx() / rhs.approx(),
            coarse,
            rational,
            self.node.fast_defined && rhs.node.fast_defined && rhs.approx().sign().is_some(),
        )
    }
    fn neg(&self) -> Self {
        Self::make(
            Expr::Neg(self.clone()),
            -self.approx(),
            self.node.coarse.as_ref().map(|b| -b),
            self.node.rational.as_ref().map(|r| -r),
            self.node.fast_defined,
        )
    }
}
binops!(Real; Add::add, Sub::sub, Mul::mul, Div::div);
negop!(Real);
impl From<Rational> for Real {
    fn from(r: Rational) -> Self {
        Self::from_rational(r)
    }
}
impl From<i64> for Real {
    fn from(v: i64) -> Self {
        Self::from_i64(v)
    }
}
impl From<Algebraic> for Real {
    fn from(a: Algebraic) -> Self {
        Self::from_algebraic(a)
    }
}
