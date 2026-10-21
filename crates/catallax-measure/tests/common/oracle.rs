//! Independent GMP rational / MPFR directed-rounding oracle. No kernel Ball calls.
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::Rational;
use rug::{
    Float, Rational as Q,
    float::{Constant, Round},
};
use std::ops::{Add, Div, Mul, Sub};

#[derive(Clone, Debug)]
pub enum Expr {
    Rational(Q),
    Pi,
    Sqrt(u32),
    Add(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
}
pub fn rat(n: i64, d: i64) -> Expr {
    Expr::Rational(Q::from((n, d)))
}
pub fn from_kernel(r: &Rational) -> Expr {
    Expr::Rational(to_q(r))
}
pub fn pi() -> Expr {
    Expr::Pi
}
pub fn sqrt(n: u32) -> Expr {
    Expr::Sqrt(n)
}
pub fn to_q(r: &Rational) -> Q {
    Q::from_str_radix(&r.to_string(), 10).unwrap()
}
pub fn to_kernel(q: &Q) -> Rational {
    Rational::parse(&q.to_string()).unwrap()
}

macro_rules! binary {
    ($trait:ident, $method:ident, $variant:ident) => {
        impl $trait for Expr {
            type Output = Self;
            fn $method(self, rhs: Self) -> Self {
                Self::$variant(Box::new(self), Box::new(rhs))
            }
        }
    };
}
binary!(Add, add, Add);
binary!(Mul, mul, Mul);
binary!(Div, div, Div);
impl Sub for Expr {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + rat(-1, 1) * rhs
    }
}

#[derive(Clone, Debug)]
pub struct Bounds {
    pub lo: Q,
    pub hi: Q,
}
impl Bounds {
    fn point(q: Q) -> Self {
        Self {
            lo: q.clone(),
            hi: q,
        }
    }
    fn product(self, rhs: Self) -> Self {
        let ends = [
            self.lo.clone() * &rhs.lo,
            self.lo * &rhs.hi,
            self.hi.clone() * rhs.lo,
            self.hi * rhs.hi,
        ];
        Self {
            lo: ends.iter().min().unwrap().clone(),
            hi: ends.iter().max().unwrap().clone(),
        }
    }
}
impl Expr {
    pub fn bounds(&self, bits: u32) -> Bounds {
        match self {
            Self::Rational(q) => Bounds::point(q.clone()),
            Self::Pi => Bounds {
                lo: Float::with_val_round(bits, Constant::Pi, Round::Down)
                    .0
                    .to_rational()
                    .unwrap(),
                hi: Float::with_val_round(bits, Constant::Pi, Round::Up)
                    .0
                    .to_rational()
                    .unwrap(),
            },
            Self::Sqrt(n) => {
                let mut lo = Float::with_val(bits, *n);
                let mut hi = lo.clone();
                lo.sqrt_round(Round::Down);
                hi.sqrt_round(Round::Up);
                Bounds {
                    lo: lo.to_rational().unwrap(),
                    hi: hi.to_rational().unwrap(),
                }
            }
            Self::Add(a, b) => {
                let a = a.bounds(bits);
                let b = b.bounds(bits);
                Bounds {
                    lo: a.lo + b.lo,
                    hi: a.hi + b.hi,
                }
            }
            Self::Mul(a, b) => a.bounds(bits).product(b.bounds(bits)),
            Self::Div(a, b) => {
                let b = b.bounds(bits);
                assert!(b.lo > 0 || b.hi < 0, "oracle divisor crosses zero");
                a.bounds(bits).product(Bounds {
                    lo: b.hi.recip(),
                    hi: b.lo.recip(),
                })
            }
        }
    }
}

pub fn width(e: &Enclosure) -> Rational {
    &e.hi - &e.lo
}
pub fn meets(e: &Enclosure, a: &Accuracy) -> bool {
    let m = if e.lo > Rational::zero() || e.hi < Rational::zero() {
        e.lo.abs().min(e.hi.abs())
    } else {
        Rational::zero()
    };
    width(e) <= (&a.relative * &m).max(a.absolute.clone())
}

pub fn contains(e: &Enclosure, expected: &Expr) {
    let lo = to_q(&e.lo);
    let hi = to_q(&e.hi);
    assert!(lo <= hi, "reversed enclosure: {e:?}");
    // Narrow the independent reference to < 1/1024 of the returned width.
    // Exact rational answers, including zero centroids, need no approximate comparison.
    let result_width = hi.clone() - &lo;
    let mut bits = 256;
    loop {
        let b = expected.bounds(bits);
        assert!(
            lo <= b.hi && b.lo <= hi,
            "closed form outside enclosure at {bits} bits: {e:?}"
        );
        if b.lo == b.hi {
            assert!(lo <= b.lo && b.hi <= hi);
            return;
        }
        if b.lo >= lo && b.hi <= hi && (b.hi - b.lo) * 1024 < result_width {
            return;
        }
        assert!(
            bits < 32768,
            "independent oracle could not certify containment: {e:?}"
        );
        bits *= 2;
    }
}

pub fn exact(e: &Enclosure, expected: &Expr) {
    contains(e, expected);
    assert_eq!(e.lo, e.hi, "rational-only computation must be exact");
}
