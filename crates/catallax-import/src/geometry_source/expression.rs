//! Differentiable parameter expressions with interval evaluation.
use super::mapping::MappingError;
use catallax_numbers::{Ball, Decider, Rational, Real, RealPoly, Refusal};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub(super) struct Expr(Arc<Node>);
#[derive(Debug)]
enum Node {
    Constant(Real),
    Pi,
    Parameter(usize),
    Add(Expr, Expr),
    Mul(Expr, Expr),
    Div(Expr, Expr),
    Sqrt(Expr),
    Sin(Expr),
    Cos(Expr),
    Exp(Expr),
}
impl Expr {
    pub(super) fn constant(value: Real) -> Self {
        Self(Arc::new(Node::Constant(value)))
    }
    pub(super) fn rational(value: Rational) -> Self {
        Self::constant(Real::from_rational(value))
    }
    pub(super) fn integer(value: i64) -> Self {
        Self::constant(Real::from_i64(value))
    }
    pub(super) fn pi() -> Self {
        Self(Arc::new(Node::Pi))
    }
    pub(super) fn parameter(axis: usize) -> Self {
        Self(Arc::new(Node::Parameter(axis)))
    }
    pub(super) fn add(&self, b: &Self) -> Self {
        Self(Arc::new(Node::Add(self.clone(), b.clone())))
    }
    pub(super) fn neg(&self) -> Self {
        self.mul(&Self::integer(-1))
    }
    pub(super) fn sub(&self, b: &Self) -> Self {
        self.add(&b.neg())
    }
    pub(super) fn mul(&self, b: &Self) -> Self {
        Self(Arc::new(Node::Mul(self.clone(), b.clone())))
    }
    pub(super) fn div(&self, b: &Self) -> Self {
        Self(Arc::new(Node::Div(self.clone(), b.clone())))
    }
    pub(super) fn sqrt(&self) -> Self {
        Self(Arc::new(Node::Sqrt(self.clone())))
    }
    pub(super) fn sin(&self) -> Self {
        Self(Arc::new(Node::Sin(self.clone())))
    }
    pub(super) fn cos(&self) -> Self {
        Self(Arc::new(Node::Cos(self.clone())))
    }
    pub(super) fn exp(&self) -> Self {
        Self(Arc::new(Node::Exp(self.clone())))
    }
    pub(super) fn pow(&self, n: u32) -> Self {
        (0..n).fold(Self::integer(1), |p, _| p.mul(self))
    }
    pub(super) fn polynomial(
        p: &RealPoly,
        vars: &[Self],
        dec: &mut Decider,
    ) -> Result<Self, MappingError> {
        let mut out = Self::integer(0);
        for (powers, c) in p.terms() {
            dec.charge_work(1)?;
            let mut term = Self::constant(c);
            for (v, n) in vars.iter().zip(powers.exponents()) {
                term = term.mul(&v.pow(*n));
            }
            out = out.add(&term);
        }
        Ok(out)
    }
    pub(super) fn derivative(&self, axis: usize, dec: &mut Decider) -> Result<Self, MappingError> {
        dec.charge_work(1)?;
        let d = |e: &Self, dec: &mut Decider| e.derivative(axis, dec);
        Ok(match self.0.as_ref() {
            Node::Constant(_) | Node::Pi => Self::integer(0),
            Node::Parameter(i) => Self::integer(i64::from(*i == axis)),
            Node::Add(a, b) => d(a, dec)?.add(&d(b, dec)?),
            Node::Mul(a, b) => d(a, dec)?.mul(b).add(&a.mul(&d(b, dec)?)),
            Node::Div(a, b) => d(a, dec)?.mul(b).sub(&a.mul(&d(b, dec)?)).div(&b.pow(2)),
            Node::Sqrt(a) => d(a, dec)?.div(&Self::integer(2).mul(self)),
            Node::Sin(a) => d(a, dec)?.mul(&a.cos()),
            Node::Cos(a) => d(a, dec)?.mul(&a.sin()).neg(),
            Node::Exp(a) => d(a, dec)?.mul(self),
        })
    }
    pub(super) fn evaluate(
        &self,
        parameters: &[Ball],
        bits: u32,
        dec: &mut Decider,
    ) -> Result<Ball, MappingError> {
        dec.charge_work(1)?;
        let e = |a: &Self, dec: &mut Decider| a.evaluate(parameters, bits, dec);
        let domain = || {
            Refusal::new(
                "import.geometry-parameter-domain",
                "parameter interval contains a singularity; subdivide the source chart",
            )
        };
        Ok(match self.0.as_ref() {
            Node::Constant(v) => dec.enclose_at(v, bits)?,
            Node::Pi => Ball::pi_checked(bits, dec)?,
            Node::Parameter(i) => parameters.get(*i).cloned().ok_or_else(|| {
                Refusal::new(
                    "import.geometry-parameter-dimension",
                    "missing source parameter",
                )
            })?,
            Node::Add(a, b) => e(a, dec)? + e(b, dec)?,
            Node::Mul(a, b) => e(a, dec)? * e(b, dec)?,
            Node::Div(a, b) => e(a, dec)?.checked_div(&e(b, dec)?).ok_or_else(domain)?,
            Node::Sqrt(a) => {
                let v = e(a, dec)?;
                if v.lo() < Rational::zero() {
                    return Err(domain().into());
                }
                v.sqrt().ok_or_else(domain)?
            }
            Node::Sin(a) => e(a, dec)?.sin(),
            Node::Cos(a) => e(a, dec)?.cos(),
            Node::Exp(a) => e(a, dec)?.exp(),
        })
    }
}
pub(super) type Triple = [Expr; 3];
pub(super) fn add(a: &Triple, b: &Triple) -> Triple {
    std::array::from_fn(|i| a[i].add(&b[i]))
}
pub(super) fn sub(a: &Triple, b: &Triple) -> Triple {
    std::array::from_fn(|i| a[i].sub(&b[i]))
}
pub(super) fn scale(a: &Triple, t: &Expr) -> Triple {
    a.each_ref().map(|v| v.mul(t))
}
pub(super) fn dot(a: &Triple, b: &Triple) -> Expr {
    (0..3).fold(Expr::integer(0), |s, i| s.add(&a[i].mul(&b[i])))
}
pub(super) fn cross(a: &Triple, b: &Triple) -> Triple {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        a[j].mul(&b[k]).sub(&a[k].mul(&b[j]))
    })
}
pub(super) fn unit(a: &Triple) -> Triple {
    scale(a, &Expr::integer(1).div(&dot(a, a).sqrt()))
}
pub(super) fn derivative(
    a: &Triple,
    axis: usize,
    dec: &mut Decider,
) -> Result<Triple, MappingError> {
    Ok([
        a[0].derivative(axis, dec)?,
        a[1].derivative(axis, dec)?,
        a[2].derivative(axis, dec)?,
    ])
}
