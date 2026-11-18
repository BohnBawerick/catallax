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
/// One step of a [`Program`]; operands are earlier steps.
#[derive(Clone, Copy, Debug)]
enum Step {
    /// A value and the bound on what rounding it to binary64 lost.
    Constant(f64, f64),
    Parameter(usize),
    Add(usize, usize),
    Mul(usize, usize),
    Div(usize, usize),
    Sqrt(usize),
    Sin(usize),
    Cos(usize),
    Exp(usize),
}

/// Three expressions as one binary64 straight-line program.
///
/// A chart's expressions share most of their nodes, and the tree walk of
/// [`Expr::evaluate`] visits a shared node once for every use, enclosing each
/// constant again at the requested precision. The program evaluates every
/// node once, in binary64, with its constants rounded when it is compiled.
///
/// Each value carries a running bound on its rounding error, to first order
/// and with rounding to nearest. It tells the caller when binary64 has lost
/// the value to cancellation. It is not an enclosure and proves nothing.
#[derive(Clone, Debug)]
pub(super) struct Program {
    steps: Vec<Step>,
    outputs: [usize; 3],
}
impl Program {
    pub(super) fn compile(coordinates: &Triple, dec: &mut Decider) -> Result<Self, MappingError> {
        let mut steps = Vec::new();
        let mut seen = std::collections::BTreeMap::new();
        let mut outputs = [0; 3];
        for (output, expression) in outputs.iter_mut().zip(coordinates) {
            *output = Self::step(expression, &mut steps, &mut seen, dec)?;
        }
        Ok(Self { steps, outputs })
    }

    fn step(
        expression: &Expr,
        steps: &mut Vec<Step>,
        seen: &mut std::collections::BTreeMap<usize, usize>,
        dec: &mut Decider,
    ) -> Result<usize, MappingError> {
        // A shared node is one allocation, so its address names it.
        let node = Arc::as_ptr(&expression.0) as usize;
        if let Some(index) = seen.get(&node) {
            return Ok(*index);
        }
        dec.charge_work(1)?;
        let mut operand = |a: &Expr, dec: &mut Decider| Self::step(a, steps, seen, dec);
        let step = match expression.0.as_ref() {
            Node::Constant(value) => {
                let rounded = dec.enclose_at(value, 64)?.mid().to_f64_nearest();
                let exact = value
                    .as_rational()
                    .is_some_and(|value| Rational::from_f64_exact(rounded) == Some(value));
                Step::Constant(
                    rounded,
                    if exact {
                        0.0
                    } else {
                        f64::EPSILON * rounded.abs()
                    },
                )
            }
            Node::Pi => Step::Constant(std::f64::consts::PI, f64::EPSILON * 4.0),
            Node::Parameter(axis) => Step::Parameter(*axis),
            Node::Add(a, b) => Step::Add(operand(a, dec)?, operand(b, dec)?),
            Node::Mul(a, b) => Step::Mul(operand(a, dec)?, operand(b, dec)?),
            Node::Div(a, b) => Step::Div(operand(a, dec)?, operand(b, dec)?),
            Node::Sqrt(a) => Step::Sqrt(operand(a, dec)?),
            Node::Sin(a) => Step::Sin(operand(a, dec)?),
            Node::Cos(a) => Step::Cos(operand(a, dec)?),
            Node::Exp(a) => Step::Exp(operand(a, dec)?),
        };
        steps.push(step);
        seen.insert(node, steps.len() - 1);
        Ok(steps.len() - 1)
    }

    /// Steps one evaluation runs.
    pub(super) fn len(&self) -> usize {
        self.steps.len()
    }

    /// The three values at `parameters`, each with its rounding bound. A
    /// missing parameter, a division by what may be zero or a root of what may
    /// be negative gives a bound that is not finite.
    pub(super) fn evaluate(&self, parameters: &[f64]) -> [(f64, f64); 3] {
        const ROUNDING: f64 = f64::EPSILON;
        let mut values: Vec<(f64, f64)> = Vec::with_capacity(self.steps.len());
        for step in &self.steps {
            let (value, bound) = match *step {
                Step::Constant(value, bound) => (value, bound),
                Step::Parameter(axis) => (parameters.get(axis).copied().unwrap_or(f64::NAN), 0.0),
                Step::Add(a, b) => {
                    let ((a, ea), (b, eb)) = (values[a], values[b]);
                    let sum = a + b;
                    (sum, ea + eb + ROUNDING * sum.abs())
                }
                Step::Mul(a, b) => {
                    let ((a, ea), (b, eb)) = (values[a], values[b]);
                    let product = a * b;
                    (
                        product,
                        a.abs() * eb + b.abs() * ea + ea * eb + ROUNDING * product.abs(),
                    )
                }
                Step::Div(a, b) => {
                    let ((a, ea), (b, eb)) = (values[a], values[b]);
                    let quotient = a / b;
                    let least = b.abs() - eb;
                    let bound = if least > 0.0 {
                        (ea + quotient.abs() * eb) / least + ROUNDING * quotient.abs()
                    } else {
                        f64::INFINITY
                    };
                    (quotient, bound)
                }
                Step::Sqrt(a) => {
                    let (a, ea) = values[a];
                    let root = a.sqrt();
                    let bound = if ea == 0.0 {
                        ROUNDING * root
                    } else if a - ea > 0.0 {
                        ea / (2.0 * (a - ea).sqrt()) + ROUNDING * root
                    } else {
                        f64::INFINITY
                    };
                    (root, bound)
                }
                Step::Sin(a) => (values[a].0.sin(), values[a].1 + 2.0 * ROUNDING),
                Step::Cos(a) => (values[a].0.cos(), values[a].1 + 2.0 * ROUNDING),
                Step::Exp(a) => {
                    let (a, ea) = values[a];
                    let power = a.exp();
                    (power, power * ea.exp_m1() + 2.0 * ROUNDING * power)
                }
            };
            values.push(if value.is_finite() && bound.is_finite() {
                (value, bound)
            } else {
                (value, f64::INFINITY)
            });
        }
        self.outputs.map(|output| values[output])
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

#[cfg(test)]
mod tests {
    use super::*;

    fn program(coordinates: &Triple) -> Program {
        Program::compile(coordinates, &mut Decider::default()).unwrap()
    }

    #[test]
    fn a_shared_node_is_one_step_of_the_binary64_program() {
        let sine = Expr::parameter(0).sin();
        let shared = sine.mul(&sine).add(&sine);
        // The parameter, its sine, the product and the sum, once for all three.
        assert_eq!(program(&[shared.clone(), shared.clone(), shared]).len(), 4);
    }

    #[test]
    fn the_binary64_program_agrees_with_the_enclosure_within_its_bound() {
        let (u, v) = (Expr::parameter(0), Expr::parameter(1));
        let radius = Expr::rational(Rational::parse("7/3").unwrap());
        let coordinates = [
            radius.mul(&u.cos()).add(&v.pow(3)),
            radius.mul(&u.sin()).div(&v.mul(&v).add(&Expr::integer(1))),
            u.mul(&u).add(&Expr::pi()).sqrt().sub(&v.exp()),
        ];
        let compiled = program(&coordinates);
        let mut dec = Decider::default();
        for parameters in [[0.0, 0.0], [0.3, -1.25], [5.9, 2.0], [-17.0, 0.001]] {
            let balls = parameters
                .map(|value| Ball::from_rational(&Rational::from_f64_exact(value).unwrap(), 96));
            let rounded = compiled.evaluate(&parameters);
            for (coordinate, (value, bound)) in coordinates.iter().zip(rounded) {
                let enclosure = coordinate.evaluate(&balls, 96, &mut dec).unwrap();
                let miss = (value - enclosure.mid().to_f64_nearest()).abs();
                assert!(bound < 1e-13, "{parameters:?}: bound {bound}");
                assert!(miss <= bound, "{parameters:?}: {miss} over {bound}");
            }
        }
    }

    #[test]
    fn the_binary64_program_reports_what_cancellation_loses() {
        let large = Expr::rational(Rational::parse("100000000000000000").unwrap());
        let lost = Expr::parameter(0).add(&large).sub(&large);
        let [(value, bound), ..] = program(&[lost.clone(), lost.clone(), lost]).evaluate(&[1.0]);
        // The exact value is 1; binary64 cannot hold 1e17 + 1.
        assert!((value - 1.0).abs() <= bound, "{value} within {bound}");
        assert!(bound > 1.0);
    }

    #[test]
    fn a_singular_point_has_no_finite_bound() {
        let one = Expr::integer(1);
        let pole = one.div(&Expr::parameter(0));
        let root = Expr::parameter(0).sub(&one).sqrt();
        let outputs = program(&[pole.clone(), root, pole]).evaluate(&[0.0]);
        assert!(outputs.iter().all(|(_, bound)| !bound.is_finite()));
    }
}
