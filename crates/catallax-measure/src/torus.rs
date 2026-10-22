//! Torus forms on lifted angular charts. The boundary correction records the
//! holonomy of noncontractible loops; the scalar area fixes the whole-torus term.

use crate::integrate::{Builder, Expr, Program};
use crate::planar::{components, direction};
use crate::surface::Property;
use catallax_numbers::Rational;
use catallax_shapes::{RPoint3, Torus};
use std::collections::BTreeMap;

pub(crate) struct Programs {
    pub(crate) boundary: Program,
    pub(crate) angular_derivative: Program,
    pub(crate) coordinates: Program,
    pub(crate) correction: Program,
    pub(crate) parallel: Program,
    pub(crate) full: Program,
    pub(crate) indices: Vec<usize>,
}

#[derive(Clone)]
struct Polynomial<'a> {
    terms: BTreeMap<[u8; 4], Expr<'a>>,
    builder: &'a Builder,
}

impl<'a> Polynomial<'a> {
    fn constant(builder: &'a Builder, value: Expr<'a>) -> Self {
        Self {
            terms: BTreeMap::from([([0; 4], value)]),
            builder,
        }
    }
    fn variable(builder: &'a Builder, i: usize) -> Self {
        let mut powers = [0; 4];
        powers[i] = 1;
        Self {
            terms: BTreeMap::from([(powers, builder.k(1))]),
            builder,
        }
    }
    fn plus(&self, rhs: &Self) -> Self {
        let mut result = self.clone();
        for (powers, coefficient) in &rhs.terms {
            result
                .terms
                .entry(*powers)
                .and_modify(|c| *c = *c + *coefficient)
                .or_insert(*coefficient);
        }
        result
    }
    fn times(&self, rhs: &Self) -> Self {
        let mut result = Self {
            terms: BTreeMap::new(),
            builder: self.builder,
        };
        for (a, ca) in &self.terms {
            for (b, cb) in &rhs.terms {
                let powers = std::array::from_fn(|i| a[i] + b[i]);
                let coefficient = *ca * *cb;
                result
                    .terms
                    .entry(powers)
                    .and_modify(|c| *c = *c + coefficient)
                    .or_insert(coefficient);
            }
        }
        result
    }
    fn scaled(&self, value: Expr<'a>) -> Self {
        Self {
            terms: self.terms.iter().map(|(p, c)| (*p, *c * value)).collect(),
            builder: self.builder,
        }
    }
}

fn power<'a>(b: &'a Builder, base: Expr<'a>, exponent: u8) -> Expr<'a> {
    (0..exponent).fold(b.k(1), |p, _| p * base)
}

fn mean(a: u8, b: u8) -> Rational {
    let (numerator, denominator) = match (a, b) {
        (0, 0) => (1, 1),
        (2, 0) | (0, 2) => (1, 2),
        (4, 0) | (0, 4) => (3, 8),
        (2, 2) => (1, 8),
        _ if a % 2 == 1 || b % 2 == 1 => (0, 1),
        _ => unreachable!("torus flux polynomial degree exceeds four"),
    };
    Rational::from_ratio(numerator, denominator).unwrap()
}
// Integral of cos(t)^a sin(t)^b, anchored at t=0. The first component
// multiplies t; the remaining terms are a periodic polynomial in cos and sin.
fn primitive(a: u8, b: u8) -> (Rational, Vec<(u8, u8, Rational)>) {
    let r = |n, d| Rational::from_ratio(n, d).unwrap();
    let terms: Vec<(u8, u8, Rational)> = match (a, b) {
        (0, 0) => vec![],
        (1, 0) => vec![(0, 1, r(1, 1))],
        (0, 1) => vec![(1, 0, r(-1, 1))],
        (2, 0) => vec![(1, 1, r(1, 2))],
        (0, 2) => vec![(1, 1, r(-1, 2))],
        (1, 1) => vec![(0, 2, r(1, 2))],
        (3, 0) => vec![(0, 1, r(1, 1)), (0, 3, r(-1, 3))],
        (0, 3) => vec![(1, 0, r(-1, 1)), (3, 0, r(1, 3))],
        (2, 1) => vec![(3, 0, r(-1, 3))],
        (1, 2) => vec![(0, 3, r(1, 3))],
        (4, 0) => vec![(1, 1, r(1, 2)), (3, 1, r(1, 8)), (1, 3, r(-1, 8))],
        (0, 4) => vec![(1, 1, r(-1, 2)), (3, 1, r(1, 8)), (1, 3, r(-1, 8))],
        (3, 1) => vec![(4, 0, r(-1, 4))],
        (1, 3) => vec![(0, 4, r(1, 4))],
        (2, 2) => vec![(3, 1, r(-1, 8)), (1, 3, r(1, 8))],
        _ => unreachable!("torus flux polynomial degree exceeds four"),
    };
    let mean = mean(a, b);
    let at_zero = terms
        .iter()
        .filter(|(_, s, _)| *s == 0)
        .fold(Rational::zero(), |sum, (_, _, c)| &sum + c);
    let mut terms = terms;
    if !at_zero.is_zero() {
        terms.push((0, 0, -at_zero));
    }
    (mean, terms)
}

fn evaluate_primitive<'a>(
    b: &'a Builder,
    exponents: [u8; 2],
    cosine: Expr<'a>,
    sine: Expr<'a>,
    angle: Expr<'a>,
) -> Expr<'a> {
    let (mean, terms) = primitive(exponents[0], exponents[1]);
    terms
        .into_iter()
        .fold(b.r(mean) * angle, |sum, (c, s, coefficient)| {
            sum + b.r(coefficient) * power(b, cosine, c) * power(b, sine, s)
        })
}

fn dot<'a>(a: [Expr<'a>; 3], b: [Expr<'a>; 3]) -> Expr<'a> {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn compile(torus: &Torus, anchor: &RPoint3, property: Property) -> Programs {
    let b = Builder::default();
    let q =
        std::array::from_fn::<_, 3, _>(|i| b.input(i) - b.r(components(&torus.center)[i].clone()));
    let dq = std::array::from_fn::<_, 3, _>(|i| b.input(i + 3));
    let unit = |v: &catallax_shapes::RVector3| {
        let norm = b.r(v.norm_squared()).sqrt();
        direction(v).map(|c| b.r(c.clone()) / norm)
    };
    let axis = unit(&torus.axis);
    let first = unit(&torus.ref_dir);
    let second = [
        axis[1] * first[2] - axis[2] * first[1],
        axis[2] * first[0] - axis[0] * first[2],
        axis[0] * first[1] - axis[1] * first[0],
    ];
    let x = dot(q, first);
    let y = dot(q, second);
    let z = dot(q, axis);
    let dx = dot(dq, first);
    let dy = dot(dq, second);
    let dz = dot(dq, axis);
    let rho = (x.square() + y.square()).sqrt();
    let major = b.r(torus.major.clone());
    let minor = b.r(torus.minor.clone());
    let cu = x / rho;
    let su = y / rho;
    let cv = (rho - major) / minor;
    let sv = z / minor;
    let du = (x * dy - y * dx) / rho.square();
    let drho = (x * dx + y * dy) / rho;
    let dv = ((rho - major) * dz - z * drho) / minor.square();
    let variables = std::array::from_fn::<_, 4, _>(|i| Polynomial::variable(&b, i));
    let constant = |e| Polynomial::constant(&b, e);
    let er = std::array::from_fn::<_, 3, _>(|i| {
        variables[0]
            .scaled(first[i])
            .plus(&variables[1].scaled(second[i]))
    });
    let radial = constant(major).plus(&variables[2].scaled(minor));
    let area = radial.scaled(minor);
    let normal = std::array::from_fn::<_, 3, _>(|i| {
        er[i]
            .times(&variables[2])
            .plus(&variables[3].scaled(axis[i]))
    });
    let origin = torus.center.minus(anchor);
    let origin = direction(&origin).map(|c| b.r(c.clone()));
    let position = std::array::from_fn::<_, 3, _>(|i| {
        constant(origin[i])
            .plus(&radial.times(&er[i]))
            .plus(&variables[3].scaled(minor * axis[i]))
    });
    let flux = normal[0]
        .times(&position[0])
        .plus(&normal[1].times(&position[1]))
        .plus(&normal[2].times(&position[2]));
    let volume = area.times(&flux).scaled(b.k(1) / b.k(3));
    let moments = std::array::from_fn::<_, 3, _>(|i| {
        area.times(&position[i])
            .times(&position[i])
            .times(&normal[i])
            .scaled(b.k(1) / b.k(2))
    });
    let polynomials = [
        area,
        volume,
        moments[0].clone(),
        moments[1].clone(),
        moments[2].clone(),
    ];
    let indices = match property {
        Property::Area => vec![0],
        Property::Volume => vec![0, 1],
        Property::Centroid => vec![0, 1, 2, 3, 4],
    };
    let mut forms = Vec::new();
    let mut corrections = Vec::new();
    let mut parallels = Vec::with_capacity(indices.len());
    for index in &indices {
        let mut h = b.k(0);
        let mut g = b.k(0);
        let mut parallel = b.k(0);
        for (powers, coefficient) in &polynomials[*index].terms {
            let u_polynomial = power(&b, cu, powers[0]) * power(&b, su, powers[1]);
            let v_primitive = evaluate_primitive(&b, [powers[2], powers[3]], cv, sv, b.input(6));
            h = h + *coefficient * u_polynomial * v_primitive;
            let v_mean = mean(powers[2], powers[3]);
            if !v_mean.is_zero() {
                let u_primitive =
                    evaluate_primitive(&b, [powers[0], powers[1]], cu, su, b.input(6));
                g = g + *coefficient * b.r(v_mean) * u_primitive;
            }
            let u_mean = mean(powers[0], powers[1]);
            if !u_mean.is_zero() {
                let v_primitive = evaluate_primitive(
                    &b,
                    [powers[2], powers[3]],
                    b.input(7),
                    b.input(8),
                    b.input(6),
                );
                parallel = parallel + *coefficient * b.r(u_mean) * v_primitive;
            }
        }
        forms.push(-h * du);
        corrections.push(g);
        parallels.push(-b.k(2) * b.input(9) * parallel);
    }
    let whole_area = b.k(4) * b.input(6).square() * major * minor;
    let whole_volume = b.k(2) * b.input(6).square() * major * minor.square();
    let full = [
        whole_area,
        whole_volume,
        whole_volume * origin[0],
        whole_volume * origin[1],
        whole_volume * origin[2],
    ];
    Programs {
        boundary: b.finish(&forms),
        angular_derivative: b.finish(&[dv]),
        coordinates: b.finish(&[cu, su, cv, sv]),
        correction: b.finish(&corrections),
        parallel: b.finish(&parallels),
        full: b.finish(&indices.iter().map(|i| full[*i]).collect::<Vec<_>>()),
        indices,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trigonometric_primitives_differentiate_to_every_flux_monomial() {
        // Nine rational circle samples determine the degree-eight numerator
        // after the stereographic substitution, including all four quadrants.
        let samples = ["-3", "-2", "-1", "-1/2", "0", "1/2", "1", "2", "3"].map(|t| {
            let t = Rational::parse(t).unwrap();
            let square = &t * &t;
            let denominator = &Rational::one() + &square;
            [
                &(&Rational::one() - &square) / &denominator,
                &(&Rational::from_i64(2) * &t) / &denominator,
            ]
        });
        for a in 0..=4u8 {
            for b in 0..=4 - a {
                let (mean, terms) = primitive(a, b);
                for [cos, sin] in &samples {
                    let mut derivative = mean.clone();
                    for (c, s, coefficient) in &terms {
                        if *c > 0 {
                            derivative = derivative
                                - (coefficient * &Rational::from_i64(i64::from(*c)))
                                    * cos.pow(i32::from(c - 1)).unwrap()
                                    * sin.pow(i32::from(s + 1)).unwrap();
                        }
                        if *s > 0 {
                            derivative = derivative
                                + (coefficient * &Rational::from_i64(i64::from(*s)))
                                    * cos.pow(i32::from(c + 1)).unwrap()
                                    * sin.pow(i32::from(s - 1)).unwrap();
                        }
                    }
                    assert_eq!(
                        derivative,
                        cos.pow(i32::from(a)).unwrap() * sin.pow(i32::from(b)).unwrap(),
                        "primitive of cos^{a} sin^{b} at ({cos}, {sin})"
                    );
                }
            }
        }
    }
}
