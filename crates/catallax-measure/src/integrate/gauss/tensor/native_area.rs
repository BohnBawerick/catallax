//! The same Gauss panel and pure-partial remainder as its parent, evaluated
//! with directed binary64 intervals for the positive area integral.

use super::*;
use catallax_numbers::{Interval, RealPoly};

pub(crate) struct NativeArea {
    square: Polynomial,
    denominator: Polynomial,
}

struct Polynomial {
    // Coefficients by increasing power in both orientations. Keeping both
    // lets the tensor rule reuse the contraction in its inner direction.
    axes: [Vec<Vec<Interval>>; 2],
}

fn point(x: f64) -> Interval {
    Interval::point(x).unwrap()
}

impl Polynomial {
    fn new(poly: &RealPoly) -> Option<Self> {
        let coefficients = poly.rational()?.terms();
        let mut sizes = [1; 2];
        for (power, _) in &coefficients {
            for (size, power) in sizes.iter_mut().zip(power.exponents()) {
                *size = (*size).max(*power as usize + 1);
            }
        }
        let mut axes = [
            vec![vec![point(0.); sizes[1]]; sizes[0]],
            vec![vec![point(0.); sizes[0]]; sizes[1]],
        ];
        for (power, value) in coefficients {
            let [u, v] = <[u32; 2]>::try_from(power.exponents())
                .unwrap()
                .map(|n| n as usize);
            let value = Interval::from_rational(&value);
            if !value.is_finite() {
                return None;
            }
            axes[0][u][v] = value;
            axes[1][v][u] = value;
        }
        Some(Self { axes })
    }

    fn contract(&self, kept: usize, x: Interval) -> Vec<Interval> {
        self.axes[kept].iter().map(|row| horner(row, x)).collect()
    }

    fn degree(&self, axis: usize) -> usize {
        self.axes[axis].len() - 1
    }

    fn contraction_units(&self, axis: usize) -> usize {
        2 * self.axes[axis].len() * self.degree(1 - axis)
    }

    fn jet_units(&self, axis: usize) -> usize {
        (1..=self.degree(axis))
            .map(|i| 2 * (i.min(ORDER) + 1))
            .sum()
    }
}

impl NativeArea {
    pub(crate) fn new(square: &RealPoly, denominator: &RealPoly) -> Option<Self> {
        Some(Self {
            square: Polynomial::new(square)?,
            denominator: Polynomial::new(denominator)?,
        })
    }

    pub(crate) fn integrate(
        &self,
        u: [&Rational; 2],
        v: [&Rational; 2],
        rule: &Rule,
        work: &mut Work,
    ) -> Checked<Vec<Ball>> {
        let (value, errors, dominates) = self.panel(u, v, rule, work)?;
        for (total, error) in work.tensor_error.iter_mut().zip(errors) {
            *total += error;
        }
        let lo = Rational::from_f64_exact(value.lo()).ok_or(Stop::Split)?;
        let hi = Rational::from_f64_exact(value.hi()).ok_or(Stop::Split)?;
        work.quadrature_dominates |= dominates;
        Ok(vec![work.interval(&lo, &hi)?])
    }

    fn panel(
        &self,
        u: [&Rational; 2],
        v: [&Rational; 2],
        rule: &Rule,
        work: &mut Work,
    ) -> Checked<(Interval, [f64; 2], bool)> {
        let axes = [u, v];
        let widths = axes.map(|[lo, hi]| hi - lo);
        let ranges = axes.map(|[lo, hi]| {
            Interval::new(
                Interval::from_rational(lo).lo(),
                Interval::from_rational(hi).hi(),
            )
            .unwrap()
        });
        let mut error = point(0.);
        let mut errors = [0.; 2];
        for axis in 0..2 {
            // Reserve each coefficient operation before evaluating the jet.
            // The three initial operations start the root and quotient series.
            let recurrence_units = 3
                + (1..SIZE)
                    .map(|i| {
                        3 * ((i - 1) / 2)
                            + 2 * usize::from(i % 2 == 0)
                            + 2 * i.min(self.denominator.degree(axis))
                            + 2
                    })
                    .sum::<usize>();
            work.charge_native(
                (self.square.contraction_units(axis)
                    + self.denominator.contraction_units(axis)
                    + self.square.jet_units(axis)
                    + self.denominator.jet_units(axis)
                    + recurrence_units
                    + 2) as u64,
            )?;
            let p = self.square.contract(axis, ranges[1 - axis]);
            let w = self.denominator.contract(axis, ranges[1 - axis]);
            let derivative = root_quotient(jet(&p, ranges[axis]), jet(&w, ranges[axis]))?;
            let correction = derivative * Interval::from_rational(&remainder_factor(&widths, axis));
            errors[axis] = correction.width();
            error = error + correction;
        }
        let count = rule.native_points().len();
        // Map both sets of nodes, then contract the travel direction once
        // per row. Each inner sample has a root, quotient, weight and sum.
        let quadrature_units = 4 * count
            + count
                * (self.square.contraction_units(0)
                    + self.denominator.contraction_units(0)
                    + count * (2 * (self.square.degree(0) + self.denominator.degree(0)) + 4)
                    + 2)
            + 2;
        work.charge_native(quadrature_units as u64)?;
        let nodes = axes.map(|[lo, hi]| {
            let center = Interval::from_rational(&(lo + hi).mul_pow2(-1));
            let half = Interval::from_rational(&(hi - lo).mul_pow2(-1));
            rule.native_points()
                .iter()
                .map(|&(x, w)| (center + x * half, w))
                .collect::<Vec<_>>()
        });
        let mut sum = point(0.);
        for &(v, wv) in &nodes[1] {
            let square = self.square.contract(0, v);
            let denominator = self.denominator.contract(0, v);
            let mut row = point(0.);
            for &(u, wu) in &nodes[0] {
                let density = divide(sqrt(horner(&square, u))?, horner(&denominator, u))?;
                row = row + wu * density;
            }
            sum = sum + wv * row;
        }
        let scale = Interval::from_rational(&(&widths[0] * &widths[1]).mul_pow2(-2));
        let quadrature = sum * scale;
        Ok((
            quadrature + error,
            errors,
            error.width() > quadrature.width(),
        ))
    }
}

fn divide(numerator: Interval, denominator: Interval) -> Checked<Interval> {
    if denominator.contains_zero() {
        return Err(Stop::Split);
    }
    Ok(numerator / denominator)
}

fn sqrt(x: Interval) -> Checked<Interval> {
    if x.lo() <= 0. {
        return Err(Stop::Split);
    }
    x.sqrt().ok_or(Stop::Split)
}

fn horner(coefficients: &[Interval], x: Interval) -> Interval {
    let (&highest, lower) = coefficients.split_last().unwrap();
    lower
        .iter()
        .rev()
        .fold(highest, |sum, &coefficient| sum * x + coefficient)
}

// Horner in x + t gives Taylor coefficients in t, for every x in the range.
fn jet(coefficients: &[Interval], x: Interval) -> [Interval; SIZE] {
    let mut jet = [point(0.); SIZE];
    let (&highest, lower) = coefficients.split_last().unwrap();
    jet[0] = highest;
    for (degree, &coefficient) in lower.iter().rev().enumerate() {
        for j in (1..=(degree + 1).min(ORDER)).rev() {
            jet[j] = jet[j] * x + jet[j - 1];
        }
        jet[0] = jet[0] * x + coefficient;
    }
    jet
}

// Solve root*root = p and quotient*w = root one coefficient at a time.
fn root_quotient(p: [Interval; SIZE], w: [Interval; SIZE]) -> Checked<Interval> {
    let zero = point(0.);
    let mut root = [zero; SIZE];
    root[0] = sqrt(p[0])?;
    let mut quotient = [zero; SIZE];
    quotient[0] = divide(root[0], w[0])?;
    let two_root = root[0] + root[0];
    let w_order = w.iter().rposition(|x| *x != zero).unwrap_or(0);
    for i in 1..SIZE {
        let mut rhs = p[i];
        for j in 1..=(i - 1) / 2 {
            let product = root[j] * root[i - j];
            rhs = rhs - product - product;
        }
        if i % 2 == 0 {
            rhs = rhs - root[i / 2] * root[i / 2];
        }
        root[i] = divide(rhs, two_root)?;
        let mut rhs = root[i];
        for j in 1..=i.min(w_order) {
            rhs = rhs - quotient[i - j] * w[j];
        }
        quotient[i] = divide(rhs, w[0])?;
    }
    Ok(quotient[ORDER])
}
