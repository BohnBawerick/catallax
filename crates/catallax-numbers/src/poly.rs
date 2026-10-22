//! Canonical rational polynomials and fraction-free Sylvester determinants.
// Checked and unbounded entry points share each exact algorithm. Checkpoints
// do not charge work; callers keep their existing deterministic reservations.

use crate::{Ball, Integer, Interval, Rational, Real, Sign};
use std::{cmp::Ordering, collections::BTreeMap};

/// Univariate polynomial, lowest degree first, with trailing zeros removed.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct UPoly {
    coefficients: Vec<Rational>,
}
impl UPoly {
    /// Constructs a polynomial, removing trailing zero coefficients.
    pub fn new(mut coefficients: Vec<Rational>) -> Self {
        while coefficients.last().is_some_and(Rational::is_zero) {
            coefficients.pop();
        }
        Self { coefficients }
    }
    /// Zero polynomial.
    pub fn zero() -> Self {
        Self::default()
    }
    /// Constant polynomial.
    pub fn constant(c: Rational) -> Self {
        Self::new(vec![c])
    }
    /// The polynomial `x`.
    pub fn x() -> Self {
        Self::new(vec![Rational::zero(), Rational::one()])
    }
    /// Monic product of the specified linear factors.
    pub fn from_roots(roots: &[Rational]) -> Self {
        roots.iter().fold(Self::constant(Rational::one()), |p, r| {
            &p * &Self::new(vec![-r, Rational::one()])
        })
    }
    /// Canonical coefficients, lowest degree first.
    pub fn coefficients(&self) -> &[Rational] {
        &self.coefficients
    }
    /// Degree, or `None` for zero.
    pub fn degree(&self) -> Option<usize> {
        self.coefficients.len().checked_sub(1)
    }
    /// Whether this is zero.
    pub fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }
    /// Leading coefficient.
    pub fn leading(&self) -> Option<&Rational> {
        self.coefficients.last()
    }
    /// Exact evaluation by Horner's rule.
    pub fn eval(&self, x: &Rational) -> Rational {
        self.coefficients
            .iter()
            .rev()
            .fold(Rational::zero(), |a, c| &(&a * x) + c)
    }
    /// Interval evaluation.
    pub fn eval_interval(&self, x: Interval) -> Interval {
        self.coefficients
            .iter()
            .rev()
            .fold(Interval::point(0.0).unwrap(), |a, c| {
                a * x + Interval::from_rational(c)
            })
    }
    /// Ball evaluation at the argument's working precision.
    pub fn eval_ball(&self, x: &Ball) -> Ball {
        self.coefficients.iter().rev().fold(
            Ball::from_rational(&Rational::zero(), x.precision()),
            |a, c| &(&a * x) + &Ball::from_rational(c, x.precision()),
        )
    }
    /// Exact expression evaluation.
    pub fn eval_real(&self, x: &Real) -> Real {
        let Some((leading, rest)) = self.coefficients.split_last() else {
            return Real::zero();
        };
        rest.iter()
            .rev()
            .fold(Real::from_rational(leading.clone()), |a, c| {
                &(&a * x) + &Real::from_rational(c.clone())
            })
    }
    /// Formal derivative.
    pub fn derivative(&self) -> Self {
        Self::new(
            self.coefficients
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, c)| c * &Rational::from_integer(Integer::from_u64(i as u64)))
                .collect(),
        )
    }
    /// Polynomial division; `None` for a zero divisor.
    pub fn div_rem(&self, divisor: &Self) -> Option<(Self, Self)> {
        self.div_rem_checked(divisor, &|| Ok(()))
            .expect("unbounded polynomial division")
    }
    pub(crate) fn div_rem_checked(
        &self,
        divisor: &Self,
        check: &impl Fn() -> Result<(), &'static str>,
    ) -> Result<Option<(Self, Self)>, &'static str> {
        let Some(d) = divisor.degree() else {
            return Ok(None);
        };
        let mut r = self.clone();
        let mut q = vec![Rational::zero(); self.coefficients.len().saturating_sub(d)];
        while let Some(n) = r.degree().filter(|n| *n >= d) {
            check()?;
            let c = r.leading().unwrap() / divisor.leading().unwrap();
            let k = n - d;
            q[k] = c.clone();
            for i in 0..=d {
                check()?;
                r.coefficients[i + k] = &r.coefficients[i + k] - &(&c * &divisor.coefficients[i]);
            }
            r = Self::new(r.coefficients);
        }
        Ok(Some((Self::new(q), r)))
    }
    /// Monic greatest common divisor.
    pub fn gcd(&self, other: &Self) -> Self {
        self.gcd_checked(other, &|| Ok(()))
            .expect("unbounded polynomial gcd")
    }
    pub(crate) fn gcd_checked(
        &self,
        other: &Self,
        check: &impl Fn() -> Result<(), &'static str>,
    ) -> Result<Self, &'static str> {
        let primitive = |p: Self| -> Result<Self, &'static str> {
            Ok(Self::new(
                p.primitive_integer_checked(check)?
                    .into_iter()
                    .map(Rational::from_integer)
                    .collect(),
            ))
        };
        let mut a = primitive(self.clone())?;
        let mut b = primitive(other.clone())?;
        while let Some(d) = b.degree() {
            let mut r = a;
            while let Some(n) = r.degree().filter(|n| *n >= d) {
                check()?;
                let c = r.leading().unwrap().clone();
                let leading = b.leading().unwrap();
                // A pseudo-division step differs from ordinary division only
                // by a nonzero rational factor. Removing its integer content
                // after each cancellation preserves the GCD over Q while
                // avoiding enormous rational denominators inside a remainder.
                for coefficient in &mut r.coefficients {
                    check()?;
                    *coefficient = &*coefficient * leading;
                }
                for i in 0..=d {
                    check()?;
                    r.coefficients[i + n - d] =
                        &r.coefficients[i + n - d] - &(&c * &b.coefficients[i]);
                }
                r = primitive(Self::new(r.coefficients))?;
            }
            a = b;
            b = r;
        }
        check()?;
        Ok(a.monic())
    }
    /// Monic square-free factors, in increasing multiplicity.
    pub fn square_free_factors(&self) -> Vec<(UPoly, u32)> {
        self.square_free_factors_checked(&|| Ok(()))
            .expect("unbounded square-free factorization")
    }
    pub(crate) fn square_free_factors_checked(
        &self,
        check: &impl Fn() -> Result<(), &'static str>,
    ) -> Result<Vec<(UPoly, u32)>, &'static str> {
        if self.degree().is_none_or(|d| d == 0) {
            return Ok(Vec::new());
        }
        let f = self.monic();
        let mut c = f.gcd_checked(&f.derivative(), check)?;
        let mut w = f.div_rem_checked(&c, check)?.unwrap().0;
        let mut result = Vec::new();
        let mut i = 1u32;
        while w.degree().is_some_and(|d| d > 0) {
            check()?;
            let y = w.gcd_checked(&c, check)?;
            let z = w.div_rem_checked(&y, check)?.unwrap().0;
            if z.degree().is_some_and(|d| d > 0) {
                result.push((z.monic(), i));
            }
            w = y;
            c = c.div_rem_checked(&w, check)?.unwrap().0;
            i += 1;
        }
        Ok(result)
    }
    /// Polynomial composition.
    pub fn compose(&self, inner: &Self) -> Self {
        self.coefficients.iter().rev().fold(Self::zero(), |a, c| {
            &(&a * inner) + &Self::constant(c.clone())
        })
    }
    /// Coprime integer coefficients with positive leading coefficient.
    pub fn primitive_integer(&self) -> Vec<Integer> {
        self.primitive_integer_checked(&|| Ok(()))
            .expect("unbounded primitive integer coefficients")
    }
    pub(crate) fn primitive_integer_checked(
        &self,
        check: &impl Fn() -> Result<(), &'static str>,
    ) -> Result<Vec<Integer>, &'static str> {
        let mut lcm = Integer::one();
        for c in &self.coefficients {
            check()?;
            lcm = (&lcm * c.denominator())
                .div_rem_floor(&lcm.gcd(c.denominator()))
                .unwrap()
                .0;
        }
        let mut result: Vec<_> = self
            .coefficients
            .iter()
            .map(|c| {
                check()?;
                Ok((&c.n * &lcm).div_rem_floor(&c.d).unwrap().0)
            })
            .collect::<Result<_, &'static str>>()?;
        let mut g = Integer::zero();
        for n in &result {
            check()?;
            g = g.gcd(n);
        }
        if g.is_zero() {
            return Ok(Vec::new());
        }
        if result.last().unwrap().sign() == Sign::Negative {
            g = -g;
        }
        for n in &mut result {
            check()?;
            *n = n.div_rem_floor(&g).unwrap().0;
        }
        Ok(result)
    }
    pub(crate) fn monic(&self) -> Self {
        match self.leading() {
            None => Self::zero(),
            Some(c) => self.scale(&c.recip().unwrap()),
        }
    }
    pub(crate) fn scale(&self, c: &Rational) -> Self {
        Self::new(self.coefficients.iter().map(|a| a * c).collect())
    }
    fn add(&self, rhs: &Self) -> Self {
        let mut c = vec![Rational::zero(); self.coefficients.len().max(rhs.coefficients.len())];
        for (i, a) in self.coefficients.iter().enumerate() {
            c[i] = a.clone();
        }
        for (i, a) in rhs.coefficients.iter().enumerate() {
            c[i] = &c[i] + a;
        }
        Self::new(c)
    }
    fn sub(&self, rhs: &Self) -> Self {
        self + &(-rhs)
    }
    fn mul(&self, rhs: &Self) -> Self {
        if self.is_zero() || rhs.is_zero() {
            return Self::zero();
        }
        let mut c = vec![Rational::zero(); self.coefficients.len() + rhs.coefficients.len() - 1];
        for (i, a) in self.coefficients.iter().enumerate() {
            for (j, b) in rhs.coefficients.iter().enumerate() {
                c[i + j] = &c[i + j] + &(a * b);
            }
        }
        Self::new(c)
    }
    fn neg(&self) -> Self {
        self.scale(&Rational::from_i64(-1))
    }
}

/// A monomial, ordered by total degree and then lexicographic exponents.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Monomial {
    exponents: Vec<u32>,
}
impl Monomial {
    /// Constructs a monomial.
    pub fn new(exponents: Vec<u32>) -> Self {
        Self { exponents }
    }
    /// Exponent vector.
    pub fn exponents(&self) -> &[u32] {
        &self.exponents
    }
    /// Total degree.
    pub fn degree(&self) -> u32 {
        self.exponents
            .iter()
            .try_fold(0u32, |a, b| a.checked_add(*b))
            .expect("monomial degree overflow")
    }
}
impl Ord for Monomial {
    fn cmp(&self, other: &Self) -> Ordering {
        self.degree()
            .cmp(&other.degree())
            .then_with(|| self.exponents.cmp(&other.exponents))
    }
}
impl PartialOrd for Monomial {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Canonical multivariate polynomial in a fixed number of variables.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct MPoly {
    nvars: usize,
    terms: BTreeMap<Monomial, Rational>,
}
impl MPoly {
    /// Zero polynomial.
    pub fn zero(nvars: usize) -> Self {
        Self {
            nvars,
            terms: BTreeMap::new(),
        }
    }
    /// Constant polynomial.
    pub fn constant(nvars: usize, c: Rational) -> Self {
        Self::from_terms(nvars, vec![(vec![0; nvars], c)])
    }
    /// A variable; panics for an out-of-range index.
    pub fn var(nvars: usize, index: usize) -> Self {
        assert!(index < nvars);
        let mut e = vec![0; nvars];
        e[index] = 1;
        Self::from_terms(nvars, vec![(e, Rational::one())])
    }
    /// Sums like terms and removes zero coefficients.
    pub fn from_terms(nvars: usize, terms: Vec<(Vec<u32>, Rational)>) -> Self {
        let mut p = Self::zero(nvars);
        for (e, c) in terms {
            assert_eq!(e.len(), nvars);
            p.accumulate(Monomial::new(e), c);
        }
        p
    }
    /// Number of variables.
    pub fn nvars(&self) -> usize {
        self.nvars
    }
    /// Whether the polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
    /// Total degree, or `None` for zero.
    pub fn total_degree(&self) -> Option<u32> {
        self.terms.last_key_value().map(|(m, _)| m.degree())
    }
    /// Terms in graded lexicographic order.
    pub fn terms(&self) -> Vec<(Monomial, Rational)> {
        self.terms
            .iter()
            .map(|(m, c)| (m.clone(), c.clone()))
            .collect()
    }
    /// Scalar multiplication.
    pub fn scale(&self, c: &Rational) -> Self {
        Self::from_terms(
            self.nvars,
            self.terms
                .iter()
                .map(|(m, a)| (m.exponents.clone(), a * c))
                .collect(),
        )
    }
    /// Exact rational evaluation.
    pub fn eval(&self, x: &[Rational]) -> Rational {
        assert_eq!(x.len(), self.nvars);
        self.terms.iter().fold(Rational::zero(), |a, (m, c)| {
            &a + &m
                .exponents
                .iter()
                .zip(x)
                .filter(|(e, _)| **e != 0)
                .fold(c.clone(), |b, (e, v)| &b * &pow_rational(v, *e))
        })
    }
    /// Interval evaluation on a box.
    pub fn eval_interval(&self, x: &[Interval]) -> Interval {
        assert_eq!(x.len(), self.nvars);
        self.terms
            .iter()
            .fold(Interval::point(0.0).unwrap(), |a, (m, c)| {
                a + m
                    .exponents
                    .iter()
                    .zip(x)
                    .filter(|(e, _)| **e != 0)
                    .fold(Interval::from_rational(c), |b, (e, v)| {
                        b * pow_interval(*v, *e)
                    })
            })
    }
    /// Ball evaluation on a box.
    pub fn eval_ball(&self, x: &[Ball]) -> Ball {
        assert_eq!(x.len(), self.nvars);
        let p = x.iter().map(Ball::precision).max().unwrap_or(1);
        self.terms
            .iter()
            .fold(Ball::from_rational(&Rational::zero(), p), |a, (m, c)| {
                &a + &m
                    .exponents
                    .iter()
                    .zip(x)
                    .filter(|(e, _)| **e != 0)
                    .fold(Ball::from_rational(c, p), |b, (e, v)| &b * &v.pow_u32(*e))
            })
    }
    /// Exact real expression evaluation.
    pub fn eval_real(&self, x: &[Real]) -> Real {
        assert_eq!(x.len(), self.nvars);
        self.terms.iter().fold(Real::zero(), |a, (m, c)| {
            &a + &m
                .exponents
                .iter()
                .zip(x)
                .filter(|(e, _)| **e != 0)
                .fold(Real::from_rational(c.clone()), |b, (e, v)| {
                    &b * &v.pow_u32(*e)
                })
        })
    }
    /// Partial derivative.
    pub fn partial(&self, var: usize) -> Self {
        assert!(var < self.nvars);
        Self::from_terms(
            self.nvars,
            self.terms
                .iter()
                .filter(|(m, _)| m.exponents[var] != 0)
                .map(|(m, c)| {
                    let mut e = m.exponents.clone();
                    let n = e[var];
                    e[var] -= 1;
                    (
                        e,
                        c * &Rational::from_integer(Integer::from_u64(u64::from(n))),
                    )
                })
                .collect(),
        )
    }
    /// Substitution of a polynomial for one variable.
    pub fn substitute(&self, var: usize, value: &MPoly) -> Self {
        assert!(var < self.nvars);
        assert_eq!(self.nvars, value.nvars);
        self.terms.iter().fold(Self::zero(self.nvars), |a, (m, c)| {
            let mut e = m.exponents.clone();
            let n = e[var];
            e[var] = 0;
            &a + &(&Self::from_terms(self.nvars, vec![(e, c.clone())]) * &value.power(n))
        })
    }
    /// Sylvester resultant, retaining the original variable count.
    pub fn resultant(&self, other: &Self, var: usize) -> Self {
        self.resultant_checked(other, var, &|| Ok(()))
            .expect("unbounded polynomial resultant")
    }
    pub(crate) fn resultant_checked(
        &self,
        other: &Self,
        var: usize,
        check: &impl Fn() -> Result<(), &'static str>,
    ) -> Result<Self, &'static str> {
        check()?;
        assert!(var < self.nvars);
        assert_eq!(self.nvars, other.nvars);
        if self.is_zero() || other.is_zero() {
            return Ok(Self::zero(self.nvars));
        }
        let a = self.coefficients_in(var);
        let b = other.coefficients_in(var);
        let m = a.len() - 1;
        let n = b.len() - 1;
        if m == 0 {
            return a[0].power_checked(n.try_into().expect("polynomial degree fits u32"), check);
        }
        if n == 0 {
            return b[0].power_checked(m.try_into().expect("polynomial degree fits u32"), check);
        }
        let size = m + n;
        let zero = Self::zero(self.nvars);
        let mut matrix = vec![vec![zero.clone(); size]; size];
        for i in 0..n {
            check()?;
            for j in 0..=m {
                matrix[i][i + j] = a[m - j].clone();
            }
        }
        for i in 0..m {
            check()?;
            for j in 0..=n {
                matrix[n + i][i + j] = b[n - j].clone();
            }
        }
        let mut prev = Self::constant(self.nvars, Rational::one());
        let mut negative = false;
        for k in 0..size - 1 {
            check()?;
            let Some(pivot_row) = (k..size).find(|i| !matrix[*i][k].is_zero()) else {
                return Ok(zero);
            };
            if pivot_row != k {
                matrix.swap(k, pivot_row);
                negative = !negative;
            }
            let pivot = matrix[k][k].clone();
            for i in k + 1..size {
                for j in k + 1..size {
                    check()?;
                    let numerator = &pivot.mul_checked(&matrix[i][j], check)?
                        - &matrix[i][k].mul_checked(&matrix[k][j], check)?;
                    matrix[i][j] = numerator
                        .exact_div_checked(&prev, check)?
                        .expect("Bareiss exact division");
                }
            }
            for row in matrix.iter_mut().skip(k + 1) {
                row[k] = zero.clone();
            }
            prev = pivot;
        }
        let det = matrix[size - 1][size - 1].clone();
        Ok(if negative { -det } else { det })
    }
    /// Converts to a univariate polynomial if all other variables are absent.
    pub fn to_upoly(&self, var: usize) -> Option<UPoly> {
        assert!(var < self.nvars);
        if self.terms.keys().any(|m| {
            m.exponents
                .iter()
                .enumerate()
                .any(|(i, e)| i != var && *e != 0)
        }) {
            return None;
        }
        let degree = self
            .terms
            .keys()
            .map(|m| m.exponents[var] as usize)
            .max()
            .unwrap_or(0);
        let mut c = vec![Rational::zero(); degree + 1];
        for (m, a) in &self.terms {
            c[m.exponents[var] as usize] = a.clone();
        }
        Some(UPoly::new(c))
    }
    pub(crate) fn coefficients_in(&self, var: usize) -> Vec<Self> {
        let degree = self
            .terms
            .keys()
            .map(|m| m.exponents[var] as usize)
            .max()
            .unwrap_or(0);
        let mut c = vec![Self::zero(self.nvars); degree + 1];
        for (m, a) in &self.terms {
            let mut e = m.exponents.clone();
            let i = e[var] as usize;
            e[var] = 0;
            c[i].accumulate(Monomial::new(e), a.clone());
        }
        c
    }
    pub(crate) fn power(&self, n: u32) -> Self {
        self.power_checked(n, &|| Ok(()))
            .expect("unbounded polynomial power")
    }
    fn power_checked(
        &self,
        mut n: u32,
        check: &impl Fn() -> Result<(), &'static str>,
    ) -> Result<Self, &'static str> {
        let mut a = self.clone();
        let mut r = Self::constant(self.nvars, Rational::one());
        while n != 0 {
            check()?;
            if n & 1 != 0 {
                r = r.mul_checked(&a, check)?;
            }
            n >>= 1;
            if n != 0 {
                a = a.mul_checked(&a, check)?;
            }
        }
        Ok(r)
    }
    fn exact_div_checked(
        &self,
        rhs: &Self,
        check: &impl Fn() -> Result<(), &'static str>,
    ) -> Result<Option<Self>, &'static str> {
        if rhs.is_zero() {
            return Ok(None);
        }
        let mut r = self.clone();
        let mut q = Self::zero(self.nvars);
        let (m, c) = rhs.terms.last_key_value().unwrap();
        while let Some((n, d)) = r.terms.last_key_value() {
            check()?;
            let e: Option<Vec<_>> = n
                .exponents
                .iter()
                .zip(&m.exponents)
                .map(|(a, b)| a.checked_sub(*b))
                .collect();
            let Some(e) = e else {
                return Ok(None);
            };
            let term = Self::from_terms(self.nvars, vec![(e, d / c)]);
            q = &q + &term;
            r = &r - &term.mul_checked(rhs, check)?;
        }
        Ok(Some(q))
    }
    pub(crate) fn construction_cmp(&self, other: &Self) -> Ordering {
        self.nvars
            .cmp(&other.nvars)
            .then_with(|| self.terms.cmp(&other.terms))
    }
    pub(crate) fn term_count(&self) -> usize {
        self.terms.len()
    }
    fn accumulate(&mut self, m: Monomial, c: Rational) {
        if c.is_zero() {
            return;
        }
        let value = self.terms.get(&m).map_or(c.clone(), |a| a + &c);
        if value.is_zero() {
            self.terms.remove(&m);
        } else {
            self.terms.insert(m, value);
        }
    }
    fn add(&self, rhs: &Self) -> Self {
        assert_eq!(self.nvars, rhs.nvars);
        let mut p = self.clone();
        for (m, c) in &rhs.terms {
            p.accumulate(m.clone(), c.clone());
        }
        p
    }
    fn sub(&self, rhs: &Self) -> Self {
        self + &(-rhs)
    }
    fn mul(&self, rhs: &Self) -> Self {
        self.mul_checked(rhs, &|| Ok(()))
            .expect("unbounded polynomial product")
    }
    fn mul_checked(
        &self,
        rhs: &Self,
        check: &impl Fn() -> Result<(), &'static str>,
    ) -> Result<Self, &'static str> {
        assert_eq!(self.nvars, rhs.nvars);
        let mut p = Self::zero(self.nvars);
        for (m, c) in &self.terms {
            for (n, d) in &rhs.terms {
                check()?;
                let e = m
                    .exponents
                    .iter()
                    .zip(&n.exponents)
                    .map(|(a, b)| a.checked_add(*b).expect("monomial exponent overflow"))
                    .collect();
                p.accumulate(Monomial::new(e), c * d);
            }
        }
        Ok(p)
    }
    fn neg(&self) -> Self {
        self.scale(&Rational::from_i64(-1))
    }
}
fn pow_rational(x: &Rational, mut e: u32) -> Rational {
    let mut a = x.clone();
    let mut r = Rational::one();
    while e != 0 {
        if e & 1 != 0 {
            r = &r * &a;
        }
        e >>= 1;
        if e != 0 {
            a = &a * &a;
        }
    }
    r
}
fn pow_interval(x: Interval, mut e: u32) -> Interval {
    let mut a = x;
    let mut r = Interval::point(1.0).unwrap();
    while e != 0 {
        if e & 1 != 0 {
            r = r * a;
        }
        e >>= 1;
        if e != 0 {
            a = a.square();
        }
    }
    r
}
binops!(UPoly; Add::add, Sub::sub, Mul::mul);
negop!(UPoly);
binops!(MPoly; Add::add, Sub::sub, Mul::mul);
negop!(MPoly);
