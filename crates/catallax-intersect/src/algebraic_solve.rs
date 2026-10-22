//! Certified coefficient-field reduction and finite coordinate projections.

use crate::{math::*, solve::Solutions};
use catallax_numbers::{Decider, Real, RealPoly, Refusal, Sign};
use catallax_shapes::Point3;
use std::collections::{BTreeMap, VecDeque};

#[derive(Clone)]
struct Polynomial(BTreeMap<[u32; 3], Real>);
impl Polynomial {
    fn zero() -> Self {
        Self(BTreeMap::new())
    }
    fn from_real(
        p: &RealPoly,
        permutation: [usize; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut result = Self::zero();
        for (m, c) in p.terms() {
            let e = std::array::from_fn(|i| m.exponents()[permutation[i]]);
            result.add(e, c, dec)?;
        }
        Ok(result)
    }
    fn add(&mut self, e: [u32; 3], c: Real, dec: &mut Decider) -> Result<(), Refusal> {
        let value = self.0.get(&e).map_or_else(|| c.clone(), |v| v + &c);
        let value = if value.rational_ref().is_none() {
            match dec.rational_value(&value)? {
                Some(value) => Real::from(value),
                None => value,
            }
        } else {
            value
        };
        // Coefficients are removed only after a proof, including those produced
        // by cancellation of independently constructed algebraic expressions.
        if sign(&value, dec)? == Sign::Zero {
            self.0.remove(&e);
        } else {
            self.0.insert(e, value);
        }
        Ok(())
    }
    fn leading(&self) -> Option<([u32; 3], &Real)> {
        self.0.last_key_value().map(|(m, c)| (*m, c))
    }
    fn monic(mut self, dec: &mut Decider) -> Result<Self, Refusal> {
        if let Some((leading, coefficient)) = self.leading() {
            let coefficient = coefficient.clone();
            if coefficient
                .rational_ref()
                .is_some_and(|value| value == &catallax_numbers::Rational::one())
            {
                return Ok(self);
            }
            for (monomial, value) in &mut self.0 {
                // add proved the divisor defined and nonzero. Its quotient by
                // itself is one; other quotients retain their exact field proof.
                if *monomial == leading {
                    *value = Real::one();
                    continue;
                }
                let quotient = &*value / &coefficient;
                *value = if quotient.rational_ref().is_none() {
                    match dec.rational_value(&quotient)? {
                        Some(value) => Real::from(value),
                        None => quotient,
                    }
                } else {
                    quotient
                };
            }
        }
        Ok(self)
    }
    fn multiple(
        &self,
        exponents: [u32; 3],
        scalar: &Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut result = Self::zero();
        for (m, c) in &self.0 {
            let e = std::array::from_fn(|i| {
                m[i].checked_add(exponents[i])
                    .expect("monomial exponent overflow")
            });
            result.add(e, c * scalar, dec)?;
        }
        Ok(result)
    }
    fn minus(&self, other: &Self, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut result = self.clone();
        for (m, c) in &other.0 {
            result.add(*m, -c, dec)?;
        }
        Ok(result)
    }
    fn univariate(&self) -> Option<RealPoly> {
        if self.0.keys().any(|m| m[0] != 0 || m[1] != 0) {
            return None;
        }
        Some(RealPoly::from_terms(
            1,
            self.0
                .iter()
                .map(|(m, c)| (vec![m[2]], c.clone()))
                .collect(),
        ))
    }
}
fn divides(a: [u32; 3], b: [u32; 3]) -> bool {
    (0..3).all(|i| a[i] <= b[i])
}
fn reduce(
    mut polynomial: Polynomial,
    basis: &[Polynomial],
    dec: &mut Decider,
) -> Result<Polynomial, Refusal> {
    let mut remainder = Polynomial::zero();
    while let Some((m, c)) = polynomial.leading() {
        dec.check_deadline()?;
        let c = c.clone();
        if let Some((divisor, (n, d))) = basis
            .iter()
            .filter_map(|p| p.leading().map(|l| (p, l)))
            .find(|(_, (n, _))| divides(*n, m))
        {
            let multiple =
                divisor.multiple(std::array::from_fn(|i| m[i] - n[i]), &(&c / d), dec)?;
            polynomial = polynomial.minus(&multiple, dec)?;
        } else {
            polynomial.0.remove(&m);
            remainder.add(m, c, dec)?;
        }
    }
    Ok(remainder)
}
fn groebner(generators: Vec<Polynomial>, dec: &mut Decider) -> Result<Vec<Polynomial>, Refusal> {
    let mut basis = Vec::new();
    let mut pairs = VecDeque::new();
    for p in generators {
        let p = reduce(p, &basis, dec)?.monic(dec)?;
        if !p.0.is_empty() {
            let j = basis.len();
            for i in 0..j {
                pairs.push_back((i, j));
            }
            basis.push(p);
        }
    }
    while let Some((i, j)) = pairs.pop_front() {
        dec.check_deadline()?;
        let a = basis[i].leading().unwrap().0;
        let b = basis[j].leading().unwrap().0;
        if (0..3).all(|k| a[k] == 0 || b[k] == 0) {
            continue;
        }
        let l: [u32; 3] = std::array::from_fn(|k| a[k].max(b[k]));
        let left = basis[i].multiple(std::array::from_fn(|k| l[k] - a[k]), &Real::one(), dec)?;
        let right = basis[j].multiple(std::array::from_fn(|k| l[k] - b[k]), &Real::one(), dec)?;
        let p = reduce(left.minus(&right, dec)?, &basis, dec)?.monic(dec)?;
        if !p.0.is_empty() {
            let j = basis.len();
            for i in 0..j {
                pairs.push_back((i, j));
            }
            basis.push(p);
        }
    }
    Ok(basis)
}

pub(super) fn intersections(
    equations: &[RealPoly],
    surface: &RealPoly,
    dec: &mut Decider,
) -> Result<Solutions, Refusal> {
    let curve = equations
        .iter()
        .map(|p| Polynomial::from_real(p, [0, 1, 2], dec))
        .collect::<Result<Vec<_>, _>>()?;
    let basis = groebner(curve.clone(), dec)?;
    let polynomial = Polynomial::from_real(surface, [0, 1, 2], dec)?;
    if reduce(polynomial, &basis, dec)?.0.is_empty() {
        return Ok(Solutions::Contained);
    }
    let all: Vec<_> = equations.iter().chain(std::iter::once(surface)).collect();
    let mut coordinates: [Vec<Real>; 3] = std::array::from_fn(|_| Vec::new());
    for (axis, output) in coordinates.iter_mut().enumerate() {
        let mut permutation = [0, 1, 2];
        permutation.swap(axis, 2);
        let generators = all
            .iter()
            .map(|p| Polynomial::from_real(p, permutation, dec))
            .collect::<Result<Vec<_>, _>>()?;
        let basis = groebner(generators, dec)?;
        if basis
            .iter()
            .any(|p| p.0.len() == 1 && p.0.contains_key(&[0; 3]))
        {
            return Ok(Solutions::Points(vec![]));
        }
        let mut projection = None;
        for p in basis.iter().filter_map(Polynomial::univariate) {
            if p.total_degree(dec)?.is_some_and(|degree| degree > 0) {
                projection = Some(p);
                break;
            }
        }
        let projection = projection.ok_or_else(|| {
            Refusal::new(
                "numbers.root-not-certified",
                "finite algebraic coordinate projection was not obtained",
            )
        })?;
        *output = projection.real_roots(dec)?;
    }
    // Complete coordinate projections contain every solution, including singular
    // ones. Exact substitution selects correlated tuples, never a Cartesian
    // product of unrelated conjugates certified by a nonsingular-root assertion.
    let mut points = Vec::new();
    for x in &coordinates[0] {
        for y in &coordinates[1] {
            for z in &coordinates[2] {
                let values = [x.clone(), y.clone(), z.clone()];
                let mut incident = true;
                for polynomial in &all {
                    if sign(&polynomial.eval_real(&values), dec)? != Sign::Zero {
                        incident = false;
                        break;
                    }
                }
                if incident {
                    points.push(Point3::new(x.clone(), y.clone(), z.clone()));
                }
            }
        }
    }
    Ok(Solutions::Points(points))
}
