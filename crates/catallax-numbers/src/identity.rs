//! Rational-function identities over shared, opaque real definitions.
//!
//! This proves polynomial identities and rational constants, never algebraic
//! independence. Callers must first establish the entire expression's domain.

use crate::{Rational, Real, decider::Work, real::Expr};
use std::{
    collections::{BTreeMap, HashMap},
    sync::Arc,
};

const LIMIT: &str = "identity-size-limit";
const MAX_TERMS: usize = 128;
const MAX_FACTORS: usize = 64;
const MAX_NODES: usize = 512;
const MAX_WORK: u64 = 10_000;
type Polynomial = BTreeMap<Vec<usize>, Rational>;

// Definition matching is also bounded and iterative. An opaque square-root
// atom can contain a deep radicand even when the expanded polynomial is small.
fn same_definition(a: &Real, b: &Real, work: &mut Work) -> Result<bool, &'static str> {
    let mut pending = vec![(a, b)];
    let mut visited = 0;
    while let Some((a, b)) = pending.pop() {
        work.charge(1)?;
        if visited >= MAX_NODES {
            return Err(LIMIT);
        }
        visited += 1;
        if a.same_node(b) {
            continue;
        }
        if let (Some(a), Some(b)) = (a.rational_ref(), b.rational_ref()) {
            work.exact(a, b, false)?;
            if a != b {
                return Ok(false);
            }
            continue;
        }
        match (&a.node.expr, &b.node.expr) {
            (Expr::Algebraic(a), Expr::Algebraic(b)) => {
                let (
                    crate::algebraic::AlgebraicRepr::Root(a),
                    crate::algebraic::AlgebraicRepr::Root(b),
                ) = (&a.repr, &b.repr)
                else {
                    return Ok(false);
                };
                if Arc::ptr_eq(a, b) {
                    continue;
                }
                let (ac, bc) = (a.polynomial.coefficients(), b.polynomial.coefficients());
                if ac.len() != bc.len() {
                    return Ok(false);
                }
                for (a, b) in ac.iter().zip(bc).chain([(&a.lo, &b.lo), (&a.hi, &b.hi)]) {
                    work.exact(a, b, false)?;
                    if a != b {
                        return Ok(false);
                    }
                }
            }
            (Expr::Neg(a), Expr::Neg(b))
            | (Expr::Sqrt(a), Expr::Sqrt(b))
            | (Expr::Square(a), Expr::Square(b)) => pending.push((a, b)),
            (Expr::Pow(a, n), Expr::Pow(b, m)) if n == m => pending.push((a, b)),
            (Expr::Add(a, b), Expr::Add(c, d))
            | (Expr::Sub(a, b), Expr::Sub(c, d))
            | (Expr::Mul(a, b), Expr::Mul(c, d))
            | (Expr::Div(a, b), Expr::Div(c, d)) => {
                pending.extend([(a, c), (b, d)]);
            }
            _ => return Ok(false),
        }
    }
    Ok(true)
}

#[derive(Clone)]
struct Fraction {
    numerator: Polynomial,
    denominator: Polynomial,
}

fn constant(value: Rational) -> Polynomial {
    if value.is_zero() {
        BTreeMap::new()
    } else {
        BTreeMap::from([(vec![], value)])
    }
}

fn add(
    mut a: Polynomial,
    b: Polynomial,
    negative: bool,
    work: &mut Work,
) -> Result<Polynomial, &'static str> {
    for (monomial, coefficient) in b {
        work.charge(1)?;
        let coefficient = if negative { -coefficient } else { coefficient };
        let previous = a.remove(&monomial).unwrap_or_default();
        work.exact(&previous, &coefficient, false)?;
        let sum = previous + coefficient;
        if !sum.is_zero() {
            a.insert(monomial, sum);
        }
        if a.len() > MAX_TERMS {
            return Err(LIMIT);
        }
    }
    Ok(a)
}

fn multiply(a: &Polynomial, b: &Polynomial, work: &mut Work) -> Result<Polynomial, &'static str> {
    let mut result = BTreeMap::new();
    for (am, ac) in a {
        for (bm, bc) in b {
            work.charge(1)?;
            if am.len() + bm.len() > MAX_FACTORS {
                return Err(LIMIT);
            }
            let mut monomial = am.clone();
            monomial.extend(bm);
            monomial.sort_unstable();
            work.exact(ac, bc, true)?;
            let product = ac * bc;
            let previous = result.remove(&monomial).unwrap_or_default();
            work.exact(&previous, &product, false)?;
            let sum = previous + product;
            if !sum.is_zero() {
                result.insert(monomial, sum);
            }
            // Products with the same monomial share one stored term. Bound
            // actual storage while charging every coefficient operation.
            if result.len() > MAX_TERMS {
                return Err(LIMIT);
            }
        }
    }
    Ok(result)
}

impl Fraction {
    fn constant(value: Rational) -> Self {
        Self {
            numerator: constant(value),
            denominator: constant(Rational::one()),
        }
    }

    fn combine(&self, rhs: &Self, operation: u8, work: &mut Work) -> Result<Self, &'static str> {
        let (numerator, denominator) = match operation {
            0 | 1 => (
                add(
                    multiply(&self.numerator, &rhs.denominator, work)?,
                    multiply(&rhs.numerator, &self.denominator, work)?,
                    operation == 1,
                    work,
                )?,
                multiply(&self.denominator, &rhs.denominator, work)?,
            ),
            2 => (
                multiply(&self.numerator, &rhs.numerator, work)?,
                multiply(&self.denominator, &rhs.denominator, work)?,
            ),
            _ => (
                multiply(&self.numerator, &rhs.denominator, work)?,
                multiply(&self.denominator, &rhs.numerator, work)?,
            ),
        };
        // Definedness is established before entering this proof method. Still
        // decline a formally zero denominator rather than use it as an identity.
        if denominator.is_empty() {
            return Err(LIMIT);
        }
        if numerator.is_empty() {
            return Ok(Self::constant(Rational::zero()));
        }
        Ok(Self {
            numerator,
            denominator,
        })
    }
}

fn ratio(value: &Fraction, work: &mut Work) -> Result<Option<Rational>, &'static str> {
    if value.numerator.is_empty() {
        return Ok(Some(Rational::zero()));
    }
    if value.numerator.len() != value.denominator.len() {
        return Ok(None);
    }
    let (monomial, denominator) = value
        .denominator
        .first_key_value()
        .expect("nonzero formal denominator");
    let Some(numerator) = value.numerator.get(monomial) else {
        return Ok(None);
    };
    work.exact(numerator, denominator, true)?;
    let ratio = numerator / denominator;
    for (monomial, coefficient) in &value.denominator {
        work.exact(&ratio, coefficient, true)?;
        if value.numerator.get(monomial) != Some(&(&ratio * coefficient)) {
            return Ok(None);
        }
    }
    Ok(Some(ratio))
}

// Reduce a whole polynomial in u by u^2=N/D. A single denominator D^k
// replaces the previous product of one denominator per expanded monomial.
fn reduce_root(
    polynomial: &Polynomial,
    atom: usize,
    relation: &Fraction,
    work: &mut Work,
) -> Result<Fraction, &'static str> {
    let count = |m: &[usize]| m.iter().filter(|&&i| i == atom).count();
    let degree = polynomial.keys().map(|m| count(m) / 2).max().unwrap_or(0);
    let powers = |base: &Polynomial, work: &mut Work| -> Result<Vec<Polynomial>, &'static str> {
        let mut values = vec![constant(Rational::one())];
        for i in 1..=degree {
            values.push(multiply(&values[i - 1], base, work)?);
        }
        Ok(values)
    };
    let numerators = powers(&relation.numerator, work)?;
    let denominators = powers(&relation.denominator, work)?;
    let mut numerator = Polynomial::new();
    for (monomial, coefficient) in polynomial {
        work.charge(1)?;
        let exponent = count(monomial);
        let mut reduced: Vec<_> = monomial.iter().copied().filter(|&i| i != atom).collect();
        if exponent % 2 != 0 {
            reduced.push(atom);
            reduced.sort_unstable();
        }
        let term = BTreeMap::from([(reduced, coefficient.clone())]);
        let term = multiply(&term, &numerators[exponent / 2], work)?;
        let term = multiply(&term, &denominators[degree - exponent / 2], work)?;
        numerator = add(numerator, term, false, work)?;
    }
    Ok(Fraction {
        numerator,
        denominator: denominators[degree].clone(),
    })
}

struct Identity<'a, 'w> {
    atoms: Vec<&'a Real>,
    cache: HashMap<*const crate::real::Node, Fraction>,
    work: &'w mut Work,
    rewrites: usize,
    nodes: usize,
}

impl<'a> Identity<'a, '_> {
    fn normalize(&mut self, mut value: Fraction) -> Result<Fraction, &'static str> {
        loop {
            let square = [(false, &value.numerator), (true, &value.denominator)]
                .into_iter()
                .find_map(|(denominator, polynomial)| {
                    polynomial.keys().find_map(|monomial| {
                        monomial.windows(2).find_map(|pair| {
                            (pair[0] == pair[1]
                                && matches!(self.atoms[pair[0]].node.expr, Expr::Sqrt(_)))
                            .then_some((denominator, pair[0]))
                        })
                    })
                });
            let Some((denominator, atom)) = square else {
                return Ok(match ratio(&value, self.work)? {
                    Some(r) => Fraction::constant(r),
                    None => value,
                });
            };
            self.work.charge(1)?;
            if self.rewrites >= MAX_NODES {
                return Err(LIMIT);
            }
            self.rewrites += 1;
            let Expr::Sqrt(radicand) = &self.atoms[atom].node.expr else {
                unreachable!()
            };
            // Replace every even power together, using one common denominator.
            // Definedness proves sqrt(r)^2=r; atom independence is never assumed.
            let relation = self.visit(radicand)?;
            let polynomial = if denominator {
                &value.denominator
            } else {
                &value.numerator
            };
            let replacement = reduce_root(polynomial, atom, &relation, self.work)?;
            value = if denominator {
                Fraction {
                    numerator: value.numerator,
                    denominator: constant(Rational::one()),
                }
                .combine(&replacement, 3, self.work)?
            } else {
                replacement.combine(
                    &Fraction {
                        numerator: value.denominator,
                        denominator: constant(Rational::one()),
                    },
                    3,
                    self.work,
                )?
            };
        }
    }

    fn visit(&mut self, value: &'a Real) -> Result<Fraction, &'static str> {
        self.work.charge(1)?;
        let key = Arc::as_ptr(&value.node);
        if let Some(result) = self.cache.get(&key) {
            return Ok(result.clone());
        }
        if self.nodes >= MAX_NODES {
            return Err(LIMIT);
        }
        // Count before descending: unfinished ancestors are not yet cached.
        // A deep chain must obey the same bound as a wide expression DAG.
        self.nodes += 1;
        let result = if let Some(rational) = value.rational_ref() {
            self.work.exact(rational, rational, false)?;
            Fraction::constant(rational.clone())
        } else {
            match &value.node.expr {
                Expr::Rational(_) => unreachable!(),
                Expr::Algebraic(_) | Expr::Coordinate(..) | Expr::Sqrt(_) => {
                    let mut index = None;
                    for (i, atom) in self.atoms.iter().enumerate() {
                        self.work.charge(1)?;
                        if same_definition(atom, value, self.work)? {
                            index = Some(i);
                            break;
                        }
                    }
                    let index = if let Some(i) = index {
                        i
                    } else {
                        let i = self.atoms.len();
                        if i >= MAX_FACTORS {
                            return Err(LIMIT);
                        }
                        self.atoms.push(value);
                        i
                    };
                    Fraction {
                        numerator: BTreeMap::from([(vec![index], Rational::one())]),
                        denominator: constant(Rational::one()),
                    }
                }
                Expr::Neg(a) => {
                    let mut a = self.visit(a)?;
                    a.numerator = a.numerator.into_iter().map(|(m, c)| (m, -c)).collect();
                    a
                }
                Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                    let a = self.visit(a)?;
                    let b = self.visit(b)?;
                    let operation = match value.node.expr {
                        Expr::Add(..) => 0,
                        Expr::Sub(..) => 1,
                        Expr::Mul(..) => 2,
                        _ => 3,
                    };
                    a.combine(&b, operation, self.work)?
                }
                Expr::Square(a) | Expr::Pow(a, _) => {
                    let mut exponent = match value.node.expr {
                        Expr::Pow(_, n) => n,
                        _ => 2,
                    };
                    let mut factor = self.visit(a)?;
                    let mut result = Fraction::constant(Rational::one());
                    while exponent != 0 {
                        if exponent & 1 != 0 {
                            result = result.combine(&factor, 2, self.work)?;
                        }
                        exponent >>= 1;
                        if exponent != 0 {
                            factor = factor.combine(&factor, 2, self.work)?;
                        }
                    }
                    result
                }
            }
        };
        let result = self.normalize(result)?;
        self.cache.insert(key, result.clone());
        Ok(result)
    }
}

pub(crate) fn eligible(value: &Real) -> bool {
    value.node.identity_atoms == 1 || value.rational_ref().is_some()
}

pub(crate) fn rational(value: &Real, work: &mut Work) -> Result<Option<Rational>, &'static str> {
    // Quadratic towers and pinned system coordinates have their own proof
    // methods. Formal expansion cannot use a system's defining relations; do
    // not spend their allowance on that inconclusive route. Immutable metadata
    // makes this dispatch independent of cache history and expression depth.
    if !eligible(value) {
        return Ok(None);
    }
    let decision_limit = work.policy.max_work;
    let allowance = MAX_WORK.min(decision_limit.saturating_sub(work.used) / 8);
    work.policy.max_work = work.used + allowance;
    let mut identity = Identity {
        atoms: vec![],
        cache: HashMap::new(),
        work,
        rewrites: 0,
        nodes: 0,
    };
    let result = identity
        .visit(value)
        .and_then(|value| ratio(&value, identity.work));
    work.policy.max_work = decision_limit;
    match result {
        Err(LIMIT | "numbers.budget-exhausted") => Ok(None),
        other => other,
    }
}
