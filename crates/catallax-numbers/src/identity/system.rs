//! Exact consequences of equations retained by a certified root constructor.

use super::*;
use std::collections::HashSet;

fn relations(value: &Real, work: &mut Work) -> Result<Vec<Real>, &'static str> {
    let mut pending = vec![value];
    let mut visited = HashSet::new();
    let mut roots: Vec<(&Arc<crate::IsolatedRoot>, &Arc<[crate::RealPoly]>)> = Vec::new();
    while let Some(value) = pending.pop() {
        if value.node.identity_atoms & 4 == 0 || !visited.insert(Arc::as_ptr(&value.node)) {
            continue;
        }
        work.charge(1)?;
        if visited.len() > MAX_NODES {
            return Err(LIMIT);
        }
        match &value.node.expr {
            Expr::Coordinate(root, _) => {
                if let Some(equations) = &root.original_equations
                    && !roots.iter().any(|(_, old)| Arc::ptr_eq(old, equations))
                {
                    if roots.len() >= MAX_FACTORS {
                        return Err(LIMIT);
                    }
                    roots.push((root, equations));
                }
            }
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                pending.extend([a, b]);
            }
            Expr::Neg(a) | Expr::Sqrt(a) | Expr::Square(a) | Expr::Pow(a, _) => pending.push(a),
            Expr::Rational(_) | Expr::Algebraic(_) => {}
        }
    }
    let mut result = Vec::new();
    for (root, equations) in roots {
        let coordinates: Vec<_> = (0..equations.len())
            .map(|axis| Real::from_shared_root_coordinate(root.clone(), axis))
            .collect();
        for equation in equations.iter() {
            work.charge(1)?;
            result.push(equation.eval_real(&coordinates));
        }
    }
    Ok(result)
}

fn leading(polynomial: &Polynomial) -> Option<(&Vec<usize>, &Rational)> {
    // Graded lexicographic order is preserved by multiplication. The map's
    // sorted-factor-list order alone is not a monomial order.
    polynomial
        .iter()
        .max_by(|(a, _), (b, _)| a.len().cmp(&b.len()).then_with(|| a.cmp(b)))
}

fn quotient(product: &[usize], factor: &[usize]) -> Option<Vec<usize>> {
    let mut index = 0;
    for atom in product {
        if factor.get(index) == Some(atom) {
            index += 1;
        }
    }
    if index != factor.len() {
        return None;
    }
    let mut result = Vec::with_capacity(product.len() - factor.len());
    index = 0;
    for atom in product {
        if factor.get(index) == Some(atom) {
            index += 1;
        } else {
            result.push(*atom);
        }
    }
    Some(result)
}

fn reduce(
    mut polynomial: Polynomial,
    equations: &[Polynomial],
    work: &mut Work,
) -> Result<(Polynomial, bool), &'static str> {
    let mut remainder = Polynomial::new();
    let mut changed = false;
    while let Some((monomial, coefficient)) = leading(&polynomial) {
        work.charge(1)?;
        let mut divisor = None;
        for equation in equations {
            work.charge(1)?;
            let Some((power, scale)) = leading(equation) else {
                continue;
            };
            if let Some(power) = quotient(monomial, power) {
                divisor = Some((equation, power, scale));
                break;
            }
        }
        if let Some((equation, power, scale)) = divisor {
            work.exact(coefficient, scale, true)?;
            let factor = Polynomial::from([(power, coefficient / scale)]);
            let product = multiply(&factor, equation, work)?;
            polynomial = add(polynomial, product, true, work)?;
            changed = true;
        } else {
            let key = monomial.clone();
            let (power, coefficient) = polynomial.remove_entry(&key).expect("leading term");
            remainder.insert(power, coefficient);
            if remainder.len() > MAX_TERMS {
                return Err(LIMIT);
            }
        }
    }
    Ok((remainder, changed))
}

pub(super) fn rational(value: &Real, work: &mut Work) -> Result<Option<Rational>, &'static str> {
    // The caller establishes the complete expression's domain before using
    // equations. In particular, a vanishing numerator cannot hide division by zero.
    let decision_limit = work.policy.max_work;
    let allowance = MAX_WORK.min(decision_limit.saturating_sub(work.used) / 8);
    work.policy.max_work = work.used + allowance;
    let result = (|| {
        let equations = relations(value, work)?;
        let mut identity = Identity {
            atoms: vec![],
            cache: HashMap::new(),
            work,
            rewrites: 0,
            nodes: 0,
        };
        let mut value = identity.visit(value)?;
        let mut basis = Vec::with_capacity(equations.len());
        for equation in &equations {
            let equation = identity.visit(equation)?;
            let (numerator, _) = reduce(equation.numerator, &basis, identity.work)?;
            let equation = identity.normalize(Fraction {
                numerator,
                denominator: constant(Rational::one()),
            })?;
            if !equation.numerator.is_empty() {
                // A nonzero constant cannot be a defining zero at a certified root.
                if equation.numerator.keys().all(Vec::is_empty) {
                    return Err(LIMIT);
                }
                basis.push(equation.numerator);
            }
        }
        for _ in 0..MAX_NODES {
            let (numerator, a) = reduce(value.numerator, &basis, identity.work)?;
            let (denominator, b) = reduce(value.denominator, &basis, identity.work)?;
            if denominator.is_empty() {
                return Ok(None);
            }
            let before = identity.rewrites;
            value = identity.normalize(Fraction {
                numerator,
                denominator,
            })?;
            if let Some(result) = ratio(&value, identity.work)? {
                return Ok(Some(result));
            }
            if !a && !b && before == identity.rewrites {
                return Ok(None);
            }
        }
        Ok(None)
    })();
    work.policy.max_work = decision_limit;
    match result {
        Err(LIMIT | "numbers.budget-exhausted") => Ok(None),
        other => other,
    }
}
