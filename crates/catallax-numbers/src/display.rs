//! Unproved binary64 reads of immutable definitions. No proof cache is used.

use std::cmp::Ordering;

use crate::{IsolatedRoot, Rational, Real, algebraic::AlgebraicRepr, real::Expr};

impl Real {
    /// Read a display value from this definition without exact arithmetic
    /// after converting stored rational literals to binary64.
    ///
    /// A wide root selector is refined with floating-point arithmetic, rather
    /// than treated as an accurate point. This is not a sign test or enclosure.
    pub fn display_value(&self) -> Option<f64> {
        read(self, 0)
    }
}

fn finite(x: f64) -> Option<f64> {
    x.is_finite().then_some(x)
}

fn read(value: &Real, depth: u32) -> Option<f64> {
    if depth == 512 {
        return None;
    }
    let expression = &value.node.expr;
    match expression {
        Expr::Rational(r) => return finite(r.to_f64_nearest()),
        Expr::Algebraic(a) => {
            if let AlgebraicRepr::Rational(r) = &a.repr {
                return finite(r.to_f64_nearest());
            }
        }
        _ => {}
    }
    let bounds = value.approx();
    let lo = bounds.lo();
    let hi = bounds.hi();
    // The constructor's immutable hardware bounds are often already tight.
    if bounds.is_finite() && hi - lo <= lo.abs().max(hi.abs()).max(1.0) * 8.0 * f64::EPSILON {
        return Some(lo.midpoint(hi));
    }
    let a = |x| read(x, depth + 1);
    finite(match expression {
        Expr::Rational(_) => unreachable!(),
        Expr::Algebraic(a) => match &a.repr {
            AlgebraicRepr::Rational(_) => unreachable!(),
            AlgebraicRepr::Root(r) => root(r.polynomial.coefficients(), [&r.lo, &r.hi])?,
        },
        Expr::Coordinate(r, i) => system_coordinate(r, *i)?,
        Expr::Add(x, y) => a(x)? + a(y)?,
        Expr::Sub(x, y) => a(x)? - a(y)?,
        Expr::Mul(x, y) => a(x)? * a(y)?,
        Expr::Div(x, y) => a(x)? / a(y)?,
        Expr::Neg(x) => -a(x)?,
        Expr::Sqrt(x) => a(x)?.sqrt(),
        Expr::Square(x) => {
            let v = a(x)?;
            v * v
        }
        Expr::Pow(x, n) => power(a(x)?, *n),
    })
}

fn power(mut x: f64, mut n: u32) -> f64 {
    let mut result = 1.0;
    while n != 0 {
        if n & 1 != 0 {
            result *= x;
        }
        n >>= 1;
        if n != 0 {
            x *= x;
        }
    }
    result
}

fn root(coefficients: &[Rational], ends: [&Rational; 2]) -> Option<f64> {
    let [mut lo, mut hi] = ends.map(Rational::to_f64_nearest);
    finite(lo)?;
    finite(hi)?;
    if lo.partial_cmp(&hi) == Some(Ordering::Equal) {
        return Some(lo);
    }
    let c = coefficients
        .iter()
        .map(|x| finite(x.to_f64_nearest()))
        .collect::<Option<Vec<_>>>()?;
    let at = |t| c.iter().rev().fold(0.0, |sum, &v| sum * t + v);
    let start = finite(at(lo))?;
    let end = finite(at(hi))?;
    if start == 0.0 {
        return Some(lo);
    }
    if end == 0.0 {
        return Some(hi);
    }
    if start.is_sign_negative() == end.is_sign_negative() {
        return None;
    }
    for _ in 0..80 {
        let mid = lo.midpoint(hi);
        if [lo, hi]
            .iter()
            .any(|end| mid.partial_cmp(end) == Some(Ordering::Equal))
        {
            return Some(mid);
        }
        let value = finite(at(mid))?;
        if value == 0.0 {
            return Some(mid);
        }
        if value.is_sign_negative() == start.is_sign_negative() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(lo.midpoint(hi))
}

fn system_coordinate(root: &IsolatedRoot, axis: usize) -> Option<f64> {
    let n = root.system().dimension();
    let mut bounds = Vec::with_capacity(n);
    let mut x = Vec::with_capacity(n);
    for (a, b) in root.root_box().lo().iter().zip(root.root_box().hi()) {
        let lo = finite(a.to_f64_nearest())?;
        let hi = finite(b.to_f64_nearest())?;
        bounds.push([lo, hi]);
        x.push(lo.midpoint(hi));
    }
    let rows = root
        .system()
        .equations()
        .iter()
        .map(|p| {
            p.definition_terms()
                .map(|(monomial, c)| Some((monomial.exponents(), finite(c.to_f64_nearest())?)))
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<Vec<_>>>()?;
    let mut matrix = vec![0.0; n.checked_mul(n.checked_add(1)?)?];
    let width = n + 1;
    for _ in 0..32 {
        matrix.fill(0.0);
        for (i, terms) in rows.iter().enumerate() {
            let row = &mut matrix[i * width..(i + 1) * width];
            for &(powers, c) in terms {
                for j in 0..=n {
                    if j < n && powers[j] == 0 {
                        continue;
                    }
                    let mut term = if j == n { -c } else { c * f64::from(powers[j]) };
                    for k in 0..n {
                        term *= power(x[k], powers[k] - u32::from(k == j));
                    }
                    row[j] += term;
                }
            }
        }
        eliminate(&mut matrix, n)?;
        let mut stable = true;
        for i in 0..n {
            let delta = finite(matrix[i * width + n])?;
            let [lo, hi] = bounds[i];
            let next = (x[i] + delta).clamp(lo, hi);
            stable &= delta.abs() <= next.abs().max(1.0) * 32.0 * f64::EPSILON;
            x[i] = next;
        }
        if stable {
            return Some(x[axis]);
        }
    }
    None
}

fn eliminate(matrix: &mut [f64], n: usize) -> Option<()> {
    let width = n + 1;
    for col in 0..n {
        let index = |row, column| row * width + column;
        let pivot = (col..n).max_by(|&a, &b| {
            matrix[index(a, col)]
                .abs()
                .total_cmp(&matrix[index(b, col)].abs())
        })?;
        let divisor = finite(matrix[index(pivot, col)])?;
        if divisor == 0.0 {
            return None;
        }
        for j in col..=n {
            matrix.swap(index(col, j), index(pivot, j));
        }
        for j in col..=n {
            matrix[index(col, j)] /= divisor;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let multiplier = matrix[index(row, col)];
            for j in col..=n {
                matrix[index(row, j)] -= multiplier * matrix[index(col, j)];
            }
        }
    }
    Some(())
}
