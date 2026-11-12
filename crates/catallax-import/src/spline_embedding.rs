//! A common positive projection of the exact rational hodograph proves embedding.

use catallax_numbers::{Decider, Rational, Sign};
use catallax_shapes::spline_tools::RationalSplineCurve;

use crate::Vec3;
use crate::incidence_arithmetic::{
    arithmetic as op, binomial, checkpoint, compare, dot, equal, sub,
};

type Point = [Rational; 4];
type Net = Vec<Point>;
const CELLS: usize = 4096;
const DEPTH: usize = 24;

// Blossoming on one nonzero knot span avoids global power-basis cancellation.
pub(crate) fn bezier(
    source: &RationalSplineCurve,
    span: usize,
    dec: &mut Decider,
) -> Result<Net, String> {
    let degree = source.degree();
    let knots = source.knots();
    let mut input = Vec::new();
    for i in span - degree..=span {
        checkpoint(dec)?;
        let w = &source.weights()[i];
        let p = &source.controls()[i];
        input.push([
            op(&p[0], w, '*', dec)?,
            op(&p[1], w, '*', dec)?,
            op(&p[2], w, '*', dec)?,
            w.clone(),
        ]);
    }
    let mut out = Vec::new();
    for j in 0..=degree {
        checkpoint(dec)?;
        let mut work = input.clone();
        for r in 1..=degree {
            let t = &knots[span + usize::from(r > degree - j)];
            for i in (r..=degree).rev() {
                checkpoint(dec)?;
                let index = span - degree + i;
                let width = op(&knots[index + degree - r + 1], &knots[index], '-', dec)?;
                if width.sign() != Sign::Positive {
                    return Err("embedding blossom has a nonpositive knot interval".into());
                }
                let delta = op(t, &knots[index], '-', dec)?;
                let alpha = op(&delta, &width, '/', dec)?;
                if alpha.sign() == Sign::Negative || compare(&alpha, &Rational::one(), dec)?.is_gt()
                {
                    return Err("embedding blossom is outside its knot span".into());
                }
                let beta = op(&Rational::one(), &alpha, '-', dec)?;
                let (previous, next) = work.split_at_mut(i);
                for (value, previous) in next[0].iter_mut().zip(&previous[i - 1]) {
                    let a = op(previous, &beta, '*', dec)?;
                    let b = op(value, &alpha, '*', dec)?;
                    *value = op(&a, &b, '+', dec)?;
                }
            }
        }
        out.push(work[degree].clone());
    }
    Ok(out)
}

fn point(p: &Point, dec: &mut Decider) -> Result<Vec3, String> {
    if p[3].sign() != Sign::Positive {
        return Err("embedding needs positive homogeneous weights".into());
    }
    Ok([
        op(&p[0], &p[3], '/', dec)?,
        op(&p[1], &p[3], '/', dec)?,
        op(&p[2], &p[3], '/', dec)?,
    ])
}
fn joined(a: &Point, b: &Point, dec: &mut Decider) -> Result<bool, String> {
    for i in 0..3 {
        let x = op(&a[i], &b[3], '*', dec)?;
        let y = op(&b[i], &a[3], '*', dec)?;
        if !equal(&x, &y, dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}
fn product(a: &[Rational], b: &[Rational], dec: &mut Decider) -> Result<Vec<Rational>, String> {
    let (m, n) = (a.len() - 1, b.len() - 1);
    let mut out = vec![Rational::zero(); m + n + 1];
    for (i, a) in a.iter().enumerate() {
        for (j, b) in b.iter().enumerate() {
            checkpoint(dec)?;
            let factor = op(&binomial(m, i, dec)?, &binomial(n, j, dec)?, '*', dec)?;
            let factor = op(&factor, &binomial(m + n, i + j, dec)?, '/', dec)?;
            let value = op(a, b, '*', dec)?;
            let value = op(&value, &factor, '*', dec)?;
            out[i + j] = op(&out[i + j], &value, '+', dec)?;
        }
    }
    Ok(out)
}
fn derivative(p: &[Rational], dec: &mut Decider) -> Result<Vec<Rational>, String> {
    let degree = Rational::from_i64((p.len() - 1) as i64);
    p.windows(2)
        .map(|p| {
            let d = op(&p[1], &p[0], '-', dec)?;
            op(&d, &degree, '*', dec)
        })
        .collect()
}
fn hodograph(
    net: &Net,
    origin: &Vec3,
    direction: &Vec3,
    dec: &mut Decider,
) -> Result<Vec<Rational>, String> {
    let mut scalar = Vec::new();
    let mut weights = Vec::new();
    for p in net {
        checkpoint(dec)?;
        if p[3].sign() != Sign::Positive {
            return Err("embedding needs positive Bernstein weights".into());
        }
        let mut relative = std::array::from_fn(|_| Rational::zero());
        for i in 0..3 {
            let offset = op(&origin[i], &p[3], '*', dec)?;
            relative[i] = op(&p[i], &offset, '-', dec)?;
        }
        scalar.push(dot(&relative, direction, dec)?);
        weights.push(p[3].clone());
    }
    let a = product(&derivative(&scalar, dec)?, &weights, dec)?;
    let b = product(&scalar, &derivative(&weights, dec)?, dec)?;
    a.iter().zip(&b).map(|(a, b)| op(a, b, '-', dec)).collect()
}
fn split(p: &[Rational], dec: &mut Decider) -> Result<(Vec<Rational>, Vec<Rational>), String> {
    let mut work = p.to_vec();
    let mut left = vec![work[0].clone()];
    let mut right = vec![work[work.len() - 1].clone()];
    let half = Rational::from_ratio(1, 2).unwrap();
    while work.len() > 1 {
        let mut next = Vec::new();
        for pair in work.windows(2) {
            let sum = op(&pair[0], &pair[1], '+', dec)?;
            next.push(op(&sum, &half, '*', dec)?);
        }
        work = next;
        left.push(work[0].clone());
        right.push(work[work.len() - 1].clone());
    }
    right.reverse();
    Ok((left, right))
}
fn positive(p: Vec<Rational>, cells: &mut usize, dec: &mut Decider) -> Result<(), String> {
    let mut stack = vec![(p, 0)];
    while let Some((p, depth)) = stack.pop() {
        checkpoint(dec)?;
        if *cells >= CELLS {
            return Err("embedding cone exceeded its 4096-cell certificate allowance".into());
        }
        *cells += 1;
        if p.iter().all(|x| x.sign() == Sign::Positive) {
            continue;
        }
        if p[0].sign() != Sign::Positive || p[p.len() - 1].sign() != Sign::Positive {
            return Err("embedding projection has a nonpositive derivative witness".into());
        }
        if depth == DEPTH {
            return Err("embedding cone exceeded its depth-24 certificate allowance".into());
        }
        let (left, right) = split(&p, dec)?;
        stack.push((right, depth + 1));
        stack.push((left, depth + 1));
    }
    Ok(())
}

pub(crate) fn prove(source: &RationalSplineCurve, dec: &mut Decider) -> Result<(), String> {
    checkpoint(dec)?;
    let mut nets = Vec::new();
    for span in source.degree()..source.controls().len() {
        checkpoint(dec)?;
        if compare(&source.knots()[span], &source.knots()[span + 1], dec)?.is_lt() {
            nets.push(bezier(source, span, dec)?);
        }
    }
    let first = nets.first().ok_or("empty embedding span set")?;
    let last = nets.last().unwrap();
    let origin = point(&first[0], dec)?;
    let end = point(last.last().unwrap(), dec)?;
    let direction = sub(&end, &origin, dec)?;
    if direction.iter().all(Rational::is_zero) {
        return Err("embedding endpoint projection is zero".into());
    }
    let mut cells = 0;
    for (i, net) in nets.iter().enumerate() {
        checkpoint(dec)?;
        if i > 0 && !joined(nets[i - 1].last().unwrap(), &net[0], dec)? {
            return Err("embedding spans have no exact C0 join".into());
        }
        positive(hodograph(net, &origin, &direction, dec)?, &mut cells, dec)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;
    fn quarter() -> RationalSplineCurve {
        let controls = [
            ["1", "0", "0"],
            ["1", "0.55228475", "0"],
            ["0.55228475", "1", "0"],
            ["0", "1", "0"],
        ]
        .map(|p| p.map(|v| Rational::parse(v).unwrap()))
        .to_vec();
        RationalSplineCurve::new(
            3,
            [0, 0, 0, 0, 1, 1, 1, 1].map(Rational::from_i64).to_vec(),
            controls,
            vec![Rational::one(); 4],
            &Decider::default(),
        )
        .unwrap()
    }
    #[test]
    fn quarter_circle_cone_is_exact_and_charged() {
        let source = quarter();
        let (result, work) =
            with_operation_budget(20_000_000, || prove(&source, &mut Decider::default()));
        assert!(result.is_ok(), "{result:?}");
        assert!(work.used > 0);
        for limit in [0, work.used / 2, work.used - 1] {
            let (result, work) =
                with_operation_budget(limit, || prove(&source, &mut Decider::default()));
            assert!(result.is_err());
            assert!(work.exhausted);
        }
        let shift = Rational::parse(&format!("1{}", "0".repeat(300))).unwrap();
        let controls = source
            .controls()
            .iter()
            .map(|p| p.clone().map(|x| &x + &shift))
            .collect();
        let larger = RationalSplineCurve::new(
            source.degree(),
            source.knots().to_vec(),
            controls,
            source.weights().to_vec(),
            &Decider::default(),
        )
        .unwrap();
        let (result, larger_work) =
            with_operation_budget(20_000_000, || prove(&larger, &mut Decider::default()));
        assert!(result.is_ok(), "{result:?}");
        assert!(larger_work.used > work.used);
    }
    #[test]
    fn exact_subdivision_proves_a_strictly_positive_hodograph_with_a_negative_coefficient() {
        let p = ["0.26", "-0.24", "0.26"]
            .map(|v| Rational::parse(v).unwrap())
            .to_vec();
        let mut cells = 0;
        positive(p, &mut cells, &mut Decider::default()).unwrap();
        assert_eq!(cells, 3);
        let p = [1, -1, 1].map(Rational::from_i64).to_vec();
        assert!(positive(p, &mut 0, &mut Decider::default()).is_err());
    }
    #[test]
    fn discontinuous_spans_cannot_share_an_embedding_certificate() {
        let source = RationalSplineCurve::new(
            1,
            [0, 0, 1, 1, 2, 2].map(Rational::from_i64).to_vec(),
            [[0, 0, 0], [2, 0, 0], [1, 0, 0], [3, 0, 0]]
                .map(|p| p.map(Rational::from_i64))
                .to_vec(),
            vec![Rational::one(); 4],
            &Decider::default(),
        )
        .unwrap();
        assert!(
            prove(&source, &mut Decider::default())
                .unwrap_err()
                .contains("C0")
        );
    }
    #[test]
    fn rational_blossoms_match_the_independent_kernel_extraction() {
        let source = RationalSplineCurve::new(
            2,
            [0, 1, 2, 3, 4, 5, 6].map(Rational::from_i64).to_vec(),
            [[0, 0, 0], [1, 1, 0], [2, 1, 0], [3, 0, 0]]
                .map(|p| p.map(Rational::from_i64))
                .to_vec(),
            [1, 2, 3, 4].map(Rational::from_i64).to_vec(),
            &Decider::default(),
        )
        .unwrap();
        let curve = source.to_curve(&mut Decider::default()).unwrap();
        let mut dec = Decider::default();
        for (span, _) in curve.spans() {
            let net = bezier(&source, span, &mut dec).unwrap();
            let reference = curve.span_bezier_controls(span, &dec).unwrap();
            let factor = &net[0][3] / &reference[0][3].as_rational().unwrap();
            assert_eq!(net.len(), reference.len());
            for (a, b) in net.iter().zip(reference) {
                for i in 0..4 {
                    assert_eq!(a[i], &b[i].as_rational().unwrap() * &factor);
                }
            }
        }
    }
    #[test]
    fn subdivision_limits_refuse_a_cone_that_needs_a_finer_certificate() {
        let epsilon = Rational::parse(&format!("1/1{}", "0".repeat(100))).unwrap();
        let p = [
            Rational::from_ratio(1, 9).unwrap(),
            Rational::from_ratio(-2, 9).unwrap(),
            Rational::from_ratio(4, 9).unwrap(),
        ]
        .map(|x| &x + &epsilon)
        .to_vec();
        let error = positive(p, &mut 0, &mut Decider::default()).unwrap_err();
        assert!(error.contains("depth-24"), "{error}");
        let mut cells = CELLS;
        let error = positive(
            vec![Rational::one(); 2],
            &mut cells,
            &mut Decider::default(),
        )
        .unwrap_err();
        assert!(error.contains("4096-cell"), "{error}");
    }
}
