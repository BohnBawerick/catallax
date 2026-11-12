//! Pointwise certificates for rational stereographic circle replacements.

use catallax_numbers::{Decider, Rational, Sign};

use crate::incidence_arithmetic::{
    arithmetic as op, binomial, checkpoint, compare, equal, maximum, same, sqrt_enclosure,
};
use crate::incidence_repair::{Net, net_bounds, split};
use crate::{OrientedEdge, Vec3};

pub(crate) type Poly = Vec<Rational>;

pub(crate) fn product(a: &Poly, b: &Poly, dec: &mut Decider) -> Result<Poly, String> {
    checkpoint(dec)?;
    let (m, n) = (a.len() - 1, b.len() - 1);
    let mut out = vec![Rational::zero(); m + n + 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            checkpoint(dec)?;
            let factor = op(&binomial(m, i, dec)?, &binomial(n, j, dec)?, '*', dec)?;
            let factor = op(&factor, &binomial(m + n, i + j, dec)?, '/', dec)?;
            let value = op(x, y, '*', dec)?;
            let value = op(&value, &factor, '*', dec)?;
            out[i + j] = op(&out[i + j], &value, '+', dec)?;
        }
    }
    Ok(out)
}

fn combine(a: &Poly, b: &Poly, operation: char, dec: &mut Decider) -> Result<Poly, String> {
    a.iter()
        .zip(b)
        .map(|(a, b)| op(a, b, operation, dec))
        .collect()
}
fn scale(a: &Poly, b: &Rational, dec: &mut Decider) -> Result<Poly, String> {
    a.iter().map(|a| op(a, b, '*', dec)).collect()
}
fn derivative(a: &Poly, dec: &mut Decider) -> Result<Poly, String> {
    if a.len() < 2 {
        return Err("circle parameter has degree zero".into());
    }
    let n = Rational::from_i64((a.len() - 1) as i64);
    a.windows(2)
        .map(|p| {
            let d = op(&p[1], &p[0], '-', dec)?;
            op(&d, &n, '*', dec)
        })
        .collect()
}

/// Every coefficient has one sign, with at least one strict coefficient.
fn direction(a: &Poly) -> Result<Sign, String> {
    let mut sign = Sign::Zero;
    for x in a {
        let s = x.sign();
        if s == Sign::Zero {
            continue;
        }
        if sign != Sign::Zero && s != sign {
            return Err("circle angular monotonicity coefficients have mixed signs".into());
        }
        sign = s;
    }
    if sign == Sign::Zero {
        return Err("circle replacement has a constant angular span".into());
    }
    Ok(sign)
}

pub(crate) fn bound(
    nets: Vec<Net>,
    cap: &Rational,
    dec: &mut Decider,
) -> Result<(Rational, Rational), String> {
    let mut cells = Vec::new();
    if nets.is_empty() || nets.len() > 4096 {
        return Err("circle displacement certificate: initial cell allowance exceeded".into());
    }
    for net in nets {
        let bounds = net_bounds(&net, dec)?;
        cells.push((net, bounds, 0));
    }
    let mut visited = cells.len();
    let mut lower = Rational::zero();
    loop {
        checkpoint(dec)?;
        let mut upper = Rational::zero();
        let mut index = 0;
        for (i, (_, (lo, hi), _)) in cells.iter().enumerate() {
            lower = maximum(lower, lo.clone(), dec)?;
            if compare(hi, &upper, dec)?.is_ge() {
                upper = hi.clone();
                index = i;
            }
        }
        let lo = sqrt_enclosure(&lower, dec)?.0;
        let hi = sqrt_enclosure(&upper, dec)?.1;
        if compare(&hi, cap, dec)?.is_le()
            || compare(&lo, cap, dec)?.is_gt()
            || equal(&lower, &upper, dec)?
            || visited + 2 > 4096
            || cells[index].2 == 24
        {
            return Ok((lo, hi));
        }
        let (net, _, depth) = cells.swap_remove(index);
        let (a, b) = split(&net, dec)?;
        for net in [a, b] {
            let bounds = net_bounds(&net, dec)?;
            cells.push((net, bounds, depth + 1));
        }
        visited += 2;
    }
}

/// Coordinates are homogeneous in an exact orthonormal circle frame. A
/// positive chart denominator and a strictly ordered stereographic parameter
/// prove the native finite arc is the whole target image, in both directions.
pub(crate) fn image(
    nets: &[Net],
    radius: &Rational,
    dec: &mut Decider,
) -> Result<(Sign, Vec<Net>), String> {
    checkpoint(dec)?;
    let mut differences = Vec::new();
    let mut orientation = None;
    let mut previous: Option<Vec3> = None;
    for net in nets {
        checkpoint(dec)?;
        if net.len() < 2 {
            return Err("circle source span has no finite parameter interval".into());
        }
        let p: [Poly; 4] = std::array::from_fn(|i| net.iter().map(|p| p[i].clone()).collect());
        for w in &p[3] {
            if w.sign() != Sign::Positive {
                return Err("circle source weights are not positive".into());
            }
        }
        let rw = scale(&p[3], radius, dec)?;
        let q = combine(&rw, &p[0], '+', dec)?;
        if q.iter().any(|x| x.sign() != Sign::Positive) {
            return Err("circle replacement needs one positive stereographic chart".into());
        }
        let dy = derivative(&p[1], dec)?;
        let dq = derivative(&q, dec)?;
        let dyq = product(&dy, &q, dec)?;
        let ydq = product(&p[1], &dq, dec)?;
        let sign = direction(&combine(&dyq, &ydq, '-', dec)?)?;
        let dx = derivative(&p[0], dec)?;
        let xdy = product(&p[0], &dy, dec)?;
        let ydx = product(&p[1], &dx, dec)?;
        if direction(&combine(&xdy, &ydx, '-', dec)?)? != sign {
            return Err(
                "circle source polar monotonicity certificate disagrees with the target".into(),
            );
        }
        if p[2].iter().any(|z| !z.is_zero()) {
            return Err("circle source cap certificate: source boundary is not exactly in its unchanged cap plane".into());
        }
        if orientation.is_some_and(|old| old != sign) {
            return Err("circle replacement reverses angular direction between spans".into());
        }
        orientation = Some(sign);
        let qq = product(&q, &q, dec)?;
        let yy = product(&p[1], &p[1], dec)?;
        let denominator = combine(&qq, &yy, '+', dec)?;
        if denominator.iter().any(|x| x.sign() != Sign::Positive) {
            return Err("circle target positive-weight certificate is unproved".into());
        }
        let x = scale(&combine(&qq, &yy, '-', dec)?, radius, dec)?;
        let y = product(&q, &p[1], dec)?;
        let twice_r = op(radius, &Rational::from_i64(2), '*', dec)?;
        let y = scale(&y, &twice_r, dec)?;
        let wd = product(&p[3], &denominator, dec)?;
        let mut numerators = Vec::new();
        for (i, target) in [x, y, vec![Rational::zero(); denominator.len()]]
            .iter()
            .enumerate()
        {
            let source = product(&p[i], &denominator, dec)?;
            let target = product(target, &p[3], dec)?;
            numerators.push(combine(&source, &target, '-', dec)?);
        }
        let diff = (0..wd.len())
            .map(|i| {
                [
                    numerators[0][i].clone(),
                    numerators[1][i].clone(),
                    numerators[2][i].clone(),
                    wd[i].clone(),
                ]
            })
            .collect();
        differences.push(diff);
        let endpoint = |point: &[Rational; 4], dec: &mut Decider| -> Result<Vec3, String> {
            Ok([
                op(&point[0], &point[3], '/', dec)?,
                op(&point[1], &point[3], '/', dec)?,
                op(&point[2], &point[3], '/', dec)?,
            ])
        };
        let first = endpoint(&net[0], dec)?;
        if let Some(old) = &previous
            && !same(old, &first, dec)?
        {
            return Err("circle source spans are not exactly continuous".into());
        }
        previous = Some(endpoint(&net[net.len() - 1], dec)?);
    }
    let sign = orientation.ok_or("empty circle span set")?;
    Ok((sign, differences))
}

pub(crate) fn source_nets(edge: &OrientedEdge, dec: &mut Decider) -> Result<Vec<Net>, String> {
    let crate::Curve::BSpline(s) = &edge.curve else {
        return Err("expected source spline".into());
    };
    let g = &s.geometry;
    let mut nets = Vec::new();
    for span in g.degree()..g.controls().len() {
        if compare(&g.knots()[span], &g.knots()[span + 1], dec)?.is_lt() {
            nets.push(crate::spline_embedding::bezier(g, span, dec)?);
        }
    }
    let first = nets
        .first()
        .and_then(|n| n.first())
        .ok_or("empty source spline interval")?;
    let last = nets
        .last()
        .and_then(|n| n.last())
        .ok_or("empty source spline interval")?;
    let ends = if edge.curve_same_sense {
        [&edge.start.point, &edge.end.point]
    } else {
        [&edge.end.point, &edge.start.point]
    };
    for (point, vertex) in [(first, ends[0]), (last, ends[1])] {
        if point[3].sign() != Sign::Positive {
            return Err("source endpoint weight is not positive".into());
        }
        let actual = [
            op(&point[0], &point[3], '/', dec)?,
            op(&point[1], &point[3], '/', dec)?,
            op(&point[2], &point[3], '/', dec)?,
        ];
        if !same(&actual, vertex, dec)? {
            return Err("curved endpoint correspondence: complete source interval does not have the authored vertices as exact endpoints".into());
        }
    }
    Ok(nets)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;
    fn quarter() -> Net {
        let r = |x: &str| Rational::parse(x).unwrap();
        [
            ["1", "0", "0", "1"],
            ["1", "0.55228475", "0", "1"],
            ["0.55228475", "1", "0", "1"],
            ["0", "1", "0", "1"],
        ]
        .map(|p| p.map(r))
        .to_vec()
    }
    fn run(
        net: Net,
        budget: u64,
    ) -> (
        Result<(Sign, Rational, Rational), String>,
        catallax_numbers::OperationWork,
    ) {
        with_operation_budget(budget, || {
            let mut dec = Decider::default();
            let (sign, nets) = image(&[net], &Rational::one(), &mut dec)?;
            let (lo, hi) = bound(nets, &Rational::parse("0.01").unwrap(), &mut dec)?;
            Ok((sign, lo, hi))
        })
    }
    #[test]
    fn circle_certificate_charges_every_operation_and_operand_size() {
        let (answer, work) = run(quarter(), 20_000_000);
        let (_, lo, hi) = answer.unwrap();
        assert!(lo.sign() == Sign::Positive);
        assert!(hi <= Rational::parse("0.01").unwrap());
        assert!(work.used > 0);
        for limit in [0, work.used / 2, work.used - 1] {
            let (result, work) = run(quarter(), limit);
            assert!(result.is_err());
            assert!(work.exhausted);
        }
        let scale = Rational::parse(&format!("1{}", "0".repeat(300))).unwrap();
        let mut net = quarter();
        for p in &mut net {
            for x in p {
                *x = &*x * &scale;
            }
        }
        let (larger, larger_work) = run(net, 20_000_000);
        assert_eq!(larger.unwrap(), (Sign::Positive, lo, hi));
        assert!(larger_work.used > work.used);
    }
    #[test]
    fn circle_certificate_preserves_reverse_direction_and_refuses_a_fold() {
        let mut net = quarter();
        net.reverse();
        assert_eq!(run(net, 20_000_000).0.unwrap().0, Sign::Negative);
        let mut net = quarter();
        net[2][1] = Rational::from_i64(-2);
        assert!(run(net, 20_000_000).0.is_err());
    }
    #[test]
    fn pointwise_displacement_encloses_independent_rational_samples() {
        let net = quarter();
        let (_, lo, hi) = run(net.clone(), 20_000_000).0.unwrap();
        let lo2 = &lo * &lo;
        let hi2 = &hi * &hi;
        let mut maximum = Rational::zero();
        for i in 0..=16 {
            let t = Rational::from_ratio(i, 16).unwrap();
            let one_minus = &Rational::one() - &t;
            let mut p = net.clone();
            while p.len() > 1 {
                p = p
                    .windows(2)
                    .map(|pair| {
                        std::array::from_fn(|i| &(&pair[0][i] * &one_minus) + &(&pair[1][i] * &t))
                    })
                    .collect();
            }
            let x = &p[0][0] / &p[0][3];
            let y = &p[0][1] / &p[0][3];
            let q = &Rational::one() + &x;
            let d = &(&q * &q) + &(&y * &y);
            let tx = &(&(&q * &q) - &(&y * &y)) / &d;
            let ty = &(&Rational::from_i64(2) * &(&q * &y)) / &d;
            let dx = &x - &tx;
            let dy = &y - &ty;
            let distance = &(&dx * &dx) + &(&dy * &dy);
            assert!(distance <= hi2);
            maximum = maximum.max(distance);
        }
        assert!(maximum >= lo2);
    }
}
