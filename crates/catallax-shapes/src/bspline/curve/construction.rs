//! Positive homogeneous construction and exact C0 concatenation.

use super::super::validate_homogeneous;
use super::*;
use crate::rational_patch::exact_sign;

fn elevate(mut net: Vec<[Real; 4]>, degree: usize) -> Vec<[Real; 4]> {
    while net.len() <= degree {
        let n = net.len();
        let mut next = Vec::with_capacity(n + 1);
        next.push(net[0].clone());
        for i in 1..n {
            let alpha = Real::from_rational(
                Rational::from_ratio(i as i64, n as i64).expect("positive degree"),
            );
            let beta = &Real::one() - &alpha;
            next.push(std::array::from_fn(|axis| {
                &net[i - 1][axis] * &alpha + &net[i][axis] * &beta
            }));
        }
        next.push(net[n - 1].clone());
        net = next;
    }
    net
}

fn normalize(net: Vec<[Real; 4]>, dec: &mut Decider) -> Result<Vec<[Real; 4]>, Refusal> {
    let first = &net[0][3];
    let last = &net[net.len() - 1][3];
    if exact_sign(&(first - last), dec)? == Sign::Zero {
        return Ok(net
            .iter()
            .map(|p| p.each_ref().map(|x| x / first))
            .collect());
    }
    let n = net.len();
    let mut result = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let alpha =
            Real::from_rational(Rational::from_ratio(i as i64, n as i64).expect("positive degree"));
        let beta = &Real::one() - &alpha;
        result.push(std::array::from_fn(|axis| {
            let low = if i < n {
                &(&net[i][axis] / first) * &beta
            } else {
                Real::zero()
            };
            let high = if i > 0 {
                &(&net[i - 1][axis] / last) * &alpha
            } else {
                Real::zero()
            };
            low + high
        }));
    }
    Ok(result)
}

impl BSplineCurve {
    /// Validates positive homogeneous controls without reweighting their exact
    /// numerators. Knot, coordinate and weight rules match `new`.
    pub fn from_homogeneous(
        degree: usize,
        knots: Vec<Rational>,
        homogeneous: Vec<[Real; 4]>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        validate_axis(degree, &knots, homogeneous.len(), dec)?;
        let data = validate_homogeneous(&homogeneous, dec)?;
        let net = Net {
            degrees: [degree, 0],
            knots: [knots.into(), Arc::from([Rational::zero(), Rational::one()])],
            count: [homogeneous.len(), 1],
            controls: homogeneous,
        };
        Ok(Self::from_net(
            net,
            data.controls,
            data.weights,
            data.lower,
            data.constant,
            None,
        ))
    }

    /// Concatenates exact complete pieces in the supplied directions. Each
    /// original span receives one unit interval. Positive linear homogeneous
    /// factors normalize joining weights without moving the rational curve.
    /// Every common endpoint is proved exactly equal before sharing its node.
    pub fn joined(pieces: &[(BSplineCurve, bool)], dec: &mut Decider) -> Result<Self, Refusal> {
        let mut nets = Vec::new();
        for (curve, forward) in pieces {
            let mut spans: Vec<_> = curve.spans().map(|(span, _)| span).collect();
            if !forward {
                spans.reverse();
            }
            for span in spans {
                let mut net = curve.span_bezier_controls(span, dec)?;
                if !forward {
                    net.reverse();
                }
                nets.push(normalize(net, dec)?);
            }
        }
        let degree = nets
            .iter()
            .map(|net| net.len() - 1)
            .max()
            .ok_or_else(|| invalid("a joined curve needs at least one piece"))?;
        let mut controls: Vec<[Real; 4]> = Vec::new();
        let mut knots = vec![Rational::zero(); degree + 1];
        let count = nets.len();
        for (index, net) in nets.into_iter().enumerate() {
            let mut net = elevate(net, degree);
            if let Some(previous) = controls.last() {
                for axis in 0..3 {
                    if exact_sign(&(&previous[axis] - &net[0][axis]), dec)? != Sign::Zero {
                        return Err(invalid("joined curve endpoints differ exactly"));
                    }
                }
                net[0] = previous.clone();
                controls.extend(net.into_iter().skip(1));
            } else {
                controls.extend(net);
            }
            knots.extend(std::iter::repeat_n(
                Rational::from_i64((index + 1) as i64),
                degree + usize::from(index + 1 == count),
            ));
        }
        // A closed chain keeps the same exact seam node, not two equivalent
        // algebraic expressions with unrelated identity.
        let last = controls.len() - 1;
        let mut closed = true;
        for (first, last) in controls[0].iter().zip(&controls[last]).take(3) {
            closed &= exact_sign(&(first - last), dec)? == Sign::Zero;
        }
        if closed {
            controls[last] = controls[0].clone();
        }
        Self::from_homogeneous(degree, knots, controls, dec)
    }
}
