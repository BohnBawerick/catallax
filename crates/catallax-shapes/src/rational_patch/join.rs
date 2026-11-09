//! C0 knot joins must prove a continuous regular tangent plane.

use super::*;

fn pair(
    carrier: &BSplineSurface,
    a: &Span,
    b: &Span,
    axis: usize,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let knot = &a.domain[axis][1];
    let knots = carrier.knots()[axis];
    let begin = knots.partition_point(|value| value < knot);
    let end = knots.partition_point(|value| value <= knot);
    if end - begin < carrier.degrees()[axis] {
        // Positive-weight B-spline quotient charts are C1 here by construction.
        return Ok(());
    }
    let value = RealPoly::constant(2, Real::from_rational(knot.clone()));
    let left = a
        .normal_polynomials()
        .each_ref()
        .map(|p| p.substitute(axis, &value));
    let right = b
        .normal_polynomials()
        .each_ref()
        .map(|p| p.substitute(axis, &value));
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        if !(&left[i] * &right[j] - &left[j] * &right[i]).is_zero(dec)? {
            return Err(unproven(
                "an internal knot has distinct tangent planes; split the crease into separate patch faces",
            ));
        }
    }
    // The shared positive projection determinant fixes the same normal side,
    // so parallel nonzero normals cannot be oppositely oriented.
    Ok(())
}

pub(super) fn check(
    carrier: &BSplineSurface,
    spans: &[Span],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let nv = carrier.spans(1).count();
    let nu = spans.len() / nv;
    for u in 0..nu {
        for v in 0..nv {
            dec.check_deadline()?;
            let current = &spans[u * nv + v];
            if u + 1 < nu {
                pair(carrier, current, &spans[(u + 1) * nv + v], 0, dec)?;
            }
            if v + 1 < nv {
                pair(carrier, current, &spans[u * nv + v + 1], 1, dec)?;
            }
        }
    }
    Ok(())
}
