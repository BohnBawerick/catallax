//! Recenter power equations before enclosing the retained graph.

use super::*;

fn shifted(poly: &RealPoly, middle: &[Rational]) -> RealPoly {
    let mut local = poly.clone();
    for (axis, value) in middle.iter().enumerate() {
        local = local.substitute(
            axis,
            &(RealPoly::var(poly.nvars(), axis)
                + RealPoly::constant(poly.nvars(), Real::from(value.clone()))),
        );
    }
    local
}

pub(super) fn eval(poly: &RealPoly, at: &[Ball], dec: &Decider) -> Result<Ball, Refusal> {
    let middle: Vec<_> = at.iter().map(Ball::mid).collect();
    let local: Vec<_> = at
        .iter()
        .zip(&middle)
        .map(|(value, centre)| value - &Ball::from_rational(centre, value.precision()))
        .collect();
    shifted(poly, &middle).eval_ball(&local, dec)
}
