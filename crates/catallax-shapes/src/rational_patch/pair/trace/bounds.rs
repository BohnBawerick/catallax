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

pub(super) fn series(
    poly: &RealPoly,
    at: &[Vec<Ball>],
    order: usize,
    bits: u32,
    dec: &Decider,
) -> Result<Vec<Ball>, Refusal> {
    let middle: Vec<_> = at.iter().map(|row| row[0].mid()).collect();
    let mut local = at.to_vec();
    for (row, centre) in local.iter_mut().zip(&middle) {
        row[0] = &row[0] - &Ball::from_rational(centre, bits);
    }
    crate::rational_patch::trace_series::series(&shifted(poly, &middle), &local, order, bits, dec)
}
