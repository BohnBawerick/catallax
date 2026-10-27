//! Positive Bernstein factors with simple endpoint roots extracted before jets.

use super::*;
use crate::pair::cone_torus::{ConeTorusFrame, Polynomial, Sector};

fn polynomial(polynomial: &Polynomial, u: &E, complement: &E) -> E {
    let [a, b, c] = &polynomial.bernstein;
    match polynomial.degree {
        0 => E::constant(a.clone()),
        1 => E::constant(a.clone()) * complement.clone() + E::constant(b.clone()) * u.clone(),
        2 => {
            let left =
                E::constant(a.clone()) * complement.clone() + E::constant(b.clone()) * u.clone();
            let right =
                E::constant(b.clone()) * complement.clone() + E::constant(c.clone()) * u.clone();
            left * complement.clone() + right * u.clone()
        }
        _ => unreachable!("quadratic Bernstein factor"),
    }
}
fn sector(frame: &ConeTorusFrame, sector: &Sector, sheet: i64) -> Chart {
    let odd = sector.orders.map(|order| order % 2 == 1);
    let t = E::t();
    let one_minus_t = E::integer(1) - t.clone();
    let (u, complement, left_root, right_root) = match odd {
        [true, true] => {
            let (cosine, sine) = trig();
            (
                sine.clone().square(),
                cosine.clone().square(),
                Some(sine),
                Some(cosine),
            )
        }
        [true, false] => (
            t.clone().square(),
            one_minus_t * (E::integer(1) + t.clone()),
            Some(t),
            None,
        ),
        [false, true] => (
            t.clone() * (E::integer(2) - t),
            one_minus_t.clone().square(),
            None,
            Some(one_minus_t),
        ),
        [false, false] => (t, one_minus_t, None, None),
    };
    let parameter = E::constant(sector.lo.clone()) + E::constant(sector.span.clone()) * u.clone();
    let denominator = E::integer(1) + parameter.clone().square();
    let cosine = (E::integer(1) - parameter.clone().square()) / denominator.clone();
    let sine = E::integer(2) * parameter / denominator.clone();
    let (cosine, sine) = match sector.quarter {
        0 => (cosine, sine),
        1 => (-sine, cosine),
        2 => (-cosine, -sine),
        3 => (sine, -cosine),
        _ => unreachable!("minor-angle quarter"),
    };
    let coordinates = &frame.coordinates;
    let rho =
        E::constant(coordinates.major.clone()) + E::constant(coordinates.minor.clone()) * cosine;
    let radius = E::constant(frame.radius_base.clone())
        + E::constant(frame.radius_sine.clone()) * sine.clone();
    let twice_distance = E::constant(Real::from_i64(2) * frame.distance.clone());
    let x = (rho.square() + E::constant(frame.distance.square()) - radius.square())
        / twice_distance.clone();
    let mut positive = polynomial(&sector.factors[0], &u, &complement);
    for factor in &sector.factors[1..] {
        positive = positive * polynomial(factor, &u, &complement);
    }
    let mut y = positive.sqrt();
    for (order, variable) in sector.orders.into_iter().zip([u, complement]) {
        if order / 2 != 0 {
            let distance = E::constant(sector.span.clone()) * variable;
            for _ in 0..order / 2 {
                y = y * distance.clone();
            }
        }
    }
    for root in [left_root, right_root].into_iter().flatten() {
        y = y * root;
    }
    match odd.iter().filter(|&&odd| odd).count() {
        0 => {}
        1 => y = y * E::constant(sector.span.sqrt()),
        2 => y = y * E::constant(sector.span.clone()),
        _ => unreachable!(),
    }
    y = E::integer(sheet) * y / (twice_distance * denominator.square());
    let z = E::constant(coordinates.minor.clone()) * sine;
    Chart {
        coords: point(
            &coordinates.origin,
            &coordinates.e1,
            &coordinates.e2,
            &coordinates.e3,
            x,
            y,
            z,
        ),
        inverse: Inverse::TorusArc {
            coordinates: Arc::clone(coordinates),
            quarter: sector.quarter,
            lo: sector.lo.clone(),
            span: sector.span.clone(),
            orders: sector.orders,
        },
        checks: vec![(lin(&coordinates.origin, &coordinates.e2), sheet)],
        reverse: false,
    }
}
fn reverse(charts: &mut [Chart]) {
    charts.reverse();
    for chart in charts {
        chart.reverse();
    }
}
pub(super) fn charts(frame: &ConeTorusFrame, component: usize) -> Vec<Chart> {
    let domain = &frame.domains[component];
    let copies = if domain.sheet == 0 { 2 } else { 1 };
    let mut charts = Vec::with_capacity(domain.count * copies);
    for offset in 0..domain.count {
        let piece = &frame.sectors[(domain.start + offset) % frame.sectors.len()];
        charts.push(sector(
            frame,
            piece,
            if domain.sheet == 0 { -1 } else { domain.sheet },
        ));
    }
    if domain.sheet == 0 {
        for offset in (0..domain.count).rev() {
            let piece = &frame.sectors[(domain.start + offset) % frame.sectors.len()];
            let mut chart = sector(frame, piece, 1);
            chart.reverse();
            charts.push(chart);
        }
    } else if domain.sheet > 0 {
        reverse(&mut charts);
    }
    // grad(cone) x grad(torus) = -16 R d y times increasing minor angle.
    if frame.sense < 0 {
        reverse(&mut charts);
    }
    charts
}
