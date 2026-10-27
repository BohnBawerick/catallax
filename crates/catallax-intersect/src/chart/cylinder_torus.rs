//! Rational half-angle charts with every vanishing endpoint factor extracted.

use super::*;
use crate::pair::CylinderTorusMode;

pub(super) fn charts(frame: &Frame, branch: i64) -> (Vec<Chart>, i64) {
    let Frame::CylinderTorus {
        origin,
        e1,
        e2,
        e3,
        radius,
        radial_shift,
        radial_slope,
        lo,
        span,
        lower_radius,
        upper_radius,
        y_ends,
        z_ends,
        mode,
        sense,
    } = frame
    else {
        unreachable!("cylinder/torus frame");
    };
    let (cos, sin) = trig();
    let x = E::constant(lo.clone()) + E::constant(span.clone()) * sin.clone().square();
    let y = match y_ends {
        [true, true] => E::constant(span.clone()) * sin.clone() * cos.clone(),
        [true, false] => {
            sin.clone()
                * (E::constant(span.clone()) * (E::constant(radius.clone()) - x.clone())).sqrt()
        }
        [false, true] => {
            cos.clone()
                * (E::constant(span.clone()) * (E::constant(radius.clone()) + x.clone())).sqrt()
        }
        [false, false] => ((E::constant(radius.clone()) - x.clone())
            * (E::constant(radius.clone()) + x.clone()))
        .sqrt(),
    };
    let radial =
        (E::constant(radial_shift.clone()) + E::constant(radial_slope.clone()) * x.clone()).sqrt();
    // rho-L = (rho²-L²)/(rho+L), and similarly for H-rho. The
    // radial endpoints are positive on a ring torus, so these denominators
    // stay positive while the extracted sine/cosine factors carry the zeros.
    let z = match z_ends {
        [true, true] => {
            E::constant(radial_slope * span) * sin * cos
                / ((radial.clone() + E::constant(lower_radius.clone()))
                    * (radial.clone() + E::constant(upper_radius.clone())))
                .sqrt()
        }
        [true, false] => {
            sin * (E::constant(radial_slope * span)
                * (E::constant(upper_radius.clone()) - radial.clone())
                / (radial.clone() + E::constant(lower_radius.clone())))
            .sqrt()
        }
        [false, true] => {
            cos * (E::constant(radial_slope * span)
                * (radial.clone() - E::constant(lower_radius.clone()))
                / (radial.clone() + E::constant(upper_radius.clone())))
            .sqrt()
        }
        [false, false] => ((E::constant(upper_radius.clone()) - radial.clone())
            * (radial - E::constant(lower_radius.clone())))
        .sqrt(),
    };
    let signs: &[(i64, i64, bool)] = match mode {
        CylinderTorusMode::Axial => &[(1, branch, false), (-1, branch, true)],
        CylinderTorusMode::Transverse => &[(branch, 1, false), (branch, -1, true)],
        CylinderTorusMode::Single if y_ends[1] => {
            &[(1, 1, false), (-1, 1, true), (-1, -1, false), (1, -1, true)]
        }
        CylinderTorusMode::Single => {
            &[(1, 1, false), (1, -1, true), (-1, -1, false), (-1, 1, true)]
        }
    };
    let mut charts = Vec::with_capacity(signs.len());
    let xcoord = lin(origin, e1);
    let ycoord = lin(origin, e2);
    let zcoord = lin(origin, e3);
    for &(ys, zs, reverse) in signs {
        let mut chart = Chart {
            coords: point(
                origin,
                e1,
                e2,
                e3,
                x.clone(),
                E::integer(ys) * y.clone(),
                E::integer(zs) * z.clone(),
            ),
            inverse: Inverse::Squared {
                value: xcoord.clone(),
                lo: lo.clone(),
                span: span.clone(),
            },
            checks: vec![(ycoord.clone(), ys), (zcoord.clone(), zs)],
            reverse: false,
        };
        if reverse {
            chart.reverse();
        }
        charts.push(chart);
    }
    // In this right-handed frame, grad(cylinder) x grad(torus) has
    // x-component with sign y*z. The first chart increases x.
    (charts, sense * signs[0].0 * signs[0].1)
}
