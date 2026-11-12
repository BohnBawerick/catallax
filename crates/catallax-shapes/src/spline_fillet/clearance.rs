//! Each contact is its support's closest point to the ball's centre.
//!
//! Each contact's parameters are first proved inside its support's domain,
//! so each rail lies on its face. For `tau` in a span and the centre `C` there, let
//! `g(q) = |C - P(q)|^2` over the support's closed domain. On a box `U`
//! round the contact parameters, `g` is proved convex by a positive definite
//! Hessian over the whole box. Its gradient vanishes at the contact, where
//! `C - P` is normal to the support, so the contact is the minimum of `g` on
//! `U`, where `g = r^2`. Off `U`, cells of `tau` and of the domain are
//! excluded by `|C(t) - P| > r + m`, where `t` is the cell's middle `tau` and
//! `m` bounds how far the centre moves over the cell. So the ball meets the
//! support's face only at the contact.
//!
//! Spine and centre over the span are second-order Taylor forms about the
//! span's middle: the value and slope there from a point solve, and the
//! span's proved second Taylor coefficient as the remainder.

use super::bound::{SpanProof, root_up};
use super::jets::{Meter, ball, magnitude};
use super::*;

/// Most halvings of one exclusion cell.
const MAX_DEPTH: u32 = 24;

/// Most quarterings of the box round a contact in the convexity proof.
const CONVEX_DEPTH: u32 = 3;

/// Proves `g` convex on a parameter box: its Hessian is positive definite
/// at every point, proved on sub-boxes, each quartered when it fails.
fn convex_on(
    support: &Support,
    centre: &[Ball; 3],
    around: &[[Rational; 2]; 2],
    meter: &mut Meter,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let two = Rational::from_i64(2);
    let mut pending = vec![(around.clone(), 0u32)];
    while let Some((cell, depth)) = pending.pop() {
        let boxed: [Ball; 2] = std::array::from_fn(|i| ball(&cell[i][0]).hull(&ball(&cell[i][1])));
        if convex(support, centre, &boxed, meter, dec)? {
            continue;
        }
        if depth >= CONVEX_DEPTH {
            return Ok(false);
        }
        let middle = cell
            .each_ref()
            .map(|side| (&side[0] + &side[1]) / two.clone());
        for a in 0..2 {
            for b in 0..2 {
                let u = if a == 0 {
                    [cell[0][0].clone(), middle[0].clone()]
                } else {
                    [middle[0].clone(), cell[0][1].clone()]
                };
                let v = if b == 0 {
                    [cell[1][0].clone(), middle[1].clone()]
                } else {
                    [middle[1].clone(), cell[1][1].clone()]
                };
                pending.push(([u, v], depth + 1));
            }
        }
    }
    Ok(true)
}

/// Proves the Hessian of `g` positive definite over a parameter box, from
/// the support's own derivatives there:
/// `H_ab = 2 (P_a . P_b - (C - P) . P_ab)`.
fn convex(
    support: &Support,
    centre: &[Ball; 3],
    around: &[Ball; 2],
    meter: &mut Meter,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let one = ball(&Rational::one());
    let along = |directions: [bool; 2], meter: &mut Meter| -> Result<[Jet; 3], Refusal> {
        let at: [Jet; 2] = std::array::from_fn(|i| {
            if directions[i] {
                Jet::linear(around[i].clone(), one.clone(), 2)
            } else {
                Jet::constant(around[i].clone(), 2)
            }
        });
        support.point(&at, meter)
    };
    let pu = along([true, false], meter)?;
    let pv = along([false, true], meter)?;
    let pd = along([true, true], meter)?;
    meter.charge(dec)?;
    let two = ball(&Rational::from_i64(2));
    let mut entries = [
        ball(&Rational::zero()),
        ball(&Rational::zero()),
        ball(&Rational::zero()),
    ];
    for i in 0..3 {
        let gap = &centre[i] - &pu[i].0[0];
        let (du, dv) = (&pu[i].0[1], &pv[i].0[1]);
        let duu = &two * &pu[i].0[2];
        let dvv = &two * &pv[i].0[2];
        let duv = &(&pd[i].0[2] - &pu[i].0[2]) - &pv[i].0[2];
        entries[0] = &entries[0] + &(&du.pow_u32(2) - &(&gap * &duu));
        entries[1] = &entries[1] + &(&dv.pow_u32(2) - &(&gap * &dvv));
        entries[2] = &entries[2] + &(&(du * dv) - &(&gap * &duv));
    }
    meter.balls(24, true);
    meter.balls(30, false);
    meter.charge(dec)?;
    let [huu, hvv, huv] = entries.map(|entry| &two * &entry);
    let (a, c) = (huu.lo(), hvv.lo());
    if a.sign() != Sign::Positive || c.sign() != Sign::Positive {
        return Ok(false);
    }
    let off = magnitude(&huv);
    Ok(a * c > off.clone() * off)
}

fn inside(cell: &[[Rational; 2]; 2], around: &[[Rational; 2]; 2]) -> bool {
    (0..2).all(|i| cell[i][0] >= around[i][0] && cell[i][1] <= around[i][1])
}

/// A second-order Taylor form about the span's middle: value and slope
/// there, and an enclosure of half the second derivative over the span.
struct Taylor {
    value: [Ball; 3],
    slope: [Ball; 3],
    curve: [Ball; 3],
}

impl Taylor {
    fn new(point: &[Jet; 3], span: &[Jet; 3]) -> Self {
        Self {
            value: point.each_ref().map(|j| j.0[0].clone()),
            slope: point.each_ref().map(|j| j.0[1].clone()),
            curve: span.each_ref().map(|j| j.0[2].clone()),
        }
    }

    /// The enclosure at offset `m` from the middle.
    fn at(&self, m: &Rational, meter: &mut Meter) -> [Ball; 3] {
        let (m, m2) = (ball(m), ball(&(m.clone() * m.clone())));
        meter.balls(6, true);
        meter.balls(6, false);
        std::array::from_fn(|i| {
            &(&self.value[i] + &(&self.slope[i] * &m)) + &(&self.curve[i] * &m2)
        })
    }

    /// The enclosure over offsets in `[-h, h]`, met with `over`, an
    /// enclosure over the span.
    fn over(&self, h: &Rational, over: &[Ball; 3], meter: &mut Meter) -> [Ball; 3] {
        let m = ball(&-h.clone()).hull(&ball(h));
        let m2 = ball(&Rational::zero()).hull(&ball(&(h.clone() * h.clone())));
        meter.balls(6, true);
        meter.balls(6, false);
        std::array::from_fn(|i| {
            let value = &(&self.value[i] + &(&self.slope[i] * &m)) + &(&self.curve[i] * &m2);
            let (lo, hi) = (value.lo().max(over[i].lo()), value.hi().min(over[i].hi()));
            if lo <= hi {
                ball(&lo).hull(&ball(&hi))
            } else {
                value
            }
        })
    }
}

/// Every cell of `tau` offsets in `[-h, h]` and of the domain off `around`
/// keeps the support farther than `r` from the ball's centre. A cell's
/// centre is the Taylor form at its middle offset, widened by `speed` times
/// its half width; a cell splits in `tau` when that widening exceeds its
/// point box, and otherwise along its longer parameter side.
#[allow(clippy::too_many_arguments)]
fn excluded(
    support: &Support,
    centre: &Taylor,
    speed: &Rational,
    h: &Rational,
    radius: &Rational,
    around: &[[Rational; 2]; 2],
    meter: &mut Meter,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let two = Rational::from_i64(2);
    let mut pending = vec![(support.domain.clone(), [-h.clone(), h.clone()], 0u32)];
    while let Some((cell, offsets, depth)) = pending.pop() {
        dec.check_deadline()?;
        if inside(&cell, around) {
            continue;
        }
        let middle = (&offsets[0] + &offsets[1]) / two.clone();
        let fixed = centre.at(&middle, meter);
        let drift = speed.clone() * (&offsets[1] - &offsets[0]) / two.clone();
        let reach = radius.clone() + drift.clone();
        let at: [Jet; 2] =
            std::array::from_fn(|i| Jet::constant(ball(&cell[i][0]).hull(&ball(&cell[i][1])), 0));
        let point = support.point(&at, meter)?;
        meter.charge(dec)?;
        let mut lower = Rational::zero();
        let mut size = Rational::zero();
        for i in 0..3 {
            let p = &point[i].0[0];
            let gap = (fixed[i].lo() - p.hi()).max(p.lo() - fixed[i].hi());
            if gap.sign() == Sign::Positive {
                lower = lower + gap.clone() * gap;
            }
            size = size.max(p.hi() - p.lo());
        }
        if lower > reach.clone() * reach {
            continue;
        }
        if depth >= MAX_DEPTH {
            return Err(Refusal::new(
                SPLINE_FILLET_CLEARANCE,
                format!(
                    "the ball was not proved clear of its support over u [{}, {}] v [{}, {}]",
                    cell[0][0], cell[0][1], cell[1][0], cell[1][1]
                ),
            ));
        }
        if drift > size {
            pending.push((
                cell.clone(),
                [middle.clone(), offsets[1].clone()],
                depth + 1,
            ));
            pending.push((cell, [offsets[0].clone(), middle], depth + 1));
            continue;
        }
        let axis = usize::from(
            cell[1][1].clone() - cell[1][0].clone() > cell[0][1].clone() - cell[0][0].clone(),
        );
        let split = (&cell[axis][0] + &cell[axis][1]) / two.clone();
        let mut low = cell.clone();
        let mut high = cell;
        low[axis][1] = split.clone();
        high[axis][0] = split;
        pending.push((high, offsets.clone(), depth + 1));
        pending.push((low, offsets, depth + 1));
    }
    Ok(())
}

fn span(
    exact: &Exact,
    ends: [&Rational; 2],
    proof: &SpanProof,
    estimates: [[f64; 3]; 2],
    index: usize,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut meter = Meter::default();
    let two = Rational::from_i64(2);
    let middle = (ends[0] + ends[1]) / two.clone();
    let h = (ends[1] - ends[0]) / two;
    let halfway: [f64; 3] = std::array::from_fn(|i| 0.5 * (estimates[0][i] + estimates[1][i]));
    let point = exact.solve([&middle, &middle], [halfway, halfway], 1, dec)?;
    let point_rows = exact.rows(&point, &mut meter)?;
    meter.charge(dec)?;
    let centre = Taylor::new(&point_rows.centre, &proof.centre);
    let spine = Taylor::new(&point.x, &proof.x);
    let centre_box = centre.over(
        &h,
        &proof.centre.each_ref().map(|j| j.0[0].clone()),
        &mut meter,
    );
    let x = spine.over(&h, &proof.x.each_ref().map(|j| j.0[0].clone()), &mut meter);
    let speed = root_up(
        &proof
            .centre
            .iter()
            .map(|j| magnitude(&j.0[1]).pow(2).expect("square"))
            .fold(Rational::zero(), |s, v| s + v),
        dec,
    )?;
    let params = exact.params(&x, &proof.tau.0[0]);
    for (k, support) in exact.supports.iter().enumerate() {
        let within = (0..2).all(|i| {
            params[k][i].lo() >= support.domain[i][0] && params[k][i].hi() <= support.domain[i][1]
        });
        if !within {
            return Err(Refusal::new(
                SPLINE_FILLET_SPINE,
                format!("support {k}'s contact was not proved inside its face on span {index}"),
            ));
        }
        let mut delta = Rational::from_ratio(1, 4).expect("nonzero");
        let around = loop {
            let around: [[Rational; 2]; 2] = std::array::from_fn(|i| {
                let [lo, hi] = &support.domain[i];
                [
                    (params[k][i].lo() - delta.clone()).max(lo.clone()),
                    (params[k][i].hi() + delta.clone()).min(hi.clone()),
                ]
            });
            if convex_on(support, &centre_box, &around, &mut meter, dec)? {
                break around;
            }
            delta = delta.mul_pow2(-1);
            if delta < Rational::one().mul_pow2(-8) {
                return Err(Refusal::new(
                    SPLINE_FILLET_CLEARANCE,
                    format!("support {k} was not proved convex round its contact on span {index}"),
                ));
            }
        };
        excluded(
            support,
            &centre,
            &speed,
            &h,
            &exact.radius,
            &around,
            &mut meter,
            dec,
        )?;
    }
    Ok(())
}

/// Proves every span's contacts closest points of their supports.
pub(super) fn check(
    exact: &Exact,
    nodes: &[Node],
    estimates: &[[f64; 3]],
    proofs: &[SpanProof],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    order::in_order(
        proofs.len(),
        |i, dec| {
            span(
                exact,
                [&nodes[i].tau, &nodes[i + 1].tau],
                &proofs[i],
                [estimates[i], estimates[i + 1]],
                i,
                dec,
            )
        },
        dec,
    )
    .map(|_| ())
}
