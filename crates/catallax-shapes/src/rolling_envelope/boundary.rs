//! Parameter tables of the four boundary rows on their analytic neighbours.
//!
//! A row keeps its own chart parameter, so its neighbour sees an angle such as
//! `2 atan(t s)`, which no polynomial represents exactly. Each table is a
//! piecewise cubic Hermite fit with a bound recomputed from an interval jet.

use super::*;

/// One boundary row of the carrier's `[-1,1] x [0,1]` chart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryRow {
    /// `v = 0`, on the cylinder; the angle is the ball's station angle.
    Rail,
    /// `v = 1`, on the mouth circle; the angle is measured from the mouth's far
    /// point, so the near point has angle zero and adds a half turn.
    Mouth,
    /// `u = -1`, the transition meridian; the angle runs from the cylinder
    /// contact to the mouth plane contact.
    Negative,
    /// `u = 1`, the other transition meridian.
    Positive,
}

/// Piecewise cubic Hermite parameter data over strictly increasing knots.
#[derive(Clone, Debug)]
pub struct ParameterTable {
    /// Chart parameters of the row.
    pub knots: Vec<Rational>,
    /// Angle and its derivative at each knot, in radians.
    pub angle: Vec<[Rational; 2]>,
    /// Rail height above the mouth plane and its derivative at each knot,
    /// exactly the represented rail. Empty on every other row.
    pub height: Vec<[Rational; 2]>,
    /// Upper bound on the angle's distance from its table, in radians.
    pub bound: Rational,
}

impl RollingEnvelope {
    /// Fits the row's angle within `tolerance` radians. A row the jet cannot
    /// resolve within 2^12 pieces per chart span refuses.
    pub fn parameter_table(
        &self,
        row: BoundaryRow,
        tolerance: &Rational,
        dec: &Decider,
    ) -> Result<ParameterTable, Refusal> {
        let bits = self.data.net.bits;
        let def = &self.data.definition;
        let k = Constants::new(def, bits, dec)?;
        let nodes = &self.data.nodes;
        let (spans, enclosed): (Vec<[Rational; 2]>, Vec<RowData>) = match row {
            BoundaryRow::Rail | BoundaryRow::Mouth => (
                nodes
                    .windows(2)
                    .map(|pair| [pair[0].s.clone(), pair[1].s.clone()])
                    .collect(),
                nodes
                    .iter()
                    .map(|node| {
                        let mut out: RowData =
                            std::array::from_fn(|_| [k.one.clone(), k.one.clone()]);
                        for ((slot, value), slope) in
                            out.iter_mut().zip(&node.values).zip(&node.slopes)
                        {
                            *slot = [dec.enclose_at(value, bits)?, b(slope, bits)];
                        }
                        Ok(out)
                    })
                    .collect::<Result<_, Refusal>>()?,
            ),
            _ => (vec![[Rational::zero(), Rational::one()]], Vec::new()),
        };
        let transition = def.transition();
        let square = dec.enclose_at(&transition.square, bits)?;
        let ratio = dec.enclose_at(&transition.mouth, bits)?;
        // Quarter-arc conic: with q = rho v / (1 - v + rho v) and k = tan(pi/8),
        // the angle from the rail contact is pi/4 + 2 atan(k (2 q - 1)).
        let rho = {
            let w2 = &k.one + &(&square * &ratio.square());
            w2.checked_div(&(&k.one + &square))
                .as_ref()
                .and_then(Ball::sqrt)
                .ok_or_else(|| unresolved("transition weight ratio is unresolved"))?
        };
        let half = b(&q(2), bits)
            .sqrt()
            .ok_or_else(|| unresolved("no enclosure of sqrt(2)"))?;
        let eighth = &half - &k.one;
        let quarter = Ball::pi_checked(bits, dec)? * b(&Rational::from_ratio(1, 4).unwrap(), bits);
        // The tangent half-angle argument of the row, as a jet.
        let argument = |s: &Jet, index: usize| -> Result<Jet, Refusal> {
            match row {
                BoundaryRow::Rail => Ok(s.scale(&k.t)),
                BoundaryRow::Mouth => {
                    let span = [&spans[index][0], &spans[index][1]];
                    hermite(
                        s,
                        span,
                        [
                            [&enclosed[index][MOUTH][0], &enclosed[index][MOUTH][1]],
                            [
                                &enclosed[index + 1][MOUTH][0],
                                &enclosed[index + 1][MOUTH][1],
                            ],
                        ],
                        bits,
                    )
                }
                _ => {
                    let scaled = s.scale(&rho);
                    let denominator = Jet::constant(k.one.clone()).add(&s.scale(&(&rho - &k.one)));
                    let share = scaled.div(&denominator)?;
                    Ok(share
                        .scale(&k.two)
                        .sub(&Jet::constant(k.one.clone()))
                        .scale(&eighth))
                }
            }
        };
        let angle = |s: &Jet, index: usize| -> Result<[Ball; 2], Refusal> {
            let f = argument(s, index)?;
            let value = (&k.two * &f.0[0].atan_checked(dec)?)
                + if matches!(row, BoundaryRow::Negative | BoundaryRow::Positive) {
                    quarter.clone()
                } else {
                    b(&Rational::zero(), bits)
                };
            let slope = slope_jet(&f, &k)?;
            Ok([value, slope.0[0].clone()])
        };
        let mut pieces = 1i64;
        loop {
            dec.check_deadline()?;
            let mut worst = Rational::zero();
            for (index, span) in spans.iter().enumerate() {
                let h = (&span[1] - &span[0]) / q(pieces);
                let h4 = h.pow(4).expect("nonzero span");
                for piece in 0..pieces {
                    dec.check_deadline()?;
                    let lo = &span[0] + &(&h * &q(piece));
                    let hi = &lo + &h;
                    let s = Jet::linear(b(&lo, bits).hull(&b(&hi, bits)), k.one.clone());
                    let fourth = slope_jet(&argument(&s, index)?, &k)?;
                    // g'''' = 6 c3 and the cubic remainder is h^4 |g''''| / 384.
                    worst = worst.max(magnitude(&fourth.0[3]) * h4.clone() / q(64));
                }
            }
            if worst <= *tolerance || pieces >= 4096 {
                if worst > *tolerance {
                    return Err(unresolved(
                        "boundary parameter table does not reach its tolerance",
                    ));
                }
                let mut table = ParameterTable {
                    knots: Vec::new(),
                    angle: Vec::new(),
                    height: Vec::new(),
                    bound: worst,
                };
                for (index, span) in spans.iter().enumerate() {
                    let h = (&span[1] - &span[0]) / q(pieces);
                    let last = index + 1 == spans.len();
                    for piece in 0..pieces + i64::from(last) {
                        dec.check_deadline()?;
                        let x = &span[0] + &(&h * &q(piece));
                        let s = Jet::linear(b(&x, bits), k.one.clone());
                        let [value, slope] = angle(&s, index)?;
                        table.bound = &table.bound + &(value.rad() + slope.rad());
                        table.angle.push([value.mid(), slope.mid()]);
                        if row == BoundaryRow::Rail {
                            let pair = [&nodes[index], &nodes[index + 1]];
                            let height = hermite(
                                &s,
                                [&pair[0].s, &pair[1].s],
                                [
                                    [&enclosed[index][HEIGHT][0], &enclosed[index][HEIGHT][1]],
                                    [
                                        &enclosed[index + 1][HEIGHT][0],
                                        &enclosed[index + 1][HEIGHT][1],
                                    ],
                                ],
                                bits,
                            )?;
                            table.height.push([height.0[0].mid(), height.0[1].mid()]);
                        }
                        table.knots.push(x);
                    }
                }
                return Ok(table);
            }
            pieces *= 2;
        }
    }
}

// d/ds 2 atan(f) = 2 f' / (1 + f^2), valid through the jet's third coefficient.
fn slope_jet(f: &Jet, k: &Constants) -> Result<Jet, Refusal> {
    let bits = k.one.precision();
    let zero = b(&Rational::zero(), bits);
    let derivative = Jet(std::array::from_fn(|i| {
        if i == 4 {
            zero.clone()
        } else {
            &f.0[i + 1] * &b(&q(i as i64 + 1), bits)
        }
    }));
    derivative
        .scale(&k.two)
        .div(&Jet::constant(k.one.clone()).add(&f.square()))
}
