//! Parameter tables of the two rails on their own cylinders.
//!
//! A rail keeps the carrier's chart parameter, so its cylinder sees an angle of
//! the form `2 atan(f)`, which no polynomial represents exactly. The rail's
//! height along its cylinder is a fitted cubic, so a table holds it exactly.
//! The narrow rail's `f` is the chart's `s`, and the wide rail's is the fitted
//! half-angle tangent. Each table is a piecewise cubic Hermite fit whose bound
//! comes from an interval jet of the angle's slope.

use super::*;
use crate::ParameterTable;

// d/ds 2 atan(f) = 2 f' / (1 + f^2), valid through the jet's third coefficient.
fn turning(f: &Jet, k: &Constants) -> Result<Jet, Refusal> {
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

fn fitted(s: &Jet, pair: &[Node], row: usize, bits: u32) -> Result<Jet, Refusal> {
    let (_, ends) = span_chart(pair);
    let data = [0, 1].map(|end| {
        [
            b(&pair[end].values[row], bits),
            b(&pair[end].slopes[row], bits),
        ]
    });
    hermite(
        s,
        [&ends[0], &ends[1]],
        [[&data[0][0], &data[0][1]], [&data[1][0], &data[1][1]]],
        bits,
    )
}

impl RollingSaddle {
    /// Fits rail `rail`'s angle within `tolerance` radians: `0` is the wide rail
    /// and `1` the narrow one. The wide angle runs from the ring's `along`
    /// direction toward `across`. The narrow angle is the law's `t`, from
    /// `across` toward the side axis, and it runs from `-pi/2` to `3 pi/2`. The
    /// height is the rail's station along that cylinder's axis, from the law's
    /// origin on the wide axis and from `foot` on the narrow one. A rail the jet
    /// cannot resolve within 2^12 pieces per span refuses.
    pub fn parameter_table(
        &self,
        rail: usize,
        tolerance: &Rational,
        dec: &Decider,
    ) -> Result<ParameterTable, Refusal> {
        let bits = self.data.net.bits;
        let k = Constants::new(&self.data.definition, bits, dec)?;
        let nodes = &self.data.nodes;
        let pi = Ball::pi_checked(bits, dec)?;
        let argument = |s: &Jet, pair: &[Node]| -> Result<Jet, Refusal> {
            if rail == 0 {
                fitted(s, pair, HALF, bits)
            } else {
                Ok(Jet(s.0.clone()))
            }
        };
        let station = if rail == 0 { LIFT } else { HEIGHT };
        let mut pieces = 1i64;
        loop {
            dec.check_deadline()?;
            let mut worst = Rational::zero();
            for pair in nodes.windows(2) {
                let (_, ends) = span_chart(pair);
                let h = (&ends[1] - &ends[0]) / q(pieces);
                let h4 = h.pow(4).expect("nonzero span");
                for piece in 0..pieces {
                    dec.check_deadline()?;
                    let lo = &ends[0] + &(&h * &q(piece));
                    let hi = &lo + &h;
                    let s = Jet::linear(b(&lo, bits).hull(&b(&hi, bits)), k.one.clone());
                    let slope = turning(&argument(&s, pair)?, &k)?;
                    // g'''' = 6 c3 and the cubic remainder is h^4 |g''''| / 384.
                    worst = worst.max(magnitude(&slope.0[3]) * h4.clone() / q(64));
                }
            }
            if worst > *tolerance && pieces < 4096 {
                pieces *= 2;
                continue;
            }
            if worst > *tolerance {
                return Err(unresolved("rail parameter table misses its tolerance"));
            }
            let mut table = ParameterTable {
                knots: Vec::new(),
                angle: Vec::new(),
                height: Vec::new(),
                bound: worst,
            };
            let spans = nodes.len() - 1;
            for (index, pair) in nodes.windows(2).enumerate() {
                let (chart, ends) = span_chart(pair);
                let h = (&ends[1] - &ends[0]) / q(pieces);
                let turn = &pi * &b(&q(chart), bits);
                let last = index + 1 == spans;
                for piece in 0..pieces + i64::from(last) {
                    dec.check_deadline()?;
                    let x = &ends[0] + &(&h * &q(piece));
                    let s = Jet::linear(b(&x, bits), k.one.clone());
                    let f = argument(&s, pair)?;
                    let mut value = &k.two * &f.0[0].atan_checked(dec)?;
                    if rail == 1 {
                        value = value + turn.clone();
                    }
                    let slope = turning(&f, &k)?;
                    table.bound = &table.bound + &(value.rad() + slope.0[0].rad());
                    table.angle.push([value.mid(), slope.0[0].mid()]);
                    let height = fitted(&s, pair, station, bits)?;
                    table.height.push([height.0[0].mid(), height.0[1].mid()]);
                    table.knots.push(&x + &q(2 * chart));
                }
            }
            return Ok(table);
        }
    }
}
