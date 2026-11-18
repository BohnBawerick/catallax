//! Evaluate one homogeneous cubic row once for every section and native rail.

use super::{Packed, ROWS, bad};
use catallax_numbers::{Ball, Rational, Refusal};

const BITS: u32 = 128;

pub(super) struct Basis {
    pub shares: [Rational; 3],
    pub section_shares: [Rational; 5],
    longitudinal: [[Ball; 4]; 3],
    transverse: [[Ball; 3]; 3],
}

impl Basis {
    pub fn new() -> Self {
        Self {
            shares: std::array::from_fn(|i| Rational::from_ratio(i as i64 + 1, 4).unwrap()),
            section_shares: std::array::from_fn(|i| Rational::from_ratio(i as i64, 4).unwrap()),
            longitudinal: [[27, 27, 9, 1], [8, 24, 24, 8], [1, 9, 27, 27]].map(|row| {
                row.map(|value| {
                    Ball::from_rational(&Rational::from_ratio(value, 64).unwrap(), BITS)
                })
            }),
            transverse: [[9, 6, 1], [4, 8, 4], [1, 6, 9]].map(|row| {
                row.map(|value| {
                    Ball::from_rational(&Rational::from_ratio(value, 16).unwrap(), BITS)
                })
            }),
        }
    }
}

pub(super) struct Span {
    controls: [[Ball; ROWS]; 4],
}

pub(super) struct Sample {
    pub points: [[Ball; 3]; 5],
    pub parameters: [[Rational; 2]; 2],
}

fn cartesian(row: &[Ball; 4]) -> Result<[Ball; 3], Refusal> {
    let coordinate = |i: usize| {
        row[i]
            .checked_div(&row[3])
            .ok_or_else(|| bad("fitted section homogeneous denominator cannot be enclosed"))
    };
    Ok([coordinate(0)?, coordinate(1)?, coordinate(2)?])
}

impl Span {
    pub fn new(controls: &[Packed; 4]) -> Self {
        Self {
            controls: controls
                .each_ref()
                .map(|row| row.each_ref().map(|value| Ball::from_rational(value, BITS))),
        }
    }

    pub fn at(&self, quarter: usize, basis: &Basis) -> Result<Sample, Refusal> {
        let weights = &basis.longitudinal[quarter];
        let packed: [Ball; ROWS] = std::array::from_fn(|axis| {
            (1..4).fold(&self.controls[0][axis] * &weights[0], |sum, i| {
                sum + &self.controls[i][axis] * &weights[i]
            })
        });
        let [
            x0,
            y0,
            z0,
            w0,
            x1,
            y1,
            z1,
            w1,
            x2,
            y2,
            z2,
            w2,
            u0,
            v0,
            u1,
            v1,
        ] = packed;
        let rows = [[x0, y0, z0, w0], [x1, y1, z1, w1], [x2, y2, z2, w2]];
        let mut points: [Option<[Ball; 3]>; 5] = [const { None }; 5];
        points[0] = Some(cartesian(&rows[0])?);
        points[4] = Some(cartesian(&rows[2])?);
        for (i, weights) in basis.transverse.iter().enumerate() {
            let row = std::array::from_fn(|axis| {
                let first = &rows[0][axis] * &weights[0];
                let second = &rows[1][axis] * &weights[1];
                first + second + &rows[2][axis] * &weights[2]
            });
            points[i + 1] = Some(cartesian(&row)?);
        }
        Ok(Sample {
            points: points.map(|point| point.expect("every section sample was assigned")),
            parameters: [[u0.mid(), v0.mid()], [u1.mid(), v1.mid()]],
        })
    }
}
