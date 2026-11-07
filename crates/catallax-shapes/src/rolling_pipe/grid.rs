//! A uniform display grid of the band with a whole-cell chord certificate.
//!
//! Each grid cell is a sub-patch of a rational Bezier span P / W. Against the
//! bilinear map B through its corners, `S - B = (P - B W) / W`, and both terms
//! of that numerator are polynomials one degree higher in each direction. Its
//! Bernstein coefficients over the smallest weight bound the distance at
//! matching parameters, and the quad's split adds its twist. Outward-rounded
//! floats keep this certificate cheap enough to run on every cell.

use catallax_numbers::Interval;

use super::*;
use crate::math;

const ROWS: usize = 3;
const MAX_PASSES: usize = 10;

type Control = [Interval; 4];

fn at(i: usize, j: usize) -> usize {
    i * ROWS + j
}
fn fraction(top: usize, bottom: usize) -> Interval {
    Interval::from_rational(&Rational::from_ratio(top as i64, bottom as i64).expect("nonzero"))
}
fn one() -> Interval {
    Interval::point(1.0).expect("finite")
}
fn mix(a: &Control, b: &Control, ratio: Interval) -> Control {
    let keep = one() - ratio;
    std::array::from_fn(|k| a[k] * keep + b[k] * ratio)
}

// De Casteljau division of a net at `ratio` along one axis, keeping one side.
fn divide(net: &[Control], axis: usize, ratio: Interval, high: bool) -> Vec<Control> {
    let degree = net.len() / ROWS - 1;
    let (lines, length) = if axis == 0 {
        (ROWS, degree + 1)
    } else {
        (degree + 1, ROWS)
    };
    let place = |line: usize, k: usize| if axis == 0 { at(k, line) } else { at(line, k) };
    let mut result = net.to_vec();
    for line in 0..lines {
        let mut row: Vec<Control> = (0..length).map(|k| net[place(line, k)]).collect();
        let mut kept = vec![row[0]; length];
        kept[if high { length - 1 } else { 0 }] = row[if high { length - 1 } else { 0 }];
        for step in 1..length {
            for k in 0..length - step {
                row[k] = mix(&row[k], &row[k + 1], ratio);
            }
            let slot = if high { length - 1 - step } else { step };
            kept[slot] = row[if high { length - 1 - step } else { 0 }];
        }
        for (k, value) in kept.into_iter().enumerate() {
            result[place(line, k)] = value;
        }
    }
    result
}

// Cut ratios for the pieces between neighbouring lines of `count`, computed once.
fn cuts(count: usize) -> Vec<[Interval; 2]> {
    (0..count)
        .map(|index| [fraction(index, count), fraction(1, count - index)])
        .collect()
}

// The piece of a net between grid lines `index` and `index + 1` of `count`.
fn piece(net: &[Control], axis: usize, index: usize, cuts: &[[Interval; 2]]) -> Vec<Control> {
    let mut part = net.to_vec();
    if index > 0 {
        part = divide(&part, axis, cuts[index][0], true);
    }
    if index + 1 < cuts.len() {
        part = divide(&part, axis, cuts[index][1], false);
    }
    part
}

fn point(h: &Control) -> Option<[Interval; 3]> {
    if h[3].lo() <= 0.0 {
        return None;
    }
    let inverse = one() / h[3];
    Some(std::array::from_fn(|k| h[k] * inverse))
}

// The Euclidean length's upper bound.
fn size(value: &[Interval; 3]) -> Option<f64> {
    let sum = value[0].square() + value[1].square() + value[2].square();
    Some(sum.sqrt()?.hi())
}

// The lifting factors of a chord numerator, once per grid.
struct Lift {
    degree: usize,
    factors: Vec<Interval>,
    quarter: Interval,
}
impl Lift {
    fn new(degree: usize) -> Self {
        // (1 - t) and t raise a Bernstein coefficient i of degree n to i and i + 1.
        let lifted = |degree: usize, step: usize, i: usize| {
            if step == 0 { degree + 1 - i } else { i + 1 }
        };
        let mut factors = Vec::with_capacity(4 * (degree + 1) * ROWS);
        for step_u in 0..2 {
            for step_v in 0..2 {
                for i in 0..=degree {
                    for j in 0..ROWS {
                        factors.push(fraction(
                            lifted(degree, step_u, i) * lifted(ROWS - 1, step_v, j),
                            (degree + 1) * ROWS,
                        ));
                    }
                }
            }
        }
        Self {
            degree,
            factors,
            quarter: fraction(1, 4),
        }
    }
    fn at(&self, step_u: usize, step_v: usize, i: usize, j: usize) -> Interval {
        self.factors[((step_u * 2 + step_v) * (self.degree + 1) + i) * ROWS + j]
    }
}

// The upper bound of the cell's chord, or `None` off a positive weight.
fn chord(net: &[Control], lift: &Lift) -> Option<f64> {
    let degree = lift.degree;
    let corners = [
        [point(&net[at(0, 0)])?, point(&net[at(0, ROWS - 1)])?],
        [
            point(&net[at(degree, 0)])?,
            point(&net[at(degree, ROWS - 1)])?,
        ],
    ];
    let mut numerator = 0.0_f64;
    for k in 0..=degree + 1 {
        for l in 0..=ROWS {
            let mut sum = [Interval::point(0.0)?; 3];
            for (step_u, along_v) in corners.iter().enumerate() {
                let Some(i) = k.checked_sub(step_u).filter(|i| *i <= degree) else {
                    continue;
                };
                for (step_v, corner) in along_v.iter().enumerate() {
                    let Some(j) = l.checked_sub(step_v).filter(|j| *j < ROWS) else {
                        continue;
                    };
                    let h = &net[at(i, j)];
                    let factor = lift.at(step_u, step_v, i, j);
                    for (m, value) in sum.iter_mut().enumerate() {
                        *value = *value + (h[m] - corner[m] * h[3]) * factor;
                    }
                }
            }
            numerator = numerator.max(size(&sum)?);
        }
    }
    let low = net.iter().map(|h| h[3].lo()).fold(f64::MAX, f64::min);
    if low <= 0.0 {
        return None;
    }
    let [[a, d], [c, e]] = corners;
    let twist: [Interval; 3] = std::array::from_fn(|k| (a[k] - c[k] - d[k] + e[k]) * lift.quarter);
    let bound = Interval::point(numerator)? / Interval::point(low)?;
    Some((bound + Interval::point(size(&twist)?)?).hi())
}

pub(super) fn fit(
    pipe: &RollingPipe,
    tolerance: &Rational,
    max_points: usize,
    dec: &Decider,
) -> Result<EnvelopeGrid, Refusal> {
    let failed = || invalid("a pipe grid cell has no chord certificate");
    let refused = || {
        Refusal::new(
            "mesh.resource-limit",
            "pipe grid cannot reach the chord tolerance within its limits",
        )
    };
    let carrier = pipe.carrier();
    let mut net = Vec::with_capacity(carrier.controls().len());
    for (control, weight) in carrier.controls().iter().zip(carrier.weights()) {
        let weight = dec.enclose_at(weight, 128)?;
        let [x, y, z] = math::point(control, 128, dec)?;
        net.push([&x * &weight, &y * &weight, &z * &weight, weight].map(|v| v.to_interval()));
    }
    let degree = carrier.degrees()[0];
    let spans: Vec<&[Control]> = (0..4)
        .map(|span| &net[span * degree * ROWS..(span * degree + degree + 1) * ROWS])
        .collect();
    let target = tolerance.to_f64_nearest() * 0.98;
    let lift = Lift::new(degree);
    let (mut per_span, mut lines) = (2usize, 2usize);
    let mut settled = None;
    for _ in 0..MAX_PASSES {
        dec.check_deadline()?;
        let mut worst = 0.0_f64;
        let mut worst_cell = None;
        let (along_u, along_v) = (cuts(per_span), cuts(lines));
        for span in &spans {
            for column in 0..per_span {
                let strip = piece(span, 0, column, &along_u);
                for line in 0..lines {
                    let cell = piece(&strip, 1, line, &along_v);
                    let bound = chord(&cell, &lift).ok_or_else(failed)?;
                    if bound > worst {
                        worst = bound;
                        worst_cell = Some(cell);
                    }
                }
            }
        }
        if worst <= target {
            settled = Some(worst);
            break;
        }
        // A cut along one axis takes its share of the bound to a quarter.
        let cell = worst_cell.ok_or_else(failed)?;
        let halved = |axis: usize| -> Result<f64, Refusal> {
            let half = fraction(1, 2);
            let low = chord(&divide(&cell, axis, half, false), &lift).ok_or_else(failed)?;
            let high = chord(&divide(&cell, axis, half, true), &lift).ok_or_else(failed)?;
            Ok(low.max(high))
        };
        let share = |after: f64| ((worst - after) * 4.0 / 3.0).max(0.0);
        let [along, across] = [share(halved(0)?), share(halved(1)?)];
        let scale = |count: usize, part: f64| -> usize {
            let factor = (2.0 * part / target).sqrt().max(1.0);
            (count as f64 * factor).ceil() as usize
        };
        let (next_span, next_lines) = (scale(per_span, along), scale(lines, across));
        if next_span == per_span && next_lines == lines {
            per_span += 1;
            lines += 1;
        } else {
            per_span = next_span;
            lines = next_lines;
        }
        if (4 * per_span + 1).saturating_mul(lines + 1) > max_points {
            return Err(refused());
        }
    }
    let chord = settled.ok_or_else(refused)?;
    let chord = Rational::from_f64_exact(chord).ok_or_else(failed)?;
    let columns = 4 * per_span + 1;
    let rows = lines + 1;
    let mut points = Vec::with_capacity(columns * rows);
    for column in 0..columns {
        let u = Rational::from_ratio(column as i64, per_span as i64).expect("nonzero");
        for row in 0..rows {
            let v = Rational::from_ratio(row as i64, lines as i64).expect("nonzero");
            points.push(math::point(&pipe.point(&u, &v), 128, dec)?);
        }
    }
    Ok(EnvelopeGrid {
        columns,
        rows,
        points,
        chord,
    })
}
