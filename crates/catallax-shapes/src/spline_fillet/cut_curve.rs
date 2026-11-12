//! A plane's cut written as the carrier composed with a chart curve.
//!
//! The cut is the graph `tau(theta)` of `g = n . H - c W = 0`, where `H` and
//! `W` are the carrier's homogeneous coordinates. The chart curve here is
//! linear in `theta` between breakpoints. Where it passes a node of `tau`, a
//! breakpoint is set exactly at the node, so each piece lies in one span.
//! Composing a span's tensor Bernstein form with a line of its chart gives a
//! rational Bezier curve of degree `D + 2` in `theta`. Its controls are
//! enclosed in balls, and the written ones are their middles.
//!
//! On a piece, the composed `g` is bounded by `M`, the largest of its own
//! Bernstein coefficients. On every span the cut crosses, the differences of
//! `g`'s coefficients along `tau` share one strict sign, so `|dg/dtau|` is at
//! least `m`, the least of `D / h` times the smallest such difference. So the
//! cut's point with the same `theta` lies within `M / m` in `tau`, and in
//! space within `M / m` times the carrier's speed along `tau` there.

use super::carrier::SpanRows;
use super::cut::coefficients;
use super::*;
use crate::math;

/// The first grid of `theta`.
const FIRST_PIECES: i64 = 16;
/// Halvings of a piece beyond the first grid. A piece still above the
/// tolerance after them is kept, and its bound reported, as the exchange
/// writer's own curve fits do.
const MAX_HALVINGS: u32 = 6;
/// Subdivisions of a piece's equation for its bound.
const BOUND_SPLITS: u32 = 2;

/// A plane's cut written as the carrier composed with a chart curve.
#[derive(Clone, Debug)]
pub struct CutCurve {
    /// Breakpoints of `theta`, increasing from 0 to 1.
    pub theta: Vec<Rational>,
    /// The chart curve's `tau` at each breakpoint; it is linear between them.
    pub tau: Vec<Rational>,
    /// The degree of each piece: the carrier's along `tau`, plus 2.
    pub degree: usize,
    /// The placed curve's controls and weights: each piece's Bezier controls,
    /// consecutive pieces sharing their common end.
    pub controls: Vec<(Point3, Rational)>,
    /// An upper bound, in millimetres, on the distance from each point of the
    /// curve to the cut's point with the same `theta`.
    pub upper: Rational,
}

fn refused(detail: impl Into<String>) -> Refusal {
    Refusal::new(SPLINE_FILLET_CUT, detail)
}

/// Precision of the composed controls.
const BITS: u32 = 256;

fn ball(value: &Rational) -> Ball {
    Ball::from_rational(value, BITS)
}

/// One de Casteljau split of Bernstein coefficients on `[0, 1]` at `t`:
/// the coefficients on `[0, t]` and on `[t, 1]`.
fn split(coefficients: &[Ball], t: &Ball) -> (Vec<Ball>, Vec<Ball>) {
    let n = coefficients.len();
    let mut work = coefficients.to_vec();
    let mut left = Vec::with_capacity(n);
    let mut right = vec![Ball::from_rational(&Rational::zero(), BITS); n];
    left.push(work[0].clone());
    right[n - 1] = work[n - 1].clone();
    for step in 1..n {
        for j in 0..n - step {
            work[j] = &work[j] + &(t * &(&work[j + 1] - &work[j]));
        }
        left.push(work[0].clone());
        right[n - 1 - step] = work[n - 1 - step].clone();
    }
    (left, right)
}

/// Bernstein coefficients of the same polynomial on `t = a + x (b - a)`,
/// for `x` in `[0, 1]`; `a` may exceed `b`.
fn restrict(coefficients: &[Ball], a: &Rational, b: &Rational) -> Vec<Ball> {
    let (lo, hi, reversed) = if a <= b { (a, b, false) } else { (b, a, true) };
    let mut out = coefficients.to_vec();
    if hi < &Rational::one() {
        out = split(&out, &ball(hi)).0;
    }
    if lo.sign() == Sign::Positive {
        let at = lo.checked_div(hi).expect("hi exceeds lo");
        out = split(&out, &ball(&at)).1;
    }
    if reversed {
        out.reverse();
    }
    out
}

fn binomial(n: usize, k: usize) -> Rational {
    Rational::from_i64((0..k).fold(1i64, |c, i| c * (n - i) as i64 / (i as i64 + 1)))
}

/// Enclosures of the homogeneous Bezier controls of one span's rows
/// composed with the chart line from `(s[0], theta[0])` to `(s[1],
/// theta[1])`, in the span's own `s` on `[0, 1]`.
fn compose(rows: &SpanRows, s: [&Rational; 2], theta: [&Rational; 2]) -> Vec<[Ball; 4]> {
    let degree = rows[0][0].len() - 1;
    let along: Vec<[Vec<Ball>; 4]> = rows
        .iter()
        .map(|row| {
            std::array::from_fn(|c| {
                let values: Vec<Ball> = row[c].iter().map(ball).collect();
                restrict(&values, s[0], s[1])
            })
        })
        .collect();
    // Across: each column of three row values, quadratic in `theta`.
    let across: Vec<[Vec<Ball>; 4]> = (0..=degree)
        .map(|i| {
            std::array::from_fn(|c| {
                let quadratic = [
                    along[0][c][i].clone(),
                    along[1][c][i].clone(),
                    along[2][c][i].clone(),
                ];
                restrict(&quadratic, theta[0], theta[1])
            })
        })
        .collect();
    // The diagonal of the tensor Bezier form.
    let total = degree + 2;
    (0..=total)
        .map(|m| {
            std::array::from_fn(|c| {
                let mut sum = ball(&Rational::zero());
                for (i, column) in across.iter().enumerate() {
                    let Some(k) = m.checked_sub(i).filter(|k| *k <= 2) else {
                        continue;
                    };
                    let factor = (binomial(degree, i) * binomial(2, k))
                        .checked_div(&binomial(total, m))
                        .expect("nonzero binomial");
                    sum = &sum + &(&ball(&factor) * &column[c][k]);
                }
                sum
            })
        })
        .collect()
}

/// An upper bound of `|f|` on `[0, 1]` for Bernstein coefficients `f`, from
/// the hulls of its pieces after `depth` halvings.
fn most(coefficients: &[Ball], depth: u32) -> Rational {
    if depth == 0 {
        return coefficients
            .iter()
            .map(|value| value.lo().abs().max(value.hi().abs()))
            .max()
            .unwrap_or_else(Rational::zero);
    }
    let (left, right) = split(
        coefficients,
        &ball(&Rational::from_ratio(1, 2).expect("nonzero")),
    );
    most(&left, depth - 1).max(most(&right, depth - 1))
}

/// An upper bound of the carrier's speed along `s` on the box `s` in
/// `[s0, s1]` (a part of one span, in its own `s`) and `theta` in
/// `[t0, t1]`. On the box the carrier is a rational Bezier patch with
/// homogeneous controls `X` and weights `W`. Along `s`, with `r` the point,
/// `r_s = (n / W) sum B_i (dX_i - r dW_i)` for each `theta`, where `dX_i` and
/// `dW_i` are convex combinations across of the controls' differences. So
/// `|r_s| <= n max |dX - r dW| / min W`, with `r` in the controls' hull.
fn speed(rows: &SpanRows, s: [&Rational; 2], theta: [&Rational; 2]) -> Option<Rational> {
    let degree = rows[0][0].len() - 1;
    // Restricted controls, across then along: `net[k][c][i]`.
    let net: Vec<[Vec<Ball>; 4]> = {
        let along: Vec<[Vec<Ball>; 4]> = rows
            .iter()
            .map(|row| {
                std::array::from_fn(|c| {
                    let values: Vec<Ball> = row[c].iter().map(ball).collect();
                    restrict(&values, s[0], s[1])
                })
            })
            .collect();
        let mut net: Vec<[Vec<Ball>; 4]> = (0..3)
            .map(|_| std::array::from_fn(|_| Vec::with_capacity(degree + 1)))
            .collect();
        for c in 0..4 {
            let columns = along[0][c].iter().zip(&along[1][c]).zip(&along[2][c]);
            for ((first, middle), last) in columns {
                let quadratic = [first.clone(), middle.clone(), last.clone()];
                for (k, value) in restrict(&quadratic, theta[0], theta[1])
                    .into_iter()
                    .enumerate()
                {
                    net[k][c].push(value);
                }
            }
        }
        net
    };
    let mut least: Option<Rational> = None;
    let mut hull: [Option<Ball>; 3] = [None, None, None];
    for row in &net {
        for (i, weight) in row[3].iter().enumerate() {
            let low = weight.lo();
            if low.sign() != Sign::Positive {
                return None;
            }
            if least.as_ref().is_none_or(|least| &low < least) {
                least = Some(low);
            }
            for (c, slot) in hull.iter_mut().enumerate() {
                let point = row[c][i].checked_div(weight)?;
                *slot = Some(match slot.take() {
                    Some(old) => old.hull(&point),
                    None => point,
                });
            }
        }
    }
    let hull = hull.map(|value| value.expect("controls"));
    let mut largest = Rational::zero();
    for row in &net {
        for i in 0..degree {
            let dw = &row[3][i + 1] - &row[3][i];
            let squared = (0..3).fold(Rational::zero(), |sum, c| {
                let value = &(&row[c][i + 1] - &row[c][i]) - &(&hull[c] * &dw);
                let magnitude = value.lo().abs().max(value.hi().abs());
                sum + &magnitude * &magnitude
            });
            largest = largest.max(squared);
        }
    }
    let scale = (s[1] - s[0]).abs();
    let bound =
        (&math::sqrt_upper(&largest) * &Rational::from_i64(degree as i64)).checked_div(&least?)?;
    bound.checked_div(&scale)
}

/// Binary64 Bernstein evaluation, for choosing breakpoints only.
fn bernstein(values: &[f64], t: f64) -> f64 {
    let mut work = values.to_vec();
    for step in 1..work.len() {
        for j in 0..work.len() - step {
            work[j] += t * (work[j + 1] - work[j]);
        }
    }
    work[0]
}

/// The cut's equation on the crossing spans, in binary64, for choosing
/// breakpoints. Every choice is then bounded exactly.
struct Guide<'a> {
    nodes: &'a [Node],
    first: usize,
    /// Per crossing span, `g`'s coefficients per row.
    g: Vec<[Vec<f64>; 3]>,
}

impl Guide<'_> {
    fn at(&self, span: usize, s: f64, theta: f64) -> f64 {
        let rows = &self.g[span - self.first];
        let values = [0, 1, 2].map(|k| bernstein(&rows[k], s));
        bernstein(&values, theta)
    }

    fn width(&self, span: usize) -> f64 {
        (&self.nodes[span + 1].tau - &self.nodes[span].tau).to_f64_nearest()
    }

    /// `tau` of the cut at `theta`, by bisection on the span where `g`
    /// changes sign.
    fn tau(&self, theta: f64) -> Rational {
        let last = self.first + self.g.len() - 1;
        let mut span = last;
        for candidate in self.first..=last {
            let ends = [
                self.at(candidate, 0.0, theta),
                self.at(candidate, 1.0, theta),
            ];
            if ends[0] == 0.0 || ends[0].signum() != ends[1].signum() {
                span = candidate;
                break;
            }
        }
        let start = self.at(span, 0.0, theta);
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            let value = self.at(span, mid, theta);
            if value == 0.0 {
                lo = mid;
                hi = mid;
                break;
            }
            if value.signum() == start.signum() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let s = 0.5 * (lo + hi);
        let tau = self.nodes[span].tau.to_f64_nearest() + s * self.width(span);
        Rational::from_f64_exact(tau).expect("finite")
    }

    /// `theta` in `(a, b)` where the cut meets `tau = nodes[node]`.
    fn theta_at_node(&self, node: usize, a: &Rational, b: &Rational) -> Rational {
        let (span, s) = if node > self.first {
            (node - 1, 1.0)
        } else {
            (node, 0.0)
        };
        let (fa, fb) = (a.to_f64_nearest(), b.to_f64_nearest());
        let start = self.at(span, s, fa);
        let (mut lo, mut hi) = (fa, fb);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            let value = self.at(span, s, mid);
            if value.signum() == start.signum() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let theta = Rational::from_f64_exact(0.5 * (lo + hi)).expect("finite");
        if &theta > a && &theta < b {
            theta
        } else {
            (a + b).mul_pow2(-1)
        }
    }
}

impl SplineFillet {
    /// The cut written as the carrier composed with a chart curve `tau(theta)`
    /// that is linear between breakpoints. Pieces are halved toward
    /// `tolerance` millimetres, at most [`MAX_HALVINGS`] times, and the curve
    /// reports the bound it reaches.
    pub fn cut_curve(
        &self,
        cut: &FilletCut,
        tolerance: &Rational,
        dec: &mut Decider,
    ) -> Result<CutCurve, Refusal> {
        let fit = &self.data.fit;
        let degree = fit.carrier.degrees()[0];
        let [first, last] = cut.spans;
        // The cut names the caller's plane; the rows live in the fit's frame.
        let (normal, offset) = surface::local_plane(&self.data.frame, &cut.normal, &cut.offset)
            .ok_or_else(|| refused("a face moved by an irrational frame has no rational cut"))?;
        // The least slope of `g` along `tau` on each crossing span.
        let mut slopes: Vec<Rational> = Vec::with_capacity(last - first + 1);
        let mut guide = Guide {
            nodes: &fit.nodes,
            first,
            g: Vec::with_capacity(last - first + 1),
        };
        for span in first..=last {
            dec.check_deadline()?;
            let g = coefficients(&fit.rows[span], &normal, &offset, dec)?;
            let width = &fit.nodes[span + 1].tau - &fit.nodes[span].tau;
            let mut slope: Option<Rational> = None;
            for row in &g {
                for pair in row.windows(2) {
                    let difference = (&pair[1] - &pair[0]).abs();
                    let bound = (&difference * &Rational::from_i64(degree as i64))
                        .checked_div(&width)
                        .expect("a span has width");
                    if slope.as_ref().is_none_or(|least| &bound < least) {
                        slope = Some(bound);
                    }
                }
            }
            slopes.push(
                slope
                    .filter(|least| least.sign() == Sign::Positive)
                    .ok_or_else(|| refused("the cut's equation has no proved slope along tau"))?,
            );
            guide
                .g
                .push(g.map(|row| row.iter().map(Rational::to_f64_nearest).collect()));
        }
        let slope = slopes.iter().min().expect("a crossing span").clone();
        // The least slope over the spans a box of `tau` touches.
        let slope_over = |lo: &Rational, hi: &Rational| {
            (first..=last)
                .filter(|&span| &fit.nodes[span + 1].tau >= lo && &fit.nodes[span].tau <= hi)
                .map(|span| slopes[span - first].clone())
                .min()
                .unwrap_or_else(|| slope.clone())
        };
        let range = [&fit.nodes[first].tau, &fit.nodes[last + 1].tau];
        // The carrier's speed along `tau` over a box of `tau` and `theta`,
        // from each crossing span's part of it.
        let speed_over =
            |taus: [&Rational; 2], thetas: [&Rational; 2]| -> Result<Rational, Refusal> {
                let mut most = Rational::zero();
                for span in first..=last {
                    let start = &fit.nodes[span].tau;
                    let end = &fit.nodes[span + 1].tau;
                    if end < taus[0] || start > taus[1] {
                        continue;
                    }
                    let width = end - start;
                    let local = |tau: &Rational| (tau - start).checked_div(&width).expect("width");
                    let lo = local(taus[0].max(start));
                    let hi = local(taus[1].min(end));
                    if lo >= hi {
                        continue;
                    }
                    let along = speed(&fit.rows[span], [&lo, &hi], thetas)
                        .ok_or_else(|| refused("a carrier weight is not proved positive"))?;
                    most = most.max(along.checked_div(&width).expect("width"));
                }
                Ok(most)
            };
        let span_of = |a: &Rational, b: &Rational| {
            let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
            (first..=last).find(|&s| &fit.nodes[s].tau <= lo && hi <= &fit.nodes[s + 1].tau)
        };
        // The first grid, with each end exactly at its rail's node.
        let mut points: Vec<(Rational, Rational)> = (0..=FIRST_PIECES)
            .map(|j| {
                let theta = Rational::from_ratio(j, FIRST_PIECES).expect("nonzero");
                let tau = match j {
                    0 => cut.ends[0].clone(),
                    j if j == FIRST_PIECES => cut.ends[1].clone(),
                    _ => guide.tau(theta.to_f64_nearest()),
                };
                (theta, tau)
            })
            .collect();
        let limbs = u64::from(BITS.div_ceil(64));
        let piece_cost = limbs * limbs * 8 * (degree as u64 + 3);
        // Room for the written controls' rounding.
        let target = tolerance * &Rational::from_ratio(3, 4).expect("nonzero");
        let mut accepted: Vec<Option<(Vec<[Ball; 4]>, Rational)>> = Vec::new();
        let mut upper = Rational::zero();
        for round in 0..=MAX_HALVINGS {
            // Breakpoints where the chart curve passes a node.
            let mut j = 0;
            while j + 1 < points.len() {
                let (a, b) = (&points[j], &points[j + 1]);
                if span_of(&a.1, &b.1).is_some() {
                    j += 1;
                    continue;
                }
                let node = if a.1 < b.1 {
                    (first + 1..=last).find(|&n| fit.nodes[n].tau > a.1)
                } else {
                    (first + 1..=last).rev().find(|&n| fit.nodes[n].tau < a.1)
                }
                .ok_or_else(|| refused("the chart curve leaves the crossing spans"))?;
                let theta = guide.theta_at_node(node, &a.0, &b.0);
                points.insert(j + 1, (theta, fit.nodes[node].tau.clone()));
                if accepted.len() > j {
                    accepted.insert(j + 1, None);
                    accepted[j] = None;
                }
            }
            accepted.resize(points.len() - 1, None);
            upper = Rational::zero();
            let mut halve = Vec::new();
            for (j, piece) in accepted.iter_mut().enumerate() {
                if let Some((_, bound)) = piece {
                    upper = upper.max(bound.clone());
                    continue;
                }
                // Each piece composes and bounds a few nets of balls, charged
                // as the curve fitter charges its spans: squared limbs per
                // control computed.
                dec.charge_work(piece_cost)?;
                let (a, b) = (&points[j], &points[j + 1]);
                let span = span_of(&a.1, &b.1).expect("breakpoints at nodes");
                let start = &fit.nodes[span].tau;
                let width = &fit.nodes[span + 1].tau - start;
                let local = |tau: &Rational| (tau - start).checked_div(&width).expect("width");
                let controls = compose(&fit.rows[span], [&local(&a.1), &local(&b.1)], [&a.0, &b.0]);
                let g: Vec<Ball> = controls
                    .iter()
                    .map(|control| {
                        let mut g = -(&ball(&offset) * &control[3]);
                        for (n, x) in normal.iter().zip(control) {
                            g = &g + &(&ball(n) * x);
                        }
                        g
                    })
                    .collect();
                let largest = most(&g, BOUND_SPLITS);
                // The cut's own `tau` is within `largest / slope` of the
                // chart curve's; the spans between them have at least `near`.
                let (lo, hi) = if a.1 <= b.1 {
                    (&a.1, &b.1)
                } else {
                    (&b.1, &a.1)
                };
                let reach = largest.checked_div(&slope).expect("positive slope");
                let near = slope_over(&(lo - &reach), &(hi + &reach));
                let shift = largest.checked_div(&near).expect("positive slope");
                let box_lo = (lo - &shift).max(range[0].clone());
                let box_hi = (hi + &shift).min(range[1].clone());
                let bound = &shift * &speed_over([&box_lo, &box_hi], [&a.0, &b.0])?;
                if bound > target && round < MAX_HALVINGS {
                    halve.push(j);
                } else {
                    upper = upper.max(bound.clone());
                    *piece = Some((controls, bound));
                }
            }
            if halve.is_empty() {
                break;
            }
            for &j in halve.iter().rev() {
                let theta = (&points[j].0 + &points[j + 1].0).mul_pow2(-1);
                let tau = guide.tau(theta.to_f64_nearest());
                points.insert(j + 1, (theta, tau));
                accepted.insert(j + 1, None);
            }
        }
        // One placed rational B-spline, consecutive pieces sharing ends:
        // each control is the middle of its enclosure, or of the hull of the
        // two that meet at a junction.
        let mut enclosures: Vec<[Ball; 4]> = Vec::new();
        for piece in &accepted {
            let (piece, _) = piece.as_ref().expect("every piece is bounded");
            for (m, control) in piece.iter().enumerate() {
                if m == 0
                    && let Some(end) = enclosures.last_mut()
                {
                    if (0..4)
                        .any(|c| end[c].hi() < control[c].lo() || control[c].hi() < end[c].lo())
                    {
                        return Err(refused("consecutive pieces do not share their end"));
                    }
                    *end = std::array::from_fn(|c| end[c].hull(&control[c]));
                    continue;
                }
                enclosures.push(control.clone());
            }
        }
        // The ends are exactly the rails' points at their nodes.
        let rail = |side: usize, tau: &Rational| -> [Rational; 4] {
            let node = fit
                .nodes
                .iter()
                .position(|node| &node.tau == tau)
                .expect("a node");
            let (span, i) = if node <= last {
                (node, 0)
            } else {
                (node - 1, degree)
            };
            std::array::from_fn(|c| fit.rows[span][2 * side][c][i].clone())
        };
        let ends = [rail(0, &cut.ends[0]), rail(1, &cut.ends[1])];
        let count = enclosures.len();
        let mut controls: Vec<(Point3, Rational)> = Vec::with_capacity(count);
        let mut position = Rational::zero();
        let mut relative = Rational::zero();
        let mut lows = [None::<Rational>, None, None];
        let mut highs = [None::<Rational>, None, None];
        for (index, enclosure) in enclosures.iter().enumerate() {
            let exact = match index {
                0 => Some(&ends[0]),
                i if i + 1 == count => Some(&ends[1]),
                _ => None,
            };
            let (point, weight) = if let Some(exact) = exact {
                let weight = exact[3].clone();
                let xyz = [0, 1, 2].map(|c| exact[c].checked_div(&weight).expect("weight"));
                (xyz, weight)
            } else {
                let weight = enclosure[3].clone();
                let low = weight.lo();
                if low.sign() != Sign::Positive {
                    return Err(refused("a composed weight is not proved positive"));
                }
                let cartesian: [Ball; 3] = std::array::from_fn(|c| {
                    enclosure[c].checked_div(&weight).expect("positive weight")
                });
                position = position.max(
                    cartesian
                        .iter()
                        .fold(Rational::zero(), |sum, value| sum + value.rad().mul_pow2(1)),
                );
                relative = relative.max(
                    weight
                        .rad()
                        .mul_pow2(1)
                        .checked_div(&low)
                        .expect("positive"),
                );
                (cartesian.map(|value| value.mid()), weight.mid())
            };
            for c in 0..3 {
                if lows[c].as_ref().is_none_or(|low| &point[c] < low) {
                    lows[c] = Some(point[c].clone());
                }
                if highs[c].as_ref().is_none_or(|high| &point[c] > high) {
                    highs[c] = Some(point[c].clone());
                }
            }
            let [x, y, z] = point;
            let placed = self
                .data
                .frame
                .apply_point(&Point3::new(x.into(), y.into(), z.into()));
            controls.push((placed, weight));
        }
        // Each written control is within twice its enclosure's radius of the
        // exact one, and a relative weight change `e` moves a convex
        // combination by at most `4 e` times the hull's size.
        let diameter = (0..3).fold(Rational::zero(), |sum, c| {
            sum + (highs[c].clone().expect("controls") - lows[c].clone().expect("controls"))
        });
        let upper = &(&upper + &position) + &(&relative.mul_pow2(2) * &diameter);
        Ok(CutCurve {
            theta: points.iter().map(|(theta, _)| theta.clone()).collect(),
            tau: points.into_iter().map(|(_, tau)| tau).collect(),
            degree: degree + 2,
            controls,
            upper,
        })
    }
}
