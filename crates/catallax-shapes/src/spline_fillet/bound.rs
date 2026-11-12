//! The recomputed two-sided position bound.
//!
//! Carrier and exact face are both rational quadratics in the section
//! parameter with rows `P_i` and positive weights `w_i`, so each is a convex
//! combination `sum l_i P_i`. With carrier rows `P~_i` and weights
//! `w~_i = w_i (1 + e_i)`, `|e_i| <= eps < 1`,
//!
//! ```text
//! |F - S| <= max_i |P~_i - P_i| + (2 eps / (1 - eps)) max_i |P_i - C|
//! ```
//!
//! where `C` is the ball's centre and `max_i |P_i - C| = r / cos(phi/2)`.
//! Each fitted row's error is its cubic Hermite remainder, from the fourth
//! Taylor coefficient over the span, plus the effect of rounded node data.
//! A rail moves by at most its support's parameter derivatives times its
//! parameters' errors, over a box holding both the exact and the fitted
//! parameters. Both surfaces share their parameters, so the bound holds from
//! each surface to the other.

use super::carrier::{fitted, node_tau};
use super::jets::{Meter, Solution, ball, magnitude};
use super::*;

/// A rational upper bound on the square root of a nonnegative rational.
pub(super) fn root_up(value: &Rational, dec: &mut Decider) -> Result<Rational, Refusal> {
    if value.sign() != Sign::Positive {
        return Ok(Rational::zero());
    }
    let guess = value.to_f64_nearest().sqrt() * (1.0 + 1.0 / f64::from(1u32 << 30));
    let mut candidate = Rational::from_f64_exact(guess)
        .filter(|_| guess.is_finite() && guess > 0.0)
        .unwrap_or_else(|| value.clone().max(Rational::one()));
    for _ in 0..64 {
        dec.reserve_exact_arithmetic(&candidate, &candidate, true)?;
        if &(&candidate * &candidate) >= value {
            return Ok(candidate);
        }
        candidate = candidate.mul_pow2(1);
    }
    Ok(value.clone().max(Rational::one()))
}

/// Upper bounds on the parameter derivatives of support `k`'s point and
/// weight over a parameter box: per parameter, the Euclidean size of the
/// point's partial and the size of the weight's partial.
fn slopes(
    support: &Support,
    at: &[Ball; 2],
    meter: &mut Meter,
    dec: &mut Decider,
) -> Result<[[Rational; 2]; 2], Refusal> {
    let mut out: [[Rational; 2]; 2] = std::array::from_fn(|_| [Rational::zero(), Rational::zero()]);
    for (axis, slot) in out.iter_mut().enumerate() {
        let jets: [Jet; 2] = std::array::from_fn(|i| {
            if i == axis {
                Jet::linear(at[i].clone(), ball(&Rational::one()), 1)
            } else {
                Jet::constant(at[i].clone(), 1)
            }
        });
        let contact = support.contact(&jets, meter)?;
        meter.charge(dec)?;
        let square = contact
            .point
            .iter()
            .map(|j| magnitude(&j.0[1]).pow(2).expect("square"))
            .fold(Rational::zero(), |s, v| s + v);
        *slot = [root_up(&square, dec)?, magnitude(&contact.weight.0[1])];
    }
    Ok(out)
}

fn widened(value: &Ball, by: &Rational) -> Ball {
    ball(&(value.lo() - by.clone())).hull(&ball(&(value.hi() + by.clone())))
}

/// One span's proof: its bound, and the second-order jets of the spine and
/// the centre over the span, which the clearance proof reuses.
pub(super) struct SpanProof {
    pub(super) upper: Rational,
    pub(super) x: [Jet; 3],
    pub(super) tau: Jet,
    pub(super) centre: [Jet; 3],
}

/// The bound over one span.
fn span_bound(
    exact: &Exact,
    pair: [&Node; 2],
    ends: [&NodeData; 2],
    estimates: [[f64; 3]; 2],
    dec: &mut Decider,
) -> Result<SpanProof, Refusal> {
    let span = [&pair[0].tau, &pair[1].tau];
    let h = span[1] - span[0];
    let solution: Solution = exact.solve(span, estimates, 4, dec)?;
    let mut meter = Meter::default();
    let rows = exact.rows(&solution, &mut meter)?;
    meter.charge(dec)?;
    let whole = fitted(&solution, &rows);
    let h4 = h.pow(4).expect("nonempty span");
    let sixteen = Rational::from_i64(16);
    let slope_factor = Rational::from_ratio(4, 27).expect("nonzero") * h.clone();
    let error: [Rational; ROWS] = std::array::from_fn(|j| {
        let remainder = magnitude(&whole[j].0[4]) * h4.clone() / sixteen.clone();
        let value = |k: usize| magnitude(&(&ball(&pair[k].values[j]) - &ends[k][j][0]));
        let slope = |k: usize| magnitude(&(&ball(&pair[k].slopes[j]) - &ends[k][j][1]));
        remainder + value(0).max(value(1)) + slope_factor.clone() * (slope(0) + slope(1))
    });
    // Rails: parameter errors through the supports' derivatives over boxes
    // that hold both the exact and the fitted parameters.
    let x: [Ball; 3] = std::array::from_fn(|i| widened(&solution.x[i].0[0], &error[i]));
    let tau = solution.tau.0[0].clone();
    let [pa, pb] = exact.params(&x, &tau);
    let [da, db] = [
        slopes(&exact.supports[0], &pa, &mut meter, dec)?,
        slopes(&exact.supports[1], &pb, &mut meter, dec)?,
    ];
    // The first support's free parameter is the unknown x0; tau is exact.
    let free = 1 - exact.axis;
    let first = da[free][0].clone() * error[0].clone();
    let second = db[0][0].clone() * error[1].clone() + db[1][0].clone() * error[2].clone();
    let middle = root_up(
        &(0..3).fold(Rational::zero(), |s, i| {
            s + error[MIDDLE + i].clone() * error[MIDDLE + i].clone()
        }),
        dec,
    )?;
    let positive = |value: &Ball, what: &str| -> Result<Rational, Refusal> {
        let lo = value.lo();
        if lo.sign() != Sign::Positive {
            return Err(spine(format!("{what} is not separated from zero")));
        }
        Ok(lo)
    };
    let weights = [
        positive(&rows.contacts[0].weight.0[0], "the first support's weight")?,
        positive(&rows.weight.0[0], "the middle weight")?,
        positive(&rows.contacts[1].weight.0[0], "the second support's weight")?,
    ];
    let relative = [
        da[free][1].clone() * error[0].clone() / weights[0].clone(),
        error[WEIGHT].clone() / weights[1].clone(),
        (db[0][1].clone() * error[1].clone() + db[1][1].clone() * error[2].clone())
            / weights[2].clone(),
    ];
    let eps = relative.into_iter().max().expect("three");
    if eps >= Rational::one() {
        return Err(Refusal::new(
            SPLINE_FILLET_DEVIATION,
            "a fitted weight's error is not below its exact value",
        ));
    }
    let half = positive(&rows.half.0[0], "the section's half-angle cosine")?;
    let reach = exact.radius.clone() / half;
    let spread = Rational::from_i64(2) * eps.clone() / (Rational::one() - eps) * reach;
    Ok(SpanProof {
        upper: first.max(middle).max(second) + spread,
        x: solution.x.each_ref().map(|jet| jet.truncated(2)),
        tau: solution.tau.truncated(1),
        centre: rows.centre.each_ref().map(|jet| jet.truncated(2)),
    })
}

/// The largest bound over every span of a retained fit.
pub(super) fn recompute(
    exact: &Exact,
    nodes: &[Node],
    estimates: &[[f64; 3]],
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let data = order::in_order(
        nodes.len(),
        |i, dec| Ok(carrier::node(exact, &nodes[i].tau, estimates[i], dec)?.1),
        dec,
    )?;
    let bounds = order::in_order(
        nodes.len() - 1,
        |i, dec| {
            span_bound(
                exact,
                [&nodes[i], &nodes[i + 1]],
                [&data[i], &data[i + 1]],
                [estimates[i], estimates[i + 1]],
                dec,
            )
        },
        dec,
    )?;
    Ok(bounds
        .into_iter()
        .map(|proof| proof.upper)
        .max()
        .expect("nonempty"))
}

/// A proved fit: its nodes, their binary64 estimates, each span's proof and
/// the largest bound.
pub(super) struct Fit {
    pub(super) nodes: Vec<Node>,
    pub(super) estimates: Vec<[f64; 3]>,
    pub(super) spans: Vec<SpanProof>,
    pub(super) upper: Rational,
}

/// A refusal a narrower span may avoid: a span solve that was not proved,
/// or a weight error that was not below its weight.
fn narrower_may_prove(refusal: &Refusal) -> bool {
    refusal.code == SPLINE_FILLET_SPINE || refusal.code == SPLINE_FILLET_DEVIATION
}

struct Entry {
    tau: Rational,
    estimate: [f64; 3],
    proved: Option<(Node, NodeData)>,
}

struct Piece {
    halvings: u32,
    outcome: Option<Result<SpanProof, Refusal>>,
}

/// Fits from eight uniform spans, halving each span that misses the target
/// or refuses in a way a narrower span may avoid, up to the limit of
/// halvings. A span's proof depends only on its two nodes, so a span
/// that passes is never proved again. New nodes and spans are proved in
/// order on every core.
pub(super) fn refine(
    exact: &Exact,
    range: &[Rational; 2],
    seed: (f64, [f64; 3]),
    limits: &Limits,
    dec: &mut Decider,
) -> Result<Fit, Refusal> {
    let first = carrier::estimates(exact, range, FIRST_SPANS, seed)?;
    let mut entries: Vec<Entry> = first
        .into_iter()
        .enumerate()
        .map(|(i, estimate)| Entry {
            tau: node_tau(range, FIRST_SPANS, i),
            estimate,
            proved: None,
        })
        .collect();
    let mut pieces: Vec<Piece> = (0..FIRST_SPANS)
        .map(|_| Piece {
            halvings: 0,
            outcome: None,
        })
        .collect();
    loop {
        let open: Vec<usize> = (0..entries.len())
            .filter(|&i| entries[i].proved.is_none())
            .collect();
        let proved = order::in_order(
            open.len(),
            |k, dec| {
                let entry = &entries[open[k]];
                carrier::node(exact, &entry.tau, entry.estimate, dec)
            },
            dec,
        )?;
        for (i, node) in open.into_iter().zip(proved) {
            entries[i].proved = Some(node);
        }
        let open: Vec<usize> = (0..pieces.len())
            .filter(|&i| pieces[i].outcome.is_none())
            .collect();
        let outcomes = order::in_order(
            open.len(),
            |k, dec| {
                let i = open[k];
                let [a, b] = [i, i + 1].map(|j| entries[j].proved.as_ref().expect("proved node"));
                match span_bound(
                    exact,
                    [&a.0, &b.0],
                    [&a.1, &b.1],
                    [entries[i].estimate, entries[i + 1].estimate],
                    dec,
                ) {
                    Err(refusal) if narrower_may_prove(&refusal) => Ok(Err(refusal)),
                    other => other.map(Ok),
                }
            },
            dec,
        )?;
        for (i, outcome) in open.into_iter().zip(outcomes) {
            pieces[i].outcome = Some(outcome);
        }
        // A span is halved when it misses the target or refused, while it
        // may be halved again. Its bound often falls much faster than the
        // fourth power of the width, as its enclosures tighten too, so one
        // halving at a time keeps the fewest spans.
        let halve: Vec<bool> = pieces
            .iter()
            .map(|piece| {
                let misses = match &piece.outcome {
                    Some(Ok(proof)) => proof.upper > limits.target,
                    _ => true,
                };
                misses && piece.halvings < limits.halvings
            })
            .collect();
        if !halve.contains(&true) {
            break;
        }
        let mut old = std::mem::take(&mut entries).into_iter();
        let mut left = old.next().expect("a first node");
        for (piece, halve) in std::mem::take(&mut pieces).into_iter().zip(halve) {
            let right = old.next().expect("a node after each span");
            if halve {
                let tau = (&left.tau + &right.tau) / Rational::from_i64(2);
                let estimate = carrier::between(exact, &tau, [left.estimate, right.estimate], 0.5)?;
                let halvings = piece.halvings + 1;
                entries.push(left);
                entries.push(Entry {
                    tau,
                    estimate,
                    proved: None,
                });
                for _ in 0..2 {
                    pieces.push(Piece {
                        halvings,
                        outcome: None,
                    });
                }
            } else {
                entries.push(left);
                pieces.push(piece);
            }
            left = right;
        }
        entries.push(left);
    }
    let mut spans = Vec::with_capacity(pieces.len());
    for (i, piece) in pieces.into_iter().enumerate() {
        match piece.outcome.expect("every span tried") {
            Ok(proof) => spans.push(proof),
            Err(refusal) => {
                return Err(Refusal::new(
                    refusal.code,
                    format!(
                        "on tau in [{}, {}] after {} halvings: {}",
                        entries[i].tau,
                        entries[i + 1].tau,
                        piece.halvings,
                        refusal.detail
                    ),
                ));
            }
        }
    }
    let upper = spans
        .iter()
        .map(|proof| proof.upper.clone())
        .max()
        .expect("nonempty");
    let (nodes, estimates) = entries
        .into_iter()
        .map(|entry| (entry.proved.expect("proved node").0, entry.estimate))
        .unzip();
    Ok(Fit {
        nodes,
        estimates,
        spans,
        upper,
    })
}
