//! Where a plane cuts the face.
//!
//! On each span the face's homogeneous coordinates are
//! `H = sum_k B_k(theta) R_k(s)`, with the quadratic Bernstein basis `B_k`
//! and the span's exact rows `R_k`. For the plane `n . x = c`, the function
//! `g = n . H - c W` is then a tensor Bernstein polynomial whose
//! coefficients are those of the rows, and its sign is the side of the
//! plane, since `W > 0`. Exact signs of those coefficients prove:
//!
//! - every span before the cut has `g` of one sign, and every span after it
//!   the other;
//! - on every span in between, the differences of consecutive coefficients
//!   along `s` have the after sign, so `g` is strictly monotone along each
//!   curve of constant `theta`.
//!
//! So each such curve crosses the plane exactly once, and the cut is a graph
//! `tau(theta)` over the whole of `[0, 1]`. Its ends are the exact zeros of
//! `g` on the two rails, which must lie at fit nodes; a rail border there
//! makes them rational.

use super::carrier::SpanRows;
use super::*;

/// A plane's cut of the face: a graph `tau(theta)` over `theta` in `[0, 1]`
/// that crosses each curve of constant `theta` once.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FilletCut {
    /// The plane's normal.
    pub normal: [Rational; 3],
    /// The plane's offset: its points `x` have `normal . x = offset`.
    pub offset: Rational,
    /// The cut's ends, as `tau` on the first rail and on the second.
    pub ends: [Rational; 2],
    /// The first and last span the cut crosses.
    pub spans: [usize; 2],
    /// The side of the plane, the sign of `normal . x - offset`, of the face
    /// before the cut along `tau`. The face after it is on the other side.
    pub before: Sign,
}

fn refused(detail: impl Into<String>) -> Refusal {
    Refusal::new(SPLINE_FILLET_CUT, detail)
}

/// The one strict sign shared by every value, if any.
fn common<'a>(values: impl IntoIterator<Item = &'a Rational>) -> Option<Sign> {
    let mut out = None;
    for value in values {
        let sign = value.sign();
        if sign == Sign::Zero || out.is_some_and(|seen| seen != sign) {
            return None;
        }
        out = Some(sign);
    }
    out
}

fn opposite(sign: Sign) -> Sign {
    match sign {
        Sign::Positive => Sign::Negative,
        Sign::Negative => Sign::Positive,
        Sign::Zero => Sign::Zero,
    }
}

/// The coefficients of `g` on one span: row, then index along `s`.
fn coefficients(
    rows: &SpanRows,
    normal: &[Rational; 3],
    offset: &Rational,
    dec: &mut Decider,
) -> Result<[Vec<Rational>; 3], Refusal> {
    let mut out: [Vec<Rational>; 3] = std::array::from_fn(|_| Vec::new());
    for (row, g) in rows.iter().zip(&mut out) {
        for (j, weight) in row[3].iter().enumerate() {
            let mut value = Rational::zero();
            for (coordinate, n) in row.iter().zip(normal) {
                if n.is_zero() {
                    continue;
                }
                dec.reserve_exact_arithmetic(n, &coordinate[j], true)?;
                value = value + n * &coordinate[j];
            }
            dec.reserve_exact_arithmetic(offset, weight, true)?;
            g.push(value - offset * weight);
        }
    }
    Ok(out)
}

pub(super) fn plane(
    rows: &[SpanRows],
    nodes: &[Node],
    normal: [Rational; 3],
    offset: Rational,
    dec: &mut Decider,
) -> Result<Option<FilletCut>, Refusal> {
    if normal.iter().all(Rational::is_zero) {
        return Err(definition("a cutting plane needs a nonzero normal"));
    }
    let mut spans = Vec::with_capacity(rows.len());
    for span in rows {
        dec.check_deadline()?;
        let g = coefficients(span, &normal, &offset, dec)?;
        let sign = common(g.iter().flatten());
        let differences: Vec<Rational> = g
            .iter()
            .flat_map(|row| row.windows(2).map(|pair| &pair[1] - &pair[0]))
            .collect();
        spans.push((g, sign, common(&differences)));
    }
    let Some(first) = spans.iter().position(|span| span.1.is_none()) else {
        return Ok(None);
    };
    let last = spans
        .iter()
        .rposition(|span| span.1.is_none())
        .expect("one crossing span");
    let Some(after) = spans[first].2 else {
        return Err(refused(format!(
            "the face is not proved monotone across the plane on span {first}"
        )));
    };
    let before = opposite(after);
    for (index, span) in spans.iter().enumerate() {
        let fits = if index < first {
            span.1 == Some(before)
        } else if index > last {
            span.1 == Some(after)
        } else {
            span.2 == Some(after)
        };
        if !fits {
            return Err(refused(format!(
                "span {index} is not proved on its side of the plane, or monotone across it"
            )));
        }
    }
    // The cut must leave the face through its rails, not its ends.
    let column = |index: usize, j: usize| spans[index].0.iter().map(move |row| &row[j]);
    if first == 0 && common(column(0, 0)) != Some(before) {
        return Err(refused("the cut reaches the face's first end"));
    }
    let degree = spans[last].0[0].len() - 1;
    if last == spans.len() - 1 && common(column(last, degree)) != Some(after) {
        return Err(refused("the cut reaches the face's last end"));
    }
    // Each rail's zero, at a node of the crossing spans.
    let mut ends: [Option<Rational>; 2] = [None, None];
    for (end, row) in ends.iter_mut().zip([0, 2]) {
        for node in first..=last + 1 {
            let value = if node <= last {
                &spans[node].0[row][0]
            } else {
                &spans[last].0[row][degree]
            };
            if value.is_zero() {
                *end = Some(nodes[node].tau.clone());
            }
        }
    }
    let [Some(start), Some(end)] = ends else {
        return Err(refused("the cut does not meet both rails at fit nodes"));
    };
    Ok(Some(FilletCut {
        normal,
        offset,
        ends: [start, end],
        spans: [first, last],
        before,
    }))
}
