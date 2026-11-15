//! Complete univariate pieces joined through owned source parameter ends.

use super::*;
use crate::rational_patch::{axis_section, exact_sign};
use catallax_numbers::{Real, RealPoly, Sign};

struct End {
    at: [Real; 2],
    physical: bool,
}

struct Piece {
    curve: BSplineCurve,
    ends: [usize; 2],
    crossing: bool,
}

fn same_at(a: &[Real; 2], b: &[Real; 2], dec: &mut Decider) -> Result<bool, Refusal> {
    for axis in 0..2 {
        if exact_sign(&(&a[axis] - &b[axis]), dec)? != Sign::Zero {
            return Ok(false);
        }
    }
    Ok(true)
}

fn endpoint(
    atlas: &RationalAtlas,
    chart: &AtlasChart,
    mut at: [Real; 2],
    ends: &mut Vec<End>,
    dec: &mut Decider,
) -> Result<usize, Refusal> {
    let mut physical = false;
    for (side, axis, end) in [
        (SplineBoundary::UMin, 0, 0),
        (SplineBoundary::UMax, 0, 1),
        (SplineBoundary::VMin, 1, 0),
        (SplineBoundary::VMax, 1, 1),
    ] {
        if chart.neighbor(side).is_none()
            && exact_sign(
                &(&at[axis] - &Real::from_rational(chart.patch().domain()[axis][end].clone())),
                dec,
            )? == Sign::Zero
        {
            physical = true;
        }
    }
    for (axis, value) in at.iter_mut().enumerate() {
        if atlas.has_periodic_join(axis)
            && exact_sign(
                &(&*value - &Real::from_rational(atlas.source().domain()[axis][1].clone())),
                dec,
            )? == Sign::Zero
        {
            *value = Real::from_rational(atlas.source().domain()[axis][0].clone());
        }
    }
    for (index, previous) in ends.iter_mut().enumerate() {
        if same_at(&previous.at, &at, dec)? {
            previous.physical |= physical;
            return Ok(index);
        }
    }
    ends.push(End { at, physical });
    Ok(ends.len() - 1)
}

// A duplicate seam row needs full homogeneous identity, not shared endpoints.
fn same_piece(
    a: &BSplineCurve,
    b: &BSplineCurve,
    forward: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut a_spans = a.spans();
    let mut b_spans = b.spans();
    let Some((a_span, _)) = a_spans.next() else {
        return Ok(false);
    };
    let Some((b_span, _)) = b_spans.next() else {
        return Ok(false);
    };
    if a_spans.next().is_some() || b_spans.next().is_some() {
        return Ok(false);
    }
    let a_poly = a.span_polynomials(a_span, dec)?;
    let b_poly = b.span_polynomials(b_span, dec)?;
    let [a_lo, a_hi] = a.domain();
    let [b_lo, b_hi] = b.domain();
    let scale = (b_hi - b_lo) / (a_hi - a_lo);
    let scale = if forward { scale } else { -scale };
    let origin = if forward { b_lo } else { b_hi };
    let parameter = RealPoly::var(1, 0).scale(&Real::from_rational(scale.clone()))
        + RealPoly::constant(1, Real::from_rational(origin - &(a_lo * &scale)));
    let b_poly = b_poly.map(|poly| poly.substitute(0, &parameter));
    for axis in 0..3 {
        if !(&a_poly[axis] * &b_poly[3] - &b_poly[axis] * &a_poly[3]).is_zero(dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn components(
    pieces: &[Piece],
    ends: &[End],
    dec: &mut Decider,
) -> Result<Option<Vec<(BSplineCurve, bool)>>, Refusal> {
    let mut incident = vec![Vec::new(); ends.len()];
    for (index, piece) in pieces.iter().enumerate() {
        for (side, &end) in piece.ends.iter().enumerate() {
            incident[end].push((index, side));
        }
    }
    for (end, incident) in ends.iter().zip(&incident) {
        if incident.len() > 2 || (incident.len() == 1 && !end.physical) {
            return Ok(None);
        }
    }
    let mut used = vec![false; pieces.len()];
    let mut result = Vec::new();
    while let Some(first) = used.iter().position(|used| !used) {
        let mut start = pieces[first].ends[0];
        let mut previous = first;
        let mut probe = pieces[first].ends[1];
        if incident[start].len() == 2 {
            let mut visited = vec![false; pieces.len()];
            visited[first] = true;
            loop {
                let Some(&(next, side)) = incident[probe]
                    .iter()
                    .find(|&&(index, _)| index != previous)
                else {
                    start = probe;
                    break;
                };
                if visited[next] {
                    break;
                }
                visited[next] = true;
                previous = next;
                probe = pieces[next].ends[1 - side];
            }
        }
        let crossing = pieces[first].crossing;
        let mut chain = Vec::new();
        let mut at = start;
        while let Some(&(index, side)) = incident[at].iter().find(|&&(index, _)| !used[index]) {
            let piece = &pieces[index];
            if piece.crossing != crossing {
                return Ok(None);
            }
            used[index] = true;
            chain.push((piece.curve.clone(), side == 0));
            at = piece.ends[1 - side];
        }
        // Multiple one-sided higher-order germs need a two-sided contact proof.
        if chain.len() > 1 && !crossing {
            return Ok(None);
        }
        let curve = if chain.len() == 1 && chain[0].1 {
            chain.pop().expect("one complete piece").0
        } else {
            BSplineCurve::joined(&chain, dec)?
        };
        if at != start {
            let [lo, hi] = curve.domain();
            let first = curve.point_at(&Real::from_rational(lo.clone()), dec)?;
            let last = curve.point_at(&Real::from_rational(hi.clone()), dec)?;
            if exact_sign(&first.minus(&last).norm_squared(), dec)? == Sign::Zero {
                return Ok(None);
            }
        }
        result.push((curve, crossing));
    }
    Ok(Some(result))
}

pub(super) fn section(
    atlas: &RationalAtlas,
    surface: &SectionSurface,
    join_rows: bool,
    dec: &mut Decider,
) -> Result<Option<Vec<(BSplineCurve, bool)>>, Refusal> {
    let mut pieces: Vec<Piece> = Vec::new();
    let mut ends = Vec::new();
    for chart in atlas.charts() {
        let Some(section) = axis_section(chart.patch(), surface, dec)? else {
            return Ok(None);
        };
        if !section.points().is_empty() || !section.branches().is_empty() {
            return Ok(None);
        }
        let curves: Vec<_> = section.curves().collect();
        if !join_rows && !curves.is_empty() && !pieces.is_empty() {
            return Ok(None);
        }
        for (curve, crossing) in curves {
            let mut nodes = [0; 2];
            for (side, bound) in curve.domain().into_iter().enumerate() {
                let point = curve.point_at(&Real::from_rational(bound.clone()), dec)?;
                let Some(at) = chart.patch().parameters_at(&point, dec)? else {
                    return Ok(None);
                };
                nodes[side] = endpoint(atlas, chart, at, &mut ends, dec)?;
                if !join_rows && !ends[nodes[side]].physical {
                    return Ok(None);
                }
            }
            let mut duplicate = false;
            for old in &pieces {
                let forward = old.ends == nodes;
                let reverse = old.ends == [nodes[1], nodes[0]];
                if (forward || reverse) && same_piece(&old.curve, curve, forward, dec)? {
                    if old.crossing != crossing {
                        return Ok(None);
                    }
                    duplicate = true;
                    break;
                }
            }
            if !duplicate {
                pieces.push(Piece {
                    curve: curve.clone(),
                    ends: nodes,
                    crossing,
                });
            }
        }
    }
    components(&pieces, &ends, dec)
}
