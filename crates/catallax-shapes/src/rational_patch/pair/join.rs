//! Joins certified arcs at their shared border points into branches.

use super::*;
use face::{Endpoint, isolate, overlap};
use std::borrow::Cow;
use std::collections::BTreeMap;

/// One accepted unit holding one arc.
pub(super) struct Resolved {
    pub(super) first: Arc<Cell>,
    pub(super) second_span: usize,
    pub(super) hull: Domain,
    pub(super) certificate: PairCertificate,
    pub(super) ends: [Endpoint; 2],
}

fn failure(detail: impl Into<String>) -> Refusal {
    Refusal::new("intersect.spline-pair-join", detail.into())
}

fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// Re-index only an exactly fixed, proved periodic row. Free root variables
/// and their equations do not change: C0 gluing identifies the same face map.
fn canonical_hit<'a>(ctx: &Context, hit: &'a face::Hit) -> Cow<'a, face::Hit> {
    let mut moved: Option<face::Hit> = None;
    for (index, (axis, value)) in hit.fixed.iter().enumerate() {
        let canonical = ctx.canonical(*axis, value, hit.spans[axis / 2]);
        if canonical != value {
            let target = moved.get_or_insert_with(|| hit.clone());
            target.fixed[index].1 = canonical.clone();
            target.bounds[*axis] = [canonical.clone(), canonical.clone()];
            target.face[*axis] = [canonical.clone(), canonical.clone()];
        }
    }
    moved.map_or(Cow::Borrowed(hit), Cow::Owned)
}

/// Whether two arc ends found from neighbouring cell pairs are the same
/// point: their boxes overlap and the face system has exactly one root in
/// the hull of both.
fn same_point(
    ctx: &Context,
    a: &Endpoint,
    b: &Endpoint,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let a = canonical_hit(ctx, &a.hit);
    let b = canonical_hit(ctx, &b.hit);
    if !overlap(&a.bounds, &b.bounds, dec)? {
        return Ok(false);
    }
    if a.fixed != b.fixed {
        return Err(failure(
            "two overlapping border points were isolated on different strata",
        ));
    }
    // Both faces have the same point map: across a knot by continuity, or
    // across an identified periodic row by the atlas's exact gluing proof.
    if inside(&a.bounds, &b.face, dec)? || inside(&b.bounds, &a.face, dec)? {
        return Ok(true);
    }
    let free = &a.free;
    let mut lo = Vec::with_capacity(free.len());
    let mut hi = Vec::with_capacity(free.len());
    for &axis in free {
        let [x, y] = if axis < 4 {
            [&a.bounds[axis], &b.bounds[axis]]
        } else {
            [&a.lengths[axis - 4], &b.lengths[axis - 4]]
        };
        lo.push(x[0].clone().min(y[0].clone()));
        hi.push(x[1].clone().max(y[1].clone()));
    }
    match isolate(&a.equations, &lo, &hi, dec)? {
        Some(Some(_)) => Ok(true),
        Some(None) => Err(failure("a certified border point vanished from its hull")),
        None => Err(failure(
            "two overlapping border points could not be proved equal or distinct",
        )),
    }
}

fn inside(
    inner: &[[Rational; 2]; 4],
    outer: &[[Rational; 2]; 4],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for axis in 0..4 {
        if net::less(&inner[axis][0], &outer[axis][0], dec)?
            || net::less(&outer[axis][1], &inner[axis][1], dec)?
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn point(ctx: &Context, end: &Endpoint) -> PairPoint {
    let span = &ctx.sides[0].spans[end.hit.spans[0]];
    PairPoint {
        parameters: end.hit.parameters(),
        enclosure: end.hit.bounds.clone(),
        first_chart: span.chart,
        first_span: span.local_span,
    }
}

/// Every branch, from the accepted arcs in search order.
pub(super) fn branches(
    ctx: &Context,
    arcs: Vec<Resolved>,
    dec: &mut Decider,
) -> Result<Vec<PatchPairBranch>, Refusal> {
    let ends: Vec<&Endpoint> = arcs.iter().flat_map(|arc| arc.ends.iter()).collect();
    let mut parent: Vec<usize> = (0..ends.len()).collect();
    let mut planes: BTreeMap<(usize, Rational), Vec<usize>> = BTreeMap::new();
    for (index, end) in ends.iter().enumerate() {
        for (axis, value) in &end.active {
            planes
                .entry((
                    *axis,
                    ctx.canonical(*axis, value, end.hit.spans[axis / 2]).clone(),
                ))
                .or_default()
                .push(index);
        }
    }
    for members in planes.values() {
        for (i, &a) in members.iter().enumerate() {
            for &b in &members[i + 1..] {
                if find(&mut parent, a) == find(&mut parent, b) {
                    continue;
                }
                if same_point(ctx, ends[a], ends[b], dec)? {
                    let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for index in 0..ends.len() {
        let root = find(&mut parent, index);
        groups.entry(root).or_default().push(index);
    }
    // next[arc] is the arc that starts where this one ends.
    let mut next = vec![None; arcs.len()];
    let mut open_start = vec![false; arcs.len()];
    for members in groups.values() {
        let physical = members.iter().any(|&m| {
            ends[m]
                .active
                .iter()
                .any(|(axis, value)| ctx.physical(*axis, value, ends[m].hit.spans[axis / 2]))
        });
        match members.as_slice() {
            [only] if physical => {
                if ends[*only].start {
                    open_start[only / 2] = true;
                }
            }
            [x, y] if !physical && ends[*x].start != ends[*y].start => {
                let (from, to) = if ends[*x].start { (*y, *x) } else { (*x, *y) };
                next[from / 2] = Some(to / 2);
            }
            _ => {
                return Err(failure(format!(
                    "a border point joins {} arc ends{}",
                    members.len(),
                    if physical { " on a patch border" } else { "" }
                )));
            }
        }
    }
    let mut used = vec![false; arcs.len()];
    let mut out = Vec::new();
    let follow = |first: usize, used: &mut Vec<bool>| -> Result<(Vec<usize>, bool), Refusal> {
        let mut order = Vec::new();
        let mut here = first;
        loop {
            if used[here] {
                return Ok((order, here == first));
            }
            used[here] = true;
            order.push(here);
            match next[here] {
                Some(after) => here = after,
                None => return Ok((order, false)),
            }
        }
    };
    let starts: Vec<usize> = (0..arcs.len()).filter(|&i| open_start[i]).collect();
    let mut chains = Vec::new();
    for first in starts {
        let (order, closed) = follow(first, &mut used)?;
        if closed {
            return Err(failure("an open branch closed on itself"));
        }
        chains.push((order, false));
    }
    for first in 0..arcs.len() {
        if used[first] {
            continue;
        }
        let (order, closed) = follow(first, &mut used)?;
        if !closed {
            return Err(failure("a branch without a border end does not close"));
        }
        chains.push((order, true));
    }
    for (order, closed) in chains {
        let mut list = Vec::with_capacity(order.len());
        for index in order {
            let arc = &arcs[index];
            list.push(PairArc {
                first_chart: ctx.sides[0].spans[arc.first.span].chart,
                first_span: ctx.sides[0].spans[arc.first.span].local_span,
                first_cell: arc.first.rect.clone(),
                second_chart: ctx.sides[1].spans[arc.second_span].chart,
                second_span: ctx.sides[1].spans[arc.second_span].local_span,
                second_cell: arc.hull.clone(),
                certificate: arc.certificate,
                ends: [point(ctx, &arc.ends[0]), point(ctx, &arc.ends[1])],
            });
        }
        out.push(PatchPairBranch {
            arcs: list.into(),
            closed,
        });
    }
    Ok(out)
}
