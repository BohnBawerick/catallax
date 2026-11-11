//! Joins certified arcs at their shared border points into branches.

use super::*;
use face::{Endpoint, isolate, overlap};
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

/// Whether two arc ends found from neighbouring cell pairs are the same
/// point: their boxes overlap and the face system has exactly one root in
/// the hull of both.
fn same_point(a: &Endpoint, b: &Endpoint, dec: &mut Decider) -> Result<bool, Refusal> {
    if !overlap(&a.hit.bounds, &b.hit.bounds, dec)? {
        return Ok(false);
    }
    if a.hit.fixed != b.hit.fixed {
        return Err(failure(
            "two overlapping border points were isolated on different strata",
        ));
    }
    // Each point is the only root of the shared face equations in its own
    // face box, and the B-spline is continuous across knots, so a point
    // inside the other's face box is that other point.
    if inside(&a.hit.bounds, &b.hit.face, dec)? || inside(&b.hit.bounds, &a.hit.face, dec)? {
        return Ok(true);
    }
    let free = &a.hit.free;
    let mut lo = Vec::with_capacity(free.len());
    let mut hi = Vec::with_capacity(free.len());
    for &axis in free {
        lo.push(
            a.hit.bounds[axis][0]
                .clone()
                .min(b.hit.bounds[axis][0].clone()),
        );
        hi.push(
            a.hit.bounds[axis][1]
                .clone()
                .max(b.hit.bounds[axis][1].clone()),
        );
    }
    // Both ends solve the same face equations on their shared stratum; a
    // knot line's two spans agree there exactly.
    match isolate(&a.hit.equations, &lo, &hi, dec)? {
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

fn point(end: &Endpoint) -> PairPoint {
    PairPoint {
        parameters: end.hit.parameters(),
        enclosure: end.hit.bounds.clone(),
        first_span: end.hit.spans[0],
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
                .entry((*axis, value.clone()))
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
                if same_point(ends[a], ends[b], dec)? {
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
                .any(|(axis, value)| ctx.physical(*axis, value))
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
                first_span: arc.first.span,
                first_cell: arc.first.rect.clone(),
                second_span: arc.second_span,
                second_cell: arc.hull.clone(),
                certificate: arc.certificate,
                ends: [point(&arc.ends[0]), point(&arc.ends[1])],
            });
        }
        out.push(PatchPairBranch {
            arcs: list.into(),
            closed,
        });
    }
    Ok(out)
}
