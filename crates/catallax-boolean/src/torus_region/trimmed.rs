//! A parallel or meridian is a coordinate line on S1 x S1. Insert every arc
//! endpoint into the opposite angular partition, then join neighboring open
//! rectangles wherever no finite boundary arc separates them. Virtual support
//! extensions are traversable, including at the parameter seam.

use super::{Phase, add_cut, circle_cut, point_phase};
use crate::geometry::{natural_ends, refusal, sign};
use catallax_numbers::{Decider, Refusal, Sign};
use catallax_shapes::{PointId, Torus};
use catallax_topology::{Loop, LoopId, Model};
use std::{cmp::Ordering, collections::BTreeSet};

// The general chart sampler also admits at most 4096 cells. This arrangement
// has its own finite cap rather than allocating a quadratic grid unchecked.
const MAX_CELLS: usize = 4096;

struct Arc {
    dimension: usize,
    cut: usize,
    extent: Option<[usize; 2]>,
}

fn budget() -> Refusal {
    refusal(
        "boolean.sampling-budget",
        "circular torus arrangement exceeded its deterministic cell allowance",
    )
}

fn ordered(cuts: &[Phase], dec: &mut Decider) -> Result<Vec<usize>, Refusal> {
    let mut order = Vec::with_capacity(cuts.len());
    for i in 0..cuts.len() {
        let mut at = order.len();
        while at > 0 && cuts[i].compare(&cuts[order[at - 1]], dec)? == Ordering::Less {
            at -= 1;
        }
        order.insert(at, i);
    }
    Ok(order)
}

fn root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

fn join(parent: &mut [usize], a: usize, b: usize) {
    let a = root(parent, a);
    let b = root(parent, b);
    parent[a.max(b)] = a.min(b);
}

fn neighbors(
    cuts: &[Phase],
    order: &[usize],
    query: &Phase,
    dec: &mut Decider,
) -> Result<[usize; 2], Refusal> {
    if order.is_empty() {
        return Ok([0, 0]);
    }
    for (i, &cut) in order.iter().enumerate() {
        let before = (i + order.len() - 1) % order.len();
        match query.compare(&cuts[cut], dec)? {
            Ordering::Less => return Ok([before, before]),
            Ordering::Equal => return Ok([before, i]),
            Ordering::Greater => {}
        }
    }
    Ok([order.len() - 1; 2])
}

pub(super) fn signature(
    model: &Model,
    torus: &Torus,
    groups: &[Vec<LoopId>],
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Vec<i64>>, Refusal> {
    let mut cuts = [Vec::new(), Vec::new()];
    let mut arcs = Vec::new();
    let mut seen = BTreeSet::new();
    for &loop_id in groups.iter().flatten() {
        let Loop::Coedges(coedges) = model.loop_(loop_id) else {
            return Ok(None);
        };
        for &coedge in coedges {
            dec.check_deadline()?;
            let edge = model.coedge(coedge).edge;
            if !seen.insert(edge) {
                continue;
            }
            let Some(cut) = circle_cut(model, torus, model.edge(edge).curve, dec)? else {
                return Ok(None);
            };
            let dimension = cut.dimension;
            let fixed = add_cut(&mut cuts[dimension], cut.phase, dec)?;
            let extent = match natural_ends(model, edge) {
                (None, None) => None,
                (Some(from), Some(to)) => {
                    let mut extent = [0; 2];
                    for (slot, endpoint) in extent.iter_mut().zip([from, to]) {
                        *slot = add_cut(
                            &mut cuts[1 - dimension],
                            point_phase(
                                torus,
                                &model.geom().position(endpoint),
                                1 - dimension,
                                dec,
                            )?,
                            dec,
                        )?;
                    }
                    if extent[0] == extent[1] {
                        // A circular edge with coincident ends covers the full ring.
                        None
                    } else {
                        let direction = if dimension == 0 {
                            cut.normal.dot(&cut.offset.cross(&torus.axis))
                        } else {
                            cut.normal.dot(&torus.axis)
                        };
                        match sign(&direction, dec, "boolean.classify-ambiguous")? {
                            Sign::Positive => {}
                            Sign::Negative => extent.swap(0, 1),
                            Sign::Zero => {
                                return Err(refusal(
                                    "boolean.classify-ambiguous",
                                    "circular torus boundary has no angular orientation",
                                ));
                            }
                        }
                        Some(extent)
                    }
                }
                _ => return Ok(None),
            };
            arcs.push(Arc {
                dimension,
                cut: fixed,
                extent,
            });
            if cuts[0]
                .len()
                .max(1)
                .checked_mul(cuts[1].len().max(1))
                .is_none_or(|n| n > MAX_CELLS)
            {
                return Err(budget());
            }
        }
    }
    let order = [ordered(&cuts[0], dec)?, ordered(&cuts[1], dec)?];
    let ranks = order.each_ref().map(|order| {
        let mut rank = vec![0; order.len()];
        for (i, &cut) in order.iter().enumerate() {
            rank[cut] = i;
        }
        rank
    });
    let [nu, nv] = order.each_ref().map(|order| order.len().max(1));
    let count = nu
        .checked_mul(nv)
        .filter(|n| *n <= MAX_CELLS)
        .ok_or_else(budget)?;
    let mut blocked = [vec![false; count], vec![false; count]];
    for arc in arcs {
        dec.check_deadline()?;
        let fixed = ranks[arc.dimension][arc.cut];
        let n = if arc.dimension == 0 { nv } else { nu };
        let (start, length) = match arc.extent {
            None => (0, n),
            Some([from, to]) => {
                let from = ranks[1 - arc.dimension][from];
                let to = ranks[1 - arc.dimension][to];
                (from, (to + n - from) % n)
            }
        };
        for offset in 0..length {
            let moving = (start + offset) % n;
            let index = if arc.dimension == 0 {
                fixed * nv + moving
            } else {
                moving * nv + fixed
            };
            blocked[arc.dimension][index] = true;
        }
    }
    let mut parent: Vec<_> = (0..count).collect();
    for u in 0..nu {
        dec.check_deadline()?;
        for v in 0..nv {
            let index = u * nv + v;
            if !blocked[0][index] {
                join(&mut parent, index, ((u + nu - 1) % nu) * nv + v);
            }
            if !blocked[1][index] {
                join(&mut parent, index, u * nv + (v + nv - 1) % nv);
            }
        }
    }
    let position = model.geom().position(point);
    let u = neighbors(
        &cuts[0],
        &order[0],
        &point_phase(torus, &position, 0, dec)?,
        dec,
    )?;
    let v = neighbors(
        &cuts[1],
        &order[1],
        &point_phase(torus, &position, 1, dec)?,
        dec,
    )?;
    let component = root(&mut parent, u[0] * nv + v[0]);
    for u in u {
        for v in v {
            if root(&mut parent, u * nv + v) != component {
                return Err(refusal(
                    "boolean.classify-ambiguous",
                    "torus cell probe lies on a finite boundary",
                ));
            }
        }
    }
    Ok(Some(vec![component as i64]))
}
