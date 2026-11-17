//! Conservative sweep-and-prune; unknown bounds always reach the narrow phase.

use catallax_numbers::{Decider, Rational, Refusal};

pub(crate) type Box3 = [[Rational; 2]; 3];

pub(crate) fn pairs(
    bounds: &[Option<Box3>],
    dec: &mut Decider,
    mut visit: impl FnMut(usize, usize, &mut Decider) -> Result<(), Refusal>,
) -> Result<(), Refusal> {
    dec.charge_work(bounds.len() as u64)?;
    let mut ordered = Vec::with_capacity(bounds.len());
    let mut unknown = Vec::new();
    for (index, bound) in bounds.iter().enumerate() {
        match bound {
            Some(bound) => ordered.push((index, bound)),
            None => unknown.push(index),
        }
    }
    if ordered.len() > 1 {
        let sorting_work =
            (ordered.len() as u64).saturating_mul(u64::from(ordered.len().ilog2()) + 1);
        dec.charge_work(ordered.len() as u64)?;
        let mut ends = ordered.clone();
        let mut best_axis = 0;
        let mut best_count = u128::MAX;
        // Count interval overlaps without enumerating them. Picking the least
        // crowded axis also handles separated parts aligned along y or z.
        for axis in 0..3 {
            dec.charge_work(sorting_work.saturating_mul(2))?;
            ordered.sort_unstable_by(|a, b| a.1[axis][0].cmp(&b.1[axis][0]).then(a.0.cmp(&b.0)));
            ends.sort_unstable_by(|a, b| a.1[axis][1].cmp(&b.1[axis][1]).then(a.0.cmp(&b.0)));
            dec.check_deadline()?;
            let mut expired = 0;
            let mut count = 0u128;
            for (position, (_, bound)) in ordered.iter().enumerate() {
                dec.charge_work(2)?;
                while expired < ends.len() && ends[expired].1[axis][1] < bound[axis][0] {
                    expired += 1;
                }
                count += (position - expired) as u128;
            }
            if count < best_count {
                best_axis = axis;
                best_count = count;
            }
        }
        if best_axis != 2 {
            dec.charge_work(sorting_work)?;
            ordered.sort_unstable_by(|a, b| {
                a.1[best_axis][0]
                    .cmp(&b.1[best_axis][0])
                    .then(a.0.cmp(&b.0))
            });
            dec.check_deadline()?;
        }
        // Reuse the end-order storage as the active interval set.
        ends.clear();
        for &(index, bound) in &ordered {
            dec.charge_work((ends.len() as u64).saturating_mul(2).saturating_add(1))?;
            ends.retain(|(_, other)| other[best_axis][1] >= bound[best_axis][0]);
            for &(other_index, other) in &ends {
                if (0..3).all(|axis| {
                    axis == best_axis
                        || (bound[axis][0] <= other[axis][1] && other[axis][0] <= bound[axis][1])
                }) {
                    visit(index.min(other_index), index.max(other_index), dec)?;
                }
            }
            ends.push((index, bound));
        }
    }
    for index in unknown {
        dec.charge_work(bounds.len() as u64)?;
        for (other, bound) in bounds.iter().enumerate() {
            if other == index || (other < index && bound.is_none()) {
                continue;
            }
            visit(index.min(other), index.max(other), dec)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn disjoint_boxes_prune_but_touching_and_unknown_boxes_survive() {
        let interval = |lo, hi| [Rational::from_i64(lo), Rational::from_i64(hi)];
        let bounds = [
            Some([interval(0, 1), interval(0, 1), interval(0, 1)]),
            Some([interval(1, 2), interval(0, 1), interval(0, 1)]),
            Some([interval(10, 11), interval(0, 1), interval(0, 1)]),
            None,
            None,
        ];
        let mut actual = BTreeSet::new();
        pairs(&bounds, &mut Decider::default(), |a, b, _| {
            assert!(actual.insert((a, b)), "duplicate candidate");
            Ok(())
        })
        .unwrap();
        assert_eq!(
            actual,
            BTreeSet::from([
                (0, 1),
                (0, 3),
                (0, 4),
                (1, 3),
                (1, 4),
                (2, 3),
                (2, 4),
                (3, 4),
            ])
        );
    }
}
