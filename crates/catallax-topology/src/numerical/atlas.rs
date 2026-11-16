//! Atlas faces whose charts cover only part of their retained source.
//!
//! An atlas face evaluates its whole source surface, which can hold other
//! components. Sampled boundary containment alone does not place a filled loop
//! on the face's charts: a ring of charts round a missing centre has a loop
//! whose samples all lie on charts while its interior does not. Admission
//! therefore requires the charts to fill their parameter box exactly. A loop
//! whose boundary lies in that box then encloses only chart points.

use super::*;
use catallax_shapes::RationalAtlas;

type Domain = [[Rational; 2]; 2];

/// Comparisons charged in advance for sorting `entries` exact values:
/// `entries * ceil(log2(entries))`, at least one per entry.
fn sort_bound(entries: usize) -> u64 {
    let depth = usize::BITS - entries.saturating_sub(1).leading_zeros();
    (entries as u64) * u64::from(depth.max(1))
}

/// Whether the closed chart rectangles cover their whole bounding box. The
/// grid of all chart sides splits the box into cells that each lie inside a
/// chart or meet no chart interior, so one exact test per cell decides it.
/// Endpoints are borrowed, and their collection and sorts are charged first.
fn fills_box(domains: &[&Domain], dec: &mut Decider) -> Result<bool, Refusal> {
    let entries = 2 * domains.len();
    dec.charge_work(2 * entries as u64)?;
    let mut cuts: [Vec<&Rational>; 2] = [Vec::with_capacity(entries), Vec::with_capacity(entries)];
    for domain in domains {
        for axis in 0..2 {
            cuts[axis].extend(domain[axis].iter());
        }
    }
    for axis in &mut cuts {
        dec.charge_work(sort_bound(entries))?;
        axis.sort();
        dec.check_deadline()?;
        axis.dedup();
    }
    for u in cuts[0].windows(2) {
        for v in cuts[1].windows(2) {
            dec.charge_work(domains.len() as u64)?;
            if !domains.iter().any(|d| {
                &d[0][0] <= u[0] && u[1] <= &d[0][1] && &d[1][0] <= v[0] && v[1] <= &d[1][1]
            }) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Refuses an atlas face whose charts leave a gap inside their box.
pub(super) fn check_filled(atlas: &RationalAtlas, dec: &mut Decider) -> Result<(), Refusal> {
    dec.charge_work(atlas.charts().len() as u64)?;
    let domains: Vec<&Domain> = atlas
        .charts()
        .iter()
        .map(|chart| chart.patch().domain())
        .collect();
    if fills_box(&domains, dec)? {
        Ok(())
    } else {
        Err(bad(
            "unsupported",
            "an atlas face whose charts leave a gap in their parameter box needs a filled-region coverage check",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(u: i64, v: i64) -> Domain {
        [[u, u + 1], [v, v + 1]].map(|side| side.map(Rational::from_i64))
    }

    #[test]
    fn a_ring_of_charts_round_a_missing_centre_does_not_fill_its_box() {
        let mut dec = Decider::default();
        let full: Vec<_> = (0..3)
            .flat_map(|u| (0..3).map(move |v| cell(u, v)))
            .collect();
        fn refs(domains: &[Domain]) -> Vec<&Domain> {
            domains.iter().collect()
        }
        assert!(fills_box(&refs(&full), &mut dec).unwrap());
        let ring: Vec<_> = full.iter().filter(|d| *d != &cell(1, 1)).cloned().collect();
        assert!(!fills_box(&refs(&ring), &mut dec).unwrap());
        // Unequal charts that tile a box, as the loft roof's 2 by 2 cover does.
        let half = Rational::from_ratio(1, 2).unwrap();
        let roof: Vec<Domain> = [[[0, 1], [0, 1]], [[0, 1], [1, 2]]]
            .into_iter()
            .flat_map(|[_, v]| {
                let v = v.map(Rational::from_i64);
                [
                    [[Rational::zero(), half.clone()], v.clone()],
                    [[half.clone(), Rational::one()], v],
                ]
            })
            .collect();
        assert!(fills_box(&refs(&roof), &mut dec).unwrap());
        assert!(!fills_box(&refs(&roof[1..]), &mut dec).unwrap());
    }

    #[test]
    fn chart_cuts_are_charged_before_they_are_collected_and_sorted() {
        let side = 30;
        let grid: Vec<_> = (0..side)
            .flat_map(|u| (0..side).map(move |v| cell(u, v)))
            .collect();
        let refs: Vec<_> = grid.iter().collect();
        let entries = 2 * grid.len();
        let upfront = 2 * entries as u64 + 2 * sort_bound(entries);
        assert_eq!(sort_bound(entries), entries as u64 * 11);
        let (result, work) = catallax_numbers::with_operation_budget(upfront - 1, || {
            fills_box(&refs, &mut Decider::default())
        });
        assert!(result.is_err() && work.exhausted, "{work:?}");
        let (result, work) = catallax_numbers::with_operation_budget(20_000_000, || {
            fills_box(&refs, &mut Decider::default())
        });
        assert!(
            result.unwrap() && !work.exhausted && work.used >= upfront,
            "{work:?}"
        );
    }
}
