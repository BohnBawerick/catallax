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

/// Whether the closed chart rectangles cover their whole bounding box. The
/// grid of all chart sides splits the box into cells that each lie inside a
/// chart or meet no chart interior, so one exact test per cell decides it.
fn fills_box(domains: &[Domain], dec: &mut Decider) -> Result<bool, Refusal> {
    let mut cuts: [Vec<Rational>; 2] = [Vec::new(), Vec::new()];
    for domain in domains {
        dec.charge_work(1)?;
        for axis in 0..2 {
            cuts[axis].extend(domain[axis].iter().cloned());
        }
    }
    for axis in &mut cuts {
        axis.sort();
        axis.dedup();
    }
    for u in cuts[0].windows(2) {
        for v in cuts[1].windows(2) {
            dec.charge_work(domains.len() as u64)?;
            if !domains
                .iter()
                .any(|d| d[0][0] <= u[0] && u[1] <= d[0][1] && d[1][0] <= v[0] && v[1] <= d[1][1])
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Refuses an atlas face whose charts leave a gap inside their box.
pub(super) fn check_filled(atlas: &RationalAtlas, dec: &mut Decider) -> Result<(), Refusal> {
    let domains: Vec<Domain> = atlas
        .charts()
        .iter()
        .map(|chart| chart.patch().domain().clone())
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
        assert!(fills_box(&full, &mut dec).unwrap());
        let ring: Vec<_> = full.iter().filter(|d| *d != &cell(1, 1)).cloned().collect();
        assert!(!fills_box(&ring, &mut dec).unwrap());
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
        assert!(fills_box(&roof, &mut dec).unwrap());
        assert!(!fills_box(&roof[1..], &mut dec).unwrap());
    }
}
