//! Nonperiodic minor-angle lifts for toroidal bands and contractible patches.

use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Loop, Model, Sense};

use crate::{
    chart::{self, Chart, Uv, align},
    edges::Boundaries,
};

fn rings(
    model: &Model,
    face: FaceId,
    chart: &Chart,
    boundaries: &Boundaries,
) -> Option<Vec<Vec<Uv>>> {
    let mut rings = Vec::new();
    for &boundary in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(boundary) else {
            return None;
        };
        let mut ring = Vec::new();
        for &coedge in coedges {
            let coedge = model.coedge(coedge);
            let vertices = &boundaries.mesh.edges[boundaries.edges[&coedge.edge]].vertices;
            for i in 0..vertices.len().saturating_sub(1) {
                let index = if coedge.sense == Sense::Same {
                    i
                } else {
                    vertices.len() - 1 - i
                };
                ring.push(
                    chart.canonical(chart.project(boundaries.mesh.positions[vertices[index]])),
                );
            }
        }
        if ring.len() < 3 {
            return None;
        }
        rings.push(ring);
    }
    (!rings.is_empty()).then_some(rings)
}

fn lift(value: f64, center: f64) -> f64 {
    value + libm::round(center - value)
}

fn cut_is_clear(rings: &[Vec<Uv>], dimension: usize, center: f64) -> bool {
    rings.iter().all(|ring| {
        (0..ring.len()).all(|i| {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            let natural = align(b, a, [true, true]);
            chart::q(lift(b[dimension], center)) - chart::q(lift(a[dimension], center))
                == chart::q(natural[dimension]) - chart::q(a[dimension])
        })
    })
}

fn contractible_center(ring: &[Uv], same: bool, dec: &Decider) -> Result<Option<f64>, Refusal> {
    let mut previous = ring[0];
    let mut bounds = previous.map(|v| [v, v]);
    let mut area = Rational::zero();
    for i in 1..=ring.len() {
        dec.check_deadline()?;
        let next = align(ring[i % ring.len()], previous, [true, true]);
        area = area + chart::q(previous[0]) * chart::q(next[1])
            - chart::q(previous[1]) * chart::q(next[0]);
        for dimension in 0..2 {
            bounds[dimension][0] = bounds[dimension][0].min(next[dimension]);
            bounds[dimension][1] = bounds[dimension][1].max(next[dimension]);
        }
        previous = next;
    }
    // A closed lift inside one period bounds a disk. Its oriented interior,
    // rather than the complementary punctured torus, must be the retained face.
    if previous.map(f64::to_bits) != ring[0].map(f64::to_bits)
        || bounds.iter().any(|[lo, hi]| hi - lo >= 1.0 || hi <= lo)
        || area.sign() != if same { Sign::Positive } else { Sign::Negative }
    {
        return Ok(None);
    }
    Ok(Some((bounds[1][0] + bounds[1][1]) * 0.5))
}

pub(crate) fn minor_angle_center(
    model: &Model,
    face: FaceId,
    original: &Chart,
    boundaries: &Boundaries,
    dec: &Decider,
) -> Result<Option<f64>, Refusal> {
    if !matches!(original.surface, Surface::Torus(_)) {
        return Ok(None);
    }
    let Some(rings) = rings(model, face, original, boundaries) else {
        return Ok(None);
    };
    dec.check_deadline()?;
    let same = model.face(face).sense == Sense::Same;
    if let [ring] = rings.as_slice() {
        let Some(center) = contractible_center(ring, same, dec)? else {
            return Ok(None);
        };
        return Ok(cut_is_clear(&rings, 1, center).then_some(center));
    }
    if rings.len() != 2 {
        return Ok(None);
    }
    let mut means = [0.0; 2];
    let mut increasing = [false; 2];
    for (index, ring) in rings.iter().enumerate() {
        let mut previous = ring[0];
        let mut winding = [0.0; 2];
        let mut direction = None;
        let mut mean = previous[1];
        for i in 1..=ring.len() {
            let next = align(ring[i % ring.len()], previous, original.periodic);
            let delta = next[0] - previous[0];
            if delta == 0.0 || direction.is_some_and(|forward| forward != (delta > 0.0)) {
                return Ok(None);
            }
            direction = Some(delta > 0.0);
            for k in 0..2 {
                winding[k] += next[k] - previous[k];
            }
            if i < ring.len() {
                mean += next[1];
            }
            previous = next;
        }
        if winding[0].abs().to_bits() != 1.0f64.to_bits() || winding[1] != 0.0 {
            return Ok(None);
        }
        increasing[index] = same == direction.unwrap();
        means[index] = mean / ring.len() as f64;
    }
    if increasing[0] == increasing[1] {
        return Ok(None);
    }
    if !increasing[0] {
        means.swap(0, 1);
    }
    means[1] += libm::floor(means[0] - means[1]) + 1.0;
    let center = (means[0] + means[1]) * 0.5;
    if !cut_is_clear(&rings, 1, center) {
        return Ok(None);
    }
    Ok(Some(center))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complements_and_noncontractible_loops_keep_periodic_domains() {
        let ring = [[0.9, 0.9], [0.1, 0.9], [0.1, 0.1], [0.9, 0.1]];
        let dec = Decider::default();
        assert!(contractible_center(&ring, false, &dec).unwrap().is_none());
        let reversed: Vec<_> = ring.into_iter().rev().collect();
        assert!(
            contractible_center(&reversed, true, &dec)
                .unwrap()
                .is_none()
        );
        let noncontractible = [[0.0, 0.25], [0.25, 0.25], [0.5, 0.25], [0.75, 0.25]];
        assert!(
            contractible_center(&noncontractible, true, &dec)
                .unwrap()
                .is_none()
        );
    }
}
