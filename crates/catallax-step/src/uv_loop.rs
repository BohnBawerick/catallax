//! Orientation and nesting of written rational Bezier chart boundaries.

use crate::{fit, number, uv};
use catallax_numbers::{Decider, Rational, Refusal, Sign};

type Point = [Rational; 2];
type Controls = [Point; fit::DEGREE + 1];

pub(crate) struct Curve<'a> {
    pub fit: &'a uv::Fit,
    pub forward: bool,
}

pub(crate) struct Group {
    pub outer: usize,
    pub holes: Vec<usize>,
}

fn controls(points: &[[number::Number; 2]]) -> Controls {
    std::array::from_fn(|i| points[i].each_ref().map(|n| n.value.clone()))
}

/// Twice the signed area; no quadrature or floating-point orientation test.
fn area(wire: &[Curve<'_>]) -> Rational {
    let mut area = Rational::zero();
    for curve in wire {
        let mut contribution = Rational::zero();
        for points in curve
            .fit
            .points
            .windows(fit::DEGREE + 1)
            .step_by(fit::DEGREE)
        {
            let powers = fit::bezier_power(controls(points));
            for i in 0..=fit::DEGREE {
                for j in i + 1..=fit::DEGREE {
                    let factor = Rational::from_ratio((j - i) as i64, (i + j) as i64).unwrap();
                    contribution = contribution
                        + (&powers[i][0] * &powers[j][1] - &powers[i][1] * &powers[j][0]) * factor;
                }
            }
        }
        area = if curve.forward {
            area + contribution
        } else {
            area - contribution
        };
    }
    area
}

fn split(mut points: Controls) -> (Controls, Controls) {
    let mut left: Controls = std::array::from_fn(|_| std::array::from_fn(|_| Rational::zero()));
    let mut right = left.clone();
    left[0] = points[0].clone();
    right[fit::DEGREE] = points[fit::DEGREE].clone();
    for level in 1..=fit::DEGREE {
        for i in 0..=fit::DEGREE - level {
            points[i] =
                std::array::from_fn(|axis| (&points[i][axis] + &points[i + 1][axis]).mul_pow2(-1));
        }
        left[level] = points[0].clone();
        right[fit::DEGREE - level] = points[fit::DEGREE - level].clone();
    }
    (left, right)
}

/// Returns `None` when a span endpoint is the query point itself.
fn winding(wire: &[Curve<'_>], query: &Point, dec: &Decider) -> Result<Option<i64>, Refusal> {
    let mut winding = 0i64;
    for curve in wire {
        let mut contribution = 0i64;
        for points in curve
            .fit
            .points
            .windows(fit::DEGREE + 1)
            .step_by(fit::DEGREE)
        {
            let points =
                controls(points).map(|point| std::array::from_fn(|i| &point[i] - &query[i]));
            let mut pending = vec![(points, 0u32)];
            while let Some((points, depth)) = pending.pop() {
                dec.check_deadline()?;
                if points.iter().all(|p| p[0].sign() == Sign::Negative)
                    || points.iter().all(|p| p[1].sign() == Sign::Positive)
                    || points.iter().all(|p| p[1].sign() == Sign::Negative)
                {
                    continue;
                }
                if points.iter().all(|p| p[0].sign() == Sign::Positive) {
                    // All crossings are on the positive ray. Their oriented
                    // sum depends only on endpoint halfspaces, even at tangency.
                    contribution += i64::from(points[fit::DEGREE][1].sign() == Sign::Positive)
                        - i64::from(points[0][1].sign() == Sign::Positive);
                    continue;
                }
                if points[0].iter().all(Rational::is_zero)
                    || points[fit::DEGREE].iter().all(Rational::is_zero)
                {
                    return Ok(None);
                }
                // Each midpoint introduces one further binary parameter bit.
                if depth >= dec.policy().max_precision_bits {
                    return Err(number::budget(
                        "UV boundary containment exhausted its precision",
                    ));
                }
                let (left, right) = split(points);
                pending.push((right, depth + 1));
                pending.push((left, depth + 1));
            }
        }
        winding += if curve.forward {
            contribution
        } else {
            -contribution
        };
    }
    Ok(Some(winding))
}

pub(crate) fn groups(wires: &[Vec<Curve<'_>>], dec: &Decider) -> Result<Vec<Group>, Refusal> {
    let areas: Vec<_> = wires.iter().map(|wire| area(wire)).collect();
    if areas.iter().any(Rational::is_zero) {
        return Err(number::budget("UV boundary has zero signed area"));
    }
    let mut groups: Vec<_> = areas
        .iter()
        .enumerate()
        .filter(|(_, a)| a.sign() == Sign::Positive)
        .map(|(outer, _)| Group {
            outer,
            holes: Vec::new(),
        })
        .collect();
    if groups.is_empty() {
        return Err(number::budget(
            "periodic chart has no positive outer boundary",
        ));
    }
    for (hole, area) in areas.iter().enumerate() {
        if area.sign() != Sign::Negative {
            continue;
        }
        // A hole may touch its outer boundary at a vertex, as a one-edge ring
        // does at its only vertex. Query the hole's span ends in walk order
        // until one lies on no outer boundary.
        let mut queries = wires[hole].iter().flat_map(|curve| {
            let points = &curve.fit.points;
            let count = points.len() / fit::DEGREE;
            (0..count).map(move |span| {
                let index = span * fit::DEGREE;
                points[if curve.forward {
                    index
                } else {
                    points.len() - 1 - index
                }]
                .each_ref()
                .map(|n| n.value.clone())
            })
        });
        let owner = loop {
            let point = queries.next().ok_or_else(|| {
                number::budget("UV boundary containment query touches a boundary")
            })?;
            let mut owner = None;
            let mut touches = false;
            for (i, group) in groups.iter().enumerate() {
                match winding(&wires[group.outer], &point, dec)? {
                    None => {
                        touches = true;
                        break;
                    }
                    Some(0) => {}
                    Some(1) => {
                        if owner.is_none_or(|previous: usize| {
                            areas[group.outer] < areas[groups[previous].outer]
                        }) {
                            owner = Some(i);
                        }
                    }
                    Some(_) => {
                        return Err(number::budget(
                            "UV outer boundary has no simple winding certificate",
                        ));
                    }
                }
            }
            if !touches {
                break owner;
            }
        };
        groups[owner.ok_or_else(|| number::budget("UV hole has no containing outer boundary"))?]
            .holes
            .push(hole);
    }
    Ok(groups)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_edge_ring_hole_touching_its_outer_boundary_is_nested() {
        // The outer boundary passes through the ring's only vertex at the
        // origin, so the ring's first point is no containment query.
        let outer = [
            uv::Fit::polyline(&[[-2, 0], [0, 0]], 1),
            uv::Fit::polyline(&[[0, 0], [2, 0]], 1),
            uv::Fit::polyline(&[[2, 0], [2, 3], [-2, 3], [-2, 0]], 1),
        ];
        let ring = uv::Fit::polyline(&[[0, 0], [8, 1], [8, 16], [-8, 16], [-8, 1], [0, 0]], 8);
        let wires = vec![
            outer
                .iter()
                .map(|fit| Curve { fit, forward: true })
                .collect(),
            vec![Curve {
                fit: &ring,
                forward: false,
            }],
        ];
        let groups = groups(&wires, &Decider::default()).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!((groups[0].outer, groups[0].holes.as_slice()), (0, &[1][..]));
    }
}
