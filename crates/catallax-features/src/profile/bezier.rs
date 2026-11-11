//! Exact polynomial boundaries and odd nesting for line and Bezier loops.

use super::*;
use catallax_numbers::{RealPoly, ZeroProof};
use catallax_shapes::{BezierExtrusion, Point3};

mod separation;

fn exact_sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::ProvenCoincident(proof) if proof.proof == ZeroProof::SeparationBound => Err(bad(
            "sketch.bezier-proof",
            "a separation bound is not exact incidence",
        )),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(error) => Err(error),
    }
}

fn carrier(controls: &[P], dec: &mut Decider) -> Result<BezierExtrusion, Refusal> {
    BezierExtrusion::new(
        controls
            .iter()
            .map(|p| RPoint3::new(p[0].clone(), p[1].clone(), Rational::zero()))
            .collect(),
        RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        dec,
    )
}

pub(super) fn contains(controls: &[P], p: &X, dec: &mut Decider) -> Result<bool, Refusal> {
    let curve = carrier(controls, dec)?.lower_edge();
    Ok(curve
        .unique_parameter_of(&Point3::new(p[0].clone(), p[1].clone(), Real::zero()), dec)?
        .is_some())
}

pub(super) fn midpoint(controls: &[P]) -> P {
    let mut row = controls.to_vec();
    for count in (1..row.len()).rev() {
        for i in 0..count {
            row[i] = std::array::from_fn(|j| (&row[i][j] + &row[i + 1][j]) / Rational::from_i64(2));
        }
    }
    row.remove(0)
}

fn controls(segment: &Arc) -> Vec<P> {
    segment
        .bezier
        .clone()
        .unwrap_or_else(|| vec![segment.start.clone(), segment.end.clone()])
}

fn power(controls: &[P]) -> Vec<P> {
    let degree = controls.len() - 1;
    let mut row = controls.to_vec();
    let mut result = vec![row[0].clone()];
    for order in 1..=degree {
        for i in 0..degree + 1 - order {
            row[i] = sub(&row[i + 1], &row[i]);
        }
        let factor = if order == 1 || degree == 3 && order == 2 {
            degree as i64
        } else {
            1
        };
        result.push(row[0].clone().map(|x| x * Rational::from_i64(factor)));
    }
    result
}

fn signed_area(segments: &[Arc]) -> Rational {
    let mut area = Rational::zero();
    for segment in segments {
        let coefficients = power(&controls(segment));
        for (i, a) in coefficients.iter().enumerate() {
            for (j, b) in coefficients.iter().enumerate().skip(1) {
                area = area
                    + cross(a, b) * Rational::from_ratio(j as i64, (2 * (i + j)) as i64).unwrap();
            }
        }
    }
    area
}

fn enclosed(segments: &[Arc], point: &P, dec: &mut Decider) -> Result<bool, Refusal> {
    let mut crossings = 0;
    for segment in segments {
        dec.check_deadline()?;
        let coefficients = power(&controls(segment));
        if coefficients.iter().enumerate().all(|(i, p)| {
            p[1] == if i == 0 {
                point[1].clone()
            } else {
                Rational::zero()
            }
        }) {
            continue;
        }
        let polynomial = |axis: usize| {
            RealPoly::from_terms(
                1,
                coefficients
                    .iter()
                    .enumerate()
                    .map(|(i, p)| (vec![i as u32], Real::from(p[axis].clone())))
                    .collect(),
            )
        };
        let x = polynomial(0);
        let y = polynomial(1) - RealPoly::constant(1, point[1].clone().into());
        for root in y.real_roots(dec)? {
            let from = exact_sign(&root, dec)?;
            let to = exact_sign(&(&root - &Real::one()), dec)?;
            if from == Sign::Negative || to == Sign::Positive {
                continue;
            }
            let mut derivative = y.clone();
            let mut first = None;
            for order in 1..coefficients.len() {
                derivative = derivative.partial(0);
                let sign = exact_sign(&derivative.eval_real(std::slice::from_ref(&root)), dec)?;
                if sign != Sign::Zero {
                    first = Some((order, sign));
                    break;
                }
            }
            let Some((order, sign)) = first else {
                return Err(bad(
                    "sketch.bezier-proof",
                    "the ray root multiplicity is unproved",
                ));
            };
            let counts = if from == Sign::Zero {
                sign == Sign::Positive
            } else if to == Sign::Zero {
                sign == if order % 2 == 0 {
                    Sign::Positive
                } else {
                    Sign::Negative
                }
            } else {
                order % 2 == 1
            };
            if counts
                && exact_sign(&(x.eval_real(&[root]) - Real::from(point[0].clone())), dec)?
                    == Sign::Positive
            {
                crossings += 1;
            }
        }
    }
    Ok(crossings % 2 == 1)
}

pub(super) fn profile(
    plane: &catallax_recipe::Pose,
    loops: &[Loop],
    dec: &mut Decider,
) -> Result<Profile, Refusal> {
    let pose = placement(plane, dec)?;
    let mut boundaries = Vec::new();
    for lp in loops.iter().filter(|lp| !lp.0.is_empty()) {
        let mut segments = Vec::new();
        let end_of = |s: &Segment| match s {
            Segment::Line { to, .. } | Segment::Bezier { to, .. } => Ok(point(to)),
            _ => Err(bad(
                "sketch.bezier-curve-mix",
                "Bezier profiles currently accept lines and polynomial spans",
            )),
        };
        for (i, s) in lp.0.iter().enumerate() {
            dec.check_deadline()?;
            let start = end_of(&lp.0[(i + lp.0.len() - 1) % lp.0.len()])?;
            let end = end_of(s)?;
            if start == end {
                return Err(bad("sketch.zero-segment", "segment endpoints must differ"));
            }
            let bezier = if let Segment::Bezier { controls, .. } = s {
                if !(1..=2).contains(&controls.len()) {
                    return Err(bad(
                        "sketch.bezier-controls",
                        "a quadratic or cubic needs one or two interior controls",
                    ));
                }
                let controls: Vec<_> = std::iter::once(start.clone())
                    .chain(controls.iter().map(point))
                    .chain(std::iter::once(end.clone()))
                    .collect();
                let chord = sub(&end, &start);
                let straight = controls
                    .iter()
                    .all(|p| cross(&sub(p, &start), &chord).is_zero())
                    && controls
                        .windows(2)
                        .all(|p| dot(&sub(&p[1], &p[0]), &chord).sign() != Sign::Negative);
                if straight {
                    None
                } else {
                    carrier(&controls, dec)?;
                    Some(controls)
                }
            } else {
                None
            };
            segments.push(Arc {
                ids: vec![s.id().clone()],
                start,
                end,
                circle: None,
                bezier,
                full: false,
            });
        }
        boundaries.push(Boundary {
            segments,
            depth: 0,
            parent: None,
        });
    }
    for i in 0..boundaries.len() {
        for j in i..boundaries.len() {
            let a = &boundaries[i].segments;
            let b = &boundaries[j].segments;
            for (k, first) in a.iter().enumerate() {
                for (l, second) in b.iter().enumerate() {
                    if i == j && l <= k {
                        continue;
                    }
                    let adjacent = i == j && ((k + 1) % a.len() == l || (l + 1) % a.len() == k);
                    let allowed: Vec<_> = [&first.start, &first.end]
                        .into_iter()
                        .filter(|p| adjacent && (*p == &second.start || *p == &second.end))
                        .cloned()
                        .collect();
                    separation::prove(&controls(first), &controls(second), &allowed, dec).map_err(
                        |e| {
                            if e.code == "sketch.bezier-separation" {
                                bad(
                                    if i == j {
                                        "sketch.self-intersection"
                                    } else {
                                        "sketch.loop-intersection"
                                    },
                                    &e.detail,
                                )
                            } else {
                                e
                            }
                        },
                    )?;
                }
            }
        }
    }
    let mut enclosing = vec![Vec::new(); boundaries.len()];
    for i in 0..boundaries.len() {
        for j in 0..boundaries.len() {
            if i != j
                && enclosed(
                    &boundaries[j].segments,
                    &boundaries[i].segments[0].start,
                    dec,
                )?
            {
                enclosing[i].push(j);
            }
        }
    }
    for i in 0..boundaries.len() {
        boundaries[i].depth = enclosing[i].len();
        boundaries[i].parent = enclosing[i]
            .iter()
            .copied()
            .max_by_key(|j| enclosing[*j].len());
        let area = signed_area(&boundaries[i].segments);
        if area.is_zero() {
            return Err(bad(
                "sketch.zero-area",
                "a profile loop must enclose nonzero area",
            ));
        }
        if (area.sign() == Sign::Positive) != (boundaries[i].depth % 2 == 0) {
            boundaries[i].segments.reverse();
            for s in &mut boundaries[i].segments {
                s.reversed();
            }
        }
    }
    Ok(Profile {
        pose,
        loops: boundaries,
        ellipse_frame: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_area_includes_the_curve_bulge() {
        let p = |x, y| [Rational::from_i64(x), Rational::from_i64(y)];
        for (interior, expected) in [(vec![p(2, 4)], "16/3"), (vec![p(0, 4), p(4, 4)], "48/5")] {
            let controls = std::iter::once(p(0, 0))
                .chain(interior)
                .chain(std::iter::once(p(4, 0)))
                .collect();
            let mut curve = Arc {
                ids: vec![],
                start: p(0, 0),
                end: p(4, 0),
                circle: None,
                bezier: Some(controls),
                full: false,
            };
            let mut base = Arc {
                ids: vec![],
                start: p(4, 0),
                end: p(0, 0),
                circle: None,
                bezier: None,
                full: false,
            };
            // This traversal is clockwise.
            assert_eq!(
                signed_area(&[curve.clone(), base.clone()]),
                -Rational::parse(expected).unwrap()
            );
            curve.reversed();
            base.reversed();
            assert_eq!(
                signed_area(&[base, curve]),
                Rational::parse(expected).unwrap()
            );
        }
    }
}
