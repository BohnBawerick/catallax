//! Authored segment allocation, exact equal-angle arcs and rational refinement.

use super::*;
use catallax_recipe::{ElementId, Loop, Segment, Turn};
use catallax_shapes::{BSplineCurve, Point3, Vector3};

#[cfg(test)]
mod tests;

pub(super) struct Chart {
    pub bounds: [Rational; 2],
    pub curve: BSplineCurve,
}

fn point(segment: &Segment) -> Result<[Rational; 2], Refusal> {
    match segment {
        Segment::Line { to, .. } | Segment::Arc { to, .. } | Segment::Bezier { to, .. } => {
            Ok(crate::profile::point(to))
        }
        _ => Err(unproven(
            "a non-endpoint segment cannot occur inside an open chart chain",
        )),
    }
}

fn world(profile: &Profile, point: [Real; 2]) -> Point3 {
    let at = profile.pose.rotation.apply(&Vector3::new(
        point[0].clone(),
        point[1].clone(),
        Real::zero(),
    ));
    Point3::new(
        &at.x + &profile.pose.translation.x,
        &at.y + &profile.pose.translation.y,
        &at.z + &profile.pose.translation.z,
    )
}

fn homogeneous(profile: &Profile, point: [Real; 2], weight: Real) -> [Real; 4] {
    let point = world(profile, point);
    [
        &point.x * &weight,
        &point.y * &weight,
        &point.z * &weight,
        weight,
    ]
}

fn curve(
    profile: &Profile,
    points: Vec<([Real; 2], Real)>,
    bounds: [Rational; 2],
    dec: &mut Decider,
) -> Result<Chart, Refusal> {
    let degree = points.len() - 1;
    let knots = std::iter::repeat_n(bounds[0].clone(), degree + 1)
        .chain(std::iter::repeat_n(bounds[1].clone(), degree + 1))
        .collect();
    let controls = points
        .into_iter()
        .map(|(p, w)| homogeneous(profile, p, w))
        .collect();
    Ok(Chart {
        bounds,
        curve: BSplineCurve::from_homogeneous(degree, knots, controls, dec)?,
    })
}

fn arc_step(
    a: &[Rational; 2],
    b: &[Rational; 2],
    turn: Turn,
    dec: &mut Decider,
) -> Result<(usize, Real, Real), Refusal> {
    let radius = &a[0] * &a[0] + &a[1] * &a[1];
    let dot = &a[0] * &b[0] + &a[1] * &b[1];
    let mut cross = &a[0] * &b[1] - &a[1] * &b[0];
    if turn == Turn::Cw {
        cross = -cross;
    }
    let count = if cross.sign() != Sign::Negative {
        if dot.sign() != Sign::Negative { 1 } else { 2 }
    } else if dot.sign() != Sign::Positive {
        3
    } else {
        4
    };
    let cosine = Real::from_rational(&dot / &radius);
    let half = Real::from_rational(Rational::from_ratio(1, 2).expect("nonzero"));
    let c = match count {
        1 => cosine,
        2 => (&(&Real::one() + &cosine) * &half).sqrt(),
        3 => {
            let equation = RealPoly::from_terms(
                1,
                vec![
                    (vec![3], Real::from_i64(4)),
                    (vec![1], Real::from_i64(-3)),
                    (vec![0], -cosine),
                ],
            );
            let mut selected = None;
            for root in equation.real_roots(dec)? {
                if sign(&root, dec)? != Sign::Negative
                    && sign(&(&root - &half), dec)? != Sign::Positive
                {
                    if selected.is_some() {
                        return Err(unproven("equal-angle cubic branch is not unique"));
                    }
                    selected = Some(root);
                }
            }
            selected.ok_or_else(|| unproven("equal-angle cubic branch was not proved"))?
        }
        4 => (&(&Real::one() - &(&(&Real::one() + &cosine) * &half).sqrt()) * &half).sqrt(),
        _ => unreachable!(),
    };
    let mut s = (Real::one() - c.square()).sqrt();
    if turn == Turn::Cw {
        s = -s;
    }
    Ok((count, c, s))
}

pub(super) fn build(
    profile: &Profile,
    lp: &Loop,
    start: &ElementId,
    dec: &mut Decider,
) -> Result<(Vec<Chart>, Vec<ElementId>, bool), Refusal> {
    let first =
        lp.0.iter()
            .position(|s| s.id() == start)
            .ok_or_else(|| unproven("the authored loft start is absent"))?;
    let ordered: Vec<_> = lp.0[first..].iter().chain(&lp.0[..first]).collect();
    let ids = ordered.iter().map(|s| s.id().clone()).collect();
    let mut result = Vec::new();
    let m = ordered.len();
    for (index, segment) in ordered.iter().enumerate() {
        let low = Rational::from_ratio(index as i64, m as i64).expect("nonempty loop");
        let high = Rational::from_ratio((index + 1) as i64, m as i64).expect("nonempty loop");
        match segment {
            Segment::Line { .. } | Segment::Bezier { .. } => {
                let p0 = point(ordered[(index + m - 1) % m])?;
                let end = point(segment)?;
                let mut controls = vec![(p0.map(Real::from), Real::one())];
                if let Segment::Bezier {
                    controls: interior, ..
                } = segment
                {
                    controls.extend(
                        interior
                            .iter()
                            .map(|p| (crate::profile::point(p).map(Real::from), Real::one())),
                    );
                }
                controls.push((end.map(Real::from), Real::one()));
                result.push(curve(profile, controls, [low, high], dec)?);
            }
            Segment::Arc { centre, .. } | Segment::Circle { centre, .. } => {
                let centre = crate::profile::point(centre);
                let (a, end, count, c, s) = if let Segment::Circle { radius, .. } = segment {
                    let p = [radius.value().clone(), Rational::zero()];
                    (p.clone(), p, 4, Real::zero(), Real::one())
                } else {
                    let start = point(ordered[(index + m - 1) % m])?;
                    let finish = point(segment)?;
                    let a: [Rational; 2] = std::array::from_fn(|axis| &start[axis] - &centre[axis]);
                    let end: [Rational; 2] =
                        std::array::from_fn(|axis| &finish[axis] - &centre[axis]);
                    let Segment::Arc { turn, .. } = segment else {
                        unreachable!()
                    };
                    let (count, c, s) = arc_step(&a, &end, *turn, dec)?;
                    (a, end, count, c, s)
                };
                let weight = ((&Real::one() + &c) / Real::from_i64(2)).sqrt();
                let centre = centre.map(Real::from);
                let mut a = a.map(Real::from);
                for span in 0..count {
                    let b = if span + 1 == count {
                        end.clone().map(Real::from)
                    } else {
                        [&a[0] * &c - &a[1] * &s, &a[0] * &s + &a[1] * &c]
                    };
                    let middle = std::array::from_fn(|axis| {
                        &centre[axis] + &((&a[axis] + &b[axis]) / (&Real::one() + &c))
                    });
                    let p0 = std::array::from_fn(|axis| &centre[axis] + &a[axis]);
                    let p2 = std::array::from_fn(|axis| &centre[axis] + &b[axis]);
                    let bounds = [span, span + 1].map(|part| {
                        Rational::from_ratio((index * count + part) as i64, (m * count) as i64)
                            .expect("nonempty chart")
                    });
                    result.push(curve(
                        profile,
                        vec![
                            (p0, Real::one()),
                            (middle, weight.clone()),
                            (p2, Real::one()),
                        ],
                        bounds,
                        dec,
                    )?);
                    a = b;
                }
            }
            _ => {
                return Err(unproven(
                    "loft charts currently cover lines, circular arcs and polynomial splines",
                ));
            }
        }
    }
    // Profile validation stores the material-left traversal. Determine whether
    // the authored traversal follows it, without using merged segment counts.
    let source = ordered[0];
    let normalized = profile
        .loops
        .iter()
        .flat_map(|b| &b.segments)
        .find(|s| s.ids.contains(source.id()))
        .ok_or_else(|| unproven("loft source segment has no validated profile boundary"))?;
    let material_forward = match source {
        Segment::Circle { .. } => normalized
            .circle
            .as_ref()
            .is_some_and(|(_, _, turn)| *turn == Turn::Ccw),
        Segment::Arc { turn, .. } => normalized
            .circle
            .as_ref()
            .is_some_and(|(_, _, normalized)| turn == normalized),
        _ => {
            let start = point(ordered[m - 1])?;
            let end = point(source)?;
            let raw = [&end[0] - &start[0], &end[1] - &start[1]];
            let normalized = [
                &normalized.end[0] - &normalized.start[0],
                &normalized.end[1] - &normalized.start[1],
            ];
            (&raw[0] * &normalized[0] + &raw[1] * &normalized[1]).sign() == Sign::Positive
        }
    };
    Ok((result, ids, material_forward))
}
