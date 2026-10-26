//! Extrema of the represented, checked finite prism, without target-law padding.

use super::{Accuracy, Budget, Extrema};
use crate::Enclosure;
use catallax_numbers::{Ball, Decider, MarkLog, Rational, Refusal};
use catallax_shapes::{BSplineCurve, Curve, Surface};
use catallax_topology::{Body, BodyKind, FaceId, Loop, Lump, Model, Shell, SplinePrism};
use std::collections::BTreeSet;

#[cfg(test)]
mod tests;

pub(super) fn has_spline_boundary(model: &Model, face: FaceId) -> bool {
    matches!(
        model.geom().surface(model.face(face).surface),
        Surface::BSpline(_)
    ) || model
        .face(face)
        .loops
        .iter()
        .any(|lp| match model.loop_(*lp) {
            Loop::Vertex(_) => false,
            Loop::Coedges(cs) => cs.iter().any(|c| {
                matches!(
                    model.geom().curve(model.edge(model.coedge(*c).edge).curve),
                    Curve::BSpline(_)
                )
            }),
        })
}

fn unsupported() -> Refusal {
    Refusal::new(
        "measure.unsupported-surface",
        "spline face trim has no checked complete finite prism domain",
    )
}

/// Recover only the three incident faces; the topology checker, rather than
/// this adjacency walk, establishes their actual complete parameter domains.
/// A face query has no body handle, so check an immutable scratch shell using
/// the original faces. No new geometric facts, seams or trim assumptions enter.
pub(super) fn prism_faces(
    model: &Model,
    face: FaceId,
    accuracy: &Accuracy,
    budget: &mut Budget,
    dec: &mut Decider,
) -> Result<[(FaceId, Extrema); 3], Refusal> {
    let original_bits = budget.bits;
    // Select working precision from the requested width, keeping guard bits.
    // This is only a starting precision; the extremum-width proof below is
    // authoritative. No caller limit or requested accuracy is relaxed.
    let tolerance = if accuracy.absolute.is_zero() {
        &accuracy.relative
    } else {
        &accuracy.absolute
    };
    if !tolerance.is_zero() {
        let bits = tolerance
            .denominator()
            .bit_length()
            .saturating_sub(tolerance.numerator().bit_length())
            .saturating_add(32);
        budget.bits = u32::try_from(bits)
            .unwrap_or(u32::MAX)
            .max(64)
            .min(budget.policy.max_precision_bits);
    }
    let mut faces = BTreeSet::from([face]);
    let mut pending = vec![face];
    while let Some(face) = pending.pop() {
        dec.check_deadline()?;
        for lp in &model.face(face).loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                return Err(unsupported());
            };
            for c in cs {
                for adjacent in model.coedges_of_edge(model.coedge(*c).edge) {
                    dec.check_deadline()?;
                    let Some((_, owner)) = model.owner_of_coedge(adjacent) else {
                        return Err(unsupported());
                    };
                    if faces.insert(owner) {
                        if faces.len() > 3 {
                            return Err(unsupported());
                        }
                        pending.push(owner);
                    }
                }
            }
        }
    }
    if faces.len() != 3 {
        return Err(unsupported());
    }
    let mut scratch = model.clone();
    let shell = scratch.add_shell(Shell {
        faces: faces.into_iter().collect(),
    });
    let lump = scratch.add_lump(Lump {
        shells: vec![shell],
    });
    let body = scratch.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    let checked = SplinePrism::check(&scratch, body, dec)?;
    let [lower, upper] = checked.rings().map(|edge| {
        let Curve::BSpline(curve) = scratch.geom().curve(scratch.edge(edge).curve) else {
            unreachable!("checked shared spline rim")
        };
        curve
    });
    let lower = curve_extrema(lower, lower.domain(), accuracy, budget, dec)?;
    let upper = curve_extrema(upper, upper.domain(), accuracy, budget, dec)?;
    // The checked support is an exact degree-(3,1) sweep. Every side coordinate
    // is affine in v, so its extrema occur on one of the two complete rims.
    // A linear coordinate on either compact planar cap attains its extrema on
    // the boundary (or is constant on the whole cap). Neither cap inherits the
    // other rim's extrema, and no original offset-law error is included here.
    let mut side = lower.clone();
    side.merge(&upper);
    let [a, b, c] = checked.faces();
    budget.bits = original_bits;
    Ok([(a, lower), (b, upper), (c, side)])
}

struct Panel {
    lo: Rational,
    hi: Rational,
    controls: Vec<[Ball; 4]>,
    extrema: Extrema,
}

fn panel(
    controls: Vec<[Ball; 4]>,
    lo: Rational,
    hi: Rational,
    budget: &mut Budget,
    dec: &Decider,
) -> Result<Panel, Refusal> {
    let points = budget.work(dec, |work| {
        controls
            .iter()
            .map(|p| {
                Ok([
                    work.div(&p[0], &p[3])?,
                    work.div(&p[1], &p[3])?,
                    work.div(&p[2], &p[3])?,
                ])
            })
            .collect::<Result<Vec<_>, _>>()
    })?;
    let a = &points[0];
    let b = points.last().unwrap();
    // Positive homogeneous weights prove that the rational Bezier curve stays
    // in the Cartesian control hull of each subdivided net. The endpoint
    // witnesses prove extremum widths; a hull alone cannot end the query.
    let min = std::array::from_fn(|i| Enclosure {
        lo: points.iter().map(|p| p[i].lo()).min().unwrap(),
        hi: a[i].hi().min(b[i].hi()),
    });
    let max = std::array::from_fn(|i| Enclosure {
        lo: a[i].lo().max(b[i].lo()),
        hi: points.iter().map(|p| p[i].hi()).max().unwrap(),
    });
    Ok(Panel {
        lo,
        hi,
        controls,
        extrema: Extrema { min, max },
    })
}

/// Homogeneous de Casteljau subdivision keeps the represented rational curve,
/// including its endpoints. Interval operations are charged to the same budget.
fn split(
    mut row: Vec<[Ball; 4]>,
    parameter: &Rational,
    budget: &mut Budget,
    dec: &Decider,
) -> Result<[Vec<[Ball; 4]>; 2], Refusal> {
    budget.work(dec, |work| {
        let t = work.rational(parameter)?;
        let one = work.rational(&Rational::one())?;
        let complement = work.sub(&one, &t)?;
        let mut left = vec![row[0].clone()];
        let mut right = vec![row.last().unwrap().clone()];
        for remaining in (1..row.len()).rev() {
            for i in 0..remaining {
                let (left, right) = row.split_at_mut(i + 1);
                for (value, next) in left[i].iter_mut().zip(&right[0]) {
                    let a = work.mul(&complement, value)?;
                    let b = work.mul(&t, next)?;
                    *value = work.add(&a, &b)?;
                }
            }
            left.push(row[0].clone());
            right.push(row[remaining - 1].clone());
        }
        right.reverse();
        Ok([left, right])
    })
}

fn curve_extrema(
    curve: &BSplineCurve,
    domain: [&Rational; 2],
    accuracy: &Accuracy,
    budget: &mut Budget,
    dec: &Decider,
) -> Result<Extrema, Refusal> {
    // A decider-aware seam jet supplies an attained finite-domain witness.
    // Subdivision carries the remaining witnesses without reevaluating the
    // algebraic defining expressions and derivative nets on every child panel.
    let jet = curve.jet_box([domain[0], domain[0]], budget.bits, dec)?;
    let point = jet.value.map(|v| Enclosure::from_ball(&v));
    let witness = Extrema {
        min: point.clone(),
        max: point,
    };
    let mut panels = Vec::new();
    for (span, [start, end]) in curve.spans() {
        let lo = start.max(domain[0]);
        let hi = end.min(domain[1]);
        if lo < hi {
            let degree = curve.degree();
            // SplinePrism proves the existing composite Bezier family. Do not
            // treat arbitrary B-spline controls as extracted Bezier controls.
            if !curve.knots()[span - degree + 1..=span]
                .iter()
                .all(|k| k == start)
                || !curve.knots()[span + 1..=span + degree]
                    .iter()
                    .all(|k| k == end)
            {
                return Err(unsupported());
            }
            let mut controls = budget.work(dec, |work| {
                (span - degree..=span)
                    .map(|i| {
                        let p = curve.control(i);
                        let w = work.real(curve.weight(i))?;
                        let mut point = [
                            work.real(&p.x)?,
                            work.real(&p.y)?,
                            work.real(&p.z)?,
                            w.clone(),
                        ];
                        for coordinate in &mut point[..3] {
                            *coordinate = work.mul(coordinate, &w)?;
                        }
                        Ok(point)
                    })
                    .collect::<Result<Vec<_>, _>>()
            })?;
            if hi < end {
                let [left, _] = split(controls, &((hi - start) / (end - start)), budget, dec)?;
                controls = left;
            }
            if lo > start {
                let [_, right] = split(controls, &((lo - start) / (hi - start)), budget, dec)?;
                controls = right;
            }
            panels.push(panel(controls, lo.clone(), hi.clone(), budget, dec)?);
        }
    }
    loop {
        dec.check_deadline()?;
        if panels.is_empty() {
            return Err(unsupported());
        }
        let mut total = witness.clone();
        for p in &panels {
            total.merge(&p.extrema);
        }
        if total.reached(accuracy) {
            return Ok(total);
        }
        let mut refined = Vec::new();
        for p in panels {
            dec.check_deadline()?;
            let relevant = (0..3).any(|i| {
                (!accuracy.reached(&total.min[i]) && p.extrema.min[i].lo <= total.min[i].hi)
                    || (!accuracy.reached(&total.max[i]) && p.extrema.max[i].hi >= total.max[i].lo)
            });
            if relevant {
                let mid = (&p.lo + &p.hi).mul_pow2(-1);
                let [left, right] = split(p.controls, &Rational::one().mul_pow2(-1), budget, dec)?;
                refined.push(panel(left, p.lo, mid.clone(), budget, dec)?);
                refined.push(panel(right, mid, p.hi, budget, dec)?);
            } else {
                refined.push(p);
            }
        }
        panels = refined;
    }
}
