//! Circular charts for planar caps whose source supports share one cylinder axis.
use super::{Chart, Constraint, Radial, q};
use crate::{Face, Plan};
use catallax_numbers::{Decider, Interval, Rational, Refusal};
use catallax_shapes::{Cylinder, Plane, RPoint3, Surface};

fn zero(v: &catallax_shapes::Vector3) -> bool {
    v.rational_components()
        .is_some_and(|v| v.iter().all(|q| q.is_zero()))
}
fn parallel(a: &catallax_shapes::Vector3, b: &catallax_shapes::Vector3) -> bool {
    zero(&a.cross(b))
}
fn coaxial(a: &Cylinder, b: &Cylinder) -> bool {
    parallel(&a.axis, &b.axis) && zero(&a.origin.minus(&b.origin).cross(&a.axis))
}

pub(super) fn charts(
    plan: &Plan,
    faces: &[usize],
    boxes: &[Option<[Interval; 3]>],
    bound: &Rational,
) -> Result<Option<Vec<Chart>>, Refusal> {
    let Surface::Plane(plane) = &plan.faces[faces[0]].surface else {
        return Ok(None);
    };
    let sources: Vec<&Face> = plan
        .faces
        .iter()
        .filter(|f| faces.iter().any(|&i| plan.faces[i].leaf == f.leaf))
        .collect();
    let Some(cylinder) = sources.iter().find_map(|f| match &f.surface {
        Surface::Cylinder(c) if parallel(&c.axis, &plane.normal) => Some(c),
        _ => None,
    }) else {
        return Ok(None);
    };
    if sources.iter().any(|f| match &f.surface {
        Surface::Plane(p) => {
            !parallel(&p.normal, &cylinder.axis)
                && !p
                    .normal
                    .dot(&cylinder.axis)
                    .rational_ref()
                    .is_some_and(|q| q.is_zero())
        }
        Surface::Cylinder(c) => !coaxial(c, cylinder),
        _ => true,
    }) {
        return Ok(None);
    }
    let Some(boxes) = faces
        .iter()
        .map(|&f| boxes[plan.faces[f].leaf].as_ref())
        .collect::<Option<Vec<_>>>()
    else {
        return Ok(None);
    };
    if boxes.iter().any(|b| b.iter().any(|v| !v.is_finite())) {
        return Ok(None);
    }
    Chart::disk(plane, cylinder, &boxes, bound).map(Some)
}

pub(super) fn radial_constraints(
    plan: &Plan,
    plane: &Plane,
    chart: &Chart,
    constraints: &[Vec<Constraint>],
    dec: &mut Decider,
) -> Result<Radial, Refusal> {
    let centre = chart.point(&[q(0), q(0)]);
    let at_one = chart.point(&[q(0), q(1)]);
    let centre_real = centre.to_real();
    let mut candidates = Vec::new();
    for face in &plan.faces {
        let Surface::Cylinder(c) = &face.surface else {
            continue;
        };
        if parallel(&c.axis, &plane.normal) && zero(&centre_real.minus(&c.origin).cross(&c.axis)) {
            candidates.push(face);
        }
    }
    let evaluate = |p: &catallax_numbers::MPoly, point: &RPoint3| {
        p.eval(&[point.x.clone(), point.y.clone(), point.z.clone()])
    };
    constraints
        .iter()
        .map(|leaf| {
            leaf.iter()
                .map(|c| {
                    let Some(polynomial) = &c.polynomial else {
                        return Ok(None);
                    };
                    for candidate in &candidates {
                        if c.support.equals(&candidate.support, dec)? {
                            // Equality to a coaxial cylinder proves A*r^2+B in
                            // the radial parameter, including nonunit frames.
                            let b = evaluate(polynomial, &centre);
                            let a = evaluate(polynomial, &at_one) - b.clone();
                            return Ok(Some((a, b)));
                        }
                    }
                    Ok(None)
                })
                .collect()
        })
        .collect()
}
