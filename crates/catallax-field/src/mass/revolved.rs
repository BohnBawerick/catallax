//! Exact meridian cells for Boolean combinations of coaxial polygon revolutions.
//!
//! Every profile vertex height and every radial edge crossing splits the axis.
//! Between consecutive heights, radial boundaries have a fixed affine order.
//! The resulting open cells meet no leaf boundary. Exact field membership at
//! one point therefore settles the entire cell, whose moments are polynomial
//! integrals. Caps, rails and the axis have zero volume. No mesh or B-rep is read.

use super::{Accuracy, Options, Properties, exact, failure, q, reached, valid};
use crate::{Instruction, Location, Plan, features::LinearRevolution, geometry::Leaf, profile};
use catallax_features::stand_in::MassTarget;
use catallax_measure::Enclosure;
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_shapes::{RPoint3, RVector3};
use std::collections::BTreeSet;

struct Edge {
    height: [Rational; 2],
    // r(h) = line[0] + line[1] * h.
    line: [Rational; 2],
}
struct Section {
    origin: RPoint3,
    axis: RVector3,
    radial: RVector3,
    edges: Vec<Edge>,
    heights: BTreeSet<Rational>,
}
fn coverage(detail: &str) -> Refusal {
    failure("field.mass-revolution-coverage", detail)
}
fn coordinates(p: &RPoint3) -> [Rational; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
fn components(v: &RVector3) -> [Rational; 3] {
    [v.x.clone(), v.y.clone(), v.z.clone()]
}
fn at(line: &[Rational; 2], height: &Rational) -> Rational {
    &line[0] + &(&line[1] * height)
}
fn definitions(plan: &Plan) -> Result<Vec<LinearRevolution>, Refusal> {
    plan.leaves
        .iter()
        .enumerate()
        .filter(|(i, _)| plan.leaf_registers[*i].is_some())
        .map(|(_, leaf)| {
            let Leaf::Feature(feature) = leaf else {
                return Err(coverage(
                    "every live leaf must be a full polygonal revolution",
                ));
            };
            feature
                .linear_revolution()
                .ok_or_else(|| coverage("curved or partial meridian, or nonrational frame"))
        })
        .collect()
}

fn section(
    definitions: &[LinearRevolution],
    options: &Options,
    dec: &mut Decider,
) -> Result<Section, Refusal> {
    let first = definitions
        .first()
        .ok_or_else(|| coverage("no revolved material leaves"))?;
    let origin = first.origin.clone();
    let norm = profile::sqrt(&first.axis.norm_squared())
        .filter(|n| !n.is_zero())
        .ok_or_else(|| coverage("axis must admit a rational unit direction"))?;
    let axis = first.axis.scaled(&q(1).checked_div(&norm).unwrap());
    let mut reference = None;
    let mut edges = Vec::new();
    let mut heights = BTreeSet::new();
    for definition in definitions {
        dec.charge_work(1)?;
        if !definition.axis.cross(&axis).is_zero()
            || !definition.origin.minus(&origin).cross(&axis).is_zero()
        {
            return Err(coverage("live revolution axes are not coaxial"));
        }
        let mut leaf_radial: Option<RVector3> = None;
        let mut loops = Vec::new();
        for points in &definition.loops {
            let mut meridian = Vec::new();
            for point in points {
                dec.charge_work(1)?;
                let offset = point.minus(&origin);
                let height = offset.dot(&axis);
                let radial = offset.minus(&axis.scaled(&height));
                let radius = profile::sqrt(&radial.norm_squared())
                    .ok_or_else(|| coverage("meridian radius must be rational"))?;
                if !radius.is_zero() {
                    if let Some(direction) = &leaf_radial {
                        if !radial.cross(direction).is_zero() || radial.dot(direction) < q(0) {
                            return Err(coverage("a meridian must stay on one side of its axis"));
                        }
                    } else {
                        leaf_radial = Some(radial.clone());
                    }
                    reference
                        .get_or_insert_with(|| radial.scaled(&q(1).checked_div(&radius).unwrap()));
                }
                heights.insert(height.clone());
                meridian.push([radius, height]);
            }
            loops.push(meridian);
        }
        for points in loops {
            for i in 0..points.len() {
                let a = &points[i];
                let b = &points[(i + 1) % points.len()];
                if a[1] == b[1] {
                    continue;
                }
                if edges.len() == options.max_cells {
                    return Err(failure(
                        "field.mass-budget",
                        format!(
                            "{} meridian edges exhaust cell limit {}",
                            edges.len(),
                            options.max_cells
                        ),
                    ));
                }
                let slope = (&b[0] - &a[0]).checked_div(&(&b[1] - &a[1])).unwrap();
                edges.push(Edge {
                    height: [
                        a[1].clone().min(b[1].clone()),
                        a[1].clone().max(b[1].clone()),
                    ],
                    line: [&a[0] - &(&slope * &a[1]), slope],
                });
            }
        }
    }
    for (i, a) in edges.iter().enumerate() {
        for b in &edges[..i] {
            dec.charge_work(1)?;
            let Some(height) = (&b.line[0] - &a.line[0]).checked_div(&(&a.line[1] - &b.line[1]))
            else {
                continue;
            };
            if height > a.height[0]
                && height < a.height[1]
                && height > b.height[0]
                && height < b.height[1]
            {
                heights.insert(height);
            }
        }
    }
    Ok(Section {
        origin,
        axis,
        radial: reference.ok_or_else(|| coverage("meridians have no positive radius"))?,
        edges,
        heights,
    })
}

fn inside(plan: &Plan, point: &RPoint3, dec: &mut Decider) -> Result<bool, Refusal> {
    let point = coordinates(point);
    let mut states = Vec::with_capacity(plan.instructions.len());
    for instruction in &plan.instructions {
        dec.charge_work(1)?;
        states.push(match *instruction {
            Instruction::Empty => Location::Outside,
            Instruction::Leaf(i) => plan.leaves[i].location(&point, dec)?,
            Instruction::Boolean { kind, a, b } => crate::combine(kind, states[a], states[b]),
        });
    }
    match states[plan.root] {
        Location::Inside => Ok(true),
        Location::Outside => Ok(false),
        _ => Err(failure(
            "field.mass-revolution-classification",
            format!("open meridian cell membership unresolved at {point:?}"),
        )),
    }
}

fn integral(
    line: &[Rational; 2],
    power: usize,
    height_power: usize,
    limits: &[Rational; 2],
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let binomial: &[i64] = match power {
        2 => &[1, 2, 1],
        4 => &[1, 4, 6, 4, 1],
        _ => unreachable!(),
    };
    let mut result = q(0);
    for (k, &coefficient) in binomial.iter().enumerate() {
        dec.charge_work(1)?;
        let n = (k + height_power + 1) as i32;
        let factor = &q(coefficient)
            * &(&line[0].pow((power - k) as i32).unwrap() * &line[1].pow(k as i32).unwrap());
        let span = limits[1].pow(n).unwrap() - limits[0].pow(n).unwrap();
        result = result + (factor * span).checked_div(&q(n as i64)).unwrap();
    }
    Ok(result)
}

/// Exact field moments of Boolean combinations of coaxial full revolutions
/// with polygonal meridians and rational unit axes/radii. Every material cell
/// is classified by the original field tape. Moments are rational multiples
/// of pi, enclosed under the caller's unchanged numeric and cell allowances.
/// Other definitions refuse explicitly; existing `properties` is unchanged.
pub fn properties(
    plan: &Plan,
    accuracy: &Accuracy,
    options: &Options,
    dec: &mut Decider,
) -> Result<Properties, Refusal> {
    if options.max_cells == 0
        || ![&accuracy.volume, &accuracy.centroid, &accuracy.inertia]
            .into_iter()
            .all(valid)
    {
        return Err(failure(
            "field.mass-options",
            "nonnegative accuracies and a positive cell limit are required",
        ));
    }
    if plan.is_proven_empty() {
        return Ok(Properties {
            volume: exact(q(0)),
            centroid: None,
            inertia: None,
            estimate: None,
            cells: 0,
        });
    }
    let definitions = definitions(plan)?;
    let Section {
        origin,
        axis,
        radial,
        edges,
        heights,
    } = section(&definitions, options, dec)?;
    let heights: Vec<_> = heights.into_iter().collect();
    let mut moments = [q(0), q(0), q(0), q(0)]; // V/pi, first/pi, axial second/pi, radial second/pi.
    let mut cells = 0;
    for limits in heights.windows(2) {
        let middle = (&limits[0] + &limits[1]).mul_pow2(-1);
        let mut lines: Vec<_> = edges
            .iter()
            .filter(|e| e.height[0] < middle && middle < e.height[1])
            .map(|e| e.line.clone())
            .collect();
        lines.push([q(0), q(0)]);
        lines.sort_by_key(|line| at(line, &middle));
        lines.dedup_by(|a, b| at(a, &middle) == at(b, &middle));
        for pair in lines.windows(2) {
            if cells == options.max_cells {
                return Err(failure(
                    "field.mass-budget",
                    format!(
                        "{cells} meridian cells exhaust cell limit {}",
                        options.max_cells
                    ),
                ));
            }
            cells += 1;
            let radius = (at(&pair[0], &middle) + at(&pair[1], &middle)).mul_pow2(-1);
            let point = origin
                .plus(&axis.scaled(&middle))
                .plus(&radial.scaled(&radius));
            if !inside(plan, &point, dec)? {
                continue;
            }
            let limits = [limits[0].clone(), limits[1].clone()];
            for (i, total) in moments[..3].iter_mut().enumerate() {
                *total = &*total
                    + &(integral(&pair[1], 2, i, &limits, dec)?
                        - integral(&pair[0], 2, i, &limits, dec)?);
            }
            moments[3] = &moments[3]
                + &(integral(&pair[1], 4, 0, &limits, dec)?
                    - integral(&pair[0], 4, 0, &limits, dec)?)
                .mul_pow2(-1);
        }
    }
    if moments[0].is_zero() {
        return Ok(Properties {
            volume: exact(q(0)),
            centroid: None,
            inertia: None,
            estimate: None,
            cells,
        });
    }
    let height = moments[1].checked_div(&moments[0]).unwrap();
    let centre = coordinates(&origin.plus(&axis.scaled(&height)));
    let axial = &moments[3];
    let perpendicular = &moments[3].mul_pow2(-1) + &(&moments[2] - &(&moments[1] * &height));
    let axis = components(&axis);
    let tensor: [[Rational; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let product = &axis[i] * &axis[j];
            &(&perpendicular * &(&q(i64::from(i == j)) - &product)) + &(axial * &product)
        })
    });
    let maximum = dec.policy().max_precision_bits;
    if maximum == 0 {
        return Err(failure(
            "numbers.precision-limit",
            "no moment enclosure precision available",
        ));
    }
    let mut bits = maximum.min(128);
    loop {
        dec.check_deadline()?;
        let pi = Ball::pi_checked(bits, dec)?;
        let enclosure = |factor: &Rational| {
            let value = &Ball::from_rational(factor, bits) * &pi;
            Enclosure {
                lo: value.lo(),
                hi: value.hi(),
            }
        };
        let volume = enclosure(&moments[0]);
        let inertia = std::array::from_fn(|i| std::array::from_fn(|j| enclosure(&tensor[i][j])));
        if reached(&volume, &accuracy.volume)
            && inertia
                .iter()
                .flatten()
                .all(|v| reached(v, &accuracy.inertia))
        {
            let proposal = (&pi.lo() + &pi.hi()).mul_pow2(-1);
            return Ok(Properties {
                volume,
                centroid: Some(centre.clone().map(exact)),
                inertia: Some(inertia),
                estimate: Some(MassTarget {
                    volume: &moments[0] * &proposal,
                    centroid: centre,
                    inertia: tensor.map(|row| row.map(|v| &v * &proposal)),
                }),
                cells,
            });
        }
        if bits == maximum {
            return Err(failure(
                "field.mass-accuracy",
                format!(
                    "{cells} exact meridian cells at {bits} bits; volume width {} mm3, centroid width 0 mm, maximum inertia width {} mm5; requested moments not reached",
                    volume.width().to_decimal_bounds(10).1,
                    inertia
                        .iter()
                        .flatten()
                        .map(Enclosure::width)
                        .max()
                        .unwrap()
                        .to_decimal_bounds(10)
                        .1
                ),
            ));
        }
        bits = bits.saturating_mul(2).min(maximum);
    }
}
