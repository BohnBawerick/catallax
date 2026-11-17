//! Vertex-local boundary assembly. Restriction uses homogeneous de Boor
//! coefficients; native/world agreement remains sampled numerical evidence.

use catallax_numbers::{Ball, Decider, NumericalCheck, NumericalSample, Rational, Refusal};
use catallax_shapes::{BSplineCurve, Point3};
use catallax_topology::{Model, Sense, VertexId};

use super::super::super::rolling::RolledChain;
use super::{
    BITS, BoundarySide, JunctionBoundary, JunctionFace, JunctionPatch, StripEnd, StripTermination,
    VertexContact, VertexJunction, bad, strip, upper_distance,
};

fn multiplicity(curve: &BSplineCurve, at: &Rational) -> usize {
    let knots = curve.knots();
    knots.partition_point(|knot| knot <= at) - knots.partition_point(|knot| knot < at)
}

fn restricted_curve(
    source: &BSplineCurve,
    range: [&Rational; 2],
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let domain = source.domain();
    if range[0] >= range[1] || range[0] < domain[0] || range[1] > domain[1] {
        return Err(bad(
            "junction boundary restriction leaves its finite curve domain",
        ));
    }
    let degree = source.degree();
    let spans = source
        .spans()
        .filter(|(_, span)| span[0] < range[1] && span[1] > range[0])
        .count();
    let mut controls = Vec::with_capacity((degree + 1) * spans);
    let mut weights = Vec::with_capacity((degree + 1) * spans);
    let mut knots = Vec::with_capacity((degree + 1) * (spans + 1));
    knots.extend(std::iter::repeat_n(Rational::zero(), degree + 1));
    let width = range[1] - range[0];
    for (span, interval) in source.spans() {
        let a = interval[0].max(range[0]);
        let b = interval[1].min(range[1]);
        if a >= b {
            continue;
        }
        let rows = source.span_bezier_controls_on(span, [a, b], dec)?;
        let skip = usize::from(!controls.is_empty() && multiplicity(source, a) <= degree);
        for [x, y, z, weight] in rows.into_iter().skip(skip) {
            controls.push(Point3::new(&x / &weight, &y / &weight, &z / &weight));
            weights.push(weight);
        }
        let end = &(b - range[0]) / &width;
        let count = if b == range[1] || multiplicity(source, b) > degree {
            degree + 1
        } else {
            degree
        };
        knots.extend(std::iter::repeat_n(end, count));
    }
    BSplineCurve::new(degree, knots, controls, weights, dec)
}

pub(super) fn restrict_boundary(
    boundary: &JunctionBoundary,
    range: [&Rational; 2],
    contacts: [usize; 2],
    dec: &mut Decider,
) -> Result<JunctionBoundary, Refusal> {
    Ok(JunctionBoundary {
        curve: restricted_curve(&boundary.curve, range, dec)?,
        contacts,
        sides: [
            BoundarySide {
                face: boundary.sides[0].face,
                parameters: restricted_curve(&boundary.sides[0].parameters, range, dec)?,
            },
            BoundarySide {
                face: boundary.sides[1].face,
                parameters: restricted_curve(&boundary.sides[1].parameters, range, dec)?,
            },
        ],
    })
}

fn same_end(a: StripEnd, b: StripEnd) -> bool {
    a.chain == b.chain && a.strip == b.strip && a.source_end == b.source_end
}

fn path(
    boundaries: &[JunctionBoundary],
    indices: &[usize],
    contacts: [usize; 2],
    end: StripEnd,
) -> Result<Vec<(usize, Sense)>, Refusal> {
    let mut result = Vec::with_capacity(indices.len());
    let mut at = contacts[0];
    for _ in 0..indices.len() {
        let mut next = None;
        for &index in indices {
            if result.iter().any(|&(used, _)| used == index) {
                continue;
            }
            let boundary = boundaries
                .get(index)
                .ok_or_else(|| bad("termination names a missing boundary"))?;
            if !boundary
                .sides
                .iter()
                .any(|side| matches!(side.face, JunctionFace::Strip(strip) if same_end(strip, end)))
            {
                return Err(bad("termination boundary belongs to a different strip end"));
            }
            let candidate = if boundary.contacts[0] == at {
                Some((Sense::Same, boundary.contacts[1]))
            } else if boundary.contacts[1] == at {
                Some((Sense::Reversed, boundary.contacts[0]))
            } else {
                None
            };
            if let Some((sense, to)) = candidate {
                if next.replace((index, sense, to)).is_some() {
                    return Err(bad("strip-end boundaries have an ambiguous continuation"));
                }
            }
        }
        let (index, sense, to) =
            next.ok_or_else(|| bad("strip-end boundaries do not form one connected path"))?;
        result.push((index, sense));
        at = to;
    }
    if at != contacts[1] {
        return Err(bad(
            "strip-end boundary path does not reach its second rail",
        ));
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn bind_termination(
    model: &Model,
    chains: &[RolledChain],
    end: StripEnd,
    tau: [Rational; 2],
    contacts: [usize; 2],
    boundary_indices: &[usize],
    junction: &mut VertexJunction,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let source = strip(model, chains, junction.source, end)?;
    let boundaries = path(&junction.boundaries, boundary_indices, contacts, end)?;
    let mut support_parameters = [0; 2];
    let mut samples = Vec::with_capacity(2);
    for side in 0..2 {
        let contact = junction
            .contacts
            .get_mut(contacts[side])
            .ok_or_else(|| bad("termination names a missing contact"))?;
        let t = &tau[side];
        let rail = source
            .fillet
            .rail(side)
            .ok_or_else(|| bad("termination has no world contact rail"))?;
        let mut error = upper_distance(&rail.evaluate_box([t, t], BITS, dec)?, &contact.position);
        let native = source
            .fillet
            .parameter_rail(side)
            .ok_or_else(|| bad("termination has no original-support pcurve"))?
            .evaluate_box([t, t], BITS, dec)?;
        let uv = [native[0].mid(), native[1].mid()];
        let enclosed = uv.each_ref().map(|value| Ball::from_rational(value, BITS));
        let face = source.source.faces[side];
        let world = model.geom().surface(model.face(face).surface).evaluate(
            &enclosed[0],
            &enclosed[1],
            dec,
        )?;
        error = error.max(upper_distance(&world, &contact.position));
        if error > *tolerance {
            return Err(bad(
                "termination contact misses its actual rail or original support",
            ));
        }
        support_parameters[side] = contact.parameters.len();
        contact.parameters.push((face, uv));
        samples.push(NumericalSample {
            parameters: vec![t.clone(), Rational::from_i64(side as i64)],
            distance_upper_mm: error,
        });
    }
    junction.terminations.push(StripTermination {
        strip: end,
        tau,
        contacts,
        support_parameters,
        boundaries,
    });
    junction.checks.push(NumericalCheck::Samples {
        name: "fillet.vertex.termination-rails-and-supports".into(),
        scope: format!(
            "{:?}; two actual rail endpoints and original support compositions",
            junction.source
        ),
        samples,
    });
    Ok(())
}

pub(super) fn check_boundaries(
    model: &Model,
    chains: &[RolledChain],
    junction: &mut VertexJunction,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for boundary in &junction.boundaries {
        if boundary.curve.domain() != [&Rational::zero(), &Rational::one()]
            || boundary
                .sides
                .iter()
                .any(|side| side.parameters.domain() != boundary.curve.domain())
        {
            return Err(bad(
                "junction world/native curves need the same normalized parameter",
            ));
        }
        let mut samples = Vec::new();
        for (_, span) in boundary.curve.spans() {
            for index in 0..=8 {
                dec.charge_work(1)?;
                let share = Rational::from_ratio(index, 8).unwrap();
                let t = span[0] + &((span[1] - span[0]) * share);
                let world = boundary.curve.evaluate_box([&t, &t], BITS, dec)?;
                let mut error = Rational::zero();
                for side in &boundary.sides {
                    let uv = side.parameters.evaluate_box([&t, &t], BITS, dec)?;
                    let composed = match side.face {
                        JunctionFace::Strip(end) => strip(model, chains, junction.source, end)?
                            .fillet
                            .surface()
                            .source()
                            .evaluate_box(
                                [&uv[0].lo(), &uv[0].hi()],
                                [&uv[1].lo(), &uv[1].hi()],
                                BITS,
                                dec,
                            )?,
                        JunctionFace::Source(face) => model
                            .geom()
                            .surface(model.face(face).surface)
                            .evaluate(&uv[0], &uv[1], dec)?,
                        JunctionFace::Patch(index) => junction
                            .patches
                            .get(index)
                            .ok_or_else(|| bad("boundary names a missing vertex patch"))?
                            .surface
                            .evaluate(&uv[0], &uv[1], dec)?,
                    };
                    error = error.max((0..3).fold(Rational::zero(), |sum, i| {
                        let delta = &world[i] - &composed[i];
                        sum + delta.lo().abs().max(delta.hi().abs())
                    }));
                }
                if t.is_zero() || t == Rational::one() {
                    let endpoint = usize::from(!t.is_zero());
                    let contact = junction
                        .contacts
                        .get(boundary.contacts[endpoint])
                        .ok_or_else(|| bad("boundary names a missing canonical contact"))?;
                    error = error.max(upper_distance(&world, &contact.position));
                }
                if error > *tolerance {
                    return Err(bad(
                        "assembled boundary misses its world/native or canonical-contact tolerance",
                    ));
                }
                samples.push(NumericalSample {
                    parameters: vec![t],
                    distance_upper_mm: error,
                });
            }
        }
        junction.checks.push(NumericalCheck::Samples {
            name: "fillet.vertex.assembled-world-and-native-boundary".into(),
            scope: "nine world/native composition samples per retained span and canonical endpoints; not a Hausdorff bound".into(),
            samples,
        });
    }
    Ok(())
}

pub(super) fn ball(
    model: &Model,
    chains: &[RolledChain],
    vertex: VertexId,
    ends: &[StripEnd; 3],
    radius: &Rational,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<VertexJunction, Refusal> {
    let candidate = super::super::corner(model, chains, vertex, ends, radius, tolerance, dec)?;
    let mut boundaries = Vec::with_capacity(3);
    for section in &candidate.sections {
        let geometry = super::super::boundaries::construct(
            &chains[section.strip.chain],
            &candidate,
            section,
            tolerance,
            dec,
        )?;
        let [strip, sphere] = geometry.parameters;
        boundaries.push(JunctionBoundary {
            curve: geometry.curve,
            contacts: geometry.contacts.map(|contact| contact.contact),
            sides: [
                BoundarySide {
                    face: JunctionFace::Strip(section.strip),
                    parameters: strip,
                },
                BoundarySide {
                    face: JunctionFace::Patch(0),
                    parameters: sphere,
                },
            ],
        });
    }
    let contacts = candidate
        .ball
        .contacts()
        .iter()
        .enumerate()
        .map(|(i, position)| VertexContact {
            position: position.clone(),
            parameters: vec![(candidate.faces[i], candidate.ball.parameters()[i].clone())],
        })
        .collect();
    let mut junction = VertexJunction {
        source: vertex,
        contacts,
        boundaries,
        terminations: Vec::with_capacity(3),
        patches: vec![JunctionPatch {
            surface: candidate.surface,
            sense: candidate.sense,
        }],
        checks: candidate.checks,
    };
    for (index, section) in candidate.sections.into_iter().enumerate() {
        bind_termination(
            model,
            chains,
            section.strip,
            [section.tau.clone(), section.tau],
            section.contacts,
            &[index],
            &mut junction,
            tolerance,
            dec,
        )?;
    }
    check_boundaries(model, chains, &mut junction, tolerance, dec)?;
    Ok(junction)
}

pub(super) fn cliff(
    model: &Model,
    chains: &[RolledChain],
    vertex: VertexId,
    section: super::cliff::CliffSection,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<VertexJunction, Refusal> {
    let super::cliff::CliffSection {
        ends,
        common,
        caps: [_, mut deep],
        cut,
    } = section;
    if common.iter().any(|side| *side > 1) || deep.source != vertex {
        return Err(bad(
            "Cliff section has inconsistent source or support identities",
        ));
    }
    let first = strip(model, chains, vertex, ends[0])?;
    let second = strip(model, chains, vertex, ends[1])?;
    if first.source.faces[common[0]] != second.source.faces[common[1]] {
        return Err(bad(
            "Cliff section strips do not identify the same common support",
        ));
    }
    if deep.terminations.len() != 1 || deep.boundaries.len() != 1 {
        return Err(bad(
            "Cliff neighbor cap does not identify one complete strip end",
        ));
    }
    let cap_term = deep
        .terminations
        .pop()
        .ok_or_else(|| bad("Cliff cap termination is missing"))?;
    if !same_end(cap_term.strip, ends[1]) {
        return Err(bad("Cliff cap termination belongs to a different strip"));
    }
    let mut contacts = Vec::with_capacity(3);
    for endpoint in 0..2 {
        let t = Rational::from_i64(endpoint as i64);
        let position = cut
            .curve()
            .evaluate_box([&t, &t], BITS, dec)?
            .map(|value| value.mid());
        contacts.push(VertexContact {
            position,
            parameters: Vec::new(),
        });
    }
    let outer = 1 - common[1];
    let cap_contact = cap_term.contacts[outer];
    if cap_contact >= deep.contacts.len() {
        return Err(bad("Cliff cap outer contact is missing"));
    }
    let q = deep.contacts.swap_remove(cap_contact);
    let middle = 1 - common[0];
    let separation = (0..3).fold(Rational::zero(), |sum, i| {
        sum + (&q.position[i] - &contacts[middle].position[i]).abs()
    });
    let unequal = separation > *tolerance;
    let q_index = if unequal {
        contacts.push(q);
        2
    } else {
        middle
    };
    let mut boundaries = Vec::with_capacity(if unequal { 2 } else { 1 });
    boundaries.push(JunctionBoundary {
        curve: cut.curve().clone(),
        contacts: [0, 1],
        sides: [
            BoundarySide {
                face: JunctionFace::Strip(ends[0]),
                parameters: cut.parameters()[0].clone(),
            },
            BoundarySide {
                face: JunctionFace::Strip(ends[1]),
                parameters: cut.parameters()[1].clone(),
            },
        ],
    });
    let t_common = Rational::from_i64(common[0] as i64);
    let common_uv = cut.parameters()[1].evaluate_box([&t_common, &t_common], BITS, dec)?;
    let [outer0, outer1] = cap_term.tau;
    let outer_tau = if outer == 0 { outer0 } else { outer1 };
    let mut deep_tau = [common_uv[0].mid(), outer_tau];
    let mut deep_contacts = [common[0], q_index];
    if common[1] == 1 {
        deep_tau.swap(0, 1);
        deep_contacts.swap(0, 1);
    }
    if unequal {
        let t_middle = Rational::from_i64(middle as i64);
        let native = cut.parameters()[1].evaluate_box([&t_middle, &t_middle], BITS, dec)?;
        let split = native[1].mid();
        let zero = Rational::zero();
        let one = Rational::one();
        let range = if outer == 0 {
            [&zero, &split]
        } else {
            [&split, &one]
        };
        let contact_order = if outer == 0 {
            [q_index, middle]
        } else {
            [middle, q_index]
        };
        let cap = deep
            .boundaries
            .pop()
            .ok_or_else(|| bad("Cliff deeper neighbor boundary is missing"))?;
        if !matches!(cap.sides[0].face, JunctionFace::Strip(end) if same_end(end, ends[1]))
            || !matches!(cap.sides[1].face, JunctionFace::Source(face) if face == first.source.faces[1 - common[0]])
        {
            return Err(bad(
                "Cliff deeper cap does not meet the shallower outer source face",
            ));
        }
        boundaries.push(restrict_boundary(&cap, range, contact_order, dec)?);
    }
    let mut junction = VertexJunction {
        source: vertex,
        contacts,
        boundaries,
        terminations: Vec::with_capacity(2),
        patches: Vec::new(),
        checks: vec![cut.check().clone()],
    };
    bind_termination(
        model,
        chains,
        ends[0],
        cut.tau_ends().clone(),
        [0, 1],
        &[0],
        &mut junction,
        tolerance,
        dec,
    )?;
    let indices: &[usize] = if unequal { &[0, 1] } else { &[0] };
    bind_termination(
        model,
        chains,
        ends[1],
        deep_tau,
        deep_contacts,
        indices,
        &mut junction,
        tolerance,
        dec,
    )?;
    check_boundaries(model, chains, &mut junction, tolerance, dec)?;
    Ok(junction)
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, with_operation_budget};
    use catallax_shapes::{Pose, Surface};

    #[test]
    fn box_ball_junction_retains_actual_strip_and_sphere_boundaries() {
        let ((), work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let mut model = Model::default();
            let size = [20, 30, 40].map(Rational::from_i64);
            let body =
                catallax_topology::block(&mut model, size.clone(), &Pose::identity()).unwrap();
            let vertex = model.vertices_of_body(body)[0];
            let selected: Vec<_> = model
                .edges_of_body(body)
                .into_iter()
                .filter(|edge| {
                    let edge = model.edge(*edge);
                    edge.start == Some(vertex) || edge.end == Some(vertex)
                })
                .collect();
            assert_eq!(selected.len(), 3);
            let radius = Rational::one();
            let tolerance = Rational::parse("0.0000001").unwrap();
            let chains = super::super::super::super::rolling::construct(
                &model, body, &selected, &radius, &tolerance, &mut dec,
            )
            .unwrap();
            assert_eq!(chains.len(), 3);
            let ends: [StripEnd; 3] = std::array::from_fn(|chain| {
                let edge = model.edge(chains[chain].strips[0].source.edge);
                StripEnd {
                    chain,
                    strip: 0,
                    source_end: usize::from(edge.start != Some(vertex)),
                }
            });
            let junction = ball(
                &model, &chains, vertex, &ends, &radius, &tolerance, &mut dec,
            )
            .unwrap();
            let patch = &junction.patches[0];
            let Surface::Sphere(sphere) = &patch.surface else {
                panic!("three convex strips need their checked ball surface")
            };
            assert_eq!(patch.sense, Sense::Same);
            let source = model.geom().position(model.vertex(vertex).point);
            let source = [&source.x, &source.y, &source.z]
                .map(|value| dec.enclose_budgeted_at(value, BITS).unwrap().mid());
            let center = [&sphere.center.x, &sphere.center.y, &sphere.center.z]
                .map(|value| dec.enclose_budgeted_at(value, BITS).unwrap());
            for axis in 0..3 {
                let expected = if source[axis].is_zero() {
                    radius.clone()
                } else {
                    &size[axis] - &radius
                };
                assert!(
                    (&center[axis].lo() - &expected)
                        .abs()
                        .max((&center[axis].hi() - &expected).abs())
                        <= tolerance
                );
            }
            for boundary in &junction.boundaries {
                for index in 0..=31 {
                    let t = Rational::from_ratio(index, 31).unwrap();
                    let world = boundary.curve.evaluate_box([&t, &t], BITS, &dec).unwrap();
                    let radial: [Ball; 3] = std::array::from_fn(|i| &world[i] - &center[i]);
                    let square = radial.iter().fold(
                        Ball::from_rational(&Rational::zero(), BITS),
                        |sum, value| &sum + &value.square(),
                    );
                    let distance = square.sqrt().unwrap();
                    assert!(
                        (&distance.lo() - &radius)
                            .abs()
                            .max((&distance.hi() - &radius).abs())
                            <= tolerance
                    );
                    for side in &boundary.sides {
                        let uv = side.parameters.evaluate_box([&t, &t], BITS, &dec).unwrap();
                        let composed = match side.face {
                            JunctionFace::Strip(end) => chains[end.chain].strips[end.strip]
                                .fillet
                                .surface()
                                .source()
                                .evaluate_box(
                                    [&uv[0].lo(), &uv[0].hi()],
                                    [&uv[1].lo(), &uv[1].hi()],
                                    BITS,
                                    &dec,
                                )
                                .unwrap(),
                            JunctionFace::Patch(index) => junction.patches[index]
                                .surface
                                .evaluate(&uv[0], &uv[1], &mut dec)
                                .unwrap(),
                            JunctionFace::Source(_) => {
                                panic!("ball boundary must join a strip to the ball")
                            }
                        };
                        assert!(super::super::gap(&world, &composed).unwrap() <= tolerance);
                    }
                }
            }
        });
        assert!(!work.exhausted);
        println!(
            "box ball VertexJunction geometry work={}; no body publication",
            work.used
        );
    }
}
