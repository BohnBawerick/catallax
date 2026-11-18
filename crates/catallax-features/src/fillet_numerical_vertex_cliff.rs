//! OCCT-style Cliff carrier sections. The actual fitted rails
//! remain authoritative; a failed finite-domain solve never extrapolates them.

use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_shapes::{NumericalSurfaceCut, SplineSupport, Surface};
use catallax_topology::{BodyId, Model, VertexId};

use super::super::super::{at_stage, rolling::RolledChain};
use super::{BITS, StripEnd, VertexJunction, bad, ends, gap, strip};

pub(crate) fn common_contact(
    model: &Model,
    chains: &[RolledChain],
    vertex: VertexId,
    ends: [StripEnd; 2],
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    if tolerance.sign() != Sign::Positive {
        return Err(bad("Cliff contact tolerance must be positive"));
    }
    let sources = [
        strip(model, chains, vertex, ends[0])?,
        strip(model, chains, vertex, ends[1])?,
    ];
    let common = sources[0]
        .source
        .faces
        .iter()
        .enumerate()
        .find_map(|(i, face)| {
            sources[1]
                .source
                .faces
                .iter()
                .position(|other| other == face)
                .map(|j| [i, j])
        })
        .ok_or_else(|| bad("Cliff strips have no common source support"))?;
    let rails = [
        sources[0]
            .fillet
            .rail(common[0])
            .ok_or_else(|| bad("first Cliff contact rail is missing"))?,
        sources[1]
            .fillet
            .rail(common[1])
            .ok_or_else(|| bad("second Cliff contact rail is missing"))?,
    ];
    let domains = rails.map(|rail| rail.domain().map(Rational::to_f64_nearest));
    if domains.iter().flatten().any(|value| !value.is_finite()) {
        return Err(bad("Cliff rail domain is outside binary64 proposal range"));
    }
    let mut parameters = ends.map(|end| end.source_end as f64);
    let limit = tolerance.mul_pow2(-3);
    for _ in 0..32 {
        dec.charge_work(1)?;
        let at = |side: usize| -> Result<Rational, Refusal> {
            let q = Rational::from_f64_exact(parameters[side])
                .ok_or_else(|| bad("non-finite Cliff contact proposal"))?;
            let domain = rails[side].domain();
            Ok(q.max(domain[0].clone()).min(domain[1].clone()))
        };
        let tau = [at(0)?, at(1)?];
        let jets = [
            rails[0].jet_box([&tau[0], &tau[0]], BITS, dec)?,
            rails[1].jet_box([&tau[1], &tau[1]], BITS, dec)?,
        ];
        let error = gap(&jets[0].value, &jets[1].value)?;
        if error <= limit {
            return Ok(tau);
        }
        let residual: [f64; 3] = std::array::from_fn(|i| {
            jets[0].value[i].mid().to_f64_nearest() - jets[1].value[i].mid().to_f64_nearest()
        });
        let a = jets[0]
            .first
            .each_ref()
            .map(|value| value.mid().to_f64_nearest());
        let b = jets[1]
            .first
            .each_ref()
            .map(|value| value.mid().to_f64_nearest());
        let dot = |a: &[f64; 3], b: &[f64; 3]| (0..3).map(|i| a[i] * b[i]).sum::<f64>();
        let aa = dot(&a, &a);
        let bb = dot(&b, &b);
        let ab = -dot(&a, &b);
        let determinant = aa * bb - ab * ab;
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err(bad(
                "Cliff contact rails have an unresolved intersection direction",
            ));
        }
        let rhs = [-dot(&a, &residual), dot(&b, &residual)];
        let step = [
            (bb * rhs[0] - ab * rhs[1]) / determinant,
            (aa * rhs[1] - ab * rhs[0]) / determinant,
        ];
        if step.iter().any(|value| !value.is_finite()) {
            return Err(bad("Cliff contact solve produced a non-finite step"));
        }
        let next: [f64; 2] =
            std::array::from_fn(|i| (parameters[i] + step[i]).clamp(domains[i][0], domains[i][1]));
        if (0..2).all(|i| next[i].to_bits() == parameters[i].to_bits()) {
            return Err(bad(format!(
                "common rail contact stalled in finite strip domains: tau={parameters:?}, domains={domains:?}, world separation={} mm, points={:?}",
                error.to_f64_nearest(),
                jets.each_ref().map(|jet| jet
                    .value
                    .each_ref()
                    .map(|value| value.mid().to_f64_nearest())),
            )));
        }
        parameters = next;
    }
    Err(bad(format!(
        "common rail contact did not converge in finite strip domains: tau={parameters:?}, domains={domains:?}"
    )))
}

pub(crate) struct CliffSection {
    /// Shallower strip first. Native curves in `cut` follow this order.
    pub ends: [StripEnd; 2],
    pub common: [usize; 2],
    pub caps: [VertexJunction; 2],
    pub cut: NumericalSurfaceCut,
}

/// Intersect the actual carriers after their outer-face cuts establish reach.
/// Unequal reach requires retaining an additional piece of the deeper cap.
pub(crate) fn shared_section(
    model: &Model,
    body: BodyId,
    chains: &[RolledChain],
    vertex: VertexId,
    mut strip_ends: [StripEnd; 2],
    mut common_tau: [Rational; 2],
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<CliffSection, Refusal> {
    let mut sources = [
        strip(model, chains, vertex, strip_ends[0])?,
        strip(model, chains, vertex, strip_ends[1])?,
    ];
    if !sources[0].convex || !sources[1].convex {
        return Err(bad(
            "this Cliff section needs two convex strips over the pocket corner",
        ));
    }
    let mut common = sources[0]
        .source
        .faces
        .iter()
        .enumerate()
        .find_map(|(i, face)| {
            sources[1]
                .source
                .faces
                .iter()
                .position(|other| other == face)
                .map(|j| [i, j])
        })
        .ok_or_else(|| bad("Cliff section has no common source support"))?;
    let mut caps = [
        ends::face_end(model, body, chains, vertex, strip_ends[0], tolerance, dec)
            .map_err(|error| at_stage("first Cliff neighbor section", error))?,
        ends::face_end(model, body, chains, vertex, strip_ends[1], tolerance, dec)
            .map_err(|error| at_stage("second Cliff neighbor section", error))?,
    ];
    let point = model.geom().position(model.vertex(vertex).point);
    let near = [
        dec.enclose_budgeted_at(&point.x, BITS)?.mid(),
        dec.enclose_budgeted_at(&point.y, BITS)?.mid(),
        dec.enclose_budgeted_at(&point.z, BITS)?.mid(),
    ];
    let reach = |index: usize| {
        let contact = &caps[index].contacts[1 - common[index]].position;
        (0..3).fold(Rational::zero(), |sum, i| {
            let delta = &contact[i] - &near[i];
            sum + (&delta * &delta)
        })
    };
    // A numerical ordering proposal. The retained curves and rail contacts
    // must still satisfy the actual world/native composition tolerance.
    if reach(1) < reach(0) {
        sources.swap(0, 1);
        strip_ends.swap(0, 1);
        common_tau.swap(0, 1);
        common.swap(0, 1);
        caps.swap(0, 1);
    }
    let seed = if common[0] == 0 {
        [
            common_tau[0].clone(),
            common_tau[1].clone(),
            Rational::from_i64(common[1] as i64),
        ]
    } else {
        [
            caps[0].terminations[0].tau[0].clone(),
            caps[1].terminations[0].tau[1 - common[1]].clone(),
            Rational::from_i64((1 - common[1]) as i64),
        ]
    };
    let other = Surface::BSpline(SplineSupport::RationalAtlas(
        sources[1].fillet.surface().clone(),
    ));
    let cut = sources[0]
        .fillet
        .surface_cut(other, seed, dec)
        .map_err(|error| at_stage("two fitted Cliff carriers", error))?;
    let at = Rational::from_i64(common[0] as i64);
    let world = cut.curve().evaluate_box([&at, &at], BITS, dec)?;
    for side in 0..2 {
        let rail = sources[side]
            .fillet
            .rail(common[side])
            .ok_or_else(|| bad("Cliff common rail is missing"))?;
        let native = cut.parameters()[side].evaluate_box([&at, &at], BITS, dec)?;
        let tau = native[0].mid();
        let contact = rail.evaluate_box([&tau, &tau], BITS, dec)?;
        if gap(&world, &contact)? > *tolerance {
            return Err(bad(
                "Cliff carrier intersection leaves a common contact rail",
            ));
        }
    }
    Ok(CliffSection {
        ends: strip_ends,
        common,
        caps,
        cut,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, with_operation_budget};
    use catallax_recipe::{Recipe, Role};

    #[test]
    fn owner_vertex_cylindrical_cliff_section_has_paired_native_curves() {
        let ((), work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let recipe =
                Recipe::parse(include_str!("../tests/fixtures/owner-curved-pocket.json")).unwrap();
            let built = crate::evaluate(&recipe, &mut dec).unwrap();
            let model = &built.model;
            let selected: Vec<_> = model.edges_of_body(built.result).into_iter().filter(|edge| {
                model.edge_name(built.result, *edge).is_some_and(|name| {
                    name.between().iter().any(|origin| origin.step.as_str() == "blank"
                        && matches!(&origin.role, Role::Side { segment, .. } if segment.as_str() == "right-round"))
                        && name.between().iter().any(|origin| origin.step.as_str() == "pocket-tool"
                            && matches!(&origin.role, Role::Side { segment, .. } if ["far", "right"].contains(&segment.as_str())))
                })
            }).collect();
            assert_eq!(selected.len(), 2);
            let tolerance = Rational::parse("0.0000001").unwrap();
            let chains = super::super::super::super::rolling::construct(
                model,
                built.result,
                &selected,
                &Rational::one(),
                &tolerance,
                &mut dec,
            )
            .unwrap();
            assert_eq!(chains.len(), 2);
            let a = model.edge(chains[0].strips[0].source.edge);
            let b = model.edge(chains[1].strips[0].source.edge);
            let vertex = [a.start, a.end]
                .into_iter()
                .flatten()
                .find(|vertex| b.start == Some(*vertex) || b.end == Some(*vertex))
                .unwrap();
            let ends = [a, b].map(|edge| usize::from(edge.start != Some(vertex)));
            let ends = std::array::from_fn(|i| StripEnd {
                chain: i,
                strip: 0,
                source_end: ends[i],
            });
            let tau = common_contact(model, &chains, vertex, ends, &tolerance, &mut dec).unwrap();
            for side in 0..2 {
                let rail = chains[side].strips[0].fillet.rail(0).unwrap();
                let point = rail
                    .evaluate_box([&tau[side], &tau[side]], BITS, &dec)
                    .unwrap();
                for (axis, expected) in [(0, "33.75"), (1, "-25.665")] {
                    let expected = Rational::parse(expected).unwrap();
                    assert!(
                        (&point[axis].lo() - &expected)
                            .abs()
                            .max((&point[axis].hi() - &expected).abs())
                            <= tolerance
                    );
                }
            }
            let section = shared_section(
                model,
                built.result,
                &chains,
                vertex,
                ends,
                tau,
                &tolerance,
                &mut dec,
            )
            .unwrap();
            for index in 0..=31 {
                let s = Rational::from_ratio(index, 31).unwrap();
                let world = section
                    .cut
                    .curve()
                    .evaluate_box([&s, &s], BITS, &dec)
                    .unwrap();
                for side in 0..2 {
                    let end = section.ends[side];
                    let uv = section.cut.parameters()[side]
                        .evaluate_box([&s, &s], BITS, &dec)
                        .unwrap();
                    let surface = chains[end.chain].strips[end.strip]
                        .fillet
                        .surface()
                        .source();
                    let composed = surface
                        .evaluate_box(
                            [&uv[0].lo(), &uv[0].hi()],
                            [&uv[1].lo(), &uv[1].hi()],
                            BITS,
                            &dec,
                        )
                        .unwrap();
                    assert!(
                        gap(&world, &composed).unwrap() <= tolerance,
                        "actual Cliff section at {s}"
                    );
                }
            }
            let deep_end = section.ends[1];
            let junction = super::super::assembly::cliff(
                model, &chains, vertex, section, &tolerance, &mut dec,
            )
            .unwrap();
            assert!(
                junction.patches.is_empty(),
                "Cliff must not substitute a sphere or Smooth patch"
            );
            let deep = junction
                .terminations
                .iter()
                .find(|termination| {
                    termination.strip.chain == deep_end.chain
                        && termination.strip.strip == deep_end.strip
                        && termination.strip.source_end == deep_end.source_end
                })
                .unwrap();
            assert_ne!(
                deep.tau[0], deep.tau[1],
                "unequal reach needs distinct rail parameters"
            );
            let mut previous = junction.contacts[deep.contacts[0]].position.clone();
            for &(index, sense) in &deep.boundaries {
                let boundary = &junction.boundaries[index];
                let zero = Rational::zero();
                let one = Rational::one();
                let (a, b) = if sense == catallax_topology::Sense::Same {
                    (&zero, &one)
                } else {
                    (&one, &zero)
                };
                let start = boundary.curve.evaluate_box([a, a], BITS, &dec).unwrap();
                assert!(super::super::upper_distance(&start, &previous) <= tolerance);
                previous = boundary
                    .curve
                    .evaluate_box([b, b], BITS, &dec)
                    .unwrap()
                    .map(|value| value.mid());
                for j in 0..=31 {
                    let t = Rational::from_ratio(j, 31).unwrap();
                    let world = boundary.curve.evaluate_box([&t, &t], BITS, &dec).unwrap();
                    for side in &boundary.sides {
                        let native = side
                            .parameters
                            .point_at(&catallax_numbers::Real::from_rational(t.clone()), &mut dec)
                            .unwrap();
                        let uv = [
                            dec.enclose_budgeted_at(&native.x, BITS).unwrap(),
                            dec.enclose_budgeted_at(&native.y, BITS).unwrap(),
                        ];
                        let composed = match side.face {
                            super::super::JunctionFace::Strip(end) => chains[end.chain].strips
                                [end.strip]
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
                            super::super::JunctionFace::Source(face) => model
                                .geom()
                                .surface(model.face(face).surface)
                                .evaluate(&uv[0], &uv[1], &mut dec)
                                .unwrap(),
                            super::super::JunctionFace::Patch(_) => {
                                panic!("Cliff termination should not invent a patch")
                            }
                        };
                        assert!(gap(&world, &composed).unwrap() <= tolerance);
                    }
                }
            }
            let target = &junction.contacts[deep.contacts[1]].position;
            assert!(
                (0..3).fold(Rational::zero(), |sum, i| sum
                    + (&previous[i] - &target[i]).abs())
                    <= tolerance
            );
        });
        assert!(!work.exhausted);
        println!(
            "owner cylindrical Cliff section work={}; actual fitted carriers and native curves, no extrapolation",
            work.used
        );
    }
    #[test]
    fn owner_rim_vertex_junctions_at_half_one_and_two_mm() {
        for radius in ["0.5", "1", "2"] {
            let ((), work) = with_operation_budget(100_000_000, || {
                let mut dec = Decider::new(BudgetPolicy {
                    max_work: 100_000_000,
                    ..Default::default()
                });
                let recipe =
                    Recipe::parse(include_str!("../tests/fixtures/owner-curved-pocket.json"))
                        .unwrap();
                let built = crate::evaluate(&recipe, &mut dec).unwrap();
                let model = &built.model;
                let selected: Vec<_> = model
                    .edges_of_body(built.result)
                    .into_iter()
                    .filter(|edge| {
                        model.edge_name(built.result, *edge).is_some_and(|name| {
                            name.between()
                                .iter()
                                .any(|origin| origin.step.as_str() == "blank")
                                && name
                                    .between()
                                    .iter()
                                    .any(|origin| origin.step.as_str() == "pocket-tool")
                        })
                    })
                    .collect();
                assert_eq!(selected.len(), 6, "the actual six-edge owner rim");
                let radius = Rational::parse(radius).unwrap();
                let tolerance = Rational::parse("0.0000001").unwrap();
                let chains = super::super::super::super::rolling::construct(
                    model,
                    built.result,
                    &selected,
                    &radius,
                    &tolerance,
                    &mut dec,
                )
                .unwrap();
                let junctions = super::super::construct(
                    model,
                    built.result,
                    &chains,
                    &radius,
                    &tolerance,
                    &mut dec,
                )
                .unwrap();
                for (chain, value) in chains.iter().enumerate() {
                    for strip in 0..value.strips.len() {
                        for source_end in 0..2 {
                            let count = junctions
                                .iter()
                                .flat_map(|junction| &junction.terminations)
                                .filter(|end| {
                                    end.strip.chain == chain
                                        && end.strip.strip == strip
                                        && end.strip.source_end == source_end
                                })
                                .count();
                            assert_eq!(
                                count, 1,
                                "R{radius}: every actual strip end needs one geometric termination"
                            );
                        }
                    }
                }
                for junction in &junctions {
                    assert!(
                        junction.patches.is_empty(),
                        "the owner Cliff rim must not acquire a substitute sphere or Smooth patch"
                    );
                    for boundary in &junction.boundaries {
                        for i in 1..7 {
                            let t = Rational::from_ratio(i, 7).unwrap();
                            let world = boundary.curve.evaluate_box([&t, &t], BITS, &dec).unwrap();
                            for side in &boundary.sides {
                                let native = side
                                    .parameters
                                    .point_at(
                                        &catallax_numbers::Real::from_rational(t.clone()),
                                        &mut dec,
                                    )
                                    .unwrap();
                                let uv = [
                                    dec.enclose_budgeted_at(&native.x, BITS).unwrap(),
                                    dec.enclose_budgeted_at(&native.y, BITS).unwrap(),
                                ];
                                let composed = match side.face {
                                    super::super::JunctionFace::Strip(end) => chains[end.chain]
                                        .strips[end.strip]
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
                                    super::super::JunctionFace::Source(face) => model
                                        .geom()
                                        .surface(model.face(face).surface)
                                        .evaluate(&uv[0], &uv[1], &mut dec)
                                        .unwrap(),
                                    super::super::JunctionFace::Patch(_) => {
                                        panic!("unexpected owner-rim patch")
                                    }
                                };
                                assert!(
                                    gap(&world, &composed).unwrap() <= tolerance,
                                    "R{radius}: retained world/native boundary mismatch"
                                );
                            }
                        }
                    }
                }
            });
            assert!(!work.exhausted);
            println!(
                "owner six-edge rim R{radius} vertex geometry work={}; no body publication",
                work.used
            );
        }
    }
}
