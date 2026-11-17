//! Actual neighboring-face cuts, following OCCT's one-stripe end case and
//! SMLib's Nx1 end-rib construction. Face-domain publication is separate.

use catallax_numbers::{Ball, Decider, NumericalCheck, NumericalSample, Rational, Refusal, Sign};
use catallax_topology::{BodyId, Loop, Model, Sense, VertexId};

use super::super::super::rolling::RolledChain;
use super::{
    BITS, BoundarySide, JunctionBoundary, JunctionFace, StripEnd, StripTermination, VertexContact,
    VertexJunction, bad, strip, upper_distance,
};

pub(crate) fn face_end(
    model: &Model,
    body: BodyId,
    chains: &[RolledChain],
    vertex: VertexId,
    end: StripEnd,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<VertexJunction, Refusal> {
    if tolerance.sign() != Sign::Positive {
        return Err(bad("end-section tolerance must be positive"));
    }
    let source = strip(model, chains, vertex, end)?;
    let mut neighbor = None;
    for face in model.faces_of_body(body) {
        if source.source.faces.contains(&face) {
            continue;
        }
        dec.charge_work(1)?;
        let incident = model
            .face(face)
            .loops
            .iter()
            .any(|boundary| match model.loop_(*boundary) {
                Loop::Vertex(at) => *at == vertex,
                Loop::Coedges(uses) => uses.iter().any(|coedge| {
                    let edge = model.edge(model.coedge(*coedge).edge);
                    edge.start == Some(vertex) || edge.end == Some(vertex)
                }),
            });
        if incident && neighbor.replace(face).is_some() {
            return Err(bad("an end section needs one non-support source face"));
        }
    }
    let neighbor = neighbor.ok_or_else(|| bad("an end section has no neighboring source face"))?;
    let surface = model.geom().surface(model.face(neighbor).surface);
    let point = model.geom().position(model.vertex(vertex).point);
    let seed = crate::numerical_embedding::surface_parameters(surface, &point, tolerance, dec)?;
    let cut = source.fillet.surface_cut(
        surface.clone(),
        [
            Rational::from_i64(end.source_end as i64),
            seed[0].clone(),
            seed[1].clone(),
        ],
        dec,
    )?;
    let mut contacts = Vec::with_capacity(2);
    let mut samples = Vec::with_capacity(2);
    for side in 0..2 {
        let s = Rational::from_i64(side as i64);
        let position = cut
            .curve()
            .evaluate_box([&s, &s], BITS, dec)?
            .map(|value| value.mid());
        let tau = &cut.tau_ends()[side];
        let rail = source
            .fillet
            .rail(side)
            .ok_or_else(|| bad("end contact has no world rail"))?;
        let mut error = upper_distance(&rail.evaluate_box([tau, tau], BITS, dec)?, &position);
        let native = source
            .fillet
            .parameter_rail(side)
            .ok_or_else(|| bad("end contact has no original-support rail"))?
            .evaluate_box([tau, tau], BITS, dec)?;
        let support_uv = [native[0].mid(), native[1].mid()];
        let uv = support_uv
            .each_ref()
            .map(|value| Ball::from_rational(value, BITS));
        let face = source.source.faces[side];
        error = error.max(upper_distance(
            &model
                .geom()
                .surface(model.face(face).surface)
                .evaluate(&uv[0], &uv[1], dec)?,
            &position,
        ));
        let native = cut.parameters()[1].evaluate_box([&s, &s], BITS, dec)?;
        let neighbor_uv = [native[0].mid(), native[1].mid()];
        let uv = neighbor_uv
            .each_ref()
            .map(|value| Ball::from_rational(value, BITS));
        error = error.max(upper_distance(
            &surface.evaluate(&uv[0], &uv[1], dec)?,
            &position,
        ));
        if error > *tolerance {
            return Err(bad(
                "end-section contact misses its rail, support or neighbor",
            ));
        }
        contacts.push(VertexContact {
            position,
            parameters: vec![(face, support_uv), (neighbor, neighbor_uv)],
        });
        samples.push(NumericalSample {
            parameters: vec![tau.clone(), s],
            distance_upper_mm: error,
        });
    }
    Ok(VertexJunction {
        source: vertex,
        contacts,
        boundaries: vec![JunctionBoundary {
            curve: cut.curve().clone(),
            contacts: [0, 1],
            sides: [
                BoundarySide {
                    face: JunctionFace::Strip(end),
                    parameters: cut.parameters()[0].clone(),
                },
                BoundarySide {
                    face: JunctionFace::Source(neighbor),
                    parameters: cut.parameters()[1].clone(),
                },
            ],
        }],
        terminations: vec![StripTermination {
            strip: end,
            tau: cut.tau_ends().clone(),
            contacts: [0, 1],
            support_parameters: [0, 0],
            boundaries: vec![(0, Sense::Same)],
        }],
        patches: Vec::new(),
        checks: vec![
            cut.check().clone(),
            NumericalCheck::Samples {
                name: "fillet.vertex.face-end-contacts".into(),
                scope: "two canonical rail contacts on their actual supports and neighbor".into(),
                samples,
            },
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{BudgetPolicy, current_operation_work, with_operation_budget};
    use catallax_recipe::{Recipe, Role};

    #[test]
    fn owner_vertex_sloping_face_end_keeps_different_rail_parameters() {
        let ((rolling_work, end_work), work) = with_operation_budget(100_000_000, || {
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
                    [ ("blank", "right-round"), ("pocket-tool", "far") ].iter().all(|(step, segment)| {
                        name.between().iter().any(|origin| origin.step.as_str() == *step
                            && matches!(&origin.role, Role::Side { segment: actual, .. } if actual.as_str() == *segment))
                    })
                })
            }).collect();
            assert_eq!(selected.len(), 1);
            let edge = model.edge(selected[0]);
            let vertex = [edge.start, edge.end]
                .into_iter()
                .flatten()
                .find(|vertex| {
                    model
                        .geom()
                        .position(model.vertex(*vertex).point)
                        .x
                        .as_rational()
                        == Some(Rational::from_i64(32))
                })
                .unwrap();
            let tolerance = Rational::parse("0.0000001").unwrap();
            let before = current_operation_work().used;
            let chains = super::super::super::super::rolling::construct(
                model,
                built.result,
                &selected,
                &Rational::one(),
                &tolerance,
                &mut dec,
            )
            .unwrap();
            let after_rolling = current_operation_work().used;
            let section = face_end(
                model,
                built.result,
                &chains,
                vertex,
                StripEnd {
                    chain: 0,
                    strip: 0,
                    source_end: usize::from(edge.start != Some(vertex)),
                },
                &tolerance,
                &mut dec,
            )
            .unwrap();
            let after_end = current_operation_work().used;
            let termination = &section.terminations[0];
            assert_ne!(
                termination.tau[0], termination.tau[1],
                "a sloping end is not a fixed-U section"
            );
            let boundary = &section.boundaries[0];
            let JunctionFace::Source(neighbor) = boundary.sides[1].face else {
                panic!("end must retain its actual neighbor")
            };
            let source = &chains[0].strips[0];
            assert!(!source.source.faces.contains(&neighbor));
            assert!(section.patches.is_empty());
            for index in 0..=31 {
                let s = Rational::from_ratio(index, 31).unwrap();
                let world = boundary.curve.evaluate_box([&s, &s], BITS, &dec).unwrap();
                let uv = boundary.sides[0]
                    .parameters
                    .evaluate_box([&s, &s], BITS, &dec)
                    .unwrap();
                let strip = source
                    .fillet
                    .surface()
                    .source()
                    .evaluate_box(
                        [&uv[0].lo(), &uv[0].hi()],
                        [&uv[1].lo(), &uv[1].hi()],
                        BITS,
                        &dec,
                    )
                    .unwrap();
                assert!(super::super::gap(&world, &strip).unwrap() <= tolerance);
                let uv = boundary.sides[1]
                    .parameters
                    .evaluate_box([&s, &s], BITS, &dec)
                    .unwrap();
                let neighbor = model
                    .geom()
                    .surface(model.face(neighbor).surface)
                    .evaluate(&uv[0], &uv[1], &mut dec)
                    .unwrap();
                assert!(super::super::gap(&world, &neighbor).unwrap() <= tolerance);
            }
            (after_rolling - before, after_end - after_rolling)
        });
        assert!(!work.exhausted);
        println!(
            "owner sloping face end rolling_work={rolling_work} end_work={end_work} total_work={}; geometry only",
            work.used
        );
    }
}
