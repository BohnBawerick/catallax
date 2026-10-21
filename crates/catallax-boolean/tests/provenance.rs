//! Stage mark order, incidence provenance and budgets across chained results.
mod common;
use common::*;
use std::time::{Duration, Instant};

fn manual(
    m: &mut Model,
    a: BodyId,
    b: BodyId,
    op: BoolOp,
    dec: &mut Decider,
) -> Result<(BodyId, MarkLog), BooleanRefusal> {
    let mut marks = MarkLog::new();
    macro_rules! stage {
        ($stage:expr,$call:expr) => {{
            let answer = $call.map_err(|refusal| BooleanRefusal {
                stage: $stage,
                refusal,
            })?;
            marks.extend(&answer.marks);
            answer.value
        }};
    }
    let graph = stage!(Stage::Intersect, intersect(m, a, b, dec));
    let imprinted = stage!(Stage::Imprint, imprint(m, &graph, dec));
    let classification = stage!(Stage::Classify, classify(m, &imprinted, dec));
    let selection = keep(&classification, &imprinted, op);
    let before = m.geom().facts().clone();
    let body = stage!(Stage::Sew, sew(m, &imprinted, &selection, dec));
    // Any incidence already justified by earlier stages, new intersection
    // construction, or proven equivalence needs no duplicate sew record.
    // All remaining result incidences must be supplied by boolean.sew.
    let mut justified = before;
    for fact in m.geom().facts().iter().copied() {
        if matches!(
            fact.origin,
            FactOrigin::Construction("boolean.sew" | "intersect.surfaces" | "shapes.definition")
        ) || matches!(fact.origin, FactOrigin::Proof(_)) && fact.kind == FactKind::Coincident
        {
            justified.record(fact);
        }
    }
    for f in m.faces_of_body(body) {
        let surface = GeomKey::Surface(m.face(f).surface);
        for &l in &m.face(f).loops {
            if let Loop::Coedges(cs) = m.loop_(l) {
                for &c in cs {
                    let edge = m.edge(m.coedge(c).edge);
                    let curve = GeomKey::Curve(edge.curve);
                    assert!(
                        justified.holds(FactKind::Incident, curve, surface),
                        "missing boolean.sew provenance"
                    );
                    for v in [edge.start, edge.end].into_iter().flatten() {
                        assert!(
                            justified.holds(
                                FactKind::Incident,
                                GeomKey::Point(m.vertex(v).point),
                                curve
                            ),
                            "missing endpoint provenance"
                        );
                    }
                }
            }
        }
    }
    marks.extend(&m.body(a).marks);
    marks.extend(&m.body(b).marks);
    Ok((body, marks))
}

blind_test!(
    stage_marks_match_whole_operation_in_order_under_both_budgets,
    {
        for &id in IDS {
            for policy in [
                BudgetPolicy::default(),
                BudgetPolicy {
                    max_precision_bits: 64,
                    max_work: 10_000,
                },
            ] {
                for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
                    let c = case(id);
                    let mut whole = Model::new();
                    let a = c.a.build(&mut whole);
                    let b = c.b.build(&mut whole);
                    let mut stages = whole.clone();
                    let complete = boolean(&mut whole, a, b, op, &mut Decider::new(policy));
                    let individual = manual(&mut stages, a, b, op, &mut Decider::new(policy));
                    match (complete, individual) {
                        (Ok(answer), Ok((body, marks))) => {
                            assert_eq!(answer.marks, marks);
                            assert_eq!(whole.body(answer.body).marks, marks);
                            assert_eq!(counts(&whole, answer.body), counts(&stages, body));
                            assert_eq!(
                                topology_snapshot(&whole, &[answer.body]),
                                topology_snapshot(&stages, &[body])
                            );
                            assert_eq!(
                                answer.report.stages.iter().map(|s| s.marks).sum::<usize>(),
                                marks.len()
                            );
                        }
                        (Err(a), Err(b)) => {
                            assert_eq!(a.stage, b.stage);
                            assert_eq!(a.refusal.code, b.refusal.code);
                        }
                        (a, b) => panic!("whole/stage disagreement for {id} {op:?}: {a:?} / {b:?}"),
                    }
                }
            }
        }
    }
);

blind_test!(
    case_59_chain_steps_and_empty_queries_honour_the_caller_budget,
    {
        for number in 41..=44 {
            for expired in [false, true] {
                let spec = chain(number);
                let mut m = Model::new();
                let mut bodies: Vec<_> = spec.operands.iter().map(|p| p.build(&mut m)).collect();
                let mut marked = false;
                for (step, (a, b, op, want)) in spec.steps.iter().enumerate() {
                    let before = topology_snapshot(&m, &bodies);
                    let mut dec = if expired {
                        Decider::default()
                    } else {
                        Decider::new(BudgetPolicy {
                            max_precision_bits: 64,
                            max_work: 10_000,
                        })
                    };
                    if expired {
                        dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
                    }
                    match boolean(&mut m, bodies[*a], bodies[*b], *op, &mut dec) {
                        Err(e) => {
                            if expired {
                                assert_eq!(e.refusal.code, "numbers.deadline");
                            } else {
                                match e.refusal.code {
                                    "boolean.imprint-ambiguous" => assert!(matches!(
                                        e.stage,
                                        Stage::Intersect | Stage::Imprint
                                    )),
                                    "boolean.merge-ambiguous" | "boolean.self-touching-shell" => {
                                        assert_eq!(e.stage, Stage::Sew)
                                    }
                                    "boolean.classify-ambiguous"
                                    | "boolean.no-rational-point"
                                    | "numbers.budget-exhausted"
                                    | "numbers.root-not-certified" => {}
                                    _ => panic!("unexpected chain refusal: {e:?}"),
                                }
                            }
                            assert_eq!(before, topology_snapshot(&m, &bodies));
                            break;
                        }
                        Ok(answer) => {
                            assert_eq!(before, topology_snapshot(&m, &bodies));
                            topology_rules(&m, answer.body);
                            marked |= !answer.marks.is_empty();
                            if !marked {
                                check_counts(&m, answer.body, want);
                                let mut rng = Rng(0x59c4a1);
                                for _ in 0..32 {
                                    let point = rng.point();
                                    if let Some(want) = spec.locate(&point, step) {
                                        query(
                                            &mut m,
                                            answer.body,
                                            &point,
                                            want,
                                            &mut Decider::default(),
                                        );
                                    }
                                }
                            }
                            bodies.push(answer.body);
                        }
                    }
                }
            }
        }
        let mut m = Model::new();
        let body = empty(&mut m);
        let mut dec = Decider::new(BudgetPolicy {
            max_precision_bits: 0,
            max_work: 0,
        });
        dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
        for point in [p("0", "0", "0"), p("1e-30", "1e30", "-1e30")] {
            query(&mut m, body, &point, Location::Outside, &mut dec);
        }
    }
);
