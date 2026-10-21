//! Independently checked section graphs and the reusable five-stage API.
mod common;
use common::*;

blind_test!(case_48_through_hole_graph_has_two_circle_rings, {
    let c = case("10");
    let mut m = Model::new();
    let a = c.a.build(&mut m);
    let b = c.b.build(&mut m);
    let before = topology_snapshot(&m, &[a, b]);
    let mut dec = Decider::default();
    let graph = exact(intersect(&mut m, a, b, &mut dec));
    assert_eq!(graph.first(), a);
    assert_eq!(graph.second(), b);
    assert!(graph.points().is_empty());
    assert_eq!(graph.edges().len(), 2);
    let mut levels = BTreeSet::new();
    for edge in graph.edges() {
        assert!(edge.start.is_none() && edge.end.is_none());
        assert_eq!(edge.faces_of_first.len(), 1);
        assert_eq!(edge.faces_of_second.len(), 1);
        let plane = m.geom().surface(m.face(edge.faces_of_first[0]).surface);
        let side = m.geom().surface(m.face(edge.faces_of_second[0]).surface);
        assert_eq!(plane.kind(), SurfaceKind::Plane);
        assert_eq!(side.kind(), SurfaceKind::Cylinder);
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = curve_form(m.geom(), edge.curve, &mut dec).unwrap()
        else {
            panic!("not a circle")
        };
        assert_eq!(radius_sq, r(9));
        assert!(normal.cross(&v(0, 0, 1)).is_zero());
        let Surface::Plane(plane) = plane else {
            unreachable!()
        };
        let z = plane.origin.z.clone();
        assert!(z == r(0) || z == r(10));
        same_real_point(&center, &RPoint3::new(r(5), r(5), z.clone()));
        levels.insert(z);
    }
    assert_eq!(levels, BTreeSet::from([r(0), r(10)]));
    assert_eq!(before, topology_snapshot(&m, &[a, b]));
});

blind_test!(case_49_tangent_spheres_graph_is_one_isolated_point, {
    let c = case("37");
    let mut m = Model::new();
    let a = c.a.build(&mut m);
    let b = c.b.build(&mut m);
    let graph = exact(intersect(&mut m, a, b, &mut Decider::default()));
    assert!(graph.edges().is_empty());
    assert_eq!(graph.points().len(), 1);
    same_real_point(&m.geom().position(graph.points()[0]), &p("5", "0", "0"));
});

blind_test!(case_50_coplanar_square_graph_has_ordered_face_owners, {
    let c = case("02");
    let mut m = Model::new();
    let a = c.a.build(&mut m);
    let b = c.b.build(&mut m);
    let mut dec = Decider::default();
    let graph = exact(intersect(&mut m, a, b, &mut dec));
    assert!(graph.points().is_empty());
    assert_eq!(graph.edges().len(), 4);
    let corners = [
        p("1", "0", "0"),
        p("1", "1", "0"),
        p("1", "1", "1"),
        p("1", "0", "1"),
    ];
    let index = |m: &Model, id: PointId| {
        let point = m.geom().position(id);
        corners
            .iter()
            .position(|p| {
                [&point.x, &point.y, &point.z]
                    .into_iter()
                    .zip([&p.x, &p.y, &p.z])
                    .all(|(x, y)| {
                        matches!(
                            Decider::default().compare(x, &Real::from(y.clone())),
                            Outcome::ProvenCoincident(_)
                        )
                    })
            })
            .expect("unexpected corner")
    };
    let mut found = BTreeSet::new();
    for edge in graph.edges() {
        assert!(matches!(
            curve_form(m.geom(), edge.curve, &mut dec).unwrap(),
            CurveForm::Line { .. }
        ));
        let i = index(&m, edge.start.unwrap());
        let j = index(&m, edge.end.unwrap());
        found.insert((i.min(j), i.max(j)));
        for (body, owners) in [(a, &edge.faces_of_first), (b, &edge.faces_of_second)] {
            let want: Vec<_> = m
                .faces_of_body(body)
                .into_iter()
                .filter(|&f| {
                    let s = m.geom().surface(m.face(f).surface);
                    implicit_normal(s, &corners[i]).0.is_zero()
                        && implicit_normal(s, &corners[j]).0.is_zero()
                })
                .collect();
            assert_eq!(want.len(), 2);
            assert_eq!(owners, &want);
        }
    }
    assert_eq!(found, BTreeSet::from([(0, 1), (1, 2), (2, 3), (0, 3)]));
});

blind_test!(cases_51_52_imprint_classify_keep_and_direct_sew, {
    for &id in IDS {
        let c = case(id);
        let mut m = Model::new();
        let a = c.a.build(&mut m);
        let b = c.b.build(&mut m);
        let unrelated = ball_at(1, [100; 3]).build(&mut m);
        let unknown = m.faces_of_body(unrelated)[0];
        let before = topology_snapshot(&m, &[a, b, unrelated]);
        let mut dec = Decider::default();
        let graph = exact(intersect(&mut m, a, b, &mut dec));
        let imprinted = exact(imprint(&mut m, &graph, &mut dec));
        let classification = exact(classify(&mut m, &imprinted, &mut dec));
        assert_eq!(classification.class_of(unknown), None);
        let sides = [
            (a, &c.a, &c.b, imprinted.pieces_of_first()),
            (b, &c.b, &c.a, imprinted.pieces_of_second()),
        ];
        for (body, own, other, pieces) in &sides {
            for original in m.faces_of_body(*body) {
                assert!(
                    pieces
                        .iter()
                        .any(|&piece| m.face(piece).sense == m.face(original).sense
                            && matches!(
                                m.geom().same_surface(
                                    m.face(piece).surface,
                                    m.face(original).surface,
                                    &mut dec
                                ),
                                Outcome::ProvenCoincident(_)
                            )),
                    "an operand face disappeared during imprint"
                );
            }
            for &piece in pieces {
                let mut point = exact(face_point(&mut m, piece, &mut dec));
                assert_eq!(
                    implicit_normal(m.geom().surface(m.face(piece).surface), &point).0,
                    r(0)
                );
                assert_eq!(own.locate(&point), Location::On);
                let is_isolated = graph.points().iter().any(|&id| {
                    let real = m.geom().position(id);
                    [&real.x, &real.y, &real.z]
                        .into_iter()
                        .zip([&point.x, &point.y, &point.z])
                        .all(|(x, y)| {
                            matches!(
                                Decider::default().compare(x, &Real::from(y.clone())),
                                Outcome::ProvenCoincident(_)
                            )
                        })
                });
                if is_isolated {
                    assert!(
                        matches!(c.number, 19 | 20 | 21 | 29 | 34 | 37 | 38),
                        "unexpected isolated contact witness"
                    );
                    point = contact_witness(&m, piece, own, other, &point);
                }
                let want = match other.locate(&point) {
                    Location::Inside => FaceClass::Inside,
                    Location::Outside => FaceClass::Outside,
                    Location::On => {
                        let n = own.outward(&point);
                        let other_n = other.outward(&point);
                        assert!(
                            n.cross(&other_n).is_zero(),
                            "face_point landed on a transverse section edge"
                        );
                        if n.dot(&other_n) > r(0) {
                            FaceClass::OnSame
                        } else {
                            FaceClass::OnOpposite
                        }
                    }
                };
                assert_eq!(
                    classification.class_of(piece),
                    Some(want),
                    "{id}: piece {piece:?} at {point:?}"
                );
            }
        }
        for edge in imprinted.common_edges() {
            for (_, _, _, pieces) in &sides {
                assert!(
                    pieces
                        .iter()
                        .any(|&f| m.face(f).loops.iter().any(|&l| match m.loop_(l) {
                            Loop::Coedges(cs) => cs.iter().any(|&c| m.coedge(c).edge == edge),
                            _ => false,
                        })),
                    "common edge missing from one side"
                );
            }
        }
        for op in [BoolOp::Unite, BoolOp::Intersect, BoolOp::Subtract] {
            let mut want = Vec::new();
            for (side, pieces) in [imprinted.pieces_of_first(), imprinted.pieces_of_second()]
                .into_iter()
                .enumerate()
            {
                for f in pieces {
                    let class = classification.class_of(f).unwrap();
                    let keep_it = matches!(
                        (op, side, class),
                        (BoolOp::Unite, 0, FaceClass::Outside | FaceClass::OnSame)
                            | (BoolOp::Unite, 1, FaceClass::Outside)
                            | (BoolOp::Intersect, 0, FaceClass::Inside | FaceClass::OnSame)
                            | (BoolOp::Intersect, 1, FaceClass::Inside)
                            | (
                                BoolOp::Subtract,
                                0,
                                FaceClass::Outside | FaceClass::OnOpposite
                            )
                            | (BoolOp::Subtract, 1, FaceClass::Inside)
                    );
                    if keep_it {
                        want.push((
                            f,
                            if side == 1 && op == BoolOp::Subtract {
                                Sense::Reversed
                            } else {
                                Sense::Same
                            },
                        ));
                    }
                }
            }
            let stats = dec.stats().clone();
            let selection = keep(&classification, &imprinted, op);
            assert_eq!(selection.faces(), want, "{id} {op:?}");
            assert_eq!(&stats, dec.stats());
            let result = sew(&mut m, &imprinted, &selection, &mut dec);
            let expected = expectations(id)
                .into_iter()
                .find(|(o, rev, _)| *o == op && !*rev)
                .map(|(_, _, e)| e);
            if let Some(Expected::Refused(code, Stage::Sew)) = &expected {
                assert_eq!(result.unwrap_err().code, *code);
                continue;
            }
            let body = exact(result);
            topology_rules(&m, body);
            if let Some(Expected::Counts(want)) = expected {
                check_counts(&m, body, &want);
            }
            check_faces(&mut m, body, a, b, &c.a, &c.b, op);
            for (point, location) in witnesses(id, op) {
                query(&mut m, body, &point, location, &mut dec);
            }
        }
        assert_eq!(before, topology_snapshot(&m, &[a, b, unrelated]));
    }
});

blind_test!(
    case_53_reports_preserve_fresh_and_seeded_caller_statistics,
    {
        for &id in IDS {
            for seeded in [false, true] {
                for (op, reverse, expected) in expectations(id) {
                    let c = case(id);
                    let mut m = Model::new();
                    let a = c.a.build(&mut m);
                    let b = c.b.build(&mut m);
                    let (a, b) = if reverse { (b, a) } else { (a, b) };
                    let mut dec = Decider::default();
                    if seeded {
                        // Caller work must survive every stage, including keep.
                        // sqrt(2) at width 10^-100 seeds a nonzero precision peak;
                        // its exact value follows from x²=2, not the boolean.
                        dec.enclose(&Real::from_i64(2).sqrt(), &q("1e-100"))
                            .unwrap();
                        assert!(matches!(dec.sign(&Real::from_i64(7)), Outcome::Proven(_)));
                        assert!(dec.stats().decisions > 0);
                        assert!(dec.stats().max_bits > 0);
                    }
                    let before = dec.stats().clone();
                    let result = boolean(&mut m, a, b, op, &mut dec);
                    match expected {
                        Expected::Counts(want) => {
                            let answer = result.unwrap();
                            check_report(&answer, &before, dec.stats());
                            check_counts(&m, answer.body, &want);
                        }
                        Expected::Refused(code, stage) => {
                            let refusal = result.unwrap_err();
                            assert_eq!(refusal.refusal.code, code);
                            assert_eq!(refusal.stage, stage);
                            // Refusals have no report, but may not clear totals.
                            stats_delta(&before, dec.stats());
                        }
                    }
                }
            }
        }
    }
);
