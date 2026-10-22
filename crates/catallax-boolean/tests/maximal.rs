//! Sound, deliberately nonmaximal operands must be normalised by sew.
mod common;
use common::*;

fn incident(m: &mut Model, a: GeomKey, b: GeomKey) {
    m.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("blind-test.exact-construction"),
    });
}
fn point_vertex(m: &mut Model, p: impl Into<catallax_shapes::Point3>, curve: CurveId) -> VertexId {
    let point = m.geom_mut().add_point(PointDef::Explicit(p.into()));
    incident(m, GeomKey::Point(point), GeomKey::Curve(curve));
    m.add_vertex(Vertex { point })
}
fn loop_from(m: &mut Model, edges: &[(EdgeId, Sense)]) -> LoopId {
    let cs = edges
        .iter()
        .map(|&(edge, sense)| m.add_coedge(Coedge { edge, sense }))
        .collect();
    m.add_loop(Loop::Coedges(cs))
}
fn solid(m: &mut Model, faces: Vec<FaceId>) -> BodyId {
    let shell = m.add_shell(Shell { faces });
    let lump = m.add_lump(Lump {
        shells: vec![shell],
    });
    m.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    })
}

fn split_sphere(m: &mut Model, slit: bool) -> BodyId {
    let raw = sphere(m, r(2), &Pose::identity()).unwrap();
    let surface = m.face(m.faces_of_body(raw)[0]).surface;
    let curve = m.geom_mut().add_curve(Curve::Circle(
        Circle::new(
            (RPoint3::origin()).into(),
            (v(0, 0, 1)).into(),
            (r(4)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    incident(m, GeomKey::Curve(curve), GeomKey::Surface(surface));
    let start = point_vertex(m, p("2", "0", "0"), curve);
    let end = if slit {
        point_vertex(m, p("-2", "0", "0"), curve)
    } else {
        start
    };
    let edge = m.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: Some(start),
        end: Some(end),
    });
    let faces = if slit {
        let loop_ = loop_from(m, &[(edge, Sense::Same), (edge, Sense::Reversed)]);
        vec![m.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![loop_],
        })]
    } else {
        [Sense::Same, Sense::Reversed]
            .into_iter()
            .map(|sense| {
                let loop_ = loop_from(m, &[(edge, sense)]);
                m.add_face(Face {
                    surface,
                    sense: Sense::Same,
                    loops: vec![loop_],
                })
            })
            .collect()
    };
    let body = solid(m, faces);
    assert!(m.check(body).is_empty(), "{:?}", m.check(body));
    body
}

fn rebuilt_with_edges(
    m: &mut Model,
    raw: BodyId,
    replacement: &BTreeMap<EdgeId, Vec<EdgeId>>,
) -> BodyId {
    let mut faces = Vec::new();
    for old in m.faces_of_body(raw) {
        let old = m.face(old).clone();
        let mut loops = Vec::new();
        for l in old.loops {
            let Loop::Coedges(cs) = m.loop_(l).clone() else {
                panic!("fixture has no vertex loops")
            };
            let mut uses = Vec::new();
            for c in cs {
                let c = m.coedge(c);
                let mut edges = replacement
                    .get(&c.edge)
                    .cloned()
                    .unwrap_or_else(|| vec![c.edge]);
                if c.sense == Sense::Reversed {
                    edges.reverse();
                }
                uses.extend(edges.into_iter().map(|e| (e, c.sense)));
            }
            loops.push(loop_from(m, &uses));
        }
        faces.push(m.add_face(Face {
            surface: old.surface,
            sense: old.sense,
            loops,
        }));
    }
    let body = solid(m, faces);
    assert!(m.check(body).is_empty(), "{:?}", m.check(body));
    body
}

blind_test!(
    sew_removes_coplanar_splits_slits_and_their_unused_vertices,
    {
        for slit in [false, true] {
            let mut m = Model::new();
            let a = split_sphere(&mut m, slit);
            let b = ball_at(1, [10, 0, 0]).build(&mut m);
            let before = topology_snapshot(&m, &[a, b]);
            let answer = boolean(&mut m, a, b, BoolOp::Subtract, &mut Decider::default()).unwrap();
            check_counts(
                &m,
                answer.body,
                &Counts {
                    lumps: 1,
                    shells: 1,
                    faces: 1,
                    edges: 0,
                    rings: 0,
                    vertices: 0,
                    surfaces: [0, 0, 1],
                    curves: [0; 4],
                },
            );
            assert_eq!(before, topology_snapshot(&m, &[a, b]));
            for (point, want) in [
                (p("0", "0", "0"), Location::Inside),
                (p("2", "0", "0"), Location::On),
                (p("3", "0", "0"), Location::Outside),
            ] {
                query(&mut m, answer.body, &point, want, &mut Decider::default());
            }
        }
    }
);

blind_test!(
    sew_restores_ring_edges_and_merges_degree_two_line_vertices,
    {
        for cylinder_fixture in [false, true] {
            let mut m = Model::new();
            let raw = if cylinder_fixture {
                cylinder(&mut m, r(2), r(5), &Pose::identity()).unwrap()
            } else {
                block(&mut m, [r(2), r(2), r(2)], &Pose::identity()).unwrap()
            };
            let mut replacements = BTreeMap::new();
            for e in m.edges_of_body(raw) {
                let edge = m.edge(e).clone();
                if cylinder_fixture {
                    let Curve::Circle(circle) = m.geom().curve(edge.curve) else {
                        panic!("rim")
                    };
                    let point = circle.center.plus(&v(2, 0, 0).into());
                    let vertex = point_vertex(&mut m, point, edge.curve);
                    let new = m.add_edge(Edge {
                        start: Some(vertex),
                        end: Some(vertex),
                        ..edge
                    });
                    replacements.insert(e, vec![new]);
                } else {
                    let explicit = |m: &Model, v: VertexId| {
                        let PointDef::Explicit(p) = m.geom().point(m.vertex(v).point) else {
                            panic!("block vertex")
                        };
                        p.clone()
                    };
                    let start = explicit(&m, edge.start.unwrap());
                    let end = explicit(&m, edge.end.unwrap());
                    let midpoint = start.plus(&end.minus(&start).scaled(&Real::from(q("1/2"))));
                    let vertex = point_vertex(&mut m, midpoint, edge.curve);
                    let first = m.add_edge(Edge {
                        end: Some(vertex),
                        ..edge.clone()
                    });
                    let second = m.add_edge(Edge {
                        start: Some(vertex),
                        ..edge
                    });
                    replacements.insert(e, vec![first, second]);
                    break;
                }
            }
            let a = rebuilt_with_edges(&mut m, raw, &replacements);
            let b = ball_at(1, [20; 3]).build(&mut m);
            let before = topology_snapshot(&m, &[a, b]);
            let answer = boolean(&mut m, a, b, BoolOp::Subtract, &mut Decider::default()).unwrap();
            assert_eq!(counts(&m, answer.body), counts(&m, raw));
            topology_rules(&m, answer.body);
            assert_eq!(before, topology_snapshot(&m, &[a, b]));
        }
    }
);

#[test]
fn nonmaximal_sphere_fixtures_are_sound_before_the_boolean_gate() {
    for slit in [false, true] {
        let mut m = Model::new();
        let b = split_sphere(&mut m, slit);
        let c = m.counts(b);
        assert!(m.check(b).is_empty());
        assert_eq!(
            (c.faces, c.edges, c.vertices),
            if slit { (1, 1, 2) } else { (2, 1, 1) }
        );
    }
}
