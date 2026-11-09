//! STEP read-back witnesses with a 3.5e-10 mm straight edge.
use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, Rational, Real};
use catallax_shapes::{GeomKey, Point3, PointDef, Pose};
use catallax_step::{StepOptions, write_step};
use catallax_topology::{self as topology, Coedge, Edge, Loop, Model, Sense, Vertex};

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn main() {
    let directory = std::path::PathBuf::from(std::env::args().nth(1).expect("output directory"));
    std::fs::create_dir_all(&directory).unwrap();
    for name in ["split", "readable-split", "thin"] {
        let mut model = Model::new();
        let mut body = topology::block(
            &mut model,
            [
                if name == "thin" {
                    q("3.5e-10")
                } else {
                    q("10")
                },
                q("10"),
                q("10"),
            ],
            &Pose::identity(),
        )
        .unwrap();
        if name.ends_with("split") {
            let old = model.edges_of_body(body)[0];
            let edge = model.edge(old).clone();
            let a = model
                .geom()
                .position(model.vertex(edge.start.unwrap()).point);
            let b = model.geom().position(model.vertex(edge.end.unwrap()).point);
            let t = Real::from(q(if name == "split" { "3.5e-11" } else { "2e-8" }));
            let point = Point3::new(
                &a.x + &(&(&b.x - &a.x) * &t),
                &a.y + &(&(&b.y - &a.y) * &t),
                &a.z + &(&(&b.z - &a.z) * &t),
            );
            let point = model.geom_mut().add_point(PointDef::Explicit(point));
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Point(point),
                b: GeomKey::Curve(edge.curve),
                at: None,
                origin: FactOrigin::Construction("sub_tolerance.split"),
            });
            let vertex = model.add_vertex(Vertex { point });
            let first = model.add_edge(Edge {
                end: Some(vertex),
                ..edge.clone()
            });
            let second = model.add_edge(Edge {
                start: Some(vertex),
                ..edge
            });
            let mut value = model.body(body).clone();
            for lump in &mut value.lumps {
                let mut l = model.lump(*lump).clone();
                for shell in &mut l.shells {
                    let mut s = model.shell(*shell).clone();
                    for face in &mut s.faces {
                        let mut f = model.face(*face).clone();
                        for wire in &mut f.loops {
                            let Loop::Coedges(cs) = model.loop_(*wire) else {
                                unreachable!()
                            };
                            let mut replacements = Vec::new();
                            for c in cs.clone() {
                                let coedge = model.coedge(c).clone();
                                if coedge.edge == old {
                                    let edges = if coedge.sense == Sense::Same {
                                        [first, second]
                                    } else {
                                        [second, first]
                                    };
                                    for edge in edges {
                                        replacements.push(model.add_coedge(Coedge {
                                            edge,
                                            sense: coedge.sense,
                                        }));
                                    }
                                } else {
                                    replacements.push(model.add_coedge(coedge));
                                }
                            }
                            *wire = model.add_loop(Loop::Coedges(replacements));
                        }
                        *face = model.add_face(f);
                    }
                    *shell = model.add_shell(s);
                }
                *lump = model.add_lump(l);
            }
            body = model.add_body(value);
        }
        assert!(model.check(body).is_empty());
        let original = model.counts(body);
        let mut bytes = Vec::new();
        match write_step(
            &model,
            &[body],
            &StepOptions::default(),
            &mut Decider::default(),
            &mut bytes,
        ) {
            Ok(answer) => {
                if name == "split" {
                    assert_eq!(
                        (
                            answer.value.faces,
                            answer.value.edges,
                            answer.value.vertices
                        ),
                        (6, 12, 8)
                    );
                    assert_eq!(answer.value.max_deviation, Rational::zero());
                    let [repair] = answer.value.regularizations.as_slice() else {
                        panic!("one subdivision repair")
                    };
                    assert_eq!((repair.removed_faces, repair.removed_edges), (0, 1));
                    assert_eq!(repair.max_deviation, Rational::zero());
                } else if name == "readable-split" {
                    assert!(answer.value.regularizations.is_empty());
                }
                std::fs::write(directory.join(format!("{name}.step")), bytes).unwrap();
                std::fs::write(
                    directory.join(format!("{name}.report")),
                    format!("{:?}", answer.value),
                )
                .unwrap();
            }
            Err(refusal) => {
                assert!(bytes.is_empty());
                let mut existing = vec![0x7f];
                assert_eq!(
                    write_step(
                        &model,
                        &[body],
                        &StepOptions::default(),
                        &mut Decider::default(),
                        &mut existing
                    )
                    .unwrap_err()
                    .code,
                    refusal.code
                );
                assert_eq!(existing, [0x7f]);
                std::fs::write(
                    directory.join(format!("{name}.refusal")),
                    format!("{}: {}", refusal.code, refusal.detail),
                )
                .unwrap();
            }
        }
        assert_eq!(model.counts(body), original);
        assert!(model.check(body).is_empty());
    }
}
