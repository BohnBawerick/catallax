//! Atomic propagation of a reachable piece 4 refusal.
mod common;
use common::*;

fn cone_with_intersection_ring() -> (Model, BodyId, CurveId) {
    let mut m = Model::new();
    let original = cone(&mut m, n(3), n(0), n(4), &Pose::identity()).unwrap();
    let faces = m.faces_of_body(original);
    let supports = [m.face(faces[0]).surface, m.face(faces[1]).surface];
    // (3,0,0) lies on the base plane and on the cone whose apex is (0,0,4)
    // and slope is 3/4. This same circle may be stored by its meeting definition.
    let witness = m
        .geom_mut()
        .add_point(PointDef::Explicit(RPoint3::new(n(3), n(0), n(0))));
    let curve = m
        .geom_mut()
        .add_curve(Curve::Intersection(IntersectionCurve {
            surfaces: supports,
            witness,
        }));
    let edge = m.add_edge(Edge {
        curve,
        // grad(bottom plane) x grad(cone) points opposite the original ring.
        sense: Sense::Reversed,
        start: None,
        end: None,
    });
    let mut new_faces = Vec::new();
    for fid in faces {
        let f = m.face(fid).clone();
        let mut loops = Vec::new();
        for lid in f.loops {
            match m.loop_(lid).clone() {
                Loop::Vertex(v) => loops.push(m.add_loop(Loop::Vertex(v))),
                Loop::Coedges(cs) => {
                    let cs = cs
                        .iter()
                        .map(|&c| {
                            let sense = m.coedge(c).sense;
                            m.add_coedge(Coedge { edge, sense })
                        })
                        .collect();
                    loops.push(m.add_loop(Loop::Coedges(cs)));
                }
            }
        }
        new_faces.push(m.add_face(Face {
            surface: f.surface,
            sense: f.sense,
            loops,
        }));
    }
    let shell = m.add_shell(Shell { faces: new_faces });
    let lump = m.add_lump(Lump {
        shells: vec![shell],
    });
    let b = body(&mut m, vec![lump]);
    (m, b, curve)
}
#[test]
fn unsupported_curve_fixture_is_topologically_valid() {
    let (m, b, _) = cone_with_intersection_ring();
    assert_eq!(m.check(b), vec![]);
}
blind_test!(
    upstream_curve_refusal_preserves_detail_and_writes_nothing,
    {
        let (m, b, c) = cone_with_intersection_ring();
        let expected =
            catallax_intersect::curve_form(m.geom(), c, &mut Decider::default()).unwrap_err();
        assert_eq!(expected.code, "intersect.unsupported-curve");
        let mut out = b"unchanged prefix".to_vec();
        let before = out.clone();
        let actual =
            write_step(&m, &[b], &options(17), &mut Decider::default(), &mut out).unwrap_err();
        assert_eq!(actual, expected);
        assert_eq!(out, before);
    }
);
