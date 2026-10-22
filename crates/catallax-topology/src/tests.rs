use super::*;
use crate::{BodyKind, block, cone, cylinder, empty, sphere, torus};
use catallax_numbers::{Fact, FactKind, FactOrigin, MarkLog, Rational};
use catallax_shapes::{
    Circle, Curve, GeomKey, Line, Plane, PointDef, Pose, RPoint3, RVector3, Rotation, Surface,
};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(r(x), r(y), r(z))
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}

fn body(m: &mut Model, faces: Vec<FaceId>) -> BodyId {
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

fn fact(m: &mut Model, a: GeomKey, b: GeomKey) {
    m.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("test"),
    });
}

fn solid_builds(m: &mut Model, scale: &Rational, pose: &Pose) -> Vec<BodyId> {
    vec![
        block(m, [scale * &r(2), scale * &r(3), scale * &r(5)], pose).unwrap(),
        cylinder(m, scale * &r(2), scale * &r(5), pose).unwrap(),
        sphere(m, scale * &r(3), pose).unwrap(),
        cone(m, scale * &r(2), scale.clone(), scale * &r(2), pose).unwrap(),
        cone(m, scale.clone(), scale * &r(3), scale * &r(2), pose).unwrap(),
        cone(m, scale * &r(3), r(0), scale * &r(4), pose).unwrap(),
        torus(m, scale * &r(3), scale.clone(), pose).unwrap(),
        empty(m),
    ]
}

fn sizes(m: &Model) -> [usize; 12] {
    [
        m.vertices.len(),
        m.edges.len(),
        m.coedges.len(),
        m.loops.len(),
        m.faces.len(),
        m.shells.len(),
        m.lumps.len(),
        m.bodies.len(),
        m.geom.surface_count(),
        m.geom.curve_count(),
        m.geom.point_count(),
        m.geom.facts().len(),
    ]
}

#[test]
fn primitive_counts_growth_and_construction_facts() {
    let mut m = Model::new();
    let ids = solid_builds(&mut m, &r(1), &Pose::identity());
    // Cases 1 through 7. Counts are independently specified by the contract.
    let expected = [
        (
            6,
            6,
            12,
            0,
            8,
            vec![(SurfaceKind::Plane, 6)],
            vec![(CurveKind::Line, 12)],
        ),
        (
            3,
            4,
            2,
            2,
            0,
            vec![(SurfaceKind::Plane, 2), (SurfaceKind::Cylinder, 1)],
            vec![(CurveKind::Circle, 2)],
        ),
        (1, 0, 0, 0, 0, vec![(SurfaceKind::Sphere, 1)], vec![]),
        (
            3,
            4,
            2,
            2,
            0,
            vec![(SurfaceKind::Plane, 2), (SurfaceKind::Cone, 1)],
            vec![(CurveKind::Circle, 2)],
        ),
        (
            3,
            4,
            2,
            2,
            0,
            vec![(SurfaceKind::Plane, 2), (SurfaceKind::Cone, 1)],
            vec![(CurveKind::Circle, 2)],
        ),
        (
            2,
            3,
            1,
            1,
            1,
            vec![(SurfaceKind::Plane, 1), (SurfaceKind::Cone, 1)],
            vec![(CurveKind::Circle, 1)],
        ),
        (1, 0, 0, 0, 0, vec![(SurfaceKind::Torus, 1)], vec![]),
    ];
    for (id, (faces, loops, edges, rings, vertices, fk, ek)) in ids.iter().zip(expected) {
        assert_eq!(
            m.counts(*id),
            EntityCounts {
                lumps: 1,
                shells: 1,
                faces,
                loops,
                edges,
                ring_edges: rings,
                vertices,
                faces_by_kind: fk.into_iter().collect(),
                edges_by_kind: ek.into_iter().collect()
            }
        );
        assert!(m.body(*id).marks.is_empty());
        assert_eq!(m.check(*id), vec![]);
        for f in m.faces_of_body(*id) {
            let face = m.face(f);
            for l in &face.loops {
                match m.loop_(*l) {
                    Loop::Vertex(vertex) => assert!(m.geom.facts().holds(
                        FactKind::Incident,
                        GeomKey::Point(m.vertex(*vertex).point),
                        GeomKey::Surface(face.surface)
                    )),
                    Loop::Coedges(cs) => {
                        for c in cs {
                            let edge = m.edge(m.coedge(*c).edge);
                            assert!(m.geom.facts().holds(
                                FactKind::Incident,
                                GeomKey::Curve(edge.curve),
                                GeomKey::Surface(face.surface)
                            ));
                            for vertex in [edge.start, edge.end].into_iter().flatten() {
                                assert!(m.geom.facts().holds(
                                    FactKind::Incident,
                                    GeomKey::Point(m.vertex(vertex).point),
                                    GeomKey::Curve(edge.curve)
                                ));
                            }
                        }
                    }
                }
            }
        }
    }
    assert_eq!(m.counts(ids[7]), EntityCounts::default());
    assert_eq!(m.check(ids[7]), vec![]);
    assert_eq!(&sizes(&m)[..11], &[9, 19, 38, 21, 19, 7, 7, 8, 19, 19, 9]);
    let block_faces = m.faces_of_body(ids[0]);
    for (i, j) in [(0, 1), (2, 3), (4, 5)] {
        assert!(m.geom.facts().holds(
            FactKind::Parallel,
            GeomKey::Surface(m.face(block_faces[i]).surface),
            GeomKey::Surface(m.face(block_faces[j]).surface)
        ));
    }
    assert!(!m.geom.facts().holds(
        FactKind::Parallel,
        GeomKey::Surface(m.face(block_faces[0]).surface),
        GeomKey::Surface(m.face(block_faces[2]).surface)
    ));
    let zero = m.vertices.handles().into_iter().find(|id| matches!(m.geom.point(m.vertex(*id).point), PointDef::Explicit(pt) if pt.as_rational().unwrap() == p(0,0,0))).unwrap();
    assert!(!m.geom.facts().holds(
        FactKind::Incident,
        GeomKey::Point(m.vertex(zero).point),
        GeomKey::Surface(m.face(block_faces[1]).surface)
    ));
    for (id, z, axis, slope) in [
        (
            ids[3],
            r(4),
            v(0, 0, -1),
            Rational::from_ratio(1, 2).unwrap(),
        ),
        (ids[4], r(-1), v(0, 0, 1), r(1)),
        (
            ids[5],
            r(4),
            v(0, 0, -1),
            Rational::from_ratio(3, 4).unwrap(),
        ),
    ] {
        let faces = m.faces_of_body(id);
        let Surface::Cone(c) = m.geom.surface(m.face(*faces.last().unwrap()).surface) else {
            panic!("cone side")
        };
        assert_eq!(c.apex.as_rational().unwrap(), RPoint3::new(r(0), r(0), z));
        assert_eq!(c.axis.as_rational().unwrap(), axis);
        assert_eq!(c.slope.rational_ref(), Some(&slope));
    }
}

#[test]
fn block_geometry_orientation_walk_and_backlinks() {
    let mut m = Model::new();
    let b = block(&mut m, [r(2), r(3), r(5)], &Pose::identity()).unwrap();
    let expected = [
        (
            p(0, 0, 0),
            v(-1, 0, 0),
            v(0, 1, 0),
            [p(0, 0, 0), p(0, 0, 5), p(0, 3, 5), p(0, 3, 0)],
        ),
        (
            p(2, 0, 0),
            v(1, 0, 0),
            v(0, 1, 0),
            [p(2, 0, 0), p(2, 3, 0), p(2, 3, 5), p(2, 0, 5)],
        ),
        (
            p(0, 0, 0),
            v(0, -1, 0),
            v(0, 0, 1),
            [p(0, 0, 0), p(2, 0, 0), p(2, 0, 5), p(0, 0, 5)],
        ),
        (
            p(0, 3, 0),
            v(0, 1, 0),
            v(0, 0, 1),
            [p(0, 3, 0), p(0, 3, 5), p(2, 3, 5), p(2, 3, 0)],
        ),
        (
            p(0, 0, 0),
            v(0, 0, -1),
            v(1, 0, 0),
            [p(0, 0, 0), p(0, 3, 0), p(2, 3, 0), p(2, 0, 0)],
        ),
        (
            p(0, 0, 5),
            v(0, 0, 1),
            v(1, 0, 0),
            [p(0, 0, 5), p(2, 0, 5), p(2, 3, 5), p(0, 3, 5)],
        ),
    ];
    let mut walk_edges = Vec::new();
    let mut walk_vertices = Vec::new();
    for (f, (origin, normal, ref_dir, corners)) in m.faces_of_body(b).into_iter().zip(expected) {
        let face = m.face(f);
        assert!(
            m.geom
                .surface(face.surface)
                .same_definition(&Surface::Plane(Plane {
                    origin: origin.into(),
                    normal: normal.clone().into(),
                    ref_dir: ref_dir.into()
                }))
        );
        let l = face.loops[0];
        let Loop::Coedges(cs) = m.loop_(l) else {
            panic!("block loop")
        };
        let points: Vec<_> = cs
            .iter()
            .map(|c| {
                let (start, end) = m.coedge_ends(*c);
                let start = start.unwrap();
                let end = end.unwrap();
                let PointDef::Explicit(a) = m.geom.point(m.vertex(start).point) else {
                    panic!("explicit")
                };
                let a = a.as_rational().expect("rational block fixture");
                let PointDef::Explicit(z) = m.geom.point(m.vertex(end).point) else {
                    panic!("explicit")
                };
                let z = z.as_rational().expect("rational block fixture");
                assert!(normal.cross(&z.minus(&a)).dot(&p(1, 1, 2).minus(&a)) > r(0));
                assert_eq!(m.owner_of_coedge(*c), Some((l, f)));
                let edge = m.coedge(*c).edge;
                if !walk_edges.contains(&edge) {
                    walk_edges.push(edge);
                }
                for vertex in [start, end] {
                    if !walk_vertices.contains(&vertex) {
                        walk_vertices.push(vertex);
                    }
                }
                a
            })
            .collect();
        let offset = corners.iter().position(|pt| *pt == points[0]).unwrap();
        for i in 0..4 {
            assert_eq!(points[i], corners[(i + offset) % 4]);
        }
    }
    assert_eq!(m.edges_of_body(b), walk_edges);
    assert_eq!(m.vertices_of_body(b), walk_vertices);
    for edge in m.edges_of_body(b) {
        let cs = m.coedges_of_edge(edge);
        assert_eq!(cs.len(), 2);
        assert!(cs[0].index() < cs[1].index());
        assert_ne!(m.coedge(cs[0]).sense, m.coedge(cs[1]).sense);
    }
    let f = m.faces_of_body(b)[0];
    let reused = m.add_face(m.face(f).clone());
    let Loop::Coedges(cs) = m.loop_(m.face(f).loops[0]) else {
        panic!()
    };
    assert_eq!(m.owner_of_coedge(cs[0]).unwrap().1, f);
    assert_ne!(reused, f);
    let e = m.edges_of_body(b)[0];
    let c = m.add_coedge(Coedge {
        edge: e,
        sense: Sense::Same,
    });
    assert_eq!(m.owner_of_coedge(c), None);
    assert_eq!(m.coedges_of_edge(e).last(), Some(&c));
    // Back-links from other bodies or unfinished construction do not alter checks.
    assert_eq!(m.check(b), vec![]);
}

#[test]
fn refusals_are_ordered_and_transactional() {
    let mut m = Model::new();
    solid_builds(&mut m, &r(1), &Pose::identity());
    let before = sizes(&m);
    let pose = Pose::identity();
    for axis in 0..3 {
        for bad in [0, -1] {
            let mut size = [r(1), r(1), r(1)];
            size[axis] = r(bad);
            assert_eq!(
                block(&mut m, size, &pose).unwrap_err().code,
                "topology.non-positive-size"
            );
            assert_eq!(sizes(&m), before);
        }
    }
    let failures = [
        cylinder(&mut m, r(0), r(1), &pose),
        cylinder(&mut m, r(1), r(-1), &pose),
        sphere(&mut m, r(0), &pose),
        cone(&mut m, r(0), r(1), r(1), &pose),
        cone(&mut m, r(1), r(0), r(0), &pose),
        cone(&mut m, r(1), r(-1), r(1), &pose),
        cone(&mut m, r(-1), r(-1), r(1), &pose),
        torus(&mut m, r(1), r(0), &pose),
    ];
    for failure in failures {
        assert_eq!(failure.unwrap_err().code, "topology.non-positive-size");
    }
    assert_eq!(
        cone(&mut m, r(2), r(2), r(1), &pose).unwrap_err().code,
        "topology.cone-is-cylinder"
    );
    for minor in [1, 2] {
        assert_eq!(
            torus(&mut m, r(1), r(minor), &pose).unwrap_err().code,
            "topology.torus-self-intersecting"
        );
    }
    assert_eq!(sizes(&m), before);
}

#[test]
fn scaling_poses_and_thread_determinism() {
    let mut local = Model::new();
    let ids = solid_builds(&mut local, &r(1), &Pose::identity());
    for exponent in [-30, 0, 30] {
        let scale = Rational::parse(&format!("1e{exponent}")).unwrap();
        let mut scaled = Model::new();
        let bs = solid_builds(&mut scaled, &scale, &Pose::identity());
        for (b, original) in bs.into_iter().zip(&ids) {
            assert_eq!(scaled.counts(b), local.counts(*original));
            assert_eq!(scaled.check(b), vec![]);
        }
    }
    for seed in 0..8 {
        let rotation = Rotation::from_quaternion(r(1), r(seed), r(seed + 1), r(2)).unwrap();
        let pose = Pose::new(
            rotation,
            (v(10 + seed, -13, 123)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap();
        let mut moved = Model::new();
        let bs = solid_builds(&mut moved, &r(1), &pose);
        for (b, original) in bs.into_iter().zip(&ids) {
            assert_eq!(moved.counts(b), local.counts(*original));
            assert_eq!(moved.check(b), vec![]);
        }
        for f in local.faces.handles() {
            assert_eq!(moved.face(f), local.face(f));
            assert!(
                moved
                    .geom
                    .surface(moved.face(f).surface)
                    .same_definition(&local.geom.surface(local.face(f).surface).transformed(&pose))
            );
        }
        for e in local.edges.handles() {
            assert_eq!(moved.edge(e), local.edge(e));
            assert!(
                moved.geom.curve(moved.edge(e).curve).same_definition(
                    &local
                        .geom
                        .curve(local.edge(e).curve)
                        .transformed(&pose)
                        .unwrap()
                )
            );
        }
        for vertex in local.vertices.handles() {
            let PointDef::Explicit(pt) = local.geom.point(local.vertex(vertex).point) else {
                panic!()
            };
            let PointDef::Explicit(actual) = moved.geom.point(moved.vertex(vertex).point) else {
                panic!()
            };
            assert_eq!(
                actual.as_rational().unwrap(),
                pose.apply_point(pt).as_rational().unwrap()
            );
        }
        assert_eq!(
            moved.geom.facts().iter().collect::<Vec<_>>(),
            local.geom.facts().iter().collect::<Vec<_>>()
        );
    }
    let handles: Vec<_> = (0..8)
        .map(|_| {
            std::thread::spawn(|| {
                let mut m = Model::new();
                let ids = solid_builds(&mut m, &r(1), &Pose::identity());
                (ids, m)
            })
        })
        .collect();
    for handle in handles {
        let (bs, m) = handle.join().unwrap();
        assert_eq!(bs, ids);
        for f in local.faces.handles() {
            assert_eq!(m.face(f), local.face(f));
        }
        for b in bs {
            assert_eq!(m.check(b), vec![]);
        }
        assert_eq!(
            m.geom.facts().iter().collect::<Vec<_>>(),
            local.geom.facts().iter().collect::<Vec<_>>()
        );
    }
}

fn ring_face(
    m: &mut Model,
    surface: catallax_shapes::SurfaceId,
    edge: EdgeId,
    senses: &[Sense],
) -> FaceId {
    let loops = senses
        .iter()
        .map(|sense| {
            let c = m.add_coedge(Coedge {
                edge,
                sense: *sense,
            });
            m.add_loop(Loop::Coedges(vec![c]))
        })
        .collect();
    m.add_face(Face {
        surface,
        sense: Sense::Same,
        loops,
    })
}

#[test]
fn sound_ring_faces_and_same_face_double_use() {
    let mut m = Model::new();
    let b = sphere(&mut m, r(2), &Pose::identity()).unwrap();
    let surface = m.face(m.faces_of_body(b)[0]).surface;
    let curve = m.geom_mut().add_curve(Curve::Circle(Circle {
        center: p(0, 0, 0).into(),
        normal: v(0, 0, 1).into(),
        radius_sq: r(4).into(),
        ref_dir: v(1, 0, 0).into(),
    }));
    fact(&mut m, GeomKey::Curve(curve), GeomKey::Surface(surface));
    let edge = m.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let north = ring_face(&mut m, surface, edge, &[Sense::Same]);
    let south = ring_face(&mut m, surface, edge, &[Sense::Reversed]);
    let split = body(&mut m, vec![north, south]);
    assert_eq!(m.check(split), vec![]);
    let b = torus(&mut m, r(3), r(1), &Pose::identity()).unwrap();
    let surface = m.face(m.faces_of_body(b)[0]).surface;
    let curve = m.geom_mut().add_curve(Curve::Circle(Circle {
        center: p(3, 0, 0).into(),
        normal: v(0, 1, 0).into(),
        radius_sq: r(1).into(),
        ref_dir: v(1, 0, 0).into(),
    }));
    fact(&mut m, GeomKey::Curve(curve), GeomKey::Surface(surface));
    let edge = m.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let face = ring_face(&mut m, surface, edge, &[Sense::Same, Sense::Reversed]);
    let cut = body(&mut m, vec![face]);
    assert_eq!(m.check(cut), vec![]);
}

fn assert_fault(m: &Model, b: BodyId, code: &str) {
    let faults = m.check(b);
    assert!(
        faults.iter().any(|f| f.code == code),
        "missing {code} in {faults:?}"
    );
    assert_eq!(m.check(b), faults);
}

#[test]
fn boundary_cycle_and_edge_use_faults() {
    let mut m = Model::new();
    let b = block(&mut m, [r(1), r(1), r(1)], &Pose::identity()).unwrap();
    let faces = m.faces_of_body(b);
    let missing_face = body(&mut m, faces[1..].to_vec());
    assert_fault(&m, missing_face, "topology.edge-use");
    assert_fault(&m, missing_face, "topology.euler");
    let original = m.face(faces[0]).clone();
    let Loop::Coedges(cs) = m.loop_(original.loops[0]).clone() else {
        panic!()
    };
    for flip in [false, true] {
        let mut changed = cs.clone();
        if flip {
            let mut c = m.coedge(changed[0]).clone();
            c.sense = if c.sense == Sense::Same {
                Sense::Reversed
            } else {
                Sense::Same
            };
            changed[0] = m.add_coedge(c);
        } else {
            changed.pop();
        }
        let l = m.add_loop(Loop::Coedges(changed));
        let f = m.add_face(Face {
            loops: vec![l],
            ..original.clone()
        });
        let mut new_faces = faces.clone();
        new_faces[0] = f;
        let broken = body(&mut m, new_faces);
        assert_fault(&m, broken, "topology.loop-not-closed");
        if flip {
            assert_fault(&m, broken, "topology.edge-use");
        }
    }
    // A duplicate face or an aliasing loop must not be hidden by query deduplication.
    let alias = m.add_face(original);
    for extra in [faces[0], alias] {
        let mut fs = faces.clone();
        fs.push(extra);
        let duplicate = body(&mut m, fs);
        assert_fault(&m, duplicate, "topology.shared-entity");
    }
    assert_eq!(m.check(b), vec![]);
}

#[test]
fn structural_and_incidence_faults() {
    let mut m = Model::new();
    let b = cylinder(&mut m, r(2), r(5), &Pose::identity()).unwrap();
    let faces = m.faces_of_body(b);
    let surface = m.face(faces[0]).surface;
    let plane = m.add_face(Face {
        surface,
        sense: Sense::Same,
        loops: vec![],
    });
    let unbounded = body(&mut m, vec![plane]);
    assert_fault(&m, unbounded, "topology.unbounded-face");
    let empty_loop = m.add_loop(Loop::Coedges(vec![]));
    let f = m.add_face(Face {
        surface,
        sense: Sense::Same,
        loops: vec![empty_loop],
    });
    let b = body(&mut m, vec![f]);
    assert_fault(&m, b, "topology.empty-part");
    let point = m
        .geom_mut()
        .add_point(PointDef::Explicit((p(0, 0, 0)).into()));
    let vertex = m.add_vertex(Vertex { point });
    let vertex_loop = m.add_loop(Loop::Vertex(vertex));
    let f = m.add_face(Face {
        surface,
        sense: Sense::Same,
        loops: vec![vertex_loop],
    });
    let b = body(&mut m, vec![f]);
    assert_fault(&m, b, "topology.vertex-loop");
    assert_fault(&m, b, "topology.missing-incidence");
    let line = m.geom_mut().add_curve(Curve::Line(Line {
        origin: p(0, 0, 0).into(),
        dir: v(1, 0, 0).into(),
    }));
    for (start, end) in [
        (None, None),
        (Some(vertex), None),
        (Some(vertex), Some(vertex)),
    ] {
        let edge = m.add_edge(Edge {
            curve: line,
            sense: Sense::Same,
            start,
            end,
        });
        let f = ring_face(&mut m, surface, edge, &[Sense::Same]);
        let b = body(&mut m, vec![f]);
        assert_fault(&m, b, "topology.edge-ends");
    }
    let edge = m.edges_of_body(m.bodies.handles()[0])[0];
    let curve = m.geom.curve(m.edge(edge).curve).clone();
    let new_curve = m.geom_mut().add_curve(curve);
    let new_edge = m.add_edge(Edge {
        curve: new_curve,
        ..m.edge(edge).clone()
    });
    let f = ring_face(&mut m, surface, new_edge, &[Sense::Same, Sense::Reversed]);
    let b = body(&mut m, vec![f]);
    assert_fault(&m, b, "topology.missing-incidence");
    let b = sphere(&mut m, r(2), &Pose::identity()).unwrap();
    let f = m.faces_of_body(b)[0];
    let second = m.add_face(m.face(f).clone());
    let disconnected = body(&mut m, vec![f, second]);
    assert_fault(&m, disconnected, "topology.disconnected-shell");
    assert_fault(&m, disconnected, "topology.euler");
    let b = m.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![],
        marks: MarkLog::new(),
    });
    assert_fault(&m, b, "topology.body-kind");
    let empty_lump = m.add_lump(Lump { shells: vec![] });
    let b = m.add_body(Body {
        kind: BodyKind::Empty,
        lumps: vec![empty_lump],
        marks: MarkLog::new(),
    });
    assert_fault(&m, b, "topology.body-kind");
    assert_fault(&m, b, "topology.empty-part");
}

#[test]
fn multiple_lumps_voids_and_cross_shell_sharing() {
    let mut m = Model::new();
    let a = block(&mut m, [r(1), r(1), r(1)], &Pose::identity()).unwrap();
    let b = block(
        &mut m,
        [r(1), r(1), r(1)],
        &Pose::new(
            Rotation::identity(),
            (v(1, 1, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let lumps = vec![m.body(a).lumps[0], m.body(b).lumps[0]];
    let touching = m.add_body(Body {
        kind: BodyKind::Solid,
        lumps: lumps.clone(),
        marks: MarkLog::new(),
    });
    assert_eq!(m.check(touching), vec![]);
    let c = m.counts(touching);
    assert_eq!(
        (c.lumps, c.shells, c.faces, c.edges, c.vertices),
        (2, 2, 12, 24, 16)
    );
    let faces = m.faces_of_body(touching);
    let joined_shell = body(&mut m, faces);
    assert_fault(&m, joined_shell, "topology.disconnected-shell");
    let duplicate = m.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lumps[0], lumps[0]],
        marks: MarkLog::new(),
    });
    assert_fault(&m, duplicate, "topology.shared-entity");
    let shell_a = m.lump(lumps[0]).shells[0];
    let alias_lump = m.add_lump(Lump {
        shells: vec![shell_a],
    });
    let duplicate = m.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lumps[0], alias_lump],
        marks: MarkLog::new(),
    });
    assert_fault(&m, duplicate, "topology.shared-entity");
    let outer = block(&mut m, [r(10), r(10), r(10)], &Pose::identity()).unwrap();
    let inner = block(
        &mut m,
        [r(2), r(2), r(2)],
        &Pose::new(
            Rotation::identity(),
            (v(4, 4, 4)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let mut void_faces = Vec::new();
    for f in m.faces_of_body(inner) {
        let face = m.face(f).clone();
        let mut loops = Vec::new();
        for l in face.loops {
            let Loop::Coedges(cs) = m.loop_(l).clone() else {
                panic!()
            };
            let cs = cs
                .into_iter()
                .rev()
                .map(|c| {
                    let mut c = m.coedge(c).clone();
                    c.sense = if c.sense == Sense::Same {
                        Sense::Reversed
                    } else {
                        Sense::Same
                    };
                    m.add_coedge(c)
                })
                .collect();
            loops.push(m.add_loop(Loop::Coedges(cs)));
        }
        void_faces.push(m.add_face(Face {
            surface: face.surface,
            sense: Sense::Reversed,
            loops,
        }));
    }
    let void = m.add_shell(Shell { faces: void_faces });
    let outer = m.lump(m.body(outer).lumps[0]).shells[0];
    let lump = m.add_lump(Lump {
        shells: vec![outer, void],
    });
    let hollow = m.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    assert_eq!(m.check(hollow), vec![]);
    assert_eq!(m.counts(hollow).shells, 2);
}

#[test]
fn block_with_cylindrical_hole_has_genus_one() {
    let mut m = Model::new();
    let b = block(&mut m, [r(4), r(4), r(4)], &Pose::identity()).unwrap();
    let mut faces = m.faces_of_body(b);
    let cutter = cylinder(
        &mut m,
        r(1),
        r(4),
        &Pose::new(
            Rotation::identity(),
            (v(2, 2, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let cutter_faces = m.faces_of_body(cutter);
    let side = m.face(cutter_faces[2]).clone();
    let mut hole_loops = Vec::new();
    for (i, l) in side.loops.iter().enumerate() {
        let Loop::Coedges(cs) = m.loop_(*l) else {
            panic!()
        };
        let edge = m.coedge(cs[0]).edge;
        let curve = m.edge(edge).curve;
        let mut plane = m.face(faces[4 + i]).clone();
        fact(
            &mut m,
            GeomKey::Curve(curve),
            GeomKey::Surface(plane.surface),
        );
        let sense = if i == 0 { Sense::Same } else { Sense::Reversed };
        let c = m.add_coedge(Coedge { edge, sense });
        plane.loops.push(m.add_loop(Loop::Coedges(vec![c])));
        faces[4 + i] = m.add_face(plane);
        let c = m.add_coedge(Coedge {
            edge,
            sense: if i == 0 { Sense::Reversed } else { Sense::Same },
        });
        hole_loops.push(m.add_loop(Loop::Coedges(vec![c])));
    }
    faces.push(m.add_face(Face {
        surface: side.surface,
        sense: Sense::Reversed,
        loops: hole_loops,
    }));
    let holed = body(&mut m, faces);
    assert_eq!(m.check(holed), vec![]);
    let c = m.counts(holed);
    assert_eq!(
        (c.faces, c.loops, c.edges, c.ring_edges, c.vertices),
        (7, 10, 14, 2, 8)
    );
}

#[test]
fn one_vertex_closed_edges_and_repeated_vertex_boundaries() {
    let mut m = Model::new();
    let sphere = sphere(&mut m, r(2), &Pose::identity()).unwrap();
    let surface = m.face(m.faces_of_body(sphere)[0]).surface;
    let curve = m.geom_mut().add_curve(Curve::Circle(Circle {
        center: p(0, 0, 0).into(),
        normal: v(0, 0, 1).into(),
        radius_sq: r(4).into(),
        ref_dir: v(1, 0, 0).into(),
    }));
    let point = m
        .geom_mut()
        .add_point(PointDef::Explicit((p(2, 0, 0)).into()));
    let vertex = m.add_vertex(Vertex { point });
    fact(&mut m, GeomKey::Curve(curve), GeomKey::Surface(surface));
    fact(&mut m, GeomKey::Point(point), GeomKey::Curve(curve));
    let edge = m.add_edge(Edge {
        curve,
        sense: Sense::Reversed,
        start: Some(vertex),
        end: Some(vertex),
    });
    let north = ring_face(&mut m, surface, edge, &[Sense::Same]);
    let south = ring_face(&mut m, surface, edge, &[Sense::Reversed]);
    let split = body(&mut m, vec![north, south]);
    assert_eq!(m.check(split), vec![]);
    assert_eq!(m.counts(split).vertices, 1);
    assert_eq!(m.counts(split).ring_edges, 0);
    let vertex_loop = m.add_loop(Loop::Vertex(vertex));
    let mut face = m.face(north).clone();
    face.loops.push(vertex_loop);
    let face = m.add_face(face);
    let invalid = body(&mut m, vec![face, south]);
    assert_fault(&m, invalid, "topology.vertex-loop");
    let shared_shell = m.add_shell(Shell {
        faces: m
            .shell(m.lump(m.body(split).lumps[0]).shells[0])
            .faces
            .clone(),
    });
    let lump = m.add_lump(Lump {
        shells: vec![m.lump(m.body(split).lumps[0]).shells[0], shared_shell],
    });
    let invalid = m.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    assert_fault(&m, invalid, "topology.shared-entity");
}
