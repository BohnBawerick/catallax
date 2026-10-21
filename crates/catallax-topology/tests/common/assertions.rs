use super::*;
use std::collections::{BTreeMap, BTreeSet};

pub fn assert_builder(m: &Model, b: BodyId, shape: &Primitive, pose: &Pose) {
    let expected = expected(shape, pose);
    let body = m.body(b);
    assert_eq!(body.kind, BodyKind::Solid);
    assert!(body.marks.is_empty());
    assert_eq!(body.lumps.len(), 1);
    assert_eq!(m.lump(body.lumps[0]).shells.len(), 1);
    let fs = &m.shell(shell_of(m, b)).faces;
    assert_eq!(fs.len(), expected.faces.len());
    let (wf, we, wv) = walk(m, b);
    assert_eq!(m.faces_of_body(b), wf);
    assert_eq!(m.edges_of_body(b), we);
    assert_eq!(m.vertices_of_body(b), wv);
    assert_eq!(wv.len(), expected.points.len());
    assert_eq!(we.len(), expected.edges.len());
    let mut vertex_indices = BTreeMap::new();
    for id in &wv {
        let index = expected
            .points
            .iter()
            .position(|p| p == explicit(m, *id))
            .expect("unexpected vertex position");
        assert!(
            !vertex_indices.values().any(|i| *i == index),
            "duplicate corner"
        );
        vertex_indices.insert(*id, index);
    }
    let mut edge_indices = BTreeMap::new();
    for &id in &we {
        let e = m.edge(id);
        assert_eq!(e.sense, Sense::Same);
        let candidate = (
            m.geom().curve(e.curve).clone(),
            e.start.map(|v| vertex_indices[&v]),
            e.end.map(|v| vertex_indices[&v]),
        );
        let index = expected
            .edges
            .iter()
            .position(|e| e == &candidate)
            .expect("unexpected curve or endpoints");
        assert!(
            !edge_indices.values().any(|i| *i == index),
            "duplicate edge"
        );
        edge_indices.insert(id, index);
    }
    for (&id, (surface, sense, loops)) in fs.iter().zip(&expected.faces) {
        let f = m.face(id);
        assert_eq!(m.geom().surface(f.surface), surface);
        assert_eq!(f.sense, *sense);
        assert_eq!(f.loops.len(), loops.len());
        assert_eq!(
            implicit(m.geom().surface(f.surface), &expected.interior).sign(),
            Sign::Negative
        );
        for (&l, expected_loop) in f.loops.iter().zip(loops) {
            match (m.loop_(l), expected_loop) {
                (Loop::Vertex(v), Boundary::Vertex(i)) => assert_eq!(vertex_indices[v], *i),
                (Loop::Coedges(cs), Boundary::Cycle(want)) => {
                    let actual: Vec<_> = cs
                        .iter()
                        .map(|c| {
                            let c = m.coedge(*c);
                            (edge_indices[&c.edge], c.sense)
                        })
                        .collect();
                    assert_eq!(actual.len(), want.len());
                    assert!(
                        (0..want.len()).any(|offset| (0..want.len())
                            .all(|i| actual[i] == want[(i + offset) % want.len()])),
                        "wrong oriented cycle: {actual:?}, expected {want:?}"
                    );
                    for &c in cs {
                        assert_eq!(m.owner_of_coedge(c), Some((l, id)));
                    }
                }
                _ => panic!("wrong boundary kind"),
            }
        }
        if let Surface::Plane(s) = surface
            && matches!(shape, Primitive::Block(_))
        {
            let points: Vec<_> = cs(m, f.loops[0])
                .into_iter()
                .map(|c| {
                    let c = m.coedge(c);
                    let e = m.edge(c.edge);
                    explicit(
                        m,
                        if c.sense == Sense::Same {
                            e.start.unwrap()
                        } else {
                            e.end.unwrap()
                        },
                    )
                    .clone()
                })
                .collect();
            let sum = points.iter().fold(RVector3::zero(), |sum, p| {
                sum.plus(&p.minus(&RPoint3::origin()))
            });
            let center = RPoint3::origin().plus(&sum.scaled(&q(1, 4)));
            for i in 0..4 {
                let d = points[(i + 1) % 4].minus(&points[i]);
                // For an outward-oriented convex face, n cross d points toward its centre.
                assert_eq!(
                    s.normal.cross(&d).dot(&center.minus(&points[i])).sign(),
                    Sign::Positive
                );
            }
        }
    }
    for &e in &we {
        let all = m.coedges_of_edge(e);
        assert_eq!(all.len(), 2);
        assert!(all.windows(2).all(|w| w[0].index() < w[1].index()));
        assert_eq!(m.coedge(all[0]).sense, flip(m.coedge(all[1]).sense));
        for c in all {
            assert_eq!(m.coedge(c).edge, e);
        }
    }
    let (nf, nl, ne, nr, nv, sk, ck) = match shape {
        Primitive::Block(_) => (
            6,
            6,
            12,
            0,
            8,
            vec![(SurfaceKind::Plane, 6)],
            vec![(CurveKind::Line, 12)],
        ),
        Primitive::Cylinder(..) => (
            3,
            4,
            2,
            2,
            0,
            vec![(SurfaceKind::Plane, 2), (SurfaceKind::Cylinder, 1)],
            vec![(CurveKind::Circle, 2)],
        ),
        Primitive::Sphere(_) => (1, 0, 0, 0, 0, vec![(SurfaceKind::Sphere, 1)], vec![]),
        Primitive::Cone(_, t, _) if t.is_zero() => (
            2,
            3,
            1,
            1,
            1,
            vec![(SurfaceKind::Plane, 1), (SurfaceKind::Cone, 1)],
            vec![(CurveKind::Circle, 1)],
        ),
        Primitive::Cone(..) => (
            3,
            4,
            2,
            2,
            0,
            vec![(SurfaceKind::Plane, 2), (SurfaceKind::Cone, 1)],
            vec![(CurveKind::Circle, 2)],
        ),
        Primitive::Torus(..) => (1, 0, 0, 0, 0, vec![(SurfaceKind::Torus, 1)], vec![]),
    };
    assert_eq!(
        m.counts(b),
        EntityCounts {
            lumps: 1,
            shells: 1,
            faces: nf,
            loops: nl,
            edges: ne,
            ring_edges: nr,
            vertices: nv,
            faces_by_kind: sk.into_iter().collect(),
            edges_by_kind: ck.into_iter().collect()
        }
    );
    assert_eq!(euler(m, shell_of(m, b)), expected.chi);
    assert_facts(m, b, shape);
    sound(m, b);
}

fn construction_fact(m: &Model, kind: FactKind, a: GeomKey, b: GeomKey, origin: &'static str) {
    assert!(
        m.geom().facts().holds(kind, a, b),
        "missing {kind:?} {a:?} {b:?}"
    );
    assert!(
        m.geom().facts().iter().any(|f| f.kind == kind
            && ((f.a == a && f.b == b) || (f.a == b && f.b == a))
            && f.origin == FactOrigin::Construction(origin)),
        "required construction provenance absent"
    );
}
pub fn assert_facts(m: &Model, b: BodyId, shape: &Primitive) {
    let (fs, es, _) = walk(m, b);
    for &f in &fs {
        let f = m.face(f);
        for &l in &f.loops {
            match m.loop_(l) {
                Loop::Vertex(v) => construction_fact(
                    m,
                    FactKind::Incident,
                    GeomKey::Point(m.vertex(*v).point),
                    GeomKey::Surface(f.surface),
                    shape.name(),
                ),
                Loop::Coedges(cs) => {
                    for c in cs {
                        construction_fact(
                            m,
                            FactKind::Incident,
                            GeomKey::Curve(m.edge(m.coedge(*c).edge).curve),
                            GeomKey::Surface(f.surface),
                            shape.name(),
                        );
                    }
                }
            }
        }
    }
    for e in es {
        let e = m.edge(e);
        for v in [e.start, e.end].into_iter().flatten() {
            construction_fact(
                m,
                FactKind::Incident,
                GeomKey::Point(m.vertex(v).point),
                GeomKey::Curve(e.curve),
                shape.name(),
            );
        }
    }
    let pairs = match shape {
        Primitive::Block(_) => vec![(0, 1), (2, 3), (4, 5)],
        Primitive::Cylinder(..) => vec![(0, 1)],
        Primitive::Cone(_, t, _) if !t.is_zero() => vec![(0, 1)],
        _ => vec![],
    };
    for (a, b) in pairs {
        construction_fact(
            m,
            FactKind::Parallel,
            GeomKey::Surface(m.face(fs[a]).surface),
            GeomKey::Surface(m.face(fs[b]).surface),
            shape.name(),
        );
    }
    if matches!(shape, Primitive::Block(_)) {
        assert!(!m.geom().facts().holds(
            FactKind::Parallel,
            GeomKey::Surface(m.face(fs[0]).surface),
            GeomKey::Surface(m.face(fs[2]).surface)
        ));
        // 000 is the sole corner belonging to the negative x, y and z faces.
        let vertices = |f: FaceId| -> BTreeSet<VertexId> {
            m.face(f)
                .loops
                .iter()
                .flat_map(|l| cs(m, *l))
                .flat_map(|c| {
                    let e = m.edge(m.coedge(c).edge);
                    [e.start.unwrap(), e.end.unwrap()]
                })
                .collect()
        };
        let x = vertices(fs[0]);
        let y = vertices(fs[2]);
        let z = vertices(fs[4]);
        let origin = *x.iter().find(|v| y.contains(v) && z.contains(v)).unwrap();
        assert!(!m.geom().facts().holds(
            FactKind::Incident,
            GeomKey::Point(m.vertex(origin).point),
            GeomKey::Surface(m.face(fs[1]).surface)
        ));
    }
}

// Structured equality includes handles, stored entity fields, geometry and fact
// insertion order. PointDef has no PartialEq; builder points must be Explicit.
#[derive(Debug, PartialEq)]
pub struct Snapshot {
    bodies: Vec<(BodyId, Body)>,
    lumps: Vec<(LumpId, Lump)>,
    shells: Vec<(ShellId, Shell)>,
    faces: Vec<(FaceId, Face)>,
    loops: Vec<(LoopId, Loop)>,
    coedges: Vec<(CoedgeId, Coedge)>,
    edges: Vec<(EdgeId, Edge)>,
    vertices: Vec<(VertexId, Vertex)>,
    surfaces: Vec<(SurfaceId, Surface)>,
    curves: Vec<(CurveId, Curve)>,
    points: Vec<(PointId, RPoint3)>,
    facts: Vec<Fact<GeomKey>>,
    lengths: [usize; 12],
}
pub fn snapshot(m: &Model, bs: &[BodyId]) -> Snapshot {
    let (mut lu, mut sh, mut fs, mut ls, mut cs, mut es, mut vs, mut ss, mut curves, mut ps) = (
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
    );
    for &b in bs {
        for &l in &m.body(b).lumps {
            lu.insert(l, m.lump(l).clone());
            for &s in &m.lump(l).shells {
                sh.insert(s, m.shell(s).clone());
            }
        }
        let (faces, edges, vertices) = walk(m, b);
        for f in faces {
            let face = m.face(f);
            fs.insert(f, face.clone());
            ss.insert(face.surface, m.geom().surface(face.surface).clone());
            for &l in &face.loops {
                ls.insert(l, m.loop_(l).clone());
                if let Loop::Coedges(co) = m.loop_(l) {
                    for &c in co {
                        cs.insert(c, m.coedge(c).clone());
                    }
                }
            }
        }
        for e in edges {
            let edge = m.edge(e);
            es.insert(e, edge.clone());
            curves.insert(edge.curve, m.geom().curve(edge.curve).clone());
        }
        for v in vertices {
            vs.insert(v, m.vertex(v).clone());
            ps.insert(m.vertex(v).point, explicit(m, v).clone());
        }
    }
    Snapshot {
        bodies: bs.iter().map(|&b| (b, m.body(b).clone())).collect(),
        lumps: lu.into_iter().collect(),
        shells: sh.into_iter().collect(),
        faces: fs.into_iter().collect(),
        loops: ls.into_iter().collect(),
        coedges: cs.into_iter().collect(),
        edges: es.into_iter().collect(),
        vertices: vs.into_iter().collect(),
        surfaces: ss.into_iter().collect(),
        curves: curves.into_iter().collect(),
        points: ps.into_iter().collect(),
        facts: m.geom().facts().iter().copied().collect(),
        lengths: lengths(m),
    }
}
