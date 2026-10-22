use super::*;
pub fn flip(s: Sense) -> Sense {
    if s == Sense::Same {
        Sense::Reversed
    } else {
        Sense::Same
    }
}
pub fn translated(x: &str, y: &str, z: &str) -> Pose {
    Pose::new(
        Rotation::identity(),
        (RVector3::new(r(x), r(y), r(z))).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap()
}
pub fn xpose(x: &str, y: &str, z: &str) -> Pose {
    Pose::new(
        Rotation::from_quaternion(n(1), n(0), n(1), n(0)).unwrap(),
        (RVector3::new(r(x), r(y), r(z))).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap()
}
#[derive(Clone, Debug)]
pub enum Primitive {
    Block([Rational; 3]),
    Cylinder(Rational, Rational),
    Sphere(Rational),
    Cone(Rational, Rational, Rational),
    Torus(Rational, Rational),
}
impl Primitive {
    pub fn build(&self, m: &mut Model, p: &Pose) -> BodyId {
        match self {
            Self::Block(s) => block(m, s.clone(), p),
            Self::Cylinder(r, h) => cylinder(m, r.clone(), h.clone(), p),
            Self::Sphere(r) => sphere(m, r.clone(), p),
            Self::Cone(a, b, h) => cone(m, a.clone(), b.clone(), h.clone(), p),
            Self::Torus(a, b) => torus(m, a.clone(), b.clone(), p),
        }
        .unwrap()
    }
    pub fn scale(&self, k: &Rational) -> Self {
        match self {
            Self::Block(s) => Self::Block(s.clone().map(|x| &x * k)),
            Self::Cylinder(r, h) => Self::Cylinder(r * k, h * k),
            Self::Sphere(r) => Self::Sphere(r * k),
            Self::Cone(a, b, h) => Self::Cone(a * k, b * k, h * k),
            Self::Torus(a, b) => Self::Torus(a * k, b * k),
        }
    }
}
pub fn primitives() -> Vec<Primitive> {
    vec![
        Primitive::Block([n(2), n(3), n(5)]),
        Primitive::Cylinder(n(2), n(5)),
        Primitive::Sphere(n(5)),
        Primitive::Cone(n(3), n(0), n(4)),
        Primitive::Cone(n(2), n(1), n(2)),
        Primitive::Cone(n(1), n(3), n(2)),
        Primitive::Torus(n(3), n(1)),
    ]
}
pub fn body(m: &mut Model, lumps: Vec<LumpId>) -> BodyId {
    m.add_body(Body {
        kind: BodyKind::Solid,
        lumps,
        marks: MarkLog::new(),
    })
}
pub fn void_block() -> (Model, BodyId) {
    let mut m = Model::new();
    let outer = block(&mut m, [n(10), n(10), n(10)], &Pose::identity()).unwrap();
    let inner = block(&mut m, [n(2), n(2), n(2)], &translated("4", "4", "4")).unwrap();
    let mut faces = Vec::new();
    for fid in m.faces_of_body(inner) {
        let f = m.face(fid).clone();
        let mut loops = Vec::new();
        for lid in f.loops {
            let Loop::Coedges(cs) = m.loop_(lid).clone() else {
                panic!()
            };
            let cs = cs
                .into_iter()
                .rev()
                .map(|c| {
                    let c = m.coedge(c).clone();
                    m.add_coedge(Coedge {
                        edge: c.edge,
                        sense: flip(c.sense),
                    })
                })
                .collect();
            loops.push(m.add_loop(Loop::Coedges(cs)));
        }
        faces.push(m.add_face(Face {
            surface: f.surface,
            sense: flip(f.sense),
            loops,
        }));
    }
    let inner = m.add_shell(Shell { faces });
    let outer = m.lump(m.body(outer).lumps[0]).shells[0];
    let lump = m.add_lump(Lump {
        shells: vec![outer, inner],
    });
    let b = body(&mut m, vec![lump]);
    (m, b)
}
pub fn split_sphere() -> (Model, BodyId) {
    let mut m = Model::new();
    let b = sphere(&mut m, n(2), &Pose::identity()).unwrap();
    let surface = m.face(m.faces_of_body(b)[0]).surface;
    let curve = m.geom_mut().add_curve(Curve::Circle(
        Circle::new(
            (RPoint3::origin()).into(),
            (RVector3::new(n(0), n(0), n(1))).into(),
            (n(4)).into(),
            (RVector3::new(n(1), n(0), n(0))).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    m.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a: GeomKey::Curve(curve),
        b: GeomKey::Surface(surface),
        at: None,
        origin: FactOrigin::Construction("step.tests.split-sphere"),
    });
    let edge = m.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let mut faces = Vec::new();
    for sense in [Sense::Same, Sense::Reversed] {
        let c = m.add_coedge(Coedge { edge, sense });
        let l = m.add_loop(Loop::Coedges(vec![c]));
        faces.push(m.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![l],
        }));
    }
    let shell = m.add_shell(Shell { faces });
    let lump = m.add_lump(Lump {
        shells: vec![shell],
    });
    let b = body(&mut m, vec![lump]);
    (m, b)
}
// Rebuild the second touching block with shared geometry but separate topology.
pub fn touching_blocks() -> (Model, BodyId) {
    let mut m = Model::new();
    let a = block(&mut m, [n(1), n(1), n(1)], &Pose::identity()).unwrap();
    let b = block(&mut m, [n(1), n(1), n(1)], &translated("1", "1", "0")).unwrap();
    let mut vertices = BTreeMap::new();
    for v in m.vertices_of_body(b) {
        let mut p = m.vertex(v).point;
        for av in m.vertices_of_body(a) {
            let ap = m.vertex(av).point;
            if let (PointDef::Explicit(x), PointDef::Explicit(y)) =
                (m.geom().point(p), m.geom().point(ap))
                && x.rational_components().unwrap() == y.rational_components().unwrap()
            {
                p = ap;
                break;
            }
        }
        vertices.insert(v, m.add_vertex(Vertex { point: p }));
    }
    let mut edges = BTreeMap::new();
    for e in m.edges_of_body(b) {
        let old = m.edge(e).clone();
        let mut curve = old.curve;
        for ae in m.edges_of_body(a) {
            let ac = m.edge(ae).curve;
            if m.geom().curve(curve).same_definition(m.geom().curve(ac)) {
                curve = ac;
                break;
            }
        }
        let start = old.start.map(|v| vertices[&v]);
        let end = old.end.map(|v| vertices[&v]);
        for v in [start, end].into_iter().flatten() {
            let p = m.vertex(v).point;
            m.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Point(p),
                b: GeomKey::Curve(curve),
                at: None,
                origin: FactOrigin::Construction("step.tests.touch"),
            });
        }
        edges.insert(
            e,
            m.add_edge(Edge {
                curve,
                sense: old.sense,
                start,
                end,
            }),
        );
    }
    let mut faces = Vec::new();
    for f in m.faces_of_body(b) {
        let f = m.face(f).clone();
        let mut loops = Vec::new();
        for l in f.loops {
            let Loop::Coedges(cs) = m.loop_(l).clone() else {
                panic!()
            };
            let cs = cs
                .into_iter()
                .map(|c| {
                    let c = m.coedge(c).clone();
                    let edge = edges[&c.edge];
                    let curve = m.edge(edge).curve;
                    m.geom_mut().record_fact(Fact {
                        kind: FactKind::Incident,
                        a: GeomKey::Curve(curve),
                        b: GeomKey::Surface(f.surface),
                        at: None,
                        origin: FactOrigin::Construction("step.tests.touch"),
                    });
                    m.add_coedge(Coedge {
                        edge,
                        sense: c.sense,
                    })
                })
                .collect();
            loops.push(m.add_loop(Loop::Coedges(cs)));
        }
        faces.push(m.add_face(Face {
            surface: f.surface,
            sense: f.sense,
            loops,
        }));
    }
    let shell = m.add_shell(Shell { faces });
    let lump = m.add_lump(Lump {
        shells: vec![shell],
    });
    let a = m.body(a).lumps[0];
    let b = body(&mut m, vec![a, lump]);
    (m, b)
}
pub fn boolean_fixture(case: u32) -> (Model, BodyId) {
    let (mut m, a, b) = boolean_pair(case);
    let answer = catallax_boolean::boolean(
        &mut m,
        a,
        b,
        catallax_boolean::BoolOp::Intersect,
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty());
    (m, answer.body)
}
// Piece 5 cases 1..40, base variants. Benchmarks also enumerate boundary variants below.
pub fn boolean_pair(case: u32) -> (Model, BodyId, BodyId) {
    let mut m = Model::new();
    let block_p = |s: [&str; 3]| Primitive::Block(s.map(r));
    let cube = |s: &str| block_p([s, s, s]);
    let cyl = |r0: &str, h: &str| Primitive::Cylinder(r(r0), r(h));
    let ball = |r0: &str| Primitive::Sphere(r(r0));
    let origin = Pose::identity();
    let (a, ap, b, bp) = match case {
        1 => (
            cube("2"),
            origin.clone(),
            cube("2"),
            translated("1", "1", "1"),
        ),
        2 => (
            cube("1"),
            origin.clone(),
            cube("1"),
            translated("1", "0", "0"),
        ),
        3 => (
            cube("2"),
            origin.clone(),
            cube("2"),
            translated("2", "1", "0"),
        ),
        4 => (
            cube("1"),
            origin.clone(),
            cube("1"),
            translated("1", "1", "0"),
        ),
        5 => (
            cube("1"),
            origin.clone(),
            cube("1"),
            translated("1", "1", "1"),
        ),
        6 => (
            cube("4"),
            origin.clone(),
            block_p(["2", "2", "6"]),
            translated("1", "1", "-1"),
        ),
        7 => (
            cube("2"),
            origin.clone(),
            cube("2"),
            translated("1", "0", "0"),
        ),
        8 => (cube("2"), origin.clone(), cube("2"), origin.clone()),
        9 => (
            cube("1"),
            origin.clone(),
            cube("1"),
            translated("1.000000000000000000000000000001", "0", "0"),
        ),
        10 => (
            cube("10"),
            origin.clone(),
            cyl("3", "12"),
            translated("5", "5", "-1"),
        ),
        11 => (
            cube("10"),
            origin.clone(),
            cyl("3", "5"),
            translated("5", "5", "6"),
        ),
        12 => (
            cube("10"),
            origin.clone(),
            cyl("2", "20"),
            xpose("-5", "5", "12"),
        ),
        13 => (
            cube("10"),
            origin.clone(),
            cyl("2", "20"),
            xpose("-5", "5", "11.999999999999999999999999999999"),
        ),
        14 => (
            cube("10"),
            origin.clone(),
            cyl("5", "12"),
            translated("5", "5", "-1"),
        ),
        15 => (
            cube("10"),
            origin.clone(),
            cyl("3", "12"),
            translated("0", "5", "-1"),
        ),
        16 => (
            cube("10"),
            origin.clone(),
            cyl("2", "4"),
            translated("5", "5", "3"),
        ),
        17 => (
            cube("10"),
            origin.clone(),
            ball("3"),
            translated("10", "10", "10"),
        ),
        18 => (
            cube("10"),
            origin.clone(),
            ball("3"),
            translated("5", "5", "10"),
        ),
        19 => (
            cube("10"),
            origin.clone(),
            ball("3"),
            translated("5", "5", "13"),
        ),
        20 => (
            cube("10"),
            origin.clone(),
            ball("3"),
            translated("5", "5", "7"),
        ),
        21 => (
            cube("10"),
            origin.clone(),
            ball("5"),
            translated("5", "5", "5"),
        ),
        22 => (
            cube("10"),
            origin.clone(),
            ball("3"),
            translated("10", "5", "5"),
        ),
        23 => (
            cyl("2", "10"),
            translated("0", "0", "-5"),
            cyl("1", "10"),
            xpose("-5", "0", "0"),
        ),
        24 => (
            cyl("1", "4"),
            translated("0", "0", "-2"),
            cyl("1", "4"),
            xpose("-2", "0", "0"),
        ),
        25 => (
            cyl("5", "10"),
            origin.clone(),
            cyl("5", "10"),
            translated("8", "0", "0"),
        ),
        26 => (
            cyl("5", "10"),
            origin.clone(),
            cyl("5", "10"),
            translated("10", "0", "0"),
        ),
        27 => (
            cyl("2", "5"),
            origin.clone(),
            cyl("2", "5"),
            translated("0", "0", "5"),
        ),
        28 => (
            cyl("2", "10"),
            translated("0", "0", "-5"),
            cyl("1", "10"),
            xpose("-5", "1", "0"),
        ),
        29 => (
            cyl("2", "10"),
            translated("0", "0", "-5"),
            cyl("1", "10"),
            xpose("-5", "3", "0"),
        ),
        30 => (
            cyl("3", "20"),
            translated("0", "0", "-10"),
            ball("6"),
            origin.clone(),
        ),
        31 => (
            cyl("5", "20"),
            translated("0", "0", "-10"),
            ball("5"),
            origin.clone(),
        ),
        32 => (
            cyl("1", "6"),
            translated("0", "0", "-3"),
            ball("2"),
            translated("1", "0", "0"),
        ),
        33 => (
            cyl("1", "6"),
            translated("0", "0", "-3"),
            ball("2"),
            translated("2", "0", "0"),
        ),
        34 => (
            cyl("1", "6"),
            translated("0", "0", "-3"),
            ball("2"),
            translated("3", "0", "0"),
        ),
        35 => (
            cyl("1", "6"),
            translated("0", "0", "-3"),
            ball("1"),
            translated("1", "0", "0"),
        ),
        36 => (
            ball("5"),
            origin.clone(),
            ball("5"),
            translated("6", "0", "0"),
        ),
        37 => (
            ball("5"),
            origin.clone(),
            ball("5"),
            translated("10", "0", "0"),
        ),
        38 => (
            ball("5"),
            origin.clone(),
            ball("2"),
            translated("3", "0", "0"),
        ),
        39 => (ball("5"), origin.clone(), ball("3"), origin.clone()),
        40 => (ball("5"), origin.clone(), ball("5"), origin.clone()),
        _ => panic!("unknown piece 5 fixture"),
    };
    let a = a.build(&mut m, &ap);
    let b = b.build(&mut m, &bp);
    (m, a, b)
}
// Finite solids whose side intersections are exactly piece 4's cases 17..26.
// Caps are beyond every branch. Cases 20,21,22,23 have no quartic branch at
// their base parameters; print zero spline counts instead of inventing one.
pub fn branch_pair(case: u32) -> (Model, BodyId, BodyId) {
    let mut m = Model::new();
    let (ar, br, bp, sphere_b) = match case {
        17 => ("2", "1", xpose("-20", "0", "0"), false),
        18 => ("2", "1", xpose("-20", "2", "0"), false),
        19 => ("2", "1", xpose("-20", "1", "0"), false),
        20 => ("2", "1", xpose("-20", "3", "0"), false),
        21 => ("1", "1", xpose("-20", "0", "0"), false),
        22 => {
            // q=(3,0,1,0) rotates z to (3/5,0,4/5).
            let rotation = Rotation::from_quaternion(n(3), n(0), n(1), n(0)).unwrap();
            (
                "1",
                "1",
                Pose::new(
                    rotation,
                    (RVector3::new(n(-12), n(0), n(-16))).into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
                false,
            )
        }
        23 => ("5", "6", translated("0", "0", "1"), true),
        24 => ("1", "2", translated("1", "0", "0"), true),
        25 => ("1", "2", translated("2", "0", "0"), true),
        26 => ("1", "3", translated("1", "0", "0"), true),
        _ => panic!(),
    };
    let a = cylinder(&mut m, r(ar), n(40), &translated("0", "0", "-20")).unwrap();
    let b = if sphere_b {
        sphere(&mut m, r(br), &bp).unwrap()
    } else {
        cylinder(&mut m, r(br), n(40), &bp).unwrap()
    };
    (m, a, b)
}
