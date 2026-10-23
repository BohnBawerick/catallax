#![allow(dead_code, unused_imports, unused_macros)]

pub use catallax_numbers::{Fact, FactKind, FactOrigin, MarkLog, Rational, Sign};
pub use catallax_shapes::*;
pub use catallax_topology::*;
use std::collections::{BTreeMap, BTreeSet};

macro_rules! blind_test {
    ($name:ident, $body:block) => {
        #[test]
        fn $name() $body
    };
}
pub(crate) use blind_test;

pub fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
pub fn q(n: i64, d: i64) -> Rational {
    Rational::from_ratio(n, d).unwrap()
}
pub fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}
pub fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(r(x), r(y), r(z))
}
pub fn flip(s: Sense) -> Sense {
    match s {
        Sense::Same => Sense::Reversed,
        Sense::Reversed => Sense::Same,
    }
}
pub fn plane(origin: RPoint3, normal: RVector3, ref_dir: RVector3) -> Surface {
    Surface::Plane(Plane {
        origin: origin.into(),
        normal: normal.into(),
        ref_dir: ref_dir.into(),
    })
}
pub fn circle(center: RPoint3, radius: Rational) -> Curve {
    Curve::Circle(Circle {
        center: center.into(),
        normal: v(0, 0, 1).into(),
        radius_sq: (&radius * &radius).into(),
        ref_dir: v(1, 0, 0).into(),
    })
}
pub fn sphere_surface(radius: Rational) -> Surface {
    Surface::Sphere(Sphere {
        center: p(0, 0, 0).into(),
        radius: radius.into(),
        axis: v(0, 0, 1).into(),
        ref_dir: v(1, 0, 0).into(),
    })
}
pub fn torus_surface() -> Surface {
    Surface::Torus(Torus {
        center: p(0, 0, 0).into(),
        major: r(3).into(),
        minor: r(1).into(),
        axis: v(0, 0, 1).into(),
        ref_dir: v(1, 0, 0).into(),
    })
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
    pub fn build(&self, m: &mut Model, pose: &Pose) -> Result<BodyId, catallax_numbers::Refusal> {
        match self {
            Self::Block(s) => block(m, s.clone(), pose),
            Self::Cylinder(a, b) => cylinder(m, a.clone(), b.clone(), pose),
            Self::Sphere(a) => sphere(m, a.clone(), pose),
            Self::Cone(a, b, c) => cone(m, a.clone(), b.clone(), c.clone(), pose),
            Self::Torus(a, b) => torus(m, a.clone(), b.clone(), pose),
        }
    }
    pub fn name(&self) -> &'static str {
        match self {
            Self::Block(_) => "topology.block",
            Self::Cylinder(..) => "topology.cylinder",
            Self::Sphere(_) => "topology.sphere",
            Self::Cone(..) => "topology.cone",
            Self::Torus(..) => "topology.torus",
        }
    }
    pub fn scaled(&self, k: &Rational) -> Self {
        match self {
            Self::Block(s) => Self::Block(std::array::from_fn(|i| &s[i] * k)),
            Self::Cylinder(a, b) => Self::Cylinder(a * k, b * k),
            Self::Sphere(a) => Self::Sphere(a * k),
            Self::Cone(a, b, c) => Self::Cone(a * k, b * k, c * k),
            Self::Torus(a, b) => Self::Torus(a * k, b * k),
        }
    }
}
pub fn primitives() -> Vec<Primitive> {
    vec![
        Primitive::Block([r(2), r(3), r(5)]),
        Primitive::Cylinder(r(2), r(5)),
        Primitive::Sphere(r(3)),
        Primitive::Cone(r(2), r(1), r(2)),
        Primitive::Cone(r(1), r(3), r(2)),
        Primitive::Cone(r(3), r(0), r(4)),
        Primitive::Torus(r(3), r(1)),
    ]
}

// Indices are fixture-local, never presumed to be an implementation's arena order.
#[derive(Clone, Debug)]
pub enum Boundary {
    Cycle(Vec<(usize, Sense)>),
    Vertex(usize),
}
#[derive(Clone, Debug)]
pub struct Expected {
    pub points: Vec<RPoint3>,
    pub edges: Vec<(Curve, Option<usize>, Option<usize>)>,
    pub faces: Vec<(Surface, Sense, Vec<Boundary>)>,
    pub interior: RPoint3,
    pub chi: i64,
}
pub const BOX_CYCLES: [[usize; 4]; 6] = [
    [0, 1, 3, 2],
    [4, 6, 7, 5],
    [0, 4, 5, 1],
    [2, 3, 7, 6],
    [0, 2, 6, 4],
    [1, 5, 7, 3],
];
pub fn expected(shape: &Primitive, pose: &Pose) -> Expected {
    let mut out = Expected {
        points: vec![],
        edges: vec![],
        faces: vec![],
        interior: p(0, 0, 0),
        chi: 2,
    };
    let cap = |z: Rational, sign| plane(RPoint3::new(r(0), r(0), z), v(0, 0, sign), v(1, 0, 0));
    match shape {
        Primitive::Block(s) => {
            out.points = (0..8)
                .map(|i| {
                    RPoint3::new(
                        if i & 4 != 0 { s[0].clone() } else { r(0) },
                        if i & 2 != 0 { s[1].clone() } else { r(0) },
                        if i & 1 != 0 { s[2].clone() } else { r(0) },
                    )
                })
                .collect();
            for i in 0..8 {
                for (bit, dir) in [(4, v(1, 0, 0)), (2, v(0, 1, 0)), (1, v(0, 0, 1))] {
                    if i & bit == 0 {
                        out.edges.push((
                            Curve::Line(Line {
                                origin: out.points[i].clone().into(),
                                dir: dir.into(),
                            }),
                            Some(i),
                            Some(i | bit),
                        ));
                    }
                }
            }
            let surfaces = [
                plane(p(0, 0, 0), v(-1, 0, 0), v(0, 1, 0)),
                plane(
                    RPoint3::new(s[0].clone(), r(0), r(0)),
                    v(1, 0, 0),
                    v(0, 1, 0),
                ),
                plane(p(0, 0, 0), v(0, -1, 0), v(0, 0, 1)),
                plane(
                    RPoint3::new(r(0), s[1].clone(), r(0)),
                    v(0, 1, 0),
                    v(0, 0, 1),
                ),
                cap(r(0), -1),
                cap(s[2].clone(), 1),
            ];
            for (surface, corners) in surfaces.into_iter().zip(BOX_CYCLES) {
                let cycle = (0..4)
                    .map(|j| {
                        let (a, b) = (corners[j], corners[(j + 1) % 4]);
                        let e = out
                            .edges
                            .iter()
                            .position(|(_, x, y)| {
                                (*x == Some(a) && *y == Some(b)) || (*x == Some(b) && *y == Some(a))
                            })
                            .unwrap();
                        (
                            e,
                            if out.edges[e].1 == Some(a) {
                                Sense::Same
                            } else {
                                Sense::Reversed
                            },
                        )
                    })
                    .collect();
                out.faces
                    .push((surface, Sense::Same, vec![Boundary::Cycle(cycle)]));
            }
            out.interior = RPoint3::new(&s[0] / &r(2), &s[1] / &r(2), &s[2] / &r(2));
        }
        Primitive::Cylinder(radius, height) | Primitive::Cone(radius, _, height) => {
            let top = match shape {
                Primitive::Cone(_, t, _) => t.clone(),
                _ => radius.clone(),
            };
            out.edges
                .push((circle(p(0, 0, 0), radius.clone()), None, None));
            out.faces.push((
                cap(r(0), -1),
                Sense::Same,
                vec![Boundary::Cycle(vec![(0, Sense::Reversed)])],
            ));
            let mut side_loops = vec![Boundary::Cycle(vec![(0, Sense::Same)])];
            if !top.is_zero() {
                out.edges.push((
                    circle(RPoint3::new(r(0), r(0), height.clone()), top.clone()),
                    None,
                    None,
                ));
                out.faces.push((
                    cap(height.clone(), 1),
                    Sense::Same,
                    vec![Boundary::Cycle(vec![(1, Sense::Same)])],
                ));
                side_loops.push(Boundary::Cycle(vec![(1, Sense::Reversed)]));
            } else {
                out.points.push(RPoint3::new(r(0), r(0), height.clone()));
                side_loops.push(Boundary::Vertex(0));
            }
            let side = if matches!(shape, Primitive::Cylinder(..)) {
                Surface::Cylinder(Cylinder {
                    origin: p(0, 0, 0).into(),
                    axis: v(0, 0, 1).into(),
                    radius_sq: catallax_numbers::Real::from(radius.clone()).square(),
                    ref_dir: v(1, 0, 0).into(),
                })
            } else {
                // r(z)=base+(top-base)z/h vanishes at h*base/(base-top).
                Surface::Cone(Cone {
                    apex: RPoint3::new(r(0), r(0), height * radius / (radius - &top)).into(),
                    axis: v(0, 0, if top < *radius { -1 } else { 1 }).into(),
                    slope: (&(radius - &top).abs() / height).into(),
                    ref_dir: v(1, 0, 0).into(),
                })
            };
            out.faces.push((side, Sense::Same, side_loops));
            out.interior = RPoint3::new(r(0), r(0), height / &r(2));
        }
        Primitive::Sphere(radius) => {
            out.faces
                .push((sphere_surface(radius.clone()), Sense::Same, vec![]))
        }
        Primitive::Torus(major, minor) => {
            out.faces.push((
                Surface::Torus(Torus {
                    center: p(0, 0, 0).into(),
                    axis: v(0, 0, 1).into(),
                    major: major.clone().into(),
                    minor: minor.clone().into(),
                    ref_dir: v(1, 0, 0).into(),
                }),
                Sense::Same,
                vec![],
            ));
            out.interior = RPoint3::new(major.clone(), r(0), r(0));
            out.chi = 0;
        }
    }
    for point in &mut out.points {
        *point = pose
            .apply_point(&point.to_real())
            .as_rational()
            .expect("rational topology fixture");
    }
    for (curve, _, _) in &mut out.edges {
        *curve = curve.transformed(pose).unwrap();
    }
    for (surface, _, _) in &mut out.faces {
        *surface = surface.transformed(pose);
    }
    out.interior = pose
        .apply_point(&out.interior.to_real())
        .as_rational()
        .expect("rational topology fixture");
    out
}

pub fn record(m: &mut Model, kind: FactKind, a: GeomKey, b: GeomKey) {
    m.geom_mut().record_fact(Fact {
        kind,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("blind.fixture"),
    });
}
pub fn vertex(m: &mut Model, point: RPoint3) -> VertexId {
    let point = m.geom_mut().add_point(PointDef::Explicit((point).into()));
    m.add_vertex(Vertex { point })
}
pub fn edge(m: &mut Model, curve: Curve, start: Option<VertexId>, end: Option<VertexId>) -> EdgeId {
    let curve = m.geom_mut().add_curve(curve);
    for vertex in [start, end].into_iter().flatten() {
        record(
            m,
            FactKind::Incident,
            GeomKey::Point(m.vertex(vertex).point),
            GeomKey::Curve(curve),
        );
    }
    m.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start,
        end,
    })
}
pub fn cycle(m: &mut Model, uses: &[(EdgeId, Sense)]) -> LoopId {
    let cs = uses
        .iter()
        .map(|&(edge, sense)| m.add_coedge(Coedge { edge, sense }))
        .collect();
    m.add_loop(Loop::Coedges(cs))
}
pub fn face(m: &mut Model, surface: Surface, sense: Sense, loops: Vec<LoopId>) -> FaceId {
    let surface = m.geom_mut().add_surface(surface);
    face_on(m, surface, sense, loops, true)
}
pub fn face_on(
    m: &mut Model,
    surface: SurfaceId,
    sense: Sense,
    loops: Vec<LoopId>,
    facts: bool,
) -> FaceId {
    if facts {
        for &l in &loops {
            match m.loop_(l).clone() {
                Loop::Vertex(v) => record(
                    m,
                    FactKind::Incident,
                    GeomKey::Point(m.vertex(v).point),
                    GeomKey::Surface(surface),
                ),
                Loop::Coedges(cs) => {
                    for c in cs {
                        record(
                            m,
                            FactKind::Incident,
                            GeomKey::Curve(m.edge(m.coedge(c).edge).curve),
                            GeomKey::Surface(surface),
                        );
                    }
                }
            }
        }
    }
    m.add_face(Face {
        surface,
        sense,
        loops,
    })
}
pub fn body_from_lumps(m: &mut Model, lumps: Vec<LumpId>) -> BodyId {
    m.add_body(Body {
        kind: BodyKind::Solid,
        lumps,
        marks: MarkLog::new(),
    })
}
pub fn body_from_shells(m: &mut Model, shells: Vec<ShellId>) -> BodyId {
    let lump = m.add_lump(Lump { shells });
    body_from_lumps(m, vec![lump])
}
pub fn body_from_faces(m: &mut Model, faces: Vec<FaceId>) -> BodyId {
    let shell = m.add_shell(Shell { faces });
    body_from_shells(m, vec![shell])
}
pub fn assemble(m: &mut Model, e: &Expected) -> BodyId {
    let vs: Vec<_> = e.points.iter().map(|p| vertex(m, p.clone())).collect();
    let es: Vec<_> = e
        .edges
        .iter()
        .map(|(c, a, b)| edge(m, c.clone(), a.map(|i| vs[i]), b.map(|i| vs[i])))
        .collect();
    let mut fs = vec![];
    for (s, sense, ls) in &e.faces {
        let loops = ls
            .iter()
            .map(|l| match l {
                Boundary::Vertex(i) => m.add_loop(Loop::Vertex(vs[*i])),
                Boundary::Cycle(cs) => {
                    cycle(m, &cs.iter().map(|(i, s)| (es[*i], *s)).collect::<Vec<_>>())
                }
            })
            .collect();
        fs.push(face(m, s.clone(), *sense, loops));
    }
    body_from_faces(m, fs)
}
pub fn explicit(m: &Model, id: VertexId) -> RPoint3 {
    match m.geom().point(m.vertex(id).point) {
        PointDef::Explicit(p) => p.as_rational().expect("rational topology fixture"),
        _ => panic!("builder must store an explicit point"),
    }
}
pub fn cs(m: &Model, l: LoopId) -> Vec<CoedgeId> {
    match m.loop_(l) {
        Loop::Coedges(c) => c.clone(),
        _ => panic!("expected coedge loop"),
    }
}
pub fn shell_of(m: &Model, b: BodyId) -> ShellId {
    m.lump(m.body(b).lumps[0]).shells[0]
}
pub fn sound(m: &Model, b: BodyId) {
    assert_eq!(m.check(b), vec![], "body {b:?}");
}
pub fn faults(m: &Model, b: BodyId, codes: &[&str]) {
    let actual = m.check(b);
    assert_eq!(actual, m.check(b), "fault order must be deterministic");
    for code in codes {
        assert!(
            actual.iter().any(|f| f.code == *code),
            "missing {code}: {actual:?}"
        );
    }
}

// Clone probes reveal the next index of all eight private topology arenas without
// changing the observed model or relying on Debug formatting/private fields.
pub fn lengths(m: &Model) -> [usize; 12] {
    let mut probe = m.clone();
    let s = probe.geom_mut().add_surface(sphere_surface(r(1)));
    let c = probe.geom_mut().add_curve(circle(p(0, 0, 0), r(1)));
    let p = probe
        .geom_mut()
        .add_point(PointDef::Explicit((p(0, 0, 0)).into()));
    let v = probe.add_vertex(Vertex { point: p });
    let e = probe.add_edge(Edge {
        curve: c,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let co = probe.add_coedge(Coedge {
        edge: e,
        sense: Sense::Same,
    });
    let l = probe.add_loop(Loop::Coedges(vec![co]));
    let f = probe.add_face(Face {
        surface: s,
        sense: Sense::Same,
        loops: vec![l],
    });
    let sh = probe.add_shell(Shell { faces: vec![f] });
    let lu = probe.add_lump(Lump { shells: vec![sh] });
    let b = body_from_lumps(&mut probe, vec![lu]);
    [
        m.geom().surface_count(),
        m.geom().curve_count(),
        m.geom().point_count(),
        v.index() as usize,
        e.index() as usize,
        co.index() as usize,
        l.index() as usize,
        f.index() as usize,
        sh.index() as usize,
        lu.index() as usize,
        b.index() as usize,
        m.geom().facts().len(),
    ]
}

pub fn walk(m: &Model, b: BodyId) -> (Vec<FaceId>, Vec<EdgeId>, Vec<VertexId>) {
    let (mut fs, mut es, mut vs) = (vec![], vec![], vec![]);
    let (mut sf, mut se, mut sv) = (BTreeSet::new(), BTreeSet::new(), BTreeSet::new());
    for &lump in &m.body(b).lumps {
        for &shell in &m.lump(lump).shells {
            for &f in &m.shell(shell).faces {
                if sf.insert(f) {
                    fs.push(f);
                }
                for &l in &m.face(f).loops {
                    match m.loop_(l) {
                        Loop::Vertex(v) => {
                            if sv.insert(*v) {
                                vs.push(*v);
                            }
                        }
                        Loop::Coedges(cs) => {
                            for &c in cs {
                                let co = m.coedge(c);
                                let e = m.edge(co.edge);
                                if se.insert(co.edge) {
                                    es.push(co.edge);
                                }
                                let ends = if co.sense == Sense::Same {
                                    [e.start, e.end]
                                } else {
                                    [e.end, e.start]
                                };
                                for v in ends.into_iter().flatten() {
                                    if sv.insert(v) {
                                        vs.push(v);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    (fs, es, vs)
}
pub fn euler(m: &Model, s: ShellId) -> i64 {
    let (mut edges, mut vertices) = (BTreeSet::new(), BTreeSet::new());
    let mut chi = 0;
    for &f in &m.shell(s).faces {
        let f = m.face(f);
        chi += if f.loops.is_empty() && m.geom().surface(f.surface).kind() == SurfaceKind::Torus {
            0
        } else {
            2 - f.loops.len() as i64
        };
        for &l in &f.loops {
            match m.loop_(l) {
                Loop::Vertex(v) => {
                    vertices.insert(*v);
                }
                Loop::Coedges(cs) => {
                    for &c in cs {
                        edges.insert(m.coedge(c).edge);
                    }
                }
            }
        }
    }
    let mut rings = 0;
    for &e in &edges {
        let e = m.edge(e);
        if e.start.is_none() && e.end.is_none() {
            rings += 1;
        }
        for v in [e.start, e.end].into_iter().flatten() {
            vertices.insert(v);
        }
    }
    vertices.len() as i64 + rings - edges.len() as i64 + chi
}

// The implicit formulas are evaluated here, not by Surface::implicit or side().
pub fn implicit(s: &Surface, p: &RPoint3) -> Rational {
    let sq = |x: &Rational| x * x;
    match s {
        Surface::Plane(s) => p
            .minus(&s.origin.as_rational().unwrap())
            .dot(&s.normal.as_rational().unwrap()),
        Surface::Cylinder(s) => {
            let d = p.minus(&s.origin.as_rational().unwrap());
            let axis = s.axis.as_rational().unwrap();
            d.cross(&axis).norm_squared() - s.radius_sq.as_rational().unwrap() * axis.norm_squared()
        }
        Surface::Cone(s) => {
            let d = p.minus(&s.apex.as_rational().unwrap());
            let axis = s.axis.as_rational().unwrap();
            d.cross(&axis).norm_squared() - sq(&s.slope.as_rational().unwrap()) * sq(&d.dot(&axis))
        }
        Surface::Sphere(s) => {
            p.minus(&s.center.as_rational().unwrap()).norm_squared()
                - sq(&s.radius.as_rational().unwrap())
        }
        Surface::Torus(s) => {
            let d = p.minus(&s.center.as_rational().unwrap());
            let axis = s.axis.as_rational().unwrap();
            let norm = d.norm_squared();
            let major = sq(&s.major.as_rational().unwrap());
            sq(&(&norm + &major - sq(&s.minor.as_rational().unwrap())))
                - r(4) * major * (norm - sq(&d.dot(&axis)) / axis.norm_squared())
        }
        _ => panic!("new surface requires an oracle"),
    }
}

pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }
    pub fn rational(&mut self) -> Rational {
        q((self.next() % 19) as i64 - 9, (self.next() % 7 + 1) as i64)
    }
    pub fn pose(&mut self) -> Pose {
        // w is strictly positive, so this quaternion can never be zero.
        let w = q((self.next() % 9 + 1) as i64, (self.next() % 7 + 1) as i64);
        Pose::new(
            Rotation::from_quaternion(w, self.rational(), self.rational(), self.rational())
                .unwrap(),
            (RVector3::new(self.rational(), self.rational(), self.rational())).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap()
    }
}
