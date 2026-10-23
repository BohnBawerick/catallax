use super::oracle::{Expr, from_kernel, pi, rat, sqrt};
use catallax_numbers::{Fact, FactKind, FactOrigin, MarkLog, Rational, Real};
use catallax_shapes::*;
use catallax_topology::*;

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
pub fn ten(k: i32) -> Rational {
    Rational::parse(&format!("1e{k}")).unwrap()
}
fn plane(origin: RPoint3, normal: RVector3, ref_dir: RVector3) -> Surface {
    Surface::Plane(Plane {
        origin: origin.into(),
        normal: normal.into(),
        ref_dir: ref_dir.into(),
    })
}
fn sphere_surface(radius: i64) -> Surface {
    Surface::Sphere(Sphere {
        center: p(0, 0, 0).into(),
        radius: r(radius).into(),
        axis: v(0, 0, 1).into(),
        ref_dir: v(1, 0, 0).into(),
    })
}
fn cylinder_surface(origin: RPoint3, axis: RVector3, radius: i64, ref_dir: RVector3) -> Surface {
    Surface::Cylinder(Cylinder {
        origin: origin.into(),
        axis: axis.into(),
        radius_sq: catallax_numbers::Real::from(r(radius)).square(),
        ref_dir: ref_dir.into(),
    })
}
fn circle(center: RPoint3, normal: RVector3, radius_sq: i64, ref_dir: RVector3) -> Curve {
    Curve::Circle(Circle {
        center: center.into(),
        normal: normal.into(),
        radius_sq: r(radius_sq).into(),
        ref_dir: ref_dir.into(),
    })
}
fn flip(s: Sense) -> Sense {
    match s {
        Sense::Same => Sense::Reversed,
        Sense::Reversed => Sense::Same,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Case {
    Block,
    FractionalBlock,
    Cylinder,
    Sphere,
    Frustum,
    Apex,
    Torus,
    Empty,
    Hole,
    Void,
    Hemisphere,
    Cap,
    Octant,
    HalfPositive,
    HalfNegative,
    Steinmetz,
}
pub const CASES: [Case; 16] = [
    Case::Block,
    Case::FractionalBlock,
    Case::Cylinder,
    Case::Sphere,
    Case::Frustum,
    Case::Apex,
    Case::Torus,
    Case::Empty,
    Case::Hole,
    Case::Void,
    Case::Hemisphere,
    Case::Cap,
    Case::Octant,
    Case::HalfPositive,
    Case::HalfNegative,
    Case::Steinmetz,
];
pub const CURVED: [Case; 11] = [
    Case::Cylinder,
    Case::Sphere,
    Case::Frustum,
    Case::Apex,
    Case::Torus,
    Case::Hemisphere,
    Case::Cap,
    Case::Octant,
    Case::HalfPositive,
    Case::HalfNegative,
    Case::Steinmetz,
];

#[derive(Clone, Debug)]
pub struct Expected {
    pub volume: Expr,
    pub area: Expr,
    pub centroid: Option<[Expr; 3]>,
    pub exact: bool,
}
impl Expected {
    pub fn moved_scaled(&self, scale: &Rational, pose: &Pose) -> Self {
        let s = from_kernel(scale);
        let centroid = self.centroid.as_ref().map(|c| {
            let t = pose
                .translation
                .rational_components()
                .expect("rational measurement fixture");
            let rotation = pose
                .rotation
                .rational_matrix()
                .expect("rational measurement fixture");
            std::array::from_fn(|i| {
                (0..3).fold(from_kernel(t[i]), |sum, j| {
                    sum + from_kernel(rotation[i][j]) * s.clone() * c[j].clone()
                })
            })
        });
        Self {
            volume: self.volume.clone() * s.clone() * s.clone() * s.clone(),
            area: self.area.clone() * s.clone() * s,
            centroid,
            exact: self.exact,
        }
    }
}
fn zero() -> [Expr; 3] {
    std::array::from_fn(|_| rat(0, 1))
}
impl Case {
    pub fn expected(self) -> Expected {
        let (volume, area, centroid, exact) = match self {
            Self::Block => (
                rat(30, 1),
                rat(62, 1),
                [rat(1, 1), rat(3, 2), rat(5, 2)],
                true,
            ),
            Self::FractionalBlock => (
                rat(7, 3000),
                rat(7022, 1500),
                [rat(1, 6), rat(7, 2), rat(1, 2000)],
                true,
            ),
            Self::Cylinder => (
                rat(20, 1) * pi(),
                rat(28, 1) * pi(),
                [rat(0, 1), rat(0, 1), rat(5, 2)],
                false,
            ),
            Self::Sphere => (rat(500, 3) * pi(), rat(100, 1) * pi(), zero(), false),
            Self::Frustum => (
                rat(14, 3) * pi(),
                (rat(5, 1) + rat(3, 1) * sqrt(5)) * pi(),
                [rat(0, 1), rat(0, 1), rat(11, 14)],
                false,
            ),
            Self::Apex => (
                rat(12, 1) * pi(),
                rat(24, 1) * pi(),
                [rat(0, 1), rat(0, 1), rat(1, 1)],
                false,
            ),
            Self::Torus => (
                rat(6, 1) * pi() * pi(),
                rat(12, 1) * pi() * pi(),
                zero(),
                false,
            ),
            Self::Empty => {
                return Expected {
                    volume: rat(0, 1),
                    area: rat(0, 1),
                    centroid: None,
                    exact: true,
                };
            }
            Self::Hole => (
                rat(64, 1) - rat(4, 1) * pi(),
                rat(96, 1) + rat(6, 1) * pi(),
                [rat(2, 1), rat(2, 1), rat(2, 1)],
                false,
            ),
            Self::Void => (
                rat(992, 1),
                rat(624, 1),
                [rat(5, 1), rat(5, 1), rat(5, 1)],
                true,
            ),
            Self::Hemisphere => (
                rat(16, 3) * pi(),
                rat(12, 1) * pi(),
                [rat(0, 1), rat(0, 1), rat(3, 4)],
                false,
            ),
            Self::Cap => (
                rat(52, 3) * pi(),
                rat(36, 1) * pi(),
                [rat(0, 1), rat(0, 1), rat(48, 13)],
                false,
            ),
            Self::Octant => (
                rat(9, 2) * pi(),
                rat(45, 4) * pi(),
                [rat(9, 8), rat(9, 8), rat(9, 8)],
                false,
            ),
            Self::HalfPositive | Self::HalfNegative => {
                let sign = if self == Self::HalfPositive { 1 } else { -1 };
                (
                    rat(6, 1) * pi(),
                    rat(10, 1) * pi() + rat(12, 1),
                    [rat(0, 1), rat(sign * 8, 3) / pi(), rat(3, 2)],
                    false,
                )
            }
            Self::Steinmetz => (rat(16, 3), rat(16, 1), zero(), false),
        };
        Expected {
            volume,
            area,
            centroid: Some(centroid),
            exact,
        }
    }
}

pub struct Fixture {
    pub model: Model,
    pub body: BodyId,
    pub expected: Expected,
    pub faces: Vec<(FaceId, Expr)>,
}

#[derive(Clone)]
struct FaceSpec {
    surface: Surface,
    sense: Sense,
    loops: Vec<Vec<(usize, Sense)>>,
    area: Expr,
}
#[derive(Default)]
struct Blueprint {
    points: Vec<RPoint3>,
    edges: Vec<(Curve, Option<usize>, Option<usize>)>,
    faces: Vec<FaceSpec>,
}
impl Blueprint {
    fn line(&mut self, a: usize, b: usize) -> usize {
        let pa = &self.points[a];
        let pb = &self.points[b];
        let curve = Curve::Line(Line {
            origin: pa.clone().into(),
            dir: RVector3::new(&pb.x - &pa.x, &pb.y - &pa.y, &pb.z - &pa.z).into(),
        });
        self.edge(curve, Some(a), Some(b))
    }
    fn edge(&mut self, curve: Curve, a: Option<usize>, b: Option<usize>) -> usize {
        let id = self.edges.len();
        self.edges.push((curve, a, b));
        id
    }
    fn face(&mut self, surface: Surface, loops: Vec<Vec<(usize, Sense)>>, area: Expr) {
        self.faces.push(FaceSpec {
            surface,
            sense: Sense::Same,
            loops,
            area,
        });
    }
    fn reverse(&mut self) {
        for f in &mut self.faces {
            f.sense = flip(f.sense);
            for l in &mut f.loops {
                l.reverse();
                for (_, s) in l {
                    *s = flip(*s);
                }
            }
        }
    }
    fn assemble(
        &self,
        m: &mut Model,
        scale: &Rational,
        pose: &Pose,
    ) -> (ShellId, Vec<(FaceId, Expr)>) {
        let point = |p: &RPoint3| pose.apply_point(&scaled_point(p, scale).into());
        let vs: Vec<_> = self
            .points
            .iter()
            .map(|p| {
                let point = m.geom_mut().add_point(PointDef::Explicit(point(p)));
                m.add_vertex(Vertex { point })
            })
            .collect();
        let es: Vec<_> = self
            .edges
            .iter()
            .map(|(c, a, b)| {
                let curve = m
                    .geom_mut()
                    .add_curve(scaled_curve(c, scale).transformed(pose).unwrap());
                for i in [a, b].into_iter().flatten() {
                    record(
                        m,
                        GeomKey::Point(m.vertex(vs[*i]).point),
                        GeomKey::Curve(curve),
                    );
                }
                m.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: a.map(|i| vs[i]),
                    end: b.map(|i| vs[i]),
                })
            })
            .collect();
        let fs: Vec<_> = self
            .faces
            .iter()
            .map(|f| {
                let surface = m
                    .geom_mut()
                    .add_surface(scaled_surface(&f.surface, scale).transformed(pose));
                let loops = f
                    .loops
                    .iter()
                    .map(|l| {
                        let cs = l
                            .iter()
                            .map(|(e, sense)| {
                                record(
                                    m,
                                    GeomKey::Curve(m.edge(es[*e]).curve),
                                    GeomKey::Surface(surface),
                                );
                                m.add_coedge(Coedge {
                                    edge: es[*e],
                                    sense: *sense,
                                })
                            })
                            .collect();
                        m.add_loop(Loop::Coedges(cs))
                    })
                    .collect();
                let id = m.add_face(Face {
                    surface,
                    sense: f.sense,
                    loops,
                });
                (id, f.area.clone() * from_kernel(scale) * from_kernel(scale))
            })
            .collect();
        (
            m.add_shell(Shell {
                faces: fs.iter().map(|(id, _)| *id).collect(),
            }),
            fs,
        )
    }
}
fn record(m: &mut Model, a: GeomKey, b: GeomKey) {
    m.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("measure.blind-fixture"),
    });
}
fn scaled_point(p: &RPoint3, s: &Rational) -> RPoint3 {
    RPoint3::new(&p.x * s, &p.y * s, &p.z * s)
}
fn scaled_definition_point(p: &Point3, s: &Real) -> Point3 {
    Point3::new(&p.x * s, &p.y * s, &p.z * s)
}
fn scaled_curve(c: &Curve, s: &Rational) -> Curve {
    let s = Real::from(s.clone());
    match c {
        Curve::Line(c) => Curve::Line(Line {
            origin: scaled_definition_point(&c.origin, &s),
            dir: c.dir.clone(),
        }),
        Curve::Circle(c) => Curve::Circle(Circle {
            center: scaled_definition_point(&c.center, &s),
            radius_sq: &(&c.radius_sq * &s) * &s,
            normal: c.normal.clone(),
            ref_dir: c.ref_dir.clone(),
        }),
        Curve::Ellipse(c) => Curve::Ellipse(Ellipse {
            center: scaled_definition_point(&c.center, &s),
            major_sq: &(&c.major_sq * &s) * &s,
            minor_sq: &(&c.minor_sq * &s) * &s,
            normal: c.normal.clone(),
            major_dir: c.major_dir.clone(),
        }),
        _ => unreachable!("fixtures only use analytic curves"),
    }
}
fn scaled_surface(f: &Surface, s: &Rational) -> Surface {
    let s = Real::from(s.clone());
    match f {
        Surface::Plane(f) => Surface::Plane(Plane {
            origin: scaled_definition_point(&f.origin, &s),
            normal: f.normal.clone(),
            ref_dir: f.ref_dir.clone(),
        }),
        Surface::Sphere(f) => Surface::Sphere(Sphere {
            center: scaled_definition_point(&f.center, &s),
            radius: &f.radius * &s,
            axis: f.axis.clone(),
            ref_dir: f.ref_dir.clone(),
        }),
        Surface::Cylinder(f) => Surface::Cylinder(Cylinder {
            origin: scaled_definition_point(&f.origin, &s),
            radius_sq: &f.radius_sq * &s.square(),
            axis: f.axis.clone(),
            ref_dir: f.ref_dir.clone(),
        }),
        _ => unreachable!("other surfaces come from primitive builders"),
    }
}
fn box_blueprint(side: i64, offset: i64) -> Blueprint {
    let mut d = Blueprint {
        points: (0..8)
            .map(|i| {
                p(
                    offset + if i & 4 != 0 { side } else { 0 },
                    offset + if i & 2 != 0 { side } else { 0 },
                    offset + if i & 1 != 0 { side } else { 0 },
                )
            })
            .collect(),
        ..Blueprint::default()
    };
    for i in 0..8 {
        for bit in [4, 2, 1] {
            if i & bit == 0 {
                d.line(i, i | bit);
            }
        }
    }
    // Each cycle is counter-clockwise as seen from outside the material.
    let cycles = [
        [0, 1, 3, 2],
        [4, 6, 7, 5],
        [0, 4, 5, 1],
        [2, 3, 7, 6],
        [0, 2, 6, 4],
        [1, 5, 7, 3],
    ];
    let normals = [
        v(-1, 0, 0),
        v(1, 0, 0),
        v(0, -1, 0),
        v(0, 1, 0),
        v(0, 0, -1),
        v(0, 0, 1),
    ];
    let refs = [
        v(0, 1, 0),
        v(0, 1, 0),
        v(0, 0, 1),
        v(0, 0, 1),
        v(1, 0, 0),
        v(1, 0, 0),
    ];
    for ((cycle, normal), ref_dir) in cycles.into_iter().zip(normals).zip(refs) {
        let uses = (0..4)
            .map(|i| {
                let (a, b) = (cycle[i], cycle[(i + 1) % 4]);
                let e = d
                    .edges
                    .iter()
                    .position(|(_, x, y)| {
                        (*x == Some(a) && *y == Some(b)) || (*x == Some(b) && *y == Some(a))
                    })
                    .unwrap();
                (
                    e,
                    if d.edges[e].1 == Some(a) {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    },
                )
            })
            .collect();
        d.face(
            plane(d.points[cycle[0]].clone(), normal, ref_dir),
            vec![uses],
            rat(side * side, 1),
        );
    }
    d
}
fn cap_blueprint(radius: i64, z: i64, ring_sq: i64) -> Blueprint {
    use Sense::{Reversed as R, Same as S};
    let mut d = Blueprint::default();
    let e = d.edge(
        circle(p(0, 0, z), v(0, 0, 1), ring_sq, v(1, 0, 0)),
        None,
        None,
    );
    // Curved area = 2*pi*R*(R-z); closing disc area = pi*(R^2-z^2).
    d.face(
        sphere_surface(radius),
        vec![vec![(e, S)]],
        rat(2 * radius * (radius - z), 1) * pi(),
    );
    d.face(
        plane(p(0, 0, z), v(0, 0, -1), v(1, 0, 0)),
        vec![vec![(e, R)]],
        rat(ring_sq, 1) * pi(),
    );
    d
}
fn octant_blueprint() -> Blueprint {
    use Sense::{Reversed as R, Same as S};
    let mut d = Blueprint {
        points: vec![p(0, 0, 0), p(3, 0, 0), p(0, 3, 0), p(0, 0, 3)],
        ..Blueprint::default()
    };
    let ox = d.line(0, 1);
    let oy = d.line(0, 2);
    let oz = d.line(0, 3);
    let xy = d.edge(
        circle(p(0, 0, 0), v(0, 0, 1), 9, v(1, 0, 0)),
        Some(1),
        Some(2),
    );
    let yz = d.edge(
        circle(p(0, 0, 0), v(1, 0, 0), 9, v(0, 1, 0)),
        Some(2),
        Some(3),
    );
    let zx = d.edge(
        circle(p(0, 0, 0), v(0, 1, 0), 9, v(0, 0, 1)),
        Some(3),
        Some(1),
    );
    d.face(
        plane(p(0, 0, 0), v(-1, 0, 0), v(0, 1, 0)),
        vec![vec![(oz, S), (yz, R), (oy, R)]],
        rat(9, 4) * pi(),
    );
    d.face(
        plane(p(0, 0, 0), v(0, -1, 0), v(1, 0, 0)),
        vec![vec![(ox, S), (zx, R), (oz, R)]],
        rat(9, 4) * pi(),
    );
    d.face(
        plane(p(0, 0, 0), v(0, 0, -1), v(1, 0, 0)),
        vec![vec![(oy, S), (xy, R), (ox, R)]],
        rat(9, 4) * pi(),
    );
    d.face(
        sphere_surface(3),
        vec![vec![(xy, S), (yz, S), (zx, S)]],
        rat(9, 2) * pi(),
    );
    d
}
fn half_blueprint(sign: i64) -> Blueprint {
    use Sense::{Reversed as R, Same as S};
    let mut d = Blueprint {
        points: vec![
            p(2 * sign, 0, 0),
            p(-2 * sign, 0, 0),
            p(2 * sign, 0, 3),
            p(-2 * sign, 0, 3),
        ],
        ..Blueprint::default()
    };
    let bottom = d.edge(
        circle(p(0, 0, 0), v(0, 0, 1), 4, v(1, 0, 0)),
        Some(0),
        Some(1),
    );
    let top = d.edge(
        circle(p(0, 0, 3), v(0, 0, 1), 4, v(1, 0, 0)),
        Some(2),
        Some(3),
    );
    let ac = d.line(0, 2);
    let bd = d.line(1, 3);
    let ba = d.line(1, 0);
    let dc = d.line(3, 2);
    d.face(
        cylinder_surface(p(0, 0, 0), v(0, 0, 1), 2, v(1, 0, 0)),
        vec![vec![(bottom, S), (bd, S), (top, R), (ac, R)]],
        rat(6, 1) * pi(),
    );
    d.face(
        plane(p(0, 0, 0), v(0, -sign, 0), v(1, 0, 0)),
        vec![vec![(ac, S), (dc, R), (bd, R), (ba, S)]],
        rat(12, 1),
    );
    d.face(
        plane(p(0, 0, 0), v(0, 0, -1), v(1, 0, 0)),
        vec![vec![(bottom, R), (ba, R)]],
        rat(2, 1) * pi(),
    );
    d.face(
        plane(p(0, 0, 3), v(0, 0, 1), v(1, 0, 0)),
        vec![vec![(top, S), (dc, S)]],
        rat(2, 1) * pi(),
    );
    d
}
fn steinmetz_blueprint() -> Blueprint {
    use Sense::{Reversed as R, Same as S};
    let mut d = Blueprint {
        points: vec![p(0, 1, 0), p(0, -1, 0)],
        ..Blueprint::default()
    };
    // z=+/-x. Semi-axes sqrt(2),1; normal=(-sign,0,1).
    // Natural arc y+ -> y- has x<0, reverse endpoint order has x>0.
    for sign in [1, -1] {
        let c = Curve::Ellipse(Ellipse {
            center: p(0, 0, 0).into(),
            normal: v(-sign, 0, 1).into(),
            major_dir: v(1, 0, sign).into(),
            major_sq: r(2).into(),
            minor_sq: r(1).into(),
        });
        d.edge(c.clone(), Some(0), Some(1));
        d.edge(c, Some(1), Some(0));
    }
    let zc = cylinder_surface(p(0, 0, 0), v(0, 0, 1), 1, v(1, 0, 0));
    let xc = cylinder_surface(p(0, 0, 0), v(1, 0, 0), 1, v(0, 0, -1));
    // Each cylindrical patch has area integral_-pi/2^pi/2 2*cos(u) du = 4.
    d.face(zc.clone(), vec![vec![(3, S), (1, R)]], rat(4, 1));
    d.face(zc, vec![vec![(0, S), (2, R)]], rat(4, 1));
    d.face(xc.clone(), vec![vec![(2, S), (1, S)]], rat(4, 1));
    d.face(xc, vec![vec![(0, R), (3, R)]], rat(4, 1));
    d
}

pub fn build(case: Case, scale: &Rational, pose: &Pose) -> Fixture {
    let mut m = Model::new();
    let expected = case.expected().moved_scaled(scale, pose);
    let (body, faces) = match case {
        Case::Hole
        | Case::Void
        | Case::Hemisphere
        | Case::Cap
        | Case::Octant
        | Case::HalfPositive
        | Case::HalfNegative
        | Case::Steinmetz => {
            let d = match case {
                Case::Hole => {
                    let mut d = box_blueprint(4, 0);
                    let bottom = d.edge(circle(p(2, 2, 0), v(0, 0, 1), 1, v(1, 0, 0)), None, None);
                    let top = d.edge(circle(p(2, 2, 4), v(0, 0, 1), 1, v(1, 0, 0)), None, None);
                    d.faces[4].loops.push(vec![(bottom, Sense::Same)]);
                    d.faces[5].loops.push(vec![(top, Sense::Reversed)]);
                    for i in [4, 5] {
                        d.faces[i].area = rat(16, 1) - pi();
                    }
                    d.faces.push(FaceSpec {
                        surface: cylinder_surface(p(2, 2, 0), v(0, 0, 1), 1, v(1, 0, 0)),
                        sense: Sense::Reversed,
                        loops: vec![vec![(bottom, Sense::Reversed)], vec![(top, Sense::Same)]],
                        area: rat(8, 1) * pi(),
                    });
                    d
                }
                Case::Void => box_blueprint(10, 0),
                Case::Hemisphere => cap_blueprint(2, 0, 4),
                Case::Cap => cap_blueprint(5, 3, 16),
                Case::Octant => octant_blueprint(),
                Case::HalfPositive => half_blueprint(1),
                Case::HalfNegative => half_blueprint(-1),
                Case::Steinmetz => steinmetz_blueprint(),
                _ => unreachable!(),
            };
            let (shell, mut faces) = d.assemble(&mut m, scale, pose);
            let mut shells = vec![shell];
            if case == Case::Void {
                let mut inner = box_blueprint(2, 4);
                inner.reverse();
                let (shell, fs) = inner.assemble(&mut m, scale, pose);
                shells.push(shell);
                faces.extend(fs);
            }
            let lump = m.add_lump(Lump { shells });
            let body = m.add_body(Body {
                kind: BodyKind::Solid,
                lumps: vec![lump],
                marks: MarkLog::new(),
            });
            (body, faces)
        }
        _ => {
            let s = |n| &r(n) * scale;
            let body = match case {
                Case::Block => block(&mut m, [s(2), s(3), s(5)], pose).unwrap(),
                Case::FractionalBlock => {
                    block(&mut m, [&q(1, 3) * scale, s(7), &q(1, 1000) * scale], pose).unwrap()
                }
                Case::Cylinder => cylinder(&mut m, s(2), s(5), pose).unwrap(),
                Case::Sphere => sphere(&mut m, s(5), pose).unwrap(),
                Case::Frustum => cone(&mut m, s(2), s(1), s(2), pose).unwrap(),
                Case::Apex => cone(&mut m, s(3), r(0), s(4), pose).unwrap(),
                Case::Torus => torus(&mut m, s(3), s(1), pose).unwrap(),
                Case::Empty => empty(&mut m),
                _ => unreachable!(),
            };
            let faces = m
                .faces_of_body(body)
                .into_iter()
                .map(|id| {
                    let surf = m.geom().surface(m.face(id).surface);
                    let a = match (case, surf) {
                        (Case::Block | Case::FractionalBlock, Surface::Plane(f)) => {
                            let n = pose
                                .rotation
                                .inverse()
                                .apply(&f.normal)
                                .as_rational()
                                .expect("rational measurement fixture");
                            let size = if case == Case::Block {
                                [rat(2, 1), rat(3, 1), rat(5, 1)]
                            } else {
                                [rat(1, 3), rat(7, 1), rat(1, 1000)]
                            };
                            if n.x != r(0) {
                                size[1].clone() * size[2].clone()
                            } else if n.y != r(0) {
                                size[0].clone() * size[2].clone()
                            } else {
                                size[0].clone() * size[1].clone()
                            }
                        }
                        (Case::Cylinder, Surface::Plane(_)) => rat(4, 1) * pi(),
                        (Case::Cylinder, Surface::Cylinder(_)) => rat(20, 1) * pi(),
                        (Case::Sphere, _) => rat(100, 1) * pi(),
                        (Case::Torus, _) => rat(12, 1) * pi() * pi(),
                        (Case::Frustum, Surface::Plane(f)) => {
                            let n = pose
                                .rotation
                                .inverse()
                                .apply(&f.normal)
                                .as_rational()
                                .expect("rational measurement fixture");
                            if n.z > r(0) { pi() } else { rat(4, 1) * pi() }
                        }
                        (Case::Frustum, Surface::Cone(_)) => rat(3, 1) * sqrt(5) * pi(),
                        (Case::Apex, Surface::Plane(_)) => rat(9, 1) * pi(),
                        (Case::Apex, Surface::Cone(_)) => rat(15, 1) * pi(),
                        _ => unreachable!("unexpected primitive face"),
                    };
                    (id, a * from_kernel(scale) * from_kernel(scale))
                })
                .collect();
            (body, faces)
        }
    };
    assert!(m.check(body).is_empty(), "{case:?}: {:?}", m.check(body));
    Fixture {
        model: m,
        body,
        expected,
        faces,
    }
}

pub fn poses() -> Vec<Pose> {
    // Fixed xorshift seed. Quaternion components lie in [-100,100], translations
    // are rational with magnitude <= 10^6. No platform RNG or floating point.
    let mut state = 0x91a7_60c4_f033_e5b1u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    (0..20)
        .map(|_| {
            let mut c: [i64; 4] = std::array::from_fn(|_| (next() % 201) as i64 - 100);
            if c == [0; 4] {
                c[0] = 1;
            }
            let t: [Rational; 3] = std::array::from_fn(|_| {
                q(
                    (next() % 2_000_001) as i64 - 1_000_000,
                    (next() % 17 + 1) as i64,
                )
            });
            Pose::new(
                Rotation::from_quaternion(r(c[0]), r(c[1]), r(c[2]), r(c[3])).unwrap(),
                (RVector3::new(t[0].clone(), t[1].clone(), t[2].clone())).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap()
        })
        .collect()
}

pub fn thin_spherical_shell(gap: &Rational) -> Fixture {
    let mut model = Model::new();
    let inner = &r(1) - gap;
    let mut shells = vec![];
    let mut faces = vec![];
    for (radius, sense) in [(r(1), Sense::Same), (inner.clone(), Sense::Reversed)] {
        let surface = model.geom_mut().add_surface(Surface::Sphere(Sphere {
            center: p(0, 0, 0).into(),
            radius: radius.clone().into(),
            axis: v(0, 0, 1).into(),
            ref_dir: v(1, 0, 0).into(),
        }));
        let face = model.add_face(Face {
            surface,
            sense,
            loops: vec![],
        });
        shells.push(model.add_shell(Shell { faces: vec![face] }));
        faces.push((
            face,
            rat(4, 1) * pi() * from_kernel(&radius) * from_kernel(&radius),
        ));
    }
    let lump = model.add_lump(Lump { shells });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    let ri = from_kernel(&inner);
    let expected = Expected {
        volume: rat(4, 3) * pi() * (rat(1, 1) - ri.clone() * ri.clone() * ri.clone()),
        area: rat(4, 1) * pi() * (rat(1, 1) + ri.clone() * ri),
        centroid: Some(zero()),
        exact: false,
    };
    assert!(model.check(body).is_empty());
    Fixture {
        model,
        body,
        expected,
        faces,
    }
}
