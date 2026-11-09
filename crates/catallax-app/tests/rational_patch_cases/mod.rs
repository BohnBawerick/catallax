//! Closed named solids with independently known geometry, not recipe substitutes.

use catallax_numbers::{Decider, Fact, FactKind, FactOrigin, MarkLog, Rational, Real};
use catallax_recipe::{ElementId, Origin, Role, StepId};
use catallax_shapes::{
    BSplineSurface, Curve, GeomKey, Line, Plane, Point3, PointDef, RationalPatch, SplineBoundary,
    Surface, Vector3,
};
use catallax_topology::{
    Body, BodyId, BodyKind, Coedge, Edge, EdgeId, Face, FaceId, Loop, Lump, Model, Sense, Shell,
    Vertex, VertexId,
};

pub fn r(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
pub fn p(x: &str, y: &str, z: &str) -> Point3 {
    Point3::new(r(x).into(), r(y).into(), r(z).into())
}

#[derive(Clone, Copy, Debug)]
pub enum Roof {
    Bilinear,
    Rational,
    NonuniformBump,
    JoinedDegreeElevated,
    JoinedNonlinear,
}

pub struct Fixture {
    pub model: Model,
    pub body: BodyId,
    pub roofs: Vec<FaceId>,
    pub patches: Vec<RationalPatch>,
}

fn record(model: &mut Model, a: GeomKey, b: GeomKey) {
    model.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("rational-patch.fixture"),
    });
}

fn edge(model: &mut Model, curve: Curve, start: VertexId, end: VertexId) -> EdgeId {
    let curve = model.geom_mut().add_curve(curve);
    for vertex in [start, end] {
        record(
            model,
            GeomKey::Point(model.vertex(vertex).point),
            GeomKey::Curve(curve),
        );
    }
    model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: Some(start),
        end: Some(end),
    })
}

fn line(a: Point3, b: Point3) -> Curve {
    let dir = b.minus(&a);
    Curve::Line(Line { origin: a, dir })
}

fn plane(origin: Point3, normal: [i64; 3], reference: [i64; 3]) -> Surface {
    let vector = |v: [i64; 3]| {
        let [x, y, z] = v.map(|q| Real::from_rational(Rational::from_i64(q)));
        Vector3::new(x, y, z)
    };
    Surface::Plane(Plane {
        origin,
        normal: vector(normal),
        ref_dir: vector(reference),
    })
}

fn face(model: &mut Model, surface: Surface, uses: &[(EdgeId, Sense)], label: &str) -> FaceId {
    let surface = model.geom_mut().add_surface(surface);
    let coedges = uses
        .iter()
        .map(|&(edge, sense)| {
            record(
                model,
                GeomKey::Curve(model.edge(edge).curve),
                GeomKey::Surface(surface),
            );
            model.add_coedge(Coedge { edge, sense })
        })
        .collect();
    let loop_ = model.add_loop(Loop::Coedges(coedges));
    let origin = Origin {
        step: StepId::new("fixture").unwrap(),
        role: Role::Side {
            segment: ElementId::new(label).unwrap(),
            path: None,
        },
    };
    model.add_named_face(
        Face {
            surface,
            sense: Sense::Same,
            loops: vec![loop_],
        },
        vec![origin],
    )
}

pub fn solid(kind: Roof) -> Fixture {
    let mut dec = Decider::default();
    let mut model = Model::new();
    let joined = matches!(kind, Roof::JoinedDegreeElevated | Roof::JoinedNonlinear);
    let bump = matches!(kind, Roof::NonuniformBump);
    let xs = if joined {
        vec![r("0"), r("1/2"), r("1")]
    } else {
        vec![r("0"), r("1")]
    };
    let n = xs.len() - 1;
    let height = |x: &Rational, y: &Rational| {
        if bump {
            r("1") + x * &(&r("1") - x)
        } else {
            r("1") + x * y
        }
    };
    let point = |i: usize, layer: usize, y: usize| {
        let y = Rational::from_i64(y as i64);
        let z = if layer == 0 {
            r("0")
        } else {
            height(&xs[i], &y)
        };
        Point3::new(xs[i].clone().into(), y.into(), z.into())
    };
    let vertices: Vec<[[VertexId; 2]; 2]> = (0..xs.len())
        .map(|i| {
            std::array::from_fn(|layer| {
                std::array::from_fn(|y| {
                    let point = model
                        .geom_mut()
                        .add_point(PointDef::Explicit(point(i, layer, y)));
                    model.add_vertex(Vertex { point })
                })
            })
        })
        .collect();
    let patches: Vec<_> = (0..n)
        .map(|i| {
            let elevated = matches!(kind, Roof::JoinedDegreeElevated) && i == 1;
            let weighted = matches!(kind, Roof::Rational | Roof::JoinedDegreeElevated)
                || matches!(kind, Roof::JoinedNonlinear) && i == 1;
            let ys = if elevated {
                vec![r("0"), r("2/3"), r("1")]
            } else {
                vec![r("0"), r("1")]
            };
            let weights = if elevated {
                vec![r("1"), r("3/2"), r("2")]
            } else if weighted {
                vec![r("1"), r("2")]
            } else {
                vec![r("1"); 2]
            };
            let cx = if bump {
                vec![r("0"), r("1/6"), r("2/3"), r("1")]
            } else {
                vec![xs[i].clone(), xs[i + 1].clone()]
            };
            let cz = [r("1"), r("7/6"), r("4/3"), r("1")];
            let mut controls = Vec::new();
            let mut all_weights = Vec::new();
            for (row, x) in cx.iter().enumerate() {
                for (y, weight) in ys.iter().zip(&weights) {
                    let z = if bump { cz[row].clone() } else { height(x, y) };
                    controls.push(Point3::new(x.clone().into(), y.clone().into(), z.into()));
                    all_weights.push(Real::from_rational(weight.clone()));
                }
            }
            let u_knots = if bump {
                ["2/3", "2/3", "2/3", "1", "5/3", "5/3", "5/3"]
                    .map(r)
                    .to_vec()
            } else {
                ["2/3", "2/3", "5/3", "5/3"].map(r).to_vec()
            };
            let v_knots = if elevated {
                ["-2", "-2", "-2", "7", "7", "7"].map(r).to_vec()
            } else {
                ["-2", "-2", "7", "7"].map(r).to_vec()
            };
            RationalPatch::new(
                BSplineSurface::new(
                    [if bump { 2 } else { 1 }, if elevated { 2 } else { 1 }],
                    [u_knots, v_knots],
                    [cx.len(), ys.len()],
                    controls,
                    all_weights,
                    &mut dec,
                )
                .unwrap(),
                &mut dec,
            )
            .unwrap_or_else(|error| panic!("{kind:?} patch {i}: {error:?}"))
        })
        .collect();
    let x_edges: Vec<[[EdgeId; 2]; 2]> = (0..n)
        .map(|i| {
            std::array::from_fn(|layer| {
                std::array::from_fn(|y| {
                    let curve = if layer == 1 {
                        Curve::BSpline(patches[i].boundary(if y == 0 {
                            SplineBoundary::VMin
                        } else {
                            SplineBoundary::VMax
                        }))
                    } else {
                        line(point(i, layer, y), point(i + 1, layer, y))
                    };
                    edge(
                        &mut model,
                        curve,
                        vertices[i][layer][y],
                        vertices[i + 1][layer][y],
                    )
                })
            })
        })
        .collect();
    let y_edges: Vec<[EdgeId; 2]> = (0..xs.len())
        .map(|i| {
            std::array::from_fn(|layer| {
                let curve = if layer == 1 && !(matches!(kind, Roof::JoinedNonlinear) && i == 1) {
                    let (patch, side) = if i == 0 {
                        (&patches[0], SplineBoundary::UMin)
                    } else {
                        (&patches[i - 1], SplineBoundary::UMax)
                    };
                    Curve::BSpline(patch.boundary(side))
                } else {
                    line(point(i, layer, 0), point(i, layer, 1))
                };
                edge(
                    &mut model,
                    curve,
                    vertices[i][layer][0],
                    vertices[i][layer][1],
                )
            })
        })
        .collect();
    let z_edges: Vec<[EdgeId; 2]> = (0..xs.len())
        .map(|i| {
            std::array::from_fn(|y| {
                edge(
                    &mut model,
                    line(point(i, 0, y), point(i, 1, y)),
                    vertices[i][0][y],
                    vertices[i][1][y],
                )
            })
        })
        .collect();
    let same = Sense::Same;
    let rev = Sense::Reversed;
    let mut faces = Vec::new();
    let mut roofs = Vec::new();
    for i in 0..n {
        faces.push(face(
            &mut model,
            plane(point(i, 0, 0), [0, 0, -1], [1, 0, 0]),
            &[
                (y_edges[i][0], same),
                (x_edges[i][0][1], same),
                (y_edges[i + 1][0], rev),
                (x_edges[i][0][0], rev),
            ],
            "floor",
        ));
        faces.push(face(
            &mut model,
            plane(point(i, 0, 0), [0, -1, 0], [1, 0, 0]),
            &[
                (x_edges[i][0][0], same),
                (z_edges[i + 1][0], same),
                (x_edges[i][1][0], rev),
                (z_edges[i][0], rev),
            ],
            "front",
        ));
        faces.push(face(
            &mut model,
            plane(point(i, 0, 1), [0, 1, 0], [1, 0, 0]),
            &[
                (z_edges[i][1], same),
                (x_edges[i][1][1], same),
                (z_edges[i + 1][1], rev),
                (x_edges[i][0][1], rev),
            ],
            "back",
        ));
        let roof = face(
            &mut model,
            Surface::BSpline(patches[i].clone().into()),
            &[
                (x_edges[i][1][0], same),
                (y_edges[i + 1][1], same),
                (x_edges[i][1][1], rev),
                (y_edges[i][1], rev),
            ],
            &format!("roof-{i}"),
        );
        faces.push(roof);
        roofs.push(roof);
    }
    faces.push(face(
        &mut model,
        plane(point(0, 0, 0), [-1, 0, 0], [0, 1, 0]),
        &[
            (z_edges[0][0], same),
            (y_edges[0][1], same),
            (z_edges[0][1], rev),
            (y_edges[0][0], rev),
        ],
        "left",
    ));
    faces.push(face(
        &mut model,
        plane(point(n, 0, 0), [1, 0, 0], [0, 1, 0]),
        &[
            (y_edges[n][0], same),
            (z_edges[n][1], same),
            (y_edges[n][1], rev),
            (z_edges[n][0], rev),
        ],
        "right",
    ));
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::default(),
    });
    Fixture {
        model,
        body,
        roofs,
        patches,
    }
}
