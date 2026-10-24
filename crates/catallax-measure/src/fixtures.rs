//! Hand-built contract cases, using only the public geometry and topology API.

use catallax_numbers::{Fact, FactKind, FactOrigin, MarkLog, Rational};
use catallax_shapes::{
    Circle, Curve, Cylinder, Ellipse, GeomKey, Line, Plane, PointDef, Pose, RPoint3, RVector3,
    Sphere, Surface,
};
use catallax_topology::{
    Body, BodyId, BodyKind, Coedge, Edge, Face, FaceId, Loop, Lump, Model, Sense, Shell, Vertex,
};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(r(x), r(y), r(z))
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(r(x), r(y), r(z))
}
fn plane(origin: RPoint3, normal: RVector3, reference: RVector3) -> Surface {
    Surface::Plane(
        Plane::new(
            (origin).into(),
            (normal).into(),
            (reference).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
fn circle(center: RPoint3, normal: RVector3, reference: RVector3, radius2: i64) -> Curve {
    Curve::Circle(
        Circle::new(
            (center).into(),
            (normal).into(),
            (r(radius2)).into(),
            (reference).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
fn line(a: &RPoint3, b: &RPoint3) -> Curve {
    Curve::Line(
        Line::new(
            (a.clone()).into(),
            (b.minus(a)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
}
fn fact(model: &mut Model, a: GeomKey, b: GeomKey) {
    model.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a,
        b,
        at: None,
        origin: FactOrigin::Construction("measure.fixture"),
    });
}
fn finish(model: &mut Model, faces: Vec<FaceId>) -> BodyId {
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    })
}

type FaceLoops = (usize, Vec<Vec<(usize, Sense)>>);
fn sketch(
    points: Vec<RPoint3>,
    surfaces: Vec<Surface>,
    curves: Vec<(Curve, usize, usize)>,
    faces: Vec<FaceLoops>,
    pose: &Pose,
) -> (Model, BodyId) {
    let mut model = Model::new();
    let vertices = points
        .into_iter()
        .map(|point| {
            let point = model
                .geom_mut()
                .add_point(PointDef::Explicit(pose.apply_point(&point.to_real())));
            model.add_vertex(Vertex { point })
        })
        .collect::<Vec<_>>();
    let surfaces = surfaces
        .into_iter()
        .map(|surface| model.geom_mut().add_surface(surface.transformed(pose)))
        .collect::<Vec<_>>();
    let edges = curves
        .into_iter()
        .map(|(curve, from, to)| {
            let curve = model.geom_mut().add_curve(
                curve
                    .transformed(pose)
                    .expect("fixture curves are analytic"),
            );
            for vertex in [vertices[from], vertices[to]] {
                let point = model.vertex(vertex).point;
                fact(&mut model, GeomKey::Point(point), GeomKey::Curve(curve));
            }
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(vertices[from]),
                end: Some(vertices[to]),
            })
        })
        .collect::<Vec<_>>();
    let faces = faces
        .into_iter()
        .map(|(surface, loops)| {
            let surface = surfaces[surface];
            let loops = loops
                .into_iter()
                .map(|coedges| {
                    let coedges = coedges
                        .into_iter()
                        .map(|(edge, sense)| {
                            let edge = edges[edge];
                            let curve = model.edge(edge).curve;
                            fact(&mut model, GeomKey::Curve(curve), GeomKey::Surface(surface));
                            model.add_coedge(Coedge { edge, sense })
                        })
                        .collect();
                    model.add_loop(Loop::Coedges(coedges))
                })
                .collect();
            model.add_face(Face {
                surface,
                sense: Sense::Same,
                loops,
            })
        })
        .collect();
    let body = finish(&mut model, faces);
    assert!(model.check(body).is_empty(), "{:?}", model.check(body));
    (model, body)
}

pub(super) fn octant(pose: &Pose) -> (Model, BodyId) {
    use Sense::{Reversed as R, Same as S};
    let points = vec![p(0, 0, 0), p(3, 0, 0), p(0, 3, 0), p(0, 0, 3)];
    let curves = vec![
        (line(&points[0], &points[1]), 0, 1),
        (line(&points[0], &points[2]), 0, 2),
        (line(&points[0], &points[3]), 0, 3),
        (circle(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0), 9), 1, 2),
        (circle(p(0, 0, 0), v(1, 0, 0), v(0, 1, 0), 9), 2, 3),
        (circle(p(0, 0, 0), v(0, 1, 0), v(0, 0, 1), 9), 3, 1),
    ];
    let surfaces = vec![
        plane(p(0, 0, 0), v(0, 0, -1), v(1, 0, 0)),
        plane(p(0, 0, 0), v(-1, 0, 0), v(0, 1, 0)),
        plane(p(0, 0, 0), v(0, -1, 0), v(0, 0, 1)),
        Surface::Sphere(
            Sphere::new(
                p(0, 0, 0).into(),
                r(3).into(),
                v(0, 0, 1).into(),
                v(1, 0, 0).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
    ];
    sketch(
        points,
        surfaces,
        curves,
        vec![
            (0, vec![vec![(1, S), (3, R), (0, R)]]),
            (1, vec![vec![(2, S), (4, R), (1, R)]]),
            (2, vec![vec![(0, S), (5, R), (2, R)]]),
            (3, vec![vec![(3, S), (4, S), (5, S)]]),
        ],
        pose,
    )
}

pub(super) fn half_cylinder(pose: &Pose) -> (Model, BodyId) {
    use Sense::{Reversed as R, Same as S};
    let points = vec![p(2, 0, 0), p(-2, 0, 0), p(2, 0, 3), p(-2, 0, 3)];
    let curves = vec![
        (circle(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0), 4), 0, 1),
        (circle(p(0, 0, 3), v(0, 0, 1), v(1, 0, 0), 4), 2, 3),
        (line(&points[0], &points[2]), 0, 2),
        (line(&points[1], &points[3]), 1, 3),
        (line(&points[1], &points[0]), 1, 0),
        (line(&points[3], &points[2]), 3, 2),
    ];
    let surfaces = vec![
        plane(p(0, 0, 0), v(0, 0, -1), v(1, 0, 0)),
        plane(p(0, 0, 3), v(0, 0, 1), v(1, 0, 0)),
        plane(p(0, 0, 0), v(0, -1, 0), v(1, 0, 0)),
        Surface::Cylinder(
            Cylinder::new(
                p(0, 0, 0).into(),
                v(0, 0, 1).into(),
                r(2).into(),
                v(1, 0, 0).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
    ];
    sketch(
        points,
        surfaces,
        curves,
        vec![
            (0, vec![vec![(0, R), (4, R)]]),
            (1, vec![vec![(1, S), (5, S)]]),
            (2, vec![vec![(4, S), (2, S), (5, R), (3, R)]]),
            (3, vec![vec![(0, S), (3, S), (1, R), (2, R)]]),
        ],
        pose,
    )
}

pub(super) fn steinmetz(pose: &Pose) -> (Model, BodyId) {
    use Sense::{Reversed as R, Same as S};
    let minus = Curve::Ellipse(
        Ellipse::new(
            (p(0, 0, 0)).into(),
            (v(-1, 0, 1)).into(),
            (v(1, 0, 1)).into(),
            (r(2)).into(),
            (r(1)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let plus = Curve::Ellipse(
        Ellipse::new(
            (p(0, 0, 0)).into(),
            (v(1, 0, 1)).into(),
            (v(-1, 0, 1)).into(),
            (r(2)).into(),
            (r(1)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let surfaces = vec![
        Surface::Cylinder(
            Cylinder::new(
                p(0, 0, 0).into(),
                v(0, 0, 1).into(),
                r(1).into(),
                v(1, 0, 0).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
        Surface::Cylinder(
            Cylinder::new(
                p(0, 0, 0).into(),
                v(1, 0, 0).into(),
                r(1).into(),
                v(0, 0, -1).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        ),
    ];
    sketch(
        vec![p(0, 1, 0), p(0, -1, 0)],
        surfaces,
        vec![
            (minus.clone(), 0, 1),
            (minus, 1, 0),
            (plus.clone(), 0, 1),
            (plus, 1, 0),
        ],
        vec![
            (0, vec![vec![(3, S), (1, R)]]),
            (0, vec![vec![(0, S), (2, R)]]),
            (1, vec![vec![(0, R), (3, R)]]),
            (1, vec![vec![(2, S), (1, S)]]),
        ],
        pose,
    )
}

pub(super) fn void_block() -> (Model, BodyId) {
    let mut model = Model::new();
    let outer =
        catallax_topology::block(&mut model, [r(10), r(10), r(10)], &Pose::identity()).unwrap();
    let inner = catallax_topology::block(
        &mut model,
        [r(2), r(2), r(2)],
        &Pose::new(
            catallax_shapes::Rotation::identity(),
            (v(4, 4, 4)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let mut inner_faces = Vec::new();
    for id in model.faces_of_body(inner) {
        let face = model.face(id).clone();
        let mut loops = Vec::new();
        for id in face.loops {
            let Loop::Coedges(coedges) = model.loop_(id).clone() else {
                unreachable!()
            };
            let coedges = coedges
                .iter()
                .rev()
                .map(|id| {
                    let coedge = model.coedge(*id);
                    model.add_coedge(Coedge {
                        edge: coedge.edge,
                        sense: if coedge.sense == Sense::Same {
                            Sense::Reversed
                        } else {
                            Sense::Same
                        },
                    })
                })
                .collect();
            loops.push(model.add_loop(Loop::Coedges(coedges)));
        }
        inner_faces.push(model.add_face(Face {
            surface: face.surface,
            sense: Sense::Reversed,
            loops,
        }));
    }
    let inner_shell = model.add_shell(Shell { faces: inner_faces });
    let outer_shell = model.lump(model.body(outer).lumps[0]).shells[0];
    let lump = model.add_lump(Lump {
        shells: vec![outer_shell, inner_shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    assert!(model.check(body).is_empty());
    (model, body)
}

pub(super) fn through_hole(scale: &Rational, pose: &Pose) -> (Model, BodyId) {
    let mut model = Model::new();
    let local_point = |x, y, z| RPoint3::new(scale * &r(x), scale * &r(y), scale * &r(z));
    let block =
        catallax_topology::block(&mut model, std::array::from_fn(|_| scale * &r(4)), pose).unwrap();
    let mut faces = model.faces_of_body(block);
    let wall = Surface::Cylinder(
        Cylinder::new(
            (local_point(2, 2, 0)).into(),
            (v(0, 0, 1)).into(),
            (scale.clone()).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
    .transformed(pose);
    let wall = model.geom_mut().add_surface(wall);
    let mut wall_loops = Vec::new();
    for (index, z, hole_sense) in [(4, 0, Sense::Same), (5, 4, Sense::Reversed)] {
        let old_face = model.face(faces[index]).clone();
        let curve = Curve::Circle(
            Circle::new(
                (local_point(2, 2, z)).into(),
                (v(0, 0, 1)).into(),
                (scale * scale).into(),
                (v(1, 0, 0)).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        )
        .transformed(pose)
        .unwrap();
        let curve = model.geom_mut().add_curve(curve);
        fact(
            &mut model,
            GeomKey::Curve(curve),
            GeomKey::Surface(old_face.surface),
        );
        fact(&mut model, GeomKey::Curve(curve), GeomKey::Surface(wall));
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        let hole_coedge = model.add_coedge(Coedge {
            edge,
            sense: hole_sense,
        });
        let hole_loop = model.add_loop(Loop::Coedges(vec![hole_coedge]));
        let wall_coedge = model.add_coedge(Coedge {
            edge,
            sense: if hole_sense == Sense::Same {
                Sense::Reversed
            } else {
                Sense::Same
            },
        });
        wall_loops.push(model.add_loop(Loop::Coedges(vec![wall_coedge])));
        let mut loops = Vec::new();
        for old_loop in old_face.loops {
            let Loop::Coedges(coedges) = model.loop_(old_loop).clone() else {
                unreachable!()
            };
            let coedges = coedges
                .into_iter()
                .map(|id| model.add_coedge(model.coedge(id).clone()))
                .collect();
            loops.push(model.add_loop(Loop::Coedges(coedges)));
        }
        loops.push(hole_loop);
        faces[index] = model.add_face(Face {
            surface: old_face.surface,
            sense: old_face.sense,
            loops,
        });
    }
    faces.push(model.add_face(Face {
        surface: wall,
        sense: Sense::Reversed,
        loops: wall_loops,
    }));
    let body = finish(&mut model, faces);
    assert!(model.check(body).is_empty(), "{:?}", model.check(body));
    (model, body)
}

pub(super) fn half_torus() -> (Model, BodyId) {
    half_torus_posed(&Pose::identity())
}

pub(super) fn half_torus_posed(pose: &Pose) -> (Model, BodyId) {
    let mut model = Model::new();
    let torus = model.geom_mut().add_surface(
        Surface::Torus(
            catallax_shapes::Torus::new(
                p(0, 0, 0).into(),
                v(0, 0, 1).into(),
                r(3).into(),
                r(1).into(),
                v(1, 0, 0).into(),
                &mut catallax_numbers::Decider::default(),
            )
            .unwrap(),
        )
        .transformed(pose),
    );
    let mut rings = Vec::new();
    let mut caps = Vec::new();
    for (x, normal, reference) in [(3, v(0, -1, 0), v(1, 0, 0)), (-3, v(0, 1, 0), v(-1, 0, 0))] {
        let plane = model
            .geom_mut()
            .add_surface(plane(p(x, 0, 0), v(0, -1, 0), v(1, 0, 0)).transformed(pose));
        let curve = model.geom_mut().add_curve(
            circle(p(x, 0, 0), normal, reference, 1)
                .transformed(pose)
                .unwrap(),
        );
        fact(&mut model, GeomKey::Curve(curve), GeomKey::Surface(torus));
        fact(&mut model, GeomKey::Curve(curve), GeomKey::Surface(plane));
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        rings.push(edge);
        caps.push(plane);
    }
    let mut loops = Vec::new();
    let mut faces = Vec::new();
    for i in 0..2 {
        let coedge = model.add_coedge(Coedge {
            edge: rings[i],
            sense: if i == 0 { Sense::Reversed } else { Sense::Same },
        });
        loops.push(model.add_loop(Loop::Coedges(vec![coedge])));
        let coedge = model.add_coedge(Coedge {
            edge: rings[i],
            sense: if i == 0 { Sense::Same } else { Sense::Reversed },
        });
        let boundary = model.add_loop(Loop::Coedges(vec![coedge]));
        faces.push(model.add_face(Face {
            surface: caps[i],
            sense: Sense::Same,
            loops: vec![boundary],
        }));
    }
    faces.insert(
        0,
        model.add_face(Face {
            surface: torus,
            sense: Sense::Same,
            loops,
        }),
    );
    let body = finish(&mut model, faces);
    assert!(model.check(body).is_empty(), "{:?}", model.check(body));
    (model, body)
}

pub(super) fn cap(radius: i64, height: i64) -> (Model, BodyId) {
    let mut model = Model::new();
    let sphere = model.geom_mut().add_surface(Surface::Sphere(
        Sphere::new(
            (p(0, 0, 0)).into(),
            (r(radius)).into(),
            (v(0, 0, 1)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let plane = model
        .geom_mut()
        .add_surface(plane(p(0, 0, height), v(0, 0, -1), v(1, 0, 0)));
    let curve = model.geom_mut().add_curve(circle(
        p(0, 0, height),
        v(0, 0, 1),
        v(1, 0, 0),
        radius * radius - height * height,
    ));
    for surface in [sphere, plane] {
        fact(&mut model, GeomKey::Curve(curve), GeomKey::Surface(surface));
    }
    let edge = model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let mut faces = Vec::new();
    for (surface, sense) in [(sphere, Sense::Same), (plane, Sense::Reversed)] {
        let coedge = model.add_coedge(Coedge { edge, sense });
        let boundary = model.add_loop(Loop::Coedges(vec![coedge]));
        faces.push(model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![boundary],
        }));
    }
    let body = finish(&mut model, faces);
    assert!(model.check(body).is_empty());
    (model, body)
}

pub(super) fn coplanar_split_block() -> (Model, BodyId) {
    let mut model = Model::new();
    let block =
        catallax_topology::block(&mut model, [r(4), r(4), r(4)], &Pose::identity()).unwrap();
    let mut faces = model.faces_of_body(block);
    let original = model.face(faces[5]).clone();
    let curve = model
        .geom_mut()
        .add_curve(circle(p(2, 2, 4), v(0, 0, 1), v(1, 0, 0), 1));
    fact(
        &mut model,
        GeomKey::Curve(curve),
        GeomKey::Surface(original.surface),
    );
    let edge = model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let inner = model.add_coedge(Coedge {
        edge,
        sense: Sense::Reversed,
    });
    let inner = model.add_loop(Loop::Coedges(vec![inner]));
    let mut loops = Vec::new();
    for boundary in original.loops {
        let Loop::Coedges(coedges) = model.loop_(boundary).clone() else {
            unreachable!()
        };
        let coedges = coedges
            .into_iter()
            .map(|id| model.add_coedge(model.coedge(id).clone()))
            .collect();
        loops.push(model.add_loop(Loop::Coedges(coedges)));
    }
    loops.push(inner);
    faces[5] = model.add_face(Face {
        surface: original.surface,
        sense: Sense::Same,
        loops,
    });
    let disc = model.add_coedge(Coedge {
        edge,
        sense: Sense::Same,
    });
    let disc = model.add_loop(Loop::Coedges(vec![disc]));
    faces.push(model.add_face(Face {
        surface: original.surface,
        sense: Sense::Same,
        loops: vec![disc],
    }));
    let body = finish(&mut model, faces);
    assert!(model.check(body).is_empty(), "{:?}", model.check(body));
    (model, body)
}

pub(super) fn torus_rectangle(long: bool, pose: &Pose) -> (Model, BodyId) {
    let sy = if long { -1 } else { 1 };
    let torus = Surface::Torus(
        catallax_shapes::Torus::new(
            p(0, 0, 0).into(),
            v(0, 0, 1).into(),
            r(3).into(),
            r(1).into(),
            v(0, 1, 0).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    );
    let surfaces = vec![
        torus,
        plane(p(0, 0, 0), v(0, 0, -1), v(1, 0, 0)),
        plane(p(0, 0, 1), v(0, 0, 1), v(1, 0, 0)),
        plane(p(0, 0, 0), v(0, -1, 0), v(1, 0, 0)),
        plane(p(0, 0, 0), v(-sy, 0, 0), v(0, 0, 1)),
    ];
    let points = vec![
        p(4, 0, 0),
        p(0, sy * 4, 0),
        p(0, sy * 3, 1),
        p(3, 0, 1),
        p(0, 0, 0),
        p(0, 0, 1),
    ];
    let curves = vec![
        (circle(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0), 16), 0, 1),
        (circle(p(0, sy * 3, 0), v(sy, 0, 0), v(0, sy, 0), 1), 1, 2),
        (circle(p(0, 0, 1), v(0, 0, 1), v(1, 0, 0), 9), 3, 2),
        (circle(p(3, 0, 0), v(0, -1, 0), v(1, 0, 0), 1), 0, 3),
        (line(&points[4], &points[0]), 4, 0),
        (line(&points[4], &points[1]), 4, 1),
        (line(&points[5], &points[3]), 5, 3),
        (line(&points[5], &points[2]), 5, 2),
        (line(&points[4], &points[5]), 4, 5),
    ];
    let (s, n) = (Sense::Same, Sense::Reversed);
    let faces = vec![
        (0, vec![vec![(0, s), (1, s), (2, n), (3, n)]]),
        (1, vec![vec![(5, s), (0, n), (4, n)]]),
        (2, vec![vec![(6, s), (2, s), (7, n)]]),
        (3, vec![vec![(4, s), (3, s), (6, n), (8, n)]]),
        (4, vec![vec![(8, s), (7, s), (1, n), (5, n)]]),
    ];
    let (mut model, body) = sketch(points, surfaces, curves, faces, pose);
    let mut faces = model.faces_of_body(body);
    let source = catallax_recipe::Origin {
        step: catallax_recipe::StepId::new("source").unwrap(),
        role: catallax_recipe::Role::End,
    };
    faces[0] = model.add_named_face(
        model.face(faces[0]).clone(),
        vec![catallax_recipe::Origin {
            step: catallax_recipe::StepId::new("offset").unwrap(),
            role: catallax_recipe::Role::Offset {
                of: Box::new(source),
            },
        }],
    );
    let body = finish(&mut model, faces);
    assert!(model.check(body).is_empty(), "{:?}", model.check(body));
    (model, body)
}
