use std::collections::BTreeMap;

use catallax_boolean::{BoolOp, boolean};
use catallax_numbers::{Decider, MarkLog, Rational};
use catallax_shapes::{Circle, Curve, Pose, RPoint3, RVector3, Rotation, Surface};
use catallax_topology::{
    self as topology, Body, BodyId, BodyKind, Coedge, Edge, Face, Loop, Lump, Model, Sense, Shell,
};

use crate::{Mesh, MeshOptions, mesh_body};

fn options() -> MeshOptions {
    MeshOptions {
        chord_tolerance: Rational::parse("0.15").unwrap(),
        ..MeshOptions::default()
    }
}
fn pose(x: i64, y: i64, z: i64) -> Pose {
    Pose::new(
        Rotation::identity(),
        (RVector3::new(x.into(), y.into(), z.into())).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap()
}
fn mesh(model: &Model, body: BodyId) -> Mesh {
    mesh_body(model, body, &options(), &mut Decider::default())
        .unwrap()
        .value
}

fn assert_closed(mesh: &Mesh) {
    let mut edges = BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();
    let mut volume = 0.0;
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices.map(|v| mesh.positions[v]);
        let ab: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
        let ac: [f64; 3] = std::array::from_fn(|i| c[i] - a[i]);
        let cross = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        assert!(
            cross.iter().any(|x| x.abs() > 1e-20),
            "zero-area triangle {triangle:?}"
        );
        volume += a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            edges.entry((a.min(b), a.max(b))).or_default().push((a, b));
        }
    }
    for (edge, uses) in edges {
        assert_eq!(uses.len(), 2, "edge {edge:?}: {uses:?}");
        assert_eq!(uses[0], (uses[1].1, uses[1].0), "edge winding {edge:?}");
    }
    assert!(
        volume > 0.0,
        "outward winding must enclose positive material volume"
    );
    for edge in &mesh.edges {
        for pair in edge.vertices.windows(2) {
            let count = mesh
                .triangles
                .iter()
                .filter(|t| {
                    (0..3).any(|i| {
                        let (a, b) = (t.vertices[i], t.vertices[(i + 1) % 3]);
                        (a == pair[0] && b == pair[1]) || (b == pair[0] && a == pair[1])
                    })
                })
                .count();
            assert_eq!(count, 2, "shared boundary {:?} {pair:?}", edge.edge);
        }
    }
}

fn surface_distance(surface: &Surface, p: [f64; 3]) -> f64 {
    let (origin, axis) = match surface {
        Surface::Plane(s) => (&s.origin, &s.normal),
        Surface::Cylinder(s) => (&s.origin, &s.axis),
        Surface::Cone(s) => (&s.apex, &s.axis),
        Surface::Sphere(s) => (&s.center, &s.axis),
        Surface::Torus(s) => (&s.center, &s.axis),
        _ => unreachable!(),
    };
    let d: [f64; 3] = std::array::from_fn(|i| p[i] - [&origin.x, &origin.y, &origin.z][i].to_f64());
    let a = [&axis.x, &axis.y, &axis.z].map(|r| r.to_f64());
    let norm = libm::sqrt(a.iter().map(|x| x * x).sum());
    let z = (0..3).map(|i| d[i] * a[i] / norm).sum::<f64>();
    let radial = libm::sqrt((d.iter().map(|x| x * x).sum::<f64>() - z * z).max(0.0));
    match surface {
        Surface::Plane(_) => z.abs(),
        Surface::Cylinder(s) => (radial - s.radius().to_f64()).abs(),
        Surface::Cone(s) => {
            let k = s.slope.to_f64();
            (radial - k * z).abs() / libm::sqrt(1.0 + k * k)
        }
        Surface::Sphere(s) => (libm::sqrt(d.iter().map(|x| x * x).sum()) - s.radius.to_f64()).abs(),
        Surface::Torus(s) => {
            let x = radial - s.major.to_f64();
            (libm::sqrt(x * x + z * z) - s.minor.to_f64()).abs()
        }
        _ => unreachable!(),
    }
}
fn assert_accuracy(model: &Model, mesh: &Mesh) {
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options().chord_tolerance);
        let surface = model.geom().surface(model.face(face.face).surface);
        for t in &mesh.triangles[face.triangles.clone()] {
            assert_eq!(t.face, face.face);
            for weights in [
                [0.5, 0.5, 0.0],
                [0.5, 0.0, 0.5],
                [0.0, 0.5, 0.5],
                [1.0 / 3.0; 3],
                [0.2, 0.3, 0.5],
            ] {
                let p = std::array::from_fn(|i| {
                    (0..3)
                        .map(|j| weights[j] * mesh.positions[t.vertices[j]][i])
                        .sum()
                });
                let measured = surface_distance(surface, p);
                assert!(
                    measured <= face.max_chord_error.to_f64_nearest() + 1e-12,
                    "face {:?}: measured {measured}, bound {}",
                    face.face,
                    face.max_chord_error
                );
            }
        }
    }
    for edge in &mesh.edges {
        assert!(edge.max_chord_error <= options().chord_tolerance);
    }
}

#[test]
fn primitives_close_at_seams_poles_and_apices_with_bounded_chords() {
    let mut model = Model::new();
    let identity = Pose::identity();
    let rotation = Rotation::from_quaternion(1.into(), 2.into(), 3.into(), 4.into()).unwrap();
    let transformed = Pose::new(
        rotation,
        (RVector3::new(7.into(), (-2).into(), 3.into())).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    let block = topology::block(&mut model, [2.into(), 3.into(), 5.into()], &transformed).unwrap();
    let cylinder = topology::cylinder(&mut model, 2.into(), 5.into(), &transformed).unwrap();
    let sphere = topology::sphere(&mut model, 2.into(), &identity).unwrap();
    let cone = topology::cone(&mut model, 3.into(), 0.into(), 4.into(), &identity).unwrap();
    let frustum = topology::cone(&mut model, 1.into(), 3.into(), 4.into(), &transformed).unwrap();
    let torus = topology::torus(&mut model, 3.into(), 1.into(), &identity).unwrap();
    for (name, body) in [
        ("block", block),
        ("cylinder", cylinder),
        ("sphere", sphere),
        ("cone", cone),
        ("frustum", frustum),
        ("torus", torus),
    ] {
        let result = mesh_body(&model, body, &options(), &mut Decider::default())
            .unwrap_or_else(|e| panic!("{name}: {e}"))
            .value;
        assert_eq!(result.faces.len(), model.faces_of_body(body).len());
        assert_closed(&result);
        assert_accuracy(&model, &result);
    }
}

fn lens(op: BoolOp) -> (Model, BodyId) {
    let mut model = Model::new();
    let a = topology::sphere(&mut model, 2.into(), &Pose::identity()).unwrap();
    let b = topology::sphere(&mut model, 2.into(), &pose(1, 0, 0)).unwrap();
    let body = boolean(&mut model, a, b, op, &mut Decider::default())
        .unwrap()
        .body;
    (model, body)
}
#[test]
fn spherical_caps_complements_and_reversed_faces_mesh_the_kept_side() {
    for op in [BoolOp::Intersect, BoolOp::Unite, BoolOp::Subtract] {
        let (model, body) = lens(op);
        let result = mesh(&model, body);
        assert_closed(&result);
        assert_accuracy(&model, &result);
        for t in &result.triangles {
            let p: [f64; 3] = std::array::from_fn(|i| {
                t.vertices
                    .iter()
                    .map(|&v| result.positions[v][i])
                    .sum::<f64>()
                    / 3.0
            });
            let a = p.iter().map(|x| x * x).sum::<f64>() <= 4.0 + 1e-12;
            let b = (p[0] - 1.0).powi(2) + p[1] * p[1] + p[2] * p[2] <= 4.0 + 1e-12;
            match op {
                BoolOp::Intersect => assert!(a && b),
                BoolOp::Unite => assert!(a || b),
                BoolOp::Subtract => {
                    assert!(a && (p[0] - 1.0).powi(2) + p[1] * p[1] + p[2] * p[2] >= 3.5)
                }
                _ => unreachable!(),
            }
        }
    }
}

#[test]
fn planar_holes_and_a_reversed_cylinder_share_every_boundary_segment() {
    let mut model = Model::new();
    let a = topology::block(
        &mut model,
        [6.into(), 6.into(), 4.into()],
        &Pose::identity(),
    )
    .unwrap();
    let b = topology::cylinder(&mut model, 1.into(), 6.into(), &pose(3, 3, -1)).unwrap();
    let body = boolean(&mut model, a, b, BoolOp::Subtract, &mut Decider::default())
        .unwrap()
        .body;
    let result = mesh(&model, body);
    assert_closed(&result);
    assert_accuracy(&model, &result);
    for t in &result.triangles {
        let p: [f64; 3] = std::array::from_fn(|i| {
            t.vertices
                .iter()
                .map(|&v| result.positions[v][i])
                .sum::<f64>()
                / 3.0
        });
        assert!(
            (p[0] - 3.0).powi(2) + (p[1] - 3.0).powi(2) > 0.9,
            "triangle fills drilled hole"
        );
    }
}

#[test]
fn traced_sphere_cylinder_intersection_is_shared_and_within_tolerance() {
    let mut model = Model::new();
    let a = topology::cylinder(&mut model, 1.into(), 6.into(), &pose(0, 0, -3)).unwrap();
    let b = topology::sphere(&mut model, 2.into(), &pose(2, 0, 0)).unwrap();
    let body = boolean(&mut model, a, b, BoolOp::Intersect, &mut Decider::default())
        .unwrap()
        .body;
    assert!(model.edges_of_body(body).iter().any(|e| matches!(
        model.geom().curve(model.edge(*e).curve),
        Curve::Intersection(_)
    )));
    let result = mesh(&model, body);
    assert_closed(&result);
    assert_accuracy(&model, &result);
}

#[test]
fn split_torus_closes_both_periodic_dimensions_across_shared_rings() {
    let mut model = Model::new();
    let original = topology::torus(&mut model, 3.into(), 1.into(), &Pose::identity()).unwrap();
    let surface = model.face(model.faces_of_body(original)[0]).surface;
    let edges: Vec<_> = [4, 2]
        .into_iter()
        .map(|r| {
            let curve = model.geom_mut().add_curve(Curve::Circle(
                Circle::new(
                    (RPoint3::origin()).into(),
                    (RVector3::new(0.into(), 0.into(), 1.into())).into(),
                    (r * r).into(),
                    (RVector3::new(1.into(), 0.into(), 0.into())).into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ));
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            })
        })
        .collect();
    let mut faces = Vec::new();
    for senses in [
        [Sense::Same, Sense::Reversed],
        [Sense::Reversed, Sense::Same],
    ] {
        let loops = edges
            .iter()
            .zip(senses)
            .map(|(&edge, sense)| {
                let c = model.add_coedge(Coedge { edge, sense });
                model.add_loop(Loop::Coedges(vec![c]))
            })
            .collect();
        faces.push(model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops,
        }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    let result = mesh(&model, body);
    assert_closed(&result);
    assert_accuracy(&model, &result);
    for face in &result.faces {
        for triangle in &result.triangles[face.triangles.clone()] {
            let z = triangle
                .vertices
                .iter()
                .map(|&v| result.positions[v][2])
                .sum::<f64>();
            if face.face == result.faces[0].face {
                assert!(z >= -1e-12);
            } else {
                assert!(z <= 1e-12);
            }
        }
    }
}

#[test]
fn obj_bytes_are_identical_with_one_two_and_four_concurrent_workers() {
    let (model, body) = lens(BoolOp::Subtract);
    let mut expected = Vec::new();
    mesh(&model, body).write_obj(&mut expected).unwrap();
    for threads in [1, 2, 4] {
        std::thread::scope(|scope| {
            let jobs: Vec<_> = (0..threads)
                .map(|_| {
                    scope.spawn(|| {
                        let mut bytes = Vec::new();
                        mesh(&model, body).write_obj(&mut bytes).unwrap();
                        bytes
                    })
                })
                .collect();
            for job in jobs {
                assert_eq!(job.join().unwrap(), expected);
            }
        });
    }
}

#[test]
fn empty_bodies_and_unattainable_tolerances_do_not_emit_false_proofs() {
    let mut model = Model::new();
    let empty = topology::empty(&mut model);
    assert_eq!(mesh(&model, empty), Mesh::default());
    let body = topology::sphere(&mut model, 1.into(), &Pose::identity()).unwrap();
    let limited = MeshOptions {
        max_vertices: 10,
        ..options()
    };
    assert!(mesh_body(&model, body, &limited, &mut Decider::default()).is_err());
    let invalid = MeshOptions {
        chord_tolerance: Rational::zero(),
        ..options()
    };
    assert_eq!(
        mesh_body(&model, body, &invalid, &mut Decider::default())
            .unwrap_err()
            .code,
        "mesh.invalid-options"
    );
}

#[test]
fn a_sphere_boundary_through_definition_poles_uses_a_regular_chart() {
    let mut model = Model::new();
    let original = topology::sphere(&mut model, 2.into(), &Pose::identity()).unwrap();
    let surface = model.face(model.faces_of_body(original)[0]).surface;
    let curve = model.geom_mut().add_curve(Curve::Circle(
        Circle::new(
            (RPoint3::origin()).into(),
            (RVector3::new(1.into(), 0.into(), 0.into())).into(),
            4.into(),
            (RVector3::new(0.into(), 1.into(), 0.into())).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let edge = model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let mut faces = Vec::new();
    for sense in [Sense::Same, Sense::Reversed] {
        let c = model.add_coedge(Coedge { edge, sense });
        let l = model.add_loop(Loop::Coedges(vec![c]));
        faces.push(model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![l],
        }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    let result = mesh(&model, body);
    assert_closed(&result);
    assert_accuracy(&model, &result);
    for face in &result.faces {
        for triangle in &result.triangles[face.triangles.clone()] {
            let x = triangle
                .vertices
                .iter()
                .map(|&v| result.positions[v][0])
                .sum::<f64>();
            if face.face == result.faces[0].face {
                assert!(x >= -1e-12);
            } else {
                assert!(x <= 1e-12);
            }
        }
    }
}

#[test]
fn a_sub_ulp_gap_between_a_planar_hole_and_its_outer_boundary_stays_connected() {
    let mut model = Model::new();
    let plate = topology::block(
        &mut model,
        [2.into(), 3.into(), 1.into()],
        &Pose::identity(),
    )
    .unwrap();
    let epsilon = Rational::parse("1e-36").unwrap();
    let flange_pose = Pose::new(
        Rotation::identity(),
        RVector3::new(Rational::one() - epsilon, 1.into(), 0.into()).into(),
        &mut Decider::default(),
    )
    .unwrap();
    let flange = topology::block(&mut model, [1.into(), 1.into(), 2.into()], &flange_pose).unwrap();
    let body = boolean(
        &mut model,
        plate,
        flange,
        BoolOp::Unite,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let options = MeshOptions::default();
    let mesh = mesh_body(&model, body, &options, &mut Decider::default())
        .unwrap()
        .value;
    assert_eq!(mesh.faces.len(), model.faces_of_body(body).len());
    assert!(
        mesh.faces
            .iter()
            .all(|f| f.max_chord_error <= options.chord_tolerance)
    );
    // Binary64 collapses the thin strip geometrically, but not its incidence.
    let mut incidence = BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let a = triangle.vertices[i];
            let b = triangle.vertices[(i + 1) % 3];
            incidence
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((a, b));
        }
    }
    for uses in incidence.values() {
        assert_eq!(uses.len(), 2);
        assert_eq!(uses[0], (uses[1].1, uses[1].0));
    }
    for edge in &mesh.edges {
        for pair in edge.vertices.windows(2) {
            assert!(incidence.contains_key(&(pair[0].min(pair[1]), pair[0].max(pair[1]))));
        }
    }
}

#[test]
fn planar_hole_bridges_preserve_multiple_bores() {
    let mut model = Model::new();
    let mut body = topology::block(
        &mut model,
        [8.into(), 4.into(), 2.into()],
        &Pose::identity(),
    )
    .unwrap();
    for x in [2, 6] {
        let bore = topology::cylinder(&mut model, 1.into(), 4.into(), &pose(x, 2, -1)).unwrap();
        body = boolean(
            &mut model,
            body,
            bore,
            BoolOp::Subtract,
            &mut Decider::default(),
        )
        .unwrap()
        .body;
    }
    let result = mesh(&model, body);
    assert_closed(&result);
    assert_accuracy(&model, &result);
    for face in &result.faces {
        if !matches!(
            model.geom().surface(model.face(face.face).surface),
            Surface::Plane(_)
        ) {
            continue;
        }
        for triangle in &result.triangles[face.triangles.clone()] {
            let center: [f64; 3] = std::array::from_fn(|i| {
                triangle
                    .vertices
                    .iter()
                    .map(|&v| result.positions[v][i])
                    .sum::<f64>()
                    / 3.0
            });
            for x in [2.0, 6.0] {
                assert!((center[0] - x).powi(2) + (center[1] - 2.0).powi(2) > 0.9);
            }
        }
    }
}

#[test]
fn requested_controls_bound_the_returned_triangles_and_polylines() {
    let mut model = Model::new();
    let pose = Pose::identity();
    let bodies = [
        topology::block(&mut model, [2.into(), 3.into(), 5.into()], &pose).unwrap(),
        topology::cylinder(&mut model, 2.into(), 3.into(), &pose).unwrap(),
        topology::cone(&mut model, 2.into(), 0.into(), 3.into(), &pose).unwrap(),
        topology::sphere(&mut model, 2.into(), &pose).unwrap(),
        topology::torus(&mut model, 3.into(), 1.into(), &pose).unwrap(),
    ];
    let options = MeshOptions {
        angular_deviation: Some(Rational::parse("0.8").unwrap()),
        max_edge_length: Some(Rational::parse("1.5").unwrap()),
        max_aspect_ratio: Some(Rational::one()),
        ..options()
    };
    for body in bodies {
        let result = mesh_body(&model, body, &options, &mut Decider::default()).unwrap();
        let mesh = result.value;
        assert_closed(&mesh);
        assert_accuracy(&model, &mesh);
        assert!(
            mesh.unavailable_controls
                .iter()
                .all(|c| c.control == crate::MeshControl::AspectRatio)
        );
        for edge in &mesh.edges {
            for pair in edge.vertices.windows(2) {
                let a = mesh.positions[pair[0]];
                let b = mesh.positions[pair[1]];
                let length = libm::sqrt((0..3).map(|i| (a[i] - b[i]).powi(2)).sum());
                assert!(length <= 1.5 + 1e-12);
            }
        }
        for face in &mesh.faces {
            assert!(
                face.max_angular_deviation.as_ref().unwrap()
                    <= options.angular_deviation.as_ref().unwrap()
            );
            assert!(
                mesh.unavailable_controls
                    .iter()
                    .any(|c| c.face == face.face)
            );
            let surface = model.geom().surface(model.face(face.face).surface);
            for triangle in &mesh.triangles[face.triangles.clone()] {
                let points = triangle.vertices.map(|v| mesh.positions[v]);
                let lengths: [f64; 3] = std::array::from_fn(|j| {
                    libm::sqrt(
                        (0..3)
                            .map(|i| (points[j][i] - points[(j + 1) % 3][i]).powi(2))
                            .sum(),
                    )
                });
                assert!(lengths.iter().all(|l| *l <= 1.5 + 1e-12));
                let ab: [f64; 3] = std::array::from_fn(|i| points[1][i] - points[0][i]);
                let ac: [f64; 3] = std::array::from_fn(|i| points[2][i] - points[0][i]);
                let normal = [
                    ab[1] * ac[2] - ab[2] * ac[1],
                    ab[2] * ac[0] - ab[0] * ac[2],
                    ab[0] * ac[1] - ab[1] * ac[0],
                ];
                let norm = libm::sqrt(normal.iter().map(|x| x * x).sum());
                let ratio = lengths.iter().copied().fold(0.0, f64::max).powi(2) / norm;
                assert!(ratio <= face.max_aspect_ratio.as_ref().unwrap().to_f64_nearest() + 1e-10);
                for p in points {
                    let exact = catallax_shapes::Point3::new(
                        crate::chart::q(p[0]).into(),
                        crate::chart::q(p[1]).into(),
                        crate::chart::q(p[2]).into(),
                    );
                    let gradient = surface.gradient(&exact);
                    let mut analytic = [&gradient.x, &gradient.y, &gradient.z].map(|x| x.to_f64());
                    if model.face(face.face).sense == Sense::Reversed {
                        analytic = analytic.map(|x| -x);
                    }
                    let anorm = libm::sqrt(analytic.iter().map(|x| x * x).sum());
                    if anorm == 0.0 {
                        continue;
                    }
                    let cosine =
                        (0..3).map(|i| normal[i] * analytic[i]).sum::<f64>() / (norm * anorm);
                    let angle = libm::acos(cosine.clamp(-1.0, 1.0));
                    assert!(angle <= 0.8 + 1e-10, "face {:?}: angle {angle}", face.face);
                }
            }
        }
    }
}

#[test]
fn collapsed_planar_facets_preserve_incidence_and_report_unavailable_quality() {
    let mut model = Model::new();
    let body = topology::block(
        &mut model,
        [Rational::parse("1e-36").unwrap(), 2.into(), 2.into()],
        &pose(1, 0, 0),
    )
    .unwrap();
    let options = MeshOptions {
        angular_deviation: Some(Rational::parse("0.8").unwrap()),
        max_edge_length: Some(Rational::one()),
        max_aspect_ratio: Some(Rational::from_i64(10)),
        ..options()
    };
    let mesh = mesh_body(&model, body, &options, &mut Decider::default())
        .unwrap()
        .value;
    let mut incidence = BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            incidence
                .entry((a.min(b), a.max(b)))
                .or_default()
                .push((a, b));
        }
    }
    for uses in incidence.values() {
        assert_eq!(uses.len(), 2);
        assert_eq!(uses[0], (uses[1].1, uses[1].0));
    }
    for edge in &mesh.edges {
        for pair in edge.vertices.windows(2) {
            assert!(incidence.contains_key(&(pair[0].min(pair[1]), pair[0].max(pair[1]))));
        }
    }
    assert!(
        mesh.faces
            .iter()
            .any(|face| face.max_aspect_ratio.is_none())
    );
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
        let mut collapsed = false;
        if face.max_aspect_ratio.is_none() {
            assert!(face.max_angular_deviation.is_none());
            for control in [
                crate::MeshControl::AspectRatio,
                crate::MeshControl::AngularDeviation,
            ] {
                assert!(
                    mesh.unavailable_controls
                        .iter()
                        .any(|c| c.face == face.face && c.control == control)
                );
            }
        }
        for triangle in &mesh.triangles[face.triangles.clone()] {
            let points = triangle.vertices.map(|i| mesh.positions[i]);
            let ab: [f64; 3] = std::array::from_fn(|i| points[1][i] - points[0][i]);
            let ac: [f64; 3] = std::array::from_fn(|i| points[2][i] - points[0][i]);
            collapsed |= (0..3).all(|i| {
                let cross = ab[(i + 1) % 3] * ac[(i + 2) % 3] - ab[(i + 2) % 3] * ac[(i + 1) % 3];
                cross.abs().to_bits() == 0
            });
            for i in 0..3 {
                let distance2: f64 = (0..3)
                    .map(|j| (points[i][j] - points[(i + 1) % 3][j]).powi(2))
                    .sum();
                assert!(distance2 <= 1.0);
            }
        }
        assert_eq!(face.max_aspect_ratio.is_none(), collapsed);
    }
}

#[test]
fn oblique_elliptic_cylinder_shares_circular_cap_rails() {
    use catallax_shapes::{EllipticCylinder, Plane};
    let mut model = Model::new();
    let mut dec = Decider::default();
    let p = |[x, y, z]: [i64; 3]| RPoint3::new(x.into(), y.into(), z.into()).into();
    let v = |[x, y, z]: [i64; 3]| RVector3::new(x.into(), y.into(), z.into()).into();
    let side = model.geom_mut().add_surface(Surface::EllipticCylinder(
        EllipticCylinder::new(
            p([0, 0, 0]),
            v([1, 0, 1]),
            v([0, 1, 0]),
            4.into(),
            2.into(),
            &mut dec,
        )
        .unwrap(),
    ));
    let mut rings = Vec::new();
    let mut faces = Vec::new();
    for (center, normal) in [([0, 0, 0], [0, 0, -1]), ([2, 0, 2], [0, 0, 1])] {
        let surface = model.geom_mut().add_surface(Surface::Plane(
            Plane::new(p(center), v(normal), v([1, 0, 0]), &mut dec).unwrap(),
        ));
        let curve = model.geom_mut().add_curve(Curve::Circle(
            Circle::new(p(center), v([0, 0, 1]), 4.into(), v([1, 0, 0]), &mut dec).unwrap(),
        ));
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        rings.push(edge);
        let coedge = model.add_coedge(Coedge {
            edge,
            sense: if normal[2] < 0 {
                Sense::Reversed
            } else {
                Sense::Same
            },
        });
        let ring = model.add_loop(Loop::Coedges(vec![coedge]));
        faces.push(model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![ring],
        }));
    }
    let loops = rings
        .into_iter()
        .zip([Sense::Same, Sense::Reversed])
        .map(|(edge, sense)| {
            let coedge = model.add_coedge(Coedge { edge, sense });
            model.add_loop(Loop::Coedges(vec![coedge]))
        })
        .collect();
    faces.push(model.add_face(Face {
        surface: side,
        sense: Sense::Same,
        loops,
    }));
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    let options = MeshOptions::default();
    let mesh = mesh_body(&model, body, &options, &mut dec).unwrap().value;
    assert_closed(&mesh);
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
        for triangle in &mesh.triangles[face.triangles.clone()] {
            for weights in [[1.0 / 3.0; 3], [0.2, 0.3, 0.5]] {
                let p: [f64; 3] = std::array::from_fn(|i| {
                    (0..3)
                        .map(|j| weights[j] * mesh.positions[triangle.vertices[j]][i])
                        .sum()
                });
                let error = if model.face(face.face).surface == side {
                    let major = p[1];
                    let minor = (p[2] - p[0]) / 2.0_f64.sqrt();
                    // Radial projection supplies an independent upper bound on
                    // distance to the ellipse, whose major semi-axis is two.
                    2.0 * ((major * major / 4.0 + minor * minor / 2.0).sqrt() - 1.0).abs()
                } else {
                    p[2].abs().min((p[2] - 2.0).abs())
                };
                assert!(error <= face.max_chord_error.to_f64_nearest() + 1e-12);
            }
        }
    }
}

#[test]
fn torus_trimmed_between_non_grid_latitudes_closes_across_both_seams() {
    let mut model = Model::new();
    let original = topology::torus(&mut model, 3.into(), 1.into(), &Pose::identity()).unwrap();
    let surface = model.face(model.faces_of_body(original)[0]).surface;
    let mut dec = Decider::default();
    let edges: Vec<_> = [18, 12]
        .into_iter()
        .map(|radius| {
            let radius = Rational::from_ratio(radius, 5).unwrap();
            let curve = model.geom_mut().add_curve(Curve::Circle(
                Circle::new(
                    RPoint3::new(0.into(), 0.into(), Rational::from_ratio(4, 5).unwrap()).into(),
                    RVector3::new(0.into(), 0.into(), 1.into()).into(),
                    (&radius * &radius).into(),
                    RVector3::new(3.into(), 4.into(), 0.into()).into(),
                    &mut dec,
                )
                .unwrap(),
            ));
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            })
        })
        .collect();
    let faces: Vec<_> = [
        [Sense::Same, Sense::Reversed],
        [Sense::Reversed, Sense::Same],
    ]
    .into_iter()
    .map(|senses| {
        let loops = edges
            .iter()
            .zip(senses)
            .map(|(&edge, sense)| {
                let coedge = model.add_coedge(Coedge { edge, sense });
                model.add_loop(Loop::Coedges(vec![coedge]))
            })
            .collect();
        model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops,
        })
    })
    .collect();
    let upper = faces[0];
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    let options = MeshOptions::default();
    let mesh = mesh_body(&model, body, &options, &mut dec).unwrap().value;
    assert_closed(&mesh);
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
        for triangle in &mesh.triangles[face.triangles.clone()] {
            let z = triangle
                .vertices
                .iter()
                .map(|&i| mesh.positions[i][2])
                .sum::<f64>()
                / 3.0;
            if face.face == upper {
                assert!(z >= 0.8 - 1e-12);
            } else {
                assert!(z <= 0.8 + 1e-12);
            }
        }
    }
}

#[test]
fn strongly_curved_shared_rails_do_not_emit_boundary_only_ears() {
    let mut model = Model::new();
    let a = topology::cylinder(&mut model, 2.into(), 8.into(), &pose(0, 0, -4)).unwrap();
    let rotation = Rotation::from_quaternion(1.into(), 1.into(), 0.into(), 0.into()).unwrap();
    let placed = Pose::new(
        rotation,
        RVector3::new(0.into(), 4.into(), 0.into()).into(),
        &mut Decider::default(),
    )
    .unwrap();
    let b = topology::cylinder(&mut model, 2.into(), 8.into(), &placed).unwrap();
    let body = boolean(&mut model, a, b, BoolOp::Intersect, &mut Decider::default())
        .unwrap()
        .body;
    let options = MeshOptions::default();
    let mesh = mesh_body(&model, body, &options, &mut Decider::default())
        .unwrap()
        .value;
    assert_closed(&mesh);
    assert_accuracy(&model, &mesh);
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
    }
    for edge in &mesh.edges {
        let rail: std::collections::BTreeSet<_> = edge.vertices.iter().copied().collect();
        for triangle in &mesh.triangles {
            assert!(
                !triangle.vertices.iter().all(|vertex| rail.contains(vertex)),
                "one rail must not become a zero-volume face ear"
            );
        }
    }
}
