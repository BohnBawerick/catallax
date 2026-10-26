//! Fine drawing meshes retain geometric bounds, shared rails and deterministic output.
use std::collections::BTreeMap;

use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::Surface;
use catallax_topology::{Model, Sense};

fn check_geometry(model: &Model, mesh: &Mesh, options: &MeshOptions) {
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    for triangle in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let uses = incidence.entry((a.min(b), a.max(b))).or_default();
            uses.0 += 1;
            uses.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(incidence.values().all(|&uses| uses == (2, 0)));
    assert!(mesh.unavailable_controls.is_empty());
    for edge in &mesh.edges {
        assert!(edge.max_chord_error <= options.chord_tolerance);
        for pair in edge.vertices.windows(2) {
            assert_eq!(
                incidence[&(pair[0].min(pair[1]), pair[0].max(pair[1]))],
                (2, 0)
            );
        }
    }
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
        assert_eq!(face.max_angular_deviation, options.angular_deviation);
        let definition = model.face(face.face);
        assert!(!model.face_origins(face.face).is_empty());
        let surface = model.geom().surface(definition.surface);
        let dec = Decider::default();
        let scalar = |value: &catallax_numbers::Real| {
            dec.enclose_at(value, 256).unwrap().mid().to_f64_nearest()
        };
        let [first, second] = match surface {
            Surface::Plane(_) => [0.0, 0.0],
            Surface::Cylinder(s) => [scalar(&s.radius()), 0.0],
            Surface::Cone(s) => [scalar(&s.slope), 0.0],
            Surface::Torus(s) => [scalar(&s.major), scalar(&s.minor)],
            _ => unreachable!(),
        };
        let (origin, axis) = match surface {
            Surface::Plane(s) => (&s.origin, &s.normal),
            Surface::Cylinder(s) => (&s.origin, &s.axis),
            Surface::Cone(s) => (&s.apex, &s.axis),
            Surface::Torus(s) => (&s.center, &s.axis),
            _ => panic!("unexpected reproduction surface"),
        };
        let origin = [&origin.x, &origin.y, &origin.z].map(scalar);
        let mut axis = [&axis.x, &axis.y, &axis.z].map(scalar);
        let norm = axis.iter().map(|x| x * x).sum::<f64>().sqrt();
        axis = axis.map(|x| x / norm);
        for triangle in &mesh.triangles[face.triangles.clone()] {
            assert_eq!(triangle.face, face.face);
            let points = triangle.vertices.map(|i| mesh.positions[i]);
            let ab: [f64; 3] = std::array::from_fn(|i| points[1][i] - points[0][i]);
            let ac: [f64; 3] = std::array::from_fn(|i| points[2][i] - points[0][i]);
            let cross: [f64; 3] = std::array::from_fn(|i| {
                ab[(i + 1) % 3] * ac[(i + 2) % 3] - ab[(i + 2) % 3] * ac[(i + 1) % 3]
            });
            let norm = cross.iter().map(|x| x * x).sum::<f64>().sqrt();
            assert!(norm > 0.0);
            for weights in [[0.5, 0.5, 0.0], [1.0 / 3.0; 3], [0.2, 0.3, 0.5]] {
                let d: [f64; 3] = std::array::from_fn(|i| {
                    (0..3).map(|j| weights[j] * points[j][i]).sum::<f64>() - origin[i]
                });
                let z = (0..3).map(|i| d[i] * axis[i]).sum::<f64>();
                let radial: [f64; 3] = std::array::from_fn(|i| d[i] - z * axis[i]);
                let r = radial.iter().map(|x| x * x).sum::<f64>().sqrt();
                let (distance, normal) = match surface {
                    Surface::Plane(_) => (z.abs(), axis),
                    Surface::Cylinder(_) => ((r - first).abs(), radial.map(|x| x / r)),
                    Surface::Cone(_) => {
                        let k = first;
                        (
                            (r - k * z).abs() / (1.0 + k * k).sqrt(),
                            std::array::from_fn(|i| radial[i] / r - k * axis[i]),
                        )
                    }
                    Surface::Torus(_) => {
                        let dr = r - first;
                        (
                            (dr.hypot(z) - second).abs(),
                            std::array::from_fn(|i| radial[i] * dr / r + axis[i] * z),
                        )
                    }
                    _ => unreachable!(),
                };
                assert!(
                    distance <= face.max_chord_error.to_f64_nearest() + 1e-10,
                    "face {}: measured {distance}, bound {}",
                    face.face.index(),
                    face.max_chord_error
                );
                let analytic_norm = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
                let sign = if definition.sense == Sense::Same {
                    1.0
                } else {
                    -1.0
                };
                let cosine = sign * (0..3).map(|i| cross[i] * normal[i]).sum::<f64>()
                    / (norm * analytic_norm);
                assert!(
                    cosine.clamp(-1.0, 1.0).acos()
                        <= options.angular_deviation.as_ref().unwrap().to_f64_nearest() + 1e-8
                );
            }
        }
    }
}

#[test]
fn fine_drawing_meshes_preserve_bounds_and_shared_edges() {
    for (input, chord) in [
        (
            include_str!("../benches/fine_mesh/countersink.json"),
            "0.000722841614740048",
        ),
        (
            include_str!("../benches/fine_mesh/fillet.json"),
            "0.0005462600113499065",
        ),
        (
            include_str!("../benches/fine_mesh/drill-point.json"),
            "0.0004358898943540674",
        ),
        (include_str!("../../../parts/P2-collar.json"), "0.00036"),
        (
            include_str!("../../../parts/P3-filleted-boss.json"),
            "0.0005462600113499065",
        ),
        (
            include_str!("../benches/fine_mesh/bracket.json"),
            "0.0016638509548634457",
        ),
        (
            include_str!("../benches/fine_mesh/trim.json"),
            "0.00011357816691600547",
        ),
    ] {
        let recipe = Recipe::parse(input).unwrap();
        let mut dec = Decider::default();
        let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        let options = MeshOptions {
            chord_tolerance: chord.parse().unwrap(),
            angular_deviation: Some(Rational::parse("0.03490658503988659").unwrap()),
            ..MeshOptions::default()
        };
        let mesh = mesh_body(&built.model, built.result, &options, &mut dec)
            .unwrap()
            .value;
        assert_eq!(
            mesh.faces
                .iter()
                .map(|face| face.face)
                .collect::<std::collections::BTreeSet<_>>(),
            built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .collect()
        );
        check_geometry(&built.model, &mesh, &options);
        assert_eq!(
            mesh,
            mesh_body(
                &built.model,
                built.result,
                &options,
                &mut Decider::default()
            )
            .unwrap()
            .value
        );
    }
}

#[test]
fn tilted_trim_preserves_fine_planar_certificates() {
    let mut input: serde_json::Value =
        serde_json::from_str(include_str!("../benches/fine_mesh/trim.json")).unwrap();
    input["steps"]
        .as_array_mut()
        .unwrap()
        .push(serde_json::json!({
            "id":"tilted", "op":"mirror", "body":"solid",
            "plane":{"point":["0","0","0"],"normal":["1","2","3"]},
            "keep_original":false
        }));
    input["result"] = serde_json::json!("tilted");
    let recipe = Recipe::parse(&input.to_string()).unwrap();
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    let options = MeshOptions {
        chord_tolerance: "0.00011357816691600547".parse().unwrap(),
        angular_deviation: Some("0.03490658503988659".parse().unwrap()),
        ..MeshOptions::default()
    };
    let mesh = mesh_body(
        &built.model,
        built.result,
        &options,
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    check_geometry(&built.model, &mesh, &options);
}

#[test]
fn rotated_cones_and_apices_preserve_bounds_without_increasing_baseline_vertices() {
    // Compact output counts measured on d0cea699 at these exact controls and
    // placements. They are not the temporary-chart MeshOptions vertex budget.
    for (base, top, height, baseline_vertices) in
        [("12", "11", "1", 269_308), ("2.5", "0", "1.5", 67_059)]
    {
        let input = serde_json::json!({
            "format":"catallax-recipe/3", "name":"Rotated fine cone",
            "steps":[{"id":"cone", "op":"cone", "base_radius":base, "top_radius":top,
                "height":height, "pose":{"rotation":["1","2","3","4"],
                    "translation":["7","-3","2"]}}],
            "result":"cone"
        });
        let recipe = Recipe::parse(&input.to_string()).unwrap();
        let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
        let options = MeshOptions {
            chord_tolerance: "0.00036".parse().unwrap(),
            angular_deviation: Some("0.03490658503988659".parse().unwrap()),
            ..MeshOptions::default()
        };
        let mesh = mesh_body(
            &built.model,
            built.result,
            &options,
            &mut Decider::default(),
        )
        .unwrap()
        .value;
        assert!(
            mesh.positions.len() <= baseline_vertices,
            "radii {base}/{top}: {} vertices exceed baseline {baseline_vertices}",
            mesh.positions.len()
        );
        check_geometry(&built.model, &mesh, &options);
    }
}
