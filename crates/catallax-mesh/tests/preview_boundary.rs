//! SolidDraft flange previews retain their complete shared boundary.

use std::collections::{BTreeMap, BTreeSet};

use catallax_mesh::{MeshOptions, PreviewOptions, Triangle, mesh_body, preview_mesh_body};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::Surface;
use catallax_topology::Sense;

#[test]
fn rounded_reference_plane_flange_has_an_early_preview() {
    flange(include_str!("fixtures/rounded-flange.json"));
}

#[test]
fn exact_face_attachment_flange_has_an_early_preview() {
    flange(include_str!("fixtures/attached-flange.json"));
}

fn flange(source: &str) {
    let recipe = Recipe::parse(source).unwrap();
    let (built, work) = with_operation_budget(100_000_000, || {
        catallax_features::evaluate(&recipe, &mut Decider::default())
    });
    let built = built.unwrap();
    eprintln!("{}: build work={}", recipe.name, work.used);
    assert!(!work.exhausted);
    let options = MeshOptions {
        chord_tolerance: "0.1".parse().unwrap(),
        ..MeshOptions::default()
    };
    let start = std::time::Instant::now();
    let (preview, work) = with_operation_budget(100_000_000, || {
        preview_mesh_body(
            &built.model,
            built.result,
            &PreviewOptions::from(&options),
            &mut Decider::default(),
        )
    });
    eprintln!("preview work: {}, time: {:?}", work.used, start.elapsed());
    let preview = preview.unwrap();
    assert!(!work.exhausted);
    assert!(!preview.triangles.is_empty());
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    let mut unique = BTreeSet::new();
    for triangle in &preview.triangles {
        let mut vertices = triangle.vertices;
        vertices.sort_unstable();
        assert!(unique.insert(vertices), "duplicate triangle {triangle:?}");
        let cross = triangle_cross(&preview.positions, triangle);
        assert!(cross.iter().any(|&x| x != 0.0), "zero-area triangle");
        let face = built.model.face(triangle.face);
        let Surface::Plane(plane) = built.model.geom().surface(face.surface) else {
            panic!("flange fixture must have planar faces")
        };
        let normal = [&plane.normal.x, &plane.normal.y, &plane.normal.z].map(|x| x.to_f64());
        let alignment = (0..3).map(|i| cross[i] * normal[i]).sum::<f64>();
        assert!(if face.sense == Sense::Same {
            alignment > 0.0
        } else {
            alignment < 0.0
        });
        let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z].map(|x| x.to_f64());
        let length = normal.iter().map(|x| x * x).sum::<f64>().sqrt();
        for &v in &triangle.vertices {
            let p = preview.positions[v];
            let distance = (0..3)
                .map(|i| (p[i] - origin[i]) * normal[i])
                .sum::<f64>()
                .abs()
                / length;
            assert!(distance <= options.chord_tolerance.to_f64_nearest());
            assert!(p.iter().all(|x| x.is_finite()));
            assert!(
                p[1] <= 10.0 + 1e-10 || p[2] <= 10.0 + 1e-10,
                "filled the bracket recess"
            );
        }
        for i in 0..3 {
            let a = triangle.vertices[i];
            let b = triangle.vertices[(i + 1) % 3];
            let entry = incidence.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }
    let bad: Vec<_> = incidence
        .iter()
        .filter(|(_, uses)| **uses != (2, 0))
        .collect();
    assert!(bad.is_empty(), "open or overlapping facets: {bad:?}");
    for edge in &preview.edges {
        for pair in edge.vertices.windows(2).filter(|pair| pair[0] != pair[1]) {
            assert_eq!(
                incidence.get(&(pair[0].min(pair[1]), pair[0].max(pair[1]))),
                Some(&(2, 0)),
                "unshared edge {:?}",
                edge.edge
            );
        }
    }
    let certified = mesh_body(
        &built.model,
        built.result,
        &options,
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert_eq!(
        preview
            .faces
            .iter()
            .map(|f| f.face)
            .collect::<BTreeSet<_>>(),
        certified.faces.iter().map(|f| f.face).collect()
    );
    assert_eq!(
        preview
            .edges
            .iter()
            .map(|e| e.edge)
            .collect::<BTreeSet<_>>(),
        certified.edges.iter().map(|e| e.edge).collect()
    );
    for face in &preview.faces {
        let reference = certified
            .faces
            .iter()
            .find(|f| f.face == face.face)
            .unwrap();
        let measured = area(
            &preview.positions,
            &preview.triangles[face.triangles.clone()],
        );
        let expected = area(
            &certified.positions,
            &certified.triangles[reference.triangles.clone()],
        );
        assert!(
            (measured - expected).abs() < 1e-8,
            "face {:?}: area {measured} vs {expected}",
            face.face
        );
    }
    for edge in &preview.edges {
        let reference = certified
            .edges
            .iter()
            .find(|e| e.edge == edge.edge)
            .unwrap();
        let a: Vec<_> = edge
            .vertices
            .iter()
            .map(|&i| preview.positions[i])
            .collect();
        let b: Vec<_> = reference
            .vertices
            .iter()
            .map(|&i| certified.positions[i])
            .collect();
        for (source, target) in [(&a, &b), (&b, &a)] {
            for &point in source {
                let distance = target
                    .windows(2)
                    .map(|pair| segment_distance(point, pair[0], pair[1]))
                    .fold(f64::INFINITY, f64::min);
                assert!(
                    distance <= options.chord_tolerance.to_f64_nearest(),
                    "edge {:?}: deviation {distance}",
                    edge.edge
                );
            }
        }
    }
    eprintln!(
        "preview triangles={}, certified triangles={}",
        preview.triangles.len(),
        certified.triangles.len()
    );
}

fn triangle_cross(positions: &[[f64; 3]], triangle: &Triangle) -> [f64; 3] {
    let [a, b, c] = triangle.vertices.map(|i| positions[i]);
    let u: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
    let v: [f64; 3] = std::array::from_fn(|i| c[i] - a[i]);
    std::array::from_fn(|i| u[(i + 1) % 3] * v[(i + 2) % 3] - u[(i + 2) % 3] * v[(i + 1) % 3])
}

fn area(positions: &[[f64; 3]], triangles: &[Triangle]) -> f64 {
    triangles
        .iter()
        .map(|t| {
            triangle_cross(positions, t)
                .iter()
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt()
                * 0.5
        })
        .sum()
}

fn segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let d: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
    let norm2 = d.iter().map(|x| x * x).sum::<f64>();
    let t = if norm2 == 0.0 {
        0.0
    } else {
        ((0..3).map(|i| (p[i] - a[i]) * d[i]).sum::<f64>() / norm2).clamp(0.0, 1.0)
    };
    (0..3)
        .map(|i| (p[i] - a[i] - t * d[i]).powi(2))
        .sum::<f64>()
        .sqrt()
}
