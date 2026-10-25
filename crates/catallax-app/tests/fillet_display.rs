//! Shared-rail display closure through the application's exact recipe entry.

use std::collections::{BTreeMap, BTreeSet};

use catallax_mesh::MeshOptions;
use catallax_numbers::{BudgetPolicy, Decider};
use catallax_recipe::{Operation, Recipe, Role};

#[test]
fn p3_fillet_display_closes_shared_rails_before_and_after_a_radius_edit() {
    let mut recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let mut previous_names = None;
    for radius in ["1", "3/2"] {
        if let Operation::Fillet { radius: value, .. } = &mut recipe.steps[4].op {
            *value = serde_json::from_value(serde_json::json!(radius)).unwrap();
        }
        let options = MeshOptions::default();
        let mut dec = Decider::new(BudgetPolicy {
            max_work: 1_000_000,
            ..Default::default()
        });
        let body = catallax_app::build(&recipe, &options, &mut dec).unwrap();
        let mesh = &body.display.value;
        let expected: BTreeSet<_> = body
            .built
            .model
            .faces_of_body(body.built.result)
            .into_iter()
            .collect();
        let actual: BTreeSet<_> = mesh.triangles.iter().map(|t| t.face).collect();
        assert_eq!(actual, expected);
        let mut incidence = BTreeMap::<_, (usize, i32)>::new();
        for triangle in &mesh.triangles {
            assert!(body.names.faces.contains_key(&triangle.face));
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                uses.0 += 1;
                uses.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(incidence.values().all(|uses| *uses == (2, 0)));
        for edge in &mesh.edges {
            assert!(edge.max_chord_error <= options.chord_tolerance);
            for pair in edge.vertices.windows(2) {
                assert_eq!(
                    incidence.get(&(pair[0].min(pair[1]), pair[0].max(pair[1]))),
                    Some(&(2, 0))
                );
            }
        }
        for face in &mesh.faces {
            assert!(face.max_chord_error <= options.chord_tolerance);
        }
        let mut names: Vec<_> = body
            .names
            .faces
            .values()
            .flatten()
            .filter(|origin| matches!(origin.role, Role::Blend { .. }))
            .cloned()
            .collect();
        names.sort();
        assert_eq!(names.len(), 5);
        if let Some(previous) = &previous_names {
            assert_eq!(&names, previous);
        }
        previous_names = Some(names);
    }
}

#[test]
fn thin_boss_rim_fillets_mesh_with_default_accuracy_and_bounded_work() {
    use std::time::{Duration, Instant};

    use catallax_numbers::with_operation_budget;
    use serde_json::json;

    for (radius, tolerance) in [
        ("1", "0.01"),
        ("0.1", "0.01"),
        ("0.05", "0.01"),
        ("0.01", "0.01"),
        ("1", "0.1"),
        ("0.1", "0.1"),
        ("0.05", "0.1"),
        ("0.01", "0.1"),
        ("0.05", "0.5"),
        ("0.01", "0.5"),
    ] {
        for base in ["block", "boss"] {
            let mut input: serde_json::Value =
                serde_json::from_str(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
            input["steps"].as_array_mut().unwrap().truncate(3);
            input["steps"].as_array_mut().unwrap().push(json!({
                "id":"rim", "op":"fillet", "body":"joined", "radius":radius,
                "edges":[{"between":[
                    {"step":base,"role":"end"},
                    {"step":"boss","role":{"side":{"segment":"wall"}}}
                ]}]
            }));
            input["result"] = json!("rim");
            if radius == "0.05" {
                input["steps"].as_array_mut().unwrap().push(json!({
                    "id":"tilted", "op":"mirror", "body":"rim",
                    "plane":{"point":["0","0","0"],"normal":["1","2","3"]},
                    "keep_original":false
                }));
                input["result"] = json!("tilted");
            }
            let recipe = Recipe::parse(&input.to_string()).unwrap();
            let mut dec = Decider::default();
            let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
            let options = MeshOptions {
                chord_tolerance: catallax_numbers::Rational::parse(tolerance).unwrap(),
                angular_deviation: Some(
                    catallax_numbers::Rational::parse("0.2617993877991494").unwrap(),
                ),
                ..MeshOptions::default()
            };
            let start = Instant::now();
            dec.set_deadline(Some(start + Duration::from_secs(5)));
            let (answer, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
                catallax_mesh::mesh_body(&built.model, built.result, &options, &mut dec)
            });
            let elapsed = start.elapsed();
            eprintln!(
                "rim base={base} radius={radius} tolerance={tolerance} mesh={elapsed:?} work={work:?}"
            );
            let mesh = answer.unwrap().value;
            assert!(!work.exhausted);
            assert!(
                elapsed < Duration::from_secs(5),
                "radius={radius}: {elapsed:?}"
            );
            assert!(mesh.positions.len() < 10_000, "radius={radius}");
            let expected: BTreeSet<_> = built
                .model
                .faces_of_body(built.result)
                .into_iter()
                .collect();
            let actual: BTreeSet<_> = mesh
                .triangles
                .iter()
                .map(|triangle| triangle.face)
                .collect();
            assert_eq!(actual, expected);
            assert!(
                mesh.faces
                    .iter()
                    .all(|face| face.max_chord_error <= options.chord_tolerance)
            );
            assert!(mesh.faces.iter().all(|face| {
                face.max_angular_deviation
                    .as_ref()
                    .is_some_and(|bound| bound <= options.angular_deviation.as_ref().unwrap())
            }));
            assert!(
                mesh.edges
                    .iter()
                    .all(|edge| edge.max_chord_error <= options.chord_tolerance)
            );
            let mut incidence = BTreeMap::<_, (usize, i32)>::new();
            for triangle in &mesh.triangles {
                for i in 0..3 {
                    let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                    let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                    uses.0 += 1;
                    uses.1 += if a < b { 1 } else { -1 };
                }
            }
            assert!(incidence.values().all(|uses| *uses == (2, 0)));
        }
    }
}
